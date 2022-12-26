module SolidMixtureMod

    use kind_parameters, only: rkind,clen
    use constants,       only: zero,epssmall,eps,one,two,third,half
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use LADMod,          only: ladobject
    use exits,           only: GracefulExit
    use EOSMod,          only: eos
    use StiffGasEOS,     only: stiffgas
    use Sep1SolidEOS,    only: sep1solid
    use SolidMod,        only: solid

    implicit none

    type :: solid_mixture

        integer :: ns
        integer :: nxp, nyp, nzp
        type(solid), dimension(:), allocatable :: material

        type(decomp_info), pointer :: decomp
        type(derivatives), pointer :: der
        type(filters),     pointer :: fil
        type(filters),     pointer :: gfil
        type(ladobject),   pointer :: LAD

        real(rkind), allocatable, dimension(:,:,:) :: rho0mix, mumix, yieldmix, solidVF

        logical :: SOSmodel = .FALSE.           ! is sound speed given by `equilibrium' model? Alternative is `frozen' model. Check Saurel et al., JCP 2009.
        logical :: PTeqb = .TRUE., pEqb = .FALSE., pRelax = .FALSE., updateEtot = .FALSE.
        logical :: use_gTg = .FALSE., useOneG = .FALSE.

    contains

        procedure :: init
        procedure :: set_material
        procedure :: relaxPressure
        procedure :: relaxPressure_os
        procedure :: equilibratePressure
        procedure :: equilibratePressureTemperature
        procedure :: getLAD
        procedure :: calculate_source
        procedure :: update_g
        procedure :: update_Ys
        procedure :: update_eh
        procedure :: update_VF
        procedure :: filter
        procedure :: get_rho
        procedure :: get_primitive
        procedure :: get_conserved
        procedure :: get_ehydro_from_p
        procedure :: get_p_from_ehydro
        procedure :: get_rhoYs_from_gVF
        procedure :: get_emix
        procedure :: get_pmix
        procedure :: get_Tmix
        procedure :: getSOS
        procedure :: get_J
        procedure :: get_q
        procedure :: get_qmix
        procedure :: get_dt
        procedure :: get_eelastic_devstress
        procedure :: get_mixture_properties
        procedure :: checkNaN
        procedure :: fnumden
        procedure :: rootfind_nr_1d
        final     :: destroy

    end type

    !interface solid_mixture
    !    module procedure init
    !end interface

contains

    !function init(decomp,der,fil,LAD,ns) result(this)
    subroutine init(this,decomp,der,fil,gfil,LAD,ns,PTeqb,pEqb,pRelax,SOSmodel,use_gTg,updateEtot,useOneG)
        !type(solid_mixture)      , intent(inout) :: this
        class(solid_mixture)      , intent(inout) :: this
        type(decomp_info), target, intent(in)    :: decomp
        type(filters),     target, intent(in)    :: fil, gfil
        type(derivatives), target, intent(in)    :: der
        type(ladobject),   target, intent(in)    :: LAD
        integer,                   intent(in)    :: ns
        logical,                   intent(in)    :: PTeqb,pEqb,pRelax,updateEtot
        logical,                   intent(in)    :: SOSmodel
        logical,                   intent(in)    :: use_gTg, useOneG

        type(solid), allocatable :: dummy
        integer :: i

        if (ns < 1) call GracefulExit("Must have at least 1 species in the problem. Check input file for errors",3457)

        this%PTeqb      = PTeqb
        this%pEqb       = pEqb
        this%pRelax     = pRelax
        this%SOSmodel   = SOSmodel
        this%use_gTg    = use_gTg
        this%updateEtot = updateEtot
        this%useOneG    = useOneG

        this%ns = ns

        this%nxp = decomp%ysz(1)
        this%nyp = decomp%ysz(2)
        this%nzp = decomp%ysz(3)

        this%decomp => decomp
        this%der  => der
        this%fil  => fil
        this%gfil => gfil
        this%LAD  => LAD

        ! Allocate array of solid objects (Use a dummy to avoid memory leaks)
        allocate(dummy)
        call dummy%init(decomp,der,fil,gfil,this%PTeqb,this%pEqb,this%pRelax,this%use_gTg,this%updateEtot)

        if (allocated(this%material)) deallocate(this%material)
        allocate(this%material(this%ns))!, source=dummy)
        do i=1,this%ns
            call this%material(i)%init(decomp,der,fil,gfil,this%PTeqb,this%pEqb,this%pRelax,this%use_gTg,this%updateEtot)
        end do
        deallocate(dummy)

        this%material(1)%Ys = one
        this%material(1)%VF = one
        do i=2,this%ns
            this%material(i)%Ys = zero
            this%material(i)%VF = zero
        end do

        if(allocated(this%rho0mix)) deallocate(this%rho0mix)
        allocate(this%rho0mix(this%nxp, this%nyp, this%nzp))

        if(allocated(this%mumix)) deallocate(this%mumix)
        allocate(this%mumix(this%nxp, this%nyp, this%nzp))

        if(allocated(this%yieldmix)) deallocate(this%yieldmix)
        allocate(this%yieldmix(this%nxp, this%nyp, this%nzp))

        if(allocated(this%solidVF)) deallocate(this%solidVF)
        allocate(this%solidVF(this%nxp, this%nyp, this%nzp))

    end subroutine
    !end function

    pure elemental subroutine destroy(this)
        type(solid_mixture), intent(inout)  :: this

        if(allocated(this%solidVF)) deallocate(this%solidVF)
        if(allocated(this%yieldmix)) deallocate(this%yieldmix)
        if(allocated(this%mumix)) deallocate(this%mumix)
        if(allocated(this%rho0mix)) deallocate(this%rho0mix)

        ! Deallocate array of solids (Destructor of solid should take care of everything else)
        if (allocated(this%material)) deallocate(this%material)

        nullify(this%LAD)
        nullify(this%fil)
        nullify(this%gfil)
        nullify(this%der)
        nullify(this%decomp)

    end subroutine

    subroutine set_material(this, imat, hydro, elastic)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: imat
        class(stiffgas ),     intent(in)    :: hydro
        class(sep1solid),     intent(in)    :: elastic

        if ((imat .GT. this%ns) .OR. (imat .LE. 0)) call GracefulExit("Cannot set material with index greater than the number of species.",4534)

        if (allocated(this%material(imat)%hydro)) deallocate(this%material(imat)%hydro)
        allocate( this%material(imat)%hydro, source=hydro )
        
        if (allocated(this%material(imat)%elastic)) deallocate(this%material(imat)%elastic)
        allocate( this%material(imat)%elastic, source=elastic )
    end subroutine

    subroutine relaxPressure_os(this,rho,u,v,w,mixE,dtsim,mixP)
        class(solid_mixture), intent(inout), target :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho,u,v,w,mixE
        real(rkind),                                        intent(in)  :: dtsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: Fsrc,dt
        type(solid), pointer :: mat1, mat2
        integer :: i, tt 
        real(rkind) :: muc = one

        mat1 => this%material(1); mat2 => this%material(2)


        do tt = 1, 1000
            !print *, tt, maxval(abs(mat1%p-mat2%p))
            Fsrc = abs(mat1%VF*mat2%VF)*(mat1%p - mat2%p)/muc
            !if(maxval(abs(Fsrc))*dtsim < 1.0d-10) exit
            if(maxval(abs((mat1%p-mat2%p)/(mat1%p+mat2%p))) < 1.0d-10) exit
            dt = 0.00000002_rkind*(mat1%VF/Fsrc)

            !do i = 1, 9
            !    mat1%g(:,:,:,i) = mat1%g(:,:,:,i) * (one - dt*Fsrc/mat1%VF)
            !    mat2%g(:,:,:,i) = mat2%g(:,:,:,i) * (one - dt*Fsrc/mat2%VF)
            !enddo

            mat1%VF = mat1%VF + Fsrc*dt
            mat2%VF = mat2%VF - Fsrc*dt
            !print *, tt, maxval((mat1%p*Fsrc*dt)), minval(mat1%p*Fsrc*dt)
            !print *, "p1 before = ", mat1%p(100,1,1)

            mat1%consrv(:,:,:,2) = mat1%consrv(:,:,:,2) - dt*mat1%p*Fsrc
            mat2%consrv(:,:,:,2) = mat2%consrv(:,:,:,2) + dt*mat1%p*Fsrc

            call mat1%get_primitive(rho,u,v,w)
            call mat2%get_primitive(rho,u,v,w)
            call this%get_p_from_ehydro(rho)   ! Get species pressures from species hydrodynamic energy 
            !print *, "p1 after = ", mat1%p(100,1,1)

            if (tt==1000) then
                print *, "Pressure relaxation OS did not converge", maxval(abs(Fsrc)), dtsim
            endif
        enddo

        call this%get_pmix(mixP)              ! Get mixture pressure from species pressures
stop
    end subroutine

    subroutine relaxPressure(this,rho,mixE,mixP)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho, mixE
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix
        real(rkind), dimension(4*this%ns), target :: fparams
        integer, dimension(1)             :: iparams
        real(rkind), dimension(:), pointer :: vf, gam, psph, pinf

        integer :: i,j,k,imat,minlc(1),thisind,othrind
        real(rkind) :: refp, peqb, peqb2, minvf

        real(rkind), dimension(1:this%ns) :: fac

        ehmix = mixE
        do imat = 1, this%ns
            ehmix = ehmix - this%material(imat)%Ys * this%material(imat)%eel
        enddo

        ! equilibrate and reset species pressures, reset volume fractions

        vf   => fparams(  1:this%ns)
        gam  => fparams(  this%ns+1:2*this%ns)
        psph => fparams(2*this%ns+1:3*this%ns)
        pinf => fparams(3*this%ns+1:4*this%ns)
        
        do imat=1,this%ns
          fparams(  this%ns+imat) = this%material(imat)%hydro%gam    ! gamma
          fparams(3*this%ns+imat) = this%material(imat)%hydro%PInf   ! PInf
        enddo

        ! Set reference pressure
        refp = maxval(fparams(3*this%ns+1:4*this%ns)) ! max over all PInfs
        if(refp < 1.0D-5) refp = 1.0D0

        do k=1,this%nzp
         do j=1,this%nyp
          do i=1,this%nxp

            do imat=1,this%ns
                ! set fparams
                fparams(          imat) = this%material(imat)%VF(i,j,k)    ! volume fractions
                fparams(2*this%ns+imat) = this%material(imat)%p(i,j,k)     ! pressure before eqb
                !print '(1(i4,1x),3(e19.12,1x))', i, fparams(2*this%ns+imat), fparams(imat), this%material(1)%VF(i,j,k) + this%material(2)%VF(i,j,k)
            end do

            !minvf = minval(vf(1:this%ns));     minlc = minloc(vf(1:this%ns))
            !if(minvf < zero) then
            !    ! do not solve non-linear problem. assume relaxed pressure is
            !    ! equal to the pressure of the dominant species
            !    thisind = minlc(1);  othrind = mod(thisind, 2) + 1
            !    psph = this%material(othrind)%p(i,j,k)
            !    do imat=1,this%ns
            !        vf(imat) = this%material(imat)%Ys(i,j,k)* this%material(imat)%eh(i,j,k)*rho(i,j,k)*(gam(imat)-one)/ &
            !                   (psph(imat)+gam(imat)*pinf(imat))
            !    end do
            !    peqb = psph(1)
            !    peqb2 = peqb
            !    
            !else
                ! solve non-linear problem for relaxed pressure

                ! set iparams
                iparams(1) = 1     !   used in fnumden; 2 for PTeqb 

                ! scale all pressures by pref
                fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)/refp

                ! set initial guess
                peqb = sum(fparams(1:this%ns)*fparams(2*this%ns+1:3*this%ns))
                ! solve non-linear equation
                call this%rootfind_nr_1d(peqb,fparams,iparams)

                ! rescale all pressures by refp
                fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)*refp
                peqb = peqb*refp

                ! update species VF, eh, ...
                fac = (psph + gam*pinf + (gam-one)*peqb)/(gam*(pinf+peqb))
                vf = vf*fac
                !this%material(1:this%ns)%g = this%material(1:this%ns)%g*fac**third        !! --- not clear if this is needed or if it works
                peqb2 = (rho(i,j,k)*ehmix(i,j,k) - sum(vf*gam*pinf/(gam-one))) / sum(vf/(gam-one))
                psph = peqb2
            !endif
                !print *, 'peqb = ', peqb
                !if(peqb2 < zero) then
                !    write(*,*) i, peqb2, peqb
                !    write(*,*) 'rho: ', rho(i,j,k)
                !    write(*,*) 'eh: ', ehmix(i,j,k)
                !    write(*,*) 'vf : ', vf
                !    write(*,*) 'gam: ', gam
                !    write(*,*) 'fac: ', fac
                !    write(*,*) 'Ys : ', this%material(1)%Ys(i,j,k), this%material(2)%Ys(i,j,k)
                !    write(*,*) 'eh : ', this%material(1)%eh(i,j,k), this%material(2)%eh(i,j,k)
                !    write(*,*) '-----------------------------------'
                !    write(*,*) '-----------------------------------'
                !endif

            mixP(i,j,k) = peqb2
            do imat=1,this%ns
              this%material(imat)%VF(i,j,k) = vf(imat)
              this%material(imat)%p(i,j,k) = psph(imat)
            enddo

            

          enddo
         enddo
        enddo
        ! get species energy from species pressure
        do imat = 1, this%ns
            call this%material(imat)%get_ehydroT_from_p(rho)
        end do

    end subroutine

    subroutine equilibratePressure(this,mixRho,mixE,mixP)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: mixRho, mixE
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix, rhom, tmp

        integer :: imat, i,j,k,a(3)
        real(rkind) :: gamfac

        ! print *, '---'
        ! subtract elastic energy to determine hydrostatic energy. Temperature
        ! is assumed a function of only hydrostatic energy. Not sure if this is
        ! correct.
        ehmix = mixE
        tmp = zero
        do imat = 1, this%ns
            ehmix = ehmix - this%material(imat)%Ys * this%material(imat)%eel

            gamfac =  this%material(imat)%hydro%gam * this%material(imat)%hydro%onebygam_m1*this%material(imat)%hydro%PInf
            call this%material(imat)%getSpeciesDensity(mixRho,rhom)
            
            ! ehmix = ehmix - gamfac*this%material(imat)%Ys/rhom
            ehmix = ehmix - gamfac*this%material(imat)%VF/mixRho
            ! tmp = tmp + this%material(imat)%hydro%onebygam_m1 * this%material(imat)%Ys/rhom
            tmp = tmp + this%material(imat)%hydro%onebygam_m1 * this%material(imat)%VF/mixRho

            !do k=1,this%nzp
            !    do j = 1,this%nyp
            !        do i = 1,this%nxp
            !            if ( this%material(imat)%Ys(i,j,k) * this%material(imat)%VF(i,j,k) < zero ) then
            !                print *, "Found negative Ys*VF in material ", imat, " at ", i, j, k
            !                print *, "  Ys = ", this%material(imat)%Ys(i,j,k)
            !                print *, "  VF = ", this%material(imat)%Vf(i,j,k)
            !                ! exit
            !            end if
            !        end do
            !    end do
            !end do

        enddo

         mixP = ehmix/tmp

         ! a = minloc(ehmix)
         ! print *, "  loc: ", minloc(ehmix)
         ! !a(1) = 122; a(2:3) = 1
         ! print *, "Mat1 Ys, VF = ", this%material(1)%Ys(a(1),a(2),a(3)), this%material(1)%VF(a(1),a(2),a(3))
         ! print *, "Mat2 Ys, VF = ", this%material(2)%Ys(a(1),a(2),a(3)), this%material(2)%VF(a(1),a(2),a(3))
         ! print *, "Mat1,2 T1   = ", this%material(1)%hydro%gam * this%material(1)%hydro%onebygam_m1 * this%material(1)%hydro%PInf   * this%material(1)%VF(a(1),a(2),a(3))/mixRho(a(1),a(2),a(3)),  this%material(2)%hydro%gam * this%material(2)%hydro%onebygam_m1 * this%material(2)%hydro%PInf   * this%material(2)%VF(a(1),a(2),a(3))/mixRho(a(1),a(2),a(3))
         ! print *, "Mat1,2 eel = ",  this%material(1)%eel(a(1),a(2),a(3)),  this%material(2)%eel(a(1),a(2),a(3))
         ! print *, "Mat1,2 g11 = ",  this%material(1)%g11(a(1),a(2),a(3)),  this%material(2)%g11(a(1),a(2),a(3))
         ! print *, "Mat1,2 g22 = ",  this%material(1)%g22(a(1),a(2),a(3)),  this%material(2)%g22(a(1),a(2),a(3))
         ! print *, "Mix eh      = ",  mixE(a(1),a(2),a(3)), mixE(a(1),a(2),a(3)) -   (this%material(1)%Ys(a(1),a(2),a(3)) * this%material(1)%eel(a(1),a(2),a(3)) + this%material(2)%Ys(a(1),a(2),a(3)) * this%material(2)%eel(a(1),a(2),a(3)) )
         ! print *, "ehmix      = ",  ehmix(a(1),a(2),a(3))
         ! print *, "tmp        = ",  tmp(a(1),a(2),a(3))
         ! print *, "mixP       = ",  mixP(a(1),a(2),a(3))
         
         !print*, "ehmix: ", minval(ehmix), maxval(ehmix)
         !print*, "tmp  : ", minval(tmp), maxval(tmp)
         !print*, "mixP : ", minval(mixP), maxval(mixP)
      
         !if (minval(mixP) < zero) then
         !    print *, "  loc: ", minloc(mixP)
         !end if

         do imat = 1, this%ns
             this%material(imat)%p = mixP
             call this%material(imat)%get_ehydroT_from_p(mixRho)
         end do

    end subroutine

    subroutine equilibratePressureTemperature(this,mixRho,mixE,mixP,mixT)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: mixRho, mixE
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP, mixT

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix
        real(rkind), dimension(4*this%ns), target :: fparams
        integer, dimension(4)             :: iparams
        ! real(rkind), dimension(:), pointer :: vf, psph, pinf

        integer :: i,j,k,imat
        real(rkind) :: maxp, peqb !, pest, pdiffmax

        ! subtract elastic energy to determine hydrostatic energy. Temperature
        ! is assumed a function of only hydrostatic energy. Not sure if this is
        ! correct.
        ehmix = mixE
        do imat = 1, this%ns
            ehmix = ehmix - this%material(imat)%Ys * this%material(imat)%eel
        enddo
!print *, 'pinpiut: ', mixE(89,1,1), ehmix(89,1,1), mixRho(89,1,1)
!print *, 'VF     : ', this%material(1)%VF(89,1,1), this%material(2)%VF(89,1,1)
!print *, 'pstart : ', this%material(1)%p(89,1,1), this%material(2)%p(89,1,1)

        do k=1,this%nzp
         do j=1,this%nyp
          do i=1,this%nxp
            ! set fparams
            fparams(1) = mixRho(i,j,k)*ehmix(i,j,k)

            ! set iparams
            iparams(1) = 2
            iparams(2) = i; iparams(3) = j; iparams(4) = k;

            maxp = zero; peqb = zero
            do imat=1,this%ns
              !! determine max over all PInfs
              !maxp = maxval(maxp, this%material(imat)%hydro%PInf)

              ! set initial guess
              peqb = peqb + this%material(imat)%VF(i,j,k)*this%material(imat)%p(i,j,k)
            end do
            !pest = peqb

            ! solve non-linear equation
            call this%rootfind_nr_1d(peqb,fparams,iparams)
            !pdiffmax = max(dabs(pest-peqb),pdiffmax)

            !! rescale all pressures by maxp
            !fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)*maxp
            mixP(i,j,k) = peqb !*maxp

          enddo
         enddo
        enddo
!print *, 'pafter: ', this%material(1)%p(89,1,1), this%material(2)%p(89,1,1)

        mixT = zero
        do i = 1, this%ns
          mixT = mixT + this%material(i)%Ys*this%material(i)%hydro%Cv * &
                (mixP + this%material(i)%hydro%gam*this%material(i)%hydro%PInf)/(mixP + this%material(i)%hydro%PInf)
        enddo
        mixT = ehmix/mixT

        do i = 1, this%ns
            this%material(i)%T = mixT
            this%material(i)%p = mixP
            this%material(i)%VF = mixRho*this%material(i)%Ys*(this%material(i)%hydro%gam-one)* &
                                         this%material(i)%hydro%Cv*mixT/(mixP + this%material(i)%hydro%PInf)

            this%material(i)%eh = this%material(i)%hydro%Cv*mixT*(mixP + this%material(i)%hydro%gam*this%material(i)%hydro%PInf) / &
                                                                 (mixP + this%material(i)%hydro%PInf)
        end do

    end subroutine

    subroutine get_dt(this, rho, delta, dtkap, dtDiff, dtplast)
        use reductions, only: P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho
        real(rkind), intent(in)  :: delta
        real(rkind), intent(out) :: dtkap, dtDiff, dtplast

        integer :: imat

        dtkap = one/epssmall; dtDiff = dtkap; dtplast = dtDiff

        do imat = 1, this%ns
          if (.NOT. this%PTeqb) then
              dtkap =   min(dtkap,   one / ( (P_MAXVAL( this%material(imat)%kap*this%material(imat)%T/(rho* delta**4)))**(third) + eps))
          end if
          dtDiff =  min(dtDiff,  one / ( (P_MAXVAL( this%material(imat)%diff/delta**2) + eps)) )
          dtplast = min(dtplast, this%material(imat)%elastic%tau0)
        enddo

        ! For now disable plastic time step limit by setting a large value
        dtplast = real(1.0D32,rkind)

    end subroutine

    ! Subroutine to get species art. conductivities and diffusivities
    subroutine getLAD(this,rho,e,sos,x_bc,y_bc,z_bc)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho,e,sos  ! Mixture density and speed of sound
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: i

        do i = 1,this%ns
            call this%material(i)%getPhysicalProperties()

            if (.NOT. this%PTeqb) then
                ! Artificial conductivity
                !call this%LAD%get_conductivity(rho, this%material(i)%eh, this%material(i)%T, sos, &
                !                                    this%material(i)%kap, x_bc, y_bc, z_bc)
                call this%LAD%get_conductivity(rho, this%material(i)%Ys*this%material(i)%eh, e, sos, &     ! --change2
                                                    this%material(i)%kap, x_bc, y_bc, z_bc)
            end if

            ! Artificial diffusivity (grad(Ys) is stored in Ji at this stage)
!print *, 'bef LAD:', this%material(1)%Ji(89,1,1,1), sos(89,1,1), this%material(1)%diff(89,1,1)
            call this%LAD%get_diffusivity(this%material(i)%Ys, this%material(i)%Ji(:,:,:,1), &
                                          this%material(i)%Ji(:,:,:,2), this%material(i)%Ji(:,:,:,3), &
                                          sos, this%material(i)%diff, x_bc, y_bc, z_bc)
!print *, 'aft LAD:', this%material(1)%Ji(89,1,1,1), sos(89,1,1), this%material(1)%diff(89,1,1)
        end do

    end subroutine

    subroutine get_p_from_ehydro(this, rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho

        integer :: imat

        ! get species pressure from species energy
        do imat = 1, this%ns
            call this%material(imat)%get_p_from_ehydro(rho)
        enddo

    end subroutine

    subroutine get_ehydro_from_p(this,rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho

        integer :: imat

        ! get species energy from species pressure
        do imat = 1, this%ns
            call this%material(imat)%get_ehydroT_from_p(rho)             ! computes species ehydro and T from species p
        enddo

    end subroutine

    subroutine get_Tmix(this,T)
        class(solid_mixture), intent(in) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: T  ! Mixture temperature

        integer :: i

        T = zero
        do i = 1,this%ns
            T = T + this%material(i)%VF * this%material(i)%T  ! Volume fraction weighted sum
        end do
    end subroutine

    subroutine get_mixture_properties(this)
        class(solid_mixture), intent(inout) :: this

        integer :: imat

        ! compute rho0mix, mumix, yieldmix as VF weighted sums
        this%rho0mix  = this%material(1)%elastic%rho0  * this%material(1)%VF
        !this%mumix    = this%material(1)%elastic%mu/this%material(1)%elastic%rho0 * this%material(1)%VF
        this%mumix    = this%material(1)%elastic%mu * this%material(1)%VF
        this%yieldmix = this%material(1)%elastic%yield * this%material(1)%VF
        do imat = 2, this%ns
            this%rho0mix  = this%rho0mix  + this%material(imat)%elastic%rho0 * this%material(imat)%VF
            !this%mumix    = this%mumix    + this%material(imat)%elastic%mu/this%material(imat)%elastic%rho0 * this%material(imat)%VF
            this%mumix    = this%mumix    + this%material(imat)%elastic%mu * this%material(imat)%VF
            this%yieldmix = this%yieldmix + this%material(imat)%elastic%yield * this%material(imat)%VF
        enddo
        !this%mumix = this%mumix * this%rho0mix

        this%solidVF = zero
        do imat = 1, this%ns
          if(this%material(imat)%elastic%mu > eps) this%solidVF = this%solidVF + this%material(imat)%VF
        enddo

    end subroutine

    subroutine get_eelastic_devstress(this,devstress)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6), intent(out) :: devstress

        integer :: imat

        devstress = zero

        if(this%useOneG) then
            call this%material(1)%get_eelastic_devstress(this%rho0mix,this%mumix)
            !print '(a,4(e21.14,1x))', 'mu/rho0_m', this%mumix(89,1,1), this%rho0mix(89,1,1), this%mumix(89,1,1)/this%rho0mix(89,1,1)
            !print '(a,4(e21.14,1x))', 'mu/rho0_1', this%material(1)%elastic%mu, this%material(1)%elastic%rho0, this%material(1)%elastic%mu/this%material(1)%elastic%rho0, this%material(1)%VF(89,1,1)
            !print '(a,4(e21.14,1x))', 'mu/rho0_2', this%material(2)%elastic%mu, this%material(2)%elastic%rho0, this%material(2)%elastic%mu/this%material(2)%elastic%rho0, this%material(2)%VF(89,1,1)
            do imat = 2, this%ns
              this%material(imat)%eel = this%material(1)%eel
              !this%material(imat)%devstress = this%material(1)%devstress
            enddo
            devstress = this%material(1)%devstress
            !print '(a,4(e21.14,1x))', 'mixt sxx', this%material(1)%devstress(89,1,1,1)
            !print '(a,4(e21.14,1x))', 'mixt eel', this%material(1)%eel(89,1,1), this%material(2)%eel(89,1,1)
            !print '(a,4(e21.14,1x))', 'mixt Ys ', this%material(1)%Ys(89,1,1), this%material(2)%Ys(89,1,1)
            !print '(a,4(e21.14,1x))', 'tot  mu ', this%mumix(89,1,1)/this%rho0mix(89,1,1), this%mumix(89,1,1), this%rho0mix(89,1,1)
            !print '(a,4(e21.14,1x))', 'eel 1   ', this%material(1)%Ys(89,1,1)*this%material(1)%eel(89,1,1)
            !print '(a,4(e21.14,1x))', 'eel 2   ', this%material(2)%Ys(89,1,1)*this%material(2)%eel(89,1,1)
            !print '(a,4(e21.14,1x))', 'tot  eel', this%material(1)%Ys(89,1,1)*this%material(1)%eel(89,1,1) + this%material(2)%Ys(89,1,1)*this%material(2)%eel(89,1,1)
            !print '(a,9(e21.14,1x))', 'mix1 g  ', this%material(1)%g(89,1,1,:)
            !print '(a,9(e21.14,1x))', 'mix2 g  ', this%material(2)%g(89,1,1,:)
        else
            do imat = 1, this%ns
              call this%material(imat)%get_eelastic_devstress()
              ! print *, "Material ", imat, " sigma:"
              ! print *, "   ",  this%material(imat)%devstress(200,1,1,1), this%material(imat)%devstress(200,1,1,2), this%material(imat)%devstress(200,1,1,3)
              ! print *, "   ",  this%material(imat)%devstress(200,1,1,2), this%material(imat)%devstress(200,1,1,4), this%material(imat)%devstress(200,1,1,5)
              ! print *, "   ",  this%material(imat)%devstress(200,1,1,3), this%material(imat)%devstress(200,1,1,5), this%material(imat)%devstress(200,1,1,6)
              devstress(:,:,:,1) = devstress(:,:,:,1) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,1)
              devstress(:,:,:,2) = devstress(:,:,:,2) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,2)
              devstress(:,:,:,3) = devstress(:,:,:,3) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,3)
              devstress(:,:,:,4) = devstress(:,:,:,4) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,4)
              devstress(:,:,:,5) = devstress(:,:,:,5) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,5)
              devstress(:,:,:,6) = devstress(:,:,:,6) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,6)
            end do
            !print '(a,4(e21.14,1x))', 'mult sxx', devstress(89,1,1,1)
            !print '(a,4(e21.14,1x))', 'mult eel', this%material(1)%eel(89,1,1), this%material(2)%eel(89,1,1)
            !print '(a,4(e21.14,1x))', 'mult Ys ', this%material(1)%Ys(89,1,1), this%material(2)%Ys(89,1,1)
            !print '(a,4(e21.14,1x))', 'tot  mu1', (this%material(1)%VF(89,1,1)*this%material(1)%elastic%mu + this%material(2)%VF(89,1,1)*this%material(2)%elastic%mu)/(this%material(1)%VF(89,1,1)*this%material(1)%elastic%rho0 + this%material(2)%VF(89,1,1)*this%material(2)%elastic%rho0), (this%material(1)%VF(89,1,1)*this%material(1)%elastic%mu + this%material(2)%VF(89,1,1)*this%material(2)%elastic%mu), (this%material(1)%VF(89,1,1)*this%material(1)%elastic%rho0 + this%material(2)%VF(89,1,1)*this%material(2)%elastic%rho0)
            !print '(a,4(e21.14,1x))', 'tot  mu2', (this%material(1)%Ys(89,1,1)*this%material(1)%elastic%mu/this%material(1)%elastic%rho0 + this%material(2)%Ys(89,1,1)*this%material(2)%elastic%mu/this%material(2)%elastic%rho0)
            !print '(a,4(e21.14,1x))', 'tot  eel', this%material(1)%Ys(89,1,1)*this%material(1)%eel(89,1,1) + this%material(2)%Ys(89,1,1)*this%material(2)%eel(89,1,1)
            !print '(a,9(e21.14,1x))', 'mat1 g  ', this%material(1)%g(89,1,1,:)
            !print '(a,9(e21.14,1x))', 'mat2 g  ', this%material(2)%g(89,1,1,:)
            !print '(a,4(e21.14,1x))', 'mat1 sxx', this%material(1)%devstress(89,1,1,1), this%material(1)%VF(89,1,1)
            !print '(a,4(e21.14,1x))', 'mat2 sxx', this%material(2)%devstress(89,1,1,1), this%material(2)%VF(89,1,1)
        endif

    end subroutine

    subroutine get_J(this,rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho  ! Mixture density

        integer :: i
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: sumJx,sumJy,sumJz
!print *, 'diff Ji:', this%material(1)%Ji(89,1,1,1), this%material(1)%diff(89,1,1)
        sumJx = zero; sumJy = zero; sumJz = zero;
        ! Get diff*gradYs (gradYs are in Ji's)
        do i=1,this%ns
            this%material(i)%Ji(:,:,:,1) = this%material(i)%Ji(:,:,:,1)*this%material(i)%diff
            sumJx = sumJx + this%material(i)%Ji(:,:,:,1)

            this%material(i)%Ji(:,:,:,2) = this%material(i)%Ji(:,:,:,2)*this%material(i)%diff
            sumJy = sumJy + this%material(i)%Ji(:,:,:,2)

            this%material(i)%Ji(:,:,:,3) = this%material(i)%Ji(:,:,:,3)*this%material(i)%diff
            sumJz = sumJz + this%material(i)%Ji(:,:,:,3)
        end do

        ! Correct Ji's so this sum becomes zero (No net diffusive flux)
        do i=1,this%ns
            this%material(i)%Ji(:,:,:,1) = -rho*( this%material(i)%Ji(:,:,:,1) - this%material(i)%Ys*sumJx )
            this%material(i)%Ji(:,:,:,2) = -rho*( this%material(i)%Ji(:,:,:,2) - this%material(i)%Ys*sumJy )
            this%material(i)%Ji(:,:,:,3) = -rho*( this%material(i)%Ji(:,:,:,3) - this%material(i)%Ys*sumJz )
        end do

    end subroutine

    pure subroutine get_conserved(this,rho,u,v,w)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%get_conserved(rho,u,v,w)
        end do

    end subroutine

    subroutine get_primitive(this,rho,u,v,w,e,devstress,p,sos)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho, u,v,w,e
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6), intent(out) :: devstress
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: p, sos

        integer :: imat
        do imat = 1, this%ns
          call this%material(imat)%get_primitive(rho,u,v,w)
        end do
  
        call this%get_eelastic_devstress(devstress)   ! Get species elastic energies, and mixture and species devstress
        if(this%ns == 1) then
          this%material(1)%eh = e - this%material(1)%eel ! Since eh equation is not exact and this is a better alternative for 1 species
        endif
        if(this%PTeqb) then    ! --- why is this condition needed here? this provides initial guess for pressure even in pRelax case -- NSG
            call this%get_p_from_ehydro(rho)   ! Get species pressures from species hydrodynamic energy 
            call this%get_pmix(p)              ! Get mixture pressure from species pressures
        endif
        
        call this%getSOS(rho,p,sos)

    end subroutine

    subroutine get_rho(this,rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rho

        integer :: imat

        rho = zero
        do imat = 1, this%ns
          rho = rho + this%material(imat)%consrv(:,:,:,1)
        end do
!print *, '--rho--', rho(89,1,1), this%material(1)%consrv(89,1,1,1), this%material(2)%consrv(89,1,1,1)
    end subroutine

    subroutine get_rhoYs_from_gVF(this,rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rho
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhom

        integer :: imat

        rho = zero
        do imat = 1, this%ns
          call this%material(imat)%getSpeciesDensity_from_g(rhom)
          rho = rho + this%material(imat)%VF * rhom
        end do

        do imat = 1,this%ns
            call this%material(imat)%getSpeciesDensity_from_g(rhom)
            this%material(imat)%Ys = this%material(imat)%VF * rhom / rho
        end do
    end subroutine

    subroutine get_q(this,x_bc,y_bc,z_bc)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y, transpose_y_to_z, transpose_z_to_y
        class(solid_mixture), intent(inout) :: this
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: i
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: tmp1_in_x, tmp2_in_x
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: tmp1_in_y
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: tmp1_in_z, tmp2_in_z

        do i = 1,this%ns
            if(.NOT. this%PTeqb) then
                ! Step 1: Get qy
                !call this%der%ddy(this%material(i)%T,tmp1_in_y,y_bc(1),y_bc(2))
                !this%material(i)%qi(:,:,:,2) = -this%material(i)%VF*this%material(i)%kap*tmp1_in_y
                call this%der%ddy(this%material(i)%Ys*this%material(i)%eh,tmp1_in_y,y_bc(1),y_bc(2)) ! --change2
                this%material(i)%qi(:,:,:,2) = -this%material(i)%kap*tmp1_in_y

                ! Step 2: Get qx
                !call transpose_y_to_x(this%material(i)%T,tmp1_in_x,this%decomp)
                !call this%der%ddx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
                !call transpose_x_to_y(tmp2_in_x,tmp1_in_y,this%decomp)
                !this%material(i)%qi(:,:,:,1) = -this%material(i)%VF*this%material(i)%kap*tmp1_in_y
                call transpose_y_to_x(this%material(i)%Ys*this%material(i)%eh,tmp1_in_x,this%decomp) ! --change2
                call this%der%ddx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
                call transpose_x_to_y(tmp2_in_x,tmp1_in_y,this%decomp)
                this%material(i)%qi(:,:,:,1) = -this%material(i)%kap*tmp1_in_y

                ! Step 3: Get qz
                !call transpose_y_to_z(this%material(i)%T,tmp1_in_z,this%decomp)
                !call this%der%ddz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
                !call transpose_z_to_y(tmp2_in_z,tmp1_in_y)
                !this%material(i)%qi(:,:,:,3) = -this%material(i)%VF*this%material(i)%kap*tmp1_in_y
                call transpose_y_to_z(this%material(i)%Ys*this%material(i)%eh,tmp1_in_z,this%decomp) ! --change2
                call this%der%ddz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
                call transpose_z_to_y(tmp2_in_z,tmp1_in_y)
                this%material(i)%qi(:,:,:,3) = -this%material(i)%kap*tmp1_in_y
            else
                this%material(i)%qi(:,:,:,:) = zero     ! calculated in sgrid
            endif

            ! If multispecies, add the inter-species enthalpy flux
            if (this%ns .GT. 1) then
                call this%material(i)%get_enthalpy(tmp1_in_y)
                this%material(i)%qi(:,:,:,1) = this%material(i)%qi(:,:,:,1) + ( tmp1_in_y * this%material(i)%Ji(:,:,:,1) )
                this%material(i)%qi(:,:,:,2) = this%material(i)%qi(:,:,:,2) + ( tmp1_in_y * this%material(i)%Ji(:,:,:,2) )
                this%material(i)%qi(:,:,:,3) = this%material(i)%qi(:,:,:,3) + ( tmp1_in_y * this%material(i)%Ji(:,:,:,3) )
            end if
        end do

        ! Done
    end subroutine

    subroutine get_qmix(this,qx,qy,qz)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: qx,qy,qz

        integer :: i

        qx = zero; qy = zero; qz = zero

        do i = 1,this%ns
            qx = qx + this%material(i)%qi(:,:,:,1)
            qy = qy + this%material(i)%qi(:,:,:,2)
            qz = qz + this%material(i)%qi(:,:,:,3)
        end do

    end subroutine

    subroutine filter(this, iflag,x_bc,y_bc,z_bc)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: iflag
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%filter(iflag,x_bc,y_bc,z_bc)
        end do
       
        if(this%ns==2) then
            this%material(2)%VF = one - this%material(1)%VF
        endif

    end subroutine

    subroutine calculate_source(this,rho,divu,u,v,w,p,Fsource,x_bc,y_bc,z_bc)
        use operators, only: divergence,gradient
        class(solid_mixture), intent(in), target :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho,divu,u,v,w,p
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: Fsource
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc
      
        type(solid), pointer :: mat1, mat2
        real(rkind), dimension(this%nxp,this%nyp,this%nzp)  :: tmp1,tmp2,tmp3,rhocsq1,rhocsq2

            mat1 => this%material(1); mat2 => this%material(2)

        if(this%pEqb) then

            !call mat1%getSpeciesDensity(rho,tmp1)
            call mat1%hydro%get_sos2(mat1%rhom,mat1%p,tmp2)
            call mat1%elastic%get_sos2(mat1%rhom,tmp2)
            rhocsq1 = mat1%rhom*tmp2

            !call mat2%getSpeciesDensity(rho,tmp1)
            call mat2%hydro%get_sos2(mat2%rhom,mat2%p,tmp3)
            call mat2%elastic%get_sos2(mat2%rhom,tmp3)
            rhocsq2 = mat2%rhom*tmp3

            Fsource = mat1%VF*tmp3*(one - mat2%rhom/mat1%rhom)

            call mat1%get_enthalpy(tmp2)
            call mat2%get_enthalpy(tmp3)

            call divergence(this%decomp,this%der,-mat1%Ji(:,:,:,1),-mat1%Ji(:,:,:,2),-mat1%Ji(:,:,:,3),tmp1,-x_bc,-y_bc,-z_bc)    ! mass fraction equation is anti-symmetric

            Fsource = -(rhocsq1 - rhocsq2)*divu*mat1%VF*mat2%VF + (Fsource + mat1%VF*(mat2%hydro%gam-one)*(tmp2-tmp3))*tmp1

            Fsource = Fsource / (mat2%VF*rhocsq1 + mat1%VF*rhocsq2)

        elseif(this%pRelax) then

            call divergence(this%decomp,this%der,mat1%Ji(:,:,:,1),mat1%Ji(:,:,:,2),mat1%Ji(:,:,:,3),tmp1,-x_bc,-y_bc,-z_bc)    ! mass fraction equation is anti-symmetric
            Fsource = -tmp1*(mat1%eh + mat1%eel + half*(u*u+v*v+w*w))

            call gradient(this%decomp,this%der,mat2%VF,tmp1,tmp2,tmp3)
            Fsource = Fsource - p*(u*tmp1 + v*tmp2 + w*tmp3)

            call gradient(this%decomp,this%der,p,tmp1,tmp2,tmp3)
            Fsource = Fsource - mat1%VF*mat2%VF/rho*(mat1%rhom-mat2%rhom)*(u*tmp1 + v*tmp2 + w*tmp3)

        endif

    end subroutine 

    subroutine update_g(this,isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc)
        use operators, only: filter3D
        use reductions, only: P_SUM
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho,u,v,w,src
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: imat

        if(this%useOneG) then
            if (this%use_gTg) then
                call this%material(1)%update_gTg(isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc,this%rho0mix,this%mumix,this%yieldmix,this%solidVF)
            else
                call this%material(1)%update_g(isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc,this%rho0mix,this%mumix,this%yieldmix,this%solidVF)
            end if
            do imat = 2, this%ns
                this%material(imat)%g = this%material(1)%g
            enddo
        else
            do imat = 1, this%ns
                if (this%use_gTg) then
                    call this%material(imat)%update_gTg(isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc)
                else
                    call this%material(imat)%update_g(isub,dt,rho,u,v,w,x,y,z,src,tsim,x_bc,y_bc,z_bc)
                end if
            end do
        endif

    end subroutine

    subroutine update_Ys(this,isub,dt,rho,u,v,w,x,y,z,tsim,x_bc,y_bc,z_bc)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%update_Ys(isub,dt,rho,u,v,w,x,y,z,tsim,x_bc,y_bc,z_bc)
        end do

    end subroutine

    subroutine get_pmix(this,p)
        class(solid_mixture), intent(in) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: p  ! Mixture pressure

        integer :: i

        p = zero
        do i = 1,this%ns
            p = p + this%material(i)%VF * this%material(i)%p  ! Volume fraction weighted sum
        end do

    end subroutine

    subroutine get_emix(this,e)
        class(solid_mixture), intent(in) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: e  ! Mixture internal energy

        integer :: i

        e = zero
        do i = 1,this%ns
            e = e + this%material(i)%Ys * ( this%material(i)%eh + this%material(i)%eel )  ! Mass fraction weighted sum
        end do

    end subroutine

    subroutine update_eh(this,isub,dt,rho,u,v,w,x,y,z,tsim,divu,viscwork,src,taustar,x_bc,y_bc,z_bc)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho,u,v,w,divu,viscwork,src
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6), intent(in) :: taustar
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: imat

        !do imat = 1, this%ns
        !  call this%material(imat)%update_eh(isub,dt,rho,u,v,w,x,y,z,tsim,divu,viscwork,src,taustar,x_bc,y_bc,z_bc)
        !end do
        call this%material(1)%update_eh(isub,dt,rho,u,v,w,x,y,z,tsim,divu,viscwork, src,taustar,x_bc,y_bc,z_bc)
        call this%material(2)%update_eh(isub,dt,rho,u,v,w,x,y,z,tsim,divu,viscwork,-src,taustar,x_bc,y_bc,z_bc)

    end subroutine

    subroutine getSOS(this,rho,p,sos)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho, p  ! Mixture density and pressure
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: sos     ! Mixture speed of sound

        integer :: i
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhom, sosm
!print *, '----In mix%getSOS----'
        sos = zero
        do i = 1,this%ns
            call this%material(i)%getSpeciesDensity(rho,rhom)
            call this%material(i)%hydro%get_sos2(rhom,p,sosm)
            !if(.not. this%useOneG) then
                call this%material(i)%elastic%get_sos2(rhom,sosm)
                !print *, 'mat  dens', i, rhom(89,1,1), p(89,1,1), this%material(i)%Ys(89,1,1), sosm(89,1,1)
            !endif
            if(this%SOSmodel) then
                ! equilibrium model
                sos = sos + this%material(i)%VF/(rhom*sosm)
            else
                ! frozen model (details in Saurel et al. 2009)
                sos = sos + this%material(i)%Ys*sosm
            endif
        end do

        if(this%SOSmodel) then
            sos = one / (sqrt(rho*sos) + epssmall)
        else
            sos = sqrt(abs(sos))
        endif

        !if(this%useOneG) then
        !    sosm = sos*sos
        !    call this%material(1)%elastic%get_sos2_mixture(rho,this%mumix,sosm)
        !    sos = sqrt(sosm)
        !    print *, 'mixt dens', rho(89,1,1), this%mumix(89,1,1)
        !endif
!print *, '----Exiting mix%getSOS----'
    end subroutine

    subroutine update_VF(this,isub,dt,rho,u,v,w,x,y,z,tsim,divu,src,x_bc,y_bc,z_bc)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w,divu,src
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: imat

        if (this%ns > 2) call GracefulExit("Figure out 2 materials first!",4356)

        !do imat = 1, this%ns
        !  call this%material(imat)%update_VF(this%material( 2-mod(imat+1,2) ),isub,dt,rho,u,v,w,x,y,z,tsim,divu,src,x_bc,y_bc,z_bc)
        !end do
        call this%material(1)%update_VF(this%material(2),isub,dt,rho,u,v,w,x,y,z,tsim,divu, src,x_bc,y_bc,z_bc)
        !call this%material(2)%update_VF(this%material(1),isub,dt,rho,u,v,w,x,y,z,tsim,divu,-src,x_bc,y_bc,z_bc)

        this%material(2)%VF = one - this%material(1)%VF

    end subroutine

    subroutine checkNaN(this)
        class(solid_mixture), intent(in)   :: this
        integer :: imat

        do imat=1,this%ns
            call this%material(imat)%checkNaN(imat)
        end do

    end subroutine

    subroutine fnumden(this,pf,fparams,iparams,num,den)
        class(solid_mixture), intent(inout)   :: this
        real(rkind), intent(in)               :: pf
        real(rkind), intent(in), dimension(:), target :: fparams
        integer, intent(in), dimension(:)     :: iparams
        real(rkind), intent(out)              :: num, den

        integer :: im, i, j, k
        real(rkind), dimension(:), pointer :: vf, gam, psph, pinf
        real(rkind) :: fac, gm1, YCv, pinfloc, rhoE

        if(iparams(1) == 1) then
            ! relaxPressure
            vf   => fparams(  1:this%ns)
            gam  => fparams(  this%ns+1:2*this%ns)
            psph => fparams(2*this%ns+1:3*this%ns)
            pinf => fparams(3*this%ns+1:4*this%ns)
            
            num = zero; den = zero;
            do im = 1, this%ns
              fac = vf(im)/gam(im)/(pinf(im)+pf)
              num = num + fac*(psph(im)-pf)
              !den = den + fac*(psph(im)-pinf(im)-two*pf)/(pinf(im)+pf)
              den = den - fac*(psph(im)+pinf(im))/(pinf(im)+pf)
              !write(*,*) im, psph(im), pf!-pinf(im)-two*pf, -(psph(im)+pinf(im))
            enddo
            nullify(vf,gam,psph,pinf)
        elseif(iparams(1)==2) then
            ! equilibratePressureTemperature
            i = iparams(2); j = iparams(3); k = iparams(4)
            num = zero; den = zero
            do im = 1, this%ns
              gm1 = this%material(im)%hydro%gam-one
              YCv = this%material(im)%Ys(i,j,k) * this%material(im)%hydro%Cv
              pinfloc = this%material(im)%hydro%PInf
              rhoE = fparams(1)

              num = num + YCv*(gm1*rhoE-(pf+this%material(im)%hydro%gam*pinfloc))/(pf+pinfloc)
              den = den + YCv*gm1*(pinfloc-rhoE) / (pf+pinfloc)**2
            enddo
        endif

    end subroutine

    subroutine rootfind_nr_1d(this,pf,fparams,iparams)
        use decomp_2d, only: nrank
        class(solid_mixture), intent(inout)   :: this
        real(rkind), intent(inout)            :: pf
        real(rkind), intent(in), dimension(:) :: fparams
        integer, intent(in), dimension(:)     :: iparams
    
        integer     :: ii, itmax = 1000
        !real(rkind) :: tol = 1.0d-8
        real(rkind) :: dpf, num, den, den_conv
    
        !pfinitguess = pf
        do ii = 1, itmax
          call this%fnumden(pf,fparams,iparams,num,den)
          ! if(dabs(den)>1.0d-12) then
          if(dabs(den)>1.0d-15) then
            dpf = num/den
          else
            write(*,*) 'den very small, please check.', num, num/den
            write(*,*) 'failure at proc ', nrank, ' at index ', iparams(2:4)
            stop
          endif
          pf = pf - dpf
          ! check for convergence
          if(dabs(pf)>1.0d-12) then
            den_conv = dabs(pf)
          else
            den_conv = one
          endif
          if(dabs(dpf)/den_conv<1.0d-8) exit
        enddo
        if(ii==itmax+1) then
          write(*,*) 'Newtons method for pf did not converge. Check details.', iparams(1)
        endif
    
    end subroutine rootfind_nr_1d


end module
