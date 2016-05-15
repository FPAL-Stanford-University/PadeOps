module SolidMixtureMod

    use kind_parameters, only: rkind,clen
    use constants,       only: zero,epssmall,eps,one,two,third
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
        type(ladobject),   pointer :: LAD

    contains

        procedure :: init
        procedure :: set_material
        procedure :: relaxPressure
        procedure :: equilibratePressureTemperature
        procedure :: getLAD
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
        procedure :: getSOS
        procedure :: get_J
        procedure :: get_q
        procedure :: get_qmix
        procedure :: get_dt
        procedure :: get_eelastic_devstress
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
    subroutine init(this,decomp,der,fil,LAD,ns)
        !type(solid_mixture)      , intent(inout) :: this
        class(solid_mixture)      , intent(inout) :: this
        type(decomp_info), target, intent(in)    :: decomp
        type(filters),     target, intent(in)    :: fil
        type(derivatives), target, intent(in)    :: der
        type(ladobject),   target, intent(in)    :: LAD
        integer,                   intent(in)    :: ns

        type(solid), allocatable :: dummy
        integer :: i

        if (ns < 1) call GracefulExit("Must have at least 1 species in the problem. Check input file for errors",3457)

        this%ns = ns
        this%nxp = decomp%ysz(1)
        this%nyp = decomp%ysz(2)
        this%nzp = decomp%ysz(3)

        this%decomp => decomp
        this%der => der
        this%fil => fil
        this%LAD => LAD

        ! Allocate array of solid objects (Use a dummy to avoid memory leaks)
        allocate(dummy)
        call dummy%init(decomp,der,fil)

        if (allocated(this%material)) deallocate(this%material)
        allocate(this%material(this%ns))!, source=dummy)
        do i=1,this%ns
            call this%material(i)%init(decomp,der,fil)
        end do
        deallocate(dummy)

        this%material(1)%Ys = one
        this%material(1)%VF = one
        do i=2,this%ns
            this%material(i)%Ys = zero
            this%material(i)%VF = zero
        end do

    end subroutine
    !end function

    pure elemental subroutine destroy(this)
        type(solid_mixture), intent(inout)  :: this

        ! Deallocate array of solids (Destructor of solid should take care of everything else)
        if (allocated(this%material)) deallocate(this%material)

        nullify(this%LAD)
        nullify(this%fil)
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

    subroutine equilibratePressureTemperature(this,mixRho,mixE,mixP,mixT)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: mixRho, mixE
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP, mixT

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix
        real(rkind), dimension(4*this%ns), target :: fparams
        integer, dimension(4)             :: iparams
        real(rkind), dimension(:), pointer :: vf, gam, psph, pinf

        integer :: i,j,k,imat
        real(rkind) :: maxp, peqb, pest, pdiffmax

        ! subtract elastic energy to determine hydrostatic energy. Temperature
        ! is assumed a function of only hydrostatic energy. Not sure if this is
        ! correct.
        ehmix = mixE
        do imat = 1, this%ns
            ehmix = ehmix - this%material(imat)%Ys * this%material(imat)%eel
        enddo

        !! identical materials
        !mixP = (this%material(1)%hydro%gam - one)*mixRho*ehmix - &
        !        this%material(1)%hydro%gam * this%material(1)%hydro%PInf

        !mixT = ehmix*(mixP + this%material(1)%hydro%PInf)/(this%material(1)%hydro%Cv *(mixP + this%material(1)%hydro%gam*this%material(1)%hydro%PInf))

        ! non-identical materials

        ! implement non-linear solve for mixP
        !vf   => fparams(  1:this%ns)
        !gam  => fparams(  this%ns+1:2*this%ns)
        !psph => fparams(2*this%ns+1:3*this%ns)
        !pinf => fparams(3*this%ns+1:4*this%ns)
        
        do k=1,this%nzp
         do j=1,this%nyp
          do i=1,this%nxp
            !do imat=1,this%ns

                !fparams(          imat) = this%material(imat)%VF(i,j,k)    ! volume fractions
                !fparams(  this%ns+imat) = this%material(imat)%hydro%gam    ! gamma
                !fparams(2*this%ns+imat) = this%material(imat)%p(i,j,k)     ! pressure before eqb
                !fparams(3*this%ns+imat) = this%material(imat)%hydro%PInf   ! PInf
            !end do

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


        mixT = zero
        do i = 1, this%ns
          mixT = mixT + this%material(i)%Ys*this%material(i)%hydro%Cv *(mixP + this%material(i)%hydro%gam*this%material(i)%hydro%PInf)/(mixP + this%material(i)%hydro%PInf)
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

    subroutine relaxPressure(this,rho,mixE,mixP)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho, mixE
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP

        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix
        real(rkind), dimension(4*this%ns), target :: fparams
        integer, dimension(1)             :: iparams
        real(rkind), dimension(:), pointer :: vf, gam, psph, pinf

        integer :: i,j,k,imat
        real(rkind) :: maxp, peqb, pest, pdiffmax

        real(rkind), dimension(1:this%ns) :: fac

        pdiffmax = zero
        ehmix = mixE
        do imat = 1, this%ns
            ehmix = ehmix - this%material(imat)%Ys * this%material(imat)%eel
            !write(*,'(a,i4,1x,2(e21.14,1x))') '-e-b-', imat, maxval(this%material(imat)%eh), minval(this%material(imat)%eh)
            !write(*,'(a,i4,1x,2(e21.14,1x))') '-p-b-', imat, maxval(this%material(imat)%p), minval(this%material(imat)%p)
        enddo

        ! equilibrate and reset species pressures, reset volume fractions

        vf   => fparams(  1:this%ns)
        gam  => fparams(  this%ns+1:2*this%ns)
        psph => fparams(2*this%ns+1:3*this%ns)
        pinf => fparams(3*this%ns+1:4*this%ns)
        
        do k=1,this%nzp
         do j=1,this%nyp
          do i=1,this%nxp
            do imat=1,this%ns
                ! set fparams
                fparams(          imat) = this%material(imat)%VF(i,j,k)    ! volume fractions
                fparams(  this%ns+imat) = this%material(imat)%hydro%gam    ! gamma
                fparams(2*this%ns+imat) = this%material(imat)%p(i,j,k)     ! pressure before eqb
                fparams(3*this%ns+imat) = this%material(imat)%hydro%PInf   ! PInf
            end do

            ! set iparams
            iparams(1) = 1     

            ! scale all pressures by max over all PInfs
            maxp = maxval(fparams(3*this%ns+1:4*this%ns))
            fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)/maxp

            ! set initial guess
            peqb = sum(fparams(1:this%ns)*fparams(2*this%ns+1:3*this%ns))
            pest = peqb

            !! solve non-linear equation
            !call this%rootfind_nr_1d(peqb,fparams,iparams)
            !pdiffmax = max(dabs(pest-peqb),pdiffmax)

            ! rescale all pressures by maxp
            fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)*maxp
            peqb = peqb*maxp

            ! update species VF, eh, ...
            fac = (psph + gam*pinf + (gam-one)*peqb)/(gam*(pinf+peqb))
            vf = vf*fac
            !this%material(1:this%ns)%g = this%material(1:this%ns)%g*fac**third        !! --- not clear if this is needed or if it works

            peqb = (rho(i,j,k)*ehmix(i,j,k) - sum(vf*gam*pinf/(gam-one))) / sum(vf/(gam-one))
            psph = peqb

            mixP(i,j,k) = peqb
            do imat=1,this%ns
              this%material(imat)%VF(i,j,k) = vf(imat)
              this%material(imat)%p(i,j,k) = psph(imat)
            enddo

          enddo
         enddo
        enddo
            !write(*,*) '--pdiff--', pdiffmax

        ! get species energy from species pressure
        do imat = 1, this%ns
            call this%material(imat)%get_ehydroT_from_p(rho)
        end do

    end subroutine

    subroutine get_dt(this, rho, delta, dtkap, dtDiff, dtplast)
        use reductions, only: P_MAXVAL
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho
        real(rkind), intent(in)  :: delta, dtkap
        real(rkind), intent(out) :: dtDiff, dtplast

        integer :: imat

        !dtkap = one/epssmall; 
        dtDiff = one/epssmall; dtplast = dtDiff

        do imat = 1, this%ns
          !dtkap =   min(dtkap,   one / ( (P_MAXVAL( this%material(imat)%kap*this%material(imat)%T/(rho* delta**4)))**(third) + eps))
          dtDiff =  min(dtDiff,  one / ( (P_MAXVAL( this%material(imat)%diff/delta**2) + eps)) )
          dtplast = min(dtplast, this%material(imat)%elastic%tau0)
        enddo

        ! For now disable plastic time step limit by setting a large value
        dtplast = real(1.0D32,rkind)

    end subroutine

    ! Subroutine to get species art. conductivities and diffusivities
    subroutine getLAD(this,rho,sos,x_bc,y_bc,z_bc)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho,sos  ! Mixture density and speed of sound
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: i

        do i = 1,this%ns
            call this%material(i)%getPhysicalProperties()

            !! Artificial conductivity
            !call this%LAD%get_conductivity(rho, this%material(i)%eh, this%material(i)%T, sos, &
            !                                    this%material(i)%kap, x_bc, y_bc, z_bc)
            this%material(i)%kap = zero
            ! Artificial diffusivity (grad(Ys) is stored in Ji at this stage)
            call this%LAD%get_diffusivity(this%material(i)%Ys, this%material(i)%Ji(:,:,:,1), &
                                          this%material(i)%Ji(:,:,:,2), this%material(i)%Ji(:,:,:,3), &
                                          sos, this%material(i)%diff, x_bc, y_bc, z_bc)
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

    subroutine get_ehydro_from_p(this,rho,p,T)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: p, T

        integer :: imat

        p = zero; T = zero

        ! get species energy from species pressure
        do imat = 1, this%ns
            call this%material(imat)%get_ehydroT_from_p(rho)             ! computes species ehydro and T from species p
            p = p + this%material(imat)%VF * this%material(imat)%p       ! computes mixture p - redundant when assuming equilibrium
            T = T + this%material(imat)%VF * this%material(imat)%T       ! computes mixture T - redundant when assuming equilibrium
        enddo

    end subroutine

    subroutine get_eelastic_devstress(this,devstress)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6), intent(out) :: devstress

        integer :: imat

        devstress = zero

        do imat = 1, this%ns
          call this%material(imat)%get_eelastic_devstress()
          devstress(:,:,:,1) = devstress(:,:,:,1) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,1)
          devstress(:,:,:,2) = devstress(:,:,:,2) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,2)
          devstress(:,:,:,3) = devstress(:,:,:,3) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,3)
          devstress(:,:,:,4) = devstress(:,:,:,4) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,4)
          devstress(:,:,:,5) = devstress(:,:,:,5) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,5)
          devstress(:,:,:,6) = devstress(:,:,:,6) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,6)
        end do

    end subroutine

    subroutine get_J(this,rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho  ! Mixture density

        integer :: i
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: sumJx,sumJy,sumJz

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

    pure subroutine get_conserved(this,rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%get_conserved(rho)
        end do

    end subroutine

    subroutine get_primitive(this,rho,devstress,p,sos,e)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho, e
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6), intent(out) :: devstress
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: p, sos

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%get_primitive(rho)
        end do
  
        call this%get_eelastic_devstress(devstress)   ! Get species elastic energies, and mixture and species devstress
        if(this%ns == 1) then
          this%material(1)%eh = e - this%material(1)%eel
        endif
        !call this%get_p_from_ehydro(rho)              ! Get species hydrodynamic energy, temperature; and mixture pressure, temperature
        !call this%get_pmix(p)
        !call this%getSOS(rho,p,sos)


    end subroutine

    subroutine get_rho(this,rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: rho

        integer :: imat

        rho = zero
        do imat = 1, this%ns
          rho = rho + this%material(imat)%consrv(:,:,:,1)
        end do

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

            !write(*,'(a,2(e21.14,1x))') '-sumJ1- ', maxval(this%material(1)%Ji(:,:,:,1) + this%material(2)%Ji(:,:,:,1)),  &
            !                                        minval(this%material(1)%Ji(:,:,:,1) + this%material(2)%Ji(:,:,:,1))
        do i = 1,this%ns
            ! Done in sgrid
            !! Step 1: Get qy
            !call this%der%ddy(this%material(i)%T,tmp1_in_y,y_bc(1),y_bc(2))
            !this%material(i)%qi(:,:,:,2) = -this%material(i)%kap*tmp1_in_y
            this%material(i)%qi(:,:,:,2) = zero

            !! Step 2: Get qx
            !call transpose_y_to_x(this%material(i)%T,tmp1_in_x,this%decomp)
            !call this%der%ddx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
            !call transpose_x_to_y(tmp2_in_x,tmp1_in_y,this%decomp)
            !this%material(i)%qi(:,:,:,1) = -this%material(i)%kap*tmp1_in_y
            this%material(i)%qi(:,:,:,1) = zero

            !! Step 3: Get qz
            !call transpose_y_to_z(this%material(i)%T,tmp1_in_z,this%decomp)
            !call this%der%ddz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
            !call transpose_z_to_y(tmp2_in_z,tmp1_in_y)
            !this%material(i)%qi(:,:,:,3) = -this%material(i)%kap*tmp1_in_y
            this%material(i)%qi(:,:,:,3) = zero

            ! If multispecies, add the inter-species enthalpy flux
            if (this%ns .GT. 1) then
                call this%material(i)%get_enthalpy(tmp1_in_y)
            !write(*,'(a,2(e21.14,1x))') '-enth-', maxval(tmp1_in_y),  minval(tmp1_in_y)
            !write(*,'(a,2(e21.14,1x))') '-Ji1- ', maxval(this%material(i)%Ji(:,:,:,1)),  minval(this%material(i)%Ji(:,:,:,1))
            !write(*,'(a,2(e21.14,1x))') '-Ji2- ', maxval(this%material(i)%Ji(:,:,:,2)),  minval(this%material(i)%Ji(:,:,:,2))
            !write(*,'(a,2(e21.14,1x))') '-Ji3- ', maxval(this%material(i)%Ji(:,:,:,3)),  minval(this%material(i)%Ji(:,:,:,3))
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
            !qx = qx + this%material(i)%VF * this%material(i)%qi(:,:,:,1)
            !qy = qy + this%material(i)%VF * this%material(i)%qi(:,:,:,2)
            !qz = qz + this%material(i)%VF * this%material(i)%qi(:,:,:,3)
            qx = qx + this%material(i)%qi(:,:,:,1)
            qy = qy + this%material(i)%qi(:,:,:,2)
            qz = qz + this%material(i)%qi(:,:,:,3)
        end do

    end subroutine

    subroutine filter(this, iflag)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: iflag

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%filter(iflag)
        end do

    end subroutine

    subroutine update_g(this,isub,dt,rho,u,v,w,x,y,z,tsim)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%update_g(isub,dt,rho,u,v,w,x,y,z,tsim)
        end do

    end subroutine

    subroutine update_Ys(this,isub,dt,rho,u,v,w,x,y,z,tsim)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%update_Ys(isub,dt,rho,u,v,w,x,y,z,tsim)
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

    subroutine update_eh(this,isub,dt,rho,u,v,w,x,y,z,tsim,divu,viscwork)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w,divu,viscwork

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%update_eh(isub,dt,rho,u,v,w,x,y,z,tsim,divu,viscwork)
          !write(*,'(a,i4,1x,2(e21.14,1x))') '--eh--', imat, maxval(this%material(imat)%eh), minval(this%material(imat)%eh)
        end do

    end subroutine

    subroutine getSOS(this,rho,p,sos)
        use constants, only: fourthird
        class(solid_mixture), intent(in) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho, p  ! Mixture density and pressure
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: sos     ! Mixture speed of sound

        integer :: i
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: mixmu, mixgam, mixgamPinf

        mixmu = zero; mixgam = zero; mixgamPinf = zero
        do i = 1,this%ns
            mixmu = mixmu + this%material(i)%VF * this%material(i)%elastic%mu
            mixgam = mixgam + this%material(i)%VF * this%material(i)%hydro%gam
            mixgamPinf = mixgamPinf + this%material(i)%VF * this%material(i)%hydro%gam * this%material(i)%hydro%Pinf
        end do

        sos = sqrt( (mixgam*p + mixgamPinf + fourthird*mixmu)/rho )

    end subroutine

    subroutine update_VF(this,isub,dt,u,v,w,x,y,z,tsim)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: u,v,w

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%update_VF(isub,dt,u,v,w,x,y,z,tsim)
        end do

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
            den = den + fac*(psph(im)-pinf(im)-two*pf)/(pinf(im)+pf)
          enddo

          nullify(vf,gam,psph,pinf)
        elseif(iparams(1) == 2) then
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
          if(dabs(den)>1.0d-12) then
            dpf = num/den
          else
            write(*,*) 'den very small, please check.', num, num/den
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
