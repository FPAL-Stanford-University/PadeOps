module SolidMixtureMod

    use kind_parameters, only: rkind,clen
    use constants,       only: zero,epssmall,eps,one,two,third
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use LADMod,          only: ladobject
    use exits,           only: GracefulExit
    use EOSMod,          only: eos
    ! use StiffGasEOS,     only: stiffgas
    ! use Sep1SolidEOS,    only: sep1solid
    use SolidMod,        only: solid
    use AbstractEOSMod,  only: abstracteos

    implicit none

    type :: solid_mixture

        integer :: ns
        integer :: nxp, nyp, nzp
        type(solid), dimension(:), allocatable :: material

        type(decomp_info), pointer :: decomp
        type(derivatives), pointer :: der
        type(filters),     pointer :: fil
        type(ladobject),   pointer :: LAD

        logical :: SOSmodel = .FALSE.           ! is sound speed given by `equilibrium' model? Alternative is `frozen' model. Check Saurel et al., JCP 2009.
        logical :: PTeqb = .TRUE.
        logical :: usegTg = .FALSE.

    contains

        procedure :: init
        procedure :: set_material
        ! procedure :: relaxPressure
        procedure :: equilibratePressureTemperature
        procedure :: getLAD
        procedure :: update_g
        procedure :: update_Ys
        !procedure :: update_eh
        !procedure :: update_VF
        procedure :: filter
        procedure :: get_rho
        procedure :: get_primitive
        procedure :: get_conserved
        ! procedure :: get_ehydro_from_p
        ! procedure :: get_p_from_ehydro
        procedure :: post_bc

        procedure :: get_rhoYs_from_gVF
        !procedure :: get_emix
        procedure :: get_pmix
        procedure :: get_Tmix
        ! procedure :: getSOS
        procedure :: get_J
        procedure :: get_q
        procedure :: get_qmix
        procedure :: get_dt
        ! procedure :: get_eelastic_devstress
        procedure :: checkNaN
        ! procedure :: fnumden
        ! procedure :: rootfind_nr_1d
        final     :: destroy

    end type

    !interface solid_mixture
    !    module procedure init
    !end interface

contains

    !function init(decomp,der,fil,LAD,ns) result(this)
    subroutine init(this,decomp,der,fil,LAD,ns,PTeqb,SOSmodel,usegTg)
        !type(solid_mixture)      , intent(inout) :: this
        class(solid_mixture)      , intent(inout) :: this
        type(decomp_info), target, intent(in)    :: decomp
        type(filters),     target, intent(in)    :: fil
        type(derivatives), target, intent(in)    :: der
        type(ladobject),   target, intent(in)    :: LAD
        integer,                   intent(in)    :: ns
        logical,                   intent(in)    :: PTeqb
        logical,                   intent(in)    :: SOSmodel
        logical,                   intent(in)    :: usegTg

        type(solid), allocatable :: dummy
        integer :: i

        if (ns < 1) call GracefulExit("Must have at least 1 species in the problem. Check input file for errors",3457)

        this%PTeqb    = PTeqb
        this%SOSmodel = SOSmodel
        this%usegTg  = usegTg

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
        call dummy%init(decomp,der,fil,this%PTeqb,this%usegTg)

        if (allocated(this%material)) deallocate(this%material)
        allocate(this%material(this%ns))!, source=dummy)
        do i=1,this%ns
            call this%material(i)%init(decomp,der,fil,this%PTeqb,this%usegTg)
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

    subroutine set_material(this, imat, eos)
        use Sep1SolidEOSMod, only: sep1solideos
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: imat
        class(abstracteos),   intent(in)    :: eos

        if ((imat .GT. this%ns) .OR. (imat .LE. 0)) call GracefulExit("Cannot set material with index greater than the number of species.",4534)
       
        if (allocated(this%material(imat)%eos)) deallocate(this%material(imat)%eos)
        allocate( this%material(imat)%eos, source=eos )
    end subroutine

    !ADD! subroutine relaxPressure(this,rho,mixE,mixP)
    !ADD!     class(solid_mixture), intent(inout) :: this
    !ADD!     real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho, mixE
    !ADD!     real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP

    !ADD!     real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix
    !ADD!     real(rkind), dimension(4*this%ns), target :: fparams
    !ADD!     integer, dimension(1)             :: iparams
    !ADD!     real(rkind), dimension(:), pointer :: vf, gam, psph, pinf

    !ADD!     integer :: i,j,k,imat
    !ADD!     real(rkind) :: maxp, peqb

    !ADD!     real(rkind), dimension(1:this%ns) :: fac

    !ADD!     ehmix = mixE
    !ADD!     do imat = 1, this%ns
    !ADD!         ehmix = ehmix - this%material(imat)%Ys * this%material(imat)%eel
    !ADD!     enddo

    !ADD!     ! equilibrate and reset species pressures, reset volume fractions

    !ADD!     vf   => fparams(  1:this%ns)
    !ADD!     gam  => fparams(  this%ns+1:2*this%ns)
    !ADD!     psph => fparams(2*this%ns+1:3*this%ns)
    !ADD!     pinf => fparams(3*this%ns+1:4*this%ns)
    !ADD!     
    !ADD!     do k=1,this%nzp
    !ADD!      do j=1,this%nyp
    !ADD!       do i=1,this%nxp
    !ADD!         do imat=1,this%ns
    !ADD!             ! set fparams
    !ADD!             fparams(          imat) = this%material(imat)%VF(i,j,k)    ! volume fractions
    !ADD!             fparams(  this%ns+imat) = this%material(imat)%hydro%gam    ! gamma
    !ADD!             fparams(2*this%ns+imat) = this%material(imat)%p(i,j,k)     ! pressure before eqb
    !ADD!             fparams(3*this%ns+imat) = this%material(imat)%hydro%PInf   ! PInf
    !ADD!         end do

    !ADD!         ! set iparams
    !ADD!         iparams(1) = 1     !   used in fnumden; 2 for PTeqb 

    !ADD!         ! scale all pressures by max over all PInfs
    !ADD!         maxp = maxval(fparams(3*this%ns+1:4*this%ns))
    !ADD!         fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)/maxp

    !ADD!         ! set initial guess
    !ADD!         peqb = sum(fparams(1:this%ns)*fparams(2*this%ns+1:3*this%ns))

    !ADD!         ! solve non-linear equation
    !ADD!         call this%rootfind_nr_1d(peqb,fparams,iparams)

    !ADD!         ! rescale all pressures by maxp
    !ADD!         fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)*maxp
    !ADD!         peqb = peqb*maxp

    !ADD!         ! update species VF, eh, ...
    !ADD!         fac = (psph + gam*pinf + (gam-one)*peqb)/(gam*(pinf+peqb))
    !ADD!         vf = vf*fac
    !ADD!         !this%material(1:this%ns)%g = this%material(1:this%ns)%g*fac**third        !! --- not clear if this is needed or if it works

    !ADD!         peqb = (rho(i,j,k)*ehmix(i,j,k) - sum(vf*gam*pinf/(gam-one))) / sum(vf/(gam-one))
    !ADD!         psph = peqb

    !ADD!         mixP(i,j,k) = peqb
    !ADD!         do imat=1,this%ns
    !ADD!           this%material(imat)%VF(i,j,k) = vf(imat)
    !ADD!           this%material(imat)%p(i,j,k) = psph(imat)
    !ADD!         enddo

    !ADD!       enddo
    !ADD!      enddo
    !ADD!     enddo

    !ADD!     ! get species energy from species pressure
    !ADD!     do imat = 1, this%ns
    !ADD!         call this%material(imat)%get_ehydroT_from_p(rho)
    !ADD!     end do

    !ADD! end subroutine

    !ADD! subroutine equilibratePressureTemperature(this,mixRho,mixE,mixP,mixT)
    !ADD!     class(solid_mixture), intent(inout) :: this
    !ADD!     real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: mixRho, mixE
    !ADD!     real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: mixP, mixT

    !ADD!     real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: ehmix
    !ADD!     real(rkind), dimension(4*this%ns), target :: fparams
    !ADD!     integer, dimension(4)             :: iparams
    !ADD!     ! real(rkind), dimension(:), pointer :: vf, psph, pinf

    !ADD!     integer :: i,j,k,imat
    !ADD!     real(rkind) :: maxp, peqb !, pest, pdiffmax

    !ADD!     ! subtract elastic energy to determine hydrostatic energy. Temperature
    !ADD!     ! is assumed a function of only hydrostatic energy. Not sure if this is
    !ADD!     ! correct.
    !ADD!     ehmix = mixE
    !ADD!     do imat = 1, this%ns
    !ADD!         ehmix = ehmix - this%material(imat)%Ys * this%material(imat)%eel
    !ADD!     enddo

    !ADD!     do k=1,this%nzp
    !ADD!      do j=1,this%nyp
    !ADD!       do i=1,this%nxp
    !ADD!         ! set fparams
    !ADD!         fparams(1) = mixRho(i,j,k)*ehmix(i,j,k)

    !ADD!         ! set iparams
    !ADD!         iparams(1) = 2
    !ADD!         iparams(2) = i; iparams(3) = j; iparams(4) = k;

    !ADD!         maxp = zero; peqb = zero
    !ADD!         do imat=1,this%ns
    !ADD!           !! determine max over all PInfs
    !ADD!           !maxp = maxval(maxp, this%material(imat)%hydro%PInf)

    !ADD!           ! set initial guess
    !ADD!           peqb = peqb + this%material(imat)%VF(i,j,k)*this%material(imat)%p(i,j,k)
    !ADD!         end do
    !ADD!         !pest = peqb

    !ADD!         ! solve non-linear equation
    !ADD!         call this%rootfind_nr_1d(peqb,fparams,iparams)
    !ADD!         !pdiffmax = max(dabs(pest-peqb),pdiffmax)

    !ADD!         !! rescale all pressures by maxp
    !ADD!         !fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)*maxp
    !ADD!         mixP(i,j,k) = peqb !*maxp

    !ADD!       enddo
    !ADD!      enddo
    !ADD!     enddo

    !ADD!     mixT = zero
    !ADD!     do i = 1, this%ns
    !ADD!       mixT = mixT + this%material(i)%Ys*this%material(i)%hydro%Cv * &
    !ADD!             (mixP + this%material(i)%hydro%gam*this%material(i)%hydro%PInf)/(mixP + this%material(i)%hydro%PInf)
    !ADD!     enddo
    !ADD!     mixT = ehmix/mixT

    !ADD!     do i = 1, this%ns
    !ADD!         this%material(i)%T = mixT
    !ADD!         this%material(i)%p = mixP
    !ADD!         this%material(i)%VF = mixRho*this%material(i)%Ys*(this%material(i)%hydro%gam-one)* &
    !ADD!                                      this%material(i)%hydro%Cv*mixT/(mixP + this%material(i)%hydro%PInf)

    !ADD!         this%material(i)%eh = this%material(i)%hydro%Cv*mixT*(mixP + this%material(i)%hydro%gam*this%material(i)%hydro%PInf) / &
    !ADD!                                                              (mixP + this%material(i)%hydro%PInf)
    !ADD!     end do

    !ADD! end subroutine

    subroutine equilibratePressureTemperature(this,mixRho,mixE)
        use eqbPTFunctionMod, only: eqbPTFunction
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: mixRho, mixE

        real(rkind), dimension(this%ns) :: VF0, Y
        real(rkind), dimension(2*this%ns) :: xvar, fvar
        real(rkind), dimension(9,this%ns) :: g0

        !real(rkind), dimension(4*this%ns), target :: fparams
        integer :: i, j, k, m
        type(eqbPTFunction) :: eqbfn

        do k=1,this%nzp
         do j=1,this%nyp
          do i=1,this%nxp
              ! create object that extends abstract class newton_functor
              do m = 1, this%ns
                  VF0(m)  = this%material(m)%VF(i,j,k)
                  Y(m)    = this%material(m)%Ys(i,j,k)
                  g0(:,m) = this%material(m)%g(i,j,k,:)
              end do
              write(*,*) '---', i
              if(minval(VF0) < 1.0D-6) then
                  !if(VF0(1) > VF0(2)) then
                  !   this%material(2)%p = this%material(1)%p
                  !   this%material(2)%T = this%material(1)%T
                  !else
                  !   this%material(1)%p = this%material(2)%p
                  !   this%material(1)%T = this%material(2)%T
                  !endif
                  cycle
              endif
              eqbfn = eqbPTFunction(this%ns, this%material, VF0, g0, Y, mixE(i,j,k))

              ! set initial guesses and call newton_solve
              do m = 1, this%ns
                  xvar(m) = this%material(m)%energy(i,j,k)
              end do
              xvar(this%ns+1:2*this%ns) = VF0(:)
              fvar = 0.0D0
                  !if(i==65 .and. j==1 .and. k==1) then
                  !    write(*,*) 'xvar before:    ', xvar
                  !    write(*,*) 'fvar before:    ', fvar
                  !endif
              call eqbfn%newton_solve(xvar, fvar)
                  !if(i==65 .and. j==1 .and. k==1) then
                      !write(*,'(a,4(e19.12,1x))') 'xvar after:    ', xvar
                  !endif

              ! store returned VFs, gs and energies
              do m = 1, this%ns
                  this%material(m)%energy(i,j,k) = xvar(m)
                  this%material(m)%VF(i,j,k)  = xvar(this%ns+m)
                  this%material(m)%g(i,j,k,:) = g0(:,m)*xvar(this%ns+m)**third
                  this%material(m)%p(i,j,k) = eqbfn%p(m)
                  this%material(m)%T(i,j,k) = eqbfn%T(m)
                  !if(i==65 .and. j==1 .and. k==1) then
                  !    write(*,*) 'pressures:    ', eqbfn%p(m)
                  !    write(*,*) 'temperatures: ', eqbfn%T(m)
                  !endif
              end do
          end do
         end do
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
          !ADD! dtplast = min(dtplast, this%material(imat)%elastic%tau0)
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

            if (.NOT. this%PTeqb) then
                ! Artificial conductivity
                call this%LAD%get_conductivity(rho, this%material(i)%energy, this%material(i)%T, sos, &
                                                    this%material(i)%kap, x_bc, y_bc, z_bc)
            end if

            ! Artificial diffusivity (grad(Ys) is stored in Ji at this stage)
            call this%LAD%get_diffusivity(this%material(i)%Ys, this%material(i)%Ji(:,:,:,1), &
                                          this%material(i)%Ji(:,:,:,2), this%material(i)%Ji(:,:,:,3), &
                                          sos, this%material(i)%diff, x_bc, y_bc, z_bc)
        end do

    end subroutine

    ! subroutine get_p_from_ehydro(this, rho)
    !     class(solid_mixture), intent(inout) :: this
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho

    !     integer :: imat

    !     ! get species pressure from species energy
    !     do imat = 1, this%ns
    !         call this%material(imat)%get_p_from_ehydro(rho)
    !     enddo

    ! end subroutine

    ! subroutine get_ehydro_from_p(this,rho)
    !     class(solid_mixture), intent(inout) :: this
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho

    !     integer :: imat

    !     ! get species energy from species pressure
    !     do imat = 1, this%ns
    !         call this%material(imat)%get_ehydroT_from_p(rho)             ! computes species ehydro and T from species p
    !     enddo

    ! end subroutine

    subroutine get_Tmix(this,T)
        class(solid_mixture), intent(in) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: T  ! Mixture temperature

        integer :: i

        T = zero
        do i = 1,this%ns
            T = T + this%material(i)%VF * this%material(i)%T  ! Volume fraction weighted sum
        end do
    end subroutine

    ! subroutine get_eelastic_devstress(this,devstress)
    !     class(solid_mixture), intent(inout) :: this
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp,6), intent(out) :: devstress

    !     integer :: imat

    !     devstress = zero

    !     do imat = 1, this%ns
    !       call this%material(imat)%get_eelastic_devstress()
    !       devstress(:,:,:,1) = devstress(:,:,:,1) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,1)
    !       devstress(:,:,:,2) = devstress(:,:,:,2) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,2)
    !       devstress(:,:,:,3) = devstress(:,:,:,3) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,3)
    !       devstress(:,:,:,4) = devstress(:,:,:,4) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,4)
    !       devstress(:,:,:,5) = devstress(:,:,:,5) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,5)
    !       devstress(:,:,:,6) = devstress(:,:,:,6) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,6)
    !     end do

    ! end subroutine

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
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in) :: rho, e
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6), intent(out) :: devstress
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(out) :: p, sos
        
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhom, sosm

        integer :: imat

        p = zero
        sos = zero

        if(this%pTeqb) then
           call this%equilibratePressureTemperature(rho, e)
        endif
         
        do imat = 1, this%ns
          call this%material(imat)%get_primitive(rho,sosm) !ADD! This is only for single species. Need to fix for multi-species
          
          ! Add contribution to mixture pressure
          p = p + this%material(imat)%VF * this%material(imat)%p
         
          ! Add contribution to mixture devstress
          devstress(:,:,:,1) = devstress(:,:,:,1) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,1)
          devstress(:,:,:,2) = devstress(:,:,:,2) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,2)
          devstress(:,:,:,3) = devstress(:,:,:,3) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,3)
          devstress(:,:,:,4) = devstress(:,:,:,4) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,4)
          devstress(:,:,:,5) = devstress(:,:,:,5) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,5)
          devstress(:,:,:,6) = devstress(:,:,:,6) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,6)

          ! Add contribution to mixture speed of sound
          sos = sos + sosm
          if(this%SOSmodel) then
              ! equilibrium model
              call this%material(imat)%getSpeciesDensity(rho,rhom)
              sos = sos + this%material(imat)%VF/(rhom*sosm)
          else
              ! frozen model (details in Saurel et al. 2009)
              sos = sos + this%material(imat)%Ys*sosm
          endif

        end do

        if(this%SOSmodel) then
            sos = one / (sqrt(rho*sos) + epssmall)
        else
            sos = sqrt(sos)
        endif

        !ADD! call this%get_eelastic_devstress(devstress)   ! Get species elastic energies, and mixture and species devstress
        !ADD! if(this%ns == 1) then
        !ADD!   this%material(1)%eh = e - this%material(1)%eel ! Since eh equation is not exact and this is a better alternative for 1 species
        !ADD! endif
        !ADD! if(this%PTeqb) then
        !ADD!     call this%get_p_from_ehydro(rho)   ! Get species pressures from species hydrodynamic energy 
        !ADD!     call this%get_pmix(p)              ! Get mixture pressure from species pressures
        !ADD! endif
        
        ! call this%getSOS(rho,p,sos)

    end subroutine

    subroutine post_bc(this, rho, e, p, devstress, T, sos)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(in)  :: rho
        real(rkind), dimension(this%nxp,this%nyp,this%nzp,6), intent(out) :: devstress
        real(rkind), dimension(this%nxp,this%nyp,this%nzp),   intent(out) :: e, p, T, sos
        
        real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhom, sosm

        integer :: imat

        e = zero
        p = zero
        T = zero
        sos = zero
        devstress = zero
         
        do imat = 1, this%ns
            call this%material(imat)%get_energy(rho) ! Update internal energy
            call this%material(imat)%get_primitive(rho,sosm) !ADD! Material density is being calculated multiple times. Add as a field to SolidMod to avoid this
        
            ! Add contribution to mixture internal energy
            e = e + this%material(imat)%Ys * this%material(imat)%energy
            
            ! Add contribution to mixture pressure
            p = p + this%material(imat)%VF * this%material(imat)%p

            ! Add contribution to mixture temperature
            T = T + this%material(imat)%Ys * this%material(imat)%T

            ! Add contribution to mixture devstress
            devstress(:,:,:,1) = devstress(:,:,:,1) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,1)
            devstress(:,:,:,2) = devstress(:,:,:,2) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,2)
            devstress(:,:,:,3) = devstress(:,:,:,3) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,3)
            devstress(:,:,:,4) = devstress(:,:,:,4) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,4)
            devstress(:,:,:,5) = devstress(:,:,:,5) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,5)
            devstress(:,:,:,6) = devstress(:,:,:,6) + this%material(imat)%VF * this%material(imat)%devstress(:,:,:,6)

            ! Add contribution to mixture speed of sound
            sos = sos + sosm
            if(this%SOSmodel) then
                ! equilibrium model
                call this%material(imat)%getSpeciesDensity(rho,rhom)
                sos = sos + this%material(imat)%VF/(rhom*sosm)
            else
                ! frozen model (details in Saurel et al. 2009)
                sos = sos + this%material(imat)%Ys*sosm
            endif

        end do

        if(this%SOSmodel) then
            sos = one / (sqrt(rho*sos) + epssmall)
        else
            sos = sqrt(sos)
        endif

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
            call this%material(imat)%get_conserved(rho)
        end do
    end subroutine

    subroutine get_q(this,rho,x_bc,y_bc,z_bc)
        use decomp_2d, only: transpose_y_to_x, transpose_x_to_y, transpose_y_to_z, transpose_z_to_y
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: i
        real(rkind), dimension(this%decomp%xsz(1),this%decomp%xsz(2),this%decomp%xsz(3)) :: tmp1_in_x, tmp2_in_x
        real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3)) :: tmp1_in_y
        real(rkind), dimension(this%decomp%zsz(1),this%decomp%zsz(2),this%decomp%zsz(3)) :: tmp1_in_z, tmp2_in_z

        do i = 1,this%ns
            if(.NOT. this%PTeqb) then
                ! Step 1: Get qy
                call this%der%ddy(this%material(i)%T,tmp1_in_y,y_bc(1),y_bc(2))
                this%material(i)%qi(:,:,:,2) = -this%material(i)%VF*this%material(i)%kap*tmp1_in_y

                ! Step 2: Get qx
                call transpose_y_to_x(this%material(i)%T,tmp1_in_x,this%decomp)
                call this%der%ddx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
                call transpose_x_to_y(tmp2_in_x,tmp1_in_y,this%decomp)
                this%material(i)%qi(:,:,:,1) = -this%material(i)%VF*this%material(i)%kap*tmp1_in_y

                ! Step 3: Get qz
                call transpose_y_to_z(this%material(i)%T,tmp1_in_z,this%decomp)
                call this%der%ddz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
                call transpose_z_to_y(tmp2_in_z,tmp1_in_y)
                this%material(i)%qi(:,:,:,3) = -this%material(i)%VF*this%material(i)%kap*tmp1_in_y
            else
                this%material(i)%qi(:,:,:,:) = zero     ! calculated in sgrid
            endif

            ! If multispecies, add the inter-species enthalpy flux
            if (this%ns .GT. 1) then
                call this%material(i)%get_enthalpy(rho,tmp1_in_y)
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

    end subroutine

    subroutine update_g(this,isub,dt,rho,u,v,w,x,y,z,tsim,x_bc,y_bc,z_bc)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt,tsim
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w
        integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

        integer :: imat

        do imat = 1, this%ns
            if (this%usegTg) then
                call this%material(imat)%update_gTg(isub,dt,rho,u,v,w,x,y,z,tsim,x_bc,y_bc,z_bc)
            else
                call this%material(imat)%update_g(isub,dt,rho,u,v,w,x,y,z,tsim,x_bc,y_bc,z_bc)
            end if
        end do

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

    !subroutine get_emix(this,e)
    !    class(solid_mixture), intent(in) :: this
    !    real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: e  ! Mixture internal energy

    !    integer :: i

    !    e = zero
    !    do i = 1,this%ns
    !        e = e + this%material(i)%Ys * ( this%material(i)%eh + this%material(i)%eel )  ! Mass fraction weighted sum
    !    end do

    !end subroutine

    !subroutine update_eh(this,isub,dt,rho,u,v,w,x,y,z,tsim,divu,viscwork,x_bc,y_bc,z_bc)
    !    class(solid_mixture), intent(inout) :: this
    !    integer,              intent(in)    :: isub
    !    real(rkind),          intent(in)    :: dt,tsim
    !    real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
    !    real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w,divu,viscwork
    !    integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

    !    integer :: imat

    !    do imat = 1, this%ns
    !      call this%material(imat)%update_eh(isub,dt,rho,u,v,w,x,y,z,tsim,divu,viscwork,x_bc,y_bc,z_bc)
    !    end do

    !end subroutine

    ! subroutine getSOS(this,rho,p,sos)
    !     class(solid_mixture), intent(in) :: this
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in)  :: rho, p  ! Mixture density and pressure
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(out) :: sos     ! Mixture speed of sound

    !     integer :: i
    !     real(rkind), dimension(this%nxp,this%nyp,this%nzp) :: rhom, sosm

    !     sos = zero
    !     do i = 1,this%ns
    !         call this%material(i)%getSpeciesDensity(rho,rhom)
    !         call this%material(i)%hydro%get_sos2(rhom,p,sosm)
    !         call this%material(i)%elastic%get_sos2(rhom,sosm)
    !         if(this%SOSmodel) then
    !             ! equilibrium model
    !             sos = sos + this%material(i)%VF/(rhom*sosm)
    !         else
    !             ! frozen model (details in Saurel et al. 2009)
    !             sos = sos + this%material(i)%Ys*sosm
    !         endif
    !     end do
    !     if(this%SOSmodel) then
    !         sos = one / (sqrt(rho*sos) + epssmall)
    !     else
    !         sos = sqrt(sos)
    !     endif

    ! end subroutine

    !subroutine update_VF(this,isub,dt,u,v,w,x,y,z,tsim,x_bc,y_bc,z_bc)
    !    class(solid_mixture), intent(inout) :: this
    !    integer,              intent(in)    :: isub
    !    real(rkind),          intent(in)    :: dt,tsim
    !    real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: x,y,z
    !    real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: u,v,w
    !    integer, dimension(2), intent(in) :: x_bc, y_bc, z_bc

    !    integer :: imat

    !    do imat = 1, this%ns
    !      call this%material(imat)%update_VF(isub,dt,u,v,w,x,y,z,tsim,x_bc,y_bc,z_bc)
    !    end do

    !end subroutine

    subroutine checkNaN(this)
        class(solid_mixture), intent(in)   :: this
        integer :: imat

        do imat=1,this%ns
            call this%material(imat)%checkNaN(imat)
        end do

    end subroutine

    !ADD! subroutine fnumden(this,pf,fparams,iparams,num,den)
    !ADD!     class(solid_mixture), intent(inout)   :: this
    !ADD!     real(rkind), intent(in)               :: pf
    !ADD!     real(rkind), intent(in), dimension(:), target :: fparams
    !ADD!     integer, intent(in), dimension(:)     :: iparams
    !ADD!     real(rkind), intent(out)              :: num, den

    !ADD!     integer :: im, i, j, k
    !ADD!     real(rkind), dimension(:), pointer :: vf, gam, psph, pinf
    !ADD!     real(rkind) :: fac, gm1, YCv, pinfloc, rhoE

    !ADD!     if(iparams(1) == 1) then
    !ADD!         ! relaxPressure
    !ADD!         vf   => fparams(  1:this%ns)
    !ADD!         gam  => fparams(  this%ns+1:2*this%ns)
    !ADD!         psph => fparams(2*this%ns+1:3*this%ns)
    !ADD!         pinf => fparams(3*this%ns+1:4*this%ns)
    !ADD!         
    !ADD!         num = zero; den = zero;
    !ADD!         do im = 1, this%ns
    !ADD!           fac = vf(im)/gam(im)/(pinf(im)+pf)
    !ADD!           num = num + fac*(psph(im)-pf)
    !ADD!           den = den + fac*(psph(im)-pinf(im)-two*pf)/(pinf(im)+pf)
    !ADD!         enddo
    !ADD!         nullify(vf,gam,psph,pinf)
    !ADD!     elseif(iparams(1)==2) then
    !ADD!         ! equilibratePressureTemperature
    !ADD!         i = iparams(2); j = iparams(3); k = iparams(4)
    !ADD!         num = zero; den = zero
    !ADD!         do im = 1, this%ns
    !ADD!           gm1 = this%material(im)%hydro%gam-one
    !ADD!           YCv = this%material(im)%Ys(i,j,k) * this%material(im)%hydro%Cv
    !ADD!           pinfloc = this%material(im)%hydro%PInf
    !ADD!           rhoE = fparams(1)

    !ADD!           num = num + YCv*(gm1*rhoE-(pf+this%material(im)%hydro%gam*pinfloc))/(pf+pinfloc)
    !ADD!           den = den + YCv*gm1*(pinfloc-rhoE) / (pf+pinfloc)**2
    !ADD!         enddo
    !ADD!     endif

    !ADD! end subroutine

    !ADD! subroutine rootfind_nr_1d(this,pf,fparams,iparams)
    !ADD!     use decomp_2d, only: nrank
    !ADD!     class(solid_mixture), intent(inout)   :: this
    !ADD!     real(rkind), intent(inout)            :: pf
    !ADD!     real(rkind), intent(in), dimension(:) :: fparams
    !ADD!     integer, intent(in), dimension(:)     :: iparams
    !ADD! 
    !ADD!     integer     :: ii, itmax = 1000
    !ADD!     !real(rkind) :: tol = 1.0d-8
    !ADD!     real(rkind) :: dpf, num, den, den_conv
    !ADD! 
    !ADD!     !pfinitguess = pf
    !ADD!     do ii = 1, itmax
    !ADD!       call this%fnumden(pf,fparams,iparams,num,den)
    !ADD!       ! if(dabs(den)>1.0d-12) then
    !ADD!       if(dabs(den)>1.0d-15) then
    !ADD!         dpf = num/den
    !ADD!       else
    !ADD!         write(*,*) 'den very small, please check.', num, num/den
    !ADD!         write(*,*) 'failure at proc ', nrank, ' at index ', iparams(2:4)
    !ADD!         stop
    !ADD!       endif
    !ADD!       pf = pf - dpf
    !ADD!       ! check for convergence
    !ADD!       if(dabs(pf)>1.0d-12) then
    !ADD!         den_conv = dabs(pf)
    !ADD!       else
    !ADD!         den_conv = one
    !ADD!       endif
    !ADD!       if(dabs(dpf)/den_conv<1.0d-8) exit
    !ADD!     enddo
    !ADD!     if(ii==itmax+1) then
    !ADD!       write(*,*) 'Newtons method for pf did not converge. Check details.', iparams(1)
    !ADD!     endif
    !ADD! 
    !ADD! end subroutine rootfind_nr_1d


end module
