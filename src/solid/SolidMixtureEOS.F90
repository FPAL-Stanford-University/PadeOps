module SolidMixtureMod

    use kind_parameters, only: rkind,clen
    use constants,       only: zero,one,two,third
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
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

    contains

        procedure :: set_material
        procedure :: relaxPressure
        procedure :: update_g
        procedure :: update_Ys
        procedure :: update_eh
        procedure :: update_VF
        procedure :: filter
        procedure :: get_rho
        procedure :: get_primitive
        procedure :: get_conserved
        final     :: destroy

    end type

    interface solid_mixture
        module procedure init
    end interface

contains

    function init(decomp,der,ns) result(this)
        type(solid_mixture)               :: this
        type(decomp_info),  intent(in)    :: decomp
        type(derivatives),  intent(in)    :: der
        integer,            intent(in)    :: ns

        type(solid), allocatable :: dummy
        integer :: i

        this%ns = ns
        this%nxp = decomp%ysz(1)
        this%nyp = decomp%ysz(2)
        this%nzp = decomp%ysz(3)

        ! Allocate array of solid objects (Use a dummy to avoid memory leaks)
        allocate(dummy, source=solid(decomp,der))
        if (allocated(this%material)) deallocate(this%material)
        allocate(this%material(this%ns), source=dummy)
        deallocate(dummy)

    end function

    pure elemental subroutine destroy(this)
        type(solid_mixture), intent(inout)  :: this
        integer :: i

        ! Deallocate array of solids (Destructor of solid should take care of everything else)
        if (allocated(this%material)) deallocate(this%material)
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

    subroutine relaxPressure(this)
        class(solid_mixture), intent(inout) :: this

        real(rkind), dimension(4*this%ns), target :: fparams
        integer, dimension(1)             :: iparams
        real(rkind), dimension(:), pointer :: vf, gam, psph, pinf

        real(rkind), dimension(1:this%ns) :: fac

        ! get species hydrodynamic energy; mixture hydrodynamic energy
        ! get species pressure from species energy

        vf   => fparams(  1:this%ns)
        gam  => fparams(  this%ns+1:2*this%ns)
        psph => fparams(2*this%ns+1:3*this%ns)
        pinf => fparams(3*this%ns+1:4*this%ns)
        
        do k=1,this%nzp
         do j=1,this%nyp
          do i=1,this%nxp
            ! set fparams
            fparams(  1:this%ns)           = this%material(1:this%ns)%VF(i,j,k)    ! volume fractions
            fparams(  this%ns+1:2*this%ns) = this%material(1:this%ns)%hydro%gam    ! gamma
            fparams(2*this%ns+1:3*this%ns) = this%material(1:this%ns)%p(i,j,k)     ! pressure before eqb
            fparams(3*this%ns+1:4*this%ns) = this%material(1:this%ns)%hydro%PInf   ! PInf

            ! set iparams
            iparams(1) = 0     ! dummy; not really used

            ! scale all pressures by max over all PInfs
            maxp = maxval(fparams(3*this%ns+1:4*this%ns))
            fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)/maxp

            ! set initial guess
            peqb = sum(fparams(1:this%ns)*fparams(2*this%ns+1:3*this%ns))

            ! solve non-linear equation
            call this%rootfind_nr_1d(peqb,fparams,iparams)

            ! rescale all pressures by maxp
            fparams(2*this%ns+1:4*this%ns) = fparams(2*this%ns+1:4*this%ns)*maxp
            peqb = peqb*maxp

            ! update species VF, eh, ...
            fac = (psph + gam*pinf + (gam-one)*peqb)/(gam*(pinf+peqb))
            vf = vf*fac
            !this%material(1:this%ns)%g = this%material(1:this%ns)%g*fac**third        !! --- not clear if this is needed or if it works

            peqb = (rho(i,j,k)*eh(i,j,k) - sum(vf*gam*pinf/(gam-one))) / sum(vf/(gam-one))
            psph = peqb

          enddo
         enddo
        enddo

        ! compute get species energy from species pressure

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

    subroutine get_conserved(this,rho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%get_conserved(rho)
        end do

    end subroutine

    subroutine get_primitive(this,onebyrho)
        class(solid_mixture), intent(inout) :: this
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: onebyrho

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%get_primitive(onebyrho)
        end do

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

    subroutine filter(this, fil, 1)
        class(solid_mixture), intent(inout) :: this
        type(filters),        intent(in)    :: fil
        integer,              intent(in)    :: iflag

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%filter(fil, iflag)
        end do

    end subroutine

    subroutine update_g(this,isub,dt,rho,u,v,w)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%update_g(isub,dt,rho,u,v,w)
        end do

    end subroutine

    subroutine update_Ys(this,isub,dt,rho,u,v,w)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%update_Ys(isub,dt,rho,u,v,w)
        end do

    end subroutine

    subroutine update_eh(this,isub,dt,rho,u,v,w,tauiiart,divu)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w,tauiiart,divu

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%update_eh(isub,dt,rho,u,v,w,tauiiart,divu)
        end do

    end subroutine

    subroutine update_VF(this,isub,dt,rho,u,v,w)
        class(solid_mixture), intent(inout) :: this
        integer,              intent(in)    :: isub
        real(rkind),          intent(in)    :: dt
        real(rkind), dimension(this%nxp,this%nyp,this%nzp), intent(in) :: rho,u,v,w

        integer :: imat

        do imat = 1, this%ns
          call this%material(imat)%update_VF(isub,dt,rho,u,v,w)
        end do

    end subroutine

    subroutine fnumden(this,pf,fparams,iparams,num,den)
        class(solid_mixture), intent(inout)   :: this
        real(rkind), intent(in)               :: pf
        real(rkind), intent(in), dimension(:), target :: fparams
        integer, intent(in), dimension(:)     :: iparams
        real(rkind), intent(out)              :: num, den

        integer :: im
        real(rkind), dimension(:), pointer :: vf, gam, psph, pinf

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

    end subroutine

    subroutine rootfind_nr_1d(this,pf,fparams,iparams)
        class(solid_mixture), intent(inout)   :: this
        real(rkind), intent(inout)            :: pf
        real(rkind), intent(in), dimension(:) :: fparams
        integer, intent(in), dimension(:)     :: iparams
    
        integer     :: ii, itmax = 1000
        real(rkind) :: tol = 1.0d-8
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
