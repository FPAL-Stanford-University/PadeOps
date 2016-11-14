module eqbPTFunctionMod
    use kind_parameters, only: rkind
    use NewtonSolverMod, only: newton_functor
    use SolidMod,        only: solid
    implicit none

    type, extends(newton_functor) :: eqbPTFunction    ! only 2 materials for now
        integer                                  :: ns
        type(solid), dimension(:),   pointer     :: material
        real(rkind), dimension(:),   allocatable :: VF0, Y, p, T
        real(rkind), dimension(:,:), allocatable :: g0
        real(rkind)                              :: emix
    contains
        procedure :: set_parameters
        procedure :: evaluate
        procedure :: gradient
        final     :: destroy
    end type

    interface eqbpTFunction
        module procedure init
    end interface
contains

    function init(ns,material) result(this)
        use exits, only: GracefulExit
        type(eqbpTFunction) :: this
        integer,                            intent(in) :: ns
        type(solid), dimension(ns), target, intent(in) :: material

        if (ns > 2) call GracefulExit("Only 2 materials are supported with PT equilibrium formulation for now.",3454)

        this%ns = ns
        this%material => material

        if (allocated(this%VF0)) deallocate(this%VF0); allocate( this%VF0(this%ns)  )
        if (allocated(this%g0) ) deallocate(this%g0);  allocate( this%g0(9,this%ns) )
        if (allocated(this%Y)  ) deallocate(this%Y);   allocate( this%Y(this%ns)    )
        if (allocated(this%p)  ) deallocate(this%p);   allocate( this%p(this%ns)    )
        if (allocated(this%T)  ) deallocate(this%T);   allocate( this%T(this%ns)    )

        this%niters = 10
        this%tolerance = real(1.D-18,rkind)
    end function
    
    pure elemental subroutine destroy(this)
        type(eqbpTFunction), intent(inout) :: this

        nullify( this%material )
        if (allocated(this%VF0)) deallocate(this%VF0);
        if (allocated(this%g0) ) deallocate(this%g0);
        if (allocated(this%Y)  ) deallocate(this%Y);
        if (allocated(this%p)  ) deallocate(this%p);
        if (allocated(this%T)  ) deallocate(this%T);
    end subroutine

    subroutine set_parameters(this,VF0,g0,Y,emix)
        class(eqbpTFunction),              intent(inout) :: this
        real(rkind), dimension(this%ns),   intent(in)    :: VF0, Y
        real(rkind), dimension(9,this%ns), intent(in)    :: g0
        real(rkind),                       intent(in)    :: emix
        
        this%VF0(:) = VF0(:)
        this%g0(:,:) = g0(:,:)
        this%Y(:)  = Y(:)
        this%emix = emix
    end subroutine

    subroutine evaluate(this, x, y)
        use constants, only: one
        class(eqbpTFunction),              intent(inout) :: this
        real(rkind), dimension(:),         intent(in)    :: x
        real(rkind), dimension(size(x,1)), intent(out)   :: y

        real(rkind), dimension(this%ns) :: energy, VF!, p, T
        integer :: m

        !write(*,*) '------In evaluate------------------------'
        !write(*,*) '  xvar: ', x
        do m = 1, this%ns
            !write(*,*) '+++++ material = ', m
            energy(m) = x(m); VF(m) = x(this%ns+m)
            call this%material(m)%eos%get_pT_from_energyVF(this%VF0(m), this%g0(:,m), energy(m), VF(m), this%p(m), this%T(m))
        enddo
        !write(*,*) '  this%p: ', this%p
        !write(*,*) '  this%T: ', this%T
        ! set elements of jacobian - obviously works for 2 materials only
        y(1) = this%p(1) - this%p(2)
        y(2) = this%T(1) - this%T(2)
        y(3) = VF(1) + VF(2) - one
        y(4) = this%Y(1)*energy(1) + this%Y(2)*energy(2) - this%emix
        !write(*,*) '  fvar: ', y
        !write(*,*) '------Done evaluate------------------------'

    end subroutine

    subroutine gradient(this, x, grad)
        use constants, only: zero, one
        class(eqbpTFunction),                        intent(in)  :: this
        real(rkind), dimension(:),                   intent(in)  :: x
        real(rkind), dimension(size(x,1),size(x,1)), intent(out) :: grad

        real(rkind), dimension(this%ns) :: energy, VF, dpde, dpdVF, dTde, dTdVF
        integer :: m

        !grad = zero

        do m = 1, this%ns
            energy(m) = x(m); VF(m) = x(this%ns+m)
            call this%material(m)%eos%get_pT_derivatives_wrt_energyVF(this%VF0(m), this%g0(:,m), energy(m), VF(m), dpde(m), dpdVF(m), dTde(m), dTdVF(m))
        enddo
        ! set elements of jacobian - obviously works for 2 materials only
        grad(1,1) = dpde(1);   grad(1,2) = -dpde(2);   grad(1,3) = dpdVF(1); grad(1,4) = -dpdVF(2)
        grad(2,1) = dTde(1);   grad(2,2) = -dTde(2);   grad(2,3) = dTdVF(1); grad(2,4) = -dTdVF(2)
        grad(3,1) = zero   ;   grad(3,2) =  zero   ;   grad(3,3) = one     ; grad(3,4) =  one   
        grad(4,1) = this%Y(1);   grad(4,2) =  this%Y(2);   grad(4,3) = zero    ; grad(4,4) =  zero   

    end subroutine

end module
