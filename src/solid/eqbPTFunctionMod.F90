module eqbPTFunctionMod
    use kind_parameters, only: rkind
    use NewtonSolverMod, only: newton_functor
    use SolidMod,        only: solid
    implicit none

    type, extends(newton_functor) :: eqbPTFunction    ! only 2 materials for now
        integer                     :: ns
        type(solid), dimension(2)   :: material
        real(rkind), dimension(2)   :: VF0, Y
        real(rkind), dimension(9,2) :: g0
        real(rkind)                 :: emix
    contains
        procedure :: evaluate
        procedure :: gradient
    end type

    interface eqbpTFunction
        module procedure init
    end interface
contains

    function init(ns,material,VF0,g0,Y,emix) result(this)
        type(eqbpTFunction) :: this
        integer,                     intent(in) :: ns
        type(solid), dimension(:),   intent(in) :: material
        real(rkind), dimension(:),   intent(in) :: VF0, Y
        real(rkind), dimension(:,:), intent(in) :: g0
        real(rkind),                 intent(in) :: emix

        this%ns = ns
        this%material = material
        this%VF0(:) = VF0(:)
        this%g0(:,:) = g0(:,:)
        this%Y(:)  = Y(:)
        this%emix = emix

        this%niters = 10
        this%tolerance = real(1.D-16,rkind)
    end function

    subroutine evaluate(this, x, y)
        use constants, only: one
        class(eqbpTFunction),              intent(in)  :: this
        real(rkind), dimension(:),         intent(in)  :: x
        real(rkind), dimension(size(x,1)), intent(out) :: y

        real(rkind), dimension(this%ns) :: energy, VF, p, T
        integer :: m

        do m = 1, this%ns
            energy(m) = x(m); VF(m) = x(this%ns+m)
            call this%material(m)%eos%get_pT_from_energyVF(this%VF0(m), this%g0(:,m), energy(m), VF(m), p(m), T(m))
        enddo
        ! set elements of jacobian - obviously works for 2 materials only
        y(1) = p(1) - p(2)
        y(2) = T(1) - T(2)
        y(3) = VF(1) + VF(2) - one
        y(4) = this%Y(1)*energy(1) + this%Y(2)*energy(2) - this%emix

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
        grad(4,1) = one    ;   grad(4,2) =  one    ;   grad(4,3) = zero    ; grad(4,4) =  zero   

    end subroutine

end module
