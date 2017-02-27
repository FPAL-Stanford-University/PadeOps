module NewtonSolverProg
    use kind_parameters, only: rkind
    use NewtonSolverMod, only: newton_functor
    implicit none

    type, extends(newton_functor) :: myfunction
        real(rkind) :: a, b, c
    contains
        procedure :: evaluate
        procedure :: gradient
    end type

    interface myfunction
        module procedure init
    end interface
contains

    function init(a,b,c) result(this)
        type(myfunction) :: this
        real(rkind), intent(in) :: a, b, c

        this%a = a
        this%b = b
        this%c = c

        this%niters = 10
        this%tolerance = real(1.D-16,rkind)
    end function

    subroutine evaluate(this, x, y)
        class(myfunction),                 intent(inout)  :: this
        real(rkind), dimension(:),         intent(in)     :: x
        real(rkind), dimension(size(x,1)), intent(out)    :: y
        
        y = this%a*x**2 + this%b*x + this%c
    end subroutine

    subroutine gradient(this, x, grad)
        use constants, only: zero
        class(myfunction),                           intent(in)  :: this
        real(rkind), dimension(:),                   intent(in)  :: x
        real(rkind), dimension(size(x,1),size(x,1)), intent(out) :: grad

        integer :: i

        grad = zero

        do i = 1,size(x,1)
            grad(i,i) = 2*this%a*x(i) + this%b
        end do
    end subroutine

end module

program test_NewtonSolver
    use kind_parameters,  only: rkind
    use constants,        only: zero,half,one,two,three
    use NewtonSolverProg, only: myfunction
    implicit none

    type(myfunction) :: f
    real(rkind), dimension(2) :: x, y

    f = myfunction(one,zero,-one) ! x**2 - 1

    x = half    ! Initial guess
    y = [zero, three]    ! Final function value

    print*, "Target value: ", y
    print*, "Exact   x: ", [one, two]
    print*, "Initial x: ", x

    call f%newton_solve(x,y)
    
    print*, "Final x: ", x
    call f%evaluate(x,y)
    print*, "Final value: ", y
end program
