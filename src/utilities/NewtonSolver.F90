module NewtonSolverMod
    use kind_parameters, only: rkind
        implicit none

        ! Unary functor y = f(x) used in the Newton Solver
        type, abstract :: newton_functor
        contains
            procedure(evaluate_interface), deferred :: evaluate
            procedure(gradient_interface), deferred :: gradient
            procedure                               :: get_residual
            procedure                               :: newton_solve
        end type

        interface 
            subroutine evaluate_interface(this, x, y)
                import :: rkind
                import :: newton_functor
                class(newton_functor),             intent(in)  :: this
                real(rkind), dimension(:),         intent(in)  :: x
                real(rkind), dimension(size(x,1)), intent(out) :: y
            end subroutine

            subroutine gradient_interface(this, x, grad)
                import :: rkind
                import :: newton_functor
                class(newton_functor),                       intent(in)  :: this
                real(rkind), dimension(:),                   intent(in)  :: x
                real(rkind), dimension(size(x,1),size(x,1)), intent(out) :: grad
            end subroutine
        end interface
contains

    subroutine newton_solve(f,x,y)
        use constants, only: eps, half, one
        use exits,     only: GracefulExit
        class(newton_functor),             intent(in)    :: f
        real(rkind), dimension(:),         intent(inout) :: x
        real(rkind), dimension(size(x,1)), intent(in)    :: y
        
        real(rkind), dimension(size(x,1))           :: dx, dx_new, x_new, f0
        real(rkind), dimension(size(x,1),size(x,1)) :: gradf
        integer,     dimension(size(x,1))           :: ipiv

        real(rkind) :: residual, residual_new, alpha = half, t = one, tol = real(1.D-16,rkind)
        integer     :: n, iters, niters = 50
        
        n = size(x,1)

        ! Get original residual
        call f%get_residual(x,y,ipiv,gradf,dx,f0,residual)

        iters = 0

        do while ( (iters < niters) .and. (abs(residual) > tol) )
            ! Backtracking line search
            t = one
            x_new = x + t*dx

            ! Get new residual
            call f%get_residual(x_new,y,ipiv,gradf,dx_new,f0,residual_new)

            do while ( (abs(residual_new) >= abs(residual)) .and. (t > eps) )
                if (iters > (niters-10)) then
                    print '(A,I0,3(A,ES15.5))', 'iters = ', iters, ', t = ', t, ', residual_new = ', residual_new, ', residual = ', residual
                end if

                t = alpha*t
                x_new = x + t*dx

                ! Get new residual
                call f%get_residual(x_new,y,ipiv,gradf,dx_new,f0,residual_new)
            end do

            x = x_new
            dx = dx_new
            residual = residual_new

            iters = iters + 1

            if ((iters >= niters) .or. (t <= eps)) then
                call GracefulExit('Newton solve did not converge',6382)
            end if
        end do

    end subroutine

    subroutine get_residual(f,x,y,ipiv,gradf,dx,f0,residual)
        use exits,     only: GracefulExit
        class(newton_functor),                       intent(in)  :: f
        real(rkind), dimension(:),                   intent(in)  :: x
        real(rkind), dimension(size(x,1)),           intent(in)  :: y
        integer,     dimension(size(x,1)),           intent(out) :: ipiv
        real(rkind), dimension(size(x,1)),           intent(out) :: dx, f0
        real(rkind), dimension(size(x,1),size(x,1)), intent(out) :: gradf
        real(rkind),                                 intent(out) :: residual

        integer :: n, info

        n = size(x,1)

        ! Get initial function value
        call f%evaluate(x,f0)

        ! Get Newton step dx
        call f%gradient(x, gradf)
        dx = y - f0
        call dgesv(n, 1, gradf, n, ipiv, dx, n, info)
        if (info /= 0) call GracefulExit("Error in getting Newton step (dgesv call)",54762)

        ! Compute residual
        residual = sum( (y-f0)*dx )

    end subroutine

end module NewtonSolverMod
