module NewtonSolverMod
    use kind_parameters, only: rkind
    use constants,       only: half
    implicit none

    ! Unary functor y = f(x) used in the Newton Solver
    type, abstract :: newton_functor
        real(rkind) :: alpha = half
        real(rkind) :: tolerance = real(1.D-14,rkind)
        integer     :: niters = 50
        logical     :: isGradPosDef = .false.
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
            class(newton_functor),             intent(inout)  :: this
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
        use constants, only: eps, one
        use exits,     only: GracefulExit
        class(newton_functor),             intent(inout)    :: f
        real(rkind), dimension(:),         intent(inout) :: x
        real(rkind), dimension(size(x,1)), intent(in)    :: y
        
        real(rkind), dimension(size(x,1))           :: dx, dx_new, x_new, f0
        real(rkind), dimension(size(x,1),size(x,1)) :: gradf
        integer,     dimension(size(x,1))           :: ipiv

        real(rkind) :: residual, residual_new, t = one
        integer     :: n, iters
        
        n = size(x,1)

        ! Get original residual
        call f%get_residual(x,y,ipiv,gradf,dx,f0,residual)

        iters = 0

        do while ( (iters < f%niters) .and. (abs(residual) > f%tolerance) )
            ! Backtracking line search
            t = one
            x_new = x + t*dx

            ! Get new residual
            call f%get_residual(x_new,y,ipiv,gradf,dx_new,f0,residual_new)

            do while ( (abs(residual_new) >= abs(residual)) .and. (t > eps) )
                if (iters > (f%niters-10)) then
                    print '(A,I0,3(A,ES15.5))', 'iters = ', iters, ', t = ', t, ', residual_new = ', residual_new, ', residual = ', residual
                end if

                t = f%alpha*t
                x_new = x + t*dx

                ! Get new residual
                call f%get_residual(x_new,y,ipiv,gradf,dx_new,f0,residual_new)
            end do

            x = x_new
            dx = dx_new
            residual = residual_new

            iters = iters + 1

            if ((iters >= f%niters) .or. (t <= eps)) then
                !print *, iters, t
                call GracefulExit('Newton solve did not converge',6382)
            end if
        end do
        !write(*,'(a,4(e19.12,1x))') '   residual = ', residual
        !write(*,'(a,4(i4,1x))') 'no. iterons = ', iters

    end subroutine

    subroutine get_residual(f,x,y,ipiv,gradf,dx,f0,residual)
        use exits,     only: GracefulExit
        class(newton_functor),                       intent(inout)  :: f
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
        !write(*,*) '   y  = ', y
        !write(*,'(a,4(e19.12,1x))') '   dx = ', dx
        !write(*,*) '   f0 = ', f0
        !write(*,'(a,4(e19.12,1x))') '   gradf = ', gradf(1,:)
        !write(*,'(a,4(e19.12,1x))') '           ', gradf(2,:)
        !write(*,'(a,4(e19.12,1x))') '           ', gradf(3,:)
        !write(*,'(a,4(e19.12,1x))') '           ', gradf(4,:)
        call dgesv(n, 1, gradf, n, ipiv, dx, n, info)
        if (info /= 0) call GracefulExit("Error in getting Newton step (dgesv call)",54762)

        ! Compute residual
        if(f%isGradPosDef) then
           residual = sum( (y-f0)*dx )
        else
           residual = sum( (y-f0)**2 )
        endif
        !write(*,'(a,4(e19.12,1x))') '   dx = ', dx
        !write(*,'(a,4(e19.12,1x))') '   y  = ', y
        !write(*,'(a,4(e19.12,1x))') '   dx = ', dx
        !write(*,'(a,4(e19.12,1x))') '   f0 = ', f0
        !write(*,'(a,4(e19.12,1x))') '   residual = ', residual

    end subroutine

end module NewtonSolverMod
