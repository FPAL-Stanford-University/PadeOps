module Tridiagsolver
    use kind_parameters, only: rkind
    use constants,      only: one, half, zero, two, three, four, five, six
    implicit none
    
    real(rkind), dimension(:,:), allocatable :: Tridiag
    integer                                  :: n

    contains

    subroutine ComputeTridiag(bc1,bcn)
        integer, intent(in) :: bc1, bcn
        integer             :: i
    
        associate (a => Tridiag(:,1), b => Tridiag(:,2), c => Tridiag(:,3), cp => Tridiag(:,4), den => Tridiag(:,5))
    
            a  = one/three
            b  = one
            c  = one/three

            select case (bc1) 
            case(0)
                a (2) = one/four
                c (2) = one/four
                
                a (1) = zero
                c (1) = two
            case(1)

            end select
            
            select case (bcn) 
            case(0)
                a (n-1) = one/four
                c (n-1) = one/four
                
                a (n  ) = two
                c (n  ) = zero
            case(1)

            end select

            cp(1) = c(1)/b(1)
            do i = 2,n-1
                cp(i) = c(i)/(b(i) - a(i)*cp(i-1))
            end do
            
            den(1) = one/b(1)
            den(2:n) = one/(b(2:n) - a(2:n)*cp(1:n-1))

        end associate  
           
    end subroutine
    

    subroutine SolveTridiag(d)
        real(rkind), dimension(n), intent(inout) :: d
        integer :: i

        !associate (a => Tridiag(:,1), b => Tridiag(:,2), c => Tridiag(:,3), cp => Tridiag(:,4), den => Tridiag(:,5))
       
        d(1) = d(1)*Tridiag(1,5)
        do i = 2,n
            d(i) = (d(i) - Tridiag(i,1)*d(i-1))*Tridiag(i,5)
        end do
        
        do i = n-1,1,-1
            d(i) = d(i) - Tridiag(i,4)*d(i+1)
        end do 
        
        !end associate
    end subroutine


    subroutine ComputeD1RHS(y,rhs,dx)
        real(rkind), dimension(n), intent(in)  :: y
        real(rkind),               intent(in)  :: dx
        real(rkind), dimension(n), intent(out) :: rhs
        real(rkind) :: a, b, c
        
        a = (14._rkind/9._rkind)/(two *dx)
        b = ( 1._rkind/9._rkind)/(four*dx)

        rhs(3:n-2) = a*(y(4:n-1) - y(2:n-3)) &
                   + b*(y(5:n  ) - y(1:n-4)) 
        
        
        a = three/(four*dx)
        rhs(2  ) = a*(y(3) - y(1  ))
        rhs(n-1) = a*(y(n) - y(n-2)) 
        

        rhs(1  ) = (one/dx)*(-(five*half)*y(1) + two*y(2  ) + half*y(3  ))
        rhs(n  ) = (one/dx)*( (five*half)*y(n) - two*y(n-1) - half*y(n-2))
        
    end subroutine

end module

program testTridiag
    use kind_parameters, only: rkind
    use constants, only : pi, one, zero, two
    use TridiagSolver, only: n, tridiag,ComputeD1RHS, computeTridiag, solveTridiag
    use timer,           only: tic, toc
    implicit none
    integer :: i
    real(rkind), allocatable, dimension(:) :: x, y, dydx, rhs
    real(rkind) :: dx, omega

    n = 64
    do while (n .le.  65536 )
       omega = 12._rkind
       dx = two*pi/(n-1)

       allocate (x(n), y(n), dydx(n), rhs(n))
       do i = 1,n
           x(i) = (i-1)*dx
       end do 
       y = sin(omega*x)
       dydx = omega*cos(omega*x) 

       allocate (Tridiag(n,5))

       call ComputeTridiag(0,0)
       call ComputeD1RHS(y,rhs,dx)
       call SolveTridiag(rhs)
       print*, "Max Global Error (n =",n,"):", maxval(abs(rhs - dydx))
       !print*, "Max Int.   Error (n =",n,"):", maxval(abs(x(8:n-7) - dydx(8:n-7)))
       print*, " ------------------"
       deallocate (x, y, dydx, rhs)
       deallocate (Tridiag)
       
        n = n*2
    end do 
end program
