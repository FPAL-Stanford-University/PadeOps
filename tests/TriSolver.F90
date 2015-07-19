module Tridiagsolver
    use kind_parameters, only: rkind
    use constants,      only: one, half, zero, two, three, four, five, six
    use cd06_constants,     only: alpha06d1,a06d1,b06d1   
    implicit none
    
    real(rkind), dimension(:,:), allocatable :: Tridiag
    integer                                  :: n
    
    ! Set the scheme for the edge nodes (Ref. for notation: Lele - JCP paper)
    real(rkind), parameter                   :: alpha =  3._rkind
    real(rkind), parameter                   :: p     = 17._rkind / 6._rkind
    real(rkind), parameter                   :: q     =  3._rkind / 2._rkind
    real(rkind), parameter                   :: r     =  3._rkind / 2._rkind
    real(rkind), parameter                   :: s     = -1._rkind / 6._rkind 


    !!  Calculate the corressponding weights 
    ! Step 1: Assign the interior scheme
    real(rkind), parameter                   :: qhat        = a06d1
    real(rkind), parameter                   :: rhat        = b06d1
    real(rkind), parameter                   :: alpha_hat   = alpha06d1

    ! Step 2: Assign the scheme at node 2 to be Standard Pade (4th Order)
    real(rkind), parameter                   :: q_p          = 3._rkind/4._rkind
    real(rkind), parameter                   :: alpha_p     = 1._rkind/4._rkind
    
    ! Step 3: Get the scheme at the node 3
    real(rkind), parameter                   :: alpha_pp = ((40*alpha_hat - 1)*q  + 7*(4*alpha_hat &
                                                         -  1)*s)/(16*(alpha_hat + 2)*q + 8*(1      &
                                                         -  4*alpha_hat)*s)
    real(rkind), parameter                   :: q_pp     = (1._rkind/3._rkind)*(alpha_pp + 2)
    real(rkind), parameter                   :: r_pp     = (1._rkind/12._rkind)*(4*alpha_pp - 1)
    real(rkind), parameter                   :: s_pp     =  0._rkind

    ! Step 4: Get the weights
    real(rkind), parameter                   :: w1 = (2*alpha_hat + 1)/(2*(q + s))
    real(rkind), parameter                   :: w2 = ((8*alpha_hat + 7)*q - 6*(2*alpha_hat + 1)*r &
                                                   + (8*alpha_hat + 7)*s)/(9*(q + s))
    real(rkind), parameter                   :: w3 = (4*(alpha_hat + 2)*q + 2*(1 - 4*alpha_hat)*s) &
                                                   / (9*(q + s))

    contains

    subroutine ComputeTridiag(bc1,bcn)
        integer, intent(in) :: bc1, bcn
        integer             :: i
    
        associate (a => Tridiag(:,1), b => Tridiag(:,2), c => Tridiag(:,3), cp => Tridiag(:,4), den => Tridiag(:,5))
    
            a  = alpha_hat
            b  = one
            c  = alpha_hat 

            select case (bc1) 
            case(0)
                a (1) = w1*zero
                a (2) = w2*alpha_p
                a (3) = w3*alpha_pp 

                b (1) = w1*one
                b (2) = w2*one
                b (3) = w3*one

                c (1) = w1*alpha
                c (2) = w2*alpha_p
                c (3) = w3*alpha_pp
                
            case(1)

            end select
            
            select case (bcn) 
            case(0)
                c (n  ) = w1*zero
                c (n-1) = w2*alpha_p
                c (n-2) = w3*alpha_pp 

                b (n  ) = w1*one
                b (n-1) = w2*one
                b (n-2) = w3*one 

                a (n  ) = w1*alpha
                a (n-1) = w2*alpha_p
                a (n-2) = w3*alpha_pp 
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
        real(rkind) :: a, b, c, d
        
        a = (qhat)/(dx)
        b = (rhat)/(dx)

        rhs(4:n-3) = a*(y(5:n-2) - y(3:n-4)) &
                   + b*(y(6:n-1) - y(2:n-5)) 
        
      
        a = w3*q_pp/dx
        b = w3*r_pp/dx 
        ! 3rd node and n-2 th node
        rhs(3) = a*(y(4) - y(2)) + b*(y(5) - y(1)) 
        rhs(n-2) = a*(y(n-1) - y(n-3)) + b*(y(n) - y(n-4)) 
        
      
        ! 2nd node and n-1st node 
        a = w2*q_p/dx
        rhs(2  ) = a*(y(3) - y(1  ))
        rhs(n-1) = a*(y(n) - y(n-2)) 
        

        ! 1st and nth node
        a = w1*(-p/dx)
        b = w1*( q/dx)
        c = w1*( r/dx)
        d = w1*( s/dx)
        rhs(1  ) =  a*y(1) + b*y(  2) + c*y(  3) + d*y(  4) 
        rhs(n  ) = -a*y(n) - b*y(n-1) - c*y(n-2) - d*y(n-3)
        
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

    n = 256
    omega = 1._rkind
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
    deallocate (x, y, dydx, rhs)
    deallocate (Tridiag)
end program
