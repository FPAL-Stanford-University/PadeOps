module PentadiagonalSolver
    use kind_parameters, only: rkind
    use constants,      only: one, half, zero, two, four, six
    implicit none
    
    real(rkind), dimension(:,:), allocatable :: penta
    integer                                  :: n

    ! 10th order first derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha10d1=  1.0_rkind /  2.0_rkind
    real(rkind), parameter :: beta10d1 =  1.0_rkind / 20.0_rkind
    real(rkind), parameter :: a10d1    =( 17.0_rkind / 12.0_rkind) / 2.0_rkind
    real(rkind), parameter :: b10d1    =(101.0_rkind /150.0_rkind) / 4.0_rkind
    real(rkind), parameter :: c10d1    =(  1.0_rkind /100.0_rkind) / 6.0_rkind
    
    ! Set the scheme for the edge nodes (Ref. for notation: Lele - JCP paper)
    real(rkind), parameter                   :: alpha           =   3._rkind
    real(rkind), parameter                   :: p               = -17._rkind / 6._rkind
    real(rkind), parameter                   :: q               =   3._rkind / 2._rkind
    real(rkind), parameter                   :: r               =   3._rkind / 2._rkind
    real(rkind), parameter                   :: s               = - 1._rkind / 6._rkind
    
   
    ! Calculate the corressponding weights
    ! Step 1: Assign the interior scheme
    real(rkind), parameter                   :: q_hat           = a10d1
    real(rkind), parameter                   :: r_hat           = b10d1
    real(rkind), parameter                   :: s_hat           = c10d1
    real(rkind), parameter                   :: alpha_hat       = alpha10d1
    real(rkind), parameter                   :: beta_hat        = beta10d1
     
    ! Step 2: Assign the scheme at node 2 to be Standard Pade (4th Order)
    real(rkind), parameter                   :: q_p             = 3._rkind/4._rkind
    real(rkind), parameter                   :: alpha_p         = 1._rkind/4._rkind
   
    ! Step 3: Get the scheme at node 4
    real(rkind), parameter                   :: alpha_ppp       = (8*r_hat - 175*s_hat)/(18*r_hat - 550*s_hat)
    real(rkind), parameter                   :: beta_ppp        = (1._rkind/20._rkind)*(-3 + 8*alpha_ppp)
    real(rkind), parameter                   :: q_ppp           = (1._rkind/12._rkind)*(12 - 7*alpha_ppp) 
    real(rkind), parameter                   :: r_ppp           = (1._rkind/600._rkind)*(568*alpha_ppp - 183)  
    real(rkind), parameter                   :: s_ppp           = (1._rkind/300._rkind)*(9*alpha_ppp - 4) 

    ! Step 4: Get the scheme at node 3
    real(rkind), parameter                   :: alpha_pp        = ((17*(s*(r_hat + 2*s_hat) - q*(q_hat + r_hat &
                                                                + s_hat)))/(72*(q + s)*(q_hat + r_hat - s_hat*(q_ppp/s_ppp &
                                                                - 1))) - 8._rkind/9._rkind)/((19*(s*(r_hat + 2*s_hat) - q*(q_hat &
                                                                + r_hat + s_hat)))/(24*(q + s)*(q_hat + r_hat - s_hat*(q_ppp/s_ppp &
                                                                - 1))) - 1._rkind/3_rkind)

    
    real(rkind), parameter                   :: beta_pp         = (1._rkind/12._rkind)*(-1 + 3*alpha_pp)
    real(rkind), parameter                   :: q_pp            = (2._rkind/18._rkind)*(8 - 3*alpha_pp) 
    real(rkind), parameter                   :: r_pp            = (1._rkind/72._rkind)*(-17 + 57*alpha_pp)
    real(rkind), parameter                   :: s_pp            = 0._rkind

    ! Step 5: Get the weights
    real(rkind), parameter                   :: w1              = (q_hat + 2*r_hat + 3*s_hat)/(q + s)
    real(rkind), parameter                   :: w2              = (1/q_p)*(r_hat + s_hat*(1 + q_ppp/s_ppp) - r*(q_hat &
                                                                + 2*r_hat + 3*s_hat)/(q + s) )
    
    real(rkind), parameter                   :: w3              = (q_hat + r_hat + s_hat*(1 - q_ppp/s_ppp))/(r_pp) 
    real(rkind), parameter                   :: w4              = s_hat/s_ppp 

    contains

    subroutine ComputePenta(bc1,bcn)
        integer, intent(in) :: bc1, bcn
        integer             :: i
    
        associate (bt => penta(:,1), b => penta(:,2), d => penta(:,3),  &
                   a => penta(:,4), at => penta(:,5),                   &
                   e => penta(:,6), obc => penta(:,7),                  &
                   f => penta(:,8), g => penta(:,9)                     &
                   )
    
            at = beta_hat 
            bt = beta_hat
            a  = alpha_hat
            b  = alpha_hat
            d  = one

            select case (bc1) 
            case(0)
                bt(1) = w1*zero
                b (1) = w1*zero
                d (1) = w1*one
                a (1) = w1*alpha
                at(1) = w1*zero

                bt(2) = w2*zero
                b (2) = w2*alpha_p
                d (2) = w2*one
                a (2) = w2*alpha_p
                at(2) = w2*zero

                bt(3) = w3*beta_pp
                b (3) = w3*alpha_pp
                d (3) = w3*one
                a (3) = w3*alpha_pp
                at(3) = w3*beta_pp

                bt(4) = w4*beta_ppp
                b (4) = w4*alpha_ppp
                d (4) = w4*one
                a (4) = w4*alpha_ppp
                at(4) = w4*beta_ppp
            case(1)

            end select
            
            select case (bcn) 
            case(0)
                bt(n  ) = w1*zero
                b (n  ) = w1*alpha
                d (n  ) = w1*one
                a (n  ) = w1*zero
                at(n  ) = w1*zero

                bt(n-1) = w2*zero
                b (n-1) = w2*alpha_p
                d (n-1) = w2*one
                a (n-1) = w2*alpha_p
                at(n-1) = w2*zero
                
                bt(n-2) = w3*beta_pp
                b (n-2) = w3*alpha_pp
                d (n-2) = w3*one
                a (n-2) = w3*alpha_pp
                at(n-2) = w3*beta_pp

                bt(n-3) = w4*beta_ppp
                b (n-3) = w4*alpha_ppp
                d (n-3) = w4*one
                a (n-3) = w4*alpha_ppp
                at(n-3) = w4*beta_ppp
            
            case(1)

            end select

            ! Step 1
            obc(1) = one/d(1)

            ! Step 2
            obc(2) = one/(d(2) - b(2)*a(1)*obc(1))

            ! Step 3
            e(1) = a(1)
            f(2) = b(2)*obc(1)
            
            do i = 3,n
                g(i) = bt(i)*obc(i-2)
                e(i-1) = a(i-1) - f(i-1)*at(i-2)
                f(i) = (b(i) - g(i)*e(i-2))*obc(i-1)
                obc(i) = one/(d(i) - f(i)*e(i-1) - g(i)*at(i-2))
            end do 

        end associate  
           
    end subroutine
    

    subroutine SolvePenta(y)
        real(rkind), dimension(n), intent(inout) :: y
        integer :: i
         

        associate (bt => penta(:,1), b => penta(:,2), d => penta(:,3),  &
                   a => penta(:,4), at => penta(:,5),                   &
                   e => penta(:,6), obc => penta(:,7),                  &
                   f => penta(:,8), g => penta(:,9)                     &
                   )
       
        ! Step 1
        y(2) = y(2) - f(2)*y(1)
        do i = 3,n
            y(i) = y(i) - g(i)*y(i-2) - f(i)*y(i-1)
        end do 

        ! Step 2
        y(n) = y(n)*obc(n)
        y(n-1) = (y(n-1) - e(n-1)*y(n))*obc(n-1)
        do i = n-2,1,-1
            y(i) = (y(i) - at(i)*y(i+2) - e(i)*y(i+1))*obc(i)
        end do 

        end associate


    end subroutine


    subroutine ComputeD1RHS(y,rhs,dx)
        real(rkind), dimension(n), intent(in)  :: y
        real(rkind),               intent(in)  :: dx
        real(rkind), dimension(n), intent(out) :: rhs
        real(rkind) :: a, b, c, d
        
        a = q_hat/dx  
        b = r_hat/dx 
        c = s_hat/dx

        rhs(5:n-4) = ( a*(y(6:n-3) - y(4:n-5)) &
                   +   b*(y(7:n-2) - y(3:n-6)) &
                   +   c*(y(8:n-1) - y(2:n-7)) ) 
       
        a = w4*q_ppp/dx  
        b = w4*r_ppp/dx 
        c = w4*s_ppp/dx
        rhs(4  ) =   a*(y(5  ) - y(3  )) &
                 +   b*(y(6  ) - y(2  )) &
                 +   c*(y(7  ) - y(1  )) 
        
        rhs(n-3) =   a*(y(n-2) - y(n-4)) &
                 +   b*(y(n-1) - y(n-5)) &
                 +   c*(y(n  ) - y(n-6))  
        

        a = w3*q_pp/dx  
        b = w3*r_pp/dx 
        rhs(3  ) = ( a*(y(4) - y(2)) &
                 +   b*(y(5) - y(1)) )
        
        rhs(n-2) = ( a*(y(n-1) - y(n-3)) &
                 +   b*(y(n  ) - y(n-4)) )
        
        a = w2*q_p/dx
        rhs(2) = a*(y(3) - y(1))
        rhs(n-1) = a*(y(n) - y(n-2))
       
        a = w1*( p/dx)
        b = w1*( q/dx)
        c = w1*( r/dx)
        d = w1*( s/dx)
        rhs(1  ) =  a*y(1) + b*y(  2) + c*y(  3) + d*y(  4) 
        rhs(n  ) = -a*y(n) - b*y(n-1) - c*y(n-2) - d*y(n-3)
        
    end subroutine

end module

program testPenta
    use kind_parameters, only: rkind
    use constants, only : pi, two
    use PentadiagonalSolver, only: n, penta,ComputeD1RHS, computePenta, SolvePenta
    implicit none
    integer :: i
    real(rkind), allocatable, dimension(:) :: x, y, dydx, rhs
    real(rkind) :: dx, omega

    n = 512
    !do while (n .le. 16 )
        omega = 1._rkind
        dx = two*pi/(n-1)

        allocate (x(n), y(n), dydx(n), rhs(n))
        do i = 1,n
            x(i) = (i-1)*dx
        end do 
        y = sin(omega*x)
        dydx = omega*cos(omega*x) 

        allocate (penta(n,9))

        call ComputePenta(0,0)

        call ComputeD1RHS(y,rhs,dx)
        call SolvePenta(rhs)

        print*, "Max Global Error (n =",n,"):", maxval(abs(rhs - dydx)), maxloc(abs(rhs - dydx))

        deallocate (x, y, dydx, rhs)
        deallocate (penta)
        
    !    n = n*2
    !end do 
end program
