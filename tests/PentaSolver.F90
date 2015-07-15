module PentadiagonalSolver
    use kind_parameters, only: rkind
    use constants,      only: one, half, zero, two, four, six
    implicit none
    
    real(rkind), dimension(:,:), allocatable :: penta
    integer                                  :: n

    contains

    subroutine ComputePenta(bc1,bcn)
        integer, intent(in) :: bc1, bcn
        integer             :: i
    
        associate (bt => penta(:,1), b => penta(:,2), d => penta(:,3),  &
                   a => penta(:,4), at => penta(:,5),                   &
                   e => penta(:,6), obc => penta(:,7),                  &
                   f => penta(:,8), g => penta(:,9)                     &
                   )
    
            at = one/20._rkind * 9._rkind
            bt = one/20._rkind * 9._rkind
            a  = one/2._rkind * 9._rkind
            b  = one/2._rkind * 9._rkind
            d  = one * 9._rkind

            select case (bc1) 
            case(0)
                at(4) = 0.451390625_rkind / 9.38146875_rkind * 9.38146875_rkind
                bt(4) = 0.451390625_rkind / 9.38146875_rkind * 9.38146875_rkind 
                a (4) = 4.632718750_rkind / 9.38146875_rkind * 9.38146875_rkind
                b (4) = 4.632718750_rkind / 9.38146875_rkind * 9.38146875_rkind
                
                at(3) = 0.2964375_rkind / 10.67175_rkind * 10.67175_rkind 
                bt(3) = 0.2964375_rkind / 10.67175_rkind * 10.67175_rkind 
                a (3) =     4.743_rkind / 10.67175_rkind * 10.67175_rkind 
                b (3) =     4.743_rkind / 10.67175_rkind * 10.67175_rkind
                
                at(2) = zero 
                bt(2) = zero
                a (2) = one/4._rkind * 7.783125_rkind
                b (2) = one/4._rkind * 7.783125_rkind
                
                at(1) = zero 
                bt(1) = zero
                a (1) = two  * 4.725_rkind
                b (1) = zero * 4.725_rkind
            case(1)

            end select
            
            select case (bcn) 
            case(0)
                at(n-3) = 0.451390625_rkind / 9.38146875_rkind * 9.38146875_rkind
                bt(n-3) = 0.451390625_rkind / 9.38146875_rkind * 9.38146875_rkind 
                a (n-3) = 4.632718750_rkind / 9.38146875_rkind * 9.38146875_rkind
                b (n-3) = 4.632718750_rkind / 9.38146875_rkind * 9.38146875_rkind
                
                at(n-2) = 0.2964375_rkind / 10.67175_rkind * 10.67175_rkind 
                bt(n-2) = 0.2964375_rkind / 10.67175_rkind * 10.67175_rkind
                a (n-2) =     4.743_rkind / 10.67175_rkind * 10.67175_rkind 
                b (n-2) =     4.743_rkind / 10.67175_rkind * 10.67175_rkind
                
                at(n-1) = zero 
                bt(n-1) = zero
                a (n-1) = one/4._rkind * 7.783125_rkind
                b (n-1) = one/4._rkind * 7.783125_rkind
                
                at(n  ) = zero
                bt(n  ) = zero
                a (n  ) = zero * 4.725_rkind
                b (n  ) = two  * 4.725_rkind
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
    

    subroutine PTRANS_II(y)
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
        real(rkind) :: a, b, c
        
        a = (17._rkind/12._rkind)/(two*dx)
        b = (101._rkind/150._rkind)/(four*dx)
        c = (one/100._rkind) / (six*dx)

        rhs(5:n-4) = ( a*(y(6:n-3) - y(4:n-5)) &
                   +   b*(y(7:n-2) - y(3:n-6)) &
                   +   c*(y(8:n-1) - y(2:n-7)) ) * 9._rkind
        
        
        a = ( 6.66984375_rkind / 9.38146875_rkind)/(dx)
        b = (       1.53_rkind / 9.38146875_rkind)/(dx)
        c = (      0.015_rkind / 9.38146875_rkind)/(dx)

        rhs(4) = ( a*(y(5) - y(3)) &
               +   b*(y(6) - y(2)) &
               +   c*(y(7) - y(1)) ) * 9.38146875_rkind
        
        rhs(n-3) = ( a*(y(n-2) - y(n-4)) &
                 +   b*(y(n-1) - y(n-5)) &
                 +   c*(y(n  ) - y(n-6)) ) * 9.38146875_rkind
        
        a = (     7.905_rkind / 10.67175_rkind)/(dx)
        b = (1.23515625_rkind / 10.67175_rkind)/(dx)

        rhs(3) = ( a*(y(4) - y(2)) &
               +   b*(y(5) - y(1)) ) * 10.67175_rkind
        
        rhs(n-2) = ( a*(y(n-1) - y(n-3)) &
                 +   b*(y(n  ) - y(n-4)) ) * 10.67175_rkind
        
        a = (1.5_rkind)/(two*dx)
        rhs(2) = a*(y(3) - y(1)) * 7.783125_rkind 
        
        rhs(n-1) = a*(y(n) - y(n-2)) * 7.783125_rkind 
       
        rhs(1) = ( (-2.5_rkind * y(1) + two*y(2  ) + half*y(3  ))/dx ) * 4.725_rkind
        rhs(n) = ( ( 2.5_rkind * y(n) - two*y(n-1) - half*y(n-2))/dx ) * 4.725_rkind

    end subroutine

subroutine Algorithm3(y)
    real(rkind), dimension(n), intent(inout) :: y


end subroutine 

end module

program testPenta
    use kind_parameters, only: rkind
    use constants, only : pi, one, zero, two
    use PentadiagonalSolver, only: n, penta,ComputeD1RHS, computePenta, PTRANS_II
    implicit none
    integer :: i
    real(rkind), allocatable, dimension(:) :: x, y, dydx, rhs
    real(rkind) :: dx, omega

    n = 16
    do while (n .le. 16 )
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

        ! rhs = zero
        ! rhs(1) = 6._rkind
        ! rhs(2) = -1._rkind

        ! dydx = one

        call ComputeD1RHS(y,rhs,dx)
        call PTRANS_II(rhs)

        print*, "Max Global Error (n =",n,"):", maxval(abs(rhs - dydx)), maxloc(abs(rhs - dydx))

        print*, rhs
        print*, "Conservative error: ", SUM(rhs*dx)

        deallocate (x, y, dydx, rhs)
        deallocate (penta)
        
        n = n*2
    end do 
end program
