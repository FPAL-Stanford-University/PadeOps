module PentadiagonalSolver
    use kind_parameters, only: rkind
    use constants,      only: one, zero, two
    implicit none
    
    real(rkind), dimension(:,:), allocatable :: penta
    integer                                  :: n

    contains

    subroutine ComputePenta(bc1,bcn)
        integer, intent(in) :: bc1, bcn
        integer             :: i
    
        associate (e => penta(:,1), c => penta(:,2), d => penta(:,3), &
                   a => penta(:,4), b => penta(:,5),                  &
                   rho => penta(:,6), obp => penta(:,7),              &
                   sig => penta(:,8), phi => penta(:,9)               &
                   )
    
            e = one/20._rkind
            b = one/20._rkind
            c = one/2._rkind
            a = one/2._rkind
            d = one

            select case (bc1) 
            case(0)
                e(4) = 0.451390625_rkind / 9.38146875_rkind
                b(4) = 0.451390625_rkind / 9.38146875_rkind 
                c(4) = 4.632718750_rkind / 9.38146875_rkind
                a(4) = 4.632718750_rkind / 9.38146875_rkind

                e(3) = 0.2964375_rkind / 10.67175_rkind 
                b(3) = 0.2964375_rkind / 10.67175_rkind 
                c(3) = 4.74375_rkind / 10.67175_rkind 
                a(3) = 4.74375_rkind / 10.67175_rkind
                
                e(2) = zero 
                b(2) = zero
                c(2) = one/4._rkind
                a(2) = one/4._rkind
                
                e(1) = zero 
                b(1) = zero
                c(1) = zero
                a(1) = two
            case(1)

            end select
            
            select case (bcn) 
            case(0)
                e(n-3) = 0.451390625_rkind / 9.38146875_rkind
                b(n-3) = 0.451390625_rkind / 9.38146875_rkind 
                c(n-3) = 4.632718750_rkind / 9.38146875_rkind
                a(n-3) = 4.632718750_rkind / 9.38146875_rkind

                e(n-2) = 0.2964375_rkind / 10.67175_rkind 
                b(n-2) = 0.2964375_rkind / 10.67175_rkind 
                c(n-2) = 4.74375_rkind / 10.67175_rkind 
                a(n-2) = 4.74375_rkind / 10.67175_rkind
                
                e(n-1) = zero 
                b(n-1) = zero
                c(n-1) = one/4._rkind
                a(n-1) = one/4._rkind
                
                e(n  ) = zero 
                b(n  ) = zero
                c(n  ) = zero
                a(n  ) = two
            case(1)

            end select

            e = 1._rkind 
            c = 2._rkind
            d = 3._rkind
            a = 4._rkind
            b = 5._rkind

            a = [2,2,1,5,-7,3,-1,4,5,0]*one
            b = [1,5,-2,1,5,2,4,-3,0,0]*one
            c = [0,3,2,1,2,1,2,1,-2,4]*one
            e = [0,0,1,3,1,5,2,2,2,-1]*one
            d = [1,2,3,-4,5,6,7,-1,1,8]*one

            ! Step 3
            obp(n) = one/d(n)
            sig(n) = c(n)*obp(n)
            phi(n) = e(n)*obp(n)

            ! Step 4
            rho(n-1) = a(n-1)
            obp(n-1) = one/(d(n-1) - sig(n)*rho(n-1))
            sig(n-1) = (c(n-1) - phi(n)*rho(n-1))*obp(n-1)
            phi(n-1) = e(n-1)*obp(n-1)


            ! Step 5
            do i = n-2,3,-1
                rho(i) = a(i) - sig(i+2)*b(i)
                obp(i) = one/(d(i) - phi(i+2)*b(i) - sig(i+1)*rho(i))
                sig(i) = (c(i) - phi(i+1)*rho(i))*obp(i)
                phi(i) = e(i)*obp(i)
            end do 
            
            rho(2) = a(2) - sig(4)*b(2)
            obp(2) = one/(d(2) - phi(4)*b(2) - sig(3)*rho(2))
            sig(2) = (c(2) - phi(4)*rho(2))*obp(2)

            rho(1) = a(1) - sig(3)*b(1)
            obp(1) = one/(d(1) - phi(3)*b(1) - sig(2)*rho(1))


        end associate  
           
    end subroutine
    

    subroutine PTRANS_II(y)
        real(rkind), dimension(n), intent(inout) :: y
        integer :: i
         

        y = [8,33,8,24,29,98,99,17,57,108]
        associate (e => penta(:,1), c => penta(:,2), d => penta(:,3), &
                   a => penta(:,4), b => penta(:,5),                  &
                   rho => penta(:,6), obp => penta(:,7),              &
                   sig => penta(:,8), phi => penta(:,9)               &
                   )
       
        print*, "PSI:", one/obp 
        print*, "RHO:", rho
        print*, "PHI:", phi
        print*, "SIG:", sig
        ! Step 3
        y(n  )  = y(n)*obp(n)
        
        ! Step 4
        y(n-1)  = (y(n-1) - y(n)*rho(n-1))*obp(n-1)

        ! Step 5
        do i = n-2,1,-1
            y(i) = (y(i) - y(i+1)*rho(i) - y(i+2)*b(i))*obp(i)
        end do 

        ! Step 6
        y(2) = y(2) - sig(2)*y(1)
        do i = 3,n
            y(i) = y(i) - sig(i)*y(i-1) - phi(i)*y(i-2)
        end do 

        end associate


    end subroutine



subroutine Algorithm3(y)
    real(rkind), dimension(n), intent(inout) :: y


end subroutine 

end module

program testPenta
    use kind_parameters, only: rkind
    use constants, only : pi, one, zero, two
    use PentadiagonalSolver, only: n, penta, computePenta, PTRANS_II
    implicit none
    integer :: i
    real(rkind), allocatable, dimension(:) :: x, y, dydx
    real(rkind) :: dx, omega

    n = 10
    omega = one
    dx = two*pi/n

    allocate (x(n), y(n), dydx(n))
    do i = 1,n
        x(i) = (i-1)*dx
    end do 
    y = sin(omega*x)
    dydx = omega*cos(omega*x) 

    y = one
    dydx = zero
    allocate (penta(n,9))

    call ComputePenta(1,1)
    call PTRANS_II(y)

    !print*, "Max Global Error:", maxval(abs(y - dydx))
    !print*, "Max interior Error:", maxval(abs(y(2:n-1) - dydx(2:n-1)))
    print*, y


end program
