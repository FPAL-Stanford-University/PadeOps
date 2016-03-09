! Routines specific to 8th order Compact Filter
! Periodic LU based on  Neossi Nguetchue, Abelman (Appl. Math. & Comp. 2008)

module cf90stuff

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two
    use exits,           only: GracefulExit

    implicit none

    private
    public :: cf90, alpha90, beta90, a90, b90, c90, d90, e90
    
    ! 8th order filter coefficients with 90% truncation (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha90= real(6.6624D-1, rkind)
    real(rkind), parameter :: beta90 = real(1.6688D-1, rkind)
    real(rkind), parameter :: a90    = real(9.9965D-1, rkind)
    real(rkind), parameter :: b90    = real(6.6652D-1, rkind) ! Already divided by factor of 2 
    real(rkind), parameter :: c90    = real(1.6674D-1, rkind) ! Already divided by factor of 2 
    real(rkind), parameter :: d90    = real( 4.0D-5  , rkind) ! Already divided by factor of 2 
    real(rkind), parameter :: e90    = real(-5.0D-6  , rkind) ! Already divided by factor of 2 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic filter !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
    ! First point is just identity

    ! Second point
    real(rkind), parameter :: b2_alpha90= real( 4.997D-1, rkind)
    real(rkind), parameter :: b2_a90    = real( 9.997D-1, rkind)
    real(rkind), parameter :: b2_b90    = real(4.9985D-1, rkind) ! Already divided by factor of 2 
    
    ! Third point
    real(rkind), parameter :: b3_alpha90= real(6.6624D-1, rkind)
    real(rkind), parameter :: b3_beta90 = real(1.6688D-1, rkind)
    real(rkind), parameter :: b3_a90    = real(9.9952D-1, rkind)
    real(rkind), parameter :: b3_b90    = real(6.6656D-1, rkind) ! Already divided by factor of 2 
    real(rkind), parameter :: b3_c90    = real( 1.668D-1, rkind) ! Already divided by factor of 2 
   
    ! Fourth point
    real(rkind), parameter :: b4_alpha90= real(6.6624D-1, rkind)
    real(rkind), parameter :: b4_beta90 = real(1.6688D-1, rkind)
    real(rkind), parameter :: b4_a90    = real(9.9968D-1, rkind)
    real(rkind), parameter :: b4_b90    = real(6.6652D-1, rkind) ! Already divided by factor of 2 
    real(rkind), parameter :: b4_c90    = real(1.6672D-1, rkind) ! Already divided by factor of 2 
    real(rkind), parameter :: b4_d90    = real( 4.0D-5  , rkind) ! Already divided by factor of 2 
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    type cf90
        
        private
        
        integer     :: n

        logical     :: periodic=.TRUE.

        real(rkind), allocatable, dimension(:,:) :: LU
        real(rkind), allocatable, dimension(:,:) :: penta

        contains

        procedure :: init
        procedure :: destroy

        procedure, private :: ComputeXRHS
        procedure, private :: ComputeYRHS
        procedure, private :: ComputeZRHS
        
        procedure, private :: SolveXLU
        procedure, private :: SolveYLU
        procedure, private :: SolveZLU
        
        procedure, private :: ComputeLU
        procedure, private :: ComputePenta

        procedure, private :: SolveXPenta
        procedure, private :: SolveYPenta
        procedure, private :: SolveZPenta
        
        procedure :: filter1
        procedure :: filter2
        procedure :: filter3
        
    end type



contains

    function init(this, n_, periodic_) result(ierr)
   
        class( cf90 ), intent(inout) :: this
        integer, intent(in) :: n_
        logical, intent(in) :: periodic_
        integer :: ierr
        
        this%n = n_

        this%periodic = periodic_

        if (periodic_) then 
            ! Allocate LU matrix
            if(allocated( this%LU )) deallocate( this%LU ); allocate( this%LU(n_,9) )
    
            ! Compute LU matrix
            if (n_ .GE. 10) then
                call this%ComputeLU(this%LU,beta90,alpha90,one,alpha90,beta90)
            else if (n_ == 1) then
                this%LU = one
            else
                ierr = 7
                return
            end if
    
        else 

            ! Allocate Penta matrix
            if(allocated( this%penta )) deallocate( this%penta ); allocate( this%penta(n_,11) )
  
            if (n_ .GE. 10) then             
                call this%ComputePenta()
            else if (n_ .EQ. 1) then
                this%Penta = one
            else
                ierr = 7
                return 
            end if 
            
        end if 

        ! If everything passes
        ierr = 0
    
    end function
    
    subroutine destroy(this)

        class( cf90 ), intent(inout) :: this

        ! Dellocate LU matrix.
        if(allocated( this%LU )) deallocate( this%LU )
    
        ! Dellocate penta matrix.
        if(allocated( this%penta )) deallocate( this%penta )
    
    end subroutine

    
    subroutine ComputeLU(this,LU,e,a,d,c,f) 
    
        class( cf90 ), intent(inout) :: this
        real(rkind), intent(in) :: d,a,c,e,f
        real(rkind), dimension(this%n,9), intent(out) :: LU
        integer :: i
    
        LU = 0.0_rkind
    
        associate( b=>LU(:,1), eg=>LU(:,2), k=>LU(:,3),&
                   l=>LU(:,4),  g=>LU(:,5), h=>LU(:,6),&
                   ff=>LU(:,7),  v=>LU(:,8), w=>LU(:,9))
            
            ! Step 1       
            g(1) = d
            b(2) = a/g(1)
            h(1) = c
            k(1) = f/g(1)
            w(1) = a
            v(1) = e
            l(1) = c/g(1)
            g(2) = d - b(2)*h(1)
            k(2) = -k(1)*h(1)/g(2)
            w(2) = e - b(2)*w(1)
            v(2) = -b(2)*v(1)
            l(2) = (f - l(1)*h(1)) / g(2)
            h(2) = c - b(2)*f
    
            ! Step 2
            do i = 3,this%n-3
                b(i) = ( a - ( e/g(i-2) )*h(i-2) ) / g(i-1)
                h(i) = c - b(i)*f
                g(i) = d - ( e/g(i-2) )*f - b(i)*h(i-1)
            end do
    
            ! Step 3
            b(this%n-2) = ( a - ( e/g(this%n-4) )*h(this%n-4) ) / g(this%n-3)
            g(this%n-2) = d - ( e/g(this%n-4) )*f - b(this%n-2)*h(this%n-3)
    
            ! Step 4
            do i = 3,this%n-4
                k(i) = -( k(i-2)*f + k(i-1)*h(i-1) )/g(i)
                v(i) = -( e/g(i-2) )*v(i-2) - b(i)*v(i-1)
            end do
    
            ! Step 5
            k(this%n-3) = ( e - k(this%n-5)*f - k(this%n-4)*h(this%n-4) ) / g(this%n-3)
            k(this%n-2) = ( a - k(this%n-4)*f - k(this%n-3)*h(this%n-3) ) / g(this%n-2)
            v(this%n-3) = f - ( e/g(this%n-5) )*v(this%n-5) - b(this%n-3)*v(this%n-4)
            v(this%n-2) = c - ( e/g(this%n-4) )*v(this%n-4) - b(this%n-2)*v(this%n-3)
            g(this%n-1) = d - SUM( k(1:this%n-2)*v(1:this%n-2) )
    
            ! Step 6
            do i = 3,this%n-3
                w(i) = -( e/g(i-2) )*w(i-2) - b(i)*w(i-1)
                l(i) = -( l(i-2)*f + l(i-1)*h(i-1) ) / g(i)
            end do
    
            ! Step 7
            w(this%n-2) = f - ( e/g(this%n-4) )*w(this%n-4) - b(this%n-2)*w(this%n-3)
            w(this%n-1) = c - SUM( k(1:this%n-2)*w(1:this%n-2) )
            l(this%n-2) = ( e - l(this%n-4)*f - l(this%n-3)*h(this%n-3) ) / g(this%n-2)
            l(this%n-1) = ( a - SUM( l(1:this%n-2)*v(1:this%n-2) ) ) / g(this%n-1)
            g(this%n)   = d - SUM( l(1:this%n-1)*w(1:this%n-1) )
    
            ! Set eg(i) = e/g(i-2)
            eg(3:this%n-2) = e/g(1:this%n-4)
    
            ! Set ff = f
            ff(1:this%n-4) = f
    
            ! Set g = 1/g
            g = 1._rkind/g
    
        end associate
    
    end subroutine

    subroutine ComputePenta(this)
        class (cf90), intent(inout) :: this
        integer             :: i
    
        associate (bt   => this%penta(:,1), b   => this%penta(:,2), d => this%penta(:,3),  &
                   a    => this%penta(:,4), at  => this%penta(:,5),                         &
                   e    => this%penta(:,6), obc => this%penta(:,7),                         &
                   f    => this%penta(:,8), g   => this%penta(:,9),                         &
                   eobc => this%penta(:,10)                                                  )
    
            at = beta90 
            bt = beta90
            a  = alpha90
            b  = alpha90
            d  = one

            ! BC at 1
            bt(1) = zero
            b (1) = zero
            d (1) = one
            a (1) = zero 
            at(1) = zero

            bt(2) = zero
            b (2) = b2_alpha90
            d (2) = one
            a (2) = b2_alpha90
            at(2) = zero

            bt(3) = b3_beta90
            b (3) = b3_alpha90
            d (3) = one
            a (3) = b3_alpha90
            at(3) = b3_beta90

            bt(4) = b4_beta90
            b (4) = b4_alpha90
            d (4) = one
            a (4) = b4_alpha90
            at(4) = b4_beta90
            
            ! BC at n    
            bt(this%n  ) = zero
            b (this%n  ) = zero
            d (this%n  ) = one
            a (this%n  ) = zero 
            at(this%n  ) = zero

            bt(this%n-1) = zero
            b (this%n-1) = b2_alpha90
            d (this%n-1) = one
            a (this%n-1) = b2_alpha90
            at(this%n-1) = zero

            bt(this%n-2) = b3_beta90
            b (this%n-2) = b3_alpha90
            d (this%n-2) = one
            a (this%n-2) = b3_alpha90
            at(this%n-2) = b3_beta90

            bt(this%n-3) = b4_beta90
            b (this%n-3) = b4_alpha90
            d (this%n-3) = one
            a (this%n-3) = b4_alpha90
            at(this%n-3) = b4_beta90
            
            ! Step 1
            obc(1) = one/d(1)

            ! Step 2
            obc(2) = one/(d(2) - b(2)*a(1)*obc(1))

            ! Step 3
            e(1) = a(1)
            f(2) = b(2)*obc(1)
            
            do i = 3,this%n
                g(i) = bt(i)*obc(i-2)
                e(i-1) = a(i-1) - f(i-1)*at(i-2)
                f(i) = (b(i) - g(i)*e(i-2))*obc(i-1)
                obc(i) = one/(d(i) - f(i)*e(i-1) - g(i)*at(i-2))
            end do 

            eobc = e*obc
        end associate  
           
    end subroutine
    

    subroutine SolveXLU(this,y,n2,n3)
    
        class( cf90 ), intent(in) :: this
        integer, intent(in) :: n2,n3
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y  ! Take in RHS and put solution into it
        integer :: i,j,k
        real(rkind) :: sum1, sum2
 
        
        do k=1,n3
            do j=1,n2
                ! Step 8 ( update y instead of creating z )
                y(2,j,k) = y(2,j,k) - this%LU(2,1)*y(1,j,k) 
                sum1 = this%LU(1,3)*y(1,j,k) + this%LU(2,3)*y(2,j,k)
                sum2 = this%LU(1,4)*y(1,j,k) + this%LU(2,4)*y(2,j,k)

                ! Step 9
                do i = 3,this%n-2
                    y(i,j,k) = y(i,j,k) - this%LU(i,1)*y(i-1,j,k) - this%LU(i,2)*y(i-2,j,k)
                    sum1 = sum1 + this%LU(i,3)*y(i,j,k)
                    sum2 = sum2 + this%LU(i,4)*y(i,j,k)
                end do
    
                ! Step 10
                y(this%n-1,j,k) = y(this%n-1,j,k) - sum1
                y(this%n,j,k)   = ( y(this%n,j,k)   - sum2 - this%LU(this%n-1,4)*y(this%n-1,j,k) ) * this%LU(this%n,5)
    
                ! Step 11
                y(this%n-1,j,k) = ( y(this%n-1,j,k) - this%LU(this%n-1,9)*y(this%n,j,k) ) * this%LU(this%n-1,5)
                y(this%n-2,j,k) = ( y(this%n-2,j,k) - this%LU(this%n-2,8)*y(this%n-1,j,k) - this%LU(this%n-2,9)*y(this%n,j,k) ) * this%LU(this%n-2,5)
                y(this%n-3,j,k) = ( y(this%n-3,j,k) - this%LU(this%n-3,6)*y(this%n-2,j,k) - this%LU(this%n-3,8)*y(this%n-1,j,k) - this%LU(this%n-3,9)*y(this%n,j,k) ) * this%LU(this%n-3,5)
                do i = this%n-4,1,-1
                    y(i,j,k) = ( y(i,j,k) - this%LU(i,6)*y(i+1,j,k) - this%LU(i,7)*y(i+2,j,k) - this%LU(i,8)*y(this%n-1,j,k) - this%LU(i,9)*y(this%n,j,k) ) * this%LU(i,5)
                end do
            end do
        end do
    
    end subroutine
    
    subroutine SolveYLU(this,y,n1,n3)
    
        class( cf90 ), intent(in) :: this
        integer, intent(in) :: n1,n3
        real(rkind), dimension(n1,this%n,n3), intent(inout) :: y  ! Take in RHS and put solution into it
        integer :: j,k
        real(rkind), dimension(n1) :: sum1, sum2
 
        
        do k=1,n3
            ! Step 8 ( update y instead of creating z )
            y(:,2,k) = y(:,2,k) - this%LU(2,1)*y(:,1,k) 
            sum1 = this%LU(1,3)*y(:,1,k) + this%LU(2,3)*y(:,2,k)
            sum2 = this%LU(1,4)*y(:,1,k) + this%LU(2,4)*y(:,2,k)

            ! Step 9
            do j = 3,this%n-2
                y(:,j,k) = y(:,j,k) - this%LU(j,1)*y(:,j-1,k) - this%LU(j,2)*y(:,j-2,k)
                sum1 = sum1 + this%LU(j,3)*y(:,j,k)
                sum2 = sum2 + this%LU(j,4)*y(:,j,k)
            end do
    
            ! Step 10
            y(:,this%n-1,k) = y(:,this%n-1,k) - sum1
            y(:,this%n,k)   = ( y(:,this%n,k)   - sum2 - this%LU(this%n-1,4)*y(:,this%n-1,k) ) * this%LU(this%n,5)
    
            ! Step 11
            y(:,this%n-1,k) = ( y(:,this%n-1,k) - this%LU(this%n-1,9)*y(:,this%n,k) ) * this%LU(this%n-1,5)
            y(:,this%n-2,k) = ( y(:,this%n-2,k) - this%LU(this%n-2,8)*y(:,this%n-1,k) - this%LU(this%n-2,9)*y(:,this%n,k) ) * this%LU(this%n-2,5)
            y(:,this%n-3,k) = ( y(:,this%n-3,k) - this%LU(this%n-3,6)*y(:,this%n-2,k) - this%LU(this%n-3,8)*y(:,this%n-1,k) - this%LU(this%n-3,9)*y(:,this%n,k) ) * this%LU(this%n-3,5)
            do j = this%n-4,1,-1
                y(:,j,k) = ( y(:,j,k) - this%LU(j,6)*y(:,j+1,k) - this%LU(j,7)*y(:,j+2,k) - this%LU(j,8)*y(:,this%n-1,k) - this%LU(j,9)*y(:,this%n,k) ) * this%LU(j,5)
            end do
        end do
    
    end subroutine
    
    subroutine SolveZLU(this,y,n1,n2)
    
        class( cf90 ), intent(in) :: this
        integer, intent(in) :: n1,n2
        real(rkind), dimension(n1,n2,this%n), intent(inout) :: y  ! Take in RHS and put solution into it
        integer :: k
        real(rkind), dimension(n1,n2) :: sum1, sum2
 
        
        ! Step 8 ( update y instead of creating z )
        y(:,:,2) = y(:,:,2) - this%LU(2,1)*y(:,:,1) 
        sum1 = this%LU(1,3)*y(:,:,1) + this%LU(2,3)*y(:,:,2)
        sum2 = this%LU(1,4)*y(:,:,1) + this%LU(2,4)*y(:,:,2)

        ! Step 9
        do k = 3,this%n-2
            y(:,:,k) = y(:,:,k) - this%LU(k,1)*y(:,:,k-1) - this%LU(k,2)*y(:,:,k-2)
            sum1 = sum1 + this%LU(k,3)*y(:,:,k)
            sum2 = sum2 + this%LU(k,4)*y(:,:,k)
        end do
    
        ! Step 10
        y(:,:,this%n-1) = y(:,:,this%n-1) - sum1
        y(:,:,this%n)   = ( y(:,:,this%n)   - sum2 - this%LU(this%n-1,4)*y(:,:,this%n-1) ) * this%LU(this%n,5)
    
        ! Step 11
        y(:,:,this%n-1) = ( y(:,:,this%n-1) - this%LU(this%n-1,9)*y(:,:,this%n) ) * this%LU(this%n-1,5)
        y(:,:,this%n-2) = ( y(:,:,this%n-2) - this%LU(this%n-2,8)*y(:,:,this%n-1) - this%LU(this%n-2,9)*y(:,:,this%n) ) * this%LU(this%n-2,5)
        y(:,:,this%n-3) = ( y(:,:,this%n-3) - this%LU(this%n-3,6)*y(:,:,this%n-2) - this%LU(this%n-3,8)*y(:,:,this%n-1) - this%LU(this%n-3,9)*y(:,:,this%n) ) * this%LU(this%n-3,5)
        do k = this%n-4,1,-1
            y(:,:,k) = ( y(:,:,k) - this%LU(k,6)*y(:,:,k+1) - this%LU(k,7)*y(:,:,k+2) - this%LU(k,8)*y(:,:,this%n-1) - this%LU(k,9)*y(:,:,this%n) ) * this%LU(k,5)
        end do
    
    end subroutine
    
    subroutine SolveXPenta(this,y,n2,n3)

        class( cf90 ), intent(in) :: this
        integer, intent(in) :: n2,n3
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y
        integer :: i, j, k

        do k = 1,n3
            do j = 1,n2
                ! Step 1
                y(2,j,k) = y(2,j,k) - this%penta(2,8)*y(1,j,k)
                do i = 3,this%n
                    y(i,j,k) = y(i,j,k) - this%penta(i,9)*y(i-2,j,k) - this%penta(i,8)*y(i-1,j,k)
                end do 

                ! Step 2
                y(this%n,j,k) = y(this%n,j,k)*this%penta(this%n,7)
                
                y(this%n-1,j,k) = y(this%n-1,j,k)*this%penta(this%n-1,7) - this%penta(this%n-1,10)*y(this%n,j,k)
                do i = this%n-2,1,-1
                    y(i,j,k) = y(i,j,k)*this%penta(i,7) - y(i+2,j,k)*this%penta(i,5)*this%penta(i,7) - y(i+1,j,k)*this%penta(i,10)
                end do 
            end do 
        end do 

    end subroutine

    subroutine SolveYPenta(this,y,n1,n3)

        class( cf90 ), intent(in) :: this
        integer, intent(in) :: n1,n3
        real(rkind), dimension(n1,this%n,n3), intent(inout) :: y
        integer :: j, k

        do k = 1,n3
            ! Step 1
            y(:,2,k) = y(:,2,k) - this%penta(2,8)*y(:,1,k)
            do j = 3,this%n
                y(:,j,k) = y(:,j,k) - this%penta(j,9)*y(:,j-2,k) - this%penta(j,8)*y(:,j-1,k)
            end do 

            ! Step 2
            y(:,this%n,k) = y(:,this%n,k)*this%penta(this%n,7)
            
            y(:,this%n-1,k) = y(:,this%n-1,k)*this%penta(this%n-1,7) - this%penta(this%n-1,10)*y(:,this%n,k)
            do j = this%n-2,1,-1
                y(:,j,k) = y(:,j,k)*this%penta(j,7) - y(:,j+2,k)*this%penta(j,5)*this%penta(j,7) - y(:,j+1,k)*this%penta(j,10)
            end do 
        end do 

    end subroutine

    subroutine SolveZPenta(this,y,n1,n2)

        class( cf90 ), intent(in) :: this
        integer, intent(in) :: n1,n2
        real(rkind), dimension(n1,n2,this%n), intent(inout) :: y
        integer :: k

        ! Step 1
        y(:,:,2) = y(:,:,2) - this%penta(2,8)*y(:,:,1)
        do k = 3,this%n
            y(:,:,k) = y(:,:,k) - this%penta(k,9)*y(:,:,k-2) - this%penta(k,8)*y(:,:,k-1)
        end do 

        ! Step 2
        y(:,:,this%n) = y(:,:,this%n)*this%penta(this%n,7)
        
        y(:,:,this%n-1) = y(:,:,this%n-1)*this%penta(this%n-1,7) - this%penta(this%n-1,10)*y(:,:,this%n)
        do k = this%n-2,1,-1
            y(:,:,k) = y(:,:,k)*this%penta(k,7) - y(:,:,k+2)*this%penta(k,5)*this%penta(k,7) - y(:,:,k+1)*this%penta(k,10)
        end do 

    end subroutine

    pure subroutine ComputeXRHS(this, f, RHS, n2, n3)
    
        class( cf90 ), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        integer :: j,k


    
        select case (this%periodic)
        case (.TRUE.)
            do k=1,n3
                do j=1,n2
                    RHS(         1,j,k) = a90 * ( f(         1,j,k) )                     &
                                        + b90 * ( f(         2,j,k) + f(    this%n,j,k) ) &
                                        + c90 * ( f(         3,j,k) + f(  this%n-1,j,k) ) &
                                        + d90 * ( f(         4,j,k) + f(  this%n-2,j,k) ) &
                                        + e90 * ( f(         5,j,k) + f(  this%n-3,j,k) )
                    RHS(         2,j,k) = a90 * ( f(         2,j,k) )                     &
                                        + b90 * ( f(         3,j,k) + f(         1,j,k) ) &
                                        + c90 * ( f(         4,j,k) + f(    this%n,j,k) ) &
                                        + d90 * ( f(         5,j,k) + f(  this%n-1,j,k) ) &
                                        + e90 * ( f(         6,j,k) + f(  this%n-2,j,k) )
                    RHS(         3,j,k) = a90 * ( f(         3,j,k) )                     &
                                        + b90 * ( f(         4,j,k) + f(         2,j,k) ) &
                                        + c90 * ( f(         5,j,k) + f(         1,j,k) ) &
                                        + d90 * ( f(         6,j,k) + f(    this%n,j,k) ) &
                                        + e90 * ( f(         7,j,k) + f(  this%n-1,j,k) )
                    RHS(         4,j,k) = a90 * ( f(         4,j,k) )                     &
                                        + b90 * ( f(         5,j,k) + f(         3,j,k) ) &
                                        + c90 * ( f(         6,j,k) + f(         2,j,k) ) &
                                        + d90 * ( f(         7,j,k) + f(         1,j,k) ) &
                                        + e90 * ( f(         8,j,k) + f(    this%n,j,k) )
                    RHS(5:this%n-4,j,k) = a90 * ( f(5:this%n-4,j,k) )                     &
                                        + b90 * ( f(6:this%n-3,j,k) + f(4:this%n-5,j,k) ) &
                                        + c90 * ( f(7:this%n-2,j,k) + f(3:this%n-6,j,k) ) &
                                        + d90 * ( f(8:this%n-1,j,k) + f(2:this%n-7,j,k) ) &
                                        + e90 * ( f(9:this%n  ,j,k) + f(1:this%n-8,j,k) )
                    RHS(  this%n-3,j,k) = a90 * ( f(  this%n-3,j,k) )                     &
                                        + b90 * ( f(  this%n-2,j,k) + f(  this%n-4,j,k) ) &
                                        + c90 * ( f(  this%n-1,j,k) + f(  this%n-5,j,k) ) &
                                        + d90 * ( f(    this%n,j,k) + f(  this%n-6,j,k) ) &
                                        + e90 * ( f(         1,j,k) + f(  this%n-7,j,k) )
                    RHS(  this%n-2,j,k) = a90 * ( f(  this%n-2,j,k) )                     &
                                        + b90 * ( f(  this%n-1,j,k) + f(  this%n-3,j,k) ) &
                                        + c90 * ( f(    this%n,j,k) + f(  this%n-4,j,k) ) &
                                        + d90 * ( f(         1,j,k) + f(  this%n-5,j,k) ) &
                                        + e90 * ( f(         2,j,k) + f(  this%n-6,j,k) )
                    RHS(  this%n-1,j,k) = a90 * ( f(  this%n-1,j,k) )                     &
                                        + b90 * ( f(    this%n,j,k) + f(  this%n-2,j,k) ) &
                                        + c90 * ( f(         1,j,k) + f(  this%n-3,j,k) ) &
                                        + d90 * ( f(         2,j,k) + f(  this%n-4,j,k) ) &
                                        + e90 * ( f(         3,j,k) + f(  this%n-5,j,k) )
                    RHS(    this%n,j,k) = a90 * ( f(    this%n,j,k) )                     &
                                        + b90 * ( f(         1,j,k) + f(  this%n-1,j,k) ) &
                                        + c90 * ( f(         2,j,k) + f(  this%n-2,j,k) ) &
                                        + d90 * ( f(         3,j,k) + f(  this%n-3,j,k) ) &
                                        + e90 * ( f(         4,j,k) + f(  this%n-4,j,k) )
                end do
            end do

        case (.FALSE.)

            do k = 1,n3
                do j = 1,n2
                    
                    RHS(         1,j,k) =    one * ( f(         1,j,k) )                     

                    RHS(         2,j,k) = b2_a90 * ( f(         2,j,k) )                     &
                                        + b2_b90 * ( f(         3,j,k) + f(         1,j,k) ) 
                    
                    RHS(         3,j,k) = b3_a90 * ( f(         3,j,k) )                     &
                                        + b3_b90 * ( f(         4,j,k) + f(         2,j,k) ) &
                                        + b3_c90 * ( f(         5,j,k) + f(         1,j,k) )

                    RHS(         4,j,k) = b4_a90 * ( f(         4,j,k) )                     &
                                        + b4_b90 * ( f(         5,j,k) + f(         3,j,k) ) &
                                        + b4_c90 * ( f(         6,j,k) + f(         2,j,k) ) &
                                        + b4_d90 * ( f(         7,j,k) + f(         1,j,k) ) 

                    RHS(5:this%n-4,j,k) =    a90 * ( f(5:this%n-4,j,k) )                     &
                                        +    b90 * ( f(6:this%n-3,j,k) + f(4:this%n-5,j,k) ) &
                                        +    c90 * ( f(7:this%n-2,j,k) + f(3:this%n-6,j,k) ) &
                                        +    d90 * ( f(8:this%n-1,j,k) + f(2:this%n-7,j,k) ) &
                                        +    e90 * ( f(9:this%n  ,j,k) + f(1:this%n-8,j,k) )

                    RHS(  this%n-3,j,k) = b4_a90 * ( f(  this%n-3,j,k) )                     &
                                        + b4_b90 * ( f(  this%n-2,j,k) + f(  this%n-4,j,k) ) &
                                        + b4_c90 * ( f(  this%n-1,j,k) + f(  this%n-5,j,k) ) &
                                        + b4_d90 * ( f(    this%n,j,k) + f(  this%n-6,j,k) ) 

                    RHS(  this%n-2,j,k) = b3_a90 * ( f(  this%n-2,j,k) )                     &
                                        + b3_b90 * ( f(  this%n-1,j,k) + f(  this%n-3,j,k) ) &
                                        + b3_c90 * ( f(    this%n,j,k) + f(  this%n-4,j,k) ) 

                    RHS(  this%n-1,j,k) = b2_a90 * ( f(  this%n-1,j,k) )                     &
                                        + b2_b90 * ( f(    this%n,j,k) + f(  this%n-2,j,k) ) 

                    RHS(    this%n,j,k) =    one * ( f(    this%n,j,k) )                     
                
               end do 
            end do 
        end select
    
    end subroutine
    
    pure subroutine ComputeYRHS(this, f, RHS, n1, n3) 
    
        class( cf90 ), intent(in) :: this
        integer, intent(in) :: n1, n3
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: RHS
        integer :: k

    
        select case (this%periodic)
        case (.TRUE.)
            do k=1,n3
                RHS(:,         1,k) = a90 * ( f(:,         1,k) )                     &
                                    + b90 * ( f(:,         2,k) + f(:,    this%n,k) ) &
                                    + c90 * ( f(:,         3,k) + f(:,  this%n-1,k) ) &
                                    + d90 * ( f(:,         4,k) + f(:,  this%n-2,k) ) &
                                    + e90 * ( f(:,         5,k) + f(:,  this%n-3,k) )
                RHS(:,         2,k) = a90 * ( f(:,         2,k) )                     &
                                    + b90 * ( f(:,         3,k) + f(:,         1,k) ) &
                                    + c90 * ( f(:,         4,k) + f(:,    this%n,k) ) &
                                    + d90 * ( f(:,         5,k) + f(:,  this%n-1,k) ) &
                                    + e90 * ( f(:,         6,k) + f(:,  this%n-2,k) )
                RHS(:,         3,k) = a90 * ( f(:,         3,k) )                     &
                                    + b90 * ( f(:,         4,k) + f(:,         2,k) ) &
                                    + c90 * ( f(:,         5,k) + f(:,         1,k) ) &
                                    + d90 * ( f(:,         6,k) + f(:,    this%n,k) ) &
                                    + e90 * ( f(:,         7,k) + f(:,  this%n-1,k) )
                RHS(:,         4,k) = a90 * ( f(:,         4,k) )                     &
                                    + b90 * ( f(:,         5,k) + f(:,         3,k) ) &
                                    + c90 * ( f(:,         6,k) + f(:,         2,k) ) &
                                    + d90 * ( f(:,         7,k) + f(:,         1,k) ) &
                                    + e90 * ( f(:,         8,k) + f(:,    this%n,k) )
                RHS(:,5:this%n-4,k) = a90 * ( f(:,5:this%n-4,k) )                     &
                                    + b90 * ( f(:,6:this%n-3,k) + f(:,4:this%n-5,k) ) &
                                    + c90 * ( f(:,7:this%n-2,k) + f(:,3:this%n-6,k) ) &
                                    + d90 * ( f(:,8:this%n-1,k) + f(:,2:this%n-7,k) ) &
                                    + e90 * ( f(:,9:this%n  ,k) + f(:,1:this%n-8,k) )
                RHS(:,  this%n-3,k) = a90 * ( f(:,  this%n-3,k) )                     &
                                    + b90 * ( f(:,  this%n-2,k) + f(:,  this%n-4,k) ) &
                                    + c90 * ( f(:,  this%n-1,k) + f(:,  this%n-5,k) ) &
                                    + d90 * ( f(:,    this%n,k) + f(:,  this%n-6,k) ) &
                                    + e90 * ( f(:,         1,k) + f(:,  this%n-7,k) )
                RHS(:,  this%n-2,k) = a90 * ( f(:,  this%n-2,k) )                     &
                                    + b90 * ( f(:,  this%n-1,k) + f(:,  this%n-3,k) ) &
                                    + c90 * ( f(:,    this%n,k) + f(:,  this%n-4,k) ) &
                                    + d90 * ( f(:,         1,k) + f(:,  this%n-5,k) ) &
                                    + e90 * ( f(:,         2,k) + f(:,  this%n-6,k) )
                RHS(:,  this%n-1,k) = a90 * ( f(:,  this%n-1,k) )                     &
                                    + b90 * ( f(:,    this%n,k) + f(:,  this%n-2,k) ) &
                                    + c90 * ( f(:,         1,k) + f(:,  this%n-3,k) ) &
                                    + d90 * ( f(:,         2,k) + f(:,  this%n-4,k) ) &
                                    + e90 * ( f(:,         3,k) + f(:,  this%n-5,k) )
                RHS(:,    this%n,k) = a90 * ( f(:,    this%n,k) )                     &
                                    + b90 * ( f(:,         1,k) + f(:,  this%n-1,k) ) &
                                    + c90 * ( f(:,         2,k) + f(:,  this%n-2,k) ) &
                                    + d90 * ( f(:,         3,k) + f(:,  this%n-3,k) ) &
                                    + e90 * ( f(:,         4,k) + f(:,  this%n-4,k) )
            end do
        case (.FALSE.)

            do k = 1,n3
                RHS(:,         1,k) =    one * ( f(:,         1,k) )                     

                RHS(:,         2,k) = b2_a90 * ( f(:,         2,k) )                     &
                                    + b2_b90 * ( f(:,         3,k) + f(:,         1,k) ) 
                
                RHS(:,         3,k) = b3_a90 * ( f(:,         3,k) )                     &
                                    + b3_b90 * ( f(:,         4,k) + f(:,         2,k) ) &
                                    + b3_c90 * ( f(:,         5,k) + f(:,         1,k) )

                RHS(:,         4,k) = b4_a90 * ( f(:,         4,k) )                     &
                                    + b4_b90 * ( f(:,         5,k) + f(:,         3,k) ) &
                                    + b4_c90 * ( f(:,         6,k) + f(:,         2,k) ) &
                                    + b4_d90 * ( f(:,         7,k) + f(:,         1,k) ) 

                RHS(:,5:this%n-4,k) =    a90 * ( f(:,5:this%n-4,k) )                     &
                                    +    b90 * ( f(:,6:this%n-3,k) + f(:,4:this%n-5,k) ) &
                                    +    c90 * ( f(:,7:this%n-2,k) + f(:,3:this%n-6,k) ) &
                                    +    d90 * ( f(:,8:this%n-1,k) + f(:,2:this%n-7,k) ) &
                                    +    e90 * ( f(:,9:this%n  ,k) + f(:,1:this%n-8,k) )

                RHS(:,  this%n-3,k) = b4_a90 * ( f(:,  this%n-3,k) )                     &
                                    + b4_b90 * ( f(:,  this%n-2,k) + f(:,  this%n-4,k) ) &
                                    + b4_c90 * ( f(:,  this%n-1,k) + f(:,  this%n-5,k) ) &
                                    + b4_d90 * ( f(:,    this%n,k) + f(:,  this%n-6,k) ) 

                RHS(:,  this%n-2,k) = b3_a90 * ( f(:,  this%n-2,k) )                     &
                                    + b3_b90 * ( f(:,  this%n-1,k) + f(:,  this%n-3,k) ) &
                                    + b3_c90 * ( f(:,    this%n,k) + f(:,  this%n-4,k) ) 

                RHS(:,  this%n-1,k) = b2_a90 * ( f(:,  this%n-1,k) )                     &
                                    + b2_b90 * ( f(:,    this%n,k) + f(:,  this%n-2,k) ) 

                RHS(:,    this%n,k) =    one * ( f(:,    this%n,k) )                     
                
            end do 
        end select
    
    end subroutine

    pure subroutine ComputeZRHS(this, f, RHS, n1, n2)
    
        class( cf90 ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in) :: f
        real(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
    
        select case (this%periodic)
        case (.TRUE.)
                RHS(:,:,         1) = a90 * ( f(:,:,         1) )                     &
                                    + b90 * ( f(:,:,         2) + f(:,:,    this%n) ) &
                                    + c90 * ( f(:,:,         3) + f(:,:,  this%n-1) ) &
                                    + d90 * ( f(:,:,         4) + f(:,:,  this%n-2) ) &
                                    + e90 * ( f(:,:,         5) + f(:,:,  this%n-3) )
                RHS(:,:,         2) = a90 * ( f(:,:,         2) )                     &
                                    + b90 * ( f(:,:,         3) + f(:,:,         1) ) &
                                    + c90 * ( f(:,:,         4) + f(:,:,    this%n) ) &
                                    + d90 * ( f(:,:,         5) + f(:,:,  this%n-1) ) &
                                    + e90 * ( f(:,:,         6) + f(:,:,  this%n-2) )
                RHS(:,:,         3) = a90 * ( f(:,:,         3) )                     &
                                    + b90 * ( f(:,:,         4) + f(:,:,         2) ) &
                                    + c90 * ( f(:,:,         5) + f(:,:,         1) ) &
                                    + d90 * ( f(:,:,         6) + f(:,:,    this%n) ) &
                                    + e90 * ( f(:,:,         7) + f(:,:,  this%n-1) )
                RHS(:,:,         4) = a90 * ( f(:,:,         4) )                     &
                                    + b90 * ( f(:,:,         5) + f(:,:,         3) ) &
                                    + c90 * ( f(:,:,         6) + f(:,:,         2) ) &
                                    + d90 * ( f(:,:,         7) + f(:,:,         1) ) &
                                    + e90 * ( f(:,:,         8) + f(:,:,    this%n) )
                RHS(:,:,5:this%n-4) = a90 * ( f(:,:,5:this%n-4) )                     &
                                    + b90 * ( f(:,:,6:this%n-3) + f(:,:,4:this%n-5) ) &
                                    + c90 * ( f(:,:,7:this%n-2) + f(:,:,3:this%n-6) ) &
                                    + d90 * ( f(:,:,8:this%n-1) + f(:,:,2:this%n-7) ) &
                                    + e90 * ( f(:,:,9:this%n  ) + f(:,:,1:this%n-8) )
                RHS(:,:,  this%n-3) = a90 * ( f(:,:,  this%n-3) )                     &
                                    + b90 * ( f(:,:,  this%n-2) + f(:,:,  this%n-4) ) &
                                    + c90 * ( f(:,:,  this%n-1) + f(:,:,  this%n-5) ) &
                                    + d90 * ( f(:,:,    this%n) + f(:,:,  this%n-6) ) &
                                    + e90 * ( f(:,:,         1) + f(:,:,  this%n-7) )
                RHS(:,:,  this%n-2) = a90 * ( f(:,:,  this%n-2) )                     &
                                    + b90 * ( f(:,:,  this%n-1) + f(:,:,  this%n-3) ) &
                                    + c90 * ( f(:,:,    this%n) + f(:,:,  this%n-4) ) &
                                    + d90 * ( f(:,:,         1) + f(:,:,  this%n-5) ) &
                                    + e90 * ( f(:,:,         2) + f(:,:,  this%n-6) )
                RHS(:,:,  this%n-1) = a90 * ( f(:,:,  this%n-1) )                     &
                                    + b90 * ( f(:,:,    this%n) + f(:,:,  this%n-2) ) &
                                    + c90 * ( f(:,:,         1) + f(:,:,  this%n-3) ) &
                                    + d90 * ( f(:,:,         2) + f(:,:,  this%n-4) ) &
                                    + e90 * ( f(:,:,         3) + f(:,:,  this%n-5) )
                RHS(:,:,    this%n) = a90 * ( f(:,:,    this%n) )                     &
                                    + b90 * ( f(:,:,         1) + f(:,:,  this%n-1) ) &
                                    + c90 * ( f(:,:,         2) + f(:,:,  this%n-2) ) &
                                    + d90 * ( f(:,:,         3) + f(:,:,  this%n-3) ) &
                                    + e90 * ( f(:,:,         4) + f(:,:,  this%n-4) )
        case (.FALSE.)
                    
            RHS(:,:,         1) =    one * ( f(:,:,         1) )                     

            RHS(:,:,         2) = b2_a90 * ( f(:,:,         2) )                     &
                                + b2_b90 * ( f(:,:,         3) + f(:,:,         1) ) 
            
            RHS(:,:,         3) = b3_a90 * ( f(:,:,         3) )                     &
                                + b3_b90 * ( f(:,:,         4) + f(:,:,         2) ) &
                                + b3_c90 * ( f(:,:,         5) + f(:,:,         1) )

            RHS(:,:,         4) = b4_a90 * ( f(:,:,         4) )                     &
                                + b4_b90 * ( f(:,:,         5) + f(:,:,         3) ) &
                                + b4_c90 * ( f(:,:,         6) + f(:,:,         2) ) &
                                + b4_d90 * ( f(:,:,         7) + f(:,:,         1) ) 

            RHS(:,:,5:this%n-4) =    a90 * ( f(:,:,5:this%n-4) )                     &
                                +    b90 * ( f(:,:,6:this%n-3) + f(:,:,4:this%n-5) ) &
                                +    c90 * ( f(:,:,7:this%n-2) + f(:,:,3:this%n-6) ) &
                                +    d90 * ( f(:,:,8:this%n-1) + f(:,:,2:this%n-7) ) &
                                +    e90 * ( f(:,:,9:this%n  ) + f(:,:,1:this%n-8) )

            RHS(:,:,  this%n-3) = b4_a90 * ( f(:,:,  this%n-3) )                     &
                                + b4_b90 * ( f(:,:,  this%n-2) + f(:,:,  this%n-4) ) &
                                + b4_c90 * ( f(:,:,  this%n-1) + f(:,:,  this%n-5) ) &
                                + b4_d90 * ( f(:,:,    this%n) + f(:,:,  this%n-6) ) 

            RHS(:,:,  this%n-2) = b3_a90 * ( f(:,:,  this%n-2) )                     &
                                + b3_b90 * ( f(:,:,  this%n-1) + f(:,:,  this%n-3) ) &
                                + b3_c90 * ( f(:,:,    this%n) + f(:,:,  this%n-4) ) 

            RHS(:,:,  this%n-1) = b2_a90 * ( f(:,:,  this%n-1) )                     &
                                + b2_b90 * ( f(:,:,    this%n) + f(:,:,  this%n-2) ) 

            RHS(:,:,    this%n) =    one * ( f(:,:,    this%n) )                     

        end select
    
    end subroutine
    
    subroutine filter1(this, f, df, na, nb)
        class( cf90 ), intent(in) :: this
        integer, intent(in) :: na, nb
        real(rkind), dimension(this%n,na,nb), intent(in) :: f
        real(rkind), dimension(this%n,na,nb), intent(out) :: df

        if(this%n == 1) then
            df = f
            return
        end if
        
        call this%ComputeXRHS(f, df, na, nb)

        select case (this%periodic)
        case(.TRUE.)
            call this%SolveXLU(df, na, nb)
        case(.FALSE.)
            call this%SolveXPenta(df, na, nb)
        end select
    
    end subroutine

    subroutine filter2(this, f, df, na, nb)
        class( cf90 ), intent(in) :: this
        integer, intent(in) :: na, nb
        real(rkind), dimension(na,this%n,nb), intent(in) :: f
        real(rkind), dimension(na,this%n,nb), intent(out) :: df

        if(this%n == 1) then
            df = f
            return
        end if
        
        call this%ComputeYRHS(f, df, na, nb)

        select case (this%periodic)
        case(.TRUE.)
            call this%SolveYLU(df, na, nb)
        case(.FALSE.)
            call this%SolveYPenta(df, na, nb)
        end select
    
    end subroutine

    subroutine filter3(this, f, df, na, nb)
        class( cf90 ), intent(in) :: this
        integer, intent(in) :: na, nb
        real(rkind), dimension(na,nb,this%n), intent(in) :: f
        real(rkind), dimension(na,nb,this%n), intent(out) :: df

        if(this%n == 1) then
            df = f
            return
        end if
        
        call this%ComputeZRHS(f, df, na, nb)

        select case (this%periodic)
        case(.TRUE.)
            call this%SolveZLU(df, na, nb)
        case(.FALSE.)
            call this%SolveZPenta(df, na, nb)
        end select
    
    end subroutine

end module
