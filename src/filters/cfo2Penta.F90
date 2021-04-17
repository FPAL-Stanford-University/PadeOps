! Routines specific to 8th order Compact Filter
! Periodic LU based on  Neossi Nguetchue, Abelman (Appl. Math. & Comp. 2008)

module cfo2Pentastuff

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two,four
    use exits,           only: GracefulExit

    implicit none

    private
    public :: cfo2Penta, alphao2Penta, betao2Penta, ao2Penta, bo2Penta, co2Penta, do2Penta, eo2Penta
    
    ! second order penta diagonal filter with adjustable coefficients -- (d^2 sf(0)) / (d omega^2) = zero -- less dissipative than tridiagonal

    
    ! !set coefficients from alphao2Penta
    ! !interior filter
    ! alphao2Penta = alphao2Penta_
    ! betao2Penta = (0.5 - alphao2Penta) / 7.0
    ! ao2Penta = (4.0 + 6.0 * alphao2Penta) / 7.0
    ! bo2Penta = (2.0 + 3.0 * alphao2Penta) / 7.0

    ! !boundary point -- tridiagonal
    ! b2_alphao2Penta = alphao2Penta !maintain alpha
    ! b2_ao2Penta     = (two*alphao2Penta+one)/two
    ! b2_bo2Penta     = (two*alphao2Penta+one)/four


    real(rkind) :: alphao2Penta= zero ! Already divided by factor of 2 
    real(rkind) :: betao2Penta = zero ! Already divided by factor of 2 
    real(rkind) :: ao2Penta    = zero
    real(rkind) :: bo2Penta    = zero ! Already divided by factor of 2 
    real(rkind) :: co2Penta    = zero ! Already divided by factor of 2 
    real(rkind) :: do2Penta    = zero ! Already divided by factor of 2 
    real(rkind) :: eo2Penta    = zero ! Already divided by factor of 2 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic filter !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
    ! First point is just identity

    ! Second point
    real(rkind) :: b2_alphao2Penta= zero ! Already divided by factor of 2 
    real(rkind) :: b2_ao2Penta    = zero
    real(rkind) :: b2_bo2Penta    = zero ! Already divided by factor of 2  
    
    !other boundary points same as interior
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    type cfo2Penta
        
        private
        
        integer     :: n

        logical     :: periodic=.TRUE.

        real(rkind), allocatable, dimension(:,:) :: LU
        real(rkind), allocatable, dimension(:,:) :: penta_nn
        real(rkind), allocatable, dimension(:,:) :: penta_ns
        real(rkind), allocatable, dimension(:,:) :: penta_na
        real(rkind), allocatable, dimension(:,:) :: penta_sn
        real(rkind), allocatable, dimension(:,:) :: penta_ss
        real(rkind), allocatable, dimension(:,:) :: penta_sa
        real(rkind), allocatable, dimension(:,:) :: penta_an
        real(rkind), allocatable, dimension(:,:) :: penta_as
        real(rkind), allocatable, dimension(:,:) :: penta_aa

        contains

        procedure :: init
        procedure :: destroy

        procedure, private :: ComputeXRHS
        procedure, private :: ComputeYRHS
        procedure, private :: ComputeZRHS_REAL
        procedure, private :: ComputeZRHS_CMPLX
        
        procedure, private :: SolveXLU
        procedure, private :: SolveYLU
        procedure, private :: SolveZLU
        
        procedure, private :: ComputeLU
        procedure, private :: ComputePenta

        procedure, private :: SolveXPenta
        procedure, private :: SolveYPenta
        procedure, private :: SolveZPenta_REAL
        procedure, private :: SolveZPenta_CMPLX
        
        procedure :: filter1
        procedure :: filter2
        procedure, private :: filter3_REAL
        procedure, private :: filter3_CMPLX
        generic :: filter3 => filter3_REAL, filter3_CMPLX 
        
    end type



contains

    function init(this, n_, periodic_, alphao2Penta_) result(ierr)
   
        class( cfo2Penta ), intent(inout) :: this
        integer, intent(in) :: n_
        logical, intent(in) :: periodic_
        real(rkind), intent(in) :: alphao2Penta_
        integer :: ierr
        
        this%n = n_

        this%periodic = periodic_

        !set coefficients from alphao2Penta
        !interior filter
        alphao2Penta = alphao2Penta_
        betao2Penta = (0.5 - alphao2Penta) / 7.0
        ao2Penta = (4.0 + 6.0 * alphao2Penta) / 7.0
        bo2Penta = (2.0 + 3.0 * alphao2Penta) / 7.0

        !boundary point -- tridiagonal
        b2_alphao2Penta = alphao2Penta !maintain alpha
        b2_ao2Penta     = (two*alphao2Penta+one)/two
        b2_bo2Penta     = (two*alphao2Penta+one)/four


        if (periodic_) then 
            ! Allocate LU matrix
            if(allocated( this%LU )) deallocate( this%LU ); allocate( this%LU(n_,9) )
    
            ! Compute LU matrix
            if (n_ .GE. 10) then
                call this%ComputeLU(this%LU,betao2Penta,alphao2Penta,one,alphao2Penta,betao2Penta)
            else if (n_ == 1) then
                this%LU = one
            else
                ierr = 7
                return
            end if
    
        else 

            ! Allocate Penta matrix
            if(allocated( this%penta_nn )) deallocate( this%penta_nn ); allocate( this%penta_nn(n_,11) )
            if(allocated( this%penta_ns )) deallocate( this%penta_ns ); allocate( this%penta_ns(n_,11) )
            if(allocated( this%penta_na )) deallocate( this%penta_na ); allocate( this%penta_na(n_,11) )
            if(allocated( this%penta_sn )) deallocate( this%penta_sn ); allocate( this%penta_sn(n_,11) )
            if(allocated( this%penta_ss )) deallocate( this%penta_ss ); allocate( this%penta_ss(n_,11) )
            if(allocated( this%penta_sa )) deallocate( this%penta_sa ); allocate( this%penta_sa(n_,11) )
            if(allocated( this%penta_an )) deallocate( this%penta_an ); allocate( this%penta_an(n_,11) )
            if(allocated( this%penta_as )) deallocate( this%penta_as ); allocate( this%penta_as(n_,11) )
            if(allocated( this%penta_aa )) deallocate( this%penta_aa ); allocate( this%penta_aa(n_,11) )
  
            if (n_ .GE. 10) then             
                call this%ComputePenta(this%penta_nn, 0, 0)
                call this%ComputePenta(this%penta_ns, 0, 1)
                call this%ComputePenta(this%penta_na, 0,-1)
                call this%ComputePenta(this%penta_sn, 1, 0)
                call this%ComputePenta(this%penta_ss, 1, 1)
                call this%ComputePenta(this%penta_sa, 1,-1)
                call this%ComputePenta(this%penta_an,-1, 0)
                call this%ComputePenta(this%penta_as,-1, 1)
                call this%ComputePenta(this%penta_aa,-1,-1)
            else if (n_ .EQ. 1) then
                this%penta_nn = one
                this%penta_ns = one
                this%penta_na = one
                this%penta_sn = one
                this%penta_ss = one
                this%penta_sa = one
                this%penta_an = one
                this%penta_as = one
                this%penta_aa = one
            else
                ierr = 7
                return 
            end if 
            
        end if 

        ! If everything passes
        ierr = 0
    
    end function
    
    subroutine destroy(this)

        class( cfo2Penta ), intent(inout) :: this

        ! Dellocate LU matrix.
        if(allocated( this%LU )) deallocate( this%LU )
    
        ! Dellocate penta matrix.
        if(allocated( this%penta_nn )) deallocate( this%penta_nn )
        if(allocated( this%penta_ns )) deallocate( this%penta_ns )
        if(allocated( this%penta_na )) deallocate( this%penta_na )
        if(allocated( this%penta_sn )) deallocate( this%penta_sn )
        if(allocated( this%penta_ss )) deallocate( this%penta_ss )
        if(allocated( this%penta_sa )) deallocate( this%penta_sa )
        if(allocated( this%penta_an )) deallocate( this%penta_an )
        if(allocated( this%penta_as )) deallocate( this%penta_as )
        if(allocated( this%penta_aa )) deallocate( this%penta_aa )
    
    end subroutine

    
    subroutine ComputeLU(this,LU,e,a,d,c,f) 
    
        class( cfo2Penta ), intent(inout) :: this
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

    subroutine ComputePenta(this,penta,bc1,bcn)
        class (cfo2Penta), intent(inout) :: this
        real(rkind), dimension(this%n,11), intent(inout) :: penta
        integer, intent(in) :: bc1,bcn
        integer             :: i
    
        associate (bt   => penta(:,1), b   => penta(:,2), d => penta(:,3),  &
                   a    => penta(:,4), at  => penta(:,5),                   &
                   e    => penta(:,6), obc => penta(:,7),                   &
                   f    => penta(:,8), g   => penta(:,9),                   &
                   eobc => penta(:,10)                                      )
    
            at = betao2Penta 
            bt = betao2Penta
            a  = alphao2Penta
            b  = alphao2Penta
            d  = one

            ! BC at 1
            select case(bc1)
            case(0)
                bt(1) = zero
                b (1) = zero
                d (1) = one
                a (1) = zero 
                at(1) = zero

                bt(2) = zero
                b (2) = b2_alphao2Penta
                d (2) = one
                a (2) = b2_alphao2Penta
                at(2) = zero

                bt(3) = betao2Penta
                b (3) = alphao2Penta
                d (3) = one
                a (3) = alphao2Penta
                at(3) = betao2Penta

                bt(4) = betao2Penta
                b (4) = alphao2Penta
                d (4) = one
                a (4) = alphao2Penta
                at(4) = betao2Penta
            case(1)
                bt(1) = zero
                b (1) = zero
                d (1) = one
                a (1) = two*alphao2Penta 
                at(1) = two*betao2Penta

                bt(2) = zero
                b (2) = alphao2Penta
                d (2) = one + betao2Penta
                a (2) = alphao2Penta
                at(2) = betao2Penta
            case(-1)
                bt(1) = zero
                b (1) = zero
                d (1) = one
                a (1) = zero
                at(1) = zero

                bt(2) = zero
                b (2) = alphao2Penta
                d (2) = one - betao2Penta
                a (2) = alphao2Penta
                at(2) = betao2Penta
            end select
            
            ! BC at n    
            select case(bcn)
            case(0)
                bt(this%n  ) = zero
                b (this%n  ) = zero
                d (this%n  ) = one
                a (this%n  ) = zero 
                at(this%n  ) = zero

                bt(this%n-1) = zero
                b (this%n-1) = b2_alphao2Penta
                d (this%n-1) = one
                a (this%n-1) = b2_alphao2Penta
                at(this%n-1) = zero

                bt(this%n-2) = betao2Penta
                b (this%n-2) = alphao2Penta
                d (this%n-2) = one
                a (this%n-2) = alphao2Penta
                at(this%n-2) = betao2Penta

                bt(this%n-3) = betao2Penta
                b (this%n-3) = alphao2Penta
                d (this%n-3) = one
                a (this%n-3) = alphao2Penta
                at(this%n-3) = betao2Penta
            case(1)
                bt(this%n  ) = two*betao2Penta
                b (this%n  ) = two*alphao2Penta
                d (this%n  ) = one
                a (this%n  ) = zero 
                at(this%n  ) = zero

                bt(this%n-1) = betao2Penta
                b (this%n-1) = alphao2Penta
                d (this%n-1) = one + betao2Penta
                a (this%n-1) = alphao2Penta
                at(this%n-1) = zero
            case(-1)
                bt(this%n  ) = zero
                b (this%n  ) = zero
                d (this%n  ) = one
                a (this%n  ) = zero 
                at(this%n  ) = zero

                bt(this%n-1) = betao2Penta
                b (this%n-1) = alphao2Penta
                d (this%n-1) = one - betao2Penta
                a (this%n-1) = alphao2Penta
                at(this%n-1) = zero
            end select
            
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
    
        class( cfo2Penta ), intent(in) :: this
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
    
        class( cfo2Penta ), intent(in) :: this
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
    
        class( cfo2Penta ), intent(in) :: this
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
    
    subroutine SolveXPenta(this,penta,y,n2,n3)

        class( cfo2Penta ), intent(in) :: this
        real(rkind), dimension(this%n,11), intent(in) :: penta
        integer, intent(in) :: n2,n3
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y
        integer :: i, j, k

        do k = 1,n3
            do j = 1,n2
                ! Step 1
                y(2,j,k) = y(2,j,k) - penta(2,8)*y(1,j,k)
                do i = 3,this%n
                    y(i,j,k) = y(i,j,k) - penta(i,9)*y(i-2,j,k) - penta(i,8)*y(i-1,j,k)
                end do 

                ! Step 2
                y(this%n,j,k) = y(this%n,j,k)*penta(this%n,7)
                
                y(this%n-1,j,k) = y(this%n-1,j,k)*penta(this%n-1,7) - penta(this%n-1,10)*y(this%n,j,k)
                do i = this%n-2,1,-1
                    y(i,j,k) = y(i,j,k)*penta(i,7) - y(i+2,j,k)*penta(i,5)*penta(i,7) - y(i+1,j,k)*penta(i,10)
                end do 
            end do 
        end do 

    end subroutine

    subroutine SolveYPenta(this,penta,y,n1,n3)

        class( cfo2Penta ), intent(in) :: this
        real(rkind), dimension(this%n,11), intent(in) :: penta
        integer, intent(in) :: n1,n3
        real(rkind), dimension(n1,this%n,n3), intent(inout) :: y
        integer :: j, k

        do k = 1,n3
            ! Step 1
            y(:,2,k) = y(:,2,k) - penta(2,8)*y(:,1,k)
            do j = 3,this%n
                y(:,j,k) = y(:,j,k) - penta(j,9)*y(:,j-2,k) - penta(j,8)*y(:,j-1,k)
            end do 

            ! Step 2
            y(:,this%n,k) = y(:,this%n,k)*penta(this%n,7)
            
            y(:,this%n-1,k) = y(:,this%n-1,k)*penta(this%n-1,7) - penta(this%n-1,10)*y(:,this%n,k)
            do j = this%n-2,1,-1
                y(:,j,k) = y(:,j,k)*penta(j,7) - y(:,j+2,k)*penta(j,5)*penta(j,7) - y(:,j+1,k)*penta(j,10)
            end do 
        end do 

    end subroutine

    subroutine SolveZPenta_REAL(this,penta,y,n1,n2)

        class( cfo2Penta ), intent(in) :: this
        real(rkind), dimension(this%n,11), intent(in) :: penta
        integer, intent(in) :: n1,n2
        real(rkind), dimension(n1,n2,this%n), intent(inout) :: y
        integer :: k

#include "CFo2Penta_files/SolveZPenta_common.F90"        

    end subroutine

    subroutine SolveZPenta_CMPLX(this,penta,y,n1,n2)

        class( cfo2Penta ), intent(in) :: this
        real(rkind), dimension(this%n,11), intent(in) :: penta
        integer, intent(in) :: n1,n2
        complex(rkind), dimension(n1,n2,this%n), intent(inout) :: y
        integer :: k

#include "CFo2Penta_files/SolveZPenta_common.F90"        

    end subroutine

    pure subroutine ComputeXRHS(this, f, RHS, n2, n3, bc1, bcn)
    
        class( cfo2Penta ), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn
        integer :: j,k

        select case (this%periodic)
        case (.TRUE.)
            do k=1,n3
                do j=1,n2
                    RHS(         1,j,k) = ao2Penta * ( f(         1,j,k) )                     &
                                        + bo2Penta * ( f(         2,j,k) + f(    this%n,j,k) ) &
                                        + co2Penta * ( f(         3,j,k) + f(  this%n-1,j,k) ) &
                                        + do2Penta * ( f(         4,j,k) + f(  this%n-2,j,k) ) &
                                        + eo2Penta * ( f(         5,j,k) + f(  this%n-3,j,k) )
                    RHS(         2,j,k) = ao2Penta * ( f(         2,j,k) )                     &
                                        + bo2Penta * ( f(         3,j,k) + f(         1,j,k) ) &
                                        + co2Penta * ( f(         4,j,k) + f(    this%n,j,k) ) &
                                        + do2Penta * ( f(         5,j,k) + f(  this%n-1,j,k) ) &
                                        + eo2Penta * ( f(         6,j,k) + f(  this%n-2,j,k) )
                    RHS(         3,j,k) = ao2Penta * ( f(         3,j,k) )                     &
                                        + bo2Penta * ( f(         4,j,k) + f(         2,j,k) ) &
                                        + co2Penta * ( f(         5,j,k) + f(         1,j,k) ) &
                                        + do2Penta * ( f(         6,j,k) + f(    this%n,j,k) ) &
                                        + eo2Penta * ( f(         7,j,k) + f(  this%n-1,j,k) )
                    RHS(         4,j,k) = ao2Penta * ( f(         4,j,k) )                     &
                                        + bo2Penta * ( f(         5,j,k) + f(         3,j,k) ) &
                                        + co2Penta * ( f(         6,j,k) + f(         2,j,k) ) &
                                        + do2Penta * ( f(         7,j,k) + f(         1,j,k) ) &
                                        + eo2Penta * ( f(         8,j,k) + f(    this%n,j,k) )
                    RHS(5:this%n-4,j,k) = ao2Penta * ( f(5:this%n-4,j,k) )                     &
                                        + bo2Penta * ( f(6:this%n-3,j,k) + f(4:this%n-5,j,k) ) &
                                        + co2Penta * ( f(7:this%n-2,j,k) + f(3:this%n-6,j,k) ) &
                                        + do2Penta * ( f(8:this%n-1,j,k) + f(2:this%n-7,j,k) ) &
                                        + eo2Penta * ( f(9:this%n  ,j,k) + f(1:this%n-8,j,k) )
                    RHS(  this%n-3,j,k) = ao2Penta * ( f(  this%n-3,j,k) )                     &
                                        + bo2Penta * ( f(  this%n-2,j,k) + f(  this%n-4,j,k) ) &
                                        + co2Penta * ( f(  this%n-1,j,k) + f(  this%n-5,j,k) ) &
                                        + do2Penta * ( f(    this%n,j,k) + f(  this%n-6,j,k) ) &
                                        + eo2Penta * ( f(         1,j,k) + f(  this%n-7,j,k) )
                    RHS(  this%n-2,j,k) = ao2Penta * ( f(  this%n-2,j,k) )                     &
                                        + bo2Penta * ( f(  this%n-1,j,k) + f(  this%n-3,j,k) ) &
                                        + co2Penta * ( f(    this%n,j,k) + f(  this%n-4,j,k) ) &
                                        + do2Penta * ( f(         1,j,k) + f(  this%n-5,j,k) ) &
                                        + eo2Penta * ( f(         2,j,k) + f(  this%n-6,j,k) )
                    RHS(  this%n-1,j,k) = ao2Penta * ( f(  this%n-1,j,k) )                     &
                                        + bo2Penta * ( f(    this%n,j,k) + f(  this%n-2,j,k) ) &
                                        + co2Penta * ( f(         1,j,k) + f(  this%n-3,j,k) ) &
                                        + do2Penta * ( f(         2,j,k) + f(  this%n-4,j,k) ) &
                                        + eo2Penta * ( f(         3,j,k) + f(  this%n-5,j,k) )
                    RHS(    this%n,j,k) = ao2Penta * ( f(    this%n,j,k) )                     &
                                        + bo2Penta * ( f(         1,j,k) + f(  this%n-1,j,k) ) &
                                        + co2Penta * ( f(         2,j,k) + f(  this%n-2,j,k) ) &
                                        + do2Penta * ( f(         3,j,k) + f(  this%n-3,j,k) ) &
                                        + eo2Penta * ( f(         4,j,k) + f(  this%n-4,j,k) )
                end do
            end do

        case (.FALSE.)

            do k = 1,n3
                do j = 1,n2
                    
                    select case(bc1)
                    case(0)
                        RHS(         1,j,k) =    one * ( f(         1,j,k) )                     

                        RHS(         2,j,k) = b2_ao2Penta * ( f(         2,j,k) )                     &
                                            + b2_bo2Penta * ( f(         3,j,k) + f(         1,j,k) ) 
                        
                        RHS(         3,j,k) = ao2Penta * ( f(         3,j,k) )                     &
                                            + bo2Penta * ( f(         4,j,k) + f(         2,j,k) ) &
                                            + co2Penta * ( f(         5,j,k) + f(         1,j,k) )

                        RHS(         4,j,k) = ao2Penta * ( f(         4,j,k) )                     &
                                            + bo2Penta * ( f(         5,j,k) + f(         3,j,k) ) &
                                            + co2Penta * ( f(         6,j,k) + f(         2,j,k) ) &
                                            + do2Penta * ( f(         7,j,k) + f(         1,j,k) ) 
                    case(1)
                        RHS(1,j,k) =    ao2Penta * ( f(1,j,k) )            &
                                   +    bo2Penta * ( f(2,j,k) + f(2,j,k) ) &
                                   +    co2Penta * ( f(3,j,k) + f(3,j,k) ) &
                                   +    do2Penta * ( f(4,j,k) + f(4,j,k) ) &
                                   +    eo2Penta * ( f(5,j,k) + f(5,j,k) )
                        RHS(2,j,k) =    ao2Penta * ( f(2,j,k) )            &
                                   +    bo2Penta * ( f(3,j,k) + f(1,j,k) ) &
                                   +    co2Penta * ( f(4,j,k) + f(2,j,k) ) &
                                   +    do2Penta * ( f(5,j,k) + f(3,j,k) ) &
                                   +    eo2Penta * ( f(6,j,k) + f(4,j,k) )
                        RHS(3,j,k) =    ao2Penta * ( f(3,j,k) )            &
                                   +    bo2Penta * ( f(4,j,k) + f(2,j,k) ) &
                                   +    co2Penta * ( f(5,j,k) + f(1,j,k) ) &
                                   +    do2Penta * ( f(6,j,k) + f(2,j,k) ) &
                                   +    eo2Penta * ( f(7,j,k) + f(3,j,k) )
                        RHS(4,j,k) =    ao2Penta * ( f(4,j,k) )            &
                                   +    bo2Penta * ( f(5,j,k) + f(3,j,k) ) &
                                   +    co2Penta * ( f(6,j,k) + f(2,j,k) ) &
                                   +    do2Penta * ( f(7,j,k) + f(1,j,k) ) &
                                   +    eo2Penta * ( f(8,j,k) + f(2,j,k) )
                    case(-1)
                        RHS(1,j,k) =    ao2Penta * ( f(1,j,k) )            &
                                   +    bo2Penta * ( f(2,j,k) - f(2,j,k) ) &
                                   +    co2Penta * ( f(3,j,k) - f(3,j,k) ) &
                                   +    do2Penta * ( f(4,j,k) - f(4,j,k) ) &
                                   +    eo2Penta * ( f(5,j,k) - f(5,j,k) )
                        RHS(2,j,k) =    ao2Penta * ( f(2,j,k) )            &
                                   +    bo2Penta * ( f(3,j,k) + f(1,j,k) ) &
                                   +    co2Penta * ( f(4,j,k) - f(2,j,k) ) &
                                   +    do2Penta * ( f(5,j,k) - f(3,j,k) ) &
                                   +    eo2Penta * ( f(6,j,k) - f(4,j,k) )
                        RHS(3,j,k) =    ao2Penta * ( f(3,j,k) )            &
                                   +    bo2Penta * ( f(4,j,k) + f(2,j,k) ) &
                                   +    co2Penta * ( f(5,j,k) + f(1,j,k) ) &
                                   +    do2Penta * ( f(6,j,k) - f(2,j,k) ) &
                                   +    eo2Penta * ( f(7,j,k) - f(3,j,k) )
                        RHS(4,j,k) =    ao2Penta * ( f(4,j,k) )            &
                                   +    bo2Penta * ( f(5,j,k) + f(3,j,k) ) &
                                   +    co2Penta * ( f(6,j,k) + f(2,j,k) ) &
                                   +    do2Penta * ( f(7,j,k) + f(1,j,k) ) &
                                   +    eo2Penta * ( f(8,j,k) - f(2,j,k) )
                    end select

                    RHS(5:this%n-4,j,k) =    ao2Penta * ( f(5:this%n-4,j,k) )                     &
                                        +    bo2Penta * ( f(6:this%n-3,j,k) + f(4:this%n-5,j,k) ) &
                                        +    co2Penta * ( f(7:this%n-2,j,k) + f(3:this%n-6,j,k) ) &
                                        +    do2Penta * ( f(8:this%n-1,j,k) + f(2:this%n-7,j,k) ) &
                                        +    eo2Penta * ( f(9:this%n  ,j,k) + f(1:this%n-8,j,k) )

                    select case(bcn)
                    case(0)
                        RHS(  this%n-3,j,k) = ao2Penta * ( f(  this%n-3,j,k) )                     &
                                            + bo2Penta * ( f(  this%n-2,j,k) + f(  this%n-4,j,k) ) &
                                            + co2Penta * ( f(  this%n-1,j,k) + f(  this%n-5,j,k) ) &
                                            + do2Penta * ( f(    this%n,j,k) + f(  this%n-6,j,k) ) 

                        RHS(  this%n-2,j,k) = ao2Penta * ( f(  this%n-2,j,k) )                     &
                                            + bo2Penta * ( f(  this%n-1,j,k) + f(  this%n-3,j,k) ) &
                                            + co2Penta * ( f(    this%n,j,k) + f(  this%n-4,j,k) ) 

                        RHS(  this%n-1,j,k) = b2_ao2Penta * ( f(  this%n-1,j,k) )                     &
                                            + b2_bo2Penta * ( f(    this%n,j,k) + f(  this%n-2,j,k) ) 

                        RHS(    this%n,j,k) =    one * ( f(    this%n,j,k) )                     
                    case(1)
                        RHS(this%n-3,j,k) =    ao2Penta * ( f(this%n-3,j,k) )                   &
                                          +    bo2Penta * ( f(this%n-2,j,k) + f(this%n-4,j,k) ) &
                                          +    co2Penta * ( f(this%n-1,j,k) + f(this%n-5,j,k) ) &
                                          +    do2Penta * ( f(this%n  ,j,k) + f(this%n-6,j,k) ) &
                                          +    eo2Penta * ( f(this%n-1,j,k) + f(this%n-7,j,k) )
                        RHS(this%n-2,j,k) =    ao2Penta * ( f(this%n-2,j,k) )                   &
                                          +    bo2Penta * ( f(this%n-1,j,k) + f(this%n-3,j,k) ) &
                                          +    co2Penta * ( f(this%n  ,j,k) + f(this%n-4,j,k) ) &
                                          +    do2Penta * ( f(this%n-1,j,k) + f(this%n-5,j,k) ) &
                                          +    eo2Penta * ( f(this%n-2,j,k) + f(this%n-6,j,k) )
                        RHS(this%n-1,j,k) =    ao2Penta * ( f(this%n-1,j,k) )                   &
                                          +    bo2Penta * ( f(this%n  ,j,k) + f(this%n-2,j,k) ) &
                                          +    co2Penta * ( f(this%n-1,j,k) + f(this%n-3,j,k) ) &
                                          +    do2Penta * ( f(this%n-2,j,k) + f(this%n-4,j,k) ) &
                                          +    eo2Penta * ( f(this%n-3,j,k) + f(this%n-5,j,k) )
                        RHS(this%n  ,j,k) =    ao2Penta * ( f(this%n  ,j,k) )                   &
                                          +    bo2Penta * ( f(this%n-1,j,k) + f(this%n-1,j,k) ) &
                                          +    co2Penta * ( f(this%n-2,j,k) + f(this%n-2,j,k) ) &
                                          +    do2Penta * ( f(this%n-3,j,k) + f(this%n-3,j,k) ) &
                                          +    eo2Penta * ( f(this%n-4,j,k) + f(this%n-4,j,k) )
                    case(-1)
                        RHS(this%n-3,j,k) =    ao2Penta * ( f(this%n-3,j,k) )                   &
                                          +    bo2Penta * ( f(this%n-2,j,k) + f(this%n-4,j,k) ) &
                                          +    co2Penta * ( f(this%n-1,j,k) + f(this%n-5,j,k) ) &
                                          +    do2Penta * ( f(this%n  ,j,k) + f(this%n-6,j,k) ) &
                                          +    eo2Penta * (-f(this%n-1,j,k) + f(this%n-7,j,k) )
                        RHS(this%n-2,j,k) =    ao2Penta * ( f(this%n-2,j,k) )                   &
                                          +    bo2Penta * ( f(this%n-1,j,k) + f(this%n-3,j,k) ) &
                                          +    co2Penta * ( f(this%n  ,j,k) + f(this%n-4,j,k) ) &
                                          +    do2Penta * (-f(this%n-1,j,k) + f(this%n-5,j,k) ) &
                                          +    eo2Penta * (-f(this%n-2,j,k) + f(this%n-6,j,k) )
                        RHS(this%n-1,j,k) =    ao2Penta * ( f(this%n-1,j,k) )                   &
                                          +    bo2Penta * ( f(this%n  ,j,k) + f(this%n-2,j,k) ) &
                                          +    co2Penta * (-f(this%n-1,j,k) + f(this%n-3,j,k) ) &
                                          +    do2Penta * (-f(this%n-2,j,k) + f(this%n-4,j,k) ) &
                                          +    eo2Penta * (-f(this%n-3,j,k) + f(this%n-5,j,k) )
                        RHS(this%n  ,j,k) =    ao2Penta * ( f(this%n  ,j,k) )                   &
                                          +    bo2Penta * (-f(this%n-1,j,k) + f(this%n-1,j,k) ) &
                                          +    co2Penta * (-f(this%n-2,j,k) + f(this%n-2,j,k) ) &
                                          +    do2Penta * (-f(this%n-3,j,k) + f(this%n-3,j,k) ) &
                                          +    eo2Penta * (-f(this%n-4,j,k) + f(this%n-4,j,k) )
                    end select
                
               end do 
            end do 
        end select
    
    end subroutine
    
    pure subroutine ComputeYRHS(this, f, RHS, n1, n3, bc1, bcn) 
    
        class( cfo2Penta ), intent(in) :: this
        integer, intent(in) :: n1, n3
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn
        integer :: k

        select case (this%periodic)
        case (.TRUE.)
            do k=1,n3
                RHS(:,         1,k) = ao2Penta * ( f(:,         1,k) )                     &
                                    + bo2Penta * ( f(:,         2,k) + f(:,    this%n,k) ) &
                                    + co2Penta * ( f(:,         3,k) + f(:,  this%n-1,k) ) &
                                    + do2Penta * ( f(:,         4,k) + f(:,  this%n-2,k) ) &
                                    + eo2Penta * ( f(:,         5,k) + f(:,  this%n-3,k) )
                RHS(:,         2,k) = ao2Penta * ( f(:,         2,k) )                     &
                                    + bo2Penta * ( f(:,         3,k) + f(:,         1,k) ) &
                                    + co2Penta * ( f(:,         4,k) + f(:,    this%n,k) ) &
                                    + do2Penta * ( f(:,         5,k) + f(:,  this%n-1,k) ) &
                                    + eo2Penta * ( f(:,         6,k) + f(:,  this%n-2,k) )
                RHS(:,         3,k) = ao2Penta * ( f(:,         3,k) )                     &
                                    + bo2Penta * ( f(:,         4,k) + f(:,         2,k) ) &
                                    + co2Penta * ( f(:,         5,k) + f(:,         1,k) ) &
                                    + do2Penta * ( f(:,         6,k) + f(:,    this%n,k) ) &
                                    + eo2Penta * ( f(:,         7,k) + f(:,  this%n-1,k) )
                RHS(:,         4,k) = ao2Penta * ( f(:,         4,k) )                     &
                                    + bo2Penta * ( f(:,         5,k) + f(:,         3,k) ) &
                                    + co2Penta * ( f(:,         6,k) + f(:,         2,k) ) &
                                    + do2Penta * ( f(:,         7,k) + f(:,         1,k) ) &
                                    + eo2Penta * ( f(:,         8,k) + f(:,    this%n,k) )
                RHS(:,5:this%n-4,k) = ao2Penta * ( f(:,5:this%n-4,k) )                     &
                                    + bo2Penta * ( f(:,6:this%n-3,k) + f(:,4:this%n-5,k) ) &
                                    + co2Penta * ( f(:,7:this%n-2,k) + f(:,3:this%n-6,k) ) &
                                    + do2Penta * ( f(:,8:this%n-1,k) + f(:,2:this%n-7,k) ) &
                                    + eo2Penta * ( f(:,9:this%n  ,k) + f(:,1:this%n-8,k) )
                RHS(:,  this%n-3,k) = ao2Penta * ( f(:,  this%n-3,k) )                     &
                                    + bo2Penta * ( f(:,  this%n-2,k) + f(:,  this%n-4,k) ) &
                                    + co2Penta * ( f(:,  this%n-1,k) + f(:,  this%n-5,k) ) &
                                    + do2Penta * ( f(:,    this%n,k) + f(:,  this%n-6,k) ) &
                                    + eo2Penta * ( f(:,         1,k) + f(:,  this%n-7,k) )
                RHS(:,  this%n-2,k) = ao2Penta * ( f(:,  this%n-2,k) )                     &
                                    + bo2Penta * ( f(:,  this%n-1,k) + f(:,  this%n-3,k) ) &
                                    + co2Penta * ( f(:,    this%n,k) + f(:,  this%n-4,k) ) &
                                    + do2Penta * ( f(:,         1,k) + f(:,  this%n-5,k) ) &
                                    + eo2Penta * ( f(:,         2,k) + f(:,  this%n-6,k) )
                RHS(:,  this%n-1,k) = ao2Penta * ( f(:,  this%n-1,k) )                     &
                                    + bo2Penta * ( f(:,    this%n,k) + f(:,  this%n-2,k) ) &
                                    + co2Penta * ( f(:,         1,k) + f(:,  this%n-3,k) ) &
                                    + do2Penta * ( f(:,         2,k) + f(:,  this%n-4,k) ) &
                                    + eo2Penta * ( f(:,         3,k) + f(:,  this%n-5,k) )
                RHS(:,    this%n,k) = ao2Penta * ( f(:,    this%n,k) )                     &
                                    + bo2Penta * ( f(:,         1,k) + f(:,  this%n-1,k) ) &
                                    + co2Penta * ( f(:,         2,k) + f(:,  this%n-2,k) ) &
                                    + do2Penta * ( f(:,         3,k) + f(:,  this%n-3,k) ) &
                                    + eo2Penta * ( f(:,         4,k) + f(:,  this%n-4,k) )
            end do
        case (.FALSE.)

            do k = 1,n3
                select case(bc1)
                case(0)
                    RHS(:,         1,k) =    one * ( f(:,         1,k) )                     

                    RHS(:,         2,k) = b2_ao2Penta * ( f(:,         2,k) )                     &
                                        + b2_bo2Penta * ( f(:,         3,k) + f(:,         1,k) ) 
                    
                    RHS(:,         3,k) = ao2Penta * ( f(:,         3,k) )                     &
                                        + bo2Penta * ( f(:,         4,k) + f(:,         2,k) ) &
                                        + co2Penta * ( f(:,         5,k) + f(:,         1,k) )

                    RHS(:,         4,k) = ao2Penta * ( f(:,         4,k) )                     &
                                        + bo2Penta * ( f(:,         5,k) + f(:,         3,k) ) &
                                        + co2Penta * ( f(:,         6,k) + f(:,         2,k) ) &
                                        + do2Penta * ( f(:,         7,k) + f(:,         1,k) ) 
                case(1)
                    RHS(:,1,k) =    ao2Penta * ( f(:,1,k) )            &
                               +    bo2Penta * ( f(:,2,k) + f(:,2,k) ) &
                               +    co2Penta * ( f(:,3,k) + f(:,3,k) ) &
                               +    do2Penta * ( f(:,4,k) + f(:,4,k) ) &
                               +    eo2Penta * ( f(:,5,k) + f(:,5,k) )
                    RHS(:,2,k) =    ao2Penta * ( f(:,2,k) )            &
                               +    bo2Penta * ( f(:,3,k) + f(:,1,k) ) &
                               +    co2Penta * ( f(:,4,k) + f(:,2,k) ) &
                               +    do2Penta * ( f(:,5,k) + f(:,3,k) ) &
                               +    eo2Penta * ( f(:,6,k) + f(:,4,k) )
                    RHS(:,3,k) =    ao2Penta * ( f(:,3,k) )            &
                               +    bo2Penta * ( f(:,4,k) + f(:,2,k) ) &
                               +    co2Penta * ( f(:,5,k) + f(:,1,k) ) &
                               +    do2Penta * ( f(:,6,k) + f(:,2,k) ) &
                               +    eo2Penta * ( f(:,7,k) + f(:,3,k) )
                    RHS(:,4,k) =    ao2Penta * ( f(:,4,k) )            &
                               +    bo2Penta * ( f(:,5,k) + f(:,3,k) ) &
                               +    co2Penta * ( f(:,6,k) + f(:,2,k) ) &
                               +    do2Penta * ( f(:,7,k) + f(:,1,k) ) &
                               +    eo2Penta * ( f(:,8,k) + f(:,2,k) )
                case(-1)  
                    RHS(:,1,k) =    ao2Penta * ( f(:,1,k) )            &
                               +    bo2Penta * ( f(:,2,k) - f(:,2,k) ) &
                               +    co2Penta * ( f(:,3,k) - f(:,3,k) ) &
                               +    do2Penta * ( f(:,4,k) - f(:,4,k) ) &
                               +    eo2Penta * ( f(:,5,k) - f(:,5,k) )
                    RHS(:,2,k) =    ao2Penta * ( f(:,2,k) )            &
                               +    bo2Penta * ( f(:,3,k) + f(:,1,k) ) &
                               +    co2Penta * ( f(:,4,k) - f(:,2,k) ) &
                               +    do2Penta * ( f(:,5,k) - f(:,3,k) ) &
                               +    eo2Penta * ( f(:,6,k) - f(:,4,k) )
                    RHS(:,3,k) =    ao2Penta * ( f(:,3,k) )            &
                               +    bo2Penta * ( f(:,4,k) + f(:,2,k) ) &
                               +    co2Penta * ( f(:,5,k) + f(:,1,k) ) &
                               +    do2Penta * ( f(:,6,k) - f(:,2,k) ) &
                               +    eo2Penta * ( f(:,7,k) - f(:,3,k) )
                    RHS(:,4,k) =    ao2Penta * ( f(:,4,k) )            &
                               +    bo2Penta * ( f(:,5,k) + f(:,3,k) ) &
                               +    co2Penta * ( f(:,6,k) + f(:,2,k) ) &
                               +    do2Penta * ( f(:,7,k) + f(:,1,k) ) &
                               +    eo2Penta * ( f(:,8,k) - f(:,2,k) )
                end select

                RHS(:,5:this%n-4,k) =    ao2Penta * ( f(:,5:this%n-4,k) )                     &
                                    +    bo2Penta * ( f(:,6:this%n-3,k) + f(:,4:this%n-5,k) ) &
                                    +    co2Penta * ( f(:,7:this%n-2,k) + f(:,3:this%n-6,k) ) &
                                    +    do2Penta * ( f(:,8:this%n-1,k) + f(:,2:this%n-7,k) ) &
                                    +    eo2Penta * ( f(:,9:this%n  ,k) + f(:,1:this%n-8,k) )

                select case(bcn)
                case(0)
                    RHS(:,  this%n-3,k) = ao2Penta * ( f(:,  this%n-3,k) )                     &
                                        + bo2Penta * ( f(:,  this%n-2,k) + f(:,  this%n-4,k) ) &
                                        + co2Penta * ( f(:,  this%n-1,k) + f(:,  this%n-5,k) ) &
                                        + do2Penta * ( f(:,    this%n,k) + f(:,  this%n-6,k) ) 

                    RHS(:,  this%n-2,k) = ao2Penta * ( f(:,  this%n-2,k) )                     &
                                        + bo2Penta * ( f(:,  this%n-1,k) + f(:,  this%n-3,k) ) &
                                        + co2Penta * ( f(:,    this%n,k) + f(:,  this%n-4,k) ) 

                    RHS(:,  this%n-1,k) = b2_ao2Penta * ( f(:,  this%n-1,k) )                     &
                                        + b2_bo2Penta * ( f(:,    this%n,k) + f(:,  this%n-2,k) ) 

                    RHS(:,    this%n,k) =    one * ( f(:,    this%n,k) )                     
                case(1)
                    RHS(:,this%n-3,k) =    ao2Penta * ( f(:,this%n-3,k) )                   &
                                      +    bo2Penta * ( f(:,this%n-2,k) + f(:,this%n-4,k) ) &
                                      +    co2Penta * ( f(:,this%n-1,k) + f(:,this%n-5,k) ) &
                                      +    do2Penta * ( f(:,this%n  ,k) + f(:,this%n-6,k) ) &
                                      +    eo2Penta * ( f(:,this%n-1,k) + f(:,this%n-7,k) )
                    RHS(:,this%n-2,k) =    ao2Penta * ( f(:,this%n-2,k) )                   &
                                      +    bo2Penta * ( f(:,this%n-1,k) + f(:,this%n-3,k) ) &
                                      +    co2Penta * ( f(:,this%n  ,k) + f(:,this%n-4,k) ) &
                                      +    do2Penta * ( f(:,this%n-1,k) + f(:,this%n-5,k) ) &
                                      +    eo2Penta * ( f(:,this%n-2,k) + f(:,this%n-6,k) )
                    RHS(:,this%n-1,k) =    ao2Penta * ( f(:,this%n-1,k) )                   &
                                      +    bo2Penta * ( f(:,this%n  ,k) + f(:,this%n-2,k) ) &
                                      +    co2Penta * ( f(:,this%n-1,k) + f(:,this%n-3,k) ) &
                                      +    do2Penta * ( f(:,this%n-2,k) + f(:,this%n-4,k) ) &
                                      +    eo2Penta * ( f(:,this%n-3,k) + f(:,this%n-5,k) )
                    RHS(:,this%n  ,k) =    ao2Penta * ( f(:,this%n  ,k) )                   &
                                      +    bo2Penta * ( f(:,this%n-1,k) + f(:,this%n-1,k) ) &
                                      +    co2Penta * ( f(:,this%n-2,k) + f(:,this%n-2,k) ) &
                                      +    do2Penta * ( f(:,this%n-3,k) + f(:,this%n-3,k) ) &
                                      +    eo2Penta * ( f(:,this%n-4,k) + f(:,this%n-4,k) )
                case(-1)  
                    RHS(:,this%n-3,k) =    ao2Penta * ( f(:,this%n-3,k) )                   &
                                      +    bo2Penta * ( f(:,this%n-2,k) + f(:,this%n-4,k) ) &
                                      +    co2Penta * ( f(:,this%n-1,k) + f(:,this%n-5,k) ) &
                                      +    do2Penta * ( f(:,this%n  ,k) + f(:,this%n-6,k) ) &
                                      +    eo2Penta * (-f(:,this%n-1,k) + f(:,this%n-7,k) )
                    RHS(:,this%n-2,k) =    ao2Penta * ( f(:,this%n-2,k) )                   &
                                      +    bo2Penta * ( f(:,this%n-1,k) + f(:,this%n-3,k) ) &
                                      +    co2Penta * ( f(:,this%n  ,k) + f(:,this%n-4,k) ) &
                                      +    do2Penta * (-f(:,this%n-1,k) + f(:,this%n-5,k) ) &
                                      +    eo2Penta * (-f(:,this%n-2,k) + f(:,this%n-6,k) )
                    RHS(:,this%n-1,k) =    ao2Penta * ( f(:,this%n-1,k) )                   &
                                      +    bo2Penta * ( f(:,this%n  ,k) + f(:,this%n-2,k) ) &
                                      +    co2Penta * (-f(:,this%n-1,k) + f(:,this%n-3,k) ) &
                                      +    do2Penta * (-f(:,this%n-2,k) + f(:,this%n-4,k) ) &
                                      +    eo2Penta * (-f(:,this%n-3,k) + f(:,this%n-5,k) )
                    RHS(:,this%n  ,k) =    ao2Penta * ( f(:,this%n  ,k) )                   &
                                      +    bo2Penta * (-f(:,this%n-1,k) + f(:,this%n-1,k) ) &
                                      +    co2Penta * (-f(:,this%n-2,k) + f(:,this%n-2,k) ) &
                                      +    do2Penta * (-f(:,this%n-3,k) + f(:,this%n-3,k) ) &
                                      +    eo2Penta * (-f(:,this%n-4,k) + f(:,this%n-4,k) )
                end select
                
            end do 
        end select
    
    end subroutine

    pure subroutine ComputeZRHS_REAL(this, f, RHS, n1, n2, bc1, bcn)
    
        class( cfo2Penta ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in) :: f
        real(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn

#include "CFo2Penta_files/ComputeZRHS_common.F90"    
    
    end subroutine
    
    pure subroutine ComputeZRHS_CMPLX(this, f, RHS, n1, n2, bc1, bcn)
    
        class( cfo2Penta ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in) :: f
        complex(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn

#include "CFo2Penta_files/ComputeZRHS_common.F90"    
    
    end subroutine

    subroutine filter1(this, f, df, na, nb, bc1_, bcn_)
        class( cfo2Penta ), intent(in) :: this
        integer, intent(in) :: na, nb
        real(rkind), dimension(this%n,na,nb), intent(in) :: f
        real(rkind), dimension(this%n,na,nb), intent(out) :: df
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn

        if(this%n == 1) then
            df = f
            return
        end if
        
        if (present(bc1_)) then
            bc1 = bc1_
            if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bc1 (should be 0, 1 or -1)", 324)
            end if
        else
            bc1 = 0
        end if

        if (present(bcn_)) then
            bcn = bcn_
            if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bcn (should be 0, 1 or -1)", 324)
            end if
        else
            bcn = 0
        end if

        call this%ComputeXRHS(f, df, na, nb, bc1, bcn)

        select case (this%periodic)
        case(.TRUE.)
            call this%SolveXLU(df, na, nb)
        case(.FALSE.)
            select case(bc1)
            case(0) ! Normal non-periodic left boundary
                select case(bcn)
                case(0)  ! Normal non-periodic right boundary
                    call this%SolveXPenta(this%penta_nn, df, na, nb)
                case(1) ! Symmetric right boundary
                    call this%SolveXPenta(this%penta_ns, df, na, nb)
                case(-1) ! Antisymmetric right boundary
                    call this%SolveXPenta(this%penta_na, df, na, nb)
                end select
            case(1) ! Symmetric left boundary
                select case(bcn)
                case(0)  ! Normal non-periodic right boundary
                    call this%SolveXPenta(this%penta_sn, df, na, nb)
                case(1)  ! Symmetric right boundary
                    call this%SolveXPenta(this%penta_ss, df, na, nb)
                case(-1) ! Antisymmetric right boundary
                    call this%SolveXPenta(this%penta_sa, df, na, nb)
                end select
            case(-1) ! Antisymmetric left boundary
                select case(bcn)
                case(0)  ! Normal non-periodic right boundary
                    call this%SolveXPenta(this%penta_an, df, na, nb)
                case(1) ! Symmetric right boundary
                    call this%SolveXPenta(this%penta_as, df, na, nb)
                case(-1) ! Antisymmetric right boundary
                    call this%SolveXPenta(this%penta_aa, df, na, nb)
                end select
            end select
        end select
    
    end subroutine

    subroutine filter2(this, f, df, na, nb, bc1_, bcn_)
        class( cfo2Penta ), intent(in) :: this
        integer, intent(in) :: na, nb
        real(rkind), dimension(na,this%n,nb), intent(in) :: f
        real(rkind), dimension(na,this%n,nb), intent(out) :: df
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn

        if(this%n == 1) then
            df = f
            return
        end if
        
        if (present(bc1_)) then
            bc1 = bc1_
            if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bc1 (should be 0, 1 or -1)", 324)
            end if
        else
            bc1 = 0
        end if

        if (present(bcn_)) then
            bcn = bcn_
            if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bcn (should be 0, 1 or -1)", 324)
            end if
        else
            bcn = 0
        end if

        call this%ComputeYRHS(f, df, na, nb, bc1, bcn)

        select case (this%periodic)
        case(.TRUE.)
            call this%SolveYLU(df, na, nb)
        case(.FALSE.)
            select case(bc1)
            case(0) ! Normal non-periodic left boundary
                select case(bcn)
                case(0)  ! Normal non-periodic right boundary
                    call this%SolveYPenta(this%penta_nn, df, na, nb)
                case(1) ! Symmetric right boundary
                    call this%SolveYPenta(this%penta_ns, df, na, nb)
                case(-1) ! Antisymmetric right boundary
                    call this%SolveYPenta(this%penta_na, df, na, nb)
                end select
            case(1) ! Symmetric left boundary
                select case(bcn)
                case(0)  ! Normal non-periodic right boundary
                    call this%SolveYPenta(this%penta_sn, df, na, nb)
                case(1)  ! Symmetric right boundary
                    call this%SolveYPenta(this%penta_ss, df, na, nb)
                case(-1) ! Antisymmetric right boundary
                    call this%SolveYPenta(this%penta_sa, df, na, nb)
                end select
            case(-1) ! Antisymmetric left boundary
                select case(bcn)
                case(0)  ! Normal non-periodic right boundary
                    call this%SolveYPenta(this%penta_an, df, na, nb)
                case(1) ! Symmetric right boundary
                    call this%SolveYPenta(this%penta_as, df, na, nb)
                case(-1) ! Antisymmetric right boundary
                    call this%SolveYPenta(this%penta_aa, df, na, nb)
                end select
            end select
        end select
    
    end subroutine

    subroutine filter3_REAL(this, f, df, na, nb, bc1_, bcn_)
        class( cfo2Penta ), intent(in) :: this
        integer, intent(in) :: na, nb
        real(rkind), dimension(na,nb,this%n), intent(in) :: f
        real(rkind), dimension(na,nb,this%n), intent(out) :: df
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn

        if(this%n == 1) then
            df = f
            return
        end if
        
        if (present(bc1_)) then
            bc1 = bc1_
            if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bc1 (should be 0, 1 or -1)", 324)
            end if
        else
            bc1 = 0
        end if

        if (present(bcn_)) then
            bcn = bcn_
            if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bcn (should be 0, 1 or -1)", 324)
            end if
        else
            bcn = 0
        end if

        call this%ComputeZRHS_REAL(f, df, na, nb, bc1, bcn)

        select case (this%periodic)
        case(.TRUE.)
            call this%SolveZLU(df, na, nb)
        case(.FALSE.)
            select case(bc1)
            case(0) ! Normal non-periodic left boundary
                select case(bcn)
                case(0)  ! Normal non-periodic right boundary
                    call this%SolveZPenta_REAL(this%penta_nn, df, na, nb)
                case(1) ! Symmetric right boundary
                    call this%SolveZPenta_REAL(this%penta_ns, df, na, nb)
                case(-1) ! Antisymmetric right boundary
                    call this%SolveZPenta_REAL(this%penta_na, df, na, nb)
                end select
            case(1) ! Symmetric left boundary
                select case(bcn)
                case(0)  ! Normal non-periodic right boundary
                    call this%SolveZPenta_REAL(this%penta_sn, df, na, nb)
                case(1)  ! Symmetric right boundary
                    call this%SolveZPenta_REAL(this%penta_ss, df, na, nb)
                case(-1) ! Antisymmetric right boundary
                    call this%SolveZPenta_REAL(this%penta_sa, df, na, nb)
                end select
            case(-1) ! Antisymmetric left boundary
                select case(bcn)
                case(0)  ! Normal non-periodic right boundary
                    call this%SolveZPenta_REAL(this%penta_an, df, na, nb)
                case(1) ! Symmetric right boundary
                    call this%SolveZPenta_REAL(this%penta_as, df, na, nb)
                case(-1) ! Antisymmetric right boundary
                    call this%SolveZPenta_REAL(this%penta_aa, df, na, nb)
                end select
            end select
        end select
    
    end subroutine

    subroutine filter3_CMPLX(this, f, df, na, nb, bc1_, bcn_)
        class( cfo2Penta ), intent(in) :: this
        integer, intent(in) :: na, nb
        complex(rkind), dimension(na,nb,this%n), intent(in) :: f
        complex(rkind), dimension(na,nb,this%n), intent(out) :: df
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn

        if(this%n == 1) then
            df = f
            return
        end if
        if (present(bcn_)) then
            bc1 = bc1_; bcn = bcn_
        else 
            bc1 = 0; bcn = 0
        end if 

        call this%ComputeZRHS_CMPLX(f, df, na, nb, bc1, bcn)
        call this%SolveZPenta_CMPLX(this%penta_nn, df, na, nb)
    
    end subroutine

end module
