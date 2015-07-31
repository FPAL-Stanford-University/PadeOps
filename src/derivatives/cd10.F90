! Routines specific to 10th order Compact Finite Differencing scheme
! Periodic LU based on  Neossi Nguetchue, Abelman (Appl. Math. & Comp. 2008)

module cd10stuff

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two
    
    implicit none

    private
    public :: cd10
    
    ! 10th order first derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha10d1=  1.0_rkind /  2.0_rkind
    real(rkind), parameter :: beta10d1 =  1.0_rkind / 20.0_rkind
    real(rkind), parameter :: a10d1    =( 17.0_rkind / 12.0_rkind) / 2.0_rkind
    real(rkind), parameter :: b10d1    =(101.0_rkind /150.0_rkind) / 4.0_rkind
    real(rkind), parameter :: c10d1    =(  1.0_rkind /100.0_rkind) / 6.0_rkind
    
    ! 10th order second derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha10d2=334.0_rkind /899.0_rkind
    real(rkind), parameter :: beta10d2 = 43.0_rkind /1798.0_rkind
    real(rkind), parameter :: a10d2    =(1065.0_rkind /1798.0_rkind) / 1.0_rkind
    real(rkind), parameter :: b10d2    =(1038.0_rkind / 899.0_rkind) / 4.0_rkind
    real(rkind), parameter :: c10d2    =( 79.0_rkind /1798.0_rkind ) / 9.0_rkind

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic 1st derivative evaluation !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
    ! Set the scheme for the edge nodes (Ref. for notation: Lele - JCP paper)
    ! 1st derivative 
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
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic 2nd derivative evaluation !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    ! Incomplete
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    type cd10
        
        private
        
        integer     :: n
        real(rkind) :: dx
        real(rkind) :: onebydx
        real(rkind) :: onebydx2

        logical     :: periodic=.TRUE.
        integer     :: bc1=0                               ! Boundary condition type. 0=Dirichlet, 1=Neumann
        integer     :: bcn=0                               ! Boundary condition type. 0=Dirichlet, 1=Neumann 

        real(rkind), allocatable, dimension(:,:) :: LU1
        real(rkind), allocatable, dimension(:,:) :: LU2
        real(rkind), allocatable, dimension(:,:) :: penta1
        real(rkind), allocatable, dimension(:,:) :: penta2

        contains

        procedure :: init
        procedure :: destroy
        procedure :: GetSize
        procedure, private :: ComputeXD1RHS
        procedure, private :: ComputeYD1RHS
        procedure, private :: ComputeZD1RHS
        
        procedure, private :: ComputeD2RHS

        procedure, private :: SolveXLU1
        procedure, private :: SolveYLU1
        procedure, private :: SolveZLU1
        
        procedure, private :: SolveLU2
       
        procedure, private :: ComputePenta1
        procedure, private :: ComputePenta2

        procedure, private :: SolveXPenta1
        procedure, private :: SolveYPenta1
        procedure, private :: SolveZPenta1
        
        procedure, private :: SolvePenta2
        
        procedure :: dd1
        procedure :: dd2
        procedure :: dd3
        
        procedure :: cd10der2

    end type



contains
    
    pure function GetSize(this) result(val)
        class(cd10), intent(in) :: this
        integer  :: val 
        val = this%n
    end function

    function init(this, n_, dx_, periodic_, bc1_, bcn_) result(ierr)
   
        class( cd10 ), intent(inout) :: this
        integer, intent(in) :: n_
        real(rkind), intent(in) :: dx_
        logical, intent(in) :: periodic_
        integer, intent(in) :: bc1_, bcn_
        integer :: ierr
        
        this%n = n_
        this%dx = dx_
        this%onebydx = one/dx_
        this%onebydx2 = this%onebydx/dx_

        this%periodic = periodic_

        this%bc1 = bc1_
        this%bcn = bcn_

        if (periodic_) then 
            ! Allocate 1st derivative LU matrix.
            if(allocated( this%LU1 )) deallocate( this%LU1 ); allocate( this%LU1(n_,9) )
    
            ! Compute 1st derivative LU matrix
            if (n_ .GE. 8) then
                call ComputeLU(this%LU1,n_,beta10d1,alpha10d1,one,alpha10d1,beta10d1)
            else if (n_ == 1) then
                this%LU1 = one
            else
                ierr = 2
                return
            end if
    
            ! Allocate 2nd derivative LU matrices.
            if(allocated( this%LU2 )) deallocate( this%LU2 ); allocate( this%LU2(n_,9) )
    
            ! Compute 2nd derivative LU matrix
            if (n_ .GE. 8) then
                call ComputeLU(this%LU2,n_,beta10d2,alpha10d2,one,alpha10d2,beta10d2)
            else if (n_ == 1) then
                this%LU2 = one
            end if
        else 
            ! Allocate 1st derivative Penta matrices.
            if(allocated( this%penta1 )) deallocate( this%penta1 ); allocate( this%penta1(n_,11) )
  
            if (n_ .GE. 8) then             
                call this%ComputePenta1(bc1_,bcn_)
            else if (n_ .EQ. 1) then
                this%Penta1 = one
            else
                ierr = 2
                return 
            end if 

            
            ! Allocate 2nd derivative Penta matrices.
            if(allocated( this%penta2 )) deallocate( this%penta2 ); allocate( this%penta2(n_,11) )
        
            if (n_ .GE. 8) then             
                call this%ComputePenta2(bc1_,bcn_)
            else if (n_ .EQ. 1) then
                this%Penta1 = one
            else
                ierr = 2
                return 
            end if 

        end if 

        ! If everything passes
        ierr = 0
    
    end function
    
    subroutine destroy(this)

        class( cd10 ), intent(inout) :: this

        ! Dellocate 1st derivative LU matrix.
        if(allocated( this%LU1 )) deallocate( this%LU1 )
    
        ! Dellocate 2nd derivative LU matrix.
        if(allocated( this%LU2 )) deallocate( this%LU2 )
    
        ! Dellocate 1st derivative penta matrix.
        if(allocated( this%penta1 )) deallocate( this%penta1 )
    
        ! Dellocate 2nd derivative penta matrix.
        if(allocated( this%penta2 )) deallocate( this%penta2 )
    
    end subroutine

    subroutine ComputeLU(LU,n,e,a,d,c,f) 
    
        integer, intent(in) :: n
        real(rkind), intent(in) :: d,a,c,e,f
        real(rkind), dimension(n,9), intent(out) :: LU
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
            do i = 3,n-3
                b(i) = ( a - ( e/g(i-2) )*h(i-2) ) / g(i-1)
                h(i) = c - b(i)*f
                g(i) = d - ( e/g(i-2) )*f - b(i)*h(i-1)
            end do
    
            ! Step 3
            b(n-2) = ( a - ( e/g(n-4) )*h(n-4) ) / g(n-3)
            g(n-2) = d - ( e/g(n-4) )*f - b(n-2)*h(n-3)
    
            ! Step 4
            do i = 3,n-4
                k(i) = -( k(i-2)*f + k(i-1)*h(i-1) )/g(i)
                v(i) = -( e/g(i-2) )*v(i-2) - b(i)*v(i-1)
            end do
    
            ! Step 5
            k(n-3) = ( e - k(n-5)*f - k(n-4)*h(n-4) ) / g(n-3)
            k(n-2) = ( a - k(n-4)*f - k(n-3)*h(n-3) ) / g(n-2)
            v(n-3) = f - ( e/g(n-5) )*v(n-5) - b(n-3)*v(n-4)
            v(n-2) = c - ( e/g(n-4) )*v(n-4) - b(n-2)*v(n-3)
            g(n-1) = d - SUM( k(1:n-2)*v(1:n-2) )
    
            ! Step 6
            do i = 3,n-3
                w(i) = -( e/g(i-2) )*w(i-2) - b(i)*w(i-1)
                l(i) = -( l(i-2)*f + l(i-1)*h(i-1) ) / g(i)
            end do
    
            ! Step 7
            w(n-2) = f - ( e/g(n-4) )*w(n-4) - b(n-2)*w(n-3)
            w(n-1) = c - SUM( k(1:n-2)*w(1:n-2) )
            l(n-2) = ( e - l(n-4)*f - l(n-3)*h(n-3) ) / g(n-2)
            l(n-1) = ( a - SUM( l(1:n-2)*v(1:n-2) ) ) / g(n-1)
            g(n)   = d - SUM( l(1:n-1)*w(1:n-1) )
    
            ! Set eg(i) = e/g(i-2)
            eg(3:n-2) = e/g(1:n-4)
    
            ! Set ff = f
            ff(1:n-4) = f
    
            ! Set g = 1/g
            g = 1._rkind/g
    
        end associate
    
    end subroutine

    subroutine ComputePenta1(this,bc1,bcn)
        class (cd10), intent(inout) :: this
        integer, intent(in) :: bc1, bcn
        integer             :: i
    
        associate (bt   => this%penta1(:,1), b   => this%penta1(:,2), d => this%penta1(:,3),  &
                   a    => this%penta1(:,4), at  => this%penta1(:,5),                         &
                   e    => this%penta1(:,6), obc => this%penta1(:,7),                         &
                   f    => this%penta1(:,8), g   => this%penta1(:,9),                         &
                   eobc => this%penta1(:,10)                                                  )
    
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
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            end select
            
            select case (bcn) 
            case(0)
                bt(this%n  ) = w1*zero
                b (this%n  ) = w1*alpha
                d (this%n  ) = w1*one
                a (this%n  ) = w1*zero
                at(this%n  ) = w1*zero

                bt(this%n-1) = w2*zero
                b (this%n-1) = w2*alpha_p
                d (this%n-1) = w2*one
                a (this%n-1) = w2*alpha_p
                at(this%n-1) = w2*zero
                
                bt(this%n-2) = w3*beta_pp
                b (this%n-2) = w3*alpha_pp
                d (this%n-2) = w3*one
                a (this%n-2) = w3*alpha_pp
                at(this%n-2) = w3*beta_pp

                bt(this%n-3) = w4*beta_ppp
                b (this%n-3) = w4*alpha_ppp
                d (this%n-3) = w4*one
                a (this%n-3) = w4*alpha_ppp
                at(this%n-3) = w4*beta_ppp
            
            case(1)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    
    subroutine ComputePenta2(this,bc1,bcn)
        class (cd10), intent(inout) :: this
        integer, intent(in) :: bc1, bcn

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if ((bc1 == 1) .and. (bcn == 1)) then
            this%penta2 = zero
        end if 

    end subroutine 


    subroutine SolveXLU1(this,y,n2,n3, y1, yn, z1, zn)
    
        class( cd10 ), intent(in) :: this
        integer, intent(in) :: n2,n3, y1, yn, z1, zn
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y  ! Take in RHS and put solution into it
        integer :: i,j,k
        real(rkind) :: sum1, sum2
 
        
        do k=z1,zn
            do j=y1,yn
                ! Step 8 ( update y instead of creating z )
                y(2,j,k) = y(2,j,k) - this%LU1(2,1)*y(1,j,k) 
                sum1 = this%LU1(1,3)*y(1,j,k) + this%LU1(2,3)*y(2,j,k)
                sum2 = this%LU1(1,4)*y(1,j,k) + this%LU1(2,4)*y(2,j,k)

                ! Step 9
                do i = 3,this%n-2
                    y(i,j,k) = y(i,j,k) - this%LU1(i,1)*y(i-1,j,k) - this%LU1(i,2)*y(i-2,j,k)
                    sum1 = sum1 + this%LU1(i,3)*y(i,j,k)
                    sum2 = sum2 + this%LU1(i,4)*y(i,j,k)
                end do
    
                ! Step 10
                y(this%n-1,j,k) = y(this%n-1,j,k) - sum1
                y(this%n,j,k)   = ( y(this%n,j,k)   - sum2 - this%LU1(this%n-1,4)*y(this%n-1,j,k) ) * this%LU1(this%n,5)
    
                ! Step 11
                y(this%n-1,j,k) = ( y(this%n-1,j,k) - this%LU1(this%n-1,9)*y(this%n,j,k) ) * this%LU1(this%n-1,5)
                y(this%n-2,j,k) = ( y(this%n-2,j,k) - this%LU1(this%n-2,8)*y(this%n-1,j,k) - this%LU1(this%n-2,9)*y(this%n,j,k) ) * this%LU1(this%n-2,5)
                y(this%n-3,j,k) = ( y(this%n-3,j,k) - this%LU1(this%n-3,6)*y(this%n-2,j,k) - this%LU1(this%n-3,8)*y(this%n-1,j,k) - this%LU1(this%n-3,9)*y(this%n,j,k) ) * this%LU1(this%n-3,5)
                do i = this%n-4,1,-1
                    y(i,j,k) = ( y(i,j,k) - this%LU1(i,6)*y(i+1,j,k) - this%LU1(i,7)*y(i+2,j,k) - this%LU1(i,8)*y(this%n-1,j,k) - this%LU1(i,9)*y(this%n,j,k) ) * this%LU1(i,5)
                end do
            end do
        end do
    
    end subroutine
    
    subroutine SolveYLU1(this,y,n1,n3, x1, xn, z1, zn)
    
        class( cd10 ), intent(in) :: this
        integer, intent(in) :: n1,n3, x1, xn, z1, zn
        real(rkind), dimension(n1,this%n,n3), intent(inout) :: y  ! Take in RHS and put solution into it
        integer :: j,k
        real(rkind), dimension(x1:xn) :: sum1, sum2
 
        
        do k=z1,zn
            ! Step 8 ( update y instead of creating z )
            y(x1:xn,2,k) = y(x1:xn,2,k) - this%LU1(2,1)*y(x1:xn,1,k) 
            sum1 = this%LU1(1,3)*y(x1:xn,1,k) + this%LU1(2,3)*y(x1:xn,2,k)
            sum2 = this%LU1(1,4)*y(x1:xn,1,k) + this%LU1(2,4)*y(x1:xn,2,k)

            ! Step 9
            do j = 3,this%n-2
                y(x1:xn,j,k) = y(x1:xn,j,k) - this%LU1(j,1)*y(x1:xn,j-1,k) - this%LU1(j,2)*y(x1:xn,j-2,k)
                sum1 = sum1 + this%LU1(j,3)*y(x1:xn,j,k)
                sum2 = sum2 + this%LU1(j,4)*y(x1:xn,j,k)
            end do
    
            ! Step 10
            y(x1:xn,this%n-1,k) = y(x1:xn,this%n-1,k) - sum1
            y(x1:xn,this%n,k)   = ( y(x1:xn,this%n,k)   - sum2 - this%LU1(this%n-1,4)*y(x1:xn,this%n-1,k) ) * this%LU1(this%n,5)
    
            ! Step 11
            y(x1:xn,this%n-1,k) = ( y(x1:xn,this%n-1,k) - this%LU1(this%n-1,9)*y(x1:xn,this%n,k) ) * this%LU1(this%n-1,5)
            y(x1:xn,this%n-2,k) = ( y(x1:xn,this%n-2,k) - this%LU1(this%n-2,8)*y(x1:xn,this%n-1,k) - this%LU1(this%n-2,9)*y(x1:xn,this%n,k) ) * this%LU1(this%n-2,5)
            y(x1:xn,this%n-3,k) = ( y(x1:xn,this%n-3,k) - this%LU1(this%n-3,6)*y(x1:xn,this%n-2,k) - this%LU1(this%n-3,8)*y(x1:xn,this%n-1,k) - this%LU1(this%n-3,9)*y(x1:xn,this%n,k) ) * this%LU1(this%n-3,5)
            do j = this%n-4,1,-1
                y(x1:xn,j,k) = ( y(x1:xn,j,k) - this%LU1(j,6)*y(x1:xn,j+1,k) - this%LU1(j,7)*y(x1:xn,j+2,k) - this%LU1(j,8)*y(x1:xn,this%n-1,k) - this%LU1(j,9)*y(x1:xn,this%n,k) ) * this%LU1(j,5)
            end do
        end do
    
    end subroutine
    
    subroutine SolveZLU1(this,y,n1,n2, x1, xn, y1, yn)
    
        class( cd10 ), intent(in) :: this
        integer, intent(in) :: n1,n2, x1, xn, y1, yn
        real(rkind), dimension(n1,n2,this%n), intent(inout) :: y  ! Take in RHS and put solution into it
        integer :: k
        real(rkind), dimension(x1:xn,y1:yn) :: sum1, sum2
 
        
        ! Step 8 ( update y instead of creating z )
        y(x1:xn,y1:yn,2) = y(x1:xn,y1:yn,2) - this%LU1(2,1)*y(x1:xn,y1:yn,1) 
        sum1 = this%LU1(1,3)*y(x1:xn,y1:yn,1) + this%LU1(2,3)*y(x1:xn,y1:yn,2)
        sum2 = this%LU1(1,4)*y(x1:xn,y1:yn,1) + this%LU1(2,4)*y(x1:xn,y1:yn,2)

        ! Step 9
        do k = 3,this%n-2
            y(x1:xn,y1:yn,k) = y(x1:xn,y1:yn,k) - this%LU1(k,1)*y(x1:xn,y1:yn,k-1) - this%LU1(k,2)*y(x1:xn,y1:yn,k-2)
            sum1 = sum1 + this%LU1(k,3)*y(x1:xn,y1:yn,k)
            sum2 = sum2 + this%LU1(k,4)*y(x1:xn,y1:yn,k)
        end do
    
        ! Step 10
        y(x1:xn,y1:yn,this%n-1) = y(x1:xn,y1:yn,this%n-1) - sum1
        y(x1:xn,y1:yn,this%n)   = ( y(x1:xn,y1:yn,this%n)   - sum2 - this%LU1(this%n-1,4)*y(x1:xn,y1:yn,this%n-1) ) * this%LU1(this%n,5)
    
        ! Step 11
        y(x1:xn,y1:yn,this%n-1) = ( y(x1:xn,y1:yn,this%n-1) - this%LU1(this%n-1,9)*y(x1:xn,y1:yn,this%n) ) * this%LU1(this%n-1,5)
        y(x1:xn,y1:yn,this%n-2) = ( y(x1:xn,y1:yn,this%n-2) - this%LU1(this%n-2,8)*y(x1:xn,y1:yn,this%n-1) - this%LU1(this%n-2,9)*y(x1:xn,y1:yn,this%n) ) * this%LU1(this%n-2,5)
        y(x1:xn,y1:yn,this%n-3) = ( y(x1:xn,y1:yn,this%n-3) - this%LU1(this%n-3,6)*y(x1:xn,y1:yn,this%n-2) - this%LU1(this%n-3,8)*y(x1:xn,y1:yn,this%n-1) - this%LU1(this%n-3,9)*y(x1:xn,y1:yn,this%n) ) * this%LU1(this%n-3,5)
        do k = this%n-4,1,-1
            y(x1:xn,y1:yn,k) = ( y(x1:xn,y1:yn,k) - this%LU1(k,6)*y(x1:xn,y1:yn,k+1) - this%LU1(k,7)*y(x1:xn,y1:yn,k+2) - this%LU1(k,8)*y(x1:xn,y1:yn,this%n-1) - this%LU1(k,9)*y(x1:xn,y1:yn,this%n) ) * this%LU1(k,5)
        end do
    
    end subroutine
    
    subroutine SolveXPenta1(this,y,n2,n3,y1,yn,z1,zn)

        class( cd10 ), intent(in) :: this
        integer, intent(in) :: n2,n3, y1,yn, z1, zn
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y
        integer :: i, j, k

        do k = z1,zn
            do j = y1,yn 
                ! Step 1
                y(2,j,k) = y(2,j,k) - this%penta1(2,8)*y(1,j,k)
                do i = 3,this%n
                    y(i,j,k) = y(i,j,k) - this%penta1(i,9)*y(i-2,j,k) - this%penta1(i,8)*y(i-1,j,k)
                end do 

                ! Step 2
                y(this%n,j,k) = y(this%n,j,k)*this%penta1(this%n,7)
                
                !y(this%n-1,j,k) = (y(this%n-1,j,k) - this%penta1(this%n-1,6)*y(this%n,j,k))*this%penta1(this%n-1,7)
                y(this%n-1,j,k) = y(this%n-1,j,k)*this%penta1(this%n-1,7) - this%penta1(this%n-1,10)*y(this%n,j,k)!*this%penta1(this%n-1,7)
                do i = this%n-2,1,-1
                    !y(i,j,k) = (y(i,j,k) - this%penta1(i,5)*y(i+2,j,k) - this%penta1(i,6)*y(i+1,j,k))*this%penta1(i,7)
                    y(i,j,k) = y(i,j,k)*this%penta1(i,7) - y(i+2,j,k)*this%penta1(i,5)*this%penta1(i,7) - y(i+1,j,k)*this%penta1(i,10)
                end do 
            end do 
        end do 

    end subroutine

    subroutine SolveYPenta1(this,y,n1,n3,x1,xn,z1,zn)

        class( cd10 ), intent(in) :: this
        integer, intent(in) :: n1,n3, x1,xn, z1, zn
        real(rkind), dimension(n1,this%n,n3), intent(inout) :: y
        integer :: j, k

        do k = z1,zn
            ! Step 1
            y(x1:xn,2,k) = y(x1:xn,2,k) - this%penta1(2,8)*y(x1:xn,1,k)
            do j = 3,this%n
                y(x1:xn,j,k) = y(x1:xn,j,k) - this%penta1(j,9)*y(x1:xn,j-2,k) - this%penta1(j,8)*y(x1:xn,j-1,k)
            end do 

            ! Step 2
            y(x1:xn,this%n,k) = y(x1:xn,this%n,k)*this%penta1(this%n,7)
            
            y(x1:xn,this%n-1,k) = y(x1:xn,this%n-1,k)*this%penta1(this%n-1,7) - this%penta1(this%n-1,10)*y(x1:xn,this%n,k)
            do j = this%n-2,1,-1
                y(x1:xn,j,k) = y(x1:xn,j,k)*this%penta1(j,7) - y(x1:xn,j+2,k)*this%penta1(j,5)*this%penta1(j,7) - y(x1:xn,j+1,k)*this%penta1(j,10)
            end do 
        end do 

    end subroutine

    subroutine SolveZPenta1(this,y,n1,n2,x1,xn,y1,yn)

        class( cd10 ), intent(in) :: this
        integer, intent(in) :: n1,n2,x1,xn,y1,yn
        real(rkind), dimension(n1,n2,this%n), intent(inout) :: y
        integer :: k

        ! Step 1
        y(x1:xn,y1:yn,2) = y(x1:xn,y1:yn,2) - this%penta1(2,8)*y(x1:xn,y1:yn,1)
        do k = 3,this%n
            y(x1:xn,y1:yn,k) = y(x1:xn,y1:yn,k) - this%penta1(k,9)*y(x1:xn,y1:yn,k-2) - this%penta1(k,8)*y(x1:xn,y1:yn,k-1)
        end do 

        ! Step 2
        y(x1:xn,y1:yn,this%n) = y(x1:xn,y1:yn,this%n)*this%penta1(this%n,7)
        
        y(x1:xn,y1:yn,this%n-1) = y(x1:xn,y1:yn,this%n-1)*this%penta1(this%n-1,7) - this%penta1(this%n-1,10)*y(x1:xn,y1:yn,this%n)
        do k = this%n-2,1,-1
            y(x1:xn,y1:yn,k) = y(x1:xn,y1:yn,k)*this%penta1(k,7) - y(x1:xn,y1:yn,k+2)*this%penta1(k,5)*this%penta1(k,7) - y(x1:xn,y1:yn,k+1)*this%penta1(k,10)
        end do 

    end subroutine

    subroutine SolveLU2(this,y)
         
        class( cd10 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(inout) :: y  ! Take in RHS and put solution into it
        integer :: i
        real(rkind) :: sum1, sum2
   
        ! Step 8 ( update y instead of creating z )
        y(2) = y(2) - this%LU2(2,1)*y(1) 
        sum1 = this%LU2(1,3)*y(1) + this%LU2(2,3)*y(2)
        sum2 = this%LU2(1,4)*y(1) + this%LU2(2,4)*y(2)

        ! Step 9
        do i = 3,this%n-2
            y(i) = y(i) - this%LU2(i,1)*y(i-1) - this%LU2(i,2)*y(i-2)
            sum1 = sum1 + this%LU2(i,3)*y(i)
            sum2 = sum2 + this%LU2(i,4)*y(i)
        end do
    
        ! Step 10
        y(this%n-1) = y(this%n-1) - sum1 !SUM( this%LU2(1:this%n-2,3)*y(1:this%n-2) )
        y(this%n)   = ( y(this%n)   - sum2 - this%LU2(this%n-1,4)*y(this%n-1) ) * this%LU2(this%n,5) !SUM( this%LU2(1:this%n-1,4)*y(1:this%n-1) )
    
        ! Step 11
        !y(this%n) = y(this%n) * this%LU2(this%n,5)
        y(this%n-1) = ( y(this%n-1) - this%LU2(this%n-1,9)*y(this%n) ) * this%LU2(this%n-1,5)
        y(this%n-2) = ( y(this%n-2) - this%LU2(this%n-2,8)*y(this%n-1) - this%LU2(this%n-2,9)*y(this%n) ) * this%LU2(this%n-2,5)
        y(this%n-3) = ( y(this%n-3) - this%LU2(this%n-3,6)*y(this%n-2) - this%LU2(this%n-3,8)*y(this%n-1) - this%LU2(this%n-3,9)*y(this%n) ) * this%LU2(this%n-3,5)
        do i = this%n-4,1,-1
            y(i) = ( y(i) - this%LU2(i,6)*y(i+1) - this%LU2(i,7)*y(i+2) - this%LU2(i,8)*y(this%n-1) - this%LU2(i,9)*y(this%n) ) * this%LU2(i,5)
        end do
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete/ Needs Updating 
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    end subroutine
    
    subroutine SolvePenta2(this,b)

        class( cd10 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(inout) :: b

        b = zero
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine

    pure subroutine ComputeXD1RHS(this, f, RHS, n2, n3, y1, yn, z1, zn)
    
        class( cd10 ), intent(in) :: this
        integer, intent(in) :: n2, n3, y1, yn, z1, zn
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        real(rkind) :: a10,b10,c10
        integer :: j,k
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4, b_np_4, c_np_4    
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1


    
        select case (this%periodic)
        case (.TRUE.)
            a10 = a10d1 * this%onebydx
            b10 = b10d1 * this%onebydx
            c10 = c10d1 * this%onebydx
            do k=z1,zn
                do j=y1,yn
                    RHS(1,j,k) = a10 * ( f(2,j,k)   - f(this%n  ,j,k) ) &
                           + b10 * ( f(3,j,k)   - f(this%n-1,j,k) ) &
                           + c10 * ( f(4,j,k)   - f(this%n-2,j,k) )
                    RHS(2,j,k) = a10 * ( f(3,j,k)   - f(1       ,j,k) ) &
                           + b10 * ( f(4,j,k)   - f(this%n  ,j,k) ) &
                           + c10 * ( f(5,j,k)   - f(this%n-1,j,k) )
                    RHS(3,j,k) = a10 * ( f(4,j,k)   - f(2       ,j,k) ) &
                           + b10 * ( f(5,j,k)   - f(1       ,j,k) ) &
                           + c10 * ( f(6,j,k)   - f(this%n  ,j,k) )
                    RHS(4:this%n-3,j,k) = a10 * ( f(5:this%n-2,j,k) - f(3:this%n-4,j,k) ) &
                                    + b10 * ( f(6:this%n-1,j,k) - f(2:this%n-5,j,k) ) &
                                    + c10 * ( f(7:this%n  ,j,k) - f(1:this%n-6,j,k) )
                    RHS(this%n-2,j,k) = a10 * ( f(this%n-1,j,k) - f(this%n-3,j,k) ) &
                                  + b10 * ( f(this%n  ,j,k) - f(this%n-4,j,k) ) &
                                  + c10 * ( f(1       ,j,k) - f(this%n-5,j,k) )
                    RHS(this%n-1,j,k) = a10 * ( f(this%n  ,j,k) - f(this%n-2,j,k) ) &
                                  + b10 * ( f(1       ,j,k) - f(this%n-3,j,k) ) &
                                  + c10 * ( f(2       ,j,k) - f(this%n-4,j,k) )
                    RHS(this%n  ,j,k) = a10 * ( f(1       ,j,k) - f(this%n-1,j,k) ) &
                                  + b10 * ( f(2       ,j,k) - f(this%n-2,j,k) ) &
                                  + c10 * ( f(3       ,j,k) - f(this%n-3,j,k) )
                end do
            end do
        case (.FALSE.)
            a10    = q_hat * this%onebydx  
            b10    = r_hat * this%onebydx 
            c10    = s_hat * this%onebydx

            a_np_4 = w4*q_ppp * this%onebydx  
            b_np_4 = w4*r_ppp * this%onebydx 
            c_np_4 = w4*s_ppp * this%onebydx
            
            a_np_3 = w3*q_pp * this%onebydx  
            b_np_3 = w3*r_pp * this%onebydx 

            a_np_2 = w2*q_p * this%onebydx
            
            a_np_1 = w1*( p * this%onebydx)
            b_np_1 = w1*( q * this%onebydx)
            c_np_1 = w1*( r * this%onebydx)
            d_np_1 = w1*( s * this%onebydx)

            do k = z1,zn
                do j = y1,yn
                    
                    RHS(1         ,j,k) =   a_np_1* f(1         ,j,k) +  b_np_1*f(2         ,j,k)   &
                                        +   c_np_1* f(3         ,j,k) +  d_np_1*f(4         ,j,k) 
                   
                    RHS(2         ,j,k) =   a_np_2*(f(3         ,j,k) -         f(1         ,j,k))
                    
                    RHS(3         ,j,k) =   a_np_3*(f(4         ,j,k) -         f(2         ,j,k)) &
                                        +   b_np_3*(f(5         ,j,k) -         f(1         ,j,k)) 
                    
                    RHS(4         ,j,k) =   a_np_4*(f(5         ,j,k) -         f(3         ,j,k)) &
                                        +   b_np_4*(f(6         ,j,k) -         f(2         ,j,k)) &
                                        +   c_np_4*(f(7         ,j,k) -         f(1         ,j,k)) 
                    
                    RHS(5:this%n-4,j,k) =   a10   *(f(6:this%n-3,j,k) -         f(4:this%n-5,j,k)) &
                                        +   b10   *(f(7:this%n-2,j,k) -         f(3:this%n-6,j,k)) &
                                        +   c10   *(f(8:this%n-1,j,k) -         f(2:this%n-7,j,k)) 
                    
                    RHS(this%n-3  ,j,k) =   a_np_4*(f(this%n-2  ,j,k) -         f(this%n-4  ,j,k)) &
                                        +   b_np_4*(f(this%n-1  ,j,k) -         f(this%n-5  ,j,k)) &
                                        +   c_np_4*(f(this%n    ,j,k) -         f(this%n-6  ,j,k))  
        
                    RHS(this%n-2  ,j,k) =   a_np_3*(f(this%n-1  ,j,k) -         f(this%n-3  ,j,k)) &
                                        +   b_np_3*(f(this%n    ,j,k) -         f(this%n-4  ,j,k)) 
                    
                    RHS(this%n-1  ,j,k) =   a_np_2*(f(this%n    ,j,k) -         f(this%n-2  ,j,k))

                    RHS(this%n    ,j,k) =  -a_np_1* f(this%n    ,j,k) -  b_np_1*f(this%n-1  ,j,k)   &
                                        -   c_np_1* f(this%n-2  ,j,k) -  d_np_1*f(this%n-3  ,j,k)
                
               end do 
            end do 
        end select
    
    end subroutine
    
    pure subroutine ComputeYD1RHS(this, f, RHS, n1, n3, x1, xn, z1, zn) 
    
        class( cd10 ), intent(in) :: this
        integer, intent(in) :: n1, n3, x1, xn, z1, zn
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: RHS
        real(rkind) :: a10,b10,c10
        integer :: k
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4, b_np_4, c_np_4    
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1


    
        select case (this%periodic)
        case (.TRUE.)
            a10 = a10d1 * this%onebydx
            b10 = b10d1 * this%onebydx
            c10 = c10d1 * this%onebydx
            do k=z1,zn
                !do j=y1,yn
                    RHS(x1:xn,1,k) = a10 * ( f(x1:xn,2,k)   - f(x1:xn,this%n  ,k) ) &
                           + b10 * ( f(x1:xn,3,k)   - f(x1:xn,this%n-1,k) ) &
                           + c10 * ( f(x1:xn,4,k)   - f(x1:xn,this%n-2,k) )
                    RHS(x1:xn,2,k) = a10 * ( f(x1:xn,3,k)   - f(x1:xn,1       ,k) ) &
                           + b10 * ( f(x1:xn,4,k)   - f(x1:xn,this%n  ,k) ) &
                           + c10 * ( f(x1:xn,5,k)   - f(x1:xn,this%n-1,k) )
                    RHS(x1:xn,3,k) = a10 * ( f(x1:xn,4,k)   - f(x1:xn,2       ,k) ) &
                           + b10 * ( f(x1:xn,5,k)   - f(x1:xn,1       ,k) ) &
                           + c10 * ( f(x1:xn,6,k)   - f(x1:xn,this%n  ,k) )
                    RHS(x1:xn,4:this%n-3,k) = a10 * ( f(x1:xn,5:this%n-2,k) - f(x1:xn,3:this%n-4,k) ) &
                                    + b10 * ( f(x1:xn,6:this%n-1,k) - f(x1:xn,2:this%n-5,k) ) &
                                    + c10 * ( f(x1:xn,7:this%n  ,k) - f(x1:xn,1:this%n-6,k) )
                    RHS(x1:xn,this%n-2,k) = a10 * ( f(x1:xn,this%n-1,k) - f(x1:xn,this%n-3,k) ) &
                                  + b10 * ( f(x1:xn,this%n  ,k) - f(x1:xn,this%n-4,k) ) &
                                  + c10 * ( f(x1:xn,1       ,k) - f(x1:xn,this%n-5,k) )
                    RHS(x1:xn,this%n-1,k) = a10 * ( f(x1:xn,this%n  ,k) - f(x1:xn,this%n-2,k) ) &
                                  + b10 * ( f(x1:xn,1       ,k) - f(x1:xn,this%n-3,k) ) &
                                  + c10 * ( f(x1:xn,2       ,k) - f(x1:xn,this%n-4,k) )
                    RHS(x1:xn,this%n  ,k) = a10 * ( f(x1:xn,1       ,k) - f(x1:xn,this%n-1,k) ) &
                                  + b10 * ( f(x1:xn,2       ,k) - f(x1:xn,this%n-2,k) ) &
                                  + c10 * ( f(x1:xn,3       ,k) - f(x1:xn,this%n-3,k) )
                !end do
            end do
        case (.FALSE.)
            a10    = q_hat * this%onebydx  
            b10    = r_hat * this%onebydx 
            c10    = s_hat * this%onebydx

            a_np_4 = w4*q_ppp * this%onebydx  
            b_np_4 = w4*r_ppp * this%onebydx 
            c_np_4 = w4*s_ppp * this%onebydx
            
            a_np_3 = w3*q_pp * this%onebydx  
            b_np_3 = w3*r_pp * this%onebydx 

            a_np_2 = w2*q_p * this%onebydx
            
            a_np_1 = w1*( p * this%onebydx)
            b_np_1 = w1*( q * this%onebydx)
            c_np_1 = w1*( r * this%onebydx)
            d_np_1 = w1*( s * this%onebydx)

            do k = z1,zn
                !do j = y1,yn
                    
                    RHS(x1:xn,1         ,k) =   a_np_1* f(x1:xn,1         ,k) +  b_np_1*f(x1:xn,2         ,k)   &
                                        +   c_np_1* f(x1:xn,3         ,k) +  d_np_1*f(x1:xn,4         ,k) 
                   
                    RHS(x1:xn,2         ,k) =   a_np_2*(f(x1:xn,3         ,k) -         f(x1:xn,1         ,k))
                    
                    RHS(x1:xn,3         ,k) =   a_np_3*(f(x1:xn,4         ,k) -         f(x1:xn,2         ,k)) &
                                        +   b_np_3*(f(x1:xn,5         ,k) -         f(x1:xn,1         ,k)) 
                    
                    RHS(x1:xn,4         ,k) =   a_np_4*(f(x1:xn,5         ,k) -         f(x1:xn,3         ,k)) &
                                        +   b_np_4*(f(x1:xn,6         ,k) -         f(x1:xn,2         ,k)) &
                                        +   c_np_4*(f(x1:xn,7         ,k) -         f(x1:xn,1         ,k)) 
                    
                    RHS(x1:xn,5:this%n-4,k) =   a10   *(f(x1:xn,6:this%n-3,k) -         f(x1:xn,4:this%n-5,k)) &
                                        +   b10   *(f(x1:xn,7:this%n-2,k) -         f(x1:xn,3:this%n-6,k)) &
                                        +   c10   *(f(x1:xn,8:this%n-1,k) -         f(x1:xn,2:this%n-7,k)) 
                    
                    RHS(x1:xn,this%n-3  ,k) =   a_np_4*(f(x1:xn,this%n-2  ,k) -         f(x1:xn,this%n-4  ,k)) &
                                        +   b_np_4*(f(x1:xn,this%n-1  ,k) -         f(x1:xn,this%n-5  ,k)) &
                                        +   c_np_4*(f(x1:xn,this%n    ,k) -         f(x1:xn,this%n-6  ,k))  
        
                    RHS(x1:xn,this%n-2  ,k) =   a_np_3*(f(x1:xn,this%n-1  ,k) -         f(x1:xn,this%n-3  ,k)) &
                                        +   b_np_3*(f(x1:xn,this%n    ,k) -         f(x1:xn,this%n-4  ,k)) 
                    
                    RHS(x1:xn,this%n-1  ,k) =   a_np_2*(f(x1:xn,this%n    ,k) -         f(x1:xn,this%n-2  ,k))

                    RHS(x1:xn,this%n    ,k) =  -a_np_1* f(x1:xn,this%n    ,k) -  b_np_1*f(x1:xn,this%n-1  ,k)   &
                                        -   c_np_1* f(x1:xn,this%n-2  ,k) -  d_np_1*f(x1:xn,this%n-3  ,k)
                
               !end do 
            end do 
        end select
    
    end subroutine

    pure subroutine ComputeZD1RHS(this, f, RHS, n1, n2, x1, xn, y1, yn)
    
        class( cd10 ), intent(in) :: this
        integer, intent(in) :: n1, n2, x1, xn, y1, yn
        real(rkind), dimension(n1,n2,this%n), intent(in) :: f
        real(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
        real(rkind) :: a10,b10,c10
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4, b_np_4, c_np_4    
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1


    
        select case (this%periodic)
        case (.TRUE.)
            a10 = a10d1 * this%onebydx
            b10 = b10d1 * this%onebydx
            c10 = c10d1 * this%onebydx
            RHS(x1:xn,y1:yn,1) = a10 * ( f(x1:xn,y1:yn,2)   - f(x1:xn,y1:yn,this%n  ) ) &
                               + b10 * ( f(x1:xn,y1:yn,3)   - f(x1:xn,y1:yn,this%n-1) ) &
                               + c10 * ( f(x1:xn,y1:yn,4)   - f(x1:xn,y1:yn,this%n-2) )
            RHS(x1:xn,y1:yn,2) = a10 * ( f(x1:xn,y1:yn,3)   - f(x1:xn,y1:yn,1       ) ) &
                               + b10 * ( f(x1:xn,y1:yn,4)   - f(x1:xn,y1:yn,this%n  ) ) &
                               + c10 * ( f(x1:xn,y1:yn,5)   - f(x1:xn,y1:yn,this%n-1) )
            RHS(x1:xn,y1:yn,3) = a10 * ( f(x1:xn,y1:yn,4)   - f(x1:xn,y1:yn,2       ) ) &
                               + b10 * ( f(x1:xn,y1:yn,5)   - f(x1:xn,y1:yn,1       ) ) &
                               + c10 * ( f(x1:xn,y1:yn,6)   - f(x1:xn,y1:yn,this%n  ) )
            RHS(x1:xn,y1:yn,4:this%n-3) = a10 * ( f(x1:xn,y1:yn,5:this%n-2) - f(x1:xn,y1:yn,3:this%n-4) ) &
                                        + b10 * ( f(x1:xn,y1:yn,6:this%n-1) - f(x1:xn,y1:yn,2:this%n-5) ) &
                                        + c10 * ( f(x1:xn,y1:yn,7:this%n  ) - f(x1:xn,y1:yn,1:this%n-6) )
            RHS(x1:xn,y1:yn,this%n-2) = a10 * ( f(x1:xn,y1:yn,this%n-1) - f(x1:xn,y1:yn,this%n-3) ) &
                                      + b10 * ( f(x1:xn,y1:yn,this%n  ) - f(x1:xn,y1:yn,this%n-4) ) &
                                      + c10 * ( f(x1:xn,y1:yn,1       ) - f(x1:xn,y1:yn,this%n-5) )
            RHS(x1:xn,y1:yn,this%n-1) = a10 * ( f(x1:xn,y1:yn,this%n  ) - f(x1:xn,y1:yn,this%n-2) ) &
                                      + b10 * ( f(x1:xn,y1:yn,1       ) - f(x1:xn,y1:yn,this%n-3) ) &
                                      + c10 * ( f(x1:xn,y1:yn,2       ) - f(x1:xn,y1:yn,this%n-4) )
            RHS(x1:xn,y1:yn,this%n  ) = a10 * ( f(x1:xn,y1:yn,1       ) - f(x1:xn,y1:yn,this%n-1) ) &
                                              + b10 * ( f(x1:xn,y1:yn,2       ) - f(x1:xn,y1:yn,this%n-2) ) &
                                              + c10 * ( f(x1:xn,y1:yn,3       ) - f(x1:xn,y1:yn,this%n-3) )
        case (.FALSE.)
            a10    = q_hat * this%onebydx  
            b10    = r_hat * this%onebydx 
            c10    = s_hat * this%onebydx

            a_np_4 = w4*q_ppp * this%onebydx  
            b_np_4 = w4*r_ppp * this%onebydx 
            c_np_4 = w4*s_ppp * this%onebydx
            
            a_np_3 = w3*q_pp * this%onebydx  
            b_np_3 = w3*r_pp * this%onebydx 

            a_np_2 = w2*q_p * this%onebydx
            
            a_np_1 = w1*( p * this%onebydx)
            b_np_1 = w1*( q * this%onebydx)
            c_np_1 = w1*( r * this%onebydx)
            d_np_1 = w1*( s * this%onebydx)

                    
            RHS(x1:xn,y1:yn,1         ) =   a_np_1* f(x1:xn,y1:yn,1         ) +  b_np_1*f(x1:xn,y1:yn,2         )   &
                                +   c_np_1* f(x1:xn,y1:yn,3         ) +  d_np_1*f(x1:xn,y1:yn,4         ) 
            
            RHS(x1:xn,y1:yn,2         ) =   a_np_2*(f(x1:xn,y1:yn,3         ) -         f(x1:xn,y1:yn,1         ))
            
            RHS(x1:xn,y1:yn,3         ) =   a_np_3*(f(x1:xn,y1:yn,4         ) -         f(x1:xn,y1:yn,2         )) &
                                +   b_np_3*(f(x1:xn,y1:yn,5         ) -         f(x1:xn,y1:yn,1         )) 
            
            RHS(x1:xn,y1:yn,4         ) =   a_np_4*(f(x1:xn,y1:yn,5         ) -         f(x1:xn,y1:yn,3         )) &
                                +   b_np_4*(f(x1:xn,y1:yn,6         ) -         f(x1:xn,y1:yn,2         )) &
                                +   c_np_4*(f(x1:xn,y1:yn,7         ) -         f(x1:xn,y1:yn,1         )) 
            
            RHS(x1:xn,y1:yn,5:this%n-4) =   a10   *(f(x1:xn,y1:yn,6:this%n-3) -         f(x1:xn,y1:yn,4:this%n-5)) &
                                +   b10   *(f(x1:xn,y1:yn,7:this%n-2) -         f(x1:xn,y1:yn,3:this%n-6)) &
                                +   c10   *(f(x1:xn,y1:yn,8:this%n-1) -         f(x1:xn,y1:yn,2:this%n-7)) 
            
            RHS(x1:xn,y1:yn,this%n-3  ) =   a_np_4*(f(x1:xn,y1:yn,this%n-2  ) -         f(x1:xn,y1:yn,this%n-4  )) &
                                +   b_np_4*(f(x1:xn,y1:yn,this%n-1  ) -         f(x1:xn,y1:yn,this%n-5  )) &
                                +   c_np_4*(f(x1:xn,y1:yn,this%n    ) -         f(x1:xn,y1:yn,this%n-6  ))  
        
            RHS(x1:xn,y1:yn,this%n-2  ) =   a_np_3*(f(x1:xn,y1:yn,this%n-1  ) -         f(x1:xn,y1:yn,this%n-3  )) &
                                +   b_np_3*(f(x1:xn,y1:yn,this%n    ) -         f(x1:xn,y1:yn,this%n-4  )) 
            
            RHS(x1:xn,y1:yn,this%n-1  ) =   a_np_2*(f(x1:xn,y1:yn,this%n    ) -         f(x1:xn,y1:yn,this%n-2  ))

            RHS(x1:xn,y1:yn,this%n    ) =  -a_np_1* f(x1:xn,y1:yn,this%n    ) -  b_np_1*f(x1:xn,y1:yn,this%n-1  )   &
                                -   c_np_1* f(x1:xn,y1:yn,this%n-2  ) -  d_np_1*f(x1:xn,y1:yn,this%n-3  )
            
        end select
    
    end subroutine
    
    pure function ComputeD2RHS(this,f) result (RHS)
    
        class( cd10 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: f
        real(rkind), dimension(this%n) :: RHS
        integer :: i
    
        select case (this%periodic)
        case (.TRUE.)
            RHS(1) = a10d1 * ( f(2)   - two*f(1) + f(this%n)   ) &
                   + b10d1 * ( f(3)   - two*f(1) + f(this%n-1) ) &
                   + c10d1 * ( f(4)   - two*f(1) + f(this%n-2) )
            RHS(2) = a10d1 * ( f(3)   - two*f(2) + f(1)        ) &
                   + b10d1 * ( f(4)   - two*f(2) + f(this%n)   ) &
                   + c10d1 * ( f(5)   - two*f(2) + f(this%n-1) )
            RHS(3) = a10d1 * ( f(4)   - two*f(3) + f(2)        ) &
                   + b10d1 * ( f(5)   - two*f(3) + f(1)        ) &
                   + c10d1 * ( f(6)   - two*f(3) + f(this%n)   )
            do i = 4,this%n-3
                RHS(i) = a10d1 * ( f(i+1) - two*f(i) + f(i-1) ) &
                       + b10d1 * ( f(i+2) - two*f(i) + f(i-2) ) &
                       + c10d1 * ( f(i+3) - two*f(i) + f(i-3) )
            end do
            RHS(this%n-2) = a10d1 * ( f(this%n-1) - two*f(this%n-2) + f(this%n-3) ) &
                          + b10d1 * ( f(this%n)   - two*f(this%n-2) + f(this%n-4) ) &
                          + c10d1 * ( f(1)        - two*f(this%n-2) + f(this%n-5) )
            RHS(this%n-1) = a10d1 * ( f(this%n)   - two*f(this%n-1) + f(this%n-2) ) &
                          + b10d1 * ( f(1)        - two*f(this%n-1) + f(this%n-3) ) &
                          + c10d1 * ( f(2)        - two*f(this%n-1) + f(this%n-4) )
            RHS(this%n)   = a10d1 * ( f(1)        - two*f(this%n)   + f(this%n-1) ) &
                          + b10d1 * ( f(2)        - two*f(this%n)   + f(this%n-2) ) &
                          + c10d1 * ( f(3)        - two*f(this%n)   + f(this%n-3) )
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete/ Needs to be changed 
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case (.FALSE.)
            RHS = zero
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end select
    
    end function

    subroutine dd1(this, f, df, na, nb, a1, an, b1, bn)
        class( cd10 ), intent(in) :: this
        integer, intent(in) :: na, nb, a1, an, b1, bn
        real(rkind), dimension(this%n,na,nb), intent(in) :: f
        real(rkind), dimension(this%n,na,nb), intent(out) :: df

        if(this%n == 1) then
            df = zero
            return
        end if
        
        call this%ComputeXD1RHS(f, df, na, nb, a1, an, b1, bn)

        select case (this%periodic)
        case(.TRUE.)
            call this%SolveXLU1(df, na, nb, a1, an, b1, bn)
        case(.FALSE.)
            call this%SolveXPenta1(df, na, nb, a1, an, b1, bn)
        end select
    
    end subroutine

    subroutine dd2(this, f, df, na, nb, a1, an, b1, bn)
        class( cd10 ), intent(in) :: this
        integer, intent(in) :: na, nb, a1, an, b1, bn
        real(rkind), dimension(na,this%n,nb), intent(in) :: f
        real(rkind), dimension(na,this%n,nb), intent(out) :: df

        if(this%n == 1) then
            df = zero
            return
        end if
        
        call this%ComputeYD1RHS(f, df, na, nb, a1, an, b1, bn)

        select case (this%periodic)
        case(.TRUE.)
            call this%SolveYLU1(df, na, nb, a1, an, b1, bn)
        case(.FALSE.)
            call this%SolveYPenta1(df, na, nb, a1, an, b1, bn)
        end select
    
    end subroutine

    subroutine dd3(this, f, df, na, nb, a1, an, b1, bn)
        class( cd10 ), intent(in) :: this
        integer, intent(in) :: na, nb, a1, an, b1, bn
        real(rkind), dimension(na,nb,this%n), intent(in) :: f
        real(rkind), dimension(na,nb,this%n), intent(out) :: df

        if(this%n == 1) then
            df = zero
            return
        end if
        
        call this%ComputeZD1RHS(f, df, na, nb, a1, an, b1, bn)

        select case (this%periodic)
        case(.TRUE.)
            call this%SolveZLU1(df, na, nb, a1, an, b1, bn)
        case(.FALSE.)
            call this%SolveZPenta1(df, na, nb, a1, an, b1, bn)
        end select
    
    end subroutine

    function cd10der2(this, f) result(df)

        class( cd10 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: f
        real(rkind), dimension(this%n) :: df


        print*, f(1)
        df = zero

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end function

end module
