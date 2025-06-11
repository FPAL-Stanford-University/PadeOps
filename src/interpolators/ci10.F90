! Routines specific to 10th order Compact Finite Differencing scheme
! Periodic LU based on  Neossi Nguetchue, Abelman (Appl. Math. & Comp. 2008)

module ci10stuff

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two
    use exits,           only: GracefulExit
    
    implicit none

    private
    public :: ci10, alpha10d1, beta10d1, a10d1, b10d1, c10d1
    ! 10th order first derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha10d1= 2.0_rkind/5.0_rkind  !10.0_rkind /  21.0_rkind
    real(rkind), parameter :: beta10d1 = (14_rkind*alpha10d1-5_rkind)/42_rkind !5.0_rkind / 126.0_rkind
    real(rkind), parameter :: a10d1    = (1.0_rkind / 8.0_rkind) * (10._rkind + 7._rkind*alpha10d1 )/2.0_rkind  !( 5.0_rkind / 3.0_rkind) / 2.0_rkind
    real(rkind), parameter :: b10d1    = 1._rkind/112._rkind*(189._rkind*alpha10d1 - 50._rkind ) /2.0_rkind !(5.0_rkind /14.0_rkind) / 2.0_rkind
    real(rkind), parameter :: c10d1    = (5._rkind*alpha10d1 -2._rkind )/ 48._rkind/2.0_rkind !(  1.0_rkind /126.0_rkind) / 2.0_rkind

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
    real(rkind), parameter                   :: q_hat           = ( 5.0_rkind / 3.0_rkind) !a10d1
    real(rkind), parameter                   :: r_hat           = (5.0_rkind /14.0_rkind) ! b10d1
    real(rkind), parameter                   :: s_hat           = (  1.0_rkind /126.0_rkind)  ! c10d1
    real(rkind), parameter                   :: alpha_hat       = 10.0_rkind / 21.0_rkind! alpha10d1
    real(rkind), parameter                   :: beta_hat        = (14_rkind*alpha10d1-5_rkind)/42_rkind  !beta10d1
     
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
   
    ! 1st point  

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type ci10
        
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
        real(rkind), allocatable, dimension(:,:) :: penta1_nn
        real(rkind), allocatable, dimension(:,:) :: penta1_ns
        real(rkind), allocatable, dimension(:,:) :: penta1_na
        real(rkind), allocatable, dimension(:,:) :: penta1_sn
        real(rkind), allocatable, dimension(:,:) :: penta1_ss
        real(rkind), allocatable, dimension(:,:) :: penta1_sa
        real(rkind), allocatable, dimension(:,:) :: penta1_an
        real(rkind), allocatable, dimension(:,:) :: penta1_as
        real(rkind), allocatable, dimension(:,:) :: penta1_aa
        

        contains

        procedure :: init
        procedure :: destroy
        procedure :: GetSize
        procedure, private :: ComputeXD1RHS
        procedure, private :: ComputeYD1RHS
        procedure, private :: ComputeZD1RHS
        

        procedure, private :: SolveXLU1
        procedure, private :: SolveYLU1
        procedure, private :: SolveZLU1
        
        procedure, private :: ComputePenta1

        procedure, private :: SolveXPenta1
        procedure, private :: SolveYPenta1
        procedure, private :: SolveZPenta1
        
        
        procedure ::iN2F1
        procedure ::iN2F2
        procedure ::iN2F3
        

    end type



contains
    
    pure function GetSize(this) result(val)
        class(ci10), intent(in) :: this
        integer  :: val 
        val = this%n
    end function

    function init(this, n_, dx_, periodic_, bc1_, bcn_) result(ierr)
   
        class( ci10 ), intent(inout) :: this
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
            if(allocated( this%LU1 )) deallocate( this%LU1 ); allocate( this%LU1(n_,9) ); this%LU1 = zero
    
            ! Compute 1st derivative LU matrix
            if (n_ .GE. 8) then
                call ComputeLU(this%LU1,n_,beta10d1,alpha10d1,one,alpha10d1,beta10d1)
            else if (n_ == 1) then
                this%LU1 = one
            else
                ierr = 2
                return
            end if
    
        else 
            ! Allocate 1st derivative Penta matrices.
            if(allocated( this%penta1_nn )) deallocate( this%penta1_nn ); allocate( this%penta1_nn(n_,11) ); this%penta1_nn = zero
            if(allocated( this%penta1_ns )) deallocate( this%penta1_ns ); allocate( this%penta1_ns(n_,11) ); this%penta1_ns = zero
            if(allocated( this%penta1_na )) deallocate( this%penta1_na ); allocate( this%penta1_na(n_,11) ); this%penta1_na = zero
            if(allocated( this%penta1_sn )) deallocate( this%penta1_sn ); allocate( this%penta1_sn(n_,11) ); this%penta1_sn = zero
            if(allocated( this%penta1_ss )) deallocate( this%penta1_ss ); allocate( this%penta1_ss(n_,11) ); this%penta1_ss = zero
            if(allocated( this%penta1_sa )) deallocate( this%penta1_sa ); allocate( this%penta1_sa(n_,11) ); this%penta1_sa = zero
            if(allocated( this%penta1_an )) deallocate( this%penta1_an ); allocate( this%penta1_an(n_,11) ); this%penta1_an = zero
            if(allocated( this%penta1_as )) deallocate( this%penta1_as ); allocate( this%penta1_as(n_,11) ); this%penta1_as = zero
            if(allocated( this%penta1_aa )) deallocate( this%penta1_aa ); allocate( this%penta1_aa(n_,11) ); this%penta1_aa = zero
  
            if (n_ .GE. 8) then             
                call this%ComputePenta1(this%penta1_nn, 0, 0)    ! Standard
                call this%ComputePenta1(this%penta1_ns, 0, 1)    ! Standard-Symm
                call this%ComputePenta1(this%penta1_na, 0,-1)    ! Standard-Asym
                call this%ComputePenta1(this%penta1_sn, 1, 0)    ! Symm-Standard
                call this%ComputePenta1(this%penta1_ss, 1, 1)    ! Symm-Symm
                call this%ComputePenta1(this%penta1_sa, 1,-1)    ! Symm-Asym
                call this%ComputePenta1(this%penta1_an,-1, 0)    ! Asym-Standard
                call this%ComputePenta1(this%penta1_as,-1, 1)    ! Asym-Symm
                call this%ComputePenta1(this%penta1_aa,-1,-1)    ! Asym-Asym
            else if (n_ .EQ. 1) then
                this%penta1_nn = one
                this%penta1_ns = one
                this%penta1_na = one
                this%penta1_sn = one
                this%penta1_ss = one
                this%penta1_sa = one
                this%penta1_an = one
                this%penta1_as = one
                this%penta1_aa = one
            else
                ierr = 2
                return 
            end if 

        end if 

        ! If everything passes
        ierr = 0
    
    end function
    
    subroutine destroy(this)

        class( ci10 ), intent(inout) :: this

        ! Dellocate 1st derivative LU matrix.
        if(allocated( this%LU1 )) deallocate( this%LU1 )
    
        ! Dellocate 2nd derivative LU matrix.
        if(allocated( this%LU2 )) deallocate( this%LU2 )
    
        ! Dellocate 1st derivative penta matrix.
        if(allocated( this%penta1_nn )) deallocate( this%penta1_nn )
        if(allocated( this%penta1_ns )) deallocate( this%penta1_ns )
        if(allocated( this%penta1_na )) deallocate( this%penta1_na )
        if(allocated( this%penta1_sn )) deallocate( this%penta1_sn )
        if(allocated( this%penta1_ss )) deallocate( this%penta1_ss )
        if(allocated( this%penta1_sa )) deallocate( this%penta1_sa )
        if(allocated( this%penta1_an )) deallocate( this%penta1_an )
        if(allocated( this%penta1_as )) deallocate( this%penta1_as )
        if(allocated( this%penta1_aa )) deallocate( this%penta1_aa )
    
    
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

    subroutine ComputePenta1(this,penta1,bc1,bcn)
        class(ci10), intent(inout) :: this
        real(rkind), dimension(this%n,11), intent(inout) :: penta1
        integer, intent(in) :: bc1, bcn
        integer             :: i
    
        associate (bt   => penta1(:,1), b   => penta1(:,2), d => penta1(:,3),  &
                   a    => penta1(:,4), at  => penta1(:,5),                    &
                   e    => penta1(:,6), obc => penta1(:,7),                    &
                   f    => penta1(:,8), g   => penta1(:,9),                    &
                   eobc => penta1(:,10)                                        )
    
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

                bt(1) = zero
                b (1) = zero
                d (1) = one
                a (1) = zero
                at(1) = zero

                bt(2) = zero
                b (2) = alpha_hat
                d (2) = one - beta_hat
                a (2) = alpha_hat
                at(2) = beta_hat

            case(-1)

                bt(1) = zero
                b (1) = zero
                d (1) = one
                a (1) = two*alpha_hat
                at(1) = two*beta_hat

                bt(2) = zero
                b (2) = alpha_hat
                d (2) = one + beta_hat
                a (2) = alpha_hat
                at(2) = beta_hat

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

                bt(this%n  ) = zero
                b (this%n  ) = zero
                d (this%n  ) = one
                a (this%n  ) = zero
                at(this%n  ) = zero

                bt(this%n-1) = beta_hat
                b (this%n-1) = alpha_hat
                d (this%n-1) = one - beta_hat
                a (this%n-1) = alpha_hat
                at(this%n-1) = zero

            case(-1)

                bt(this%n  ) = two*beta_hat
                b (this%n  ) = two*alpha_hat
                d (this%n  ) = one
                a (this%n  ) = zero
                at(this%n  ) = zero
                           
                bt(this%n-1) = beta_hat
                b (this%n-1) = alpha_hat
                d (this%n-1) = one + beta_hat
                a (this%n-1) = alpha_hat
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
    

    subroutine SolveXLU1(this,y,n2,n3)
    
        class( ci10 ), intent(in) :: this
        integer, intent(in) :: n2,n3
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y  ! Take in RHS and put solution into it
        integer :: i,j,k
        real(rkind) :: sum1, sum2
 
        
        do k=1,n3
            do j=1,n2
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
    
    subroutine SolveYLU1(this,y,n1,n3)
    
        class( ci10 ), intent(in) :: this
        integer, intent(in) :: n1,n3
        real(rkind), dimension(n1,this%n,n3), intent(inout) :: y  ! Take in RHS and put solution into it
        integer :: j,k
        real(rkind), dimension(n1) :: sum1, sum2
 
        
        do k=1,n3
            ! Step 8 ( update y instead of creating z )
            y(:,2,k) = y(:,2,k) - this%LU1(2,1)*y(:,1,k) 
            sum1 = this%LU1(1,3)*y(:,1,k) + this%LU1(2,3)*y(:,2,k)
            sum2 = this%LU1(1,4)*y(:,1,k) + this%LU1(2,4)*y(:,2,k)

            ! Step 9
            do j = 3,this%n-2
                y(:,j,k) = y(:,j,k) - this%LU1(j,1)*y(:,j-1,k) - this%LU1(j,2)*y(:,j-2,k)
                sum1 = sum1 + this%LU1(j,3)*y(:,j,k)
                sum2 = sum2 + this%LU1(j,4)*y(:,j,k)
            end do
    
            ! Step 10
            y(:,this%n-1,k) = y(:,this%n-1,k) - sum1
            y(:,this%n,k)   = ( y(:,this%n,k)   - sum2 - this%LU1(this%n-1,4)*y(:,this%n-1,k) ) * this%LU1(this%n,5)
    
            ! Step 11
            y(:,this%n-1,k) = ( y(:,this%n-1,k) - this%LU1(this%n-1,9)*y(:,this%n,k) ) * this%LU1(this%n-1,5)
            y(:,this%n-2,k) = ( y(:,this%n-2,k) - this%LU1(this%n-2,8)*y(:,this%n-1,k) - this%LU1(this%n-2,9)*y(:,this%n,k) ) * this%LU1(this%n-2,5)
            y(:,this%n-3,k) = ( y(:,this%n-3,k) - this%LU1(this%n-3,6)*y(:,this%n-2,k) - this%LU1(this%n-3,8)*y(:,this%n-1,k) - this%LU1(this%n-3,9)*y(:,this%n,k) ) * this%LU1(this%n-3,5)
            do j = this%n-4,1,-1
                y(:,j,k) = ( y(:,j,k) - this%LU1(j,6)*y(:,j+1,k) - this%LU1(j,7)*y(:,j+2,k) - this%LU1(j,8)*y(:,this%n-1,k) - this%LU1(j,9)*y(:,this%n,k) ) * this%LU1(j,5)
            end do
        end do
    
    end subroutine
    
    subroutine SolveZLU1(this,y,n1,n2)
    
        class( ci10 ), intent(in) :: this
        integer, intent(in) :: n1,n2
        real(rkind), dimension(n1,n2,this%n), intent(inout) :: y  ! Take in RHS and put solution into it
        ! integer :: k
        integer :: j, k
        ! real(rkind), dimension(n1,n2) :: sum1, sum2
        real(rkind), dimension(n1) :: sum1, sum2
 
        
        ! ! Step 8 ( update y instead of creating z )
        ! y(:,:,2) = y(:,:,2) - this%LU1(2,1)*y(:,:,1) 
        ! sum1 = this%LU1(1,3)*y(:,:,1) + this%LU1(2,3)*y(:,:,2)
        ! sum2 = this%LU1(1,4)*y(:,:,1) + this%LU1(2,4)*y(:,:,2)

        ! ! Step 9
        ! do k = 3,this%n-2
        !     y(:,:,k) = y(:,:,k) - this%LU1(k,1)*y(:,:,k-1) - this%LU1(k,2)*y(:,:,k-2)
        !     sum1 = sum1 + this%LU1(k,3)*y(:,:,k)
        !     sum2 = sum2 + this%LU1(k,4)*y(:,:,k)
        ! end do
    
        ! ! Step 10
        ! y(:,:,this%n-1) = y(:,:,this%n-1) - sum1
        ! y(:,:,this%n)   = ( y(:,:,this%n)   - sum2 - this%LU1(this%n-1,4)*y(:,:,this%n-1) ) * this%LU1(this%n,5)
    
        ! ! Step 11
        ! y(:,:,this%n-1) = ( y(:,:,this%n-1) - this%LU1(this%n-1,9)*y(:,:,this%n) ) * this%LU1(this%n-1,5)
        ! y(:,:,this%n-2) = ( y(:,:,this%n-2) - this%LU1(this%n-2,8)*y(:,:,this%n-1) - this%LU1(this%n-2,9)*y(:,:,this%n) ) * this%LU1(this%n-2,5)
        ! y(:,:,this%n-3) = ( y(:,:,this%n-3) - this%LU1(this%n-3,6)*y(:,:,this%n-2) - this%LU1(this%n-3,8)*y(:,:,this%n-1) - this%LU1(this%n-3,9)*y(:,:,this%n) ) * this%LU1(this%n-3,5)
        ! do k = this%n-4,1,-1
        !     y(:,:,k) = ( y(:,:,k) - this%LU1(k,6)*y(:,:,k+1) - this%LU1(k,7)*y(:,:,k+2) - this%LU1(k,8)*y(:,:,this%n-1) - this%LU1(k,9)*y(:,:,this%n) ) * this%LU1(k,5)
        ! end do
    
        do j=1,n2
            ! Step 8 ( update y instead of creating z )
            y(:,j,2) = y(:,j,2) - this%LU1(2,1)*y(:,j,1) 
            sum1 = this%LU1(1,3)*y(:,j,1) + this%LU1(2,3)*y(:,j,2)
            sum2 = this%LU1(1,4)*y(:,j,1) + this%LU1(2,4)*y(:,j,2)

            ! Step 9
            do k = 3,this%n-2
                y(:,j,k) = y(:,j,k) - this%LU1(k,1)*y(:,j,k-1) - this%LU1(k,2)*y(:,j,k-2)
                sum1 = sum1 + this%LU1(k,3)*y(:,j,k)
                sum2 = sum2 + this%LU1(k,4)*y(:,j,k)
            end do
    
            ! Step 10
            y(:,j,this%n-1) = y(:,j,this%n-1) - sum1
            y(:,j,this%n)   = ( y(:,j,this%n)   - sum2 - this%LU1(this%n-1,4)*y(:,j,this%n-1) ) * this%LU1(this%n,5)
    
            ! Step 11
            y(:,j,this%n-1) = ( y(:,j,this%n-1) - this%LU1(this%n-1,9)*y(:,j,this%n) ) * this%LU1(this%n-1,5)
            y(:,j,this%n-2) = ( y(:,j,this%n-2) - this%LU1(this%n-2,8)*y(:,j,this%n-1) - this%LU1(this%n-2,9)*y(:,j,this%n) ) * this%LU1(this%n-2,5)
            y(:,j,this%n-3) = ( y(:,j,this%n-3) - this%LU1(this%n-3,6)*y(:,j,this%n-2) - this%LU1(this%n-3,8)*y(:,j,this%n-1) - this%LU1(this%n-3,9)*y(:,j,this%n) ) * this%LU1(this%n-3,5)
            do k = this%n-4,1,-1
                y(:,j,k) = ( y(:,j,k) - this%LU1(k,6)*y(:,j,k+1) - this%LU1(k,7)*y(:,j,k+2) - this%LU1(k,8)*y(:,j,this%n-1) - this%LU1(k,9)*y(:,j,this%n) ) * this%LU1(k,5)
            end do
        end do
    
    end subroutine
    
    subroutine SolveXPenta1(this,penta1,y,n2,n3)

        class( ci10 ), intent(in) :: this
        real(rkind), dimension(this%n,11), intent(in) :: penta1
        integer, intent(in) :: n2,n3
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y
        integer :: i, j, k

        do k = 1,n3
            do j = 1,n2
                ! Step 1
                y(2,j,k) = y(2,j,k) - penta1(2,8)*y(1,j,k)
                do i = 3,this%n
                    y(i,j,k) = y(i,j,k) - penta1(i,9)*y(i-2,j,k) - penta1(i,8)*y(i-1,j,k)
                end do 

                ! Step 2
                y(this%n,j,k) = y(this%n,j,k)*penta1(this%n,7)
                
                !y(this%n-1,j,k) = (y(this%n-1,j,k) - penta1(this%n-1,6)*y(this%n,j,k))*penta1(this%n-1,7)
                y(this%n-1,j,k) = y(this%n-1,j,k)*penta1(this%n-1,7) - penta1(this%n-1,10)*y(this%n,j,k)!*penta1(this%n-1,7)
                do i = this%n-2,1,-1
                    !y(i,j,k) = (y(i,j,k) - penta1(i,5)*y(i+2,j,k) - penta1(i,6)*y(i+1,j,k))*penta1(i,7)
                    y(i,j,k) = y(i,j,k)*penta1(i,7) - y(i+2,j,k)*penta1(i,5)*penta1(i,7) - y(i+1,j,k)*penta1(i,10)
                end do 
            end do 
        end do 

    end subroutine

    subroutine SolveYPenta1(this,penta1,y,n1,n3)

        class( ci10 ), intent(in) :: this
        real(rkind), dimension(this%n,11), intent(in) :: penta1
        integer, intent(in) :: n1,n3
        real(rkind), dimension(n1,this%n,n3), intent(inout) :: y
        integer :: j, k

        do k = 1,n3
            ! Step 1
            y(:,2,k) = y(:,2,k) - penta1(2,8)*y(:,1,k)
            do j = 3,this%n
                y(:,j,k) = y(:,j,k) - penta1(j,9)*y(:,j-2,k) - penta1(j,8)*y(:,j-1,k)
            end do 

            ! Step 2
            y(:,this%n,k) = y(:,this%n,k)*penta1(this%n,7)
            
            y(:,this%n-1,k) = y(:,this%n-1,k)*penta1(this%n-1,7) - penta1(this%n-1,10)*y(:,this%n,k)
            do j = this%n-2,1,-1
                y(:,j,k) = y(:,j,k)*penta1(j,7) - y(:,j+2,k)*penta1(j,5)*penta1(j,7) - y(:,j+1,k)*penta1(j,10)
            end do 
        end do 

    end subroutine

    subroutine SolveZPenta1(this,penta1,y,n1,n2)

        class( ci10 ), intent(in) :: this
        real(rkind), dimension(this%n,11), intent(in) :: penta1
        integer, intent(in) :: n1,n2
        real(rkind), dimension(n1,n2,this%n), intent(inout) :: y
        integer :: k

        ! Step 1
        y(:,:,2) = y(:,:,2) - penta1(2,8)*y(:,:,1)
        do k = 3,this%n
            y(:,:,k) = y(:,:,k) - penta1(k,9)*y(:,:,k-2) - penta1(k,8)*y(:,:,k-1)
        end do 

        ! Step 2
        y(:,:,this%n) = y(:,:,this%n)*penta1(this%n,7)
        
        y(:,:,this%n-1) = y(:,:,this%n-1)*penta1(this%n-1,7) - penta1(this%n-1,10)*y(:,:,this%n)
        do k = this%n-2,1,-1
            y(:,:,k) = y(:,:,k)*penta1(k,7) - y(:,:,k+2)*penta1(k,5)*penta1(k,7) - y(:,:,k+1)*penta1(k,10)
        end do 

    end subroutine


    pure subroutine ComputeXD1RHS(this, f, RHS,dir, n2, n3, bc1, bcn)
    
        class( ci10 ), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10,b10,c10
        integer :: j,k
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4, b_np_4, c_np_4    
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1

        select case (this%periodic)
        case (.TRUE.)
            a10 = a10d1 
            b10 = b10d1
            c10 = c10d1 
            select case(dir)
              case("N2F")

               do k = 1,n3
                  do j = 1,n2
                     RHS(1         ,j,k) = a10 * ( f(2,j,k)         + f(1,j,k) ) &
                                        + b10 * ( f(3,j,k)          + f(this%n,j,k) ) &
                                        + c10 *( f(4,j,k)           + f(this%n-1,j,k) )

                    RHS(2         ,j,k) = a10 * ( f(3,j,k)          + f(2 ,j,k) ) &
                                        + b10 * ( f(4,j,k)          + f(1 ,j,k)) &
                                        + c10 *( f(5,j,k)           + f(this%n,j,k) )

                    RHS(3:this%n-3,j,k) = a10 * ( f(4:this%n-2,j,k) + f(3:this%n-3,j,k) ) &
                                        + b10 * ( f(5:this%n-1  ,j,k) + f(2:this%n-4,j,k) ) &
                                        + c10 * ( f(6:this%n, j, k )  + f(1:this%n-5,j,k) )

                    RHS(this%n-2  ,j,k) = a10 * ( f(this%n-1,j,k)     + f(this%n-2,j,k) ) &
                                        + b10 * ( f(this%n,j,k)       + f(this%n-3,j,k) ) &
                                        + c10 * ( f(1,j,k)            + f(this%n-4,j,k) )

                    RHS(this%n-1  ,j,k) = a10 * ( f(this%n,j,k)       + f(this%n-1,j,k) ) &
                                        + b10 * ( f(1,j,k)            + f(this%n-2,j,k) ) &
                                        + c10 * ( f(2,j,k)            + f(this%n-3,j,k) )

                    RHS(this%n    ,j,k) = a10 * ( f(1,j,k)          + f(this%n,j,k) ) &
                                        + b10 * ( f(2,j,k)          + f(this%n-1,j,k) ) &
                                        + c10 * ( f(3,j,k)          + f(this%n-2,j,k) )
                end do
              end do
          case("F2N")

           !!!!!!!!!!!!!! TO DO
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          end select
  
        case (.FALSE.)
        end select
    
    end subroutine
    
    pure subroutine ComputeYD1RHS(this, f, RHS,dir, n1, n3, bc1, bcn) 
    
        class( ci10 ), intent(in) :: this
        integer, intent(in) :: n1, n3
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10,b10,c10
        integer :: k
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4, b_np_4, c_np_4    
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1


        select case (this%periodic)
        case (.TRUE.)
            a10 = a10d1 
            b10 = b10d1 
            c10 = c10d1 

            select case(dir)
               case("N2F")
                 do k = 1,n3


                     RHS(:,1         ,k) = a10 * ( f(:,2,k)          + f(:,1,k) ) & 
                                        + b10 * ( f(:,3,k)           + f(:,this%n,k) ) &
                                        + c10 *( f(:,4,k)            + f(:,this%n-1,k) )

                     RHS(:,2         ,k) = a10 * ( f(:,3,k)          + f(:,2 ,k) ) &
                                            + b10 * ( f(:,4,k)       + f(:,1,k) ) &
                                            + c10 *( f(:,5,k)        + f(:,this%n,k) )

                     RHS(:,3:this%n-3,k) = a10 * ( f(:,4:this%n-2,k)   + f(:,3:this%n-3,k) ) &
                                           + b10 * ( f(:,5:this%n-1,k) + f(:,2:this%n-4,k) ) &
                                           + c10 * ( f(:,6:this%n, k ) + f(:,1:this%n-5,k) )

                     RHS(:,this%n-2,k) = a10 * ( f(:,this%n-1,k)     + f(:,this%n-2,k) ) &
                                         + b10 * ( f(:,this%n,k)     + f(:,this%n-3,k) ) &
                                         + c10 * ( f(:,1,k)          + f(:,this%n-4,k) )

                     RHS(:,this%n-1,k) = a10 * ( f(:,this%n,k)       + f(:,this%n-1,k) ) &
                                          + b10 * ( f(:,1,k)         + f(:,this%n-2,k) ) &
                                          + c10 * ( f(:,2,k)         + f(:,this%n-3,k) )

                     RHS(:,this%n    ,k) = a10 * ( f(:,1,k)          + f(:,this%n,k) ) &
                                            + b10 * ( f(:,2,k)       + f(:,this%n-1,k) ) &
                                            + c10 * ( f(:,3,k)       + f(:,this%n-2,k) )


                !end do
            end do
            case("F2N")

            !!!!!!!!!!!!!!!!!!!!!!!!!!   TO DO
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         end select

        case (.FALSE.)
        end select
    
    end subroutine

    !pure subroutine ComputeZD1RHS(this, f, RHS, n1, n2, bc1, bcn)
    subroutine ComputeZD1RHS(this, f, RHS,dir, n1, n2, bc1, bcn)
    
        class( ci10 ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in) :: f
        real(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10,b10,c10
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4, b_np_4, c_np_4    
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1
        integer :: j,i,k


        select case (this%periodic)
        case (.TRUE.)
            a10 = a10d1 
            b10 = b10d1 
            c10 = c10d1 
          select case(dir)
            case("N2F")

              
             RHS(:,:         ,1) = a10 * ( f(:,:,2)          + f(:,:,1) ) &
                                   + b10 * ( f(:,:,3)        + f(:,:,this%n) ) &
                                   + c10 *( f(:,:,4)         +f(:,:,this%n-1) )

             RHS(:,:         ,2) = a10 * ( f(:,:,3)          + f(:,:,2) ) &
                                   + b10 * ( f(:,:,4)        + f(:,:,1) ) &
                                   + c10 *( f(:,:,5)         + f(:,:,this%n) )

             RHS(:,:,3:this%n-3) = a10 * ( f(:,:,4:this%n-2) + f(:,:,3:this%n-3) ) &
                                   + b10 * ( f(:,:,5:this%n-1) + f(:,:,2:this%n-4) ) &
                                   + c10 * ( f(:,:,6:this%n )  + f(:,:,1:this%n-5) )

             RHS(:,:,this%n-2) = a10 * ( f(:,:,this%n-1)     + f(:,:,this%n-2) )&
                                  + b10 * ( f(:,:,this%n)    + f(:,:,this%n-3) )&
                                  + c10 * ( f(:,:,1)         + f(:,:,this%n-4) )

             RHS(:,:,this%n-1) = a10 * ( f(:,:,this%n)       + f(:,:,this%n-1) ) &
                                 + b10 * ( f(:,:,1)         + f(:,:,this%n-2) ) &
                                 + c10 * ( f(:,:,2)         + f(:,:,this%n-3) )

             RHS(:,:,this%n    ) = a10 * ( f(:,:,1)          + f(:,:,this%n) ) &
                                   + b10 * ( f(:,:,2)       + f(:,:,this%n-1) ) &
                                   + c10 * ( f(:,:,3)       + f(:,:,this%n-2) )

            case("F2N")


           !!!!!!!!!!!!!! TO DO
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        end select

        case (.FALSE.)
            
        end select
    
    end subroutine
    
    subroutine iN2F1(this, fN, fF, na, nb, bc1_, bcn_)
        class( ci10 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(this%n,na,nb), intent(in)  :: fN
        real(rkind), dimension(this%n,na,nb), intent(out) :: fF

        call this%ComputeXD1RHS(fN, fF,"N2F", na, nb, bc1, bcn)

        select case (this%periodic)
        case(.TRUE.)
            call this%SolveXLU1(fF, na, nb)!   
        case(.FALSE.)
!            select case(bc1)
!            case(0) ! Normal non-periodic left boundary
!                select case(bcn)
!                case(0)  ! Normal non-periodic right boundary
!                    call this%SolveXPenta1(this%penta1_nn, df, na, nb)
!                case(1) ! Symmetric right boundary
!                    call this%SolveXPenta1(this%penta1_ns, df, na, nb)
!                case(-1) ! Antisymmetric right boundary
!                    call this%SolveXPenta1(this%penta1_na, df, na, nb)
!                end select
!            case(1) ! Symmetric left !boundary
!                select case(bcn)
!                case(0)  ! Normal non-periodic right boundary
!                    call this%SolveXPenta1(this%penta1_sn, df, na, nb)
!                case(1)  ! Symmetric right boundary
!                    call this%SolveXPenta1(this%penta1_ss, df, na, nb)
!                case(-1) ! Antisymmetric right boundary
!                    call this%SolveXPenta1(this%penta1_sa, df, na, nb)
!                end select
!            case(-1) ! Antisymmetric left boundary
!                select case(bcn)
!                case(0)  ! Normal non-periodic right boundary
!                    call this%SolveXPenta1(this%penta1_an, df, na, nb)
!                case(1) ! Symmetric right boundary
!                    call this%SolveXPenta1(this%penta1_as, df, na, nb)
!                case(-1) ! Antisymmetric right boundary
!                    call this%SolveXPenta1(this%penta1_aa, df, na, nb)
!                end select
!            end select
        end select
        
    end subroutine

    subroutine iN2F2(this, fN, fF, na, nb, bc1_, bcn_)
        class( ci10 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(na,this%n,nb), intent(in) :: fN
        real(rkind), dimension(na,this%n,nb), intent(out) :: fF

        if(this%n == 1) then
            fF = zero
            return
        end if
        

        call this%ComputeYD1RHS(fN, fF,"N2F", na, nb, bc1, bcn)

        select case (this%periodic)
        case(.TRUE.)
            call this%SolveYLU1(fF, na, nb)
        case(.FALSE.)
!           select case(bc1)
!           case(0) ! Normal non-periodic left boundary
!               select case(bcn)
!               case(0)  ! Normal non-periodic right boundary
!                   call this%SolveYPenta1(this%penta1_nn, df, na, nb)
!               case(1) ! Symmetric right boundary
!                   call this%SolveYPenta1(this%penta1_ns, df, na, nb)
!               case(-1) ! Antisymmetric right boundary
!                   call this%SolveYPenta1(this%penta1_na, df, na, nb)
!               end select
!           case(1) ! Symmetric left boundary
!               select case(bcn)
!               case(0)  ! Normal non-periodic right boundary
!                   call this%SolveYPenta1(this%penta1_sn, df, na, nb)
!               case(1)  ! Symmetric right boundary
!                   call this%SolveYPenta1(this%penta1_ss, df, na, nb)
!               case(-1) ! Antisymmetric right boundary
!                   call this%SolveYPenta1(this%penta1_sa, df, na, nb)
!               end select
!           case(-1) ! Antisymmetric left boundary
!               select case(bcn)
!               case(0)  ! Normal non-periodic right boundary
!                   call this%SolveYPenta1(this%penta1_an, df, na, nb)
!               case(1) ! Symmetric right boundary
!                   call this%SolveYPenta1(this%penta1_as, df, na, nb)
!               case(-1) ! Antisymmetric right boundary
!                   call this%SolveYPenta1(this%penta1_aa, df, na, nb)
!               end select
!           end select
        end select
    
    end subroutine

    subroutine iN2F3(this, fN, fF, na, nb, bc1_, bcn_)
        class( ci10 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(na,nb,this%n), intent(in) :: fN
        real(rkind), dimension(na,nb,this%n), intent(out) :: fF

        if(this%n == 1) then
            fF = zero
            return
        end if


        call this%ComputeZD1RHS(fN, fF,"N2F", na, nb, bc1, bcn)

        select case (this%periodic)
        case(.TRUE.)
            call this%SolveZLU1(fF, na, nb)
        case(.FALSE.)
!            select case(bc1)
!            case(0) ! Normal non-periodic left boundary
!                select case(bcn)
!                case(0)  ! Normal non-periodic right boundary
!                    call this%SolveZPenta1(this%penta1_nn, df, na, nb)
!                case(1) ! Symmetric right boundary
!                    call this%SolveZPenta1(this%penta1_ns, df, na, nb)
!                case(-1) ! Antisymmetric right boundary
!                    call this%SolveZPenta1(this%penta1_na, df, na, nb)
!                end select
!            case(1) ! Symmetric left boundary
!                select case(bcn)
!                case(0)  ! Normal non-periodic right boundary
!                    call this%SolveZPenta1(this%penta1_sn, df, na, nb)
!                case(1)  ! Symmetric right boundary
!                    call this%SolveZPenta1(this%penta1_ss, df, na, nb)
!                case(-1) ! Antisymmetric right boundary
!                    call this%SolveZPenta1(this%penta1_sa, df, na, nb)
!                end select
!            case(-1) ! Antisymmetric left boundary
!                select case(bcn)
!                case(0)  ! Normal non-periodic right boundary
!                    call this%SolveZPenta1(this%penta1_an, df, na, nb)
!                case(1) ! Symmetric right boundary
!                    call this%SolveZPenta1(this%penta1_as, df, na, nb)
!                case(-1) ! Antisymmetric right boundary
!                    call this%SolveZPenta1(this%penta1_aa, df, na, nb)
!                end select
!            end select
        end select
    
    end subroutine


end module
