! Routines specific to 6th order Compact Finite Differencing scheme
! Periodic LU based on  Neossi Nguetchue, Abelman (Appl. Math. & Comp. 2008)

module ci06stuff

    use kind_parameters, only: rkind
    use constants,       only : zero,one,two,half
    use exits,    only: gracefulExit, message

    implicit none

    private
    public :: ci06, alpha06d1, a06d1, b06d1
    
    ! 6th order first derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha06d1=  3.0_rkind / 10.0_rkind
    real(rkind), parameter :: a06d1    = (3.0_rkind / 2.0_rkind) / 2.0_rkind
    real(rkind), parameter :: b06d1    = ( 1.0_rkind / 10.0_rkind) / 2.0_rkind
     ! 2nd order first derivative explicit centeral difference coefficients
    real(rkind), parameter :: aI02     =  150.0/256.0 !75.0d0/64.0d0
    real(rkind), parameter :: bI02     =  -25.0/256.0!-25.0d0/128.0d0
    real(rkind), parameter :: cI02     =  3.0/256.0 ! 3.0d0/128.0d0

    !real(rkind),  parameter :: a06d1    = (10.0*alpha06d1 + 9.0) / 8.0 / 2.0;
    !const Real rhs_coef_b = ( 6.0*alpha - 1.0) / 8.0 / 2.0;
    !!!!! TODO nonperiodic implementation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic 1st derivative evaluation !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
    ! Set the scheme for the edge nodes (Ref. for notation: Lele - JCP paper)
    ! 1st derivative 
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
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic 2nd derivative evaluation !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    ! Incomplete
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    type ci06

        private

        integer     :: n
        real(rkind) :: dx
        real(rkind) :: onebydx
        real(rkind) :: onebydx2

        logical     :: periodic = .TRUE.
        integer     :: bc1 = 0                             ! Boundary condition type. 0=Dirichlet, 1=Neumann
        integer     :: bcn = 0                             ! Boundary condition type. 0=Dirichlet, 1=Neumann

        real(rkind), allocatable, dimension(:,:) :: LU1
        real(rkind), allocatable, dimension(:,:) :: LU2
        real(rkind), allocatable, dimension(:,:) :: LU1_BC
        real(rkind), allocatable, dimension(:,:) :: Tri1
        real(rkind), allocatable, dimension(:,:) :: Tri2

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
        procedure, private :: SolveXLU1_BC
        procedure, private :: SolveYLU1_BC
        procedure, private :: SolveZLU1_BC
        
       
        procedure, private :: ComputeTri1

        procedure, private :: SolveXTri1
        procedure, private :: SolveYTri1
        procedure, private :: SolveZTri1
        
        
        procedure :: iN2F1
        procedure :: iN2F2
        procedure :: iN2F3

        procedure :: iF2N1
        procedure :: iF2N2
        procedure :: iF2N3


    end type

contains

    pure function GetSize(this) result(val)
        class(ci06), intent(in) :: this
        integer  :: val 
        val = this%n
    end function
    
    function init(this, n_, dx_, periodic_, bc1_, bcn_) result(ierr)
    
        class( ci06 ), intent(inout) :: this
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
            if(allocated( this%LU1 )) deallocate( this%LU1 ); allocate( this%LU1(n_,5) )
    
            ! Compute 1st derivative LU matrix
            if (n_ .GE. 6) then
                call ComputeLU(this%LU1,n_,alpha06d1,one,alpha06d1)
            else if (n_ == 1) then
                this%LU1 = one 
            else
                ierr = 3
                return
            end if
    
        
        else
            ! Allocate 1st derivative Tri matrices.
            if(allocated( this%tri1 )) deallocate( this%tri1 ); allocate( this%tri1(n_,3) )
  
            if (n_ .GE. 8) then             
            call ComputeLU(this%LU1_BC,n_-3,alpha06d1,one,alpha06d1)
            !   call this%ComputeTri1(bc1_,bcn_)
            else if (n_ .EQ. 1) then
                this%Tri1 = one
            else
                ierr = 2
                return 
            end if 

              
        end if 

    
        ! If everything passes
        ierr = 0
    
    end function

    subroutine destroy(this)

        class( ci06 ), intent(inout) :: this

        ! Dellocate 1st derivative LU matrix.
        if(allocated( this%LU1 )) deallocate( this%LU1 )
    
        ! Dellocate 1st derivative LU matrix.
        if(allocated( this%LU1_BC )) deallocate( this%LU1_BC )

        ! Dellocate 1st derivative tri matrix.
        if(allocated( this%tri1 )) deallocate( this%tri1 )

    end subroutine
    
    subroutine ComputeLU(LU,n,b,d,a)
    
        integer, intent(in) :: n
        real(rkind), intent(in) :: d,a,b
        real(rkind), dimension(n,5), intent(out) :: LU
        integer :: i
    
        LU = 0.0_rkind
    
        associate ( bc=>LU(:,1), h=>LU(:,2), c=>LU(:,3), aa=>LU(:,4), v=>LU(:,5) )
    
            ! Step 0
            c(1) = d
            v(1) = b
            h(1) = a/c(1)
    
            ! Step 1
            do i = 2,n-1
                c(i) = d - ( b/c(i-1) )*a
            end do
            do i = 2,n-2
                v(i) = -( b/c(i-1) )*v(i-1)
                h(i) = -( a/c(i) )*h(i-1)
            end do
            v(n-1) = a - ( b/c(n-2) )*v(n-2)
            h(n-1) = ( b - h(n-2)*a ) / c(n-1)
            c(n)   = d - SUM( h(1:n-1)*v(1:n-1) )
    
            bc(2:n-1) = b / c(1:n-2)
            aa(1:n-2) = a
    
            ! Set c = 1/c
            c = 1._rkind/c
   
            ! Overwrite aa by aa*c
            aa = aa*c
            
            ! Overwrite v by v*c
            v = v*c 
        end associate
    
    end subroutine
    
    subroutine ComputeTri1(this,bc1,bcn)
        class (ci06), intent(inout) :: this
        integer, intent(in) :: bc1, bcn
        integer             :: i
        real(rkind), dimension(this%n) :: a, b, c, cp, den

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
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            ! Incomplete
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end select

        select case (bcn) 
        case(0)
            c (this%n  ) = w1*zero
            c (this%n-1) = w2*alpha_p
            c (this%n-2) = w3*alpha_pp 

            b (this%n  ) = w1*one
            b (this%n-1) = w2*one
            b (this%n-2) = w3*one 

            a (this%n  ) = w1*alpha
            a (this%n-1) = w2*alpha_p
            a (this%n-2) = w3*alpha_pp 
        case(1)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            ! Incomplete
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end select

        cp(1) = c(1)/b(1)
        do i = 2,this%n-1
            cp(i) = c(i)/(b(i) - a(i)*cp(i-1))
        end do
        
        den(1) = one/b(1)
        den(2:this%n) = one/(b(2:this%n) - a(2:this%n)*cp(1:this%n-1))

        this%Tri1(:,1) = a*den
        this%Tri1(:,2) = den
        this%Tri1(:,3) = cp
           
    end subroutine
    
    subroutine SolveXLU1_BC(this,y,n2,n3)
        class (ci06), intent(in) :: this
        integer, intent(in) :: n2,n3
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y  ! Take inRHSand put solution into it
        integer :: i, j, k
        real(rkind) :: sum1

        do k = 1,n3
            do j = 1,n2
                ! Step 2
                sum1 = this%LU1_BC(1,2)*y(2,j,k)
                do i = 3,this%n-3
                    y(i,j,k) = y(i,j,k) - this%LU1_BC(i-1,1)*y(i-1,j,k)
                    sum1 = sum1 + this%LU1_BC(i-1,2)*y(i,j,k)
                end do
                y(this%n-2,j,k) = y(this%n-2,j,k) - sum1

                ! Step 3
                y(this%n-2,j,k)   = y(this%n-2,j,k) * this%LU1_BC(this%n-3,3)

                y(this%n-3,j,k) =  y(this%n-3,j,k) * this%LU1_BC(this%n-4,3)-y(this%n-2,j,k) * this%LU1_BC(this%n-4,5)
                do i = this%n-4,2,-1
                    y(i,j,k) =  y(i,j,k) * this%LU1_BC(i-1,3)- y(i+1,j,k)*this%LU1_BC(i-1,4)- y(this%n-2,j,k) * this%LU1_BC(i-1,5)
                end do
            end do
        end do

    end subroutine

    
    subroutine SolveXLU1(this,y,n2,n3)
        
        class (ci06), intent(in) :: this
        integer, intent(in) :: n2,n3
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y  ! Take in RHS and put solution into it
        integer :: i, j, k
        real(rkind) :: sum1 

        do k = 1,n3
            do j = 1,n2
                ! Step 2
                sum1 = this%LU1(1,2)*y(1,j,k)
                do i = 2,this%n-1
                    y(i,j,k) = y(i,j,k) - this%LU1(i,1)*y(i-1,j,k)
                    sum1 = sum1 + this%LU1(i,2)*y(i,j,k)
                end do
                y(this%n,j,k) = y(this%n,j,k) - sum1
    
                ! Step 3
                y(this%n,j,k)   = y(this%n,j,k) * this%LU1(this%n,3)

                y(this%n-1,j,k) =  y(this%n-1,j,k) * this%LU1(this%n-1,3) - y(this%n,j,k) * this%LU1(this%n-1,5) 
                do i = this%n-2,1,-1
                    y(i,j,k) =  y(i,j,k) * this%LU1(i,3)- y(i+1,j,k) * this%LU1(i,4)- y(this%n,j,k) * this%LU1(i,5)
                end do
            end do 
        end do 
    
    end subroutine

    subroutine SolveYLU1_BC(this,y,n1,n3)

        class (ci06), intent(in) :: this
        integer, intent(in) :: n1,n3
        real(rkind), dimension(n1,this%n,n3), intent(inout) :: y  ! Take inRHSand put solution into it
        integer ::  j, k
        real(rkind), dimension(n1) :: sum1

        do k = 1,n3
                ! Step 2
                sum1 = this%LU1_BC(1,2)*y(:,2,k)
                do j = 3,this%n-3
                    y(:,j,k) = y(:,j,k) - this%LU1_BC(j-1,1)*y(:,j-1,k)
                    sum1 = sum1 + this%LU1_BC(j-1,2)*y(:,j,k)
                end do
                y(:,this%n-2,k) = y(:,this%n-2,k) - sum1

                ! Step 3
                y(:,this%n-2,k)   = y(:,this%n-2,k) * this%LU1_BC(this%n-3,3)

                y(:,this%n-3,k) =  y(:,this%n-3,k) * this%LU1_BC(this%n-4,3)-y(:,this%n-2,k) * this%LU1_BC(this%n-4,5)
                do j = this%n-4,2,-1
                    y(:,j,k) =  y(:,j,k) * this%LU1_BC(j-1,3)- y(:,j+1,k)*this%LU1_BC(j-1,4)- y(:,this%n-2,k) * this%LU1_BC(j-1,5)
                end do
        end do

    end subroutine



    
    subroutine SolveYLU1(this,y,n1,n3)
        
        class (ci06), intent(in) :: this
        integer, intent(in) :: n1,n3
        real(rkind), dimension(n1,this%n,n3), intent(inout) :: y  ! Take in RHS and put solution into it
        integer ::  j, k
        real(rkind), dimension(n1) :: sum1 

        do k = 1,n3
                ! Step 2
                sum1 = this%LU1(1,2)*y(:,1,k)
                do j = 2,this%n-1
                    y(:,j,k) = y(:,j,k) - this%LU1(j,1)*y(:,j-1,k)
                    sum1 = sum1 + this%LU1(j,2)*y(:,j,k)
                end do
                y(:,this%n,k) = y(:,this%n,k) - sum1
    
                ! Step 3
                y(:,this%n,k)   = y(:,this%n,k) * this%LU1(this%n,3)

                y(:,this%n-1,k) =  y(:,this%n-1,k) * this%LU1(this%n-1,3) - y(:,this%n,k) * this%LU1(this%n-1,5) 
                do j = this%n-2,1,-1
                    y(:,j,k) =  y(:,j,k) * this%LU1(j,3)- y(:,j+1,k) * this%LU1(j,4)- y(:,this%n,k) * this%LU1(j,5)
                end do
        end do 
    
    end subroutine

    
    subroutine SolveZLU1(this,y,n1,n2)
        
        class (ci06), intent(in) :: this
        integer, intent(in) :: n1,n2
        real(rkind), dimension(n1,n2,this%n), intent(inout) :: y  ! Take in RHS and put solution into it
        integer ::  k
        real(rkind), dimension(n1,n2) :: sum1 

        ! Step 2
        sum1 = this%LU1(1,2)*y(:,:,1)
        do k = 2,this%n-1
            y(:,:,k) = y(:,:,k) - this%LU1(k,1)*y(:,:,k-1)
            sum1 = sum1 + this%LU1(k,2)*y(:,:,k)
        end do
        y(:,:,this%n) = y(:,:,this%n) - sum1
    
        ! Step 3
        y(:,:,this%n)   = y(:,:,this%n) * this%LU1(this%n,3)

        y(:,:,this%n-1) =  y(:,:,this%n-1) * this%LU1(this%n-1,3) - y(:,:,this%n) * this%LU1(this%n-1,5) 
        do k = this%n-2,1,-1
            y(:,:,k) =  y(:,:,k) * this%LU1(k,3)- y(:,:,k+1) * this%LU1(k,4)- y(:,:,this%n) * this%LU1(k,5)
        end do
    
    end subroutine

   
    subroutine SolveZLU1_BC(this,y,n1,n2)

        class (ci06), intent(in) :: this
        integer, intent(in) :: n1,n2
        real(rkind), dimension(n1,n2,this%n), intent(inout) :: y  ! Take inRHSand put solution into it
        integer ::  k
        real(rkind), dimension(n1,n2) :: sum1

        ! Step 2
        sum1 = this%LU1_BC(1,2)*y(:,:,2)
        do k = 3,this%n-3
            y(:,:,k) = y(:,:,k) - this%LU1_BC(k-1,1)*y(:,:,k-1)
            sum1 = sum1 + this%LU1_BC(k-1,2)*y(:,:,k)
        end do
        y(:,:,this%n-2) = y(:,:,this%n-2) - sum1

        ! Step 3
        y(:,:,this%n-2)   = y(:,:,this%n-2) * this%LU1_BC(this%n-3,3)

        y(:,:,this%n-3) =  y(:,:,this%n-3) * this%LU1_BC(this%n-4,3)-y(:,:,this%n-2) * this%LU1_BC(this%n-4,5)
        do k = this%n-4,2,-1
            y(:,:,k) =  y(:,:,k) * this%LU1_BC(k-1,3)- y(:,:,k+1) *this%LU1_BC(k-1,4)-y(:,:,this%n-2) * this%LU1_BC(k-1,5)
        end do

    end subroutine
 
    subroutine SolveXTri1(this,y,n2,n3)
    
        class (ci06), intent(in) :: this
        integer, intent(in) :: n2,n3
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y  
        integer :: i, j, k
  
        do k = 1,n3
            do j = 1,n2
                y(1,j,k) = y(1,j,k)*this%Tri1(1,2)
                do i = 2,this%n
                    y(i,j,k) = y(i,j,k)*this%Tri1(i,2) - y(i-1,j,k)*this%Tri1(i,1)
                end do
                
                do i = this%n-1,1,-1
                    y(i,j,k) = y(i,j,k) - this%Tri1(i,3)*y(i+1,j,k)
                end do 
            end do 
        end do 
        
    end subroutine
   
    subroutine SolveYTri1(this,y,n1,n3)
    
        class (ci06), intent(in) :: this
        integer, intent(in) :: n1,n3
        real(rkind), dimension(n1,this%n,n3), intent(inout) :: y  
        integer :: j, k
  
        do k = 1,n3 
                y(:,1,k) = y(:,1,k)*this%Tri1(1,2)
                do j = 2,this%n
                    y(:,j,k) = y(:,j,k)*this%Tri1(j,2) - y(:,j-1,k)*this%Tri1(j,1)
                end do
                
                do j = this%n-1,1,-1
                    y(:,j,k) = y(:,j,k) - this%Tri1(j,3)*y(:,j+1,k)
                end do 
        end do 
        
    end subroutine
    
    subroutine SolveZTri1(this,y,n1,n2)
    
        class (ci06), intent(in) :: this
        integer, intent(in) :: n1,n2
        real(rkind), dimension(n1,n2,this%n), intent(inout) :: y  
        integer :: k
  
                y(:,:,1) = y(:,:,1)*this%Tri1(1,2)
                do k = 2,this%n
                    y(:,:,k) = y(:,:,k)*this%Tri1(k,2) - y(:,:,k-1)*this%Tri1(k,1)
                end do
                
                do k = this%n-1,1,-1
                    y(:,:,k) = y(:,:,k) - this%Tri1(k,3)*y(:,:,k+1)
                end do 
        
    end subroutine

    
    subroutine ComputeXD1RHS(this,f,RHS, dir, n2, n3,bc1,bcn) 
         
        class( ci06 ), intent(in) :: this
        integer, intent(in) :: n2, n3
        integer,optional, intent(in) ::bc1,bcn
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer ::  j, k
        real(rkind) :: a06, b06,a10,a104,a102,b10,b102,b104,c10,c102,c104,a101
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1

        select case (this%periodic)
        case (.TRUE.)
            a06 = a06d1 
            b06 = b06d1 
            RHS = 0.0d0
         select case(dir)
           case("N2F")
   
            do k = 1,n3
                do j = 1,n2
                    RHS(1         ,j,k) = a06 * ( f(2,j,k)          + f(1       ,j,k) ) &
                                        + b06 * ( f(3,j,k)          + f(this%n  ,j,k) ) 
                    

                    RHS(2:this%n-2,j,k) = a06 * ( f(3:this%n-1,j,k) + f(2:this%n-2,j,k) ) &
                                        + b06 * ( f(4:this%n  ,j,k) + f(1:this%n-3,j,k) ) 
                    
                    RHS(this%n-1  ,j,k) = a06 * ( f(this%n,j,k)          + f(this%n-1,j,k) ) &
                                        + b06 * ( f(1,j,k)          + f(this%n-2,j,k) )

                    
                    RHS(this%n    ,j,k) = a06 * ( f(1,j,k)          + f(this%n  ,j,k) ) &
                                        + b06 * ( f(2,j,k)          + f(this%n-1,j,k) )
                end do 
            end do 
          case("F2N")
    
           !!!!!!!!!!!!!! TO DO         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
        end select
        case (.FALSE.)
        
            a10 = aI02 * half
            b10 = bI02 * half
            c10 = cI02 * half
            a101 =  one
            a102 =  half
            a104 =  9.0d0/16.0d0
            b104 = -1.0d0/16.0d0
            RHS = 0.0d0
            select case (dir)
                case ("N2F")!TODO: implement better non-periodic BC: currently2466...6642
                    !interior    
                    
                     RHS(2:this%n-2,j,k) = a06 * ( f(3:this%n-1,j,k) +f(2:this%n-2,j,k) ) &
                                        + b06 * ( f(4:this%n  ,j,k) +f(1:this%n-3,j,k) )

                     select case(bc1)
                        !left boundary (1:2)
                        case(1) !symm
                            RHS(1,:,:) = RHS(1,:,:) + a10 * (f(2,:,:) +f(1,:,:)) &
                                                    + b10 * (f(3,:,:) +f(2,:,:)) &
                                                    + c10 * (f(4,:,:) +f(3,:,:))
                        case(-1) !anti-symm
                            RHS(1,:,:) = RHS(1,:,:) + a10 * (f(2,:,:) +f(1,:,:)) &
                                                    + b10 * (f(3,:,:) -f(2,:,:)) &
                                                    + c10 * (f(4,:,:) -f(3,:,:))
                        case(0)
                            RHS(1,:,:) = RHS(1,:,:) + a102 * (f(2,:,:) +f(1,:,:))   !2nd order
                    end select

                    select case(bcn)
                        !right boundary (n-2:n-1)
                        case(1)
                            RHS(this%n-1,:,:) = RHS(this%n-1,:,:) + a10 *(f(this%n  ,:,:) + f(this%n-1,:,:)) &
                                                                  + b10 *(f(this%n-1,:,:) + f(this%n-2,:,:)) &
                                                                  + c10 *(f(this%n-2,:,:) + f(this%n-3,:,:))
                        case(-1)
                            RHS(this%n-1,:,:) = RHS(this%n-1,:,:) + a10 * (f(this%n  ,:,:) + f(this%n-1,:,:)) &
                                                                  + b10 *(-f(this%n-1,:,:) + f(this%n-2,:,:)) &
                                                                  + c10 *(-f(this%n-2,:,:) + f(this%n-3,:,:))
                        case(0)
                            RHS(this%n-1,:,:) = RHS(this%n-1,:,:) + a102 *(f(this%n,  :,:) + f(this%n-1,  :,:))   !2nd order
                    end select

                case ("F2N")  

                end select
        
        end select
    
    end subroutine
   


    subroutine ComputeYD1RHS(this,f, RHS,dir, n1, n3,bc1,bcn) 
         
        class( ci06 ), intent(in) :: this
        integer, intent(in) :: n1, n3
        integer,optional, intent(in) ::bc1,bcn
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer ::  k
        real(rkind) :: a06, b06, a10,a104,a102,b10,b102,b104,c10,c102,c104,a101
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1

        select case (this%periodic)
        case (.TRUE.)
            a06 = a06d1 
            b06 = b06d1 
            RHS = 0.0d0
           
          select case(dir)
            case("N2F")
            do k = 1,n3

                    RHS(:, 1    ,k) = a06 * ( f(:,2,  k) +  f(:,1,k) ) &
                                        + b06 * ( f(:,3,k)          + f(:,this%n,k) )


                    RHS(:,2:this%n-2,k) = a06 * ( f(:,3:this%n-1,k) + f(:,2:this%n-2,k) ) &
                                        + b06 * ( f(:,4:this%n ,k) + f(:,1:this%n-3,k) ) 

                    RHS(:,this%n-1 ,k) = a06 * ( f(:,this%n,k)          +f(:,this%n-1,k) ) &
                                        + b06 * ( f(:,1,k)          + f(:,this%n-2,k) )


                    RHS(:,this%n    ,k) = a06 * ( f(:,1,k)          + f(:,this%n,k) ) &
                                        + b06 * ( f(:,2,k)          + f(:,this%n-1,k) )

                !end do 
            end do 
            case("F2N")
              
            !!!!!!!!!!!!!!!!!!!!!!!!!!   TO DO           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         end select
        case (.FALSE.)
            a10 = aI02 * half
            b10 = bI02 * half
            c10 = cI02 * half
            a101 =  one
            a102 =  half
            a104 =  9.0d0/16.0d0
            b104 = -1.0d0/16.0d0
            RHS = 0.0d0
            select case (dir)
                case ("N2F")

                        a10 = aI02 * half
            b10 = bI02 * half
            c10 = cI02 * half
            a101 =  one
            a102 =  half
            a104 =  9.0d0/16.0d0
            b104 = -1.0d0/16.0d0
            RHS = 0.0d0
            select case (dir)
                case ("N2F")!TODO: implement better non-periodic BC:currently2466...6642
                    !interior   
                    RHS(:,2:this%n-2,k) = a06 * ( f(:,3:this%n-1,k) +f(:,2:this%n-2,k) ) &
                                        + b06 * ( f(:,4:this%n ,k) + f(:,1:this%n-3,k) )

                    select case(bc1)
                        !left boundary (1:2)
                        case(1) !symm
                            RHS(:,1,:) = RHS(:,1,:) + a10 * (f(:,2,:) +f(:,1,:))&
                                                    + b10 * (f(:,3,:) +f(:,2,:))&
                                                    + c10 * (f(:,4,:) +f(:,3,:))
                        case(-1) !anti-symm
                            RHS(:,1,:) = RHS(:,1,:) + a10 * (f(:,2,:) +f(:,1,:))&
                                                    + b10 * (f(:,3,:) -f(:,2,:))&
                                                    + c10 * (f(:,4,:) -f(:,3,:))
                        case(0)
                            RHS(:,1,:) = RHS(:,1,:) + a102 * (f(:,2,:)+f(:,1,:))   !2nd order
                    end select

                    select case(bcn)
                        !right boundary (n-2:n-1)
                        case(1)
                            RHS(:,this%n-1,:) = RHS(:,this%n-1,:) + a10*(f(:,this%n  ,:) + f(:,this%n-1,:)) &
                                                                  + b10*(f(:,this%n-1,:) + f(:,this%n-2,:)) &
                                                                  + c10*(f(:,this%n-2,:) + f(:,this%n-3,:))
                        case(-1)
                            RHS(:,this%n-1,:) = RHS(:,this%n-1,:) + a10 *(f(:,this%n  ,:) + f(:,this%n-1,:)) &
                                                                  + b10*(-f(:,this%n-1,:) + f(:,this%n-2,:)) &
                                                                  + c10*(-f(:,this%n-2,:) + f(:,this%n-3,:))
                        case(0)
                            RHS(:,this%n-1,:) = RHS(:,this%n-1,:) + a102*(f(:,this%n,  :) + f(:,this%n-1,  :))   !2nd order
                    end select

                case ("F2N")

                end select
            end select 
        end select
    
    end subroutine


    subroutine ComputeZD1RHS(this,f, RHS,dir, n1, n2,bc1,bcn) 
         
        class( ci06 ), intent(in) :: this
        integer, intent(in) :: n1, n2
        integer,optional, intent(in) ::bc1,bcn
        real(rkind), dimension(n1,n2,this%n), intent(in) :: f
        real(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        real(rkind) :: a06, b06, a10,a104,a102,b10,b102,b104,c10,c102,c104,a101
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1

        select case (this%periodic)
        case (.TRUE.)
            a06 = a06d1 
            b06 = b06d1 
            RHS = 0.0d0
        select case(dir)
           case("N2F")

            
            RHS(:,:,        1)   = a06 * ( f(:,:,2)     + f(:,:,    1) ) &
                                        + b06 * ( f(:,:,3)    + f(:,:,this%n) )


            RHS(:,:,2:this%n-2)  = a06 * ( f(:,:,3:this%n-1) + f(:,:,2:this%n-2) ) &
                                        + b06 * ( f(:,:,4:this%n) + f(:,:,1:this%n-3) ) 

            RHS(:,:,this%n-1)    = a06 * ( f(:,:,this%n)     + f(:,:,this%n-1) ) &
                                        + b06 * ( f(:,:,1)  + f(:,:,this%n-2) )


            RHS(:,:,this%n)      = a06 * ( f(:,:,1)          + f(:,:,this%n) ) &
                                        + b06 * ( f(:,:,2)  + f(:,:,this%n-1) )


           case("F2N")
           !!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!!!
           end select
        case (.FALSE.)
          
            a10 = aI02 * half
            b10 = bI02 * half
            c10 = cI02 * half
            a101 =  one
            a102 =  half
            a104 =  9.0d0/16.0d0
            b104 = -1.0d0/16.0d0
            RHS = 0.0d0
            select case (dir)
                case ("N2F")!TODO: implement better non-periodic BC: currently2466...6642
                    !interior   
                     RHS(:,:,2:this%n-2)  = a06 * ( f(:,:,3:this%n-1) +f(:,:,2:this%n-2) ) &
                                        + b06 * ( f(:,:,4:this%n) +f(:,:,1:this%n-3) )

                                         select case(bc1)
                        !left boundary (1:2)
                        case(1) !symm
                            RHS(:,:,1) = RHS(:,:,1) + a10 * (f(:,:,2) +f(:,:,1)) &
                                                    + b10 * (f(:,:,3) +f(:,:,2)) &
                                                    + c10 * (f(:,:,4) +f(:,:,3))
                        case(-1) !anti-symm
                            RHS(:,:,1) = RHS(:,:,1) + a10 * (f(:,:,2) +f(:,:,1)) &
                                                    + b10 * (f(:,:,3) -f(:,:,2)) &
                                                    + c10 * (f(:,:,4) -f(:,:,3))
                        case(0)
                            RHS(:,:,1) = RHS(:,:,1) + a102 * (f(:,:,2) +f(:,:,1))   !2nd order
                    end select

                    select case(bcn)
                        !right boundary (n-2:n-1)
                        case(1)
                            RHS(:,:,this%n-1) = RHS(:,:,this%n-1) + a10 *(f(:,:,this%n  ) + f(:,:,this%n-1)) &
                                                                  + b10 *(f(:,:,this%n-1) + f(:,:,this%n-2)) &
                                                                  + c10 *(f(:,:,this%n-2) + f(:,:,this%n-3))
                        case(-1)
                            RHS(:,:,this%n-1) = RHS(:,:,this%n-1) + a10 * (f(:,:,this%n  ) + f(:,:,this%n-1)) &
                                                                  + b10 *(-f(:,:,this%n-1) + f(:,:,this%n-2)) &
                                                                  + c10 *(-f(:,:,this%n-2) + f(:,:,this%n-3))
                        case(0)
                            RHS(:,:,this%n-1) = RHS(:,:,this%n-1) + a102 *(f(:,:,this%n) + f(:,:,this%n-1))   !2nd order
                    end select



                case ("F2N")

                end select
        
        end select
    
    end subroutine
   

    subroutine iN2F1(this, fN, fF, na, nb,bc1_,bcn_)
        class( ci06 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(this%n,na,nb), intent(in)  :: fN
        real(rkind), dimension(this%n,na,nb), intent(out) :: fF

        if (this%n == 1) then
            fF = zero
            return
        end if

        if (present(bc1_)) then
            bc1 = bc1_
            if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bc1(should be 0, 1 or -1)", 324)
            end if
        else
            bc1 = 0
        end if

        if (present(bcn_)) then
            bcn = bcn_
            if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bcn(should be 0, 1 or -1)", 324)
            end if
        else
            bcn = 0
        end if


        call this%ComputeXD1RHS(fN, fF,"N2F", na, nb,bc1,bcn)
        select case (this%periodic) 
        case (.TRUE.)
           call this%SolveXLU1(fF,na,nb)
        case (.FALSE.) 
           call this%SolveXLU1_BC(fF,na,nb)
        end select


    end subroutine

    subroutine iF2N1(this, fF, fN, na, nb)
        class( ci06 ), intent(in) :: this
        integer, intent(in) :: na, nb
        real(rkind), dimension(this%n,na,nb), intent(in)  :: fF
        real(rkind), dimension(this%n,na,nb), intent(out) :: fN

        if (this%n == 1) then
            fN = zero
            return
        end if

        call this%ComputeXD1RHS(fF, fN,"F2N", na, nb)
        select case (this%periodic)
        case (.TRUE.)
           call this%SolveXLU1(fN,na,nb)
        case (.FALSE.)
           call this%SolveXTri1(fN,na,nb)
        end select


    end subroutine


    subroutine iN2F2(this, fN, fF, na, nb,bc1_,bcn_)
        class( ci06 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(na,this%n,nb), intent(in)  :: fN
        real(rkind), dimension(na,this%n,nb), intent(out) :: fF

        if (this%n == 1) then
            fF = zero
            return
        end if

        if (present(bc1_)) then
            bc1 = bc1_
            if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bc1(should be 0, 1 or -1)", 324)
            end if
        else
            bc1 = 0
        end if

        if (present(bcn_)) then
            bcn = bcn_
            if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bcn(should be 0, 1 or -1)", 324)
            end if
        else
            bcn = 0
        end if

        call this%ComputeYD1RHS(fN, fF,"N2F", na, nb,bc1,bcn)
        select case (this%periodic) 
        case (.TRUE.)
           call this%SolveYLU1(fF,na,nb)
        case (.FALSE.) 
           call this%SolveYLU1_BC(fF,na,nb)
        end select


    end subroutine
    

    subroutine iF2N2(this, fF, fN, na, nb)
        class( ci06 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer :: bc1,bcn
        real(rkind), dimension(na,this%n,nb), intent(in)  :: fF
        real(rkind), dimension(na,this%n,nb), intent(out) :: fN

        if (this%n == 1) then
            fN = zero
            return
        end if

        call this%ComputeYD1RHS(fF, fN,"F2N", na, nb)
        select case (this%periodic)
        case (.TRUE.)
           call this%SolveYLU1(fN,na,nb)
        case (.FALSE.)
           call this%SolveYTri1(fN,na,nb)
        end select


    end subroutine

    
    subroutine iN2F3(this, fN, fF, na, nb,bc1_,bcn_)
        class( ci06 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(na,nb,this%n), intent(in)  :: fN
        real(rkind), dimension(na,nb,this%n), intent(out) :: fF

        if (this%n == 1) then
            fF = zero
            return
        end if
        
        if (present(bc1_)) then
            bc1 = bc1_
            if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bc1(should be 0, 1 or -1)", 324)
            end if
        else
            bc1 = 0
        end if

        if (present(bcn_)) then
            bcn = bcn_
            if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bcn(should be 0, 1 or -1)", 324)
            end if
        else
            bcn = 0
        end if

        call this%ComputeZD1RHS(fN, fF,"N2F", na, nb,bc1,bcn)
        select case (this%periodic) 
        case (.TRUE.)
           call this%SolveZLU1(fF,na,nb)
        case (.FALSE.) 
           call this%SolveZLU1_BC(fF,na,nb)
        end select


    end subroutine
    
    subroutine iF2N3(this, fF, fN, na, nb)
        class( ci06 ), intent(in) :: this
        integer, intent(in) :: na, nb
        real(rkind), dimension(na,nb,this%n), intent(in)  :: fF
        real(rkind), dimension(na,nb,this%n), intent(out) :: fN

        if (this%n == 1) then
            fN = zero
            return
        end if

        call this%ComputeZD1RHS(fF, fN,"F2N", na, nb)
        select case (this%periodic)
        case (.TRUE.)
           call this%SolveZLU1(fN,na,nb)
        case (.FALSE.)
           call this%SolveZTri1(fN,na,nb)
        end select


    end subroutine


end module
