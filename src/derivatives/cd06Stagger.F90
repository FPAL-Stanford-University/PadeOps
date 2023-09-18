! Routines specific to 6th order Compact Finite Differencing scheme
! Periodic LU based on  Neossi Nguetchue, Abelman (Appl. Math. & Comp. 2008)

module cd06Staggerstuff

    use kind_parameters, only: rkind
    use constants,       only : zero,one,two,third,fifth
    use exits,    only: gracefulExit, message
    implicit none

    private
    public :: cd06Stagger, alpha06d1, a06d1, b06d1
    
    ! 6th order first derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha06d1=  9.0_rkind / 62.0_rkind
    real(rkind), parameter :: a06d1    = (63.0_rkind / 62.0_rkind) 
    real(rkind), parameter :: b06d1    = ( 17.0_rkind / 62.0_rkind) / 3.0_rkind

    ! 6th order first derivative explicit centeral difference coefficients
    real(rkind), parameter :: aD06d1     = 3.0_rkind / 4.0_rkind
    real(rkind), parameter :: bD06d1     = -3.0_rkind / 20.0_rkind
    real(rkind), parameter :: cD06d1     = 1.0_rkind / 60.0_rkind
    
    ! 6th order second derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha06d2=  2.0_rkind / 11.0_rkind
    real(rkind), parameter :: a06d2    = (12.0_rkind / 11.0_rkind) / 1.0_rkind
    real(rkind), parameter :: b06d2    = ( 3.0_rkind / 11.0_rkind) / 4.0_rkind
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic 1st derivative evaluation !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        ! Set the scheme for the edge nodes (Ref. for notation: Lele - JCP
        ! paper)
    ! 1st derivative
    real(rkind), parameter                   :: alpha =  0._rkind !3._rkind
    real(rkind), parameter                   :: p     =  -71._rkind / 24._rkind !17._rkind / 6._rkind
    real(rkind), parameter                   :: q     =  47._rkind / 8._rkind   !3._rkind / 2._rkind
    real(rkind), parameter                   :: r     =  -31._rkind / 8._rkind  !3._rkind / 2._rkind
    real(rkind), parameter                   :: s     =  23._rkind  / 24._rkind !-1._rkind / 6._rkind

    real(rkind), parameter                   :: alpha_N2F =  0._rkind !3._rkind
    real(rkind), parameter                   :: p_N2F     =  -23._rkind / 24._rkind !17._rkind / 6._rkind
    real(rkind), parameter                   :: q_N2F     =  7._rkind / 8._rkind  !3._rkind / 2._rkind
    real(rkind), parameter                   :: r_N2F     =  1._rkind / 8._rkind !3._rkind / 2._rkind
    real(rkind), parameter                   :: s_N2F     =  -1._rkind  / 24._rkind !-1._rkind / 6._rkind

    !!  Calculate the corressponding weights
    ! Step 1: Assign the interior scheme
    real(rkind), parameter                   :: qhat        = a06d1
    real(rkind), parameter                   :: rhat        = b06d1
    real(rkind), parameter                   :: alpha_hat   = alpha06d1

    ! Step 2: Assign the scheme at node 2 to be Standard Pade (4th Order)
    real(rkind), parameter                   :: p_p         = 3._rkind/8._rkind*(3._rkind - 2._rkind / 22._rkind)
    real(rkind), parameter                   :: q_p         = 1._rkind/8._rkind*( -1._rkind + 1._rkind)
    real(rkind), parameter                   :: alpha_p     = 1._rkind /22._rkind     !1._rkind/4._rkind

    ! Step 3: Get the scheme at the node 3
    real(rkind), parameter                   :: alpha_pp = 9089._rkind /69564._rkind    !((40*alpha_hat - 1)*q  + 7*(4*alpha_hat &                                                           
    real(rkind), parameter                   :: q_pp     = 95257._rkind /92752._rkind    !(1._rkind/3._rkind)*(alpha_pp + 2)
    real(rkind), parameter                   :: r_pp     = 5927._rkind  /25296._rkind / 3._rkind   !(1._rkind/12._rkind)*(4*alpha_pp - 1)
    real(rkind), parameter                   :: s_pp     =  0._rkind

    ! Step 4: Get the weights
    real(rkind), parameter                   :: w1 = 181._rkind / 31620._rkind   !(2*alpha_hat + 1)/(2*(q + s))
    real(rkind), parameter                   :: w2 = 11473._rkind / 10540._rkind !((8*alpha_hat + 7)*q - 6*(2*alpha_hat + 1)*r &
    real(rkind), parameter                   :: w3 = 1.1
    
    real(rkind), parameter                   :: w1_N2F = 223._rkind /186._rkind
    real(rkind), parameter                   :: w2_N2F = 61._rkind/62._rkind
    real(rkind), parameter                   :: w3_N2F = 1._rkind
    real(rkind), parameter                   :: alpha_p_N2F = 37._rkind / 183._rkind
    real(rkind), parameter                   :: p_p_N2F = 3._rkind/8._rkind*(3._rkind - 2._rkind *alpha_p_N2F)
    real(rkind), parameter                   :: q_p_N2F = 1._rkind/8._rkind*(-1._rkind + alpha_p_N2F * 22._rkind) / 3._rkind
    real(rkind), parameter                   :: alpha_pp_N2F = alpha06d1    !((40*alpha_hat - 1)*q  + 7*(4*alpha_hat &                                            
    real(rkind), parameter                   :: q_pp_N2F     = a06d1        !(1._rkind/3._rkind)*(alpha_pp + 2)
    real(rkind), parameter                   :: r_pp_N2F     = b06d1        !(1._rkind/12._rkind)*(4*alpha_pp - 1)
    real(rkind), parameter                   :: s_pp_N2F     =  0._rkind

 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic 2nd derivative evaluation !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    ! Incomplete
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    type cd06Stagger

        private

        integer     :: n
        real(rkind) :: dx
        real(rkind) :: onebydx
        real(rkind) :: onebydx2

        logical     :: periodic = .TRUE.
        integer     :: bc1 = 0                             ! Boundary condition type. 0=Dirichlet, 1=Neumann
        integer     :: bcn = 0                             ! Boundary condition type. 0=Dirichlet, 1=Neumann

        real(rkind), allocatable, dimension(:,:) :: LU1, LU1_BC
        real(rkind), allocatable, dimension(:,:) :: LU2
        real(rkind), allocatable, dimension(:,:) :: Tri1, Tri1_N2F
        real(rkind), allocatable, dimension(:,:) :: Tri2

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
       
        procedure, private :: ComputeTri1
        procedure, private :: ComputeTri2
        procedure, private :: ComputeTri1_N2F

        procedure, private :: SolveXTri1
        procedure, private :: SolveYTri1
        procedure, private :: SolveZTri1

        procedure, private :: SolveXTri1_N2F
        procedure, private :: SolveYTri1_N2F
        procedure, private :: SolveZTri1_N2F
        
        procedure, private :: SolveTri2
        
        procedure :: dd1F2N
        procedure :: dd2F2N
        procedure :: dd3F2N
        procedure :: dd1N2F
        procedure :: dd2N2F
        procedure :: dd3N2F
        procedure :: cd06Staggerder2

    end type

contains

    pure function GetSize(this) result(val)
        class(cd06Stagger), intent(in) :: this
        integer  :: val 
        val = this%n
    end function
    
    function init(this, n_, dx_, periodic_, bc1_, bcn_) result(ierr)
    
        class( cd06Stagger ), intent(inout) :: this
        integer, intent(in) :: n_
        real(rkind), intent(in) :: dx_
        logical, intent(in) :: periodic_
        integer, intent(in) :: bc1_, bcn_
        integer :: ierr, n1
    
        this%n = n_
        this%dx = dx_
        this%onebydx = one/dx_
        this%onebydx2 = this%onebydx/dx_
        n1 = this%n - 1
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
    
            ! Allocate 2nd derivative LU matrix.
            if(allocated( this%LU2 )) deallocate( this%LU2 ); allocate( this%LU2(n_,5) )
    
            ! Compute 2nd derivative LU matrix
            if (n_ .GE. 6) then
                call ComputeLU(this%LU2,n_,alpha06d2,one,alpha06d2)
            else if (n_ == 1) then
                this%LU2 = one
            end if
        
        else
            ! Allocate 1st derivative Tri matrices.
            if(allocated( this%tri1 )) deallocate( this%tri1 ); allocate( this%tri1(n_,3) )
            if(allocated( this%tri1_N2F )) deallocate( this%tri1_N2F ); allocate(this%tri1_N2F(n1,3) )
 
            if (n_ .GE. 8) then             

                call this%ComputeTri1(bc1_,bcn_)
                call this%ComputeTri1_N2F(bc1_,bcn_)
                
            else if (n_ .EQ. 1) then
                this%Tri1 = one
                this%Tri1_N2F = one
            else
                ierr = 2
                return 
            end if 

            
            ! Allocate 2nd derivative Tri matrices.
            if(allocated( this%tri2 )) deallocate( this%tri2 ); allocate( this%tri2(n_,3) )
        
            if (n_ .GE. 8) then             
                call this%ComputeTri2(bc1_,bcn_)
            else if (n_ .EQ. 1) then
                this%Tri2 = one
            end if 
              
        end if 

    
        ! If everything passes
        ierr = 0
    
    end function

    subroutine destroy(this)

        class( cd06Stagger ), intent(inout) :: this

        ! Dellocate 1st derivative LU matrix.
        if(allocated( this%LU1 )) deallocate( this%LU1 )
    
        ! Dellocate 2nd derivative LU matrix.
        if(allocated( this%LU2 )) deallocate( this%LU2 )

        if(allocated( this%LU1_BC )) deallocate( this%LU1_BC )

        ! Dellocate 1st derivative tri matrix.
        if(allocated( this%tri1 )) deallocate( this%tri1 )

        ! Dellocate 2nd derivative tri matrix.
        if(allocated( this%tri1_N2F )) deallocate( this%tri1_N2F )
  
        ! Dellocate 2nd derivative tri matrix.
        if(allocated( this%tri2 )) deallocate( this%tri2 )

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
        class (cd06Stagger), intent(inout) :: this
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
    
   subroutine ComputeTri1_N2F(this,bc1,bcn)
        class (cd06Stagger), intent(inout) :: this
        integer, intent(in) :: bc1, bcn
        integer             :: i
        real(rkind), dimension(this%n-1) :: a, b, c, cp, den

        a  = alpha_hat
        b  = one
        c  = alpha_hat
        
        
        select case (bc1)
        case(0)
            a (1) = w1_N2F*zero
            a (2) = w2_N2F*alpha_p_N2F
            a (3) = w3_N2F*alpha_pp_N2F

            b (1) = w1_N2F*one
            b (2) = w2_N2F*one
            b (3) = w3_N2F*one

            c (1) = w1_N2F*alpha_N2F
            c (2) = w2_N2F*alpha_p_N2F
            c (3) = w3_N2F*alpha_pp_N2F

        case(1)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Incomplete
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end select

        select case (bcn)
        case(0)
            c (this%n-1  ) = w1_N2F*zero
            c (this%n-2) = w2_N2F*alpha_p_N2F
            c (this%n-3) = w3_N2F*alpha_pp_N2F

            b (this%n-1) = w1_N2F*one
            b (this%n-2) = w2_N2F*one
            b (this%n-3) = w3_N2F*one

            a (this%n-1) = w1_N2F*alpha_N2F
            a (this%n-2) = w2_N2F*alpha_p_N2F
            a (this%n-3) = w3_N2F*alpha_pp_N2F
        case(1)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Incomplete
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end select

        cp(1) = c(1)/b(1)
        do i = 2,this%n-2
            cp(i) = c(i)/(b(i) - a(i)*cp(i-1))
        end do

        den(1) = one/b(1)
        den(2:this%n-1) = one/(b(2:this%n-1) - a(2:this%n-1)*cp(1:this%n-2))

        this%Tri1_N2F(:,1) = a*den
        this%Tri1_N2F(:,2) = den
        this%Tri1_N2F(:,3) = cp
    end subroutine
 
    subroutine ComputeTri2(this,bc1,bcn)
        class (cd06Stagger), intent(inout) :: this
        integer, intent(in) :: bc1, bcn
   
        

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        if ((bc1 == 1) .and. (bcn == 1)) then
            this%tri2 = zero
        end if 
    end subroutine 
    
    subroutine SolveXLU1(this,y,n2,n3)
        
        class (cd06Stagger), intent(in) :: this
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

    subroutine SolveYLU1(this,y,n1,n3)
        
        class (cd06Stagger), intent(in) :: this
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
        
        class (cd06Stagger), intent(in) :: this
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

    subroutine SolveXTri1(this,y,n2,n3)

        class (cd06Stagger), intent(in) :: this
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

    subroutine SolveXTri1_N2F(this,y,n2,n3)
    
        class (cd06Stagger), intent(in) :: this
        integer, intent(in) :: n2,n3
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y  
        integer :: i, j, k
  
        do k = 1,n3
            do j = 1,n2
                y(1,j,k) = y(1,j,k)*this%Tri1_N2F(1,2)
                do i = 2,this%n-1
                    y(i,j,k) = y(i,j,k)*this%Tri1_N2F(i,2) - y(i-1,j,k)*this%Tri1_N2F(i,1)
                end do
                
                do i = this%n-2,1,-1
                    y(i,j,k) = y(i,j,k) - this%Tri1_N2F(i,3)*y(i+1,j,k)
                end do 
            end do 

        end do 

        y(this%n,:,:) = 0        
    end subroutine
   
    subroutine SolveYTri1(this,y,n1,n3)
    
        class (cd06Stagger), intent(in) :: this
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
  
    subroutine SolveYTri1_N2F(this,y,n1,n3)

        class (cd06Stagger), intent(in) :: this
        integer, intent(in) :: n1,n3
        real(rkind), dimension(n1,this%n,n3), intent(inout) :: y
        integer :: j, k

        do k = 1,n3
                y(:,1,k) = y(:,1,k)*this%Tri1_N2F(1,2)
                do j = 2,this%n-1
                    y(:,j,k) = y(:,j,k)*this%Tri1_N2F(j,2) - y(:,j-1,k)*this%Tri1_N2F(j,1)
                end do

                do j = this%n-2,1,-1
                    y(:,j,k) = y(:,j,k) - this%Tri1_N2F(j,3)*y(:,j+1,k)
                end do
        end do
         y(:,this%n,:) = 0
    end subroutine
  
    subroutine SolveZTri1(this,y,n1,n2)
    
        class (cd06Stagger), intent(in) :: this
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

    subroutine SolveZTri1_N2F(this,y,n1,n2)

        class (cd06Stagger), intent(in) :: this
        integer, intent(in) :: n1,n2
        real(rkind), dimension(n1,n2,this%n), intent(inout) :: y
        integer :: k

                y(:,:,1) = y(:,:,1)*this%Tri1_N2F(1,2)
                do k = 2,this%n-1
                    y(:,:,k) = y(:,:,k)*this%Tri1_N2F(k,2) - y(:,:,k-1)*this%Tri1_N2F(k,1)
                end do

                do k = this%n-2,1,-1
                    y(:,:,k) = y(:,:,k) - this%Tri1_N2F(k,3)*y(:,:,k+1)
                end do
                y(:,:,this%n) = 0
    end subroutine


    subroutine SolveLU2(this,y,n2,n3)
        
        class (cd06Stagger), intent(in) :: this
        integer, intent(in) :: n2,n3
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y  ! Take in RHS and put solution into it

        
        y = zero
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine 
     
    subroutine SolveTri2(this,y,n2,n3)
    
        class (cd06Stagger), intent(in) :: this
        integer, intent(in) :: n2,n3
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y  
  
        y = zero 
    end subroutine
    
    subroutine ComputeXD1RHS(this,f, RHS, dir,n2, n3,bc1,bcn) 
         
        class( cd06Stagger ), intent(in) :: this
        integer, intent(in) :: n2, n3,bc1,bcn
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer ::  j, k
        real(rkind) :: a06, b06,a10,a102,a104,b10,b104,b102,c10,c102,c104,a101
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2, b_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1

        select case (this%periodic)
        case (.TRUE.)
            a06 = a06d1 * this%onebydx
            b06 = b06d1 * this%onebydx
            RHS = 0.0d0
        select case(dir)
           case("F2N")  
            do k = 1,n3
                do j = 1,n2
                    RHS(1         ,j,k) = a06 * ( f(1,j,k)          - f(this%n  ,j,k) ) &
                                        + b06 * ( f(2,j,k)          - f(this%n-1,j,k) ) 
                    
                    RHS(2         ,j,k) = a06 * ( f(2,j,k)          - f(1       ,j,k) ) &
                                        + b06 * ( f(3,j,k)          - f(this%n  ,j,k) )

                    RHS(3:this%n-2,j,k) = a06 * ( f(3:this%n-2,j,k) - f(2:this%n-3,j,k) ) &
                                        + b06 * ( f(4:this%n-1,j,k) - f(1:this%n-4,j,k) ) 
                    
                    RHS(this%n-1  ,j,k) = a06 * ( f(this%n-1,j,k)         - f(this%n-2,j,k) ) &
                                        + b06 * ( f(this%n,j,k)          - f(this%n-3,j,k) ) 
                    
                    RHS(this%n    ,j,k) = a06 * ( f(this%n,j,k)          - f(this%n-1,j,k) ) &
                                        + b06 * ( f(1,j,k)          - f(this%n-2,j,k) )
                end do 
            end do 
          case("N2F")
          !!!!!!!!!!!!!!!!!!!!! To do     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do k = 1,n3
                do j = 1,n2
                    RHS(1         ,j,k) = a06 * ( f(2,j,k)          - f(1,j,k) ) &
                                        + b06 * ( f(3,j,k)          - f(this%n,j,k) )


                    RHS(2:this%n-2,j,k) = a06 * ( f(3:this%n-1,j,k) - f(2:this%n-2,j,k) ) &
                                        + b06 * ( f(4:this%n,j,k)   - f(1:this%n-3,j,k) )

                    RHS(this%n-1  ,j,k) = a06 * ( f(this%n,j,k)     - f(this%n-1,j,k) ) &
                                        + b06 * ( f(1,j,k)          - f(this%n-2,j,k) )

                    RHS(this%n    ,j,k) = a06 * ( f(1,j,k)          -  f(this%n,j,k) ) &
                                        + b06 * ( f(2,j,k)          -  f(this%n-1,j,k) )
                end do
            end do

        end select

        case (.FALSE.)
       
            RHS = 0.0d0
            a06 = a06d1 * this%onebydx
            b06 = b06d1 * this%onebydx

            a_np_3 = w3*q_pp * this%onebydx
            b_np_3 = w3*r_pp * this%onebydx

            a_np_2 = w2*p_p * this%onebydx

            a_np_1 = w1*(  p * this%onebydx)
            b_np_1 = w1*(  q * this%onebydx)
            c_np_1 = w1*(  r * this%onebydx)
            d_np_1 = w1*(  s * this%onebydx)


            select case (dir)

               case("F2N")
               
                 do k = 1,n3
                   do j = 1,n2
              
                      RHS(4:this%n-3,j,k) = a06 * ( f(4:this%n-3,j,k) -f(3:this%n-4,j,k) ) &
                                        + b06 * ( f(5:this%n-2,j,k) -f(2:this%n-5,j,k) )
                    end do
                 end do

               select case(bc1)
                        ! left boundary (1:3)
                        case(1)

                        !!!!!!!!!!!!!!!!!!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!!
                        case(-1)

                        !!!!!!!!!!!!!!!!!!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!!!!!!!!!
                        case(0)
         
                           do k = 1, n3
                             do j = 1, n2
                                RHS(1,j,k) =   a_np_1* f(1,j,k) + b_np_1*f(2,j,k)   &
                                        +   c_np_1* f(3,j,k) + d_np_1*f(4,j,k)

                                RHS(2,j,k) =   a_np_2*(f(2,j,k) - f(1,j,k))

                                RHS(3,j,k) =   a_np_3*(f(3,j,k) - f(2,j,k)) &
                                        +   b_np_3*(f(4,j,k) - f(1,j,k))
                              end do
                           end do
                end select

                 select case(bcn)
                        !right boundary (n-2:n)
                        case(1)

                        !!!!!!!!!!!!!!!!!!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!!!
 
                        case(-1)

                        !!!!!!!!!!!!!!!!!!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!!!!!!!!!

                        case(0)
                         
                          do k = 1,n3
                             do j = 1, n2
                                 RHS(this%n-2  ,j,k) =   a_np_3*(f(this%n-2  ,j,k) - f(this%n-3  ,j,k)) &
                                        +   b_np_3*(f(this%n-1    ,j,k) - f(this%n-4  ,j,k))

                                 RHS(this%n-1  ,j,k) =   a_np_2*(f(this%n-1    ,j,k) - f(this%n-2  ,j,k))

                                 RHS(this%n    ,j,k) = -a_np_1* f(this%n-1    ,j,k) - b_np_1*f(this%n-2  ,j,k)   &
                                        -   c_np_1* f(this%n-3  ,j,k) - d_np_1*f(this%n-4  ,j,k)
                             end do
                          end do
                    end select

                case("N2F")
                   RHS = 0.0d0
                   a06 = a06d1 * this%onebydx
                   b06 = b06d1 * this%onebydx
                                       
                   a_np_2 = w2_N2F*p_p_N2F * this%onebydx
                   b_np_2 = w2_N2F*q_p_N2F * this%onebydx

                   a_np_1 = w1_N2F*(  p_N2F * this%onebydx)
                   b_np_1 = w1_N2F*(  q_N2F * this%onebydx)
                   c_np_1 = w1_N2F*(  r_N2F * this%onebydx)
                   d_np_1 = w1_N2F*(  s_N2F * this%onebydx)

                   do k = 1,n3
                      do j = 1,n2


                        RHS(3:this%n-3,j,k) = a06 * ( f(4:this%n-2,j,k) - f(3:this%n-3,j,k) ) &
                                        + b06 * ( f(5:this%n-1,j,k)   - f(2:this%n-4,j,k) )
                     end do
                   end do
  
                   select case(bc1)
                       case(1)
                       case(-1)
                       case(0)

                          do k = 1,n3
                              do j = 1,n2


                                RHS(1,j,k) = a_np_1*f(1,j,k) + b_np_1*f(2,j,k) + c_np_1*f(3,j,k) + d_np_1*f(4,j,k)
                                RHS(2,j,k) = a_np_2*( f(3,j,k) - f(2,j,k) ) + b_np_2*( f(4,j,k) - f(1,j,k))
                              end do
                             end do
                  end select

                  
                  select case(bcn)
                       case(1)
                       case(-1)
                       case(0)

                          do k = 1,n3
                              do j = 1,n2


                                RHS(this%n-1,j,k) = -a_np_1*f(this%n,j,k) - b_np_1*f(this%n-1,j,k) - c_np_1*f(this%n-2,j,k) - d_np_1*f(this%n-3,j,k)
                                RHS(this%n-2,j,k) = a_np_2*( f(this%n-1,j,k) - f(this%n-2,j,k) ) + b_np_2*( f(this%n,j,k) - f(this%n-3,j,k) )
                              end do
                          end do
                  end select


            end select           
        end select
    
    end subroutine
   


    subroutine ComputeYD1RHS(this,f, RHS,dir, n1, n3,bc1,bcn) 
         
        class( cd06Stagger ), intent(in) :: this
        integer, intent(in) :: n1, n3,bc1,bcn
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer ::  k
        real(rkind) :: a06, b06, a10,a102,a104,b10,b104,b102,c10,c102,c104,a101
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2, b_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1

        select case (this%periodic)
        case (.TRUE.)
            a06 = a06d1 * this%onebydx
            b06 = b06d1 * this%onebydx
            RHS = 0.0D0
         select case (dir)
           case("F2N")        
   
            do k = 1,n3
                    RHS(:,1         ,k) = a06 * ( f(:,1,k)          - f(:,this%n  ,k) ) &
                                        + b06 * ( f(:,2,k)          - f(:,this%n-1,k) ) 
                    
                    RHS(:,2         ,k) = a06 * ( f(:,2,k)          - f(:,1       ,k) ) &
                                        + b06 * ( f(:,3,k)          - f(:,this%n  ,k) )

                    RHS(:,3:this%n-2,k) = a06 * ( f(:,3:this%n-2,k) - f(:,2:this%n-3,k) ) &
                                        + b06 * ( f(:,4:this%n-1  ,k) - f(:,1:this%n-4,k) ) 
                    
                    RHS(:,this%n-1  ,k) = a06 * ( f(:,this%n-1,k)         - f(:,this%n-2,k) ) &
                                        + b06 * ( f(:,this%n,k)          - f(:,this%n-3,k) ) 
                    
                    RHS(:,this%n    ,k) = a06 * ( f(:,this%n,k)          - f(:,this%n-1,k) ) &
                                        + b06 * ( f(:,1,k)          - f(:,this%n-2,k) )
                !end do 
            end do 
           case("N2F")
         !!!!!!!!!!!!!!!!!!!!!!!!!!   TO DO        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do k = 1,n3

               RHS(:,1, k         )  = a06 * ( f(:,2,k)          - f(:,1,k) ) &
                                + b06 * ( f(:,3,k)          - f(:,this%n,k) )

               RHS(:,2:this%n-2,k)   = a06 * ( f(:,3:this%n-1,k) - f(:,2:this%n-2,k)) &
                               + b06 * ( f(:,4:this%n,k) - f(:,1:this%n-3,k) )

               RHS(:, this%n-1, k  ) = a06 * ( f(:,this%n,k)         -   f(:,this%n-1,k) ) &
                               + b06 * ( f(:,1,k)            - f(:,this%n-2,k) )

               RHS(:,this%n, k    ) = a06 * ( f(:,1,k)          - f(:,this%n,k) ) &
                                + b06 * ( f(:,2,k)          - f(:,this%n-1,k ) )

                !end do
            end do

          end select
        case (.FALSE.)
            RHS = 0.0d0
              
            a06 = a06d1 * this%onebydx
            b06 = b06d1 * this%onebydx

            a_np_3 = w3*q_pp * this%onebydx
            b_np_3 = w3*r_pp * this%onebydx

            a_np_2 = w2*p_p * this%onebydx

            a_np_1 = w1*(  p * this%onebydx)
            b_np_1 = w1*(  q * this%onebydx)
            c_np_1 = w1*(  r * this%onebydx)
            d_np_1 = w1*(  s * this%onebydx)

    
          select case(dir)

            case("N2F")
               RHS = 0.0d0
               a06 = a06d1 * this%onebydx
               b06 = b06d1 * this%onebydx

               a_np_2 = w2_N2F*p_p_N2F * this%onebydx
               b_np_2 = w2_N2F*q_p_N2F * this%onebydx

               a_np_1 = w1_N2F*(  p_N2F * this%onebydx)
               b_np_1 = w1_N2F*(  q_N2F * this%onebydx)
               c_np_1 = w1_N2F*(  r_N2F * this%onebydx)
               d_np_1 = w1_N2F*(  s_N2F * this%onebydx)

               do k = 1,n3

                     RHS(:,3:this%n-3,k) = a06 * ( f(:,4:this%n-2,k) - f(:,3:this%n-3,k) ) &
                                        + b06 * ( f(:,5:this%n-1,k)   - f(:,2:this%n-4,k) )
               end do

                select case(bc1)
                    case(1)
                    case(-1)
                    case(0)

                        do k = 1,n3

                             RHS(:,1,k) = a_np_1*f(:,1,k) + b_np_1*f(:,2,k) + c_np_1*f(:,3,k) + d_np_1*f(:,4,k)
                             RHS(:,2,k) = a_np_2*( f(:,3,k) - f(:,2,k) ) +    b_np_2*( f(:,4,k) - f(:,1,k))

                        end do

                  end select


                  select case(bcn)
                       case(1)
                       case(-1)
                       case(0)
 
                         
                          do k = 1,n3
           
                             RHS(:,this%n-1,k) = -a_np_1*f(:,this%n,k) - b_np_1*f(:,this%n-1,k) - c_np_1*f(:,this%n-2,k) - d_np_1*f(:,this%n-3,k)
                             RHS(:,this%n-2,k) = a_np_2*( f(:,this%n-1,k) - f(:,this%n-2,k) ) + b_np_2*( f(:,this%n,k) - f(:,this%n-3,k) )

                          end do
                  end select
 

            case("F2N")

              do k = 1, n3

                 RHS(:,4:this%n-3,k) = a06 * ( f(:,4:this%n-3,k) -f(:,3:this%n-4,k) ) &
                                        + b06 * ( f(:,5:this%n-2  ,k) -f(:,2:this%n-5,k) )
              end do

              select case(bc1)
                        ! left boundary (1:3)
                        case(1)
 
                        !!!!!!!!!!!!!!!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!!!!!!
                        case(-1)
                        !!!!!!!!!!!!!!!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!!!!!!
                        case(0)

                           do k =  1,n3
                               RHS(:,1         ,k) =   a_np_1* f(:,1         ,k) + b_np_1*f(:,2         ,k)   &
                                                       +   c_np_1* f(:,3         ,k) + d_np_1*f(:,4         ,k)

                               RHS(:,2         ,k) =   a_np_2*(f(:,2         ,k) - f(:,1         ,k))

                               RHS(:,3         ,k) =   a_np_3*(f(:,3         ,k) - f(:,2         ,k)) &
                                                       +   b_np_3*(f(:,4         ,k) - f(:,1         ,k))
                           end do
                    end select

                    select case(bcn)
                        !right boundary (n-2:n)
                        case(1)
                        !!!!!!!!!!!!!!!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!!!!!!
                        case(-1)
                        !!!!!!!!!!!!!!!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!!!!!!
                        case(0)

                           do k = 1, n3
                              RHS(:,this%n-2  ,k) =   a_np_3*(f(:,this%n-2  ,k) - f(:,this%n-3  ,k)) &
                                                      +   b_np_3*(f(:,this%n-1  ,k) - f(:,this%n-4  ,k))

                              RHS(:,this%n-1  ,k) =   a_np_2*(f(:,this%n - 1,k) - f(:,this%n-2  ,k))

                              RHS(:,this%n    ,k) =  -a_np_1* f(:,this%n-1  ,k) - b_np_1*f(:,this%n-2  ,k)   &
                                                     -   c_np_1* f(:,this%n-3  ,k) - d_np_1*f(:,this%n-4  ,k)
                           end do

                    end select


          end select
              
        end select
    
    end subroutine


    subroutine ComputeZD1RHS(this,f, RHS, dir,n1, n2,bc1,bcn) 
         
        class( cd06Stagger ), intent(in) :: this
        integer, intent(in) :: n1, n2,bc1,bcn
        real(rkind), dimension(n1,n2,this%n), intent(in) :: f
        real(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        real(rkind) :: a06, b06, a10,a102,a104,b10,b104,b102,c10,c102,c104,a101
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2, b_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1

        select case (this%periodic)
        case (.TRUE.)
            a06 = a06d1 * this%onebydx
            b06 = b06d1 * this%onebydx
            RHS = 0.0d0
         select case(dir)
            case("F2N")
            RHS(:,:,1         ) = a06 * ( f(:,:,1)          - f(:,:,this%n  ) ) &
                                + b06 * ( f(:,:,2)          - f(:,:,this%n-1) ) 
            
            RHS(:,:,2         ) = a06 * ( f(:,:,2)          - f(:,:,1       ) ) &
                               + b06 * ( f(:,:,3)          - f(:,:,this%n  ) )

            RHS(:,:,3:this%n-2) = a06 * ( f(:,:,3:this%n-2) - f(:,:,2:this%n-3) ) &
                               + b06 * ( f(:,:,4:this%n-1) - f(:,:,1:this%n-4) ) 
            
            RHS(:,:,this%n-1  ) = a06 * ( f(:,:,this%n-1)         - f(:,:,this%n-2) ) &
                               + b06 * ( f(:,:,this%n)          - f(:,:,this%n-3) ) 
            
            RHS(:,:,this%n    ) = a06 * ( f(:,:,this%n)          - f(:,:,this%n-1) ) &
                                + b06 * ( f(:,:,1)          - f(:,:,this%n-2) )
            case("N2F")

            !!!!! To Do      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            RHS(:,:,1         ) = a06 * ( f(:,:,2)          - f(:,:,1) ) &
                                + b06 * ( f(:,:,3)          - f(:,:,this%n) )

            RHS(:,:,2:this%n-2) = a06 * ( f(:,:,3:this%n-1) - f(:,:,2:this%n-2)) &
                               + b06 * ( f(:,:,4:this%n) - f(:,:,1:this%n-3) )

            RHS(:,:, this%n-1  ) = a06 * ( f(:,:, this%n)         - f(:,:,this%n-1) ) &
                               + b06 * ( f(:,:,1)            - f(:,:, this%n-2) )

            RHS(:,:,this%n    ) = a06 * ( f(:,:,1)          - f(:,:, this%n) ) &
                                + b06 * ( f(:,:,2)          - f(:,:, this%n-1 ) )
 
         end select
        case (.FALSE.)

            RHS = 0.0d0
            a06 = a06d1 * this%onebydx
            b06 = b06d1 * this%onebydx

            a_np_3 = w3*q_pp * this%onebydx
            b_np_3 = w3*r_pp * this%onebydx

            a_np_2 = w2*p_p * this%onebydx

            a_np_1 = w1*(  p * this%onebydx)
            b_np_1 = w1*(  q * this%onebydx)
            c_np_1 = w1*(  r * this%onebydx)
            d_np_1 = w1*(  s * this%onebydx)


            select case (dir)

               case ("N2F")
                 RHS = 0.0d0
                 a06 = a06d1 * this%onebydx
                 b06 = b06d1 * this%onebydx
                 a_np_2 = w2_N2F*p_p_N2F * this%onebydx
                 b_np_2 = w2_N2F*q_p_N2F * this%onebydx

                 a_np_1 = w1_N2F*(  p_N2F * this%onebydx)
                 b_np_1 = w1_N2F*(  q_N2F * this%onebydx)
                 c_np_1 = w1_N2F*(  r_N2F * this%onebydx)
                 d_np_1 = w1_N2F*(  s_N2F * this%onebydx)


                 RHS(:,:,3:this%n-3) = a06 * ( f(:,:,4:this%n-2) - f(:,:,3:this%n-3) ) &
                                        + b06 * ( f(:,:,5:this%n-1)   - f(:,:,2:this%n-4) )

                select case(bc1)
                    case(1)
                    case(-1)
                    case(0)

                        RHS(:,:,1) = a_np_1*f(:,:,1) + b_np_1*f(:,:,2) + c_np_1*f(:,:,3) + d_np_1*f(:,:,4)
                        RHS(:,:,2) = a_np_2*( f(:,:,3) - f(:,:,2) ) + b_np_2*( f(:,:,4) - f(:,:,1) )

                  end select

                  select case(bcn)
                       case(1)
                       case(-1)
                       case(0)

                         RHS(:,:,this%n-1) = -a_np_1*f(:,:,this%n) - b_np_1*f(:,:,this%n-1) - c_np_1*f(:,:,this%n-2) - d_np_1*f(:,:,this%n-3)
                         RHS(:,:,this%n-2) = a_np_2*( f(:,:,this%n-1) - f(:,:,this%n-2) ) + b_np_2*( f(:,:,this%n) - f(:,:,this%n-3) )

                  end select




               case ("F2N")!TODO: implement better non-periodic BC: currently12466...66421
                    !interior    
                    RHS(:,:,4:this%n-3) = a06 * ( f(:,:,4:this%n-3) -f(:,:,3:this%n-4) ) &
                               + b06 * ( f(:,:,5:this%n-2) - f(:,:,2:this%n-5) )

                    select case(bc1)
                        ! left boundary (1:3)
                        case(1)
                        !!!!!!!!!!!!!!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!
                        case(-1)
                        !!!!!!!!!!!!!!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!
                        case(0)
                             RHS(:,:,1) =   a_np_1* f(:,:,1) + b_np_1*f(:,:,2)   &
                                                     +   c_np_1* f(:,:,3) +  d_np_1*f(:,:,4)

                             RHS(:,:,2) =   a_np_2*(f(:,:,2) -         f(:,:,1))

                             RHS(:,:,3) =   a_np_3*(f(:,:,3) -         f(:,:,2)) &
                                +   b_np_3*(f(:,:,4) -         f(:,:,1))

                    end select

                    select case(bcn)
                        !right boundary (n-2:n)
                        case(1)
                        !!!!!!!!!!!!!!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!
                        case(-1)
                        !!!!!!!!!!!!!!!!!!!!!!! TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        case(0)
                             RHS(:,:,this%n-2  ) =   a_np_3*(f(:,:,this%n-2  )- f(:,:,this%n-3  )) &
                                +   b_np_3*(f(:,:,this%n-1  ) -f(:,:,this%n-4  ))

                             RHS(:,:,this%n-1  ) =   a_np_2*(f(:,:,this%n-1  ) -f(:,:,this%n-2  ))

                             RHS(:,:,this%n    ) =  -a_np_1* f(:,:,this%n-1  ) - b_np_1*f(:,:,this%n-2  )   &
                                                    -   c_np_1* f(:,:,this%n-3  ) - d_np_1*f(:,:,this%n-4  )

                    end select


            end select        
        end select
    
    end subroutine
   

    pure function ComputeD2RHS(this,f) result (RHS)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete / Needs to be updated
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        class( cd06Stagger ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: f
        real(rkind), dimension(this%n) :: RHS
        integer :: i
    
        select case (this%periodic)
        case (.TRUE.)
            RHS(1) = a06d2 * ( f(2) - two*f(1) + f(this%n)   ) &
                   + b06d2 * ( f(3) - two*f(1) + f(this%n-1) ) 
            RHS(2) = a06d2 * ( f(3) - two*f(2) + f(1)        ) &
                   + b06d2 * ( f(4) - two*f(2) + f(this%n)   ) 
            do i = 3,this%n-2
                RHS(i) = a06d2 * ( f(i+1) - two*f(i) + f(i-1) ) &
                       + b06d2 * ( f(i+2) - two*f(i) + f(i-2) ) 
            end do
            RHS(this%n-1) = a06d2 * ( f(this%n) - two*f(this%n-1) + f(this%n-2) ) &
                          + b06d2 * ( f(1)      - two*f(this%n-1) + f(this%n-3) ) 
            RHS(this%n)   = a06d2 * ( f(1)      - two*f(this%n)   + f(this%n-1) ) &
                          + b06d2 * ( f(2)      - two*f(this%n)   + f(this%n-2) )
        case (.FALSE.)
            RHS = zero
        end select
    
    end function

    subroutine dd1F2N(this, f, df, na, nb,bc1_,bcn_)
        class( cd06Stagger ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer  :: bc1,bcn
        real(rkind), dimension(this%n,na,nb), intent(in)  :: f
        real(rkind), dimension(this%n,na,nb), intent(out) :: df
        integer, optional, intent(in) :: bc1_, bcn_

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
                call GracefulExit("Incorrect boundary specification for bcn (should be 0, 1 or -1)", 324)
            end if
        else
            bcn = 0
        end if

        if (this%n == 1) then
            df = zero
            return
        end if

        call this%ComputeXD1RHS(f, df,"F2N", na, nb,bc1,bcn)
        select case (this%periodic) 
        case (.TRUE.)
           call this%SolveXLU1(df,na,nb)
        case (.FALSE.) 
           call this%SolveXTri1(df,na,nb)
        end select


    end subroutine

     subroutine dd1N2F(this, f, df, na, nb,bc1_,bcn_)
        class( cd06Stagger ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer  :: bc1,bcn
        real(rkind), dimension(this%n,na,nb), intent(in)  :: f
        real(rkind), dimension(this%n,na,nb), intent(out) :: df
        integer, optional, intent(in) :: bc1_, bcn_
        df = 0.0
        if(this%n == 1) then
            df = zero
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


        if (this%n == 1) then
            df = zero
            return
        end if
        call this%ComputeXD1RHS(f, df,"N2F", na, nb,bc1,bcn)
        select case (this%periodic)
        case (.TRUE.)
           call this%SolveXLU1(df,na,nb)
        case (.FALSE.)
           call this%SolveXTri1_N2F(df,na,nb)
        end select


    end subroutine

    subroutine dd2F2N(this, f, df, na, nb,bc1_,bcn_)
        class( cd06Stagger ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer  :: bc1,bcn
        real(rkind), dimension(na,this%n,nb), intent(in)  :: f
        real(rkind), dimension(na,this%n,nb), intent(out) :: df
        integer, optional, intent(in) :: bc1_, bcn_

        if(this%n == 1) then
            df = zero
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


        if (this%n == 1) then
            df = zero
            return
        end if

        call this%ComputeYD1RHS(f, df, "F2N", na, nb,bc1,bcn)
        select case (this%periodic) 
        case (.TRUE.)
           call this%SolveYLU1(df,na,nb)
        case (.FALSE.) 
           call this%SolveYTri1(df,na,nb)
        end select


    end subroutine
    
    subroutine dd2N2F(this, f, df, na, nb, bc1_, bcn_)
        class( cd06Stagger ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer  :: bc1,bcn
        real(rkind), dimension(na,this%n,nb), intent(in)  :: f
        real(rkind), dimension(na,this%n,nb), intent(out) :: df
        integer, optional, intent(in) :: bc1_, bcn_
        df = 0.0
        if(this%n == 1) then
            df = zero
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


        if (this%n == 1) then
            df = zero
            return
        end if

        call this%ComputeYD1RHS(f, df, "N2F", na, nb,bc1,bcn)
        select case (this%periodic)
        case (.TRUE.)
           call this%SolveYLU1(df,na,nb)
        case (.FALSE.)
           call this%SolveYTri1_N2F(df,na,nb)
        end select


    end subroutine

    subroutine dd3F2N(this, f, df, na, nb,bc1_,bcn_)
        class( cd06Stagger ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer  :: bc1,bcn
        real(rkind), dimension(na,nb,this%n), intent(in)  :: f
        real(rkind), dimension(na,nb,this%n), intent(out) :: df
        integer, optional, intent(in) :: bc1_, bcn_
        df = 0.0
        if(this%n == 1) then
            df = zero
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

        if (this%n == 1) then
            df = zero
            return
        end if

        call this%ComputeZD1RHS(f, df,"F2N", na, nb,bc1,bcn)
        select case (this%periodic) 
        case (.TRUE.)
           call this%SolveZLU1(df,na,nb)
        case (.FALSE.) 
           call this%SolveZTri1(df,na,nb)
        end select


    end subroutine

    subroutine dd3N2F(this, f, df, na, nb,bc1_,bcn_)
        class( cd06Stagger ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer  :: bc1,bcn
        real(rkind), dimension(na,nb,this%n), intent(in)  :: f
        real(rkind), dimension(na,nb,this%n), intent(out) :: df
        integer, optional, intent(in) :: bc1_, bcn_
        df = 0.0
        if(this%n == 1) then
            df = zero
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


        if (this%n == 1) then
            df = zero
            return
        end if

        call this%ComputeZD1RHS(f, df,"N2F", na, nb,bc1,bcn)
        select case (this%periodic)
        case (.TRUE.)
           call this%SolveZLU1(df,na,nb)
        case (.FALSE.)
           call this%SolveZTri1_N2F(df,na,nb)
        end select


    end subroutine

    
    function cd06Staggerder2(this, f) result(df)
        
        class( cd06Stagger ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: f
        real(rkind), dimension(this%n) :: df

        print*, f(1)
            df = zero
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete / Needs to be updated
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end function

end module
