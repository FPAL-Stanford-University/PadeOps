! Routines specific to 6th order Compact Finite Differencing scheme
! Periodic LU based on  Neossi Nguetchue, Abelman (Appl. Math. & Comp. 2008)

module cd06stuff

    use kind_parameters, only: rkind
    use constants,       only : zero,one,two
    implicit none

    private
    public :: cd06
    
    ! 6th order first derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha06d1=  1.0_rkind / 3.0_rkind
    real(rkind), parameter :: a06d1    = (14.0_rkind / 9.0_rkind) / 2.0_rkind
    real(rkind), parameter :: b06d1    = ( 1.0_rkind / 9.0_rkind) / 4.0_rkind
    
    ! 6th order second derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha06d2=  2.0_rkind / 11.0_rkind
    real(rkind), parameter :: a06d2    = (12.0_rkind / 11.0_rkind) / 1.0_rkind
    real(rkind), parameter :: b06d2    = ( 3.0_rkind / 11.0_rkind) / 4.0_rkind
    
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
    
    
    type cd06

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
        real(rkind), allocatable, dimension(:,:) :: Tri1
        real(rkind), allocatable, dimension(:,:) :: Tri2

        contains

        procedure :: init
        procedure :: destroy

        procedure, private :: ComputeD1RHS
        procedure, private :: ComputeD2RHS

        procedure, private :: SolveLU1
        procedure, private :: SolveLU2
       
        procedure, private :: ComputeTri1
        procedure, private :: ComputeTri2

        procedure, private :: SolveTri1
        procedure, private :: SolveTri2
        
        procedure :: cd06der1
        procedure :: cd06der2

    end type

contains

    function init(this, n_, dx_, periodic_, bc1_, bcn_) result(ierr)
    
        class( cd06 ), intent(inout) :: this
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
  
            if (n_ .GE. 8) then             
                call this%ComputeTri1(bc1_,bcn_)
            else if (n_ .EQ. 1) then
                this%Tri1 = one
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

        class( cd06 ), intent(inout) :: this

        ! Dellocate 1st derivative LU matrix.
        if(allocated( this%LU1 )) deallocate( this%LU1 )
    
        ! Dellocate 2nd derivative LU matrix.
        if(allocated( this%LU2 )) deallocate( this%LU2 )

        ! Dellocate 1st derivative tri matrix.
        if(allocated( this%tri1 )) deallocate( this%tri1 )
    
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
        class (cd06), intent(inout) :: this
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
    
    
    subroutine ComputeTri2(this,bc1,bcn)
        class (cd06), intent(inout) :: this
        integer, intent(in) :: bc1, bcn
   
        

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        if ((bc1 == 1) .and. (bcn == 1)) then
            this%tri2 = zero
        end if 
    end subroutine 
    
    subroutine SolveLU1(this,y,n2,n3, y1, yn, z1, zn)
        
        class (cd06), intent(in) :: this
        integer, intent(in) :: n2,n3,y1,yn,z1,zn
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y  ! Take in RHS and put solution into it
        integer :: i, j, k
        real(rkind) :: sum1 

        do k = z1,zn
            do j = y1,yn 
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

    subroutine SolveTri1(this,y,n2,n3,y1,yn,z1,zn)
    
        class (cd06), intent(in) :: this
        integer, intent(in) :: n2,n3,y1,yn,z1,zn
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y  
        integer :: i, j, k
  
        do k = z1,zn 
            do j = y1,yn
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
   

    subroutine SolveLU2(this,y,n2,n3, y1, yn, z1, zn)
        
        class (cd06), intent(in) :: this
        integer, intent(in) :: n2,n3,y1,yn,z1,zn
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y  ! Take in RHS and put solution into it

        print*, n2,n3,y1,yn,z1,zn
        y = zero
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine 
     
    subroutine SolveTri2(this,y,n2,n3,y1,yn,z1,zn)
    
        class (cd06), intent(in) :: this
        integer, intent(in) :: n2,n3,y1,yn,z1,zn
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: y  
  
        print*, n2,n3,y1,yn,z1,zn 
        y = zero 
    end subroutine
    
    subroutine ComputeD1RHS(this,f, RHS, n2, n3, y1, yn, z1, zn) 
         
        class( cd06 ), intent(in) :: this
        integer, intent(in) :: n2, n3, y1, yn, z1, zn
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        integer ::  j, k
        real(rkind) :: a06, b06
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1

        select case (this%periodic)
        case (.TRUE.)
            a06 = a06d1 * this%onebydx
            b06 = b06d1 * this%onebydx
           
            do k = z1,zn
                do j = y1,yn 
                    RHS(1         ,j,k) = a06 * ( f(2,j,k)          - f(this%n  ,j,k) ) &
                                        + b06 * ( f(3,j,k)          - f(this%n-1,j,k) ) 
                    
                    RHS(2         ,j,k) = a06 * ( f(3,j,k)          - f(1       ,j,k) ) &
                                        + b06 * ( f(4,j,k)          - f(this%n  ,j,k) )

                    RHS(3:this%n-2,j,k) = a06 * ( f(4:this%n-1,j,k) - f(2:this%n-3,j,k) ) &
                                        + b06 * ( f(5:this%n  ,j,k) - f(1:this%n-4,j,k) ) 
                    
                    RHS(this%n-1  ,j,k) = a06 * ( f(this%n,j,k)         - f(this%n-2,j,k) ) &
                                        + b06 * ( f(1,j,k)          - f(this%n-3,j,k) ) 
                    
                    RHS(this%n    ,j,k) = a06 * ( f(1,j,k)          - f(this%n-1,j,k) ) &
                                        + b06 * ( f(2,j,k)          - f(this%n-2,j,k) )
                end do 
            end do 
        case (.FALSE.)
            a06 = a06d1 * this%onebydx
            b06 = b06d1 * this%onebydx
        
            a_np_3 = w3*q_pp * this%onebydx  
            b_np_3 = w3*r_pp * this%onebydx 

            a_np_2 = w2*q_p * this%onebydx
            
            a_np_1 = w1*( -p * this%onebydx)
            b_np_1 = w1*(  q * this%onebydx)
            c_np_1 = w1*(  r * this%onebydx)
            d_np_1 = w1*(  s * this%onebydx)
       
            do k = z1,zn
                do j = y1,yn
                    RHS(1         ,j,k) =   a_np_1* f(1         ,j,k) +  b_np_1*f(2         ,j,k)   &
                                        +   c_np_1* f(3         ,j,k) +  d_np_1*f(4         ,j,k) 
                   
                    RHS(2         ,j,k) =   a_np_2*(f(3         ,j,k) -         f(1         ,j,k))
                    
                    RHS(3         ,j,k) =   a_np_3*(f(4         ,j,k) -         f(2         ,j,k)) &
                                        +   b_np_3*(f(5         ,j,k) -         f(1         ,j,k)) 

                    RHS(4:this%n-3,j,k) =   b06   *(f(6:this%n-1,j,k) -         f(2:this%n-5,j,k)) &
                                        +   a06   *(f(5:this%n-2,j,k) -         f(3:this%n-4,j,k)) 
                                          

                    RHS(this%n-2  ,j,k) =   a_np_3*(f(this%n-1  ,j,k) -         f(this%n-3  ,j,k)) &
                                        +   b_np_3*(f(this%n    ,j,k) -         f(this%n-4  ,j,k)) 
                    
                    RHS(this%n-1  ,j,k) =   a_np_2*(f(this%n    ,j,k) -         f(this%n-2  ,j,k))

                    RHS(this%n    ,j,k) =  -a_np_1* f(this%n    ,j,k) -  b_np_1*f(this%n-1  ,j,k)   &
                                        -   c_np_1* f(this%n-2  ,j,k) -  d_np_1*f(this%n-3  ,j,k)

                end do 
            end do  
        
        end select
    
    end subroutine
    
    pure function ComputeD2RHS(this,f) result (RHS)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete / Needs to be updated
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        class( cd06 ), intent(in) :: this
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

    subroutine cd06der1(this, f, df, n2, n3, y1, yn, z1, zn)
        class( cd06 ), intent(in) :: this
        integer, intent(in) :: n2, n3, y1, yn, z1, zn
        real(rkind), dimension(this%n), intent(in)  :: f
        real(rkind), dimension(this%n), intent(out) :: df

        if (this%n == 1) then
            df = zero
            return
        end if

        call this%ComputeD1RHS(f, df, n2, n3, y1, yn, z1, zn)

        select case (this%periodic) 
        case (.TRUE.)
           call this%SolveLU1(df,n2,n3,y1,yn,z1,zn)
        case (.FALSE.) 
           call this%SolveTri1(df,n2,n3, y1, yn, z1, zn)
        end select

    end subroutine

    function cd06der2(this, f) result(df)
        
        class( cd06 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: f
        real(rkind), dimension(this%n) :: df

        print*, f(1)
            df = zero
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                ! Incomplete / Needs to be updated
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end function

end module
