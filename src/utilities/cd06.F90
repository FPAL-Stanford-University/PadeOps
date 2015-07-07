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
        real(rkind), allocatable, dimension(:,:) :: tri1
        real(rkind), allocatable, dimension(:,:) :: tri2

        contains

        procedure :: init

        procedure, private :: ComputeD1RHS
        procedure, private :: ComputeD2RHS

        procedure, private :: SolveD1
        procedure, private :: SolveD2

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

         ! Allocate 1st derivative Tri matrices.
         if(allocated( this%tri1 )) deallocate( this%tri1 ); allocate( this%tri1(n_,3) )
    
         ! Allocate 2nd derivative Tri matrices.
         if(allocated( this%tri2 )) deallocate( this%tri2 ); allocate( this%tri2(n_,3) )
    
        ! If everything passes
        ierr = 0
    
    end function
    
    subroutine ComputeLU(LU,n,b,d,a)
    
        integer, intent(in) :: n
        real(rkind), intent(in) :: d,a,b
        real(rkind), dimension(n,9), intent(out) :: LU
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
    
        end associate
    
    end subroutine
    
    subroutine SolveLU(LU,k,n)
    
        integer, intent(in) :: n
        real(rkind), dimension(n,5), intent(in)  :: LU
        real(rkind), dimension(n), intent(inout) :: k  ! Take in RHS and put solution into it
        integer :: i
    
        associate ( bc=>LU(:,1), h=>LU(:,2), onebyc=>LU(:,3), a=>LU(:,4), v=>LU(:,5) )
            
            ! Step 2
            do i = 2,n-1
                k(i) = k(i) - bc(i)*k(i-1)
            end do
            k(n) = k(n) - SUM( h(1:n-1)*k(1:n-1) )
    
            ! Step 3
            k(n) = k(n) * onebyc(n)
            k(n-1) = ( k(n-1) - v(n-1)*k(n) ) * onebyc(n-1)
            do i = n-2,1,-1
                k(i) = ( k(i) - a(i)*k(i+1) - v(i)*k(n) ) * onebyc(i)
            end do
    
        end associate
    
    end subroutine
    
    subroutine SolveTri(A,b,n)
    
        integer, intent(in) :: n
        real(rkind), dimension(n,3), intent(in) :: A
        real(rkind), dimension(n), intent(inout) :: b
    
    end subroutine
    
    pure function ComputeD1RHS(this,f) result (RHS)
    
        class( cd06 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: f
        real(rkind), dimension(this%n) :: RHS
        integer :: i
    
        select case (this%periodic)
        case (.TRUE.)
            RHS(1) = a06d1 * ( f(2)   - f(this%n)   ) &
                   + b06d1 * ( f(3)   - f(this%n-1) ) 
            RHS(2) = a06d1 * ( f(3)   - f(1)        ) &
                   + b06d1 * ( f(4)   - f(this%n)   ) 
            do i = 3,this%n-2
                RHS(i) = a06d1 * ( f(i+1) - f(i-1) ) &
                       + b06d1 * ( f(i+2) - f(i-2) ) 
            end do
            RHS(this%n-1) = a06d1 * ( f(this%n)   - f(this%n-2) ) &
                          + b06d1 * ( f(1)        - f(this%n-3) ) 
            RHS(this%n)   = a06d1 * ( f(1)        - f(this%n-1) ) &
                          + b06d1 * ( f(2)        - f(this%n-2) )
        case (.FALSE.)
            RHS = zero
        end select
    
    end function
    
    pure function ComputeD2RHS(this,f) result (RHS)
    
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

    subroutine SolveD1(this, rhs)

        class( cd06 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(inout) :: rhs
        
        select case (this%periodic)
        case(.TRUE.)
            call SolveLU(this%LU1,rhs,this%n)
        case(.FALSE.)
            call SolveTri(this%tri1,rhs,this%n)
        end select
    end subroutine

    subroutine SolveD2(this, rhs)

        class( cd06 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(inout) :: rhs
        
        select case (this%periodic)
        case(.TRUE.)
            call SolveLU(this%LU2,rhs,this%n)
        case(.FALSE.)
            call SolveTri(this%tri2,rhs,this%n)
        end select
    end subroutine

    function cd06der1(this, f) result(df)
        
        class( cd06 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: f
        real(rkind), dimension(this%n) :: df

        df = this%ComputeD1RHS(f)
        call this%SolveD1(df)
        df = df * this%onebydx

    end function

    function cd06der2(this, f) result(df)
        
        class( cd06 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: f
        real(rkind), dimension(this%n) :: df

        df = this%ComputeD2RHS(f)
        call this%SolveD2(df)
        df = df * this%onebydx2

    end function

end module
