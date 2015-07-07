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

        procedure, private :: ComputeD1RHS
        procedure, private :: ComputeD2RHS

        procedure, private :: SolveD1
        procedure, private :: SolveD2

        procedure :: cd10der1
        procedure :: cd10der2

    end type

contains

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
        
        ! Allocate 1st derivative Penta matrices.
        if(allocated( this%penta1 )) deallocate( this%penta1 ); allocate( this%penta1(n_,5) )
    
        ! Allocate 2nd derivative Penta matrices.
        if(allocated( this%penta2 )) deallocate( this%penta2 ); allocate( this%penta2(n_,5) )
    
        ! If everything passes
        ierr = 0
    
    end function
    
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
    
    subroutine SolveLU(LU,y,n)
    
        integer, intent(in) :: n
        real(rkind), dimension(n,9), intent(in)  :: LU
        real(rkind), dimension(n), intent(inout) :: y  ! Take in RHS and put solution into it
        integer :: i
    
        associate( b=>LU(:,1), eg=>LU(:,2), k=>LU(:,3),&
                   l=>LU(:,4), onebyg=>LU(:,5), h=>LU(:,6),&
                   f=>LU(:,7),  v=>LU(:,8), w=>LU(:,9))
    
            ! Step 8 ( update y instead of creating z )
            y(2) = y(2) - b(2)*y(1)
    
            ! Step 9
            do i = 3,n-2
                y(i) = y(i) - b(i)*y(i-1) - eg(i)*y(i-2)
            end do
    
            ! Step 10
            y(n-1) = y(n-1) - SUM( k(1:n-2)*y(1:n-2) )
            y(n)   = y(n)   - SUM( l(1:n-1)*y(1:n-1) )
    
            ! Step 11
            y(n) = y(n) * onebyg(n)
            y(n-1) = ( y(n-1) - w(n-1)*y(n) ) * onebyg(n-1)
            y(n-2) = ( y(n-2) - v(n-2)*y(n-1) - w(n-2)*y(n) ) * onebyg(n-2)
            y(n-3) = ( y(n-3) - h(n-3)*y(n-2) - v(n-3)*y(n-1) - w(n-3)*y(n) ) * onebyg(n-3)
            do i = n-4,1,-1
                y(i) = ( y(i) - h(i)*y(i+1) - f(i)*y(i+2) - v(i)*y(n-1) - w(i)*y(n) ) * onebyg(i)
            end do
    
        end associate
            
    end subroutine
    
    subroutine SolvePenta(A,b,n)

        integer, intent(in) :: n
        real(rkind), dimension(n,5), intent(in) :: A
        real(rkind), dimension(n), intent(inout) :: b

    end subroutine

    pure function ComputeD1RHS(this, f) result (RHS)
    
        class( cd10 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: f
        real(rkind), dimension(this%n) :: RHS
        integer :: i

        if(this%n == 1) then
            RHS = zero
            return
        end if
    
        select case (this%periodic)
        case (.TRUE.)
            RHS(1) = a10d1 * ( f(2)   - f(this%n)   ) &
                   + b10d1 * ( f(3)   - f(this%n-1) ) &
                   + c10d1 * ( f(4)   - f(this%n-2) )
            RHS(2) = a10d1 * ( f(3)   - f(1)        ) &
                   + b10d1 * ( f(4)   - f(this%n)   ) &
                   + c10d1 * ( f(5)   - f(this%n-1) )
            RHS(3) = a10d1 * ( f(4)   - f(2)        ) &
                   + b10d1 * ( f(5)   - f(1)        ) &
                   + c10d1 * ( f(6)   - f(this%n)   )
            do i = 4,this%n-3
                RHS(i) = a10d1 * ( f(i+1) - f(i-1) ) &
                       + b10d1 * ( f(i+2) - f(i-2) ) &
                       + c10d1 * ( f(i+3) - f(i-3) )
            end do
            RHS(this%n-2) = a10d1 * ( f(this%n-1) - f(this%n-3) ) &
                          + b10d1 * ( f(this%n)   - f(this%n-4) ) &
                          + c10d1 * ( f(1)        - f(this%n-5) )
            RHS(this%n-1) = a10d1 * ( f(this%n)   - f(this%n-2) ) &
                          + b10d1 * ( f(1)        - f(this%n-3) ) &
                          + c10d1 * ( f(2)        - f(this%n-4) )
            RHS(this%n)   = a10d1 * ( f(1)        - f(this%n-1) ) &
                          + b10d1 * ( f(2)        - f(this%n-2) ) &
                          + c10d1 * ( f(3)        - f(this%n-3) )
        case (.FALSE.)
            RHS = zero
        end select
    
    end function
    
    pure function ComputeD2RHS(this,f) result (RHS)
    
        class( cd10 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: f
        real(rkind), dimension(this%n) :: RHS
        integer :: i
    
        if(this%n == 1) then
            RHS = zero
            return
        end if
    
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
        case (.FALSE.)
            RHS = zero
        end select
    
    end function

    subroutine SolveD1(this, rhs)

        class( cd10 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(inout) :: rhs

        select case (this%periodic)
        case(.TRUE.)
            call SolveLU(this%LU1,rhs,this%n)
        case(.FALSE.)
            call SolvePenta(this%penta1,rhs,this%n)
        end select

    end subroutine 

    subroutine SolveD2(this, rhs)

        class( cd10 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(inout) :: rhs

        select case (this%periodic)
        case(.TRUE.)
            call SolveLU(this%LU2,rhs,this%n)
        case(.FALSE.)
            call SolvePenta(this%penta2,rhs,this%n)
        end select

    end subroutine

    function cd10der1(this, f) result(df)

        class( cd10 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: f
        real(rkind), dimension(this%n) :: df

        df = this%ComputeD1RHS(f)
        call this%SolveD1(df)
        df = df * this%onebydx

    end function

    function cd10der2(this, f) result(df)

        class( cd10 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: f
        real(rkind), dimension(this%n) :: df

        df = this%ComputeD2RHS(f)
        call this%SolveD2(df)
        df = df * this%onebydx2

    end function

end module
