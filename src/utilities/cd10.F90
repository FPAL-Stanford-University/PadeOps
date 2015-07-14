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
        procedure :: destroy

        procedure, private :: ComputeD1RHS
        procedure, private :: ComputeD2RHS

        procedure, private :: SolveLU1
        procedure, private :: SolveLU2
        
        procedure, private :: SolvePenta1
        procedure, private :: SolvePenta2
        
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
    
    subroutine SolveLU1(this,y,n2,n3, y1, yn, z1, zn)
    
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
                y(this%n-1,j,k) = y(this%n-1,j,k) - sum1 !SUM( this%LU1(1:this%n-2,3)*y(1:this%n-2) )
                y(this%n,j,k)   = ( y(this%n,j,k)   - sum2 - this%LU1(this%n-1,4)*y(this%n-1,j,k) ) * this%LU1(this%n,5) !SUM( this%LU1(1:this%n-1,4)*y(1:this%n-1) )
    
                ! Step 11
                !y(this%n) = y(this%n) * this%LU1(this%n,5)
                y(this%n-1,j,k) = ( y(this%n-1,j,k) - this%LU1(this%n-1,9)*y(this%n,j,k) ) * this%LU1(this%n-1,5)
                y(this%n-2,j,k) = ( y(this%n-2,j,k) - this%LU1(this%n-2,8)*y(this%n-1,j,k) - this%LU1(this%n-2,9)*y(this%n,j,k) ) * this%LU1(this%n-2,5)
                y(this%n-3,j,k) = ( y(this%n-3,j,k) - this%LU1(this%n-3,6)*y(this%n-2,j,k) - this%LU1(this%n-3,8)*y(this%n-1,j,k) - this%LU1(this%n-3,9)*y(this%n,j,k) ) * this%LU1(this%n-3,5)
                do i = this%n-4,1,-1
                    y(i,j,k) = ( y(i,j,k) - this%LU1(i,6)*y(i+1,j,k) - this%LU1(i,7)*y(i+2,j,k) - this%LU1(i,8)*y(this%n-1,j,k) - this%LU1(i,9)*y(this%n,j,k) ) * this%LU1(i,5)
                end do
            end do
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
    
    end subroutine
    
    subroutine SolvePenta1(this,b,n2,n3)

        class( cd10 ), intent(in) :: this
        integer, intent(in) :: n2,n3
        real(rkind), dimension(this%n,n2,n3), intent(inout) :: b

    end subroutine

    subroutine SolvePenta2(this,b)

        class( cd10 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(inout) :: b

    end subroutine

    pure subroutine ComputeD1RHS(this, f, RHS, n2, n3, y1, yn, z1, zn) !result (RHS)
    
        class( cd10 ), intent(in) :: this
        integer, intent(in) :: n2, n3, y1, yn, z1, zn
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        real(rkind) :: a10,b10,c10
        integer :: i,j,k

        if(this%n == 1) then
            RHS = zero
            return
        end if

        a10 = a10d1 * this%onebydx
        b10 = b10d1 * this%onebydx
        c10 = c10d1 * this%onebydx
    
        select case (this%periodic)
        case (.TRUE.)
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
                    ! do i = 4,this%n-3
                    !     RHS(i,j,k) = a10 * ( f(i+1,j,k) - f(i-1,j,k) ) &
                    !            + b10 * ( f(i+2,j,k) - f(i-2,j,k) ) &
                    !            + c10 * ( f(i+3,j,k) - f(i-3,j,k) )
                    ! end do
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
            RHS = zero
        end select
    
    end subroutine
    
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

    subroutine SolveD2(this, rhs)

        class( cd10 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(inout) :: rhs

        select case (this%periodic)
        case(.TRUE.)
            call this%SolveLU2(rhs)
        case(.FALSE.)
            call this%SolvePenta2(rhs)
        end select

    end subroutine

    subroutine cd10der1(this, f, df, n2, n3, y1, yn, z1, zn) !result(df)

        class( cd10 ), intent(in) :: this
        integer, intent(in) :: n2, n3, y1, yn, z1, zn
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: df

        call this%ComputeD1RHS(f, df, n2, n3, y1, yn, z1, zn)
        
        select case (this%periodic)
        case(.TRUE.)
            call this%SolveLU1(df,n2,n3, y1, yn, z1, zn)
        case(.FALSE.)
            call this%SolvePenta1(df,n2,n3)
        end select

    end subroutine

    function cd10der2(this, f) result(df)

        class( cd10 ), intent(in) :: this
        real(rkind), dimension(this%n), intent(in) :: f
        real(rkind), dimension(this%n) :: df

        df = this%ComputeD2RHS(f)
        call this%SolveD2(df)
        df = df * this%onebydx2

    end function

end module
