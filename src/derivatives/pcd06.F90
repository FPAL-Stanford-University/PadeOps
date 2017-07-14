! Routines specific to 6th order Compact Finite Differencing scheme
! Periodic LU based on  Neossi Nguetchue, Abelman (Appl. Math. & Comp. 2008)

module pcd06stuff

    use mpi
    use kind_parameters, only: rkind, mpirkind
    use constants,       only: zero,one,two
    use t3dMod,          only: t3d
    use exits,           only: GracefulExit
    implicit none

    private
    public :: pcd06, alpha06d1, a06d1, b06d1
    
    ! 6th order first derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha06d1=  1.0_rkind / 3.0_rkind
    real(rkind), parameter :: a06d1    = (14.0_rkind / 9.0_rkind) / 2.0_rkind
    real(rkind), parameter :: b06d1    = ( 1.0_rkind / 9.0_rkind) / 4.0_rkind
    
    ! 6th order second derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha06d2=  2.0_rkind / 11.0_rkind
    real(rkind), parameter :: a06d2    = (12.0_rkind / 11.0_rkind) / 1.0_rkind
    real(rkind), parameter :: b06d2    = ( 3.0_rkind / 11.0_rkind) / 4.0_rkind
    
    type pcd06

        private

        type(t3d), pointer :: gp

        integer     :: n
        real(rkind) :: dx
        real(rkind) :: onebydx
        real(rkind) :: onebydx2

        logical     :: periodic = .TRUE.
        integer     :: bc1 = 0                             ! Boundary condition type. 0=Dirichlet, 1=Neumann
        integer     :: bcn = 0                             ! Boundary condition type. 0=Dirichlet, 1=Neumann

        real(rkind), allocatable, dimension(:,:) :: LU
        real(rkind), allocatable, dimension(:)   :: l_ii, u_ii
        real(rkind)                              :: l_iip1, u_iip1, e, f, g

        contains

        procedure :: init
        procedure :: destroy
        procedure :: GetSize
        procedure, private :: ComputeXD1RHS

        procedure, private :: ComputeLU
        procedure, private :: get_u_ii
        procedure, private :: get_l_ii

        procedure :: SolveLU
        procedure, private :: jacobi_iteration
        
    end type

contains

    pure function GetSize(this) result(val)
        class(pcd06), intent(in) :: this
        integer  :: val 
        val = this%n
    end function
    
    function init(this, gp_, dx_, periodic_, bc1_, bcn_) result(ierr)
    
        class( pcd06 ), intent(inout) :: this
        type(t3d), target, intent(in) :: gp_
        real(rkind), intent(in) :: dx_
        logical, intent(in) :: periodic_
        integer, intent(in) :: bc1_, bcn_
        integer :: ierr
    
        this%gp => gp_
        this%n  = this%gp%sz3D(1)
        this%dx = dx_
        this%onebydx = one/dx_
        this%onebydx2 = this%onebydx/dx_
    
        this%periodic = periodic_
    
        this%bc1 = bc1_
        this%bcn = bcn_
   
        if (periodic_) then 
            ! Allocate 1st derivative LU matrix.
            if(allocated( this%LU )) deallocate( this%LU ); allocate( this%LU(this%n-1,3) )
    
            ! Compute 1st derivative LU matrix
            call this%ComputeLU(this%LU,alpha06d1,one,alpha06d1)

            ! Allocate l_ii and u_ii vectors.
            if(allocated( this%l_ii )) deallocate( this%l_ii ); allocate( this%l_ii(this%n-1) )
            if(allocated( this%u_ii )) deallocate( this%u_ii ); allocate( this%u_ii(this%n-1) )
   
            ! Compute u_ii and l_ii
            call this%get_u_ii(this%LU, this%u_ii, alpha06d1)
            call this%get_l_ii(this%LU, this%l_ii, alpha06d1)

            ! Compute u_iip1 and l_iip1
            this%u_iip1 = alpha06d1                        ! Store only the single non-zero element
            this%l_iip1 = alpha06d1 / this%LU(this%n-1,2)  ! Store only the single non-zero element

            this%e = one - sum(this%l_ii*this%u_ii) - this%l_iip1*this%u_iip1
            this%f = - this%l_ii(this%n-1) * this%u_iip1
            this%g = - this%u_ii(this%n-1) * this%l_iip1
        else
            call GracefulExit("Non-periodic pcd06 yet to be implemented.",435)
        end if 
        
        ! if (this%gp%rank3D == 0) then
        !     print*, "LU: "
        !     print*, "    ", this%LU(:,1)
        !     print*, "    ", this%LU(:,2)
        !     print*, "    ", this%LU(:,3)
        !     print*, ""
        !     print*, "u_ii: "
        !     print*, "    ", this%u_ii
        !     print*, ""
        !     print*, "l_ii: "
        !     print*, "    ", this%l_ii
        !     print*, ""
        !     print*, "u_iip1: "
        !     print*, "    ", this%u_iip1
        !     print*, ""
        !     print*, "l_iip1: "
        !     print*, "    ", this%l_iip1
        !     print*, ""
        !     print*, "e, f, g: "
        !     print*, "    ", this%e, this%f, this%g
        !     print*, ""
        ! end if

        ! If everything passes
        ierr = 0
    
    end function

    subroutine destroy(this)
        class( pcd06 ), intent(inout) :: this

    end subroutine
    
    subroutine ComputeLU(this, LU, b, a, c)
        class (pcd06), intent(inout) :: this
        real(rkind), dimension(this%n-1,3), intent(out) :: LU
        real(rkind), intent(in) :: a, b, c
        integer :: i

        LU = 0._rkind

        LU(1,2) = a
        do i = 2,this%n-1
            LU(i,2)   = a - b*c/LU(i-1,2)
            LU(i-1,3) = c
            LU(i,1)   = b/LU(i-1,2)
        end do
    end subroutine
    
    subroutine get_u_ii(this, LU, u_ii, b)
        class (pcd06), intent(inout) :: this
        real(rkind), dimension(this%n-1,3), intent(in) :: LU
        real(rkind), dimension(this%n-1), intent(out) :: u_ii
        real(rkind), intent(in) :: b
        integer :: i

        ! Solve L u_ii = [b, 0, ..., 0]^T
        u_ii(1) = b
        do i = 2,this%n-1
            u_ii(i) = - LU(i,1) * u_ii(i-1)
        end do
    end subroutine

    subroutine get_l_ii(this, LU, l_ii, c)
        class (pcd06), intent(inout) :: this
        real(rkind), dimension(this%n-1,3), intent(in) :: LU
        real(rkind), dimension(this%n-1), intent(out) :: l_ii
        real(rkind), intent(in) :: c
        integer :: i

        ! Solve U^T l_ii = [c, 0, ..., 0]^T
        l_ii(1) = c
        do i = 2,this%n-1
            l_ii(i) = - LU(i-1,3) * l_ii(i-1) / LU(i,2)
        end do
    end subroutine
    
    subroutine ComputeXD1RHS(this,f, RHS, n2, n3) 
        class( pcd06 ), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        integer ::  j, k
        real(rkind) :: a06, b06

        a06 = a06d1 * this%onebydx
        b06 = b06d1 * this%onebydx
        
        do k = 1,n3
            do j = 1,n2
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
    
    end subroutine
   
    subroutine jacobi_iteration(this, df, n2, n3, niters, recvbuf_l, recvbuf_r)
        class(pcd06), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(0:this%n+1,n2,n3), intent(inout) :: df
        real(rkind), dimension(n2,n3), intent(out) :: recvbuf_l, recvbuf_r
        integer, intent(in) ::  niters
        integer ::  iters, ierr
        real(rkind), dimension(n2,n3) :: sendbuf
        integer :: recv_request_left, recv_request_right
        integer :: send_request_left, send_request_right
        integer, dimension(MPI_STATUS_SIZE) :: status

        df(1,:,:) = df(1,:,:) / this%e
        sendbuf = df(1,:,:)

        do iters = 1,niters
            ! Communicate updated iterations of the reduced system
            call mpi_irecv( recvbuf_l, n2*n3, mpirkind, this%gp%xleft,  0, this%gp%commX, recv_request_left,  ierr)
            call mpi_irecv( recvbuf_r, n2*n3, mpirkind, this%gp%xright, 1, this%gp%commX, recv_request_right, ierr)

            call mpi_isend( sendbuf, n2*n3, mpirkind, this%gp%xright, 0, this%gp%commX, send_request_right, ierr)
            call mpi_isend( sendbuf, n2*n3, mpirkind, this%gp%xleft,  1, this%gp%commX, send_request_left,  ierr)

            call mpi_wait(recv_request_left,  status, ierr)
            call mpi_wait(recv_request_right, status, ierr)
            call mpi_wait(send_request_left,  status, ierr)
            call mpi_wait(send_request_right, status, ierr)

            sendbuf = -(recvbuf_l*this%g + recvbuf_r*this%f)/this%e + df(1,:,:)
        end do

        ! Communicate updated iterations of the reduced system
        call mpi_irecv( recvbuf_l, n2*n3, mpirkind, this%gp%xleft,  0, this%gp%commX, recv_request_left,  ierr)
        call mpi_irecv( recvbuf_r, n2*n3, mpirkind, this%gp%xright, 1, this%gp%commX, recv_request_right, ierr)

        call mpi_isend( sendbuf, n2*n3, mpirkind, this%gp%xright, 0, this%gp%commX, send_request_right, ierr)
        call mpi_isend( sendbuf, n2*n3, mpirkind, this%gp%xleft,  1, this%gp%commX, send_request_left,  ierr)

        call mpi_wait(recv_request_left,  status, ierr)
        call mpi_wait(recv_request_right, status, ierr)
        call mpi_wait(send_request_left,  status, ierr)
        call mpi_wait(send_request_right, status, ierr)

        df(1,:,:) = sendbuf

    end subroutine

    subroutine SolveLU(this, df, n2, n3, niters)
        class( pcd06 ), intent(in) :: this
        integer, intent(in) :: n2, n3, niters
        real(rkind), dimension(0:this%n+1,n2,n3), intent(inout) :: df
        real(rkind), dimension(n2,n3) :: recvbuf_l, recvbuf_r
        integer ::  i, j, k, ierr

        real(rkind), dimension(n2,n3) :: sendbuf
        integer :: recv_request_left
        integer :: send_request_right
        integer, dimension(MPI_STATUS_SIZE) :: status

        ! L y_i = rhs_i
        do k = 1,n3
            do j = 1,n2
                do i = 3,this%n
                    df(i,j,k) = df(i,j,k) - this%LU(i-1,1)*df(i-1,j,k)
                end do
            end do
        end do

        ! Communicate last cell data to the processor to the right and recv from left
        call mpi_irecv( recvbuf_l, n2*n3, mpirkind, this%gp%xleft,  0, this%gp%commX, recv_request_left,  ierr)
        sendbuf = df(this%n,:,:)
        call mpi_isend( sendbuf, n2*n3, mpirkind, this%gp%xright, 0, this%gp%commX, send_request_right, ierr)
        call mpi_wait(recv_request_left,  status, ierr)
        call mpi_wait(send_request_right, status, ierr)

        ! get y_r
        do k = 1,n3
            do j = 1,n2
                df(1,j,k) = df(i,j,k) - sum(this%l_ii*df(2:this%n,j,k)) - this%l_iip1*recvbuf_l(j,k)
            end do
        end do

        ! Solve the reduced system
        call this%jacobi_iteration( df, n2, n3, niters, recvbuf_l, recvbuf_r )


        ! Solve for the rest of the solution vector
        do k = 1,n3
            do j = 1,n2
                df(2:this%n,j,k) = df(2:this%n,j,k) - df(1,j,k)*this%u_ii 
                df(this%n,j,k) = df(this%n,j,k) - recvbuf_r(j,k)*this%u_iip1

                df(this%n,j,k) = df(this%n,j,k) / this%LU(this%n-1,2)
                do i = this%n-1,2,-1
                    df(i,j,k) = (df(i,j,k) - this%LU(i-1,3)*df(i+1,j,k)) / this%LU(i-1,2)
                end do

            end do
        end do

    end subroutine

end module
