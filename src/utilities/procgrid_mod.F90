module procgrid_mod
    use decomp_2d
    use mpi
    use kind_parameters, only: rkind
    implicit none

    private
    public :: procgrid, PAD_SPURIOUS, num_pad

    integer , parameter :: num_pad = 4 
    logical, parameter :: PAD_SPURIOUS = .false.
    
    type :: procgrid
        type(decomp_info) :: decomp 
        integer, dimension(3) :: xsz, xst, xen
        integer :: nrow, ncol, nx, ny, nz
        integer :: YneighLeft, YneighRight, ZneighDown, ZneighUp
        real(rkind), dimension(:,:,:), allocatable :: halo_sendbuff_yleft, halo_sendbuff_yright
        real(rkind), dimension(:,:,:), allocatable :: halo_sendbuff_zup  , halo_sendbuff_zdown
        real(rkind), dimension(:,:,:), allocatable :: halo_recvbuff_yleft, halo_recvbuff_yright
        real(rkind), dimension(:,:,:), allocatable :: halo_recvbuff_zup  , halo_recvbuff_zdown
        contains 
            procedure :: init
            procedure :: halo_exchange  
            procedure :: destroy
            procedure, private :: alloc_array3
            procedure, private :: alloc_array4
            procedure, private :: alloc_array5
            generic :: alloc_array => alloc_array3, alloc_array4, alloc_array5
    end type

contains
    
    subroutine alloc_array3(this, arr)
        class(procgrid), intent(in) :: this
        real(rkind), dimension(:,:,:), allocatable, intent(out) :: arr

        if (allocated(arr)) deallocate(arr)
        allocate(arr(-num_pad+1:this%decomp%xsz(1)+num_pad,-num_pad+1:this%decomp%xsz(2)+num_pad,-num_pad+1:this%decomp%xsz(3)+num_pad))
    end subroutine 

    subroutine alloc_array4(this, arr,nfields)
        class(procgrid), intent(in) :: this
        integer, intent(in) :: nfields
        real(rkind), dimension(:,:,:,:), allocatable, intent(out) :: arr

        if (allocated(arr)) deallocate(arr)
        allocate(arr(-num_pad+1:this%decomp%xsz(1)+num_pad,-num_pad+1:this%decomp%xsz(2)+num_pad,-num_pad+1:this%decomp%xsz(3)+num_pad, nfields))
    end subroutine 

    subroutine alloc_array5(this, arr,n1,n2)
        class(procgrid), intent(in) :: this
        integer, intent(in) :: n1,n2
        real(rkind), dimension(:,:,:,:,:), allocatable, intent(out) :: arr

        if (allocated(arr)) deallocate(arr)
        allocate(arr(-num_pad+1:this%decomp%xsz(1)+num_pad,-num_pad+1:this%decomp%xsz(2)+num_pad,-num_pad+1:this%decomp%xsz(3)+num_pad, n1, n2))
    end subroutine 
    
    subroutine init(this, prow, pcol, nx, ny, nz)
        use exits, only: message
        class(procgrid), intent(inout) :: this
        logical, dimension(3) :: periodicbcs
        integer, intent(in) :: prow, pcol, nx, ny, nz
        integer :: nrow, ncol, ntasks, ierr 
    
        integer :: DECOMP_CART_X
    
        ! Get processor decomposition
        call MPI_COMM_SIZE(MPI_COMM_WORLD,ntasks,ierr) 

        this%nx = nx
        this%ny = ny
        this%nz = nz
        if ((prow == 0) .and. (pcol == 0)) then
            call get_proc_factors(ntasks,nrow,ncol)
        else
            nrow = prow
            ncol = pcol
        end if 
        if (nrow*ncol .ne. ntasks) then
            print*, nrow, ncol, ntasks
            !call GracefulExit("Processor decomposition failed.", 123)
            print*, "Processor decomposition failed inside procgrid."
            call mpi_abort(MPI_COMM_WORLD, 0, ierr) 
        end if 

        if (.not. check_stencil_validity(this%ny,this%nz,nrow,ncol)) then
            !call message(0,"Default proc decomposition incompatible with stencil.")
            !call message(0,"Will try to swap order.")
            call swap_proc_order(nrow,ncol)
            if (.not. check_stencil_validity(this%ny,this%nz,nrow,ncol)) then
                print*, nrow, ncol
                print*, "Need a better processor decomposition algorithm for this domain. Current algorithm fails"
                call mpi_abort(MPI_COMM_WORLD, 0, ierr) 
            end if 
        end if         
    

        ! Set all BCs to true (will handle non-periodicity separately)
        periodicbcs(1) = .true.; periodicbcs(2) = .true.; periodicbcs(3) = .true.
        call decomp_info_init(this%nx, this%ny, this%nz, this%decomp)
        this%xsz = this%decomp%xsz
        this%xst = this%decomp%xst
        this%xen = this%decomp%xen

        ! Get neighbors
        call MPI_CART_CREATE(MPI_COMM_WORLD,2,[nrow,ncol],[.true.,.true.],.false.,DECOMP_CART_X,ierr)
        call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 0, 1, this%YneighLeft, this%YneighRight, ierr) 
        call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 1, 1, this%ZneighDown, this%ZneighUp, ierr) 

        this%nrow = nrow
        this%ncol = ncol 

        
        allocate(this%halo_sendbuff_yleft (this%decomp%xsz(1),num_pad,this%decomp%xsz(3)))
        allocate(this%halo_sendbuff_yright(this%decomp%xsz(1),num_pad,this%decomp%xsz(3)))
        allocate(this%halo_sendbuff_zup  (this%decomp%xsz(1),this%decomp%xsz(2)+2*num_pad,num_pad))
        allocate(this%halo_sendbuff_zdown(this%decomp%xsz(1),this%decomp%xsz(2)+2*num_pad,num_pad))
        
        allocate(this%halo_recvbuff_yleft (this%decomp%xsz(1),num_pad,this%decomp%xsz(3)))
        allocate(this%halo_recvbuff_yright(this%decomp%xsz(1),num_pad,this%decomp%xsz(3)))
        allocate(this%halo_recvbuff_zup  (this%decomp%xsz(1),this%decomp%xsz(2)+2*num_pad,num_pad))
        allocate(this%halo_recvbuff_zdown(this%decomp%xsz(1),this%decomp%xsz(2)+2*num_pad,num_pad))
        
        call message(0,"Processor decomposition information:")
        call message(1,"Nrows:", nrow)
        call message(1,"Ncols:", ncol)

    end subroutine 

    subroutine destroy(this)
        class(procgrid), intent(inout) :: this

        deallocate(this%halo_sendbuff_yleft)
        deallocate(this%halo_sendbuff_yright)
        deallocate(this%halo_recvbuff_yleft)
        deallocate(this%halo_recvbuff_yright)
        
        deallocate(this%halo_sendbuff_zup)
        deallocate(this%halo_sendbuff_zdown)
        deallocate(this%halo_recvbuff_zup)
        deallocate(this%halo_recvbuff_zdown)

    end subroutine 

    subroutine halo_exchange(this, ufield)
        use kind_parameters, only: mpirkind
        class(procgrid), intent(inout) :: this
        real(rkind), dimension(-num_pad+1:,-num_pad+1:,-num_pad+1:), intent(inout) :: ufield
        integer :: send_req1, recv_req1, status(MPI_STATUS_SIZE)
        integer :: send_req2, recv_req2, ierr, tag1 = 123, tag2 = 321
        integer :: nx, ny, nz, idx
        
        nx = this%decomp%xsz(1)
        ny = this%decomp%xsz(2)
        nz = this%decomp%xsz(3)

        ! Y - halo exchange
        
        call MPI_IRECV(this%halo_recvbuff_yleft ,nx*nz*num_pad,mpirkind,this%YneighLeft ,tag1,MPI_COMM_WORLD,recv_req1,ierr)
        call MPI_IRECV(this%halo_recvbuff_yright,nx*nz*num_pad,mpirkind,this%YneighRight,tag2,MPI_COMM_WORLD,recv_req2,ierr)
        
        this%halo_sendbuff_yleft  = ufield(1:nx,1:num_pad,1:nz)
        this%halo_sendbuff_yright = ufield(1:nx,ny-num_pad+1:ny,1:nz)
  
        call MPI_ISEND(this%halo_sendbuff_yleft ,nx*nz*num_pad,mpirkind,this%YneighLeft ,tag2,MPI_COMM_WORLD,send_req1,ierr)
        call MPI_ISEND(this%halo_sendbuff_yright,nx*nz*num_pad,mpirkind,this%YneighRight,tag1,MPI_COMM_WORLD,send_req2,ierr)

        call MPI_WAIT(recv_req1,status,ierr)
        call MPI_WAIT(recv_req2,status,ierr)

        ufield(1:nx,-num_pad+1:0,1:nz) = this%halo_recvbuff_yleft
        ufield(1:nx,ny+1:ny+num_pad,1:nz) = this%halo_recvbuff_yright
        
        call MPI_WAIT(send_req1,status,ierr)
        call MPI_WAIT(send_req2,status,ierr)

        ! Z - halo exchange 

        call MPI_IRECV(this%halo_recvbuff_zdown,(ny+2*num_pad)*nx*num_pad,mpirkind,this%ZneighDown,tag1,MPI_COMM_WORLD,recv_req1,ierr)
        call MPI_IRECV(this%halo_recvbuff_zup  ,(ny+2*num_pad)*nx*num_pad,mpirkind,this%ZneighUp  ,tag2,MPI_COMM_WORLD,recv_req2,ierr)
   
        this%halo_sendbuff_zdown  = ufield(1:nx,-num_pad+1:ny+num_pad,1:num_pad      )
        this%halo_sendbuff_zup    = ufield(1:nx,-num_pad+1:ny+num_pad,nz-num_pad+1:nz)
        
        call MPI_ISEND(this%halo_sendbuff_zdown,(ny+2*num_pad)*nx*num_pad,mpirkind,this%ZneighDown,tag2,MPI_COMM_WORLD,send_req1,ierr)
        call MPI_ISEND(this%halo_sendbuff_zup  ,(ny+2*num_pad)*nx*num_pad,mpirkind,this%ZneighUp  ,tag1,MPI_COMM_WORLD,send_req2,ierr)
        
        call MPI_WAIT(recv_req1,status,ierr)
        call MPI_WAIT(recv_req2,status,ierr)

        ufield(1:nx,-num_pad+1:ny+num_pad,-num_pad+1:0) = this%halo_recvbuff_zdown
        ufield(1:nx,-num_pad+1:ny+num_pad,nz+1:nz+num_pad) = this%halo_recvbuff_zup 
        
        call MPI_WAIT(send_req1,status,ierr)
        call MPI_WAIT(send_req2,status,ierr)
    
        ! X - exchange (swap)
        ufield(-num_pad+1:0   ,-num_pad+1:ny+num_pad,-num_pad+1:nz+num_pad) = ufield(nx-num_pad+1:nx,-num_pad+1:ny+num_pad,-num_pad+1:nz+num_pad )
        ufield(nx+1:nx+num_pad,-num_pad+1:ny+num_pad,-num_pad+1:nz+num_pad) = ufield(1:num_pad      ,-num_pad+1:ny+num_pad,-num_pad+1:nz+num_pad )
   
        if (PAD_SPURIOUS) then
            ufield(-num_pad+1,:,:) = 9.d99
            ufield(nx+num_pad,:,:) = 9.d99
            ufield(:,-num_pad+1,:) = 9.d99
            ufield(:,ny+num_pad,:) = 9.d99
            ufield(:,:,-num_pad+1) = 9.d99
            ufield(:,:,nz+num_pad) = 9.d99
        end if 
    
    
    end subroutine 


    pure function check_stencil_validity(nx,ny,nrow,ncol) result(isValid)
        integer, intent(in) :: nx, ny, nrow, ncol
        logical :: isValid

        if (((ny/ncol) < num_pad) .or. ((nx/nrow) < num_pad)) then
            isValid = .false. 
        else
            isValid = .true. 
        end if  
    end function

    subroutine swap_proc_order(nrow,ncol)
        integer, intent(inout) :: nrow, ncol
        integer :: tmp
        
        tmp = ncol
        ncol = nrow
        nrow = tmp
    end subroutine

    pure subroutine get_proc_factors(nproc, nrow, ncol)
        integer, intent(in) :: nproc
        integer, intent(out) :: ncol, nrow

        integer :: res

        nrow = floor(sqrt(real(nproc)))
        res = mod(nproc,nrow)

        do while ((res .ne. 0) .and. (nrow>1)) 
            nrow = nrow - 1
            res = mod(nproc,nrow)
        end do

        ncol = nproc/nrow

    end subroutine 
    
end module 
