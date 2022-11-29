module procgrid_mod
    use decomp_2d
    use mpi
    use kind_parameters, only: rkind
    implicit none

    private
    public :: procgrid, PAD_SPURIOUS

    integer :: num_pad = 1 
    logical, parameter :: PAD_SPURIOUS = .false.
    
    type :: procgrid
        type(decomp_info) :: decomp 
        integer, dimension(3) :: zsz, zst, zen
        integer :: nrow, ncol, nx, ny, nz
        integer :: XneighLeft, XneighRight, YneighDown, YneighUp
        real(rkind), dimension(:,:,:), allocatable :: halo_sendbuff_xleft, halo_sendbuff_xright
        real(rkind), dimension(:,:,:), allocatable :: halo_sendbuff_yup  , halo_sendbuff_ydown
        real(rkind), dimension(:,:,:), allocatable :: halo_recvbuff_xleft, halo_recvbuff_xright
        real(rkind), dimension(:,:,:), allocatable :: halo_recvbuff_yup  , halo_recvbuff_ydown
        logical :: have_x_loEdge
        logical :: have_y_loEdge
        logical :: have_z_loEdge
        logical :: have_x_hiEdge
        logical :: have_y_hiEdge
        logical :: have_z_hiEdge
        contains 
            procedure :: init
            procedure :: halo_exchange  
            procedure :: halo_exchange_x 
            procedure :: halo_exchange_y 
            procedure :: halo_exchange_z 
            procedure :: destroy
            procedure, private :: alloc_array3
            procedure, private :: alloc_array4
            procedure, private :: alloc_array5
            generic :: alloc_array => alloc_array3, alloc_array4, alloc_array5
            procedure :: get_domain_avg
            procedure :: get_domain_avg_intarray 
            procedure :: print_error_fields
            procedure :: print_error_fields_interior 
            procedure :: get_min_interior 
            procedure :: get_max_interior 
    end type

contains
    
    function get_domain_avg(this, arr) result(avg)
        use reductions, only: p_sum 
        class(procgrid), intent(in) :: this
        real(rkind), dimension(-num_pad+1:,-num_pad+1:,-num_pad+1:), intent(inout) :: arr
        real(rkind) :: domsum, avg
        
        domsum = p_sum(sum(arr(1:this%zsz(1),1:this%zsz(2),1:this%zsz(3))))
        avg = domsum/(this%nx*this%ny*this%nz)

    end function 
    
    function get_domain_avg_intarray(this, arr) result(avg)
        use reductions, only: p_sum 
        class(procgrid), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(inout) :: arr
        real(rkind) :: domsum, avg
        
        domsum = p_sum(sum(arr(1:this%zsz(1),1:this%zsz(2),1:this%zsz(3))))
        avg = domsum/(this%nx*this%ny*this%nz)

    end function 

    subroutine print_error_fields(this, f1,f2, label, errf) 
        use reductions, only: p_maxval 
        use decomp_2d, only: nrank 
        class(procgrid), intent(in) :: this
        real(rkind), dimension(-num_pad+1:,-num_pad+1:,-num_pad+1:), intent(in) :: f1, f2
        character(len=*), intent(in) :: label
        real(rkind), intent(out), optional :: errf
        real(rkind) :: err

        err = p_maxval(maxval(abs(f1 - f2)))
        if (nrank == 0) then
            print*, "====================="
            print*, "Error in field: ", label
            print*, ">> ", err
            print*, "====================="
        end if 

        if (present(errf)) errf = err

    end subroutine
    
    subroutine print_error_fields_interior(this, f1,f2, label, errf) 
        use reductions, only: p_maxval 
        use decomp_2d, only: nrank 
        class(procgrid), intent(in) :: this
        real(rkind), dimension(-num_pad+1:,-num_pad+1:,-num_pad+1:), intent(in) :: f1, f2
        character(len=*), intent(in) :: label
        real(rkind), intent(out), optional :: errf
        real(rkind) :: err

        err = p_maxval(maxval(abs(f1(1:this%zsz(1),1:this%zsz(2),1:this%zsz(3)) - f2(1:this%zsz(1),1:this%zsz(2),1:this%zsz(3)))))
        if (nrank == 0) then
            print*, "====================="
            print*, "Error in field: ", label
            print*, ">> ", err
            print*, "====================="
        end if 

        if (present(errf)) errf = err

    end subroutine
    
    function get_min_interior(this, arr) result(min_interior)
        use reductions, only: p_minval
        class(procgrid), intent(in) :: this
        real(rkind), dimension(-num_pad+1:,-num_pad+1:,-num_pad+1:), intent(in) :: arr
        real(rkind) :: min_interior

        min_interior = p_minval(minval(arr(1:this%zsz(1),1:this%zsz(2),1:this%zsz(3))))

    end function 
    
    function get_max_interior(this, arr) result(max_interior)
        use reductions, only: p_maxval
        class(procgrid), intent(in) :: this
        real(rkind), dimension(-num_pad+1:,-num_pad+1:,-num_pad+1:), intent(in) :: arr
        real(rkind) :: max_interior

        max_interior = p_maxval(maxval(arr(1:this%zsz(1),1:this%zsz(2),1:this%zsz(3))))

    end function 

    subroutine alloc_array3(this, arr)
        class(procgrid), intent(in) :: this
        real(rkind), dimension(:,:,:), allocatable, intent(out) :: arr

        if (allocated(arr)) deallocate(arr)
        allocate(arr(-num_pad+1:this%decomp%zsz(1)+num_pad,-num_pad+1:this%decomp%zsz(2)+num_pad,-num_pad+1:this%decomp%zsz(3)+num_pad))
    end subroutine 

    subroutine alloc_array4(this, arr,nfields)
        class(procgrid), intent(in) :: this
        integer, intent(in) :: nfields
        real(rkind), dimension(:,:,:,:), allocatable, intent(out) :: arr

        if (allocated(arr)) deallocate(arr)
        allocate(arr(-num_pad+1:this%decomp%zsz(1)+num_pad,-num_pad+1:this%decomp%zsz(2)+num_pad,-num_pad+1:this%decomp%zsz(3)+num_pad, nfields))
    end subroutine 

    subroutine alloc_array5(this, arr,n1,n2)
        class(procgrid), intent(in) :: this
        integer, intent(in) :: n1,n2
        real(rkind), dimension(:,:,:,:,:), allocatable, intent(out) :: arr

        if (allocated(arr)) deallocate(arr)
        allocate(arr(-num_pad+1:this%decomp%zsz(1)+num_pad,-num_pad+1:this%decomp%zsz(2)+num_pad,-num_pad+1:this%decomp%zsz(3)+num_pad, n1, n2))
    end subroutine 
    
    subroutine init(this, prow, pcol, nx, ny, nz)
        use exits, only: message
        class(procgrid), intent(inout) :: this
        logical, dimension(3) :: periodicbcs
        integer, intent(in) :: prow, pcol, nx, ny, nz
        integer :: nrow, ncol, ntasks, ierr 
    
        integer :: DECOMP_CART_Z
    
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
            print*, "Processor decomposition failed."
            call mpi_abort(MPI_COMM_WORLD, 0, ierr) 
        end if 

        if (.not. check_stencil_validity(this%nx,this%ny,nrow,ncol)) then
            !call message(0,"Default proc decomposition incompatible with stencil.")
            !call message(0,"Will try to swap order.")
            call swap_proc_order(nrow,ncol)
            if (.not. check_stencil_validity(this%nx,this%ny,nrow,ncol)) then
                print*, nrow, ncol
                print*, "Need a better processor decomposition algorithm for this domain. Current algorithm fails"
                call mpi_abort(MPI_COMM_WORLD, 0, ierr) 
            end if 
        end if         
    

        ! Set all BCs to true (will handle non-periodicity separately)
        periodicbcs(1) = .true.; periodicbcs(2) = .true.; periodicbcs(3) = .true.
        call decomp_2d_init(this%nx, this%ny, this%nz, nrow, ncol, periodicbcs)
        call get_decomp_info(this%decomp)
        this%zsz = this%decomp%zsz
        this%zst = this%decomp%zst
        this%zen = this%decomp%zen

        ! Get neighbors
        call MPI_CART_CREATE(MPI_COMM_WORLD,2,[nrow,ncol],[.true.,.true.],.false.,DECOMP_CART_Z,ierr)
        call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 0, 1, this%XneighLeft, this%XneighRight, ierr) ! east & west
        call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 1, 1, this%YneighDown, this%YneighUp, ierr) ! east & west

        this%nrow = nrow
        this%ncol = ncol 
        allocate(this%halo_sendbuff_xleft (num_pad,this%decomp%zsz(2),this%decomp%zsz(3)))
        allocate(this%halo_sendbuff_xright(num_pad,this%decomp%zsz(2),this%decomp%zsz(3)))
        allocate(this%halo_sendbuff_yup  (this%decomp%zsz(1)+2*num_pad,num_pad,this%decomp%zsz(3)))
        allocate(this%halo_sendbuff_ydown(this%decomp%zsz(1)+2*num_pad,num_pad,this%decomp%zsz(3)))
        
        allocate(this%halo_recvbuff_xleft (num_pad,this%decomp%zsz(2),this%decomp%zsz(3)))
        allocate(this%halo_recvbuff_xright(num_pad,this%decomp%zsz(2),this%decomp%zsz(3)))
        allocate(this%halo_recvbuff_yup  (this%decomp%zsz(1)+2*num_pad,num_pad,this%decomp%zsz(3)))
        allocate(this%halo_recvbuff_ydown(this%decomp%zsz(1)+2*num_pad,num_pad,this%decomp%zsz(3)))
        call message(0,"Processor decomposition information:")
        call message(1,"Nrows:", nrow)
        call message(1,"Ncols:", ncol)

        this%have_x_loEdge = .false.  
        this%have_y_loEdge = .false.  
        this%have_z_loEdge = .false.  
        
        this%have_x_hiEdge = .false.  
        this%have_y_hiEdge = .false.  
        this%have_z_hiEdge = .false.  

        if (this%zst(1)  .eq.       1) this%have_x_loEdge = .true.
        if (this%zen(1) .eq. this%nx) this%have_x_hiEdge = .true.
        
        if (this%zst(2)  .eq.       1) this%have_y_loEdge = .true.
        if (this%zen(2) .eq. this%ny) this%have_y_hiEdge = .true.
            
        if (this%zst(3)  .eq.       1) this%have_z_loEdge = .true.
        if (this%zen(3) .eq. this%nz) this%have_z_hiEdge = .true.
        

    end subroutine 

    subroutine destroy(this)
        class(procgrid), intent(inout) :: this

        deallocate(this%halo_sendbuff_xleft)
        deallocate(this%halo_sendbuff_xright)
        deallocate(this%halo_recvbuff_xleft)
        deallocate(this%halo_recvbuff_xright)
        
        deallocate(this%halo_sendbuff_yup)
        deallocate(this%halo_sendbuff_ydown)
        deallocate(this%halo_recvbuff_yup)
        deallocate(this%halo_recvbuff_ydown)

    end subroutine 

    subroutine halo_exchange(this, ufield)
        use kind_parameters, only: mpirkind
        class(procgrid), intent(inout) :: this
        real(rkind), dimension(-num_pad+1:,-num_pad+1:,-num_pad+1:), intent(inout) :: ufield
        integer :: send_req1, recv_req1, status(MPI_STATUS_SIZE)
        integer :: send_req2, recv_req2, ierr, tag1 = 123, tag2 = 321
        integer :: nx, ny, nz, idx
        
        nx = this%decomp%zsz(1)
        ny = this%decomp%zsz(2)
        nz = this%decomp%zsz(3)

        ! X - halo exchange
        
        call MPI_IRECV(this%halo_recvbuff_xleft ,ny*nz*num_pad,mpirkind,this%XneighLeft ,tag1,MPI_COMM_WORLD,recv_req1,ierr)
        call MPI_IRECV(this%halo_recvbuff_xright,ny*nz*num_pad,mpirkind,this%XneighRight,tag2,MPI_COMM_WORLD,recv_req2,ierr)
        
        this%halo_sendbuff_xleft  = ufield(1:num_pad      ,1:ny,1:nz)
        this%halo_sendbuff_xright = ufield(nx-num_pad+1:nx,1:ny,1:nz)
  
        call MPI_ISEND(this%halo_sendbuff_xleft ,ny*nz*num_pad,mpirkind,this%XneighLeft ,tag2,MPI_COMM_WORLD,send_req1,ierr)
        call MPI_ISEND(this%halo_sendbuff_xright,ny*nz*num_pad,mpirkind,this%XneighRight,tag1,MPI_COMM_WORLD,send_req2,ierr)

        call MPI_WAIT(recv_req1,status,ierr)
        call MPI_WAIT(recv_req2,status,ierr)

        ufield(-num_pad+1:0   ,1:ny,1:nz) = this%halo_recvbuff_xleft
        ufield(nx+1:nx+num_pad,1:ny,1:nz) = this%halo_recvbuff_xright
        
        call MPI_WAIT(send_req1,status,ierr)
        call MPI_WAIT(send_req2,status,ierr)

        ! Y - halo exchange 

        call MPI_IRECV(this%halo_recvbuff_ydown,(nx+2*num_pad)*nz*num_pad,mpirkind,this%YneighDown,tag1,MPI_COMM_WORLD,recv_req1,ierr)
        call MPI_IRECV(this%halo_recvbuff_yup  ,(nx+2*num_pad)*nz*num_pad,mpirkind,this%YneighUp  ,tag2,MPI_COMM_WORLD,recv_req2,ierr)
   
        this%halo_sendbuff_ydown  = ufield(-num_pad+1:nx+num_pad,1:num_pad      ,1:nz)
        this%halo_sendbuff_yup    = ufield(-num_pad+1:nx+num_pad,ny-num_pad+1:ny,1:nz)
        
        call MPI_ISEND(this%halo_sendbuff_ydown,(nx+2*num_pad)*nz*num_pad,mpirkind,this%YneighDown,tag2,MPI_COMM_WORLD,send_req1,ierr)
        call MPI_ISEND(this%halo_sendbuff_yup  ,(nx+2*num_pad)*nz*num_pad,mpirkind,this%YneighUp  ,tag1,MPI_COMM_WORLD,send_req2,ierr)
        
        call MPI_WAIT(recv_req1,status,ierr)
        call MPI_WAIT(recv_req2,status,ierr)

        ufield(-num_pad+1:nx+num_pad,-num_pad+1:0,1:nz) = this%halo_recvbuff_ydown
        ufield(-num_pad+1:nx+num_pad,ny+1:ny+num_pad,1:nz) = this%halo_recvbuff_yup 
        
        call MPI_WAIT(send_req1,status,ierr)
        call MPI_WAIT(send_req2,status,ierr)
    
        ! Z - exchange (swap)
        ufield(-num_pad+1:nx+num_pad,-num_pad+1:ny+num_pad,-num_pad+1:0)    = ufield(-num_pad+1:nx+num_pad,-num_pad+1:ny+num_pad,nz-num_pad+1:nz)
        ufield(-num_pad+1:nx+num_pad,-num_pad+1:ny+num_pad,nz+1:nz+num_pad) = ufield(-num_pad+1:nx+num_pad,-num_pad+1:ny+num_pad,1:num_pad)
   
        if (PAD_SPURIOUS) then
            ufield(-num_pad+1,:,:) = 9.d99
            ufield(nx+num_pad,:,:) = 9.d99
            ufield(:,-num_pad+1,:) = 9.d99
            ufield(:,ny+num_pad,:) = 9.d99
            ufield(:,:,-num_pad+1) = 9.d99
            ufield(:,:,nz+num_pad) = 9.d99
        end if 
        
        !if (present(bcs) .and. present(bcval)) then
        !    if (this%have_x_loEdge) then
        !        select case (bcs(1))
        !        case(1)
        !            do idx = -num_pad+1,0
        !                ufield(idx,:,:) = ufield(1,:,:) 
        !            end do 
        !        case(2)
        !            ufield(-num_pad+1:0,:,:) = bcval(1)
        !        case default
        !            ! do nothing 
        !        end select 
        !    end if 
        !    
        !    if (this%have_x_hiEdge) then
        !        select case (bcs(2)) 
        !        case(1)
        !            do idx = this%zsz(1)+1,this%zsz(1)+num_pad
        !                ufield(idx,:,:) = ufield(this%zsz(1),:,:) 
        !            end do 
        !        case(2)
        !            ufield(this%zsz(1)+1:this%zsz(1)+num_pad,:,:) = bcval(2)
        !        case default
        !            ! do nothing 
        !        end select
        !    end if 


        !    if (this%have_y_loEdge) then
        !        select case (bcs(3))
        !        case(1)
        !            do idx = -num_pad+1,0
        !                ufield(:,idx,:) = ufield(:,1,:) 
        !            end do 
        !        case(2)
        !            ufield(:,-num_pad+1:0,:) = bcval(3)
        !        case default
        !            ! do nothing 
        !        end select 
        !    end if 
        !    
        !    if (this%have_y_hiEdge) then
        !        select case (bcs(4)) 
        !        case(1)
        !            do idx = this%zsz(2)+1,this%zsz(2)+num_pad
        !                ufield(:,idx,:) = ufield(:,this%zsz(2),:) 
        !            end do 
        !        case(2)
        !            ufield(:,this%zsz(2)+1:this%zsz(2)+num_pad,:) = bcval(4)
        !        case default
        !            ! do nothing 
        !        end select
        !    end if 
        !    
        !    if (this%have_z_loEdge) then
        !        select case (bcs(5))
        !        case(1)
        !            do idx = -num_pad+1,0
        !                ufield(:,:,idx) = ufield(:,:,1) 
        !            end do 
        !        case(2)
        !            ufield(:,:,-num_pad+1:0) = bcval(5)
        !        case default
        !            ! do nothing 
        !        end select 
        !    end if 
        !    
        !    if (this%have_z_hiEdge) then
        !        select case (bcs(6)) 
        !        case(1)
        !            do idx = this%zsz(3)+1,this%zsz(3)+num_pad
        !                ufield(:,:,idx) = ufield(:,:,this%zsz(3)) 
        !            end do 
        !        case(2)
        !            ufield(:,:,this%zsz(3)+1:this%zsz(3)+num_pad) = bcval(6)
        !        case default
        !            ! do nothing 
        !        end select
        !    end if 

        !end if 
    
    
    end subroutine 




    subroutine halo_exchange_x(this, ufield)
        use kind_parameters, only: mpirkind
        class(procgrid), intent(inout) :: this
        real(rkind), dimension(-num_pad+1:,-num_pad+1:,-num_pad+1:), intent(inout) :: ufield
        integer :: send_req1, recv_req1, status(MPI_STATUS_SIZE)
        integer :: send_req2, recv_req2, ierr, tag1 = 123, tag2 = 321
        integer :: nx, ny, nz
        
        nx = this%decomp%zsz(1)
        ny = this%decomp%zsz(2)
        nz = this%decomp%zsz(3)

        ! X - halo exchange
        call MPI_IRECV(this%halo_recvbuff_xleft ,ny*nz*num_pad,mpirkind,this%XneighLeft ,tag1,MPI_COMM_WORLD,recv_req1,ierr)
        call MPI_IRECV(this%halo_recvbuff_xright,ny*nz*num_pad,mpirkind,this%XneighRight,tag2,MPI_COMM_WORLD,recv_req2,ierr)
        
        this%halo_sendbuff_xleft  = ufield(1:num_pad      ,1:ny,1:nz)
        this%halo_sendbuff_xright = ufield(nx-num_pad+1:nx,1:ny,1:nz)
  
        call MPI_ISEND(this%halo_sendbuff_xleft ,ny*nz*num_pad,mpirkind,this%XneighLeft ,tag2,MPI_COMM_WORLD,send_req1,ierr)
        call MPI_ISEND(this%halo_sendbuff_xright,ny*nz*num_pad,mpirkind,this%XneighRight,tag1,MPI_COMM_WORLD,send_req2,ierr)

        call MPI_WAIT(recv_req1,status,ierr)
        call MPI_WAIT(recv_req2,status,ierr)

        ufield(-num_pad+1:0   ,1:ny,1:nz) = this%halo_recvbuff_xleft
        ufield(nx+1:nx+num_pad,1:ny,1:nz) = this%halo_recvbuff_xright
        
        call MPI_WAIT(send_req1,status,ierr)
        call MPI_WAIT(send_req2,status,ierr)

    end subroutine 

    subroutine halo_exchange_y(this, ufield)
        use kind_parameters, only: mpirkind
        class(procgrid), intent(inout) :: this
        real(rkind), dimension(-num_pad+1:,-num_pad+1:,-num_pad+1:), intent(inout) :: ufield
        integer :: send_req1, recv_req1, status(MPI_STATUS_SIZE)
        integer :: send_req2, recv_req2, ierr, tag1 = 123, tag2 = 321
        integer :: nx, ny, nz
        
        nx = this%decomp%zsz(1)
        ny = this%decomp%zsz(2)
        nz = this%decomp%zsz(3)

        call MPI_IRECV(this%halo_recvbuff_ydown,(nx+2*num_pad)*nz*num_pad,mpirkind,this%YneighDown,tag1,MPI_COMM_WORLD,recv_req1,ierr)
        call MPI_IRECV(this%halo_recvbuff_yup  ,(nx+2*num_pad)*nz*num_pad,mpirkind,this%YneighUp  ,tag2,MPI_COMM_WORLD,recv_req2,ierr)
   
        this%halo_sendbuff_ydown  = ufield(-num_pad+1:nx+num_pad,1:num_pad      ,1:nz)
        this%halo_sendbuff_yup    = ufield(-num_pad+1:nx+num_pad,ny-num_pad+1:ny,1:nz)
        
        call MPI_ISEND(this%halo_sendbuff_ydown,(nx+2*num_pad)*nz*num_pad,mpirkind,this%YneighDown,tag2,MPI_COMM_WORLD,send_req1,ierr)
        call MPI_ISEND(this%halo_sendbuff_yup  ,(nx+2*num_pad)*nz*num_pad,mpirkind,this%YneighUp  ,tag1,MPI_COMM_WORLD,send_req2,ierr)
        
        call MPI_WAIT(recv_req1,status,ierr)
        call MPI_WAIT(recv_req2,status,ierr)

        ufield(-num_pad+1:nx+num_pad,-num_pad+1:0,1:nz) = this%halo_recvbuff_ydown
        ufield(-num_pad+1:nx+num_pad,ny+1:ny+num_pad,1:nz) = this%halo_recvbuff_yup 
        
        call MPI_WAIT(send_req1,status,ierr)
        call MPI_WAIT(send_req2,status,ierr)

    end subroutine 
    
    subroutine halo_exchange_z(this, ufield)
        class(procgrid), intent(inout) :: this
        real(rkind), dimension(-num_pad+1:,-num_pad+1:,-num_pad+1:), intent(inout) :: ufield
        integer :: nx, ny, nz

        nx = this%decomp%zsz(1)
        ny = this%decomp%zsz(2)
        nz = this%decomp%zsz(3)
        
        ufield(-num_pad+1:nx+num_pad,-num_pad+1:ny+num_pad,-num_pad+1:0)    = ufield(-num_pad+1:nx+num_pad,-num_pad+1:ny+num_pad,nz-num_pad+1:nz)
        ufield(-num_pad+1:nx+num_pad,-num_pad+1:ny+num_pad,nz+1:nz+num_pad) = ufield(-num_pad+1:nx+num_pad,-num_pad+1:ny+num_pad,1:num_pad)

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
