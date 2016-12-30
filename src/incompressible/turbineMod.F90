module turbineMod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa
    use decomp_2d
    use StaggOpsMod, only: staggOps  
    use actuatorDiskMod, only: actuatorDisk
    use actuatorLineMod, only: actuatorLine
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral  
    use mpi 
    use reductions, only: p_maxval, p_sum
    use timer, only: tic, toc

    implicit none

    private
    public :: TurbineArray

    ! default initializations
    integer :: num_turbines = 1
    logical :: ADM = .TRUE. ! .FALSE. implies ALM
    character(len=clen) :: turbInfoDir

    integer :: ioUnit

    real(rkind), parameter :: degrees_to_radians = pi/180.0_rkind
    complex(rkind), parameter :: zeroC = zero + imi*zero

    type :: TurbineArray
        integer :: nTurbines
        integer :: myProc
        type(actuatorDisk), allocatable, dimension(:) :: turbArrayADM
        type(actuatorLine), allocatable, dimension(:) :: turbArrayALM

        type(decomp_info), pointer :: gpC, sp_gpC, gpE, sp_gpE
        type(spectral), pointer :: spectC, spectE
        type(staggOps), allocatable :: OpsNU
 
        real(rkind), dimension(:,:,:), allocatable :: local_rbuffxC
        real(rkind), dimension(:,:,:), pointer :: fx, fy, fz
        complex(rkind), dimension(:,:,:), pointer :: fChat, fEhat, zbuffC, zbuffE

        ! variables needed for halo communication
        integer :: neighbour(6), coord(2), dims(2), tag_s, tag_n, tag_b, tag_t
        real(rkind), allocatable, dimension(:,:,:) :: ySendBuf, zSendBuf, yRightHalo, zRightHalo, zLeftHalo

        logical :: dumpTurbField = .false.
        integer :: step = 0
    contains

        procedure :: init
        procedure :: destroy
        procedure :: init_halo_communication
        procedure :: halo_communication
        procedure :: destroy_halo_communication
        procedure :: getForceRHS 
        procedure :: dumpFullField

    end type

contains

subroutine init(this, inputFile, gpC, gpE, spectC, spectE, rbuffxC, cbuffyC, cbuffYE, cbuffzC, cbuffzE, mesh, dx, dy, dz)
    class(TurbineArray), intent(inout), target :: this
    character(len=*),    intent(in)            :: inputFile
    type(spectral), target :: spectC, spectE
    type(decomp_info), target :: gpC, gpE!, sp_gpC, sp_gpE
    real(rkind), dimension(:,:,:,:), target :: rbuffxC   ! actually 3 are required
    complex(rkind), dimension(:,:,:,:), target :: cbuffyC, cbuffyE, cbuffzC, cbuffzE
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), intent(in) :: dx, dy, dz
    logical :: useWindTurbines = .TRUE. ! .FALSE. implies ALM
    real(rkind) :: xyzPads(6)

    integer :: i, ierr

    namelist /WINDTURBINES/ useWindTurbines, num_turbines, ADM, turbInfoDir

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=WINDTURBINES)
    close(ioUnit)

    this%gpC => gpC
    this%spectC => spectC
    this%sp_gpC => this%spectC%spectdecomp

    this%gpE => gpE
    this%spectE => spectE
    this%sp_gpE => this%spectE%spectdecomp

    allocate(this%local_rbuffxC(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))

    this%fx => rbuffxC(:,:,:,1); this%fy => rbuffxC(:,:,:,2);  this%fz => this%local_rbuffxC(:,:,:)
    this%fChat => cbuffyC(:,:,:,1); this%fEhat => cbuffyE(:,:,:,1)
    this%zbuffC => cbuffzC(:,:,:,1); this%zbuffE => cbuffzE(:,:,:,1)

    allocate(this%OpsNU)
    call this%OpsNU%init(this%gpC,this%gpE,0,dx,dy,dz,this%spectC%spectdecomp,this%spectE%spectdecomp,.true.,.true.)

    ! set number of turbines
    this%nTurbines = num_turbines;

    if(ADM) then
      allocate (this%turbArrayADM(this%nTurbines))
      do i = 1, this%nTurbines
        call this%turbArrayADM(i)%init(turbInfoDir, i, mesh(:,:,:,1), mesh(:,:,:,2), mesh(:,:,:,3))
      end do
      call message(1,"WIND TURBINE ADM model initialized")
    else
       call mpi_barrier(mpi_comm_world, ierr); call message(1,"Initializing WIND TURBINE ALM model")
      call mpi_barrier(mpi_comm_world, ierr); call message(1,"Setting up for halo communication")
      call this%init_halo_communication(mesh(:,:,:,1), mesh(:,:,:,2), mesh(:,:,:,3), dx, dy, dz, xyzPads)
      call mpi_barrier(mpi_comm_world, ierr); call message(1,"Done setting up for halo communication")
      allocate (this%turbArrayALM(this%nTurbines))
      do i = 1, this%nTurbines
        call this%turbArrayALM(i)%init(turbInfoDir, i, mesh(:,:,:,1), mesh(:,:,:,2), mesh(:,:,:,3), xyzPads)
      end do
      call message(1,"WIND TURBINE ALM model initialized")
    endif

    if(this%dumpTurbField) then
        this%fx = zero; this%fy = zero; this%fz = zero
        call this%dumpFullField(this%fx, "wtfx")
        call this%dumpFullField(this%fy, "wtfy")
        call this%dumpFullField(this%fz, "wtfz")
        this%dumpTurbField = .false.
    endif

end subroutine


subroutine destroy(this)
    class(TurbineArray), intent(inout) :: this
    integer :: i, ierr

    nullify(this%gpC, this%gpE, this%spectC, this%sp_gpC, this%fx, this%fy, this%fz)
    nullify(this%zbuffC, this%zbuffE, this%fChat, this%fEhat)
    deallocate(this%local_rbuffxC)

    deallocate(this%OpsNU)

    if(ADM) then
      do i = 1, this%nTurbines
        call this%turbArrayADM(i)%destroy()
      end do
      deallocate(this%turbArrayADM)
    else
      call this%destroy_halo_communication()
      do i = 1, this%nTurbines
        call this%turbArrayALM(i)%destroy()
      end do
      deallocate(this%turbArrayALM)
    endif

end subroutine

subroutine destroy_halo_communication(this)
    class(TurbineArray), intent(inout) :: this

    deallocate(this%yRightHalo, this%zRightHalo, this%zLeftHalo, this%ySendBuf, this%zSendBuf)
end subroutine

subroutine init_halo_communication(this, xG, yG, zG, dx, dy, dz ,xyzPads)
  use decomp_2d
  class(TurbineArray),           intent(inout) :: this
  real(rkind), dimension(:,:,:), intent(in)    :: xG, yG, zG
  real(rkind),                   intent(in)    :: dx, dy, dz
  real(rkind), dimension(:),     intent(out)   :: xyzPads

  integer :: ierror
  logical :: periodic(2)

  allocate( this%yRightHalo(this%gpC%xsz(1), this%gpC%xsz(3),   3) )
  allocate( this%zRightHalo(this%gpC%xsz(1), this%gpC%xsz(2)+1, 3) )
  allocate( this%zLeftHalo (this%gpC%xsz(1), this%gpC%xsz(2)+1, 3) )
  allocate( this%ySendBuf  (this%gpC%xsz(1), this%gpC%xsz(3),   3) )
  allocate( this%zSendBuf  (this%gpC%xsz(1), this%gpC%xsz(2)+1, 3) )


  ! Initialize neighbour information - copied from halo.f90 in 2decomp&fft
  ! For X-pencil 
  this%neighbour(1) = MPI_PROC_NULL               ! east
  this%neighbour(2) = MPI_PROC_NULL               ! west
  call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 0, 1, this%neighbour(4), this%neighbour(3), ierror) ! north & south
  call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 1, 1, this%neighbour(6), this%neighbour(5), ierror) ! top & bottom

  ! set coords
  call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, this%dims, periodic, this%coord, ierror)

  ! now set tags from coords -- copied from halo_common.f90
     ! *** north/south *** 
  this%tag_s = this%coord(1)
  if (this%coord(1)==this%dims(1)-1 .AND. periodic(1)) then
     this%tag_n = 0
  else
     this%tag_n = this%coord(1) + 1
  end if

     ! *** top/bottom *** 
  this%tag_b = this%coord(2)
  if (this%coord(2)==this%dims(2)-1 .AND. periodic(2)) then
     this%tag_t = 0
  else
     this%tag_t = this%coord(2) + 1
  end if

  !write(*,'(5(i5,1x),2(L4,1x))') nrank, this%coord(1:2), this%dims(1:2), periodic(1:2)
  !write(*,'(13(i5,1x))') nrank, this%neighbour(1:6), this%tag_s, this%tag_n, this%coord(1:2), this%dims(1:2)

  ! populate xyzPads
  ! we are in x decomp
  xyzPads(1) = xG(1,1,1)                  !- this is not needed so far
  xyzPads(2) = xG(this%gpC%xsz(1),1,1) + dx 
  
  if(this%coord(1)==0) then
     xyzPads(3) = yG(1,1,1)
     if(periodic(1)) xyzPads(3) = yG(1,1,1) + ny_global*dy
  else
     xyzPads(3) = yG(1,1,1) - dy
  endif

  if(this%coord(1)==this%dims(1)-1) then
     xyzPads(4) = yG(1,this%gpC%xsz(2),1)
     if(periodic(1)) xyzPads(4) = yG(1,this%gpc%xsz(2),1) + dy
  else
     xyzPads(4) = yG(1,this%gpC%xsz(2),1) + dy
  endif

  !if(this%coord(2)==0) then
     xyzPads(5) = zG(1,1,1) - 0.5D0*dz
  !else
  !   xyzPads(5) = zG(1,1,1) - dz
  !endif
  
  !if(this%coord(2)==this%dims(2)-1) then
     xyzPads(6) = zG(1,1,this%gpC%xsz(3)) + 0.5D0*dz
  !else
  !   xyzPads(6) = zG(1,1,this%gpC%xsz(3)) + dz
  !endif

end subroutine

subroutine halo_communication(this, u, v, wC)
  use decomp_2d
  use kind_parameters, only: mpirkind
  class(TurbineArray), intent(inout) :: this
  real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)),intent(in) :: u, v, wC

  integer :: k, icount, requests(2), ierror, ierr
  integer, dimension(MPI_STATUS_SIZE,2) :: status

  ! fill this%yRightHalo
  ! pack data
  do k = 1, this%gpC%xsz(3)
    this%ySendBuf(:,k,1) = u(:, 1, k)
  enddo
  do k = 1, this%gpC%xsz(3)
    this%ySendBuf(:,k,2) = v(:, 1, k)
  enddo
  do k = 1, this%gpC%xsz(3)
    this%ySendBuf(:,k,3) = wC(:, 1, k)
  enddo
  ! receive YRight from right
  icount = this%gpC%xsz(1) * this%gpC%xsz(3) * 3
  call MPI_IRECV(this%yRightHalo, icount, mpirkind, this%neighbour(3), this%tag_n, DECOMP_2D_COMM_CART_X, requests(1), ierror)
  ! send to YRight to left
  call MPI_ISSEND(this%ySendBuf, icount, mpirkind, this%neighbour(4), this%tag_s, DECOMP_2D_COMM_CART_X, requests(2), ierror)
  call MPI_WAITALL(2, requests, status, ierror)


  ! fill this%zLeftHalo
  ! pack data
  this%zSendBuf(:,1:this%gpC%xsz(2),1) = u(:,:,this%gpC%xsz(3)); this%zSendBuf(:,this%gpC%xsz(2)+1,1) = this%yRightHalo(:,this%gpC%xsz(3),1)
  this%zSendBuf(:,1:this%gpC%xsz(2),2) = v(:,:,this%gpC%xsz(3)); this%zSendBuf(:,this%gpC%xsz(2)+1,2) = this%yRightHalo(:,this%gpC%xsz(3),2)
  this%zSendBuf(:,1:this%gpC%xsz(2),3) = wC(:,:,this%gpC%xsz(3)); this%zSendBuf(:,this%gpC%xsz(2)+1,3) = this%yRightHalo(:,this%gpC%xsz(3),3)
  ! receive ZLeft from Left
  icount = this%gpC%xsz(1) * (this%gpC%xsz(2) + 1) * 3
  call MPI_IRECV(this%zLeftHalo, icount, mpirkind, this%neighbour(6), this%tag_b, DECOMP_2D_COMM_CART_X, requests(1), ierror)
  ! send ZLeft to Right
  call MPI_ISSEND(this%zSendBuf,  icount, mpirkind, this%neighbour(5), this%tag_t, DECOMP_2D_COMM_CART_X, requests(2), ierror)
  call MPI_WAITALL(2, requests, status, ierror)


  ! fill this%zRightHalo
  ! pack data
  this%zSendBuf(:,1:this%gpC%xsz(2),1) = u(:,:,1); this%zSendBuf(:,this%gpC%xsz(2)+1,1) = this%yRightHalo(:,1,1)
  this%zSendBuf(:,1:this%gpC%xsz(2),2) = v(:,:,1); this%zSendBuf(:,this%gpC%xsz(2)+1,2) = this%yRightHalo(:,1,2)
  this%zSendBuf(:,1:this%gpC%xsz(2),3) = wC(:,:,1); this%zSendBuf(:,this%gpC%xsz(2)+1,3) = this%yRightHalo(:,1,3)
  ! receive ZRight from Right
  icount = this%gpC%xsz(1) * (this%gpC%xsz(2) + 1) * 3
  call MPI_IRECV(this%zRightHalo, icount, mpirkind, this%neighbour(5), this%tag_t, DECOMP_2D_COMM_CART_X, requests(1), ierror)
  ! send ZRight to Left
  call MPI_ISSEND(this%zSendBuf, icount, mpirkind, this%neighbour(6), this%tag_b, DECOMP_2D_COMM_CART_X, requests(2), ierror)
  call MPI_WAITALL(2, requests, status, ierror)

  if(this%coord(2)==0) then
      this%zLeftHalo(:,1:this%gpC%xsz(2),1) = two*u(:,:,1) - u(:,:,2);   this%zLeftHalo(:,this%gpC%xsz(2)+1,1) = two*this%yRightHalo(:,1,1) - this%yRightHalo(:,2,1)
      this%zLeftHalo(:,1:this%gpC%xsz(2),2) = two*v(:,:,1) - v(:,:,2);   this%zLeftHalo(:,this%gpC%xsz(2)+1,2) = two*this%yRightHalo(:,1,2) - this%yRightHalo(:,2,2)
      this%zLeftHalo(:,1:this%gpC%xsz(2),3) = two*wC(:,:,1) - wC(:,:,2); this%zLeftHalo(:,this%gpC%xsz(2)+1,3) = two*this%yRightHalo(:,1,3) - this%yRightHalo(:,2,3)
  elseif(this%coord(2)==this%dims(2)-1) then
      this%zRightHalo(:,1:this%gpC%xsz(2),1) = two*u(:,:,this%gpC%xsz(3)) - u(:,:,this%gpC%xsz(3)-1);   this%zRightHalo(:,this%gpC%xsz(2)+1,1) = two*this%yRightHalo(:,this%gpC%xsz(3),1) - this%yRightHalo(:,this%gpC%xsz(3)-1,1)
      this%zRightHalo(:,1:this%gpC%xsz(2),2) = two*v(:,:,this%gpC%xsz(3)) - v(:,:,this%gpC%xsz(3)-1);   this%zRightHalo(:,this%gpC%xsz(2)+1,2) = two*this%yRightHalo(:,this%gpC%xsz(3),2) - this%yRightHalo(:,this%gpC%xsz(3)-1,2)
      this%zRightHalo(:,1:this%gpC%xsz(2),3) = two*wC(:,:,this%gpC%xsz(3)) - wC(:,:,this%gpC%xsz(3)-1); this%zRightHalo(:,this%gpC%xsz(2)+1,3) = two*this%yRightHalo(:,this%gpC%xsz(3),3) - this%yRightHalo(:,this%gpC%xsz(3)-1,3)
  endif

end subroutine

subroutine getForceRHS(this, dt, u, v, wC, urhs, vrhs, wrhs, inst_horz_avg)
    class(TurbineArray), intent(inout), target :: this
    real(rkind),                                                                         intent(in) :: dt
    real(rkind),    dimension(this%gpC%xsz(1),   this%gpC%xsz(2),   this%gpC%xsz(3)),    intent(in) :: u, v, wC
    complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout) :: urhs, vrhs
    complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: wrhs 
    real(rkind),    dimension(:),                                                        intent(out), optional   :: inst_horz_avg
    integer :: i, ierr

    this%fx = zero; this%fy = zero; this%fz = zero
    if(ADM) then
      do i = 1, this%nTurbines
        ! CHANGED to allow avoiding inst_horz_avg calculations - useful for
        ! testing/debugging
        if (present(inst_horz_avg)) then
            call this%turbArrayADM(i)%get_RHS(u,v,wC,this%fx,this%fy,this%fz,inst_horz_avg(8*i-7:8*i))
        else
            call this%turbArrayADM(i)%get_RHS(u,v,wC,this%fx,this%fy,this%fz)   
        end if 
      end do
    else
       !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Starting halo communication"); call mpi_barrier(mpi_comm_world, ierr)
      call this%halo_communication(u, v, wC)
       !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Done halo communication"); call mpi_barrier(mpi_comm_world, ierr)
      do i = 1, this%nTurbines
        ! CHANGED to allow avoiding inst_horz_avg calculations - useful for
        ! testing/debugging
        if (present(inst_horz_avg)) then
            call this%turbArrayALM(i)%get_RHS(dt, u, v, wC, this%yRightHalo, this%zLeftHalo, this%zRightHalo, this%fx, this%fy, this%fz, inst_horz_avg(8*i-7:8*i))
        else
            call this%turbArrayALM(i)%get_RHS(dt, u, v, wC, this%yRightHalo, this%zLeftHalo, this%zRightHalo, this%fx, this%fy, this%fz)   
        end if 
      end do
    endif
       !call mpi_barrier(mpi_comm_world, ierr)
       !write(*,'(i5,6(e19.12,1x))') nrank, maxval(this%fx), minval(this%fx), maxval(this%fy), minval(this%fy), maxval(this%fz), minval(this%fz) 
       !call mpi_barrier(mpi_comm_world, ierr)

    ! add forces to rhs
    call this%spectC%fft(this%fx,this%fChat)
    urhs = urhs + this%fChat

    call this%spectC%fft(this%fy,this%fChat)
    vrhs = vrhs + this%fChat

    call this%spectC%fft(this%fz,this%fChat)
    ! interpolate fz to fzE
    call transpose_y_to_z(this%fChat,this%zbuffC,this%sp_gpC)
    call this%OpsNU%InterpZ_Cell2Edge(this%zbuffC,this%zbuffE,zeroC,zeroC)
    call transpose_z_to_y(this%zbuffE,this%fEhat,this%sp_gpE)
    wrhs = wrhs + this%fEhat

    if(this%dumpTurbField) then
      !if(ADM) then
      !else
      !  do i = 1, this%nTurbines
      !    call this%turbArrayALM(i)%dumpTurbineField()
      !  enddo
      !endif
      this%step = this%step + 1
      call this%dumpFullField(this%fx, 'wtfx')
      call this%dumpFullField(this%fy, 'wtfy')
      call this%dumpFullField(this%fz, 'wtfz')
      this%dumpTurbField = .false.
    endif

end subroutine 

    subroutine dumpFullField(this,arr,label)
        use decomp_2d_io
        use mpi
        use exits, only: message
        class(TurbineArray), intent(in) :: this
        character(len=clen) :: tempname!, fname
        real(rkind), dimension(:,:,:), intent(in) :: arr
        character(len=4), intent(in) :: label

        write(tempname,"(A6,A4,A2,I6.6,A4)") "Run02_",label,"_t",this%step,".out"
        !fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,arr,tempname)

    end subroutine


end module
