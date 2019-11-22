module turbineMod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa
    use decomp_2d
    use StaggOpsMod, only: staggOps  
    use actuatorDiskMod, only: actuatorDisk
    use actuatorDisk_T2Mod, only: actuatorDisk_T2
    use actuatorDisk_RotMod, only: actuatorDisk_Rot
    use actuatorLineMod, only: actuatorLine
    use actuatorDisk_YawMod, only: actuatorDisk_yaw
    use dynamicYawMod, only: dynamicYaw
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
    character(len=clen) :: turbInfoDir

    integer :: ioUnit

    real(rkind), parameter :: degrees_to_radians = pi/180.0_rkind
    complex(rkind), parameter :: zeroC = zero + imi*zero

    type :: TurbineArray
        integer :: nTurbines
        integer :: myProc
        type(actuatorDisk), allocatable, dimension(:) :: turbArrayADM
        type(actuatorDisk_T2), allocatable, dimension(:) :: turbArrayADM_T2
        type(actuatorDisk_Rot), allocatable, dimension(:) :: turbArrayADM_Rot
        type(actuatorLine), allocatable, dimension(:) :: turbArrayALM
        type(actuatorDisk_yaw), allocatable, dimension(:) :: turbArrayADM_Tyaw
        type(dynamicYaw) :: dyaw

        type(decomp_info), pointer :: gpC, sp_gpC, gpE, sp_gpE
        type(spectral), pointer :: spectC, spectE
        type(staggOps), allocatable :: OpsNU
 
        real(rkind), dimension(:,:,:), allocatable :: fx, fy, fz
        complex(rkind), dimension(:,:,:), pointer :: fChat, fEhat, zbuffC, zbuffE
        real(rkind), dimension(:), allocatable :: gamma, theta, meanP, gamma_nm1, meanWs, meanPbaseline, stdP
        real(rkind), dimension(:), allocatable :: power_minus_n, ws_minus_n, pb_minus_n, hubDirection
        integer :: n_moving_average, timeStep, updateCounter 
        real(rkind), dimension(:,:), allocatable :: powerUpdate
        logical :: fixedYaw = .false., considerAdvection = .true., lookup = .false.
        integer :: dynamicStart = 1, hubIndex, dirType, advectionTime
        real(rkind) :: umAngle, vmAngle, windAngle, windAngle_old

        ! variables needed for halo communication
        integer :: neighbour(6), coord(2), dims(2), tag_s, tag_n, tag_b, tag_t
        real(rkind), allocatable, dimension(:,:,:) :: ySendBuf, zSendBuf, yRightHalo, zRightHalo, zLeftHalo

        real(rkind), dimension(:,:,:), allocatable :: rbuff, blanks, speed, scalarSource
        logical :: dumpTurbField = .false., useDynamicYaw, firstStep
        integer :: step = 0, ADM_Type, yawUpdateInterval 
        character(len=clen)                           :: powerDumpDir
        ! Variables to link domain and control
        real(rkind), dimension(:,:,:), pointer :: u_ref_sim, v_ref_sim
        real(rkind), dimension(:,:,:,:), pointer :: rbuffyC_ref_sim, rbuffzC_ref_sim
        type(decomp_info), pointer :: gpC_ref_sim
        logical :: ref_domain_link = .false.
    contains

        procedure :: init
        procedure :: destroy
        procedure :: init_halo_communication
        procedure :: halo_communication
        procedure :: destroy_halo_communication
        procedure :: getForceRHS 
        procedure :: reset_turbArray 
        !procedure :: dumpFullField
        procedure :: link_pointers
        procedure :: update_wind_angle
        procedure :: link_reference_domain_for_control
        procedure :: write_turbine_power

    end type

contains

subroutine link_reference_domain_for_control(this, u, v, rbuffyC, rbuffzC, gpC)
    class(TurbineArray), intent(inout) :: this
    real(rkind), dimension(:,:,:)  , intent(in), target :: u, v 
    real(rkind), dimension(:,:,:,:)  , intent(in), target :: rbuffyC, rbuffzC
    type(decomp_info), intent(in), target :: gpC

    this%u_ref_sim => u
    this%v_ref_sim => v
    this%rbuffyC_ref_sim => rbuffyC
    this%rbuffzC_ref_sim => rbuffzC
    this%gpC_ref_sim => gpC
    this%ref_domain_link = .true.


end subroutine

subroutine update_wind_angle(this)
    class(TurbineArray), intent(inout) :: this

    ! Added by Mike 10/16/19 to compute the average wind direction in the
    ! domain which is needed for the yaw misaligned turbines        igp%rbuffxC(:,:,:,1) = atan2(igp%v, igp%u) !* 180.d0 / 3.14d0
    if (this%ref_domain_link) then
        call transpose_x_to_y(this%u_ref_sim,this%rbuffyC_ref_sim(:,:,:,1),this%gpC_ref_sim)
        call transpose_y_to_z(this%rbuffyC_ref_sim(:,:,:,1),this%rbuffzC_ref_sim(:,:,:,1),this%gpC_ref_sim)
        call transpose_x_to_y(this%v_ref_sim,this%rbuffyC_ref_sim(:,:,:,1),this%gpC_ref_sim)
        call transpose_y_to_z(this%rbuffyC_ref_sim(:,:,:,1),this%rbuffzC_ref_sim(:,:,:,2),this%gpC_ref_sim)
        this%umAngle = p_sum(sum(this%rbuffzC_ref_sim(:,:,this%hubIndex,1))) / &
              (real(this%gpC_ref_sim%xsz(1),rkind) * real(this%gpC_ref_sim%ysz(2),rkind))
        this%vmAngle = p_sum(sum(this%rbuffzC_ref_sim(:,:,this%hubIndex,2))) / &
              (real(this%gpC_ref_sim%xsz(1),rkind) * real(this%gpC_ref_sim%ysz(2),rkind))
        this%windAngle = atan2(this%vmAngle,this%umAngle) * 180.d0 / pi
    else
        this%windAngle = 0.d0
    end if

end subroutine

subroutine link_pointers(this,fxturb, fyturb, fzturb)
    class(TurbineArray), intent(in), target :: this
    real(rkind), dimension(:,:,:)  , pointer, intent(inout) :: fxturb, fyturb, fzturb

    fxturb => this%fx
    fyturb => this%fy
    fzturb => this%fz

end subroutine 


subroutine reset_turbArray(this)
    class(TurbineArray), intent(inout), target :: this
    integer :: i

    !if(ADM) then
      do i = 1, this%nTurbines
        ! CHANGED to allow avoiding inst_horz_avg calculations - useful for
        ! testing/debugging
        !call this%turbArrayADM(i)%reset_turbine()
      end do 
    !endif

end subroutine


subroutine init(this, inputFile, gpC, gpE, spectC, spectE, cbuffyC, cbuffYE, cbuffzC, cbuffzE, mesh, dx, dy, dz)
    class(TurbineArray), intent(inout), target :: this
    character(len=*),    intent(in)            :: inputFile
    character(len=clen)                           :: powerDumpDir
    type(spectral), target, intent(in) :: spectC, spectE
    type(decomp_info), target, intent(in) :: gpC, gpE!, sp_gpC, sp_gpE
    !real(rkind), dimension(:,:,:,:), target, intent(inout) :: rbuffxC   ! actually 3 are required
    complex(rkind), dimension(:,:,:,:), target, intent(inout) :: cbuffyC, cbuffyE, cbuffzC, cbuffzE
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), intent(in) :: dx, dy, dz
    logical :: useWindTurbines = .TRUE., useDynamicYaw = .FALSE. ! .FALSE. implies ALM
    real(rkind) :: xyzPads(6)
    logical :: ADM = .TRUE., WriteTurbineForce  ! .FALSE. implies ALM
    ! Dynamic yaw stuff
    character(len=clen) :: inputDirDyaw = "/home1/05294/mhowland/dynamicYawFiles/dynamicYaw.inp"
    real(rkind), dimension(:), allocatable :: xLoc, yLoc
    integer :: yawUpdateInterval = 1000

    integer :: i, ierr, ADM_Type = 2

    namelist /WINDTURBINES/ useWindTurbines, num_turbines, ADM, turbInfoDir, ADM_Type, & 
                            WriteTurbineForce, powerDumpDir, useDynamicYaw, yawUpdateInterval, inputDirDyaw

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

    allocate(this%fx(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
    allocate(this%fy(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
    allocate(this%fz(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))

    this%fChat => cbuffyC(:,:,:,1); this%fEhat => cbuffyE(:,:,:,1)
    this%zbuffC => cbuffzC(:,:,:,1); this%zbuffE => cbuffzE(:,:,:,1)

    allocate(this%OpsNU)
    call this%OpsNU%init(this%gpC,this%gpE,0,dx,dy,dz,this%spectC%spectdecomp,this%spectE%spectdecomp,.true.,.true.)

    ! set number of turbines
    this%nTurbines = num_turbines;
    this%powerDumpDir = powerDumpDir
    this%useDynamicYaw = useDynamicYaw
    this%yawUpdateInterval = yawUpdateInterval
    this%hubIndex = 1

    ! Initialize the yaw and tilf
    allocate(this%gamma(this%nTurbines))
    allocate(this%gamma_nm1(this%nTurbines))
    allocate(this%theta(this%nTurbines))
    allocate(xLoc(this%nTurbines))
    allocate(yLoc(this%nTurbines))

    if(ADM) then
      this%ADM_Type = ADM_Type
      select case (ADM_Type)
      case (1)
         allocate (this%turbArrayADM(this%nTurbines))
         do i = 1, this%nTurbines
            call this%turbArrayADM(i)%init(turbInfoDir, i, mesh(:,:,:,1), mesh(:,:,:,2), mesh(:,:,:,3),this%gpC)
         end do
         call message(0,"WIND TURBINE ADM (Type 1) array initialized")
      case (2)
         allocate (this%turbArrayADM_T2(this%nTurbines))
         do i = 1, this%nTurbines
            call this%turbArrayADM_T2(i)%init(turbInfoDir, i, mesh(:,:,:,1), mesh(:,:,:,2), mesh(:,:,:,3))
         end do
         call message(0,"WIND TURBINE ADM (Type 2) array initialized")
      case (3)
         allocate (this%turbArrayADM_Rot(this%nTurbines))
         do i = 1, this%nTurbines
            call this%turbArrayADM_Rot(i)%init(turbInfoDir, i, mesh(:,:,:,1), mesh(:,:,:,2), mesh(:,:,:,3))
         end do
         call message(0,"WIND TURBINE ROT ADM (Type 3) array initialized")
      case (4)
         allocate (this%turbArrayADM_Tyaw(this%nTurbines))
         allocate (this%rbuff(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
         allocate (this%blanks(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
         allocate (this%speed(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
         allocate (this%scalarSource(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
         if (this%useDynamicYaw) then
             allocate(this%meanP(this%nTurbines))
             allocate(this%stdP(this%nTurbines))
             allocate(this%meanPbaseline(this%nTurbines))
             allocate(this%meanWs(this%nTurbines))
             allocate(this%hubDirection(this%nTurbines))
             allocate(this%powerUpdate(this%yawUpdateInterval, this%nTurbines))
             this%powerDumpDir = powerDumpDir
             this%timeStep = 1
             this%step = 0
             this%updateCounter = 1
             this%meanP = 0.d0
             this%stdP = 0.d0
             this%meanPbaseline = 0.d0
             this%meanWs = 0.d0
             this%powerUpdate = 0.d0
             this%hubDirection = 0.d0
             this%firstStep = .TRUE.
         end if
         do i = 1, this%nTurbines
             call this%turbArrayADM_Tyaw(i)%init(turbInfoDir, i, mesh(:,:,:,1), mesh(:,:,:,2), mesh(:,:,:,3))
             this%gamma(i) = this%turbArrayADM_Tyaw(i)%yaw
             this%theta(i) = 0.d0
             call this%turbArrayADM_Tyaw(i)%link_memory_buffers(this%rbuff, this%blanks, this%speed,  & 
                   this%scalarSource)
             xLoc(i) = this%turbArrayADM_Tyaw(i)%xLoc
             yLoc(i) = this%turbArrayADM_Tyaw(i)%yLoc
         end do
         this%gamma_nm1 = this%gamma
         this%hubIndex = nint(this%turbArrayADM_Tyaw(1)%zLoc / dz)
         this%windAngle = 0.d0
         call message(0,"YAWING WIND TURBINE (Type 4) array initialized")
      end select 
    else
      call GracefulExit("Actuator Line implementation temporarily disabled. Talk to Aditya if you want to know why.",423)
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

    if (this%useDynamicYaw) then
        call this%dyaw%init(inputDirDyaw, xLoc, yLoc, &
                            this%turbArrayADM_Tyaw(1)%diam, this%nTurbines, &
                            this%fixedYaw, this%dynamicStart, this%dirType, this%considerAdvection, this%lookup)
    end if

    !if(this%dumpTurbField) then
    !    this%fx = zero; this%fy = zero; this%fz = zero
    !    call this%dumpFullField(this%fx, "wtfx")
    !    call this%dumpFullField(this%fy, "wtfy")
    !    call this%dumpFullField(this%fz, "wtfz")
    !    this%dumpTurbField = .false.
    !endif

end subroutine


subroutine destroy(this)
    class(TurbineArray), intent(inout) :: this
    integer :: i

    nullify(this%gpC, this%gpE, this%spectC, this%sp_gpC)
    nullify(this%zbuffC, this%zbuffE, this%fChat, this%fEhat)
    deallocate(this%fx, this%fy, this%fz)

    deallocate(this%OpsNU)

    !if(ADM) then
    select case (this%ADM_Type)
    case (1)
      do i = 1, this%nTurbines
        call this%turbArrayADM(i)%destroy()
      end do
    case (2)
      do i = 1, this%nTurbines
        call this%turbArrayADM_T2(i)%destroy()
      end do
    case (3)
      do i = 1, this%nTurbines
        call this%turbArrayADM_Rot(i)%destroy()
      end do
    case (4)
      do i = 1, this%nTurbines
        call this%turbArrayADM_Tyaw(i)%destroy()
      end do
    end select
      !deallocate(this%turbArrayADM)
    !else
    !  call this%destroy_halo_communication()
    !  do i = 1, this%nTurbines
    !    call this%turbArrayALM(i)%destroy()
    !  end do
    !  deallocate(this%turbArrayALM)
    !endif

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

  xyzPads(5) = zG(1,1,1) - 0.5D0*dz
  xyzPads(6) = zG(1,1,this%gpC%xsz(3)) + 0.5D0*dz

end subroutine

subroutine halo_communication(this, u, v, wC)
  use decomp_2d
  use kind_parameters, only: mpirkind
  class(TurbineArray), intent(inout) :: this
  real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)),intent(in) :: u, v, wC

  integer :: k, icount, requests(2), ierror
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

subroutine getForceRHS(this, dt, u, v, wC, urhs, vrhs, wrhs, newTimeStep, inst_horz_avg, uturb, vturb, wturb)
    class(TurbineArray), intent(inout), target :: this
    real(rkind),                                                                         intent(in) :: dt
    real(rkind),    dimension(this%gpC%xsz(1),   this%gpC%xsz(2),   this%gpC%xsz(3)),    intent(in) :: u, v, wC
    complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout) :: urhs, vrhs
    complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: wrhs 
    logical,                                                                             intent(in)    :: newTimestep
    real(rkind),    dimension(:),                                                        intent(out)   :: inst_horz_avg
    complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout), optional :: uturb, vturb
    complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout), optional :: wturb
    integer :: i
    character(len=clen) :: tempname
    real(rkind) :: tmp
    real(rkind), dimension(this%nTurbines) :: angleIn

    if (newTimeStep) then
         this%fx = zero; this%fy = zero; this%fz = zero
           select case (this%ADM_Type)
           case(1)
              do i = 1, this%nTurbines
                    call this%turbArrayADM(i)%get_RHS(u,v,wC,this%fx,this%fy,this%fz,inst_horz_avg(8*i-7:8*i))
              end do
           case(2)
               do i = 1, this%nTurbines
                    call this%turbArrayADM_T2(i)%get_RHS(u,v,wC,this%fx,this%fy,this%fz, inst_horz_avg(8*i-7:8*i))
               end do
           case(3)
               do i = 1, this%nTurbines
                    call this%turbArrayADM_Rot(i)%get_RHS(u,v,wC,this%fx,this%fy,this%fz, inst_horz_avg(8*i-7:8*i))
               end do
           case (4)
               do i = 1, this%nTurbines
                   ! If it is the first time step, implement the wind direction
                   ! based yaw alignment
                   if (this%step==1) then
                      call this%update_wind_angle()
                      this%gamma = this%windAngle*pi/180.d0
                   end if
                   ! Get RHS
                   call this%turbArrayADM_Tyaw(i)%get_RHS_withPower(u,v,wC,this%fx,this%fy,this%fz, this%gamma(i), this%theta(i), this%windAngle, this%dirType)
                   ! Proceed with dynamic yaw operations
                   if (this%useDynamicYaw) then
                       if (this%step==1) then
                           if (this%considerAdvection) then
                               this%advectionTime = nint(2*(this%turbArrayADM_Tyaw(this%nTurbines)%xLoc - this%turbArrayADM_Tyaw(1)%xLoc) / dt)
                           else
                               this%advectionTime = 1
                           end if
                       end if
                       this%powerUpdate(this%timeStep, i) = & 
                                         this%turbArrayADM_Tyaw(i)%get_power()
                       if (this%timeStep >= this%advectionTime) then
                           call this%dyaw%simpleMovingAverage(this%meanP(i), &
                                this%turbArrayADM_Tyaw(i)%get_power(), this%meanWs(i), & 
                                this%turbArrayADM_Tyaw(i)%ut, &
                                this%meanPbaseline(i), this%turbArrayADM_Tyaw(i)%powerBaseline, &
                                this%hubDirection(i), this%turbArrayADM_Tyaw(i)%hubDirection, this%stdP(i), & 
                                this%timeStep - this%advectionTime, i)
                       end if
                   end if
               end do    
               if (this%useDynamicYaw) then
                   ! Update the yaw misalignment for each turbine
                   if (mod(this%timeStep, this%yawUpdateInterval) == 0 .and. this%timeStep /= 1) then
                       ! Update the wind angle measurement
                       this%windAngle_old = this%windAngle
                       call this%update_wind_angle()
                       if (this%dirType==1 .or. this%updateCounter==1) then
                           angleIn = this%gamma - this%windAngle*pi/180.d0
                           call this%dyaw%update_and_yaw(angleIn, this%meanWs(1), & 
                                                     this%windAngle, this%meanP, this%step, this%meanPbaseline)
                           this%gamma = angleIn
                       else
                           angleIn = this%gamma - this%hubDirection*pi/180.d0
                           call this%dyaw%update_and_yaw(angleIn, this%meanWs(1), & 
                                                     this%hubDirection(1), this%meanP, this%step, this%meanPbaseline)
                           if (this%lookup==.false.) then
                               this%gamma = angleIn
                           else
                               ! These angles were taken after one step of
                               ! online yaw, this is meant to simulate the
                               ! lookup table approach
                               this%gamma(1) = 0.2935d0; this%gamma(2) = 0.2973d0; this%gamma(3) = 0.2949d0;
                               this%gamma(4) = 0.2576d0; this%gamma(5) = 0.0490d0; this%gamma(6) = 0.0000d0;
                           end if
                       endif
                       ! Add the hub height wind direction to the yaw
                       ! misalignments
                       if (this%dirType==1) then
                           this%gamma = this%gamma + this%windAngle * pi / 180.d0
                       elseif (this%dirType==2) then
                           this%gamma = this%gamma + this%hubDirection * pi / 180.d0
                       endif
                       if ((this%fixedYaw) .or. (this%dynamicStart>this%step)) then
                           if (this%dirType==1) then
                               this%gamma = this%windAngle * pi / 180.d0
                           elseif (this%dirType==2) then
                               this%gamma = this%hubDirection * pi / 180.d0
                           endif
                       end if
                       do i=1,this%nTurbines
                           write(tempname,"(A12,I3.3,A8,I3.3,A4)") "powerUpdate_",i,"_update_",this%updateCounter,".txt"
                           call this%turbArrayADM_Tyaw(i)%dumpPowerUpdate(this%powerDumpDir,& 
                                tempname, this%powerUpdate(:,i), this%dyaw%Phat, & 
                                this%gamma, this%gamma_nm1, this%meanP, &
                                this%dyaw%kw, this%dyaw%sigma_0, &
                                this%dyaw%Phat_yaw, this%updateCounter, this%meanPbaseline, this%hubDirection, this%dyaw%Popti, this%stdP)
                       end do
                       this%timeStep = 0
                       this%updateCounter=this%updateCounter+1
                       this%powerUpdate = 0.d0
                   end if
                   ! Update time step
                   if (this%firstStep == .FALSE.) then
                       this%timeStep = this%timeStep+1 
                   else
                       ! Don't count the initialization step
                       this%firstStep = .FALSE.
                   end if
                   this%gamma_nm1 = this%gamma
                   this%step=this%step+1
               end if
           end select 
    end if 

    ! add forces to rhs
    call this%spectC%fft(this%fx,this%fChat)
    urhs = urhs + this%fChat
    if (present(uturb)) uturb = this%fChat

    call this%spectC%fft(this%fy,this%fChat)
    vrhs = vrhs + this%fChat
    if (present(vturb)) vturb = this%fChat

    call this%spectC%fft(this%fz,this%fChat)
    ! interpolate fz to fzE
    call transpose_y_to_z(this%fChat,this%zbuffC,this%sp_gpC)
    call this%OpsNU%InterpZ_Cell2Edge(this%zbuffC,this%zbuffE,zeroC,zeroC)
    call transpose_z_to_y(this%zbuffE,this%fEhat,this%sp_gpE)
    wrhs = wrhs + this%fEhat
    if (present(wturb)) wturb = this%fEhat

end subroutine 

subroutine write_turbine_power(this, TID, outputdir, runID)
    use basic_io
    class(TurbineArray), intent(inout), target :: this
    integer :: i, runID, TID
    character(len=*), intent(in) :: outputdir
    character(len=clen) :: filename, tempname

    do i = 1,this%nTurbines
        if (this%ADM_Type==2) then
           if (allocated(this%turbArrayADM_T2(i)%powerTime)) then
               write(tempname,"(A3,I2.2,A2,I6.6,A6,I2.2,A4)") "Run",runID,"_t", TID, "_turbP",i,".pow"
               filename = outputDir(:len_trim(outputDir))//"/"//trim(tempname)
               call write_2d_ascii(this%turbArrayADM_T2(i)%powerTime(1:this%turbArrayADM_T2(i)%tInd-1,:),filename)  
               this%turbArrayADM_T2(i)%tInd = 1
           end if
        end if
    end do
    
end subroutine

end module
