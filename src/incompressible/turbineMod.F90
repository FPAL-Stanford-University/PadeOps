module turbineMod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa
    use decomp_2d
    use StaggOpsMod, only: staggOps  
    use actuatorDiskMod, only: actuatorDisk
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

        !type(actuatorLine), allocatable :: turbArrayALM

        !integer, dimension(:), allocatable  :: xst, xen, yst, yen
        type(decomp_info), pointer :: gpC, sp_gpC, gpE, sp_gpE
        type(spectral), pointer :: spectC, spectE
        type(staggOps), allocatable :: OpsNU
        !integer :: myLeftNeigh, myRightNeigh, myTopNeigh, myBotNeigh
 
        !integer, dimension(:),   allocatable :: num_cells_cloud                                     ! total number of cells in the cubic cloud around a turbine on this processor
        !integer, dimension(:),   allocatable :: num_blades  ! number of blades
        !integer, dimension(:), allocatable :: num_blade_points  ! number of actuator points on each blade
        !real(rkind), dimension(:), allocatable :: yaw_angle, blade_azimuth, nacelle_width, hub_radius, tip_radius, turb_thrust, turb_torque, rotspeed
        !real(rkind), dimension(:,:), allocatable :: rotor_center, turbLoc, rotor_shaft
        !real(rkind), dimension(:,:,:,:), allocatable :: blade_points  ! number of actuator points on each blade
        !real(rkind), dimension(:,:,:,:), allocatable :: blade_forces  ! forces at actuator points
        !logical, dimension(:), allocatable :: clockwise_rotation

        real(rkind), dimension(:,:,:), allocatable :: local_rbuffxC
        real(rkind), dimension(:,:,:), pointer :: fx, fy, fz
        complex(rkind), dimension(:,:,:), pointer :: fChat, fEhat, zbuffC, zbuffE

    contains

        procedure :: init
        procedure :: destroy
        procedure :: getForceRHS 

    end type

contains

subroutine init(this, inputFile, gpC, gpE, spectC, spectE, rbuffxC, cbuffyC, cbuffYE, cbuffzC, cbuffzE, mesh, dx, dy, dz)
    class(TurbineArray), intent(inout), target :: this
    character(len=*), intent(in) :: inputFile
    type(spectral), target :: spectC, spectE
    type(decomp_info), target :: gpC, gpE!, sp_gpC, sp_gpE
    real(rkind), dimension(:,:,:,:), target :: rbuffxC   ! actually 3 are required
    complex(rkind), dimension(:,:,:,:), target :: cbuffyC, cbuffyE, cbuffzC, cbuffzE
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), intent(in) :: dx, dy, dz
    logical :: useWindTurbines = .TRUE. ! .FALSE. implies ALM

    integer :: i

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

    !allocate(this%cbuffC(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), this%sp_gpC%ysz(3), 1))
    !allocate(this%cbuffE(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), this%sp_gpE%ysz(3), 1))

    allocate(this%local_rbuffxC(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))
    !allocate(this%rbuffE(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3), 1))

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
      !call this%turbArrayALM%init(this%nTurbines)
      !call message(1,"WIND TURBINE ALM model initialized")
      call GracefulExit("Wind Turbine ALM stuff is incomplete", 423)
    endif

end subroutine


subroutine destroy(this)
    class(TurbineArray), intent(inout) :: this
    integer :: i

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
      !call this%turbArrayALM%destroy()
    endif

end subroutine

subroutine getForceRHS(this, dt, u, v, wC, urhs, vrhs, wrhs, inst_horz_avg)
    class(TurbineArray), intent(inout), target :: this
    real(rkind),                                                                         intent(in) :: dt
    real(rkind),    dimension(this%gpC%xsz(1),   this%gpC%xsz(2),   this%gpC%xsz(3)),    intent(in) :: u, v, wC
    complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout) :: urhs, vrhs
    complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: wrhs 
    real(rkind),    dimension(:),                                                        intent(out), optional   :: inst_horz_avg
    integer :: i

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
      !call this%turbArrayALM%get_RHS(dt, u, v, wC, this%fx, this%fy, this%fz)
        print*, dt ! Temporary placeholder to avoid REMARK messages
    endif

    ! add forces to rhs
    call this%spectC%fft(this%fx,this%fChat)
    urhs = urhs - this%fChat

    call this%spectC%fft(this%fy,this%fChat)
    vrhs = vrhs - this%fChat

    call this%spectC%fft(this%fz,this%fChat)
    ! interpolate fz to fzE
    call transpose_y_to_z(this%fChat,this%zbuffC,this%sp_gpC)
    call this%OpsNU%InterpZ_Cell2Edge(this%zbuffC,this%zbuffE,zeroC,zeroC)
    call transpose_z_to_y(this%zbuffE,this%fEhat,this%sp_gpE)
    wrhs = wrhs - this%fEhat

end subroutine 

end module
