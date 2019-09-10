module dynamicYawMod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa
    use decomp_2d
    use exits, only: GracefulExit, message
    use mpi 
    use reductions, only: p_maxval, p_sum
    use timer, only: tic, toc
    use Gridtools, only: linspace

    implicit none

    private
    public :: dynamicYaw
    
    real(rkind), parameter :: alpha_Smooth = 1.d0 ! 0.9d0 ! Exonential smoothing constant
    integer, parameter :: xReg = 8, yReg = 8, zReg = 8

    type :: dynamicYaw
        ! Optimization parameters
        real(rkind) :: beta1 = 0.9d0, beta2 = 0.999d0, learning_rate_yaw = 10D-4
        integer :: epochsYaw = 5000, stateEstimationEpochs = 500, Ne = 20
        real(rkind) :: var_p = 0.04d0, var_k = 9D-6, var_sig = 9D-6
        ! Turbine parameters
        integer :: Nt = 1
        real(rkind), dimension(:,:), allocatable :: turbCenter
        real(rkind) :: diam=0.315d0
        real(rkind) :: yaw
        ! Wake model parameters
        real(rkind), dimension(:), allocatable :: kw, sigma_0
        real(rkind) :: powerExp, eps=1.0D-12
        integer :: Nx
        ! Wind conditions
        real(rkind) :: wind_speed, wind_direction

    contains
        procedure :: init
        procedure :: destroy
        procedure, private :: forward
        procedure, private :: EnKF_step
        procedure, private :: backward
        procedure :: yawOptimize 
    end type


contains

subroutine init(this, inputDir)
    class(actuatorDisk_yaw), intent(inout) :: this
    character(len=*), intent(in) :: inputDir
    character(len=clen) :: tempname, fname
    integer :: ioUnit
    real(rkind) :: beta1 = 0.9d0, beta2 = 0.999d0, learning_rate_yaw = 10D-4
    integer :: epochsYaw = 5000, stateEstimationEpochs = 500, Ne = 20, Nt = 2
    real(rkind) :: var_p = 0.04d0, var_k = 9D-6, var_sig = 9D-6
 
    ! Read input file for this turbine    
    namelist /DYNAMIC_YAW/ Nt, var_p, var_k, var_sig, epochsYaw, stateEstimationEpochs, Ne
    ioUnit = 534
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
    read(unit=ioUnit, NML=DYNAMIC_YAW)
    close(ioUnit)  

    this%Nt = Nt; this%var_p = var_p; this%var_k = var_k; this%var_sig = var_sig;
    this%epochsYaw = epochsYaw; this%stateEstimationEpochs = stateEstimationEpochs;
    this%Ne = Ne
 
end subroutine 

subroutine destroy(this)
    class(dynamicYaw), intent(inout) :: this

end subroutine 

subroutine yawOptimize(this)
    class(dynamicYaw), intent(inout) :: this

    ! Field data observation
    call this%observeField()
    ! Rotate domain
    [ XR, indSorted, unsort ] = rotate( X, a, conditionTurb );
    params.P = params.P(indSorted) / params.P(conditionTurb);
    turbine.turbCenter = XR;
    % Online control update
    [yaws, P_opti, P_baseline, params, errorOut{t}] = onlineUpdate(turbine, atm, params, params_opti); 
    % Store
    yaw_store(:,t) = yaws(unsort);
    P_opti_store(:,t) = P_opti(unsort);
    P_baseline_store(:,t) = P_baseline(unsort);

end

end subroutine

subroutine observeField(this)
    class(dynamicYaw), intent(inout) :: this
    this%wind_speed = 8.d0
    this%wind_direction = 0.d0
    this%powerObservation = 1.d0 ! Get the power production from ADM code 

end module 
