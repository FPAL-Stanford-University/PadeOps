module dynamicYawMod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa
    use decomp_2d
    use exits, only: GracefulExit, message
    use mpi 
    use reductions, only: p_maxval, p_sum
    use timer, only: tic, toc
    use Gridtools, only: linspace
    use sorting_mod, only: Qsort

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
        real(rkind), dimension(:), allocatable :: powerObservation
        integer :: conditionTurb
        ! Other stuff
        integer, dimension(:), allocatable :: indSorted, unsort

    contains
        procedure :: init
        procedure :: destroy
        procedure, private :: forward
        procedure, private :: EnKF_update
        procedure, private :: backward
        procedure, private :: observeField
        procedure, private :: rotate
        procedure, private :: yawOptimize 
        procedure, private :: onlineUpdate 
        procedure :: update_and_yaw
    end type


contains

subroutine init(this, inputfile)
    class(dynamicYaw), intent(inout) :: this
    character(len=*), intent(in) :: inputfile
    character(len=clen) :: tempname, fname
    integer :: ioUnit, conditionTurb, ierr
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
    this%conditionTurb = 1
    
    ! Allocate
    allocate(this%indSorted(this%Nt))
    allocate(this%unsort(this%Nt))
    allocate(this%powerObservation(this%Nt))

 
end subroutine 

subroutine destroy(this)
    class(dynamicYaw), intent(inout) :: this

end subroutine 

subroutine update_and_yaw(this)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(:,:), allocatable :: X

    ! allocate
    allocate(X(this%Nt,2))

    ! Field data observation
    call this%observeField()
    ! Rotate domain
    X = this%turbCenter
    call this%rotate(X)
    ! Normalize the power production
    this%powerObservation = this%powerObservation(this%indSorted) / this%powerObservation(this%conditionTurb)
    ! Online control update
    call this%onlineUpdate(X)

    deallocate(X)

end subroutine

subroutine onlineUpdate(this, X)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(:,:), intent(inout) :: X
    real(rkind), dimension(:), allocatable :: kw, sigma_0, Phat, kwBest, sigmaBest
    real(rkind), dimension(:,:), allocatable :: psi_kp1
    real(rkind), dimension(this%stateEstimationEpochs) :: error
    real(rkind) :: lowestError
    integer :: nt2, i, t

    ! Allocate
    allocate(kw(this%Nt))
    allocate(kwBest(this%Nt))
    allocate(sigma_0(this%Nt))
    allocate(sigmaBest(this%Nt))
    allocate(Phat(this%Nt))
    nt2 = this%Nt*2
    allocate(psi_kp1(nt2,this%Ne))
 
    ! State estimation
    error = 0.d0; lowestError = 100.d0
    do i=1, this%Ne
        psi_kp1(1:this%Nt,i) = reshape(this%kw, (/this%Nt,1/))
        psi_kp1(this%Nt:nt2,i) = reshape(this%sigma_0, (/this%Nt,1/))
    end do
    bestParams = params; lowestError = 10; bestStep = 1;
    do t=1, this%stateEstimationEpochs;
        call this%EnKF_update( psi_kp1, P_kp1, kw, sigma_0, psi_kp1 )
        ! Run forward model pass
        call this%forward(kw, sigma_0, Phat)
        ! Normalized
        Phat = Phat / Phat(1)
        error(t) = p_sum(abs(Phat-P_kp1)/Phat(1)) / real(this%Nt)
        if (error(t)<lowestError) then
            lowestError=error(t); bestStep = t;
            kwBest = kw; sigmaBest = sigma_0
        end if
    end do

    ! Generate optimal yaw angles
    call this%yawOptimize(kwBest, sigmaBest)
    ! Store the unsorted parameters
    this%kw = kwBest(this%unsort); this%sigma_0 = sigmaBest(this%unsort)   

    ! Deallocate
    deallocate(kw)
    deallocate(kwBest)
    deallocate(sigmaBest)
    deallocate(sigma_0)
    deallocate(Phat)
 
end subroutine

subroutine EnKF_update(this)
    class(dynamicYaw), intent(inout) :: this

end subroutine

subroutine yawOptimize(this)
    class(dynamicYaw), intent(inout) :: this

end subroutine

subroutine observeField(this)
    class(dynamicYaw), intent(inout) :: this
    this%wind_speed = 8.d0
    this%wind_direction = 0.d0
    this%powerObservation = 1.d0 ! Get the power production from ADM code 

end subroutine

subroutine rotate(this, X)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(2,2) ::  R
    real(rkind), dimension(:,:), intent(inout) :: X
    integer, dimension(:), allocatable :: indSorted
    type group
        integer :: order    ! original order of unsorted data
        real :: value       ! values to be sorted by
    end type group
    type(group), dimension(1) :: x_for_sorting
    integer :: i
 
    a = (this%wind_direction - 270.d0) * pi / 180.d0
    R = reshape((/cos(a), sin(a), -sin(a), cos(a)/), shape(R))
    do i=1,this%Nt
        X(i,:) = R * transpose(this%turbCenter(i,:))
        x_for_sorting(i)%value = X(i,1)
        x_for_sorting(i)%order = i
        this%unsort(i) = i
    end do

    ! Sort the turbines by upwind location
    call Qsort(x_for_sorting,1)
    this%indSorted = x_for_sorting%order 
    this%unsort(indSorted) = this%unsort      

end subroutine

subroutine forward(this)
    class(dynamicYaw), intent(inout) :: this

end subroutine

subroutine backward(this)
    class(dynamicYaw), intent(inout) :: this

end subroutine

end module 
