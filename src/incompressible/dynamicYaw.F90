module dynamicYawMod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa
    use decomp_2d
    use exits, only: GracefulExit, message
    use mpi 
    use reductions, only: p_maxval, p_sum
    use timer, only: tic, toc
    use Gridtools, only: linspace, mytrapz
    use sorting_mod
    use random,             only: gaussian_random
    use actuatorDisk_YawMod, only: actuatorDisk_yaw

    implicit none

    private
    public :: dynamicYaw

    type :: dynamicYaw
        ! Optimization parameters
        real(rkind) :: beta1 = 0.9d0, beta2 = 0.999d0, learning_rate_yaw = 10.D-4
        integer :: epochsYaw = 5000, stateEstimationEpochs = 500, Ne = 20
        real(rkind) :: var_p = 0.04d0, var_k = 9D-6, var_sig = 9D-6
        ! Turbine parameters
        integer :: Nt = 1, ts
        real(rkind), dimension(:,:), allocatable :: turbCenter, turbCenterStore
        real(rkind) :: D=0.315d0
        real(rkind), dimension(:), allocatable :: yaw
        ! Wake model parameters
        real(rkind), dimension(:), allocatable :: kw, sigma_0, kw_initial, sigma_initial
        real(rkind) :: kw_init, sigma_init
        real(rkind) :: powerExp, eps=1.0D-12, Ct, eta ! Power exponent is 3 for this AD implementation
        integer :: Nx
        ! Uncertain wake model parameter matrices
        real(rkind), dimension(:,:), allocatable :: kw_u, sigma_u, kw_u_initial, sigma_u_initial ! size: [Num_turbines, power probability bins]
        real(rkind), dimension(:,:), allocatable :: Phat_fit_u ! size: [Num_turbines, power probability bins]
        integer :: p_bins, ws_bins, dir_bins, gamma_bins
        real(rkind), dimension(:), allocatable :: p_steps
        real(rkind) :: dir_std, gamma_std, ws_std, leading_Pstd, MaxModelError
        ! Wind conditions
        real(rkind) :: wind_speed, wind_direction, wind_direction_next
        real(rkind), dimension(:), allocatable :: powerObservation, powerBaseline, Popti, Pstd
        integer :: conditionTurb
        ! Other stuff
        integer, dimension(:), allocatable :: indSorted, unsort
        real(rkind) :: Ptot_initial, Ptot_final
        ! Cache stuff
        real(rkind), dimension(:), allocatable :: u_eff, Cp, a_mom
        real(rkind), dimension(:), allocatable :: ap, delta_u_face_store
        real(rkind), dimension(:,:), allocatable :: gaussianStore, y_c, yCenter, erfStore
        real(rkind), dimension(:,:), allocatable :: deltaUIndividual, sigmaEnd
        logical, dimension(:,:), allocatable :: turbinesInFront
        real(rkind) :: A, rho
        ! Backward stuff
        real(rkind), dimension(:), allocatable :: dp_dgamma 
        ! Moving average stuff
        integer :: n_moving_average
        real(rkind), dimension(:,:), allocatable :: power_minus_n, ws_minus_n, pb_minus_n, dir_minus_n
        real(rkind), dimension(:), allocatable :: Phat, Phat_yaw, Phat_fit
        logical :: check = .true., useInitialParams = .false., secondary = .false., uncertain = .false.
        logical :: ref_turbine
        ! Superposition stuff
        ! Superpostion: 1) Linear 2) Momentum conserving 3) Modified linear
        integer :: superposition = 1, ucMaxIt = 10, Ny = 100
        real(rkind) :: epsUc = 1.0D-3
        real(rkind), dimension(:,:), allocatable :: ucr
        real(rkind), dimension(:), allocatable :: v
        ! Wind direction statistics
        integer :: Tf_init, Tf_min
        logical :: use_alpha_check
        real(rkind) :: dalpha_max

    contains
        procedure :: init
        procedure :: destroy
        procedure, private :: forward
        procedure, private :: EnKF_update
        procedure, private :: observeField
        procedure, private :: rotate
        procedure, private :: rotate_local
        procedure, private :: yawOptimize 
        procedure, private :: yawOptimize_uncertain 
        procedure, private :: onlineUpdate 
        procedure, private :: backward
        procedure :: update_and_yaw
        procedure :: simpleMovingAverage
        procedure :: alpha_check
    end type


contains

subroutine init(this, inputfile, xLoc, yLoc, diam, Nt, fixedYaw, dynamicStart, dirType, considerAdvection, lookup)
    class(dynamicYaw), intent(inout) :: this
    character(len=*), intent(in) :: inputfile
    integer :: ioUnit, conditionTurb, ierr, i, n_moving_average, superposition
    real(rkind) :: beta1 = 0.9d0, beta2 = 0.999d0, learning_rate_yaw = 10D-4
    integer :: epochsYaw = 5000, stateEstimationEpochs = 500, Ne = 20, Nx = 500
    integer :: p_bins = 0, dir_bins = 0, gamma_bins = 0, ws_bins = 0
    real(rkind) :: var_p = 0.04d0, var_k = 9D-6, var_sig = 9D-6, D = 0.315
    real(rkind) :: Ct = 0.75, eta = 0.7, powerExp = 3.d0
    real(rkind), dimension(:), intent(in) :: xLoc, yLoc
    real(rkind), intent(in) :: diam
    integer, intent(in) :: Nt
    logical, intent(out) :: fixedYaw, considerAdvection, lookup
    integer, intent(out) :: dynamicStart, dirType
    logical :: useInitialParams, secondary, uncertain, ref_turbine, use_alpha_check
    real(rkind) :: dir_std, gamma_std, ws_std, MaxModelError, dalpha_max

    ! Read input file for this turbine    
    namelist /DYNAMIC_YAW/ var_p, var_k, var_sig, epochsYaw, stateEstimationEpochs, & 
                           Ne, Ct, eta, beta1, beta2, conditionTurb, n_moving_average, &
                           fixedYaw, dynamicStart, dirType, considerAdvection, lookup, &
                           useInitialParams, powerExp, superposition, secondary, uncertain, &
                           p_bins, dir_bins, gamma_bins, ws_bins, ref_turbine, &
                           dir_std, gamma_std, ws_std, MaxModelError, use_alpha_check, dalpha_max
    ioUnit = 534
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
    read(unit=ioUnit, NML=DYNAMIC_YAW)
    close(ioUnit)  

    this%Nt = Nt; this%var_p = var_p; this%var_k = var_k; this%var_sig = var_sig;
    this%epochsYaw = epochsYaw; this%stateEstimationEpochs = stateEstimationEpochs;
    this%Ne = Ne
    this%conditionTurb = conditionTurb
    this%Nx = Nx ! for spatial integration of the spanwise velocity
    this%D = diam
    this%Ct = Ct
    this%eta = eta; this%beta1 = beta1; this%beta2= beta2;
    this%learning_rate_yaw = learning_rate_yaw   
    this%n_moving_average = n_moving_average
    this%useInitialParams = useInitialParams 
    this%powerExp = powerExp
    this%superposition = superposition ! 1) linear, 2) convective (Zong & Porte-Agel 2020), 3) mod linear
    this%secondary = secondary ! Howland & Dabiri Energies 2020
    this%uncertain = uncertain ! Howland ACC 2020 (hopefully... :)
    this%p_bins = p_bins
    this%dir_bins = dir_bins
    this%gamma_bins = gamma_bins
    this%ws_bins = ws_bins
    this%ref_turbine = ref_turbine
    this%gamma_std = gamma_std
    this%dir_std = dir_std
    this%ws_std = ws_std
    this%MaxModelError = MaxModelError
    this%use_alpha_check = use_alpha_check
    this%dalpha_max = dalpha_max ! degrees

    ! Is there an adjacent reference turbine which we are not yawing?
    if (this%ref_turbine) then
        this%Nt = this%Nt - 1
    end if

    ! Allocate
    allocate(this%yaw(this%Nt))
    allocate(this%dp_dgamma(this%Nt))
    allocate(this%indSorted(this%Nt))
    allocate(this%Cp(this%Nt))
    allocate(this%a_mom(this%Nt))
    allocate(this%unsort(this%Nt))
    allocate(this%powerObservation(this%Nt))
    allocate(this%powerBaseline(this%Nt))
    allocate(this%u_eff(this%Nt))
    allocate(this%sigmaEnd(this%Nt, this%Nt))
    allocate(this%ap(this%Nt))
    allocate(this%delta_u_face_store(this%Nt))
    allocate(this%gaussianStore(this%Nt,this%Nt))
    allocate(this%erfStore(this%Nt,this%Nt))
    allocate(this%deltaUIndividual(this%Nt,this%Nt))
    allocate(this%turbinesInFront(this%Nt,this%Nt))
    allocate(this%ucr(this%Nt,this%Nt))
    allocate(this%y_c(this%Nt,this%Nt))
    allocate(this%yCenter(this%Nt,this%Nt))
    allocate(this%turbCenter(this%Nt,2))
    allocate(this%turbCenterStore(this%Nt,2))
    allocate(this%kw(this%Nt))
    allocate(this%sigma_0(this%Nt))
    allocate(this%kw_initial(this%Nt))
    allocate(this%sigma_initial(this%Nt))
    allocate(this%Phat(this%Nt))
    allocate(this%Phat_fit(this%Nt))
    allocate(this%Phat_yaw(this%Nt))
    allocate(this%Popti(this%Nt))
    allocate(this%v(this%Nt))
    ! Uncertain allocations
    allocate(this%kw_u(this%Nt, this%p_bins))
    allocate(this%sigma_u(this%Nt, this%p_bins))
    allocate(this%Phat_fit_u(this%Nt, this%p_bins))
    allocate(this%p_steps(this%p_bins))

    ! Define
    this%yaw = 0.d0
    this%kw_init = 0.1d0;
    this%sigma_init = 0.25;
    this%kw = this%kw_init; 
    this%sigma_0 = this%sigma_init; 
    this%kw_u = this%kw_init;
    this%sigma_u = this%sigma_init;
    this%power_minus_n = 0.d0
    this%ws_minus_n = 0.d0
    this%pb_minus_n = 0.d0
    ! Wind direction statistics
    this%Tf_init = this%n_moving_average; 
    this%Tf_min = this%n_moving_average/4.d0

    ! Uncertain
    ! Discretize between -1 STD and +1 STD
    if (this%p_bins>1) then
        this%p_steps = linspace(-1.d0, 1.d0, this%p_bins)
    else
        this%p_steps = 0.d0
    end if

    ! Get the wind turbine locations
    if (this%ref_turbine) then
        ! It is assumed that the reference turbine is the first actuator disk
        ! model initiated
        do i=1,this%Nt
            this%turbCenter(i,1) = xLoc(i+1) / this%D
            this%turbCenter(i,2) = yLoc(i+1) / this%D
        end do
        this%turbCenterStore = this%turbCenter
        ! Handle modified allocations for change to Nt
        allocate(this%power_minus_n(this%n_moving_average,this%Nt+1))
        allocate(this%ws_minus_n(this%n_moving_average,this%Nt+1))
        allocate(this%pb_minus_n(this%n_moving_average,this%Nt+1))
        allocate(this%dir_minus_n(this%n_moving_average,this%Nt+1))
    else    
        do i=1,this%Nt
            this%turbCenter(i,1) = xLoc(i) / this%D
            this%turbCenter(i,2) = yLoc(i) / this%D
        end do
        this%turbCenterStore = this%turbCenter
        ! Allocate
        allocate(this%power_minus_n(this%n_moving_average,this%Nt))
        allocate(this%ws_minus_n(this%n_moving_average,this%Nt))
        allocate(this%pb_minus_n(this%n_moving_average,this%Nt))
        allocate(this%dir_minus_n(this%n_moving_average,this%Nt))
    end if

end subroutine 

subroutine destroy(this)
    class(dynamicYaw), intent(inout) :: this

end subroutine 

subroutine update_and_yaw(this, yaw, wind_speed, wind_direction, wind_direction_next, powerObservation, ts, powerBaseline, Pstd, dirStd)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(:), intent(inout) :: yaw
    real(rkind), dimension(:), intent(in) :: powerObservation, powerBaseline, Pstd
    real(rkind), intent(in) :: wind_speed, wind_direction, dirStd, wind_direction_next
    integer, intent(in) :: ts

    ! Field data observation
    this%wind_speed = wind_speed ! Get this from the data
    this%wind_direction = 270.d0-wind_direction ! Get this from the data
    this%wind_direction_next = 270.d0-wind_direction_next ! Get this from the data
    ! Use a locally linear interpolation to get the wind direction at the end of
    ! the step
    ! If one of the turbines is a reference turbine modify the power observation
    ! to only consider turbines 2:Nt+1
    if (this%ref_turbine) then    
        this%powerObservation = powerObservation(2:this%Nt+1)
        this%Pstd = Pstd(2:this%Nt+1)
        this%yaw = yaw(2:this%Nt+1) 
        this%conditionTurb = 1 ! This method assumes that ADM 1 (the first actuator disk) is the reference
    else
        this%powerObservation = powerObservation ! Get the power production from ADM code
        this%Pstd = Pstd; ! Get power std from ADM code
        this%yaw = yaw  
    end if
    this%powerBaseline = powerBaseline ! Get the power production from ADM code
    this%ts = ts
    this%dir_std = dirStd

    ! Use initial fit parameters
    if (.not. (this%check)) then
        ! For WES P1 CNBL, this was this%kw_initial and this%sigma_initial
        ! (meaning initial fit not initialization values)
        ! Modified for diurnal case, leave as init condition (not after first
        ! fit)
        this%kw = this%kw_init
        this%sigma_0 = this%sigma_init
        this%kw_u = this%kw_init
        this%sigma_u = this%sigma_init
    end if

    ! Rotate domain
    call this%rotate()

    ! Normalize the power production by baseline power (computed using ADM)
    this%powerObservation = this%powerObservation(this%indSorted) / this%powerBaseline(this%conditionTurb)
    this%Pstd = this%Pstd(this%indSorted) / this%powerBaseline(this%conditionTurb);
    this%leading_Pstd = this%Pstd(1)
    this%Pstd(1) = 0.d0

    ! Online control update
    call this%onlineUpdate()
    if (this%ref_turbine) then    
        yaw(2:this%Nt+1) = this%yaw;
        yaw(1) = 0.d0; ! reference turbine needs zero yaw misalignment
    else
        yaw = this%yaw
    end if

    ! Restore original layout
    this%turbCenter = this%turbCenterStore

    ! Save the initial fit for kw and sigma
    if (this%check) then
        this%kw_initial = this%kw;
        this%sigma_initial = this%sigma_0;
        this%kw_u_initial = this%kw_u;
        this%sigma_u_initial = this%sigma_u;
        this%check = .false.
    end if

end subroutine

subroutine onlineUpdate(this)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(:), allocatable :: kw, sigma_0, Phat, kwBest, sigmaBest, yaw
    real(rkind), dimension(:), allocatable :: Phat_zeroYaw, p_steps
    real(rkind), dimension(:,:), allocatable :: psi_kp1, kwBest_u, sigmaBest_u
    real(rkind), dimension(this%stateEstimationEpochs) :: error
    real(rkind) :: lowestError, leading_turbine_err
    integer :: nt2, i, t, bestStep, pint

    ! Allocate
    allocate(kw(this%Nt))
    allocate(yaw(this%Nt))
    allocate(kwBest(this%Nt))
    allocate(sigma_0(this%Nt))
    allocate(sigmaBest(this%Nt))
    allocate(Phat(this%Nt))
    allocate(Phat_zeroYaw(this%Nt))
    nt2 = this%Nt*2
    allocate(psi_kp1(nt2,this%Ne))

    if (.not.(this%uncertain)) then ! Certain wake model parameters
        ! Store local parameter and rotate
        kw = this%kw(this%indSorted); sigma_0 = this%sigma_0(this%indSorted)
        yaw = this%yaw(this%indSorted)
        call this%forward(kw, sigma_0, Phat_zeroYaw, yaw*0.d0)
        ! State estimation
        error = 0.d0; lowestError = 100.d0
        do i=1, this%Ne
            do t=1,this%Nt
                psi_kp1(t,i) = kw(t)
                psi_kp1(this%Nt+t,i) = sigma_0(t)
            end do
        end do
        bestStep = 1; kwBest = kw; sigmaBest = sigma_0;
        do t=1, this%stateEstimationEpochs;
            call this%EnKF_update( psi_kp1, this%powerObservation, kw, sigma_0, yaw, t)
            ! Run forward model pass
            call this%forward(kw, sigma_0, Phat, yaw)
            ! Normalized
            Phat = Phat / Phat_zeroYaw(1)
            error(t) = sum(abs(Phat-this%powerObservation)) / real(this%Nt,rkind)
            if (error(t)<lowestError) then
                lowestError=error(t); bestStep = t;
                kwBest = kw; sigmaBest = sigma_0
                this%Ptot_initial = sum(Phat)
            end if
        end do
        this%Phat_fit = Phat(this%unsort)

        ! Rotate to future wind direction
        kwBest = kwBest(this%unsort); ! unsort
        sigmaBest = sigmaBest(this%unsort) ! unsort
        this%wind_direction = this%wind_direction_next ! set next step wind direction forecast
        this%turbCenter = this%turbCenterStore ! restore original layout
        call this%rotate() ! rotate
        kwBest = kwBest(this%indSorted); ! sort
        sigmaBest = sigmaBest(this%indSorted); ! sort

        ! Generate optimal yaw angles
        this%kw = kwBest; this%sigma_0 = sigmaBest;
        ! If useInitialParams only use initial fit to optimize model
        if (this%useInitialParams .and. (.not. this%check)) then
            kwBest = this%kw_initial(this%indSorted)
            sigmaBest = this%sigma_initial(this%indSorted)
        end if
        call this%yawOptimize(kwBest, sigmaBest, yaw)
        ! Store the unsorted parameters
        this%kw = kwBest(this%unsort); this%sigma_0 = sigmaBest(this%unsort)   
        this%yaw = yaw(this%unsort) 
        this%Phat = this%Phat(this%unsort); this%Phat_yaw = this%Phat_yaw(this%unsort)
    

    else ! Wake model parameter uncertainty
        ! Initilize parameter matrix
        allocate(kwBest_u(this%Nt, this%p_bins))
        allocate(sigmaBest_u(this%Nt, this%p_bins))
        kwBest_u = this%kw_u; sigmaBest_u = this%sigma_u;
        ! Check wake model error to determinite if the parameters need to be
        ! updated
        kw = this%kw_u(this%indSorted, (this%p_bins+1)/2); 
        sigma_0 = this%sigma_u(this%indSorted,(this%p_bins+1)/2)
        call this%forward(kw, sigma_0, Phat, yaw)
        call this%forward(kw, sigma_0, Phat_zeroYaw, yaw*0.d0)
        error(1) = sum(abs(Phat-this%powerObservation)) / real(this%Nt,rkind)
        leading_turbine_err = abs(Phat(1)/Phat_zeroYaw(1) - this%powerObservation(1))

        ! Check the fitting error and the quality of the input data
        ! If the leading turbine normalized power differs significantly from the
        ! wake model, don't use the data for updating wake model parameters
        if (error(1)>this%MaxModelError) then ! .and. leading_turbine_err<this%leading_Pstd) then
            ! Update wake model parameters
            this%Phat_fit_u = 0.d0
            do pint = 1, this%p_bins
                ! Store local parameter and rotate
                kw = this%kw_u(this%indSorted,pint); sigma_0 = this%sigma_u(this%indSorted,pint)
                yaw = this%yaw(this%indSorted)
                call this%forward(kw, sigma_0, Phat_zeroYaw, yaw*0.d0)
                ! State estimation
                error = 0.d0; lowestError = 100.d0
                do i=1, this%Ne
                    do t=1,this%Nt
                        psi_kp1(t,i) = kw(t)
                        psi_kp1(this%Nt+t,i) = sigma_0(t)
                    end do
                end do
                bestStep = 1; 
                kwBest_u(:,pint) = kw; sigmaBest_u(:,pint) = sigma_0;
                do t=1, this%stateEstimationEpochs;
                    ! Ensemble Kalman update
                    call this%EnKF_update( psi_kp1, this%powerObservation+this%Pstd*this%p_steps(pint), kw, sigma_0, yaw, t)
                    ! Run forward model pass
                    call this%forward(kw, sigma_0, Phat, yaw)
                    ! Normalized
                    Phat = Phat / Phat_zeroYaw(1)
                    error(t) = sum(abs(Phat-this%powerObservation-this%Pstd*this%p_steps(pint))) / real(this%Nt,rkind)
                    if (error(t)<lowestError) then
                        lowestError=error(t); bestStep = t;
                        kwBest_u(:,pint) = kw; sigmaBest_u(:,pint) = sigma_0
                        this%Phat_fit_u(:,pint) = Phat;
                    end if
                end do
            end do
        else
            kwBest_u = this%kw_u(this%indSorted,:); 
            sigmaBest_u = this%sigma_u(this%indSorted,:)
        end if
        this%Phat_fit = this%Phat_fit_u(this%unsort,(this%p_bins+1)/2)

        ! Rotate to future wind direction
        kwBest_u = kwBest_u(this%unsort,:); ! unsort
        sigmaBest_u = sigmaBest_u(this%unsort,:) ! unsort
        this%wind_direction = this%wind_direction_next ! set next step wind direction forecast
        this%turbCenter = this%turbCenterStore ! restore original layout
        call this%rotate() ! rotate
        kwBest_u = kwBest_u(this%indSorted,:); ! sort
        sigmaBest_u = sigmaBest_u(this%indSorted,:); ! sort

        ! Generate optimal yaw angles
        call this%yawOptimize_uncertain(kwBest_u, sigmaBest_u, yaw)
        ! Terms are already unsorted in the uncertain case due to local array
        ! rotations (wind direction uncertainty)
        this%kw_u = kwBest_u(this%unsort,:); this%sigma_u = sigmaBest_u(this%unsort,:)
        this%yaw = yaw(this%unsort)
        this%Phat = this%Phat(this%unsort); this%Phat_yaw = this%Phat_yaw(this%unsort)
        this%kw = kwBest_u(this%unsort,(this%p_bins+1)/2); this%sigma_0 = sigmaBest_u(this%unsort,(this%p_bins+1)/2) 
        ! Deallocate
        deallocate(kwBest_u)
        deallocate(sigmaBest_u)
    end if

    ! Deallocate
    deallocate(kw)
    deallocate(kwBest)
    deallocate(sigmaBest)
    deallocate(sigma_0)
    deallocate(Phat)
    deallocate(Phat_zeroYaw)
    deallocate(yaw)
    deallocate(psi_kp1)
 
end subroutine

subroutine EnKF_update(this, psi_k, P_kp1, kw, sigma_0, yaw, step)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(:,:), intent(inout) :: psi_k
    real(rkind), dimension(:), intent(in) :: P_kp1, yaw
    real(rkind), dimension(:), intent(inout) :: kw, sigma_0
    integer :: NN, i, j
    real(rkind), dimension(:,:), allocatable :: randArr, randArr2, chi, rbuff
    real(rkind), dimension(:,:), allocatable :: psi_kp, psi_kp1, psiHat_kp, ones_Ne, psi_kpp, psiHat_kpp
    real(rkind), dimension(:), allocatable :: Phat, Phat_zeroYaw
    integer, intent(in) :: step

    ! Define relevant parameters
    NN = this%Nt * 2
    allocate(psi_kp(size(psi_k,1), size(psi_k,2)))
    allocate(psi_kp1(size(psi_kp,1), size(psi_kp,2)))
    allocate(psi_kpp(size(psi_kp,1), size(psi_kp,2)))
    allocate(chi(size(psi_kp,1), size(psi_kp,2)))
    allocate(psiHat_kp(this%Nt, size(psi_kp,2)))
    allocate(psiHat_kpp(this%Nt, size(psi_kp,2)))
    allocate(Phat(this%Nt))
    allocate(Phat_zeroYaw(this%Nt))
    allocate(ones_Ne(this%Ne, this%Ne))
    allocate(rbuff(size(psi_kp,1), size(psi_kp,2)))
    allocate(randArr(this%Nt, this%Ne))
    allocate(randArr2(this%Nt, this%Ne))
 
    ! Random noise
    call gaussian_random(randArr,zero,this%var_k**0.5,this%ts+step)
    call gaussian_random(randArr2,zero,this%var_sig**0.5,this%ts*5+step*5)
    chi = 0.d0
    chi(1:this%Nt, :) = randArr
    chi(this%Nt+1:NN, :) = randArr2
    psi_kp = psi_k + chi
    randArr = psi_kp(1:this%Nt,:)
    randArr2 = psi_kp(this%Nt+1:NN, :)

    ! Intermediate forcast step
    call this%forward(kw, sigma_0, Phat_zeroYaw, yaw*0.d0)
    do i=1,this%Ne
        ! Lifting line model
        call this%forward(randArr(:,i), randArr2(:,i), Phat, yaw)
        ! Normalize by first turbine
        psiHat_kp(:,i) = Phat/Phat_zeroYaw(1)
        !psiHat_kp(:,i) = Phat/Phat(1)
    end do


    ! Perturbation
    ones_Ne = 1.d0 / real(this%Ne,rkind)
    rbuff = matmul(psi_kp, ones_Ne)
    psi_kpp = psi_kp - rbuff
    rbuff = matmul(psiHat_kp, ones_Ne)
    psiHat_kpp = psiHat_kp - rbuff
    
    ! Perturbation matrices
    randArr = 0.d0; randArr2 = 0.d0
    call gaussian_random(randArr,zero,this%var_p**0.5,this%ts*7+7*step)
    do i = 1, this%Ne
        randArr2(:,i) = P_kp1
    end do
    randArr2 = randArr2 + randArr
    
    
    ! Measurement analysis step
    psi_kp1 = psi_kp + matmul(matmul(matmul(psi_kpp, transpose(psiHat_kpp)), & 
            inv(matmul(psiHat_kpp,transpose(psiHat_kpp)) + &
            matmul(randArr, transpose(randArr)))), (randArr2 - psiHat_kp))

    ! Ouput the final values
    psi_kp1 = matmul(psi_kp1,ones_Ne)
    kw = psi_kp1(1:this%Nt,1);
    sigma_0 = psi_kp1(this%Nt+1:NN,1);
    !error = p_sum(abs(rbuff(:,1)-P_kp1)) / real(this%Nt)
    psi_k = psi_kp1
    
    ! Deallocate
    deallocate(psi_kp,psi_kp1,psi_kpp,chi,psiHat_kp,psiHat_kpp,Phat,ones_Ne)
    deallocate(rbuff,randArr,randArr2, Phat_zeroYaw) 

end subroutine

subroutine yawOptimize_uncertain(this, kw_u, sigma_u, yaw)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(:,:), intent(in) :: kw_u, sigma_u
    real(rkind), dimension(:), intent(inout) :: yaw
    real(rkind), dimension(this%Nt) :: P, P_baseline, bestYaw, Ptotal, Ptotaly, Pn
    real(rkind), dimension(this%Nt) :: m, v, yaw_orig, kw, sigma_0
    real(rkind), dimension(this%epochsYaw) :: Ptot
    real(rkind) :: Ptot_baseline, bestPower, dir_orig, ws_orig
    real(rkind), dimension(this%Nt,2) :: layout_orig, turbCenter
    integer :: k, pint, i, j, mi
    logical :: check
    real(rkind), dimension(this%epochsYaw, this%Nt) :: P_time, yawTime
    real(rkind), dimension(this%Nt) :: bestPowerOut, yaw_local, grads_total
    real(rkind), dimension(this%ws_bins) :: ws_bin, rho_ws
    real(rkind), dimension(this%dir_bins) :: dir_bin, rho_dir
    real(rkind), dimension(this%gamma_bins) :: gamma_bin, rho_gamma
    real(rkind), dimension(this%p_bins) :: rho_params
    integer, dimension(this%Nt) :: indSorted, unsort


    ! Store original information
    dir_orig = this%wind_direction
    layout_orig = this%turbCenter    
    ws_orig = this%wind_speed
    yaw_orig = yaw
    Ptotal = 0.d0; Ptotaly = 0.d0
    indSorted = this%indSorted; unsort = this%unsort; turbCenter = layout_orig

    ! Define the wind condition ranges
    if (this%ws_bins>1) then
        ws_bin = linspace(-this%ws_std, this%ws_std, this%ws_bins)
    else
        ws_bin = 0.d0
    end if
    if (this%dir_bins>1) then
        dir_bin = linspace(-this%dir_std, this%dir_std, this%dir_bins)
    else
        dir_bin = 0.d0
    end if
    if (this%gamma_bins>1) then
        gamma_bin = linspace(-this%gamma_std, this%gamma_std, this%gamma_bins) * pi / 180.d0
    else
        gamma_bin = 0.d0
    end if

    ! Define the probability distributions (for now, uniform distributions)
    rho_ws = 1.d0; rho_ws = rho_ws / sum(rho_ws)
    rho_dir = 1.d0; rho_dir = rho_dir / sum(rho_dir)
    rho_gamma = 1.d0; rho_gamma = rho_gamma / sum(rho_gamma)
    rho_params = 1.d0; rho_params = rho_params / sum(rho_params)

    ! Model
    ! Baseline power
    do pint = 1, this%p_bins
        do j = 1, this%dir_bins
            ! Rotate the turbine array
            turbCenter = layout_orig
            call this%rotate_local(dir_bin(j), turbCenter, indSorted, unsort)
            this%turbCenter = turbCenter
            do i = 1, this%ws_bins
                this%wind_speed = ws_orig + ws_bin(i)
                do mi = 1, this%gamma_bins
                    yaw = yaw_orig(indSorted)
                    yaw(1) = yaw(1) + gamma_bin(mi)
                    ! Run forward pass
                    call this%forward(kw_u(indSorted,pint), sigma_u(indSorted,pint), P, yaw*0.d0)
                    Ptot_baseline = sum(P); P_baseline = P; this%Phat = P_baseline;
                    call this%forward(kw_u(indSorted,pint), sigma_u(indSorted,pint), this%Phat_yaw, yaw)
                    ! Update total power
                    Ptotaly = Ptotaly + this%Phat_yaw(unsort)/this%Phat(1) * rho_ws(i) * rho_dir(j) * & 
                             rho_params(pint) * rho_gamma(mi);
                    Ptotal = Ptotal + this%Phat(unsort)/this%Phat(1) * rho_ws(i) * rho_dir(j) * & 
                             rho_params(pint) * rho_gamma(mi);
                end do
            end do
        end do
    end do
    Ptot_baseline = sum(Ptotal); P_baseline = Ptotal; 
    ! Output model forward passes
    this%Phat = Ptotal
    this%Phat_yaw = Ptotaly
    this%Ptot_initial = Ptot_baseline

    ! Initialize the yaw to zero for re-optimization
    yaw = 0.d0
 
    ! eps also determines the termination condition
    k=1; check = .false.; 
    bestYaw = 0.d0; Ptot = 0.d0; P_time = 0.d0; P_time = 0.d0; yawTime = 0.d0
    m=0.d0; v=0.d0;
    bestPower = 0.d0; bestYaw = 0.d0; bestPowerOut = 0.d0;
    do while (k < this%epochsYaw .and. (.not. check)) 
    
        Ptotal = 0.d0; grads_total = 0.d0
        do pint = 1, this%p_bins
            do j = 1, this%dir_bins
                ! Rotate the turbine array
                turbCenter = layout_orig
                call this%rotate_local(dir_bin(j), turbCenter, indSorted, unsort)
                this%turbCenter = turbCenter
                do i = 1, this%ws_bins
                    this%wind_speed = ws_orig + ws_bin(i)
                    do mi = 1, this%gamma_bins
                        ! Forward prop
                        yaw_local = yaw(indSorted)
                        yaw_local(1) = yaw_local(1) + gamma_bin(mi)
                        this%dp_dgamma = 0.d0
                        call this%forward(kw_u(indSorted,pint), sigma_u(indSorted,pint), Pn, yaw_local*0.d0)
                        call this%forward(kw_u(indSorted,pint), sigma_u(indSorted,pint), P, yaw_local)
                        yawTime(k+1,:) = yaw;
                        ! Update total power
                        Ptotal = Ptotal + (P(unsort)/Pn(1)) * rho_ws(i) * rho_dir(j) * & 
                                 rho_params(pint) * rho_gamma(mi);
        
                        ! Total gradients
                        call this%backward(kw_u(indSorted,pint), sigma_u(indSorted,pint), yaw_local)    
                        grads_total = grads_total + this%dp_dgamma(unsort) * &
                            rho_ws(i) * rho_dir(j) * &
                            rho_params(pint) * rho_gamma(mi);                     
                        end do
                    end do
                end do
            end do
            Ptot(k) = sum(Ptotal);
            P_time(k,:) = Ptotal;

            ! Gradient ascent update yaw
            ! Adam
            m = this%beta1*m + (1.d0-this%beta1)*grads_total;
            v = this%beta2*v + (1.d0-this%beta2)*(grads_total**2);
            yaw = yaw + this%learning_rate_yaw * m / & 
                      (sqrt(v) + this%eps)
            if (Ptot(k) > bestPower) then
                bestPower = Ptot(k)
                bestPowerOut = Ptotal
                bestYaw = yaw
            end if
            k = k+1
            if (k > 10) then
                if (abs(Ptot(k-1)-Ptot(k-2))/abs(Ptot(k-1)) < 10D-9) then
                    check = .true.
                end if
            end if
    end do
    
    ! Final results
    yaw = bestYaw
    this%Ptot_final = bestPower
    this%Popti = bestPowerOut

    ! Restore 
    this%turbCenter = layout_orig
    this%wind_speed = ws_orig
    this%wind_direction = dir_orig

end subroutine


subroutine yawOptimize(this, kw, sigma_0, yaw)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(:), intent(in) :: kw, sigma_0
    real(rkind), dimension(:), intent(inout) :: yaw
    real(rkind), dimension(this%Nt) :: P, P_baseline, bestYaw
    real(rkind), dimension(this%Nt) :: m, v
    real(rkind), dimension(this%epochsYaw) :: Ptot
    real(rkind) :: Ptot_baseline, bestPower
    integer :: k
    logical :: check
    real(rkind), dimension(this%epochsYaw, this%Nt) :: P_time, yawTime
    real(rkind), dimension(this%Nt) :: bestPowerOut

    ! Model
    ! Load model inputs
    call this%forward(kw, sigma_0, P, yaw*0.d0)
    Ptot_baseline = sum(P); P_baseline = P; this%Phat = P_baseline;
    call this%forward(kw, sigma_0, this%Phat_yaw, yaw)
    this%Phat_yaw = this%Phat_yaw / this%Phat(1)
    this%Phat = this%Phat / this%Phat(1)

    ! Initialize the yaw to zero for re-optimization
    yaw = 0.d0
 
    ! eps also determines the termination condition
    k=1; check = .false.; 
    bestYaw = 0.d0; Ptot = 0.d0; P_time = 0.d0; P_time = 0.d0; yawTime = 0.d0
    m=0.d0; v=0.d0;
    bestPower = 0.d0; bestYaw = 0.d0; bestPowerOut = 0.d0;
    do while (k < this%epochsYaw .and. (.not. check)) 
    
        ! Forward prop
        this%dp_dgamma = 0.d0
        call this%forward(kw, sigma_0, P, yaw)
        yawTime(k+1,:) = yaw;
        Ptot(k) = sum(P);
        P_time(k,:) = P;
        call this%backward(kw, sigma_0, yaw)    
    
        ! Gradient ascent update yaw
        ! Adam
        m = this%beta1*m + (1.d0-this%beta1)*this%dp_dgamma;
        v = this%beta2*v + (1.d0-this%beta2)*(this%dp_dgamma**2);
        yaw = yaw + this%learning_rate_yaw * m / & 
                   (sqrt(v) + this%eps)
        if (Ptot(k) > bestPower) then
           bestPower = Ptot(k)
           bestPowerOut = P
           bestYaw = yaw
        end if
        k = k+1
        if (k > 10) then
            if (abs(Ptot(k-1)-Ptot(k-2))/abs(Ptot(k-1)) < 10D-9) then
                check = .true.
            end if
        end if
    end do
    
    ! Final results
    yaw = bestYaw
    this%Ptot_final = bestPower / P_baseline(1)
    this%Popti = bestPowerOut / P_baseline(1)

end subroutine

subroutine forward(this, kw, sigma_0, Phat, yaw)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(:), intent(in) :: kw, sigma_0, yaw
    real(rkind), dimension(this%Nt), intent(out) :: Phat
    real(rkind) :: D, R, dx, boundLow, boundHigh, edgeLow, edgeHigh, uc_error, dvi, dv
    integer :: Nx, i, j, k, ucCount
    real(rkind), dimension(this%Nt) :: delta_v0, Uc, Uc_old
    real(rkind) :: dw, aPrev, du, gaussian, delta_u_face, L, temp
    logical, dimension(this%Nt, this%Nt) :: turbinesInFront
    real(rkind), dimension(this%Nx) :: xpFront, delta_v, dwVect, dvec
    logical :: check
    real(rkind), dimension(this%Ny) :: duSum, us, ylocal, uw
    real(rkind), dimension(this%Nt, this%Nt) :: uci
 
    ! Definitions
    ! Set the relevant turbine length scale
    L = 1.d0
    turbinesInFront = .false.
    this%A = (this%D*L) ** 2. * pi / 4.d0
    D = this%D / this%D
    R = D/2.d0
    ! Atmospheric conditions
    this%rho = 1.d0
    ! Convective superposition
    Uc = this%wind_speed; uci = this%wind_speed
    Uc_old = this%wind_speed

    !! Forward prop
    ! While loop for convective velocity convergences
    uc_error = 1.0D10; ucCount = 1;
    do while (uc_error>this%epsUc .and. ucCount<this%ucMaxIt)
        ! Initialize
        this%v = 0.d0
        ! Loop over turbines
        do j = 1,this%Nt
            ! Find the leading turbine
            do i = 1,this%Nt
                dx = this%turbCenter(j, 1) - this%turbCenter(i, 1);
                if (dx > 0.d0) then
                   dw = 1.d0 + kw(i) * log(1.d0 + exp((dx - 2.0d0 * R) / R));
                   boundLow = this%turbCenter(i,2)-dw/2;
                   boundHigh = this%turbCenter(i,2)+dw/2;
                   edgeLow = this%turbCenter(j,2)-D/2;
                   edgeHigh = this%turbCenter(j,2)+D/2;
                   check = .false.
                   if (edgeLow>=boundLow .and. edgeLow<=boundHigh) then
                       check = .true.
                   end if
                   if (edgeHigh>=boundLow .and. edgeHigh<=boundHigh) then
                       check = .true.
                   end if
                   if (check) then
                       turbinesInFront(j,i) = .true.
                   end if
                end if
            end do
            ! Loop over turbines and superposition
            delta_u_face = 0.d0
            duSum = 0.d0; us = 0.d0; duSum = 0.d0; uw = 0.d0
            ylocal = linspace(this%turbCenter(j, 2)-D/2, this%turbCenter(j, 2)+D/2, this%Ny);
            do k = 1, this%Nt
                if (turbinesInFront(j,k)) then
                    ! Individual velocity deficits
                    aPrev = 0.5*(1.d0-sqrt(1.d0-this%Ct*cos(yaw(k))**2))
                    dx = this%turbCenter(j, 1) - this%turbCenter(k, 1)
                    xpFront = linspace(0.d0, dx, this%Nx)
                    dw = 1.d0 + kw(k) * log(1.d0 + exp((dx - 2.0 * R) / R));
                    if (this%superposition == 1) then
                        du = this%wind_speed * aPrev/(dw**2+this%eps) * & 
                                (1.d0+erf(dx/(R*sqrt(2.d0))));
                    elseif (this%superposition == 2) then
                        du = this%u_eff(k) * aPrev/(dw**2+this%eps) * & 
                                (1.d0+erf(dx/(R*sqrt(2.d0))));
                    elseif (this%superposition == 3) then
                        du = this%u_eff(k) * aPrev/(dw**2+this%eps) * & 
                                (1.d0+erf(dx/(R*sqrt(2.d0))));
                    end if
                    this%deltaUIndividual(j,k) = du

                    ! Deflection
                    dwVect = 1.d0 + kw(k) * log(1.d0 + exp((xpFront - 2.0 * R) / R));
                    if (this%superposition == 1) then
                        dvi = this%wind_speed * 0.25 * this%Ct * cos(yaw(k))**2 * sin(yaw(k))
                    elseif (this%superposition == 2 .or. this%superposition == 3) then
                        dvi = this%u_eff(k) * 0.25 * this%Ct * cos(yaw(k))**2 * sin(yaw(k))
                    end if
                    if (this%secondary) then
                        delta_v0(j) = dvi + this%v(k)
                    else
                        delta_v0(j) = dvi
                    end if
                    ! Velocity deficit
                    delta_v = (delta_v0(j) / (dwVect**2+this%eps)) * 0.5 * (1.d0 + erf(xpFront & 
                              / (R*sqrt(2.d0))))
                    where (xpFront<0.d0)
                        delta_v = 0.d0
                    end where
                    dv = (1.d0 / (dw**2+this%eps)) * 0.5 * (1.d0 + erf(dx & 
                              / (R*sqrt(2.d0))))
                    dvec = (1.d0 / (dwVect**2+this%eps)) * 0.5 * (1.d0 + erf(xpFront & 
                              / (R*sqrt(2.d0))))
                    
                    ! y_c
                    if (this%superposition == 1) then
                        this%y_c(k, j) = mytrapz(xpFront, -delta_v/this%wind_speed)
                    elseif (this%superposition == 2 .or. this%superposition == 3) then
                        this%y_c(k, j) = mytrapz(xpFront, -delta_v/this%u_eff(k))
                    end if
                    ! Local y frame
                    this%sigmaEnd(k, j) = sigma_0(k) * dw
                    ! take the last value of sigma which is at the next turbine
                    this%yCenter(k, j) = this%y_c(k, j) + this%turbCenter(k, 2);
                    ! Effective velocity at turbine face
                    gaussian = erf((this%turbCenter(j, 2)+D/2.d0-this%yCenter(k,j))/ & 
                                 sqrt(2.d0*this%sigmaEnd(k, j)**2+this%eps)) - &
                                 erf((this%turbCenter(j, 2)-D/2.d0-this%yCenter(k, j))/ & 
                                 sqrt(2.d0*this%sigmaEnd(k,j)**2+this%eps));
                    gaussian = gaussian * (D**2)/(16.d0*sigma_0(k)**2) & 
                               * this%sigmaEnd(k,j)*sqrt(2.d0*pi)
                    if (this%superposition == 1) then
                        delta_u_face = delta_u_face + du * gaussian;
                    elseif (this%superposition == 2 .or. this%superposition == 3) then
                        ! compute local convective velocity
                        us = du * D**2 / (8*sigma_0(k)**2) * &
                             exp(-(ylocal-this%yCenter(k,j))**2 / (2*this%sigmaEnd(k, j)**2));
                        uw = this%u_eff(k) - us
                        uci(k, j) = mytrapz(ylocal, uw*us) / mytrapz(ylocal, us);
                        if (this%superposition==3) then
                            uci(k, j) = this%wind_speed
                        end if
                        ! Update cache stuff
                        delta_u_face = delta_u_face + du * gaussian * uci(k, j) / Uc(j);
                        this%deltaUIndividual(j,k) = du * uci(k, j) / Uc(j);
                        duSum = duSum + us * uci(k, j) / Uc(j);
                        this%erfStore(k, j) = (du/this%u_eff(k)) * gaussian * uci(k, j) / Uc(j); 
                    end if
                    this%gaussianStore(k, j) = gaussian;
                    if (this%secondary) then
                        this%ucr(k, j) = (uci(k, j) / Uc(j)) * dv * gaussian;
                        this%v(j) = this%v(j) + dvi * this%ucr(k, j);
                    end if
                end if
            end do
            ! Superposition
            if (this%superposition == 2 .or. this%superposition == 3) then
                ! Global convection velocity
                Uc_old(j) = Uc(j); us = duSum;
                uw = this%wind_speed - us;
                temp = mytrapz(ylocal, us)
                if (temp < this%eps) then
                    Uc(j) = this%wind_speed
                else
                    Uc(j) = mytrapz(ylocal, uw*us) / mytrapz(ylocal, us)
                end if
                if (this%superposition==3) then
                    Uc(j) = this%wind_speed
                end if
            end if
            this%delta_u_face_store(j) = delta_u_face;
            this%u_eff(j) = (1.d0/D) * (D*this%wind_speed - delta_u_face);
            ! Calculate the induction factor at the end
            this%a_mom(j) = 0.5*(1.d0-sqrt(1.d0-this%Ct*cos(yaw(j))**2))
            ! Power
            this%ap(j) = 0.5d0*(1.d0 - sqrt(1.d0-this%Ct))
            this%Cp(j) = 4.d0*this%eta*this%ap(j)*(1.d0-this%ap(j))**2 * & 
                           cos(yaw(j))**this%powerExp
            Phat(j) = 0.5 * this%rho * this%A * this%Cp(j) * this%u_eff(j)**3;  
        end do
        ! Error update
        if (this%superposition == 2 .or. this%superposition == 3) then
            uc_error = sum(abs(Uc_old-Uc))/real(this%Nt); 
            ucCount = ucCount+1; this%turbinesInFront = turbinesInFront
        else
            uc_error=0.d0;
        end if
    end do

end subroutine


subroutine backward(this, kw, sigma_0, yaw)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(:), intent(in) :: kw, sigma_0, yaw
    integer :: i, j, k, m
    real(rkind), dimension(this%Nt) :: dp_du_eff, dp_dcp, dcp_dap
    real(rkind), dimension(this%Nt) :: du_eff_dyc, dy_c_dgamma
    real(rkind), dimension(this%Nt) :: dcp_dgamma, da_dgamma
    real(rkind) :: dx, R=0.5, D = 1.d0, du_eff_dueffUp
    logical :: check
    real(rkind) :: boundLow, boundHigh, edgeLow, edgeHigh, dw, du_eff_da
    logical, dimension(this%Nt) :: turbinesInBack
    real(rkind), dimension(this%Nx) :: dwVect, xpFront
    logical, dimension(this%Nt, this%Nt) :: turbinesInBackMat

    !! Backprop
    this%dp_dgamma = 0.d0
    ! Loop over turbines
    do i = this%Nt, 1, -1
        ! Same for all turbines
        ! d a / d gamma
        da_dgamma(i) = -0.5*this%Ct*sin(yaw(i))*cos(yaw(i))/ & 
                       sqrt(1.d0-this%Ct*cos(yaw(i))**2 + this%eps)
        ! dP / d u_eff
        dp_du_eff(i) = 0.5*this%rho*this%A*this%Cp(i)*this%u_eff(i)**2 * 3.d0;
        ! dP / dcp
        dp_dcp(i) = 0.5*this%rho*this%A*this%u_eff(i)**3
        ! d cp / da
        dcp_dap(i) = 4.d0*this%eta*(1.d0 - 4.d0*this%ap(i) + &
                     3.d0*this%ap(i)**2) * cos(yaw(i))**2
        dcp_dgamma(i) = -4.d0*this%eta * this%ap(i) * (1.d0-this%ap(i))**2*sin(yaw(i)) * &
                         cos(yaw(i))**(this%powerExp-1.d0) * this%powerExp
        this%dp_dgamma(i) = this%dp_dgamma(i) + dp_dcp(i)*dcp_dgamma(i);
    end do
    do i = this%Nt, 1, -1
        turbinesInBack = .false.; turbinesInBackMat = .false.
        do j = 1, this%Nt
            dx = this%turbCenter(j, 1) - this%turbCenter(i, 1);
            if (dx > 0) then
               dw = 1.d0 + kw(i) * log(1.d0 + exp((dx - 2.0 * R) / R));
               boundLow = this%turbCenter(i,2)-dw/2.d0;
               boundHigh = this%turbCenter(i,2)+dw/2.d0;
               edgeLow = this%turbCenter(j,2)-D/2.d0;
               edgeHigh = this%turbCenter(j,2)+D/2.d0;
               check = .false.
               if (edgeLow>=boundLow .and. edgeLow<=boundHigh) then
                   check = .true.
               end if
               if (edgeHigh>=boundLow .and. edgeHigh<=boundHigh) then
                   check = .true.
               end if
               if (check) then
                   turbinesInBack(j) = .true.
                   turbinesInBackMat(i,j) = .true.
               end if
            end if
        end do
        ! dP / dgamma
        ! Flip over whether the turbine sees a wake or not
        du_eff_da = 0;
        do k = 1, this%Nt
            if (turbinesInBack(k)) then
                du_eff_da = -this%deltaUIndividual(k,i) / (this%a_mom(i)+this%eps) * & 
                            this%gaussianStore(i,k)
                this%dp_dgamma(i) = this%dp_dgamma(i) + & 
                                    dp_du_eff(k)*du_eff_da*da_dgamma(i)
                ! d u_eff / d y_c
                du_eff_dyc(i) = -(1.d0/D) * (this%delta_u_face_store(k) * & 
                                D**2 / (8.d0*sigma_0(i)**2) * &
                                (-exp(-(this%turbCenter(k, 2)+D/2.d0-&
                                this%yCenter(i, k))**2/(2.d0*this%sigmaEnd(i,k)**2))+& 
                                exp(-(this%turbCenter(k, 2)-D/2.d0-&
                                this%yCenter(i, k))**2/(2.d0*this%sigmaEnd(i,k)**2))) )
                !!!
                ! d y_c / d gamma
                dx = this%turbCenter(k, 1) - this%turbCenter(i, 1)
                xpFront = linspace(0.d0, dx, this%Nx)
                dw = 1.d0 + kw(i) * log(1.d0 + exp((dx - 2.0d0 * R) / R))
                dwVect = 1.d0 + kw(i) * log(1.d0 + exp((xpFront - 2.0d0 * R) / R))
                ! Changed this when using the local velocity to compute y_c
                dy_c_dgamma(i) = (cos(yaw(i))**3 - 2.d0*sin(yaw(i))**2 * & 
                                 cos(yaw(i))) * &
                                 mytrapz(xpFront, -0.25d0*this%Ct*0.5d0*(1.d0+erf(xpFront/ & 
                                 (R*sqrt(2.d0))))/dwVect**2)
                ! Sum together to get dP / dgamma
                this%dp_dgamma(i) = this%dp_dgamma(i) + & 
                                    dp_du_eff(k)*du_eff_dyc(i)*dy_c_dgamma(i) 
                ! Extra terms from secondary steering
                if ((this%superposition==2 .or. this%superposition==3) .and. this%secondary) then;
                    do m = 1, this%Nt
                        if (this%turbinesInFront(i, m)) then
                            ! d y_c / d gamma
                            dy_c_dgamma(m) = (this%u_eff(m)/this%u_eff(i)) * this%ucr(m, i) &
                                             * (cos(yaw(m))**3-2*sin(yaw(m))**2*cos(yaw(m)))* &
                                             mytrapz(xpFront, -0.25*this%Ct*0.5*(1.d0+erf(xpFront &
                                             / (R*sqrt(2.d0)))) / dwVect**2);
                            ! Sum together to get dP / dgamma
                            this%dp_dgamma(m) = this%dp_dgamma(m) + &
                                           dp_du_eff(k)*du_eff_dyc(i)*dy_c_dgamma(m);
                        end if 
                    end do
                end if
            end if
        end do    
    end do
    ! Extra terms from convective superposition
    if (this%superposition == 2 .or. this%superposition == 3) then
        do i = this%Nt, 1, -1
            do m = 1, this%Nt
            if (this%turbinesInFront(i, m)) then
                do k = 1, this%Nt
                if (turbinesInBackMat(i, k)) then
                    du_eff_dueffUp = -this%erfStore(i, k);
                    ! 1
                    du_eff_da = -this%deltaUIndividual(i,m) / (this%a_mom(m)+this%eps) * this%gaussianStore(m, i);
                    this%dp_dgamma(m) = this%dp_dgamma(m) + &
                        dp_du_eff(k)*du_eff_dueffUp*du_eff_da*da_dgamma(m); 
                    ! 2
                    ! d u_eff / d y_c
                    du_eff_dyc(m) = -(1/D) * (this%delta_u_face_store(i) * & 
                                    D**2 / (8*sigma_0(m)**2) * &
                                    (-exp(-(this%turbCenter(i, 2)+D/2 - &
                                    this%yCenter(m,i))**2/(2*this%sigmaEnd(m, i)**2)) + &
                                    exp(-(this%turbCenter(i, 2)-D/2 - &
                                    this%yCenter(m,i))**2/(2*this%sigmaEnd(m, i)**2))));
                    ! d y_c / d gamma
                    dx = this%turbCenter(i, 1) - this%turbCenter(m, 1)
                    xpFront = linspace(0.d0, dx, this%Nx)
                    dwVect = 1.d0 + kw(m) * log(1.d0 + exp((xpFront - 2.0d0 * R) / R))
                    dy_c_dgamma(m) = (cos(yaw(m))**3-2*sin(yaw(m))**2*cos(yaw(m))) * & 
                                     mytrapz(xpFront, -0.25*this%Ct*0.5 * &
                                     (1.d0+erf(xpFront / (R*sqrt(2.d0))))/dwVect**2);
                    this%dp_dgamma(m) = this%dp_dgamma(m) + &
                                   dp_du_eff(k)*du_eff_dueffUp*du_eff_dyc(m)*dy_c_dgamma(m); 
                end if
                end do
            end if
            end do
        end do
    end if

 
end subroutine


!! alpha_check() is used to check the stationarity of the wind direction and
!! to produce the statistics of the wind direction relevant for wake
!! steering
subroutine alpha_check(this, alpha, t, alpha_m, alpha_std, Tf_out)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(:), intent(in) :: t, alpha ! should be time and alpha from t-2*Tf_init:t where t is current timestep
    real(rkind), intent(out) :: alpha_m, alpha_std
    integer, intent(out) :: Tf_out
    real(rkind), dimension(:), allocatable :: Rsq, x, y, yp, sig 
    real(rkind), dimension(:,:), allocatable :: cov, p
    integer, dimension(:), allocatable :: p_in
    real(rkind) :: dalpha, eps_f, eps_m, ym, t_future, dt, Tint, SStot, SSres, Rsq_min = 0.2d0
    integer :: tstart, Tf
    logical :: check, stationary
    ! Allocate
    allocate(p_in(2))
    allocate(Rsq(2))
    allocate(p(2,2))
    allocate(cov(2,2))
    allocate(x(this%Tf_init))
    allocate(y(this%Tf_init))
    allocate(yp(this%Tf_init))
    allocate(sig(this%Tf_init))

    ! Initailize
    tstart=size(t)
    check=.true.; dalpha=this%dalpha_max+1.d0; Tf=this%Tf_init;
    p_in = (/ 0, 1 /)
    dt = mean(t(2:2*this%Tf_init)-t(1:2*this%Tf_init-1)) 
    sig = 1.d0

    ! While loop over Tfit integration length
    do while (dalpha>=this%dalpha_max .and. check)
        ! Initialize loop
        ! Least squares fit, Step 1
        ! t-2*Tf:t-Tf
        x = t(tstart - 2*this%Tf_init+1 : tstart - 1*this%Tf_init); 
        y = alpha(tstart - 2*this%Tf_init+1 : tstart - 1*this%Tf_init); ym = mean(y);
        call pfit(x,y,sig,p_in,p(:,1),cov)
        yp = p(2,1)*x + p(1,1);
        SStot = sum((y - mean(y))**2.d0);
        SSres = sum((y - yp)**2.d0);
        Rsq(1) = 1.d0 - SSres/SStot;
        ! t-Tf:t
        x = t(tstart - 1*this%Tf_init+1 : tstart); 
        y = alpha(tstart - 1*this%Tf_init+1 : tstart); !ym(2) = mean(y);
        call pfit(x,y,sig,p_in,p(:,2),cov)
        yp = p(2,2)*x + p(1,2);
        SStot = sum((y - mean(y))**2.d0);
        SSres = sum((y - yp)**2.d0);
        Rsq(2) = 1.d0 - SSres/SStot;

        ! Step 2
        stationary = .true.; dalpha = 0.d0; x=0.d0; y = 0.d0;
        if (Rsq(1)>Rsq_min .and. Rsq(2)>Rsq_min) then
            ! Step 3
            y(1:Tf) = p(2,1)*t(tstart - 1*Tf+1 : tstart) + p(1,1);
            eps_f = mean( (alpha(tstart - 1*Tf+1 : tstart) - y(1:Tf) )**2.d0 );
            ! Step 4
            eps_m = mean( (alpha(tstart - 1*Tf+1 : tstart) - ym)**2.d0 );
            ! Step 5
            if (eps_f < eps_m) then
                stationary = .false.;
                ! Steps 6 & 7
                y(1:Tf) = p(2,1)*t(tstart - 1*Tf+1 : tstart) + p(1,1);
                alpha_std = std( alpha(tstart - 1*Tf+1 : tstart) - y(1:Tf) );
                ! Add previous predictive error to the STD
                alpha_std = alpha_std + sqrt(eps_f);
                ! Time
                t_future = t(tstart - 1*Tf+1) + real(Tf/2,rkind) * dt
                Tint = t(tstart) - t(tstart - 1*Tf+1)
                ! Predicted mean
                alpha_m = p(2,2)*t_future + p(1,2);
                dalpha = Tint*abs(p(2,2));
            else 
                stationary = .true.;
                ! Steps 9 & 10
                alpha_std = std( alpha(tstart - 1*Tf+1 : tstart) );
                ! Add previous predictive error to the STD
                alpha_std = alpha_std + sqrt(eps_m);
                ! Mean
                alpha_m = mean( alpha(tstart - 1*Tf+1 : tstart) );
                dalpha = abs(alpha_m - mean( alpha(tstart - 2*Tf+1 : tstart - 1*Tf) ));
            end if
        else ! step 11
            stationary = .true.;
            eps_m = mean( (alpha(tstart - 1*Tf+1 : tstart) - ym)**2.d0 );
            ! Steps 12 & 13
            alpha_std = std( alpha(tstart - 1*Tf+1 : tstart) );
            ! Add previous predictive error to the STD
            alpha_std = alpha_std + sqrt(eps_m);
            ! Mean
            alpha_m = mean( alpha(tstart - 1*Tf+1 : tstart) );
            dalpha = abs(alpha_m - mean( alpha(tstart - 2*Tf+1 : tstart - 1*Tf) ));
        end if

        ! Update Tf
        if (dalpha>this%dalpha_max) then;
            if (Tf/2>=this%Tf_min) then;
                Tf = Tf/2; 
            else
                check=.false.;
            end if
        end if
        Tf_out = Tf;

    end do
    write(*,*) 'Stationary'
    write(*,*) stationary

    ! Deallocate
    deallocate(p_in,Rsq,p,cov,x,y,sig,yp)

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Utilities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine observeField(this)
    class(dynamicYaw), intent(inout) :: this
    this%wind_speed = 8.d0 ! Get this from the data
    this%wind_direction = 270.d0 ! Get this from the data
    this%powerObservation = (/1.d0, 0.6d0/) ! Get the power production from ADM code 

end subroutine

subroutine rotate(this)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(2,2) ::  R
    real(rkind), dimension(this%Nt,2) :: X
    integer, dimension(:), allocatable :: indSorted
    type(sortgroup), dimension(this%Nt) :: x_for_sorting
    integer :: i
    real(rkind) :: a
 
    a = (this%wind_direction - 270.d0) * pi / 180.d0
    R = reshape((/cos(a), sin(a), -sin(a), cos(a)/), shape(R))
    ! shift
    this%turbCenter(:,1) = this%turbCenter(:,1) - this%turbCenter(this%conditionTurb,1)
    this%turbCenter(:,2) = this%turbCenter(:,2) - this%turbCenter(this%conditionTurb,2)
    ! rotate
    do i=1,this%Nt
        X(i,:) = matmul(this%turbCenter(i,:), transpose(R))
        x_for_sorting(i)%value = X(i,1)
        x_for_sorting(i)%zpos = i
        this%unsort(i) = i
    end do

    ! Sort the turbines by upwind location
    call Qsort(x_for_sorting,this%Nt)
    this%indSorted = x_for_sorting%zpos
    this%unsort(indSorted) = this%unsort      
    this%turbCenter = X(this%indSorted,:)

end subroutine

subroutine rotate_local(this, dir, turbCenter, indSorted, unsort)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), intent(in) :: dir
    real(rkind), dimension(:,:), intent(inout) :: turbCenter
    integer, dimension(:), intent(inout) :: indSorted, unsort
    real(rkind), dimension(2,2) ::  R
    real(rkind), dimension(this%Nt,2) :: X
    type(sortgroup), dimension(this%Nt) :: x_for_sorting
    integer :: i
    real(rkind) :: a
    ! Rotation matrix 
    a = dir * pi / 180.d0
    R = reshape((/cos(a), sin(a), -sin(a), cos(a)/), shape(R))
    ! rotate
    do i=1,this%Nt
        X(i,:) = matmul(turbCenter(i,:), transpose(R))
        x_for_sorting(i)%value = X(i,1)
        x_for_sorting(i)%zpos = i
        unsort(i) = i
    end do
    ! Sort the turbines by upwind location
    call Qsort(x_for_sorting,this%Nt)
    indSorted = x_for_sorting%zpos
    unsort(indSorted) = unsort      
    turbCenter = X(indSorted,:)

end subroutine



! Average statistics
subroutine simpleMovingAverage(this, meanP, power, meanWs, ws, meanPbaseline, powerBaseline, meanDir, windDir, stdP, i, t)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), intent(inout) :: meanP, meanWs, meanPbaseline, meanDir, stdP
    real(rkind), intent(in) :: power, ws, powerBaseline, windDir
    integer, intent(in) :: i, t

    if (i>this%n_moving_average) then
        ! Power
        meanP = meanP + (1.d0/real(this%n_moving_average,rkind)) * (power - this%power_minus_n(1,t))
        this%power_minus_n(1:this%n_moving_average-1,t) = this%power_minus_n(2:this%n_moving_average,t)
        this%power_minus_n(this%n_moving_average,t) = power
        ! Wind speed
        meanWs = meanWs + (1.d0/real(this%n_moving_average,rkind)) * (ws - this%ws_minus_n(1,t))
        this%ws_minus_n(1:this%n_moving_average-1,t) = this%ws_minus_n(2:this%n_moving_average,t)
        this%ws_minus_n(this%n_moving_average,t) = ws
        ! Power baseline
        meanPbaseline = meanPbaseline + (1.d0/real(this%n_moving_average,rkind)) * (powerBaseline - this%pb_minus_n(1,t))
        this%pb_minus_n(1:this%n_moving_average-1,t) = this%pb_minus_n(2:this%n_moving_average,t)
        ! Wind direction
        meanDir = meanDir + (1.d0/real(this%n_moving_average,rkind)) * (windDir - this%dir_minus_n(1,t))
        this%dir_minus_n(1:this%n_moving_average-1,t) = this%dir_minus_n(2:this%n_moving_average,t)
    else
        ! Power
        this%power_minus_n(i,t) = power
        meanP = sum(this%power_minus_n(1:i,t)) / real(i,rkind)
        ! Wind speed
        this%ws_minus_n(i,t) = ws
        meanWs = sum(this%ws_minus_n(1:i,t)) / real(i,rkind)
        ! Power baseline
        this%pb_minus_n(i,t) = powerBaseline
        meanPbaseline = sum(this%pb_minus_n(1:i,t)) / real(i,rkind)
        ! Wind direction
        this%dir_minus_n(i,t) = windDir
        meanDir = sum(this%dir_minus_n(1:i,t)) / real(i,rkind)
        ! STD power
        stdP = sqrt( (1.d0/real(i,rkind)) * sum((this%power_minus_n(1:i,t) - meanP) ** 2.d0) )
    end if


end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LAPACK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function mean(x) result(xm)
    real(rkind), dimension(:), intent(in) :: x
    real(rkind) :: xm
    xm = sum(x) / real(size(x),rkind)
end function mean

function std(x) result(xstd)
    real(rkind), dimension(:), intent(in) :: x
    real(rkind) :: xm, xstd
    xm = sum(x) / real(size(x),rkind)
    xstd = sqrt( (1.d0/real(size(x),rkind)) * sum((x-xm)**2.d0) )
end function std

function inv(A) result(Ainv)
  real(rkind), dimension(:,:), intent(in) :: A
  real(rkind), dimension(size(A,1),size(A,2)) :: Ainv

  real(rkind), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv

  pure function integrate(x, y) result(r)
    !! Calculates the integral of an array y with respect to x using the
    !trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be
    !uniform.
    real(rkind), intent(in)  :: x(:)         !! Variable x
    real(rkind), intent(in)  :: y(size(x))   !! Function y(x)
    real(rkind)              :: r            !! Integral ∫y(x)·dx

    ! Integrate using the trapezoidal rule
    associate(n => size(x))
      r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
    end associate
  end function

subroutine pfit(x,y,sig,p,a,cov,coeff,chi)
    ! Fit data to a polynomial a_0 + a_1 x + ... + a_d x**d
    ! Original call: subroutine pfit(x,y,sig,p,a,cov,coeff,chi)
    ! Inputs:
    !   x(1:npt)           - abscissas
    !   y(1:npt)           - data values
    !   sig(1:npt)         - data errors
    ! Outputs:
    !   a(1:np)            - max. likelihood parameters
    !   coeff(1:npt,1:np)  - coefficients giving the max. likelihood parameters
    !                        in terms of the data:
    !                           a(i) = \Sum_{j} coeff(j,i) * y(j)
    !   cov(1:n,1:np)      - Covariance matrix, cov(i,j) = Cov(a(i),a(j))
    !                        The estimated error in a(i) is sqrt(Cov(a(i),a(i)))
    !   chi                - Reduced chi value,
    !                          chi = sqrt(chi**2/(npt - np))
    ! Notes:
    !   This routine uses a QR decomposition method, which should be more
    !   numerically stable than solving the normal equations.
    implicit none
    real(rkind), intent(in)  :: x(:), y(:), sig(:)
    integer,  intent(in)  :: p(:)
    real(rkind), intent(out) :: a(:)
    real(rkind), intent(out), optional :: cov(:,:), coeff(:,:), chi

    real(rkind), allocatable :: work(:), C(:,:), Q(:,:), R(:,:), b(:)
    integer :: ipiv(size(a)), lwork
    integer :: k,npt,ifail,np,ip,jp
    real(rkind) :: coeff1(size(x),size(a)), val

    npt = size(x) ! Number of data points
    np = size(p) ! Number of polynomial terms
    if (size(a) .ne. np) stop "Error 0 in pfit"
    if (size(y) .ne. npt) stop "Error 1 in pfit"
    if (size(sig) .ne. npt) stop "Error 2 in pfit"
    if (np .gt. npt) stop "Error 4 in pfit"
    if (present(coeff)) then
       if (size(coeff,1) .ne. npt) stop "Error 6 in pfit"
       if (size(coeff,2) .ne. np) stop "Error 5 in pfit"
    end if
    if (present(cov)) then
       if (size(cov,1) .ne. np) stop "Error 7 in pfit"
       if (size(cov,2) .ne. np) stop "Error 8 in pfit"
    end if

    allocate(C(npt,np), Q(npt,np), R(np,np), b(np))

    ! Vandermonde matrix
    do jp=1,np
       do k=1,npt
          C(k,jp) = x(k)**(p(jp))/sig(k)
       end do
    end do

    ! QR decomposition
    call DQRF(C,Q,R,work)

    ! Inversion of R factor
    R = inv(R)

    ! Compute max-likelihood parameters
    ! a = R^-1 Q^T y/σ
    b = 0.d0
    do jp=1,np
       do k=1,npt
          b(jp) = b(jp) + Q(k,jp) * y(k) / sig(k)
       end do
    end do

    a = 0.d0
    do ip=1,np
       do jp=1,np
          a(ip) = a(ip) + R(ip,jp) * b(jp)
       end do
    end do

    ! Compute coefficient matrix coeff such that a(i) = Σ_j coeff(j,i) y(j)
    ! Here a(i) = R^{-1}(i,j) Q(k,j) y(k)/σ(k)
    ! So coeff(k,i) = R^{-1}(i,j) Q(k,j) / σ(k)
    if (present(coeff) .or. present(cov)) then
       coeff1 = 0.d0
       do ip=1,np
          do jp=1,np
             do k=1,npt
                coeff1(k,ip) = coeff1(k,ip) + R(ip,jp) * Q(k,jp) / sig(k)
             end do
          end do
       end do
       if (present(coeff)) coeff = coeff1
    end if

    ! Compute covariance matrix Cov(a(i),a(j)) = Σ_k C(k,i) C(k,j) σ(k)^2
    if (present(cov)) then
       cov = 0.d0
       do jp=1,np
          do ip=1,np
             do k=1,npt
                cov(ip,jp) = cov(ip,jp) + coeff1(k,ip) * coeff1(k,jp) * sig(k)**2
             end do
          end do
       end do
    end if

    ! Compute sqrt(chi^2/ndf)
    if (present(chi)) then
       chi = 0.d0
       do k=1,npt
          val = 0.d0
          do ip=1,np
             val = val + a(ip) * x(k)**p(ip)
          end do
          chi = chi + (val - y(k))**2/sig(k)**2
       end do
       chi = sqrt(chi/(npt - np))
    end if
  end subroutine pfit


  subroutine DQRF(A,Q,R,work)
    ! Compute the QR factorization of a general real matrix A:
    !   A = Q R
    ! where Q is unitary and R is upper triangular, using the LAPACK routine
    ! zgeqrf.
    ! Inputs:
    !   A:     Matrix to be factorized, m x n
    ! Ouputs:
    !   Q:     Unitary matrix, m x m
    !   R:     Upper triangular, n x n
    ! Input/output:
    !   work:  real(8) allocatable workspace array. If unallocated, this
    !          routine will allocate it to an appropriate size. If allocated,
    !          it is assumed to be the correct size for this problem.
    implicit none
    real(8), intent(in)  :: A(:,:)
    real(8), intent(out) :: Q(:,:), R(:,:)
    real(8), allocatable :: work(:)

    integer :: m, n, lwork, ierr, i, j
    real(8) :: tau(size(A,2)), qwork(1)
    !real(8) :: tau(size(A,2))
    !integer :: qwork(1) ! this line throws an error in this code
    real(8) :: A1(size(A,1),size(A,2))

    m = size(A,1)
    n = size(A,2)
    if (m .lt. n) stop "Error in DQRF: m < n"
    if (size(Q,1) .ne. m) stop "Error in DQRF (2)"
    if (size(Q,2) .ne. n) stop "Error in DQRF (3)"
    if (size(R,1) .ne. n) stop "Error in DQRF (4)"
    if (size(R,2) .ne. n) stop "Error in DQRF (5)"

    A1 = A
    if (.not. allocated(work)) then
       ! Compute size of workspace
       lwork = -1
       call DGEQRF(m, n, A1, m, TAU, qwork, LWORK, ierr)
       if (ierr .ne. 0) stop "Error calling DGEQRF (1)"
       lwork = qwork(1)
       allocate(work(lwork))
    end if

    lwork = size(work)
    call dgeqrf(m,n,A1,m,tau,work,lwork,ierr)
    if (ierr .ne. 0) stop "Error calling DGEQRF (2)"
    R = 0.d0
    do j=1,n
       do i=1,j
          R(i,j) = A1(i,j)
       end do
    end do
    Q(:,1:n) = A1
    call dorgqr(m,n,n,Q,m,tau,work,lwork,ierr)
    if (ierr .ne. 0) stop "Error calling DORGQR"
  end subroutine DQRF


end module 
