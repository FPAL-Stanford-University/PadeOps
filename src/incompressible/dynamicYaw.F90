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
        real(rkind), dimension(:), allocatable :: kw, sigma_0
        real(rkind) :: powerExp=3.d0, eps=1.0D-12, Ct, eta ! Power exponent is 3 for this AD implementation
        integer :: Nx
        ! Wind conditions
        real(rkind) :: wind_speed, wind_direction
        real(rkind), dimension(:), allocatable :: powerObservation, powerBaseline, Popti
        integer :: conditionTurb
        ! Other stuff
        integer, dimension(:), allocatable :: indSorted, unsort
        real(rkind) :: Ptot_initial, Ptot_final
        ! Cache stuff
        real(rkind), dimension(:), allocatable :: u_eff, Cp, a_mom
        real(rkind), dimension(:), allocatable :: ap, delta_u_face_store
        real(rkind), dimension(:,:), allocatable :: gaussianStore, y_c, yCenter
        real(rkind), dimension(:,:), allocatable :: deltaUIndividual, sigmaEnd
        real(rkind) :: A, rho
        ! Backward stuff
        real(rkind), dimension(:), allocatable :: dp_dgamma 
        ! Moving average stuff
        integer :: n_moving_average
        real(rkind), dimension(:,:), allocatable :: power_minus_n, ws_minus_n, pb_minus_n, dir_minus_n
        real(rkind), dimension(:), allocatable :: Phat, Phat_yaw

    contains
        procedure :: init
        procedure :: destroy
        procedure, private :: forward
        procedure, private :: EnKF_update
        procedure, private :: observeField
        procedure, private :: rotate
        procedure, private :: yawOptimize 
        procedure, private :: onlineUpdate 
        procedure, private :: backward
        procedure :: update_and_yaw
        procedure :: simpleMovingAverage
    end type


contains

subroutine init(this, inputfile, xLoc, yLoc, diam, Nt, fixedYaw, dynamicStart, dirType)
    class(dynamicYaw), intent(inout) :: this
    character(len=*), intent(in) :: inputfile
    integer :: ioUnit, conditionTurb, ierr, i, n_moving_average
    real(rkind) :: beta1 = 0.9d0, beta2 = 0.999d0, learning_rate_yaw = 10D-4
    integer :: epochsYaw = 5000, stateEstimationEpochs = 500, Ne = 20, Nx = 500
    real(rkind) :: var_p = 0.04d0, var_k = 9D-6, var_sig = 9D-6, D = 0.315
    real(rkind) :: Ct = 0.75, eta = 0.7
    real(rkind), dimension(:), intent(in) :: xLoc, yLoc
    real(rkind), intent(in) :: diam
    integer, intent(in) :: Nt
    logical, intent(out) :: fixedYaw
    integer, intent(out) :: dynamicStart, dirType
 
    ! Read input file for this turbine    
    namelist /DYNAMIC_YAW/ var_p, var_k, var_sig, epochsYaw, stateEstimationEpochs, & 
                           Ne, Ct, eta, beta1, beta2, conditionTurb, n_moving_average, &
                           fixedYaw, dynamicStart, dirType
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
    allocate(this%deltaUIndividual(this%Nt,this%Nt))
    allocate(this%y_c(this%Nt,this%Nt))
    allocate(this%yCenter(this%Nt,this%Nt))
    allocate(this%turbCenter(this%Nt,2))
    allocate(this%turbCenterStore(this%Nt,2))
    allocate(this%kw(this%Nt))
    allocate(this%sigma_0(this%Nt))
    allocate(this%Phat(this%Nt))
    allocate(this%Phat_yaw(this%Nt))
    allocate(this%Popti(this%Nt))
    allocate(this%power_minus_n(this%n_moving_average,this%Nt))
    allocate(this%ws_minus_n(this%n_moving_average,this%Nt))
    allocate(this%pb_minus_n(this%n_moving_average,this%Nt))
    allocate(this%dir_minus_n(this%n_moving_average,this%Nt))

    ! Define
    this%yaw = 0.d0
    this%kw = 0.1d0
    this%sigma_0 = 0.25
    this%power_minus_n = 0.d0
    this%ws_minus_n = 0.d0
    this%pb_minus_n = 0.d0

    ! Get the wind turbine locations
    do i=1,this%Nt
        this%turbCenter(i,1) = xLoc(i) / this%D
        this%turbCenter(i,2) = yLoc(i) / this%D
    end do
    this%turbCenterStore = this%turbCenter

end subroutine 

subroutine destroy(this)
    class(dynamicYaw), intent(inout) :: this

end subroutine 

subroutine update_and_yaw(this, yaw, wind_speed, wind_direction, powerObservation, ts, powerBaseline, wind_direction_end)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(:), intent(inout) :: yaw
    real(rkind), dimension(:), intent(in) :: powerObservation, powerBaseline
    real(rkind), intent(in) :: wind_speed, wind_direction, wind_direction_end
    integer, intent(in) :: ts

    ! Field data observation
    this%wind_speed = wind_speed ! Get this from the data
    ! Use a locally linear interpolation to get the wind direction at the end of
    ! the step
    this%wind_direction = 270.d0-(wind_direction+wind_direction_end)/2.d0 ! Get this from the data
    this%powerObservation = powerObservation ! Get the power production from ADM code
    this%powerBaseline = powerBaseline ! Get the power production from ADM code
    this%ts = ts
    this%yaw = yaw !- wind_direction*pi/180.d0
    !call this%observeField()
    ! Rotate domain
    call this%rotate()
    ! Normalize the power production
    this%powerObservation = this%powerObservation(this%indSorted) / this%powerBaseline(this%conditionTurb)
    ! Online control update
    call this%onlineUpdate()
    yaw = this%yaw
    ! Restore original layout
    this%turbCenter = this%turbCenterStore

end subroutine

subroutine onlineUpdate(this)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(:), allocatable :: kw, sigma_0, Phat, kwBest, sigmaBest, yaw
    real(rkind), dimension(:), allocatable :: Phat_zeroYaw
    real(rkind), dimension(:,:), allocatable :: psi_kp1
    real(rkind), dimension(this%stateEstimationEpochs) :: error
    real(rkind) :: lowestError
    integer :: nt2, i, t, bestStep

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
        error(t) = sum(abs(Phat-this%powerObservation)) / real(this%Nt)
        if (error(t)<lowestError) then
            lowestError=error(t); bestStep = t;
            kwBest = kw; sigmaBest = sigma_0
            this%Ptot_initial = sum(Phat)
        end if
    end do

    ! Generate optimal yaw angles
    this%kw = kwBest; this%sigma_0 = sigmaBest;
    call this%yawOptimize(kwBest, sigmaBest, yaw)
    ! Store the unsorted parameters
    this%kw = kwBest(this%unsort); this%sigma_0 = sigmaBest(this%unsort)   
    this%yaw = yaw(this%unsort) 
    this%Phat = this%Phat(this%unsort); this%Phat_yaw = this%Phat_yaw(this%unsort)

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
    real(rkind), dimension(:,:), intent(in) :: psi_k
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
    ones_Ne = 1.d0 / real(this%Ne)
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
    
    ! Deallocate
    deallocate(psi_kp,psi_kp1,psi_kpp,chi,psiHat_kp,psiHat_kpp,Phat,ones_Ne)
    deallocate(rbuff,randArr,randArr2, Phat_zeroYaw) 

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
 
    ! eps also determines the termination condition
    k=1; check = 0; 
    bestYaw = 0.d0; Ptot = 0.d0; P_time = 0.d0; P_time = 0.d0; yawTime = 0.d0
    m=0.d0; v=0.d0;
    bestPower = 0.d0; bestYaw = 0.d0; bestPowerOut = 0.d0;
    do while (k < this%epochsYaw .and. check == .false.) 
    
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
    call Qsort(x_for_sorting,1)
    this%indSorted = x_for_sorting%zpos
    this%unsort(indSorted) = this%unsort      
    this%turbCenter = X(this%indSorted,:)

end subroutine

subroutine forward(this, kw, sigma_0, Phat, yaw)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(:), intent(in) :: kw, sigma_0, yaw
    real(rkind), dimension(this%Nt), intent(out) :: Phat
    real(rkind) :: D, R, dx, boundLow, boundHigh, edgeLow, edgeHigh
    integer :: Nx, i, j, k
    real(rkind), dimension(this%Nt) :: delta_v0
    real(rkind) :: dw, aPrev, du, gaussian, delta_u_face, L
    logical, dimension(this%Nt, this%Nt) :: turbinesInFront
    real(rkind), dimension(this%Nx) :: xpFront, delta_v, dwVect
    logical :: check
    
    ! Definitions
    ! Set the relevant turbine length scale
    L = 1.d0
    turbinesInFront = .false.
    this%A = (this%D*L) ** 2. * pi / 4.d0
    D = this%D / this%D
    R = D/2.d0
    ! Atmospheric conditions
    this%rho = 1.d0
    ! Model parameters

    !! Forward prop
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
               check = 0
               if (edgeLow>=boundLow .and. edgeLow<=boundHigh) then
                   check = 1
               end if
               if (edgeHigh>=boundLow .and. edgeHigh<=boundHigh) then
                   check = 1
               end if
               if (check == .TRUE.) then
                   turbinesInFront(j,i) = .true.
               end if
            end if
        end do
        ! sum of squared deficits, loop over front turbines
        delta_u_face = 0.d0
        do k = 1, this%Nt
            if (turbinesInFront(j,k) == .TRUE.) then
                ! sum of squared deficits
                aPrev = 0.5*(1.d0-sqrt(1.d0-this%Ct*cos(yaw(k))**2))
                dx = this%turbCenter(j, 1) - this%turbCenter(k, 1)
                xpFront = linspace(0.d0, dx, this%Nx)
                dw = 1.d0 + kw(k) * log(1.d0 + exp((dx - 2.0 * R) / R));

                du = this%wind_speed * aPrev/(dw**2+this%eps) * & 
                        (1.d0+erf(dx/(R*sqrt(2.d0))));
                this%deltaUIndividual(j,k) = du

                ! Deflection
                dwVect = 1.d0 + kw(k) * log(1.d0 + exp((xpFront - 2.0 * R) / R));
                delta_v0(j) = 0.25 * this%Ct * cos(yaw(k))**2 * sin(yaw(k))
                ! Velocity deficit
                delta_v = (delta_v0(j) / (dwVect**2+this%eps)) * 0.5 * (1.d0 + erf(xpFront & 
                          / (R*sqrt(2.d0))))
                where (xpFront<0.d0)
                    delta_v = 0.d0
                end where
                ! y_c
                !this%y_c(k, j) = integrate(xpFront, -delta_v) !mytrapz(xpFront, -delta_v)
                this%y_c(k, j) = mytrapz(xpFront, -delta_v)
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
                delta_u_face = delta_u_face + du * gaussian;
                this%gaussianStore(k, j) = gaussian;
            end if
        end do
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

end subroutine


subroutine backward(this, kw, sigma_0, yaw)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), dimension(:), intent(in) :: kw, sigma_0, yaw
    integer :: i, j, k
    real(rkind), dimension(this%Nt) :: dp_du_eff, dp_dcp, dcp_dap
    real(rkind), dimension(this%Nt) :: du_eff_dyc, dy_c_dgamma
    real(rkind), dimension(this%Nt) :: dcp_dgamma, da_dgamma
    real(rkind) :: dx, R=0.5, D = 1.d0
    logical :: check
    real(rkind) :: boundLow, boundHigh, edgeLow, edgeHigh, dw, du_eff_da
    logical, dimension(this%Nt) :: turbinesInBack
    real(rkind), dimension(this%Nx) :: dwVect, xpFront

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
        turbinesInBack = .false.
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
               if (check == .true.) then
                   turbinesInBack(j) = .true.
               end if
            end if
        end do
        ! dP / dgamma
        ! Flip over whether the turbine sees a wake or not
        du_eff_da = 0;
        do k = 1, this%Nt
            if (turbinesInBack(k) == .true.) then
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
            end if
        end do    
    end do
    
end subroutine


subroutine simpleMovingAverage(this, meanP, power, meanWs, ws, meanPbaseline, powerBaseline, meanDir, windDir, i, t)
    class(dynamicYaw), intent(inout) :: this
    real(rkind), intent(inout) :: meanP, meanWs, meanPbaseline, meanDir
    real(rkind), intent(in) :: power, ws, powerBaseline, windDir
    integer, intent(in) :: i, t

    if (i>this%n_moving_average) then
        ! Power
        meanP = meanP + (1.d0/real(this%n_moving_average)) * (power - this%power_minus_n(1,t))
        this%power_minus_n(1:this%n_moving_average-1,t) = this%power_minus_n(2:this%n_moving_average,t)
        this%power_minus_n(this%n_moving_average,t) = power
        ! Wind speed
        meanWs = meanWs + (1.d0/real(this%n_moving_average)) * (ws - this%ws_minus_n(1,t))
        this%ws_minus_n(1:this%n_moving_average-1,t) = this%ws_minus_n(2:this%n_moving_average,t)
        this%ws_minus_n(this%n_moving_average,t) = ws
        ! Power baseline
        meanPbaseline = meanPbaseline + (1.d0/real(this%n_moving_average)) * (powerBaseline - this%pb_minus_n(1,t))
        this%pb_minus_n(1:this%n_moving_average-1,t) = this%pb_minus_n(2:this%n_moving_average,t)
        ! Wind direction
        meanDir = meanDir + (1.d0/real(this%n_moving_average)) * (windDir - this%dir_minus_n(1,t))
        this%dir_minus_n(1:this%n_moving_average-1,t) = this%dir_minus_n(2:this%n_moving_average,t)
    else
        ! Power
        this%power_minus_n(i,t) = power
        meanP = sum(this%power_minus_n(1:i,t)) / real(i)
        ! Wind speed
        this%ws_minus_n(i,t) = ws
        meanWs = sum(this%ws_minus_n(1:i,t)) / real(i)
        ! Power baseline
        this%pb_minus_n(i,t) = powerBaseline
        meanPbaseline = sum(this%pb_minus_n(1:i,t)) / real(i)
        ! Wind direction
        this%dir_minus_n(i,t) = windDir
        meanDir = sum(this%dir_minus_n(1:i,t)) / real(i)
    end if


end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LAPACK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

end module 
