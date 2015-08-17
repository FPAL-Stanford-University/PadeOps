! 1D shock tube problem
module ShocktubeMod
    
    use kind_parameters, only: rkind
    use constants,       only: zero,half,one,two,eight,pi
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters

    implicit none

    real(rkind) :: Rgas = one           ! Gas constant non-dimensionalized to one
    real(rkind) :: gam = 1.4_rkind      ! Ratio of specific heats for the gas

    real(rkind) :: xshk = -4._rkind     ! Location of initial discontinuity

    real(rkind) :: rhoL = 3.857143_rkind! Left density
    real(rkind) :: rhoR                 ! Right density

    real(rkind) :: uL = 2.629369_rkind  ! Left velocity
    real(rkind) :: uR = zero            ! Right velocity

    real(rkind) :: pL = 10.33333_rkind  ! Left pressure
    real(rkind) :: pR = one             ! Right pressure

    real(rkind) :: EL                   ! Left total energy
    real(rkind) :: ER                   ! Right total energy

    real(rkind) :: tstop = 1.8_rkind    ! Stop time for simulation
    real(rkind) :: dt    = 0.0001_rkind ! Time step to use for the simulation
    
    integer :: nx = 1601, ny = 1, nz = 1 ! Number of points to use for the simulation (ny and nz have to be 1)
    real(rkind) :: dx

    character(len=*), parameter :: dermethod = "cd10"    ! Use 10th order Pade for derivatives
    character(len=*), parameter :: filmethod = "cf90"    ! Use 8th order 90% filter
    
    type( derivatives ) :: der
    type( filters     ) :: fil
    type( filters     ) :: gfil
    
    real(rkind) :: Cmu = 0.002_rkind                 ! Coefficient for artificial shear viscosity
    real(rkind) :: Cbeta = 1.75_rkind                ! Coefficient for artificial bulk viscosity
    real(rkind) :: Ckap = 0.01_rkind                 ! Coefficient for artificial heat conductivity
    logical, parameter :: UseExpl4thDer = .FALSE.    ! Use explicit 4th order derivatives for art. props?
    logical, parameter :: conservative = .FALSE.     ! Use conservative formulation for the viscous terms?

contains

    subroutine expl4_d4dx4(f,df,dx)
        real(rkind), dimension(:,:,:), intent(in) :: f
        real(rkind), dimension(size(f,1),size(f,2),size(f,3)), intent(out)  :: df
        real(rkind), intent(in) :: dx
        integer :: nx
    
        real(rkind), parameter :: a = 28._rkind/3._rkind
        real(rkind), parameter :: b =-13._rkind/2._rkind
        real(rkind), parameter :: c =  2._rkind/1._rkind
        real(rkind), parameter :: d = -1._rkind/6._rkind
    
        nx = size(f,1)
    
        df(    1:3,:,:) = zero
        df( 4:nx-3,:,:) = a * ( f(4:nx-3,:,:)                 ) &
                        + b * ( f(5:nx-2,:,:) + f(3:nx-4,:,:) ) &
                        + c * ( f(6:nx-1,:,:) + f(2:nx-5,:,:) ) &
                        + d * ( f(7:nx-0,:,:) + f(1:nx-6,:,:) )
        df(nx-2:nx,:,:) = zero
    
        df = df / dx**4
    
    end subroutine

    pure subroutine GetPressure(u,p)

        real(rkind), dimension(:,:,:,:), intent(in) :: u
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: p

        p = (gam-one)*( u(:,:,:,3) - half*u(:,:,:,2)*u(:,:,:,2)/u(:,:,:,1) )

    end subroutine

    pure subroutine GetInternalEnergy(u,e)

        real(rkind), dimension(:,:,:,:), intent(in) :: u
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: e

        e = ( u(:,:,:,3) - half*u(:,:,:,2)*u(:,:,:,2)/u(:,:,:,1) )/u(:,:,:,1)

    end subroutine

    subroutine GetAdvection(u,adv,der)

        real(rkind), dimension(:,:,:,:), intent(in) :: u
        class( derivatives ),            intent(in) :: der
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3),SIZE(u,4)), intent(out) :: adv
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3),SIZE(u,4)) :: flux
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3))           :: p

        call GetPressure(u,p)

        flux(:,:,:,1) = u(:,:,:,2)                                  ! rho*u
        flux(:,:,:,2) = u(:,:,:,2)*u(:,:,:,2)/u(:,:,:,1) + p        ! rho*u^2 + p
        flux(:,:,:,3) = u(:,:,:,2) * (u(:,:,:,3) + p) / u(:,:,:,1)  ! u*(E+p)

        call der%ddx(flux(:,:,:,1),adv(:,:,:,1))
        call der%ddx(flux(:,:,:,2),adv(:,:,:,2))
        call der%ddx(flux(:,:,:,3),adv(:,:,:,3))

    end subroutine

    subroutine GetSGS(u,dudx,dTdx,mu,bulk,kap,der)
        real(rkind), dimension(:,:,:,:), intent(in) :: u
        class( derivatives ),          intent(in) :: der
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: dudx
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: dTdx
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: mu
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: bulk
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: kap
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: tmp
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: e
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: T
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: cs

        ! Get the derivatives
        tmp = u(:,:,:,2)/u(:,:,:,1)
        call der%ddx(tmp, dudx)
        
        call GetInternalEnergy(u,e)

        call GetPressure(u,cs)
        T = cs/( u(:,:,:,1) * Rgas )   ! T = p/(rho*R)
        cs = gam*cs/u(:,:,:,1)  ! Speed of sound = gamma*p/rho

        call der%ddx(T, dTdx)
        
        ! Get artificial shear viscosity
        if (.NOT. UseExpl4thDer) then
            !call fil%filterx(abs(dudx),mu)
            call der%d2dx2(abs(dudx),tmp)
            call der%d2dx2(tmp,mu)
        else
            call expl4_d4dx4(abs(dudx),mu,dx)
        end if
        tmp = Cmu*u(:,:,:,1)*abs(mu*dx**6)
        call gfil%filterx(tmp,mu)
        if (.NOT. UseExpl4thDer) then
            call gfil%filterx(mu,tmp)
            mu = tmp
        end if

        ! Get artificial bulk viscosity
        if (.NOT. UseExpl4thDer) then
            !call fil%filterx(dudx,bulk)
            call der%d2dx2(dudx,tmp)
            call der%d2dx2(tmp,bulk)
        else
            call expl4_d4dx4(dudx,bulk,dx)
        end if
        tmp = Cbeta*u(:,:,:,1)*abs(bulk*dx**6) * (MIN(dudx,zero)/(dudx+1.0D-30))
        call gfil%filterx(tmp,bulk)
        if (.NOT. UseExpl4thDer) then
            call gfil%filterx(bulk,tmp)
            bulk = tmp
        end if

        ! Get artificial conductivity
        if (.NOT. UseExpl4thDer) then
            !call fil%filterx(e,kap)
            call der%d2dx2(e,tmp)
            call der%d2dx2(tmp,kap)        
        else
            call expl4_d4dx4(e,kap,dx)
        end if
        tmp = Ckap * (u(:,:,:,1)*cs/T) * abs(kap*dx**5)
        call gfil%filterx(tmp,kap)
        if (.NOT. UseExpl4thDer) then
            call gfil%filterx(kap,tmp)
            kap = tmp
        end if

    end subroutine

    subroutine GetViscous(u,visc,der)
        
        real(rkind), dimension(:,:,:,:), intent(in) :: u
        class( derivatives ),          intent(in) :: der
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3),SIZE(u,4)), intent(out) :: visc
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: tmp
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: dudx
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: d2udx2
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: dTdx
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: mu
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: bulk
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: kap

        ! Get SGS properties
        call GetSGS(u,dudx,dTdx,mu,bulk,kap,der)
        
        if ( .NOT. conservative ) then
            visc(:,:,:,1) = zero
            
            ! Get non-conservative d(tau)/dx by expanding out the derivatives
            tmp = ((4._rkind/3._rkind)*mu + bulk ) !
            call der%ddx(tmp,visc(:,:,:,2))
            visc(:,:,:,2) = visc(:,:,:,2)*dudx     ! d(tmp)/dx * du/dx
            call der%d2dx2(u(:,:,:,2)/u(:,:,:,1),d2udx2)
            visc(:,:,:,2) = visc(:,:,:,2) + tmp*d2udx2

            ! Get non-conservative d(tau*u)/dx by expanding out the derivatives ans using above d(tau)/dx

            visc(:,:,:,3) = u(:,:,:,2)/u(:,:,:,1) * visc(:,:,:,2) ! u*d(tau)/dx
            visc(:,:,:,3) = visc(:,:,:,3) + (tmp*dudx)*dudx + kap*dTdx   ! u*d(tau)/dx + tau*du/dx + kap*dT/dx

        else
            visc(:,:,:,1) = zero

            tmp = (4._rkind/3._rkind)*mu*dudx + bulk*dudx
            call der%ddx(tmp,visc(:,:,:,2))

            tmp = tmp*u(:,:,:,2)/u(:,:,:,1) + kap*dTdx
            call der%ddx(tmp,visc(:,:,:,3))

        end if
       
    end subroutine

    subroutine calcRHS(u,RHS,der)
        real(rkind), dimension(:,:,:,:), intent(in)    :: u
        class( derivatives ),          intent(in)    :: der
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3),SIZE(u,4)), intent(out) :: RHS
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3),SIZE(u,4))              :: tmp
        logical, parameter :: viscous = .TRUE.

        call GetAdvection(u,RHS,der)
        if (viscous) then
            call GetViscous(u,tmp,der)
        else
            tmp = zero
        end if

        RHS = tmp - RHS ! viscous term - advection term

    end subroutine

    subroutine RK45(u,dt,der,fil)

        real(rkind), dimension(:,:,:,:), intent(inout) :: u
        class( derivatives ),            intent(in)    :: der
        class( filters     ),            intent(in)    :: fil
        real(rkind),                     intent(in)    :: dt
        real(rkind), dimension(5) ::  A, B
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3),SIZE(u,4)) :: RHS, Q
        integer :: step
        logical, parameter :: filteron = .TRUE.


        A(1) = 0.0_rkind;
        A(2) =  -6234157559845.0_rkind/12983515589748.0_rkind;
        A(3) =  -6194124222391.0_rkind/4410992767914.0_rkind;
        A(4) = -31623096876824.0_rkind/15682348800105.0_rkind;
        A(5) = -12251185447671.0_rkind/11596622555746.0_rkind;

        B(1) =  494393426753.0_rkind/4806282396855.0_rkind;
        B(2) = 4047970641027.0_rkind/5463924506627.0_rkind;
        B(3) = 9795748752853.0_rkind/13190207949281.0_rkind;
        B(4) = 4009051133189.0_rkind/8539092990294.0_rkind;
        B(5) = 1348533437543.0_rkind/7166442652324.0_rkind;


        Q = zero

        do step = 1,5
            call calcRHS(u,RHS,der)

            Q = dt*RHS + A(step)*Q;
            u = u + B(step)*Q

            ! Filter the solution for partial de-aliasing
            if (filteron) then
                call fil%filterx(u(:,:,:,1),RHS(:,:,:,1))
                call fil%filterx(u(:,:,:,2),RHS(:,:,:,2))
                call fil%filterx(u(:,:,:,3),RHS(:,:,:,3))
                u = RHS
            end if

            call setBC(u)

        end do

    end subroutine

    subroutine setBC(u)
        real(rkind), dimension(:,:,:,:), intent(inout) :: u
        integer :: nx

        nx = size(u,1)

        u(1,1,1,1) = rhoL;    u(nx,1,1,1) = rhoR
        u(1,1,1,2) = rhoL*uL; u(nx,1,1,2) = rhoR*uR
        u(1,1,1,3) = EL;      u(nx,1,1,3) = ER

    end subroutine

end module



program ShuOsher

    use kind_parameters, only: rkind
    use constants,       only: zero,half,one
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use ShocktubeMod,    only: nx,ny,nz,dx,tstop,dt,xshk,rhoL,rhoR,uL,uR,pL,pR,EL,ER,gam,der,fil,gfil,dermethod,filmethod, &
                               RK45,GetPressure,GetInternalEnergy,GetSGS

    implicit none


    real(rkind), dimension(:,:,:), allocatable :: x, dum, tmp, mu, bulk, kap, dudx, dTdx
    real(rkind), dimension(:,:,:,:), allocatable :: u

    real(rkind) :: t
    integer :: step

    logical :: FilterInit = .FALSE.

    integer :: i, iounit=17

    allocate(   x(nx,ny,nz)   )
    allocate( dum(nx,ny,nz)   )
    allocate( tmp(nx,ny,nz)   )
    allocate(   mu(nx,ny,nz)   )
    allocate( bulk(nx,ny,nz)   )
    allocate(  kap(nx,ny,nz)   )
    allocate( dudx(nx,ny,nz)   )
    allocate( dTdx(nx,ny,nz)   )
    allocate(   u(nx,ny,nz,3) )
   
    dx = 10._rkind / real(nx-1,rkind) 

    ! Create the solution grid on [-0.5,0.5] domain
    do i=1,nx
        x(i,1,1) = real(i-1,rkind)*dx - 5._rkind
    end do
    
    dum = half * ( one+tanh( (x-xshk)/(0.1_rkind*dx)) )

    tmp = one + 0.2_rkind*sin(5._rkind*x)   ! Right density
    rhoR = tmp(nx,1,1)                      ! Set rhoR for use in setBC

    ! Initialize conserved variables
    u(:,:,:,1) = (one-dum)*rhoL + dum*tmp               ! rho
    u(:,:,:,2) = u(:,:,:,1)*( (one-dum)*uL + dum*uR )   ! rho*u
    
    ! Calculate internal energy states from pressures
    EL = pL/(gam-one) 
    ER = pR/(gam-one)

    u(:,:,:,3) = (one-dum)*EL + dum*ER + half*u(:,:,:,2)*u(:,:,:,2)/u(:,:,:,1)  ! E

    EL = u(1,1,1,3); ER = u(nx,1,1,3)       ! Set EL and ER for use in setBC

    ! Initialize derivatives and filters objects
    call der%init(          nx,       ny,      nz, &
                            dx,      one,     one, &
                        .FALSE., .FALSE., .FALSE., &
                     dermethod,   "cd10",  "cd10"  )

    call fil%init(          nx,       ny,      nz, &
                        .FALSE., .FALSE., .FALSE., &
                     filmethod,   "cf90",  "cf90"  )

    call gfil%init(         nx,       ny,      nz, &
                        .FALSE., .FALSE., .FALSE., &
                     "gaussian",   "cf90",  "cf90"  )

    if (FilterInit) then
        call fil%filterx(u(:,:,:,1),dum)
        u(:,:,:,1) = dum
        call fil%filterx(u(:,:,:,2),dum)
        u(:,:,:,2) = dum
        call fil%filterx(u(:,:,:,3),dum)
        u(:,:,:,3) = dum
    end if

    ! Integrate in time
    step = 0
    t = zero
    print '(A6,4(A2,A13))', 'step', '|', 'time', '|', 'min density', '|', 'max velocity', '|', 'max energy'
    do while ( t .LT. (tstop - dt) )
        call RK45(u,dt,der,fil)
        t = t + dt
        step = step + 1
        call GetInternalEnergy(u,dum)
        print '(I6,4(A2,ES13.5))', step, '|', t, '|', MINVAL(u(:,:,:,1)), '|', MAXVAL(ABS(u(:,:,:,2)/u(:,:,:,1))), '|', MAXVAL(dum)
    end do
    ! Special case for the last time step
    if ( t .LT. tstop ) then
        dt = tstop - t
        call RK45(u,dt,der,fil)
        t = t + dt
        step = step + 1
        call GetInternalEnergy(u,dum)
        print '(I6,4(A2,ES13.5))', step, '|', t, '|', MINVAL(u(:,:,:,1)), '|', MAXVAL(ABS(u(:,:,:,2)/u(:,:,:,1))), '|', MAXVAL(dum)
    end if

    print*, "Done."
    print*, "Writing out visualization file 'ShuOsher.dat'"

    ! Get the SGS terms for output    
    call GetSGS(u,dudx,dTdx,mu,bulk,kap,der)

    call GetPressure(u,dum)  ! Put pressure in dum for output
    OPEN(UNIT=iounit, FILE="ShuOsher.dat", FORM='FORMATTED')
    WRITE(iounit,'(A, ES24.16)') "Final time = ", t
    WRITE(iounit,'(7A24)') "X", "Density", "Velocity", "Pressure", "Shear Visc.", "Bulk Visc.", "Conductivity"
    do i=1,nx
        WRITE(iounit,'(7ES24.16)') x(i,1,1), u(i,1,1,1), u(i,1,1,2)/u(i,1,1,1), dum(i,1,1), mu(i,1,1), bulk(i,1,1), kap(i,1,1)
    end do
    CLOSE(iounit)

    deallocate( x )
    deallocate( dum )
    deallocate( tmp )
    deallocate(   mu )
    deallocate( bulk )
    deallocate(  kap )
    deallocate( dudx )
    deallocate( dTdx )
    deallocate( u )

end program
