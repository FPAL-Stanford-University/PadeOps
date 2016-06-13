! 1D shock tube problem
module StiffShockMod
    
    use kind_parameters, only: rkind
    use constants,       only: zero,eps,half,one,two,three,eight,pi
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit,nancheck

    implicit none

    ! EOS parameters
    real(rkind) :: Rgas = one           ! Gas constant non-dimensionalized to one
    real(rkind) :: gam = 2.84_rkind     ! Ratio of specific heats for the gas
    real(rkind) :: pInf                 ! P_infty for the stiffened gas
    
    ! Problem parameters
    real(rkind) :: pRatio = 4.5_rkind   ! Pressure ratio P2/P1
    real(rkind) :: pInfbyP1 = 1000._rkind ! P_infty/P1 ratio for the stiffened gas
    
    real(rkind) :: rho1 = one           ! Right density
    real(rkind) :: rho2                 ! Left density
    
    real(rkind) :: u1                   ! Right velocity
    real(rkind) :: u2                   ! Left velocity

    real(rkind) :: p1 = one             ! Right pressure
    real(rkind) :: p2                   ! Left pressure
    
    real(rkind) :: E1                   ! Right total energy
    real(rkind) :: E2                   ! Left total energy

    real(rkind) :: tstop = 20._rkind    ! Stop time for simulation
    real(rkind) :: dt    = 0.0001_rkind ! Time step to use for the simulation
    real(rkind) :: dt_fixed = real(1.0D-6, rkind) ! Time step to use for the simulation
    real(rkind) :: CFL   = half         ! CFL number to use for the simulation
    real(rkind) :: t     = zero         ! Simulation time
    
    integer :: nx = 201, ny = 1, nz = 1 ! Number of points to use for the simulation (ny and nz have to be 1)
    real(rkind) :: dx                   ! Grid spacing

    character(len=*), parameter :: dermethod = "cd10"    ! Use 10th order Pade for derivatives
    character(len=*), parameter :: filmethod = "cf90"    ! Use 8th order 90% filter
    
    type( filters ) :: gfil
   
    ! Parameters to control LAD method
    logical, parameter :: UseExpl4thDer = .FALSE.
    logical, parameter :: conservative = .TRUE. 
    real(rkind) :: Cmu = 0.002_rkind
    real(rkind) :: Cbeta = 1.75_rkind
    real(rkind) :: Ckap = 0.00_rkind

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

        p = (gam-one)*( u(:,:,:,3) - half*u(:,:,:,2)*u(:,:,:,2)/u(:,:,:,1) ) - gam*pInf ! (gam-1)*rho*e - gam*P_infty

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
        T = e / (Rgas / (gam-one))   ! T = e / Cv
        
        call GetPressure(u,cs)  ! cs houses pressure now
        cs = sqrt( gam*(cs + pInf)/u(:,:,:,1) )  ! Speed of sound = sqrt( gamma*(p+pInf)/rho )

        call der%ddx(T, dTdx)
        
        ! Get artificial shear viscosity
        if (.NOT. UseExpl4thDer) then
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

    subroutine GetViscous(u,visc,mu,bulk,kap,der)
        
        real(rkind), dimension(:,:,:,:), intent(in) :: u
        class( derivatives ),            intent(in) :: der
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3),SIZE(u,4)), intent(out) :: visc
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: tmp
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: dudx
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: d2udx2
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: dTdx
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: mu
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: bulk
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: kap

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

    subroutine calcRHS(u,RHS,mu,bulk,kap,der)
        real(rkind), dimension(:,:,:,:), intent(in)    :: u
        class( derivatives ),          intent(in)    :: der
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3),SIZE(u,4)), intent(out) :: RHS
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3),SIZE(u,4))              :: tmp
        logical, parameter :: viscous = .TRUE.
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: mu 
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: bulk
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: kap

        call GetAdvection(u,RHS,der)
        if (viscous) then
            call GetViscous(u,tmp,mu,bulk,kap,der)
        else
            tmp = zero
        end if

        RHS = tmp - RHS ! viscous term - advection term

    end subroutine

    subroutine setBC(u)
        real(rkind), dimension(:,:,:,:), intent(inout) :: u

        u(1,:,:,1) = rho2;    u(nx,:,:,1) = rho1      ! Mass
        u(1,:,:,2) = rho2*u2; u(nx,:,:,2) = rho1*u1   ! Momentum
        u(1,:,:,3) = E2;      u(nx,:,:,3) = E1        ! Energy

    end subroutine

    subroutine get_dt(u,mu,bulk,kap,dx,dt,stability)
        real(rkind), dimension(:,:,:,:), intent(in)  :: u
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(in) :: mu
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(in) :: bulk
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(in) :: kap
        real(rkind),                     intent(in)  :: dx
        real(rkind),                     intent(out) :: dt
        character(len=*),                intent(out) :: stability

        real(rkind), dimension(nx,1,1)               :: cs, T
        real(rkind) :: dtCFL, dtmu, dtbulk, dtkap

        call GetInternalEnergy(u,cs) ! cs houses e now
        T = cs / (Rgas / (gam-one))   ! T = e / Cv
        call GetPressure(u,cs)  ! cs houses pressure now
        cs = sqrt(gam*(cs + pInf)/u(:,:,:,1))  ! Speed of sound = sqrt( gamma*(p+pInf)/rho )
        
        dtCFL  = CFL * dx / (MAXVAL( ABS(u(:,:,:,2)/u(:,:,:,1)) + ABS(cs) ))
        dtmu   = 0.2_rkind * dx**2 / (MAXVAL( mu/u(:,:,:,1) ) + eps)
        dtbulk = 0.2_rkind * dx**2 / (MAXVAL( bulk/u(:,:,:,1) ) + eps)
        dtkap  = 0.2_rkind * one / (MAXVAL( kap*T/(u(:,:,:,1)* (cs**2) * (dx**2)) ) + eps)

        if ( CFL .LE. zero ) then
            dt = dt_fixed
            stability = 'fixed'
        else
            stability = 'convective'
            dt = dtCFL
            if ( dt > dtmu ) then
                dt = dtmu
                stability = 'shear'
            else if ( dt > dtbulk ) then
                dt = dtbulk
                stability = 'bulk'
            else if ( dt > dtkap ) then
                dt = dtkap
                stability = 'conductive'
            end if
        end if

    end subroutine

    subroutine RK45(u,dt,stability,der,fil)

        real(rkind), dimension(:,:,:,:), intent(inout) :: u
        class( derivatives ),            intent(in)    :: der
        class( filters     ),            intent(in)    :: fil
        real(rkind),                     intent(inout) :: dt
        character(len=*),                intent(out)   :: stability
        real(rkind), dimension(5) ::  A, B
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3),SIZE(u,4)) :: RHS, Q
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: mu
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: bulk
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: kap
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
            ! Get the RHS
            call calcRHS(u,RHS,mu,bulk,kap,der)
    
            if (step == 1) then
                call get_dt(u,mu,bulk,kap,dx,dt,stability)
                if ( (t .LE. tstop) .AND. ( t .GE. (tstop - dt) ) ) then
                    dt = tstop - t
                    stability = 'final'
                end if
            end if

            ! Update the solution
            Q = dt*RHS + A(step)*Q;
            u = u + B(step)*Q

            ! Filter the solution for partial de-aliasing
            if (filteron) then
                call fil%filterx(u(:,:,:,1),RHS(:,:,:,1))
                call fil%filterx(u(:,:,:,2),RHS(:,:,:,2))
                call fil%filterx(u(:,:,:,3),RHS(:,:,:,3))
                u = RHS
            end if

            ! Set boundary conditions
            call setBC(u)
        end do

        if ( nancheck(u) ) then
            call GracefulExit("Oh the devil! NaN encountered. Exiting...",666)
        end if

    end subroutine

end module



program StiffShock

    use kind_parameters, only: rkind,clen
    use constants,       only: zero,half,one,two
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use StiffShockMod,   only: nx,ny,nz,dx,tstop,dt,t,rho1,rho2,p1,p2,u1,u2,E1,E2,gam,pInf,pInfbyP1,pRatio,gfil,dermethod,filmethod, &
                               RK45,GetPressure,GetInternalEnergy,GetSGS

    implicit none


    type( derivatives ) :: der
    type( filters     ) :: fil

    real(rkind), dimension(:,:,:), allocatable :: x, dum, mu, bulk, kap, dudx, dTdx
    real(rkind), dimension(:,:,:,:), allocatable :: u

    real(rkind) :: rho2rho1, rhoe1, rhoe2
    real(rkind) :: small
    
    integer :: step
    character(len=clen) :: stability

    integer :: i, iounit=17

    allocate(    x(nx,ny,nz)   )
    allocate(  dum(nx,ny,nz)   )
    allocate(   mu(nx,ny,nz)   )
    allocate( bulk(nx,ny,nz)   )
    allocate(  kap(nx,ny,nz)   )
    allocate( dudx(nx,ny,nz)   )
    allocate( dTdx(nx,ny,nz)   )
    allocate(    u(nx,ny,nz,3) )
   
    dx = one / real(nx-1,rkind) 

    ! Create the solution grid on [-0.5,0.5] domain
    do i=1,nx
        x(i,1,1) = real(i-1,rkind)*dx - half
    end do

    ! Get pInf
    pInf = pInfbyP1 * p1

    ! Calculate density ratio
    rho2rho1 = ( (gam+one)*pRatio + (two*gam)*pInfbyP1 + (gam-one) ) / ( (gam-one)*pRatio + (two*gam)*pInfbyP1 + (gam+one) )
    rho2 = rho1 * rho2rho1

    ! Calculate p2
    p2 = p1*pRatio

    ! Calculate u1
    u1 = - sqrt(rho2rho1 * (p2 - p1) / (rho2 - rho1)) ! minus since fluid moving to the left
    u2 = rho1*u1 / rho2

    ! Calculate total energy states
    rhoe1 = (p1 + gam*pInf)/(gam-one)
    rhoe2 = (p2 + gam*pInf)/(gam-one)
    E1 = rhoe1 + half * rho1 * u1 * u1 ! E = rho*e + 0.5*rho*u*u
    E2 = rhoe2 + half * rho2 * u2 * u2 ! E = rho*e + 0.5*rho*u*u

    small = one
    dum = half * ( one+tanh( x/(small*dx)) )

    ! Initialize conserved variables
    u(:,:,:,1) = (one-dum)*rho2 + dum*rho1        ! rho
    u(:,:,:,2) = (one-dum)*rho2*u2 + dum*rho1*u1  ! rho*u
    u(:,:,:,3) = (one-dum)*E2 + dum*E1            ! E

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

    ! call fil%filterx(u(:,:,:,1),dum)
    ! u(:,:,:,1) = dum
    ! call fil%filterx(u(:,:,:,2),dum)
    ! u(:,:,:,2) = dum
    ! call fil%filterx(u(:,:,:,3),dum)
    ! u(:,:,:,3) = dum

    call GetPressure(u,dum)  ! dum houses pressure now
    dum = sqrt( gam*(dum + pInf)/u(:,:,:,1) )  ! Speed of sound = sqrt( gamma*(p+pInf)/rho )
    tstop = tstop / minval(dum)

    ! Integrate in time
    step = 0
    t = zero
    print '(A6,6(A2,A13))', 'step', '|', 'time', '|', 'timestep', '|', 'stability', '|', 'min density', '|', 'max velocity', '|', 'max energy'
    do while ( t .LT. (tstop - dt) )
        call RK45(u,dt,stability,der,fil)
        t = t + dt
        step = step + 1
        call GetInternalEnergy(u,dum)
        print '(I6,2(A2,ES13.5),1(A2,A13),3(A2,ES13.5))', step, '|', t, '|', dt, '|', trim(stability), '|', MINVAL(u(:,:,:,1)), '|', MAXVAL(ABS(u(:,:,:,2)/u(:,:,:,1))), '|', MAXVAL(dum)
    end do
    ! Special case for the last time step
    if ( t .LT. tstop ) then
        dt = tstop - t
        call RK45(u,dt,stability,der,fil)
        t = t + dt
        step = step + 1
        call GetInternalEnergy(u,dum)
        print '(I6,2(A2,ES13.5),1(A2,A13),3(A2,ES13.5))', step, '|', t, '|', dt, '|', trim(stability), '|', MINVAL(u(:,:,:,1)), '|', MAXVAL(ABS(u(:,:,:,2)/u(:,:,:,1))), '|', MAXVAL(dum)
    end if
    
    print*, "Done."
    print*, "Writing out visualization file 'StiffShock.dat'"

    ! Get the SGS terms for output    
    call GetSGS(u,dudx,dTdx,mu,bulk,kap,der)

    call GetPressure(u,dum)  ! Put pressure in dum for output
    OPEN(UNIT=iounit, FILE="StiffShock.dat", FORM='FORMATTED')
    WRITE(iounit,'(A, ES24.16)') "Final time = ", t
    WRITE(iounit,'(8A24)') "X", "Density", "Velocity", "Pressure", "Energy", "Shear Visc.", "Bulk Visc.", "Conductivity"
    do i=1,nx
        WRITE(iounit,'(8ES24.16)') x(i,1,1), u(i,1,1,1), u(i,1,1,2)/u(i,1,1,1), dum(i,1,1), &
                                  ( u(i,1,1,3) - half*u(i,1,1,2)*u(i,1,1,2)/u(i,1,1,1) )/u(i,1,1,1), &
                                   mu(i,1,1), bulk(i,1,1), kap(i,1,1) ! x, density, velocity, pressure, shear visc., bulk visc., conductivity
    end do
    CLOSE(iounit)

    call der%destroy()
    call fil%destroy()
    call gfil%destroy()
    deallocate( x )
    deallocate( dum )
    deallocate(   mu )
    deallocate( bulk )
    deallocate(  kap )
    deallocate( dudx )
    deallocate( dTdx )
    deallocate( u )

end program
