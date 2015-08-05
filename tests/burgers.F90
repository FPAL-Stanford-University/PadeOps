! 1D periodic viscous burgers equation
program burgers

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two,pi
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    
    implicit none

    integer :: nx = 64, ny = 1, nz = 1

    type( derivatives ) :: der
    type( filters     ) :: fil

    real(rkind), dimension(:,:,:), allocatable :: x,u,u0

    real(rkind) :: nu = 0.01_rkind
    real(rkind) :: dx,dt,t,tend

    character(len=*), parameter :: dermethod = "cd10"    ! Use 10th order Pade
    character(len=*), parameter :: filmethod = "cf90"    ! Use 8th order filter

    integer :: i, iounit=17

    allocate( x(nx,ny,nz) )
    allocate( u(nx,ny,nz) )
   
    dx = one / real(nx,rkind)

    do i=1,nx
        x(i,1,1) = real(i-1,rkind)*dx
        u(i,1,1) = sin( two*pi*x(i,1,1) )
    end do

    u0 = u

    call der%init(          nx,     ny,     nz, &
                            dx,    one,    one, &
                        .TRUE., .TRUE., .TRUE., &
                     dermethod, "cd10", "cd10"  )

    call fil%init(          nx,     ny,     nz, &
                        .TRUE., .TRUE., .TRUE., &
                     filmethod, "cf90", "cf90"  )

    t = zero
    dt = 0.001_rkind
    tend = 0.25_rkind
    do while ( t .LT. tend )
        call RK45(u,dt,nu,der)
        t = t + dt
    end do
    if ( t .LT. tend ) then
        dt = tend - t
        call RK45(u,dt,nu,der)
        t = t + dt
    end if

    OPEN(UNIT=iounit, FILE="burgets.txt", FORM='FORMATTED')
    WRITE(iounit,'(ES24.16)') t
    do i=1,nx
        WRITE(iounit,'(3ES24.16)') x(i,1,1), u0(i,1,1), u(i,1,1)
    end do
    CLOSE(iounit)

    deallocate( x )
    deallocate( u )

contains

    subroutine GetAdvection(u,adv,der)

        real(rkind), dimension(:,:,:), intent(in) :: u
        class( derivatives ),          intent(in) :: der
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: adv
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: flux

        flux = u*u
        call der%ddx(flux,adv)

    end subroutine

    subroutine GetViscous(u,nu,visc,der)
        
        real(rkind), dimension(:,:,:), intent(in) :: u
        class( derivatives ),          intent(in) :: der
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: visc
        real(rkind)                               :: nu

        call der%d2dx2(u,visc)
        visc = nu*visc
       
    end subroutine

    subroutine calcRHS(u,nu,RHS,der)
        real(rkind), dimension(:,:,:), intent(in)    :: u
        class( derivatives ),          intent(in)    :: der
        real(rkind),                   intent(in)    :: nu
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)), intent(out) :: RHS
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3))              :: tmp

        call GetAdvection(u,tmp,der)
        call GetViscous(u,nu,RHS,der)

        RHS = RHS - tmp ! viscous term - advection term

    end subroutine

    subroutine RK45(u,dt,nu,der)

        real(rkind), dimension(:,:,:), intent(inout) :: u
        class( derivatives ),          intent(in)    :: der
        real(rkind),                   intent(in)    :: dt
        real(rkind),                   intent(in)    :: nu
        real(rkind), dimension(5) ::  A, B
        real(rkind), dimension(SIZE(u,1),SIZE(u,2),SIZE(u,3)) :: RHS, Q
        integer :: step


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
            call calcRHS(u,nu,RHS,der)

            Q = dt*RHS + A(step)*Q;
            u = u + B(step)*Q
        end do

    end subroutine

end program
