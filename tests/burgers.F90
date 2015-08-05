! 1D periodic viscous burgers equation
program burgers

    use kind_parameters, only: rkind
    use constants,       only: one,two,pi
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    
    implicit none

    integer :: nx = 64, ny = 1, nz = 1

    type( derivatives ) :: der
    type( filters     ) :: fil

    real(rkind), dimension(:,:,:), allocatable :: x,u

    real(rkind) :: nu = 0.01_rkind
    real(rkind) :: dx,dt,t,tend

    character(len=*) :: dermethod = "cd10"    ! Use 10th order Pade
    character(len=*) :: filmethod = "cf90"    ! Use 8th order filter

    allocate( x(nx,ny,nz) )
    allocate( u(nx,ny,nz) )
   
    dx = one / real(nx,rkind)

    do i=1,nx
        x(i,1,1) = real(i-1,rkind)*dx
        u(i,1,1) = one !sin( two*pi*x(i,1,1) )
    end do

    call der%init(          nx,     ny,     nz, &
                            dx,    one,    one, &
                        .TRUE., .TRUE., .TRUE., &
                     dermethod, "cd10", "cd10"  )

    call fil%init(          nx,     ny,     nz, &
                        .TRUE., .TRUE., .TRUE., &
                     filmethod, "cd10", "cd10"  )

    t = zero
    dt = 0.001_rkind
    tend = 0.5_rkind
    do while ( t .LT. tend )
        call RK45(u,dt,der)
        t = t + dt
    end do
    if ( t .LT. tend ) then
        dt = tend - t
        call RK45(u,dt,der)
        t = t + dt
    end if

    print*, u(1), one/(one-t)

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

        call der%d2dx2(u,visc)
        visc = nu*visc
       
    end subroutine

    subroutine RK45(u,dt,der)

        real(rkind), dimension(:,:,:), intent(inout) :: u
        class( derivatives ),          intent(in)    :: der
        real(rkind)                                  :: dt
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
            ! call calcRHS(u,nu,RHS,der)

            RHS = u*u

            Q = dt*RHS + A(step)*Q;
            u = u + B(step)*Q
        end do

    end subroutine

end program
