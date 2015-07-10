program test_derivatives

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two,pi
    use timer,           only: tic, toc
    use derivativestuff, only: derivatives
    implicit none

    character(len=4), parameter :: method = "CD10"
    integer, parameter :: nx=512, ny=256, nz=128
    logical, parameter :: ywrkarr=.FALSE., zwrkarr=.TRUE.
    
    type( derivatives ) :: myder

    real(rkind), dimension(:,:,:), allocatable :: f

    real(rkind), dimension(:,:,:), allocatable :: x,df,df_exact
    real(rkind) :: dx, dy, dz

    integer :: i, ierr

    allocate(f(nx ,ny,nz))
    allocate(x(nx ,ny,nz))
    allocate(df(nx ,ny,nz))
    allocate(df_exact(nx ,ny ,nz))

    if ( method .NE. "CHEB" ) then
        do i=1,nx
            x(i,:,:) = real(i-1,rkind)*two*pi/real(nx,rkind)
            f(i,:,:) = sin(4._rkind * x(i,:,:))
            df_exact(i,:,:) = 4._rkind * cos( 4._rkind * x(i,:,:))
        end do
        dx = two*pi/real(nx,rkind)
        dy = two*pi/real(ny,rkind)
        dz = two*pi/real(nz,rkind)
    else
        do i=1,nx
            x(i,:,:) = - cos( real((i-1),rkind) * pi / real((nx-1),rkind) )
            f(i,1:ny,1:nz) = cos( x(i,:,:) )
            df_exact(i,:,:) = -sin( x(i,:,:) )
        end do
    end if

    ierr = myder % init( nx, ny, nz, dx, dy, dz, .TRUE., .TRUE., .TRUE., method, method, method )

    ! call tic()
    !     df = myder % ddz(f)
    ! call toc("Time to get z derivatives")
    ! print*, "Maximum error = ", MAXVAL( ABS(df) )

    ! call tic()
    !     df = myder % ddy(f)
    ! call toc("Time to get y derivatives")
    ! print*, "Maximum error = ", MAXVAL( ABS(df) )

    ! X derivatives
    call tic()
    df = myder % ddx(f)
    call toc("Time to get x derivatives")
    print*, "Maximum error = ", MAXVAL( ABS(df - df_exact) )

    deallocate(x,f,df,df_exact)
end program
