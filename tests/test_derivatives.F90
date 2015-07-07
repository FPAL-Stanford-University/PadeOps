program test_derivatives

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two,pi
    use timer,           only: tic, toc
    use derivativestuff, only: derivatives
    implicit none

    character(len=4), parameter :: method = "CD10"
    integer, parameter :: nx=512, ny=128, nz=128
    integer, parameter :: xpadding=16, ypadding=0

    type( derivatives ) :: myder

    real(rkind), dimension(nx+xpadding,ny+ypadding,nz), target :: f_data
    real(rkind), dimension(:,:,:), pointer :: f

    real(rkind), dimension(nx,ny,nz) :: x,df,df_exact
    real(rkind) :: dx, dy, dz

    integer :: i, ierr

    f => f_data(1:nx,1:ny,1:nz)

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
            f(i,:,:) = cos( x(i,:,:) )
            df_exact(i,:,:) = -sin( x(i,:,:) )
        end do
    end if

    ierr = myder % init( nx, ny, nz, dx, dy, dz, .TRUE., .TRUE., .TRUE., method, method, method )

    call tic()
    df = myder % ddz(f)
    call toc("Time to get z derivatives")
    print*, "Maximum error = ", MAXVAL( ABS(df) )

    call tic()
    df = myder % ddy(f)
    call toc("Time to get y derivatives")
    print*, "Maximum error = ", MAXVAL( ABS(df) )

    call tic()
    df = myder % ddx(f)
    call toc("Time to get x derivatives")
    print*, "Maximum error = ", MAXVAL( ABS(df - df_exact) )

    nullify(f)

end program
