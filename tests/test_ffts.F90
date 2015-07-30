program test_cd10

    use kind_parameters, only: rkind
    use constants,       only: two,pi
    use timer,           only: tic,toc
    use fftstuff,            only: ffts
    implicit none

    integer :: nx = 128, ny=128, nz=128

    logical, parameter :: periodic = .TRUE.

    type( ffts ) :: xfour, yfour, zfour
    real(rkind), dimension(:,:,:), allocatable :: x,y,z,f,df,dfdx_exact,dfdy_exact,dfdz_exact
    real(rkind) :: dx, dy, dz
    real(rkind), parameter :: omega = 1._rkind

    integer :: i,j,k,ierr

    allocate( x(nx,ny,nz) )
    allocate( y(nx,ny,nz) )
    allocate( z(nx,ny,nz) )
    allocate( f(nx,ny,nz) )
    allocate( df(nx,ny,nz) )
    allocate( dfdx_exact(nx,ny,nz) )
    allocate( dfdy_exact(nx,ny,nz) )
    allocate( dfdz_exact(nx,ny,nz) )

    if (periodic) then
        dx = two*pi/real(nx,rkind)
    else
        stop "ffts type does not support non-peridic bc"
    end if
    dy = two*pi/real(ny,rkind)
    dz = two*pi/real(nz,rkind)

    do k=1,nz
    do j=1,ny
    do i=1,nx
        x(i,j,k) = real(i-1,rkind)*dx
        y(i,j,k) = real(j-1,rkind)*dy
        z(i,j,k) = real(k-1,rkind)*dz
        f(i,j,k) = sin(omega * x(i,j,k)) + sin(omega * y(i,j,k)) + sin(omega * z(i,j,k))
        dfdx_exact(i,j,k) = omega * cos( omega * x(i,j,k))
        dfdy_exact(i,j,k) = omega * cos( omega * y(i,j,k)) 
        dfdz_exact(i,j,k) = omega * cos( omega * z(i,j,k)) 
    end do
    end do 
    end do 

    print*, "Created initial data"

    ierr = xfour%init( nx, "x", ny, nz, dx)
    ierr = yfour%init( ny, "y", nx, nz, dy, .false.)
    ierr = zfour%init( nz, "z", nx, ny, dz, .true.)

    print*, "Initialized"
    call tic() 
    call xfour % dd1(f,df)
    call toc ("Time to get the x derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdx_exact))

    call tic() 
    call yfour % dd2(f,df)
    call toc ("Time to get the y derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))

    call tic() 
    call zfour % dd3(f,df)
    call toc ("Time to get the z derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdz_exact))

    call xfour%destroy
    call yfour%destroy
    call zfour%destroy
    deallocate( x )
    deallocate( f )
    deallocate( df )
    deallocate( dfdx_exact )
    deallocate( dfdy_exact )
    

end program
