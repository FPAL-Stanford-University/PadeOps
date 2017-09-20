program test_ed02

    use kind_parameters, only: rkind
    use constants,       only: two,pi
    use timer,           only: tic,toc
    use ed02stuff,       only: ed02
    implicit none

    integer :: nx = 64, ny=64, nz=64

    logical, parameter :: periodic = .TRUE.

    type( ed02 ) :: xed02, yed02, zed02
    real(rkind), dimension(:,:,:), allocatable :: x,y,z,f,df,df_exact,d2f,d2f_exact
    real(rkind) :: dx, dy, dz
    real(rkind), parameter :: omega = 1._rkind
    integer, dimension(0:1) :: x_bc = [1,1], y_bc = [1,1], z_bc = [1,1]

    integer :: i,ierr, j, k

    allocate( x(nx,ny,nz) )
    allocate( y(nx,ny,nz) )
    allocate( z(nx,ny,nz) )
    allocate( f(nx,ny,nz) )
    allocate( df(nx,ny,nz) )
    allocate( df_exact(nx,ny,nz) )
    allocate( d2f(nx,ny,nz) )
    allocate( d2f_exact(nx,ny,nz) )

    if (periodic) then
        dx = two*pi/real(nx,rkind)
        dy = two*pi/real(ny,rkind)
        dz = two*pi/real(nz,rkind)
    else
        dx = two*pi/real(nx-1,rkind)
        dy = two*pi/real(ny-1,rkind)
        dz = two*pi/real(nz-1,rkind)
    end if
    

    do k=1,nz
        do j=1,ny
            do i=1,nx
                x(i,j,k) = real(i-1,rkind)*dx
                y(i,j,k) = real(j-1,rkind)*dy
                z(i,j,k) = real(k-1,rkind)*dz
            end do
        end do 
    end do 
    print*, "Created initial data"

    f = cos(omega * x) * cos(omega * y) * cos(omega * z)
    ierr = xed02%init( nx, dx, periodic, 0, 0)
    ierr = yed02%init( ny, dy, periodic, 0, 0)
    ierr = zed02%init( nz, dz, periodic, 0, 0)

    call tic() 
    call xed02 % dd1(f,df,ny,nz,x_bc(0),x_bc(1))
    call xed02 % d2d1(f,d2f,ny,nz,x_bc(0),x_bc(1))
    call toc ("Time to get the x derivative:")
    df_exact =-omega* sin(omega * x) * cos(omega * y) * cos(omega * z)
    d2f_exact = -omega* omega* cos(omega * x) * cos(omega * y) * cos(omega * z)
    print*, "Maximum error (first  der) = ", MAXVAL( ABS(df - df_exact))
    print*, "Maximum error (second der) = ", MAXVAL( ABS(d2f - d2f_exact))

    call tic() 
    call yed02 % dd2(f,df,nx,nz,y_bc(0),y_bc(1))
    call yed02 % d2d2(f,d2f,nx,nz,y_bc(0),y_bc(1))
    call toc ("Time to get the y derivative:")
    df_exact =-omega* sin(omega * y) * cos(omega * x) * cos(omega * z)
    d2f_exact = -omega* omega* cos(omega * y) * cos(omega * x) * cos(omega * z)
    print*, "Maximum error (first  der) = ", MAXVAL( ABS(df - df_exact))
    print*, "Maximum error (second der) = ", MAXVAL( ABS(d2f - d2f_exact))

    call tic() 
    call zed02 % dd3(f,df,nx,ny,z_bc(0),z_bc(1))
    call zed02 % d2d3(f,d2f,nx,ny,z_bc(0),z_bc(1))
    call toc ("Time to get the z derivative:")
    df_exact =-omega* sin(omega * z) * cos(omega * x) * cos(omega * y)
    d2f_exact = -omega* omega* cos(omega * z) * cos(omega * x) * cos(omega * y)
    print*, "Maximum error (first  der) = ", MAXVAL( ABS(df - df_exact))
    print*, "Maximum error (second der) = ", MAXVAL( ABS(d2f - d2f_exact))

    call xed02%destroy()
    call yed02%destroy()
    call zed02%destroy()

    deallocate( x )
    deallocate( f )
    deallocate( df )
    deallocate( df_exact )

end program
