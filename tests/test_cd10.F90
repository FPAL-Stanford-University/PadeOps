program test_cd10

    use kind_parameters, only: rkind
    use constants,       only: two,pi
    use timer,           only: tic,toc
    use cd10stuff,       only: cd10
    implicit none

    integer :: nx = 1024, ny=128, nz=128

    logical, parameter :: periodic = .TRUE.

    type( cd10 ) :: xcd10, ycd10, zcd10
    real(rkind), dimension(:,:,:), allocatable :: x,f,df,df_exact
    real(rkind) :: dx, dy, dz
    real(rkind), parameter :: omega = 1._rkind

    integer :: x1,xn,y1,yn,z1,zn
    integer :: i,ierr

    allocate( x(nx,ny,nz) )
    allocate( f(nx,ny,nz) )
    allocate( df(nx,ny,nz) )
    allocate( df_exact(nx,ny,nz) )

    if (periodic) then
        dx = two*pi/real(nx,rkind)
    else
        dx = two*pi/real(nx-1,rkind)
    end if
    dy = dx
    dz = dx

    do i=1,nx
        x(i,:,:) = real(i-1,rkind)*dx
        f(i,:,:) = sin(omega * x(i,:,:))
        df_exact(i,:,:) = omega * cos( omega * x(i,:,:))
    end do
    print*, "Created initial data"

    ierr = xcd10%init( nx, dx, periodic, 0, 0)
    ierr = ycd10%init( ny, dy, periodic, 0, 0)
    ierr = zcd10%init( nz, dz, periodic, 0, 0)

    call tic() 
    y1 = 1;yn = ny; z1 = 1; zn = nz
    call xcd10 % dd1(f,df,ny,nz,y1,yn,z1,zn)
    call toc ("Time to get the x derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - df_exact))

    call tic() 
    x1 = 1;xn = nx; z1 = 1; zn = nz
    call ycd10 % dd2(f,df,nx,nz,x1,xn,z1,zn)
    call toc ("Time to get the y derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df))

    call tic() 
    x1 = 1;xn = nx; y1 = 1; yn = ny
    call zcd10 % dd3(f,df,nx,ny,x1,xn,y1,yn)
    call toc ("Time to get the z derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df))

    call xcd10%destroy()
    call ycd10%destroy()
    call zcd10%destroy()

    deallocate( x )
    deallocate( f )
    deallocate( df )
    deallocate( df_exact )

end program
