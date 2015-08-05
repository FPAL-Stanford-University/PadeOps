program test_cd10

    use kind_parameters, only: rkind
    use constants,       only: two,pi
    use timer,           only: tic,toc
    use cd10stuff,       only: cd10
    implicit none

    integer :: nx = 32, ny=32, nz=32

    logical, parameter :: periodic = .TRUE.

    type( cd10 ) :: xcd10, ycd10, zcd10
    real(rkind), dimension(:,:,:), allocatable :: x,y,z,f,df,df_exact,d2f,d2f_exact
    real(rkind) :: dx, dy, dz
    real(rkind), parameter :: omega = 1._rkind

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

    f = sin(omega * x) + sin(omega * y) + sin(omega * z)
    ierr = xcd10%init( nx, dx, periodic, 0, 0)
    ierr = ycd10%init( ny, dy, periodic, 0, 0)
    ierr = zcd10%init( nz, dz, periodic, 0, 0)

    call tic() 
    call xcd10 % dd1(f,df,ny,nz)
    call xcd10 % d2d1(f,d2f,ny,nz)
    call toc ("Time to get the x derivative:")
    df_exact = omega* cos(omega * x) 
    df_exact = -omega* omega* sin(omega * x) 
    print*, "Maximum error (first  der) = ", MAXVAL( ABS(df - df_exact))
    print*, "Maximum error (second der) = ", MAXVAL( ABS(d2f - d2f_exact))

    call tic() 
    call ycd10 % dd2(f,df,nx,nz)
    call ycd10 % d2d2(f,d2f,nx,nz)
    call toc ("Time to get the y derivative:")
    df_exact = omega* cos(omega * y) 
    df_exact = -omega* omega* sin(omega * y) 
    print*, "Maximum error (first  der) = ", MAXVAL( ABS(df - df_exact))
    print*, "Maximum error (second der) = ", MAXVAL( ABS(d2f - d2f_exact))

    call tic() 
    call zcd10 % dd3(f,df,nx,ny)
    call zcd10 % d2d3(f,d2f,nx,ny)
    df_exact = omega* cos(omega * z) 
    df_exact = -omega* omega* sin(omega * z) 
    print*, "Maximum error (first  der) = ", MAXVAL( ABS(df - df_exact))
    call toc ("Time to get the z derivative:")
    print*, "Maximum error (first  der) = ", MAXVAL( ABS(df))
    print*, "Maximum error (second der) = ", MAXVAL( ABS(d2f - d2f_exact))

    call xcd10%destroy()
    call ycd10%destroy()
    call zcd10%destroy()

    deallocate( x )
    deallocate( f )
    deallocate( df )
    deallocate( df_exact )

end program
