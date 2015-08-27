program test_filters

    use kind_parameters, only: rkind
    use constants,       only: two,pi
    use timer,           only: tic,toc
    use FiltersMod,      only: filters
    implicit none

    integer :: nx = 256, ny=256, nz=256


    type( filters ) :: method1, method2, method3
    real(rkind), dimension(:,:,:), allocatable :: x,y,z,f,df
    real(rkind) :: dx, dy, dz
    real(rkind), parameter :: omega = 1._rkind

    integer :: i,j,k

    allocate( x(nx,ny,nz) )
    allocate( y(nx,ny,nz) )
    allocate( z(nx,ny,nz) )
    allocate( f(nx,ny,nz) )
    allocate( df(nx,ny,nz) )

    dx = two*pi/real(nx,rkind)
    dy = two*pi/real(ny,rkind)
    dz = two*pi/real(nz,rkind)

    do k=1,nz
        do j=1,ny
            do i=1,nx
                x(i,j,k) = real(i-1,rkind)*dx
                y(i,j,k) = real(j-1,rkind)*dy
                z(i,j,k) = real(k-1,rkind)*dz
                f(i,j,k) = cos(omega * x(i,j,k)) + cos(omega * y(i,j,k)) + cos(omega * z(i,j,k))
            end do
        end do 
    end do 

    print*, "Created initial data"

    call method1%init(          nx,     ny,    nz, &
                            .TRUE., .TRUE., .TRUE., &
                            "cf90", "cf90", "cf90" )

    call method2%init(          nx,     ny,    nz, &
                            .TRUE., .TRUE.,  .TRUE., &
                            "gaussian", "gaussian", "gaussian" )
    
    call method3%init(          nx,     ny,    nz, &
                            .TRUE., .TRUE.,  .TRUE., &
                            "lstsq", "lstsq", "lstsq" )

    print*, "Initialized all methods"
   
   
    print*, "==========================================="
    print*, "Now trying METHOD 1: CF90"
    print*, "==========================================="
    call tic() 
    call method1 % filterx(f,df)
    call toc ("Time to get the x filter:")

    call tic() 
    call method1 % filtery(f,df)
    call toc ("Time to get the y filter:")

    call tic() 
    call method1 % filterz(f,df)
    call toc ("Time to get the z filter:")


    print*, "==========================================="
    print*, "Now trying METHOD 2: GAUSSIAN"
    print*, "==========================================="
    call tic() 
    call method2 % filterx(f,df)
    call toc ("Time to get the x filter:")

    call tic() 
    call method2 % filtery(f,df)
    call toc ("Time to get the y filter:")

    call tic() 
    call method2 % filterz(f,df)
    call toc ("Time to get the z filter:")


    print*, "==========================================="
    print*, "Now trying METHOD 3: LSTSQ"
    print*, "==========================================="
    call tic() 
    call method3 % filterx(f,df)
    call toc ("Time to get the x filter:")

    call tic() 
    call method3 % filtery(f,df)
    call toc ("Time to get the y filter:")

    call tic() 
    call method3 % filterz(f,df)
    call toc ("Time to get the z filter:")

    print*, "=========================================="
    print*, "Destroying everything"
    print*, "=========================================="

    call method1%destroy
    call method2%destroy
    call method3%destroy
    deallocate( x )
    deallocate( y )
    deallocate( z )
    deallocate( f )
    deallocate( df )
    

end program
