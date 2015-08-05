program test_derivatives

    use kind_parameters, only: rkind
    use constants,       only: two,pi
    use timer,           only: tic,toc
    use derivativesWrapper, only: derivatives
    implicit none

    integer :: nx = 64, ny=64, nz=64


    type( derivatives ) :: method1, method2, method3
    real(rkind), dimension(:,:,:), allocatable :: x,y,z,f,df,dfdx_exact,dfdy_exact,dfdz_exact,d2fdx2_exact,d2fdy2_exact,d2fdz2_exact
    real(rkind) :: dx, dy, dz
    real(rkind), parameter :: omega = 10._rkind

    integer :: i,j,k

    allocate( x(nx,ny,nz) )
    allocate( y(nx,ny,nz) )
    allocate( z(nx,ny,nz) )
    allocate( f(nx,ny,nz) )
    allocate( df(nx,ny,nz) )
    allocate( dfdx_exact(nx,ny,nz) )
    allocate( dfdy_exact(nx,ny,nz) )
    allocate( dfdz_exact(nx,ny,nz) )
    allocate( d2fdx2_exact(nx,ny,nz) )
    allocate( d2fdy2_exact(nx,ny,nz) )
    allocate( d2fdz2_exact(nx,ny,nz) )

    dx = two*pi/real(nx,rkind)
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
        d2fdx2_exact(i,j,k) = -omega * omega*sin( omega * x(i,j,k))
        d2fdy2_exact(i,j,k) = -omega * omega*sin( omega * y(i,j,k)) 
        d2fdz2_exact(i,j,k) = -omega * omega*sin( omega * z(i,j,k)) 
    end do
    end do 
    end do 

    print*, "Created initial data"

    call method1%init(          nx,     ny,    nz, &
                                dx,     dy,    dz, &
                            .TRUE., .TRUE., .TRUE., &
                            "cd10", "cd10", "cd10" )

    call method2%init(          nx,     ny,    nz, &
                                dx,     dy,    dz, &
                            .TRUE., .TRUE.,  .TRUE., &
                            "cd06", "cd06", "cd06" )
    
    call method3%init(          nx,     ny,    nz, &
                                dx,     dy,    dz, &
                            .TRUE., .TRUE.,  .TRUE., &
                            "four", "four", "four" )

    print*, "Initialized all methods"
  
    print*, "FIRST DERIVATIVE TESTS" 
   
    print*, "==========================================="
    print*, "Now trying METHOD 1: CD10"
    print*, "==========================================="
    call tic() 
    call method1 % ddx(f,df)
    call toc ("Time to get the x derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdx_exact))

    call tic() 
    call method1 % ddy(f,df)
    call toc ("Time to get the y derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))

    call tic() 
    call method1 % ddz(f,df)
    call toc ("Time to get the z derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdz_exact))


    print*, "==========================================="
    print*, "Now trying METHOD 2: CD06"
    print*, "==========================================="
    call tic() 
    call method2 % ddx(f,df)
    call toc ("Time to get the x derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdx_exact))

    call tic() 
    call method2 % ddy(f,df)
    call toc ("Time to get the y derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))

    call tic() 
    call method2 % ddz(f,df)
    call toc ("Time to get the z derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdz_exact))


    print*, "==========================================="
    print*, "Now trying METHOD 3: FOUR"
    print*, "==========================================="
    call tic() 
    call method3 % ddx(f,df)
    call toc ("Time to get the x derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdx_exact))

    call tic() 
    call method3 % ddy(f,df)
    call toc ("Time to get the y derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))

    call tic() 
    call method3 % ddz(f,df)
    call toc ("Time to get the z derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdz_exact))



    print*, "---------------------------------------------"


    print*, "SECOND DERIVATIVE TESTS" 
    
    print*, "==========================================="
    print*, "Now trying METHOD 1: CD10"
    print*, "==========================================="
    call tic() 
    call method1 % d2dx2(f,df)
    call toc ("Time to get the x derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - d2fdx2_exact))
   

    call tic() 
    call method1 % d2dy2(f,df)
    call toc ("Time to get the y derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - d2fdy2_exact))

    call tic() 
    call method1 % d2dz2(f,df)
    call toc ("Time to get the z derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - d2fdz2_exact))


    !print*, "==========================================="
    !print*, "Now trying METHOD 2: CD06"
    !print*, "==========================================="
    !call tic() 
    !call method2 % d2dx2(f,df)
    !call toc ("Time to get the x derivative:")
    !print*, "Maximum error = ", MAXVAL( ABS(df - d2fdx2_exact))

    !call tic() 
    !call method2 % d2dy2(f,df)
    !call toc ("Time to get the y derivative:")
    !print*, "Maximum error = ", MAXVAL( ABS(df - d2fdy2_exact))

    !call tic() 
    !call method2 % d2dz2(f,df)
    !call toc ("Time to get the z derivative:")
    !print*, "Maximum error = ", MAXVAL( ABS(df - d2fdz2_exact))


    print*, "==========================================="
    print*, "Now trying METHOD 3: FOUR"
    print*, "==========================================="
    call tic() 
    call method3 % d2dx2(f,df)
    call toc ("Time to get the x derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - d2fdx2_exact))

    call tic() 
    call method3 % d2dy2(f,df)
    call toc ("Time to get the y derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - d2fdy2_exact))

    call tic() 
    call method3 % d2dz2(f,df)
    call toc ("Time to get the z derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - d2fdz2_exact))
    
    
    print*, "=========================================="
    print*, " Destroying everything"
    print*, "=========================================="

    deallocate( x )
    deallocate( y )
    deallocate( z )
    deallocate( f )
    deallocate( df )
    deallocate( dfdx_exact )
    deallocate( dfdy_exact )
    deallocate( dfdz_exact )
    deallocate( d2fdx2_exact )
    deallocate( d2fdy2_exact )
    deallocate( d2fdz2_exact )
    call method1%destroy
    call method2%destroy
    call method3%destroy
    

end program
