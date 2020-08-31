program test_derivatives

    use kind_parameters, only: rkind
    use constants,       only: two,pi
    use DerivativesMod,  only: derivatives
    implicit none

    integer :: nx, ny, nz
    integer :: n_vec(4) 
    real(rkind) :: error
    real(rkind) :: cd10_1_error(4,3), cd10_2_error(4,3)
    real(rkind) :: cd06_1_error(4,3), cd06_2_error(4,3)
    real(rkind) :: d02_1_error(4,3), d02_2_error(4,3)
    real(rkind) :: d04_1_error(4,3), d04_2_error(4,3)
    real(rkind) :: d06_1_error(4,3), d06_2_error(4,3)
    real(rkind) :: ratio(3)

    type( derivatives ) :: method1, method2, method3, method4, method5
    real(rkind), dimension(:,:,:), allocatable :: x,y,z,f,df,dfdx_exact,dfdy_exact,dfdz_exact,d2fdx2_exact,d2fdy2_exact,d2fdz2_exact
    real(rkind) :: dx, dy, dz
    real(rkind), parameter :: omega = 2._rkind

    integer :: i,j,k, ind
    
    n_vec = (/16, 32, 64, 128/)
    !n_vec = (/256, 512, 1024, 2048/)
    !n_vec = (/512, 1024, 2048, 4096/)

    do ind=1,4

        print *, "==========================================="
        print *, "==========================================="
        print *, "Testing with nx/ny/nz =", n_vec(ind)
        print *, "==========================================="
        print *, "==========================================="

        nx = n_vec(ind)
        ny = n_vec(ind)
        nz = n_vec(ind)

        allocate( x(nx,ny,nz) )
        allocate( y(nx,ny,nz) )
        allocate( z(nx,ny,nz) )
        allocate( f(nx,ny,nz) )
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

        deallocate( x )
        deallocate( y )
        deallocate( z )
        
        allocate( df(nx,ny,nz) )

        print*, "Created initial data"


        !Note: periodicity is controlled by true/false statements
        call method1%init(          nx,     ny,    nz, &
                                    dx,     dy,    dz, &
                                .true., .true., .true., &
                                "cd10", "cd10", "cd10" )

        call method2%init(          nx,     ny,    nz, &
                                    dx,     dy,    dz, &
                                .true., .true., .true., &
                                "cd06", "cd06", "cd06" )
        
        call method3%init(          nx,     ny,    nz, &
                                    dx,     dy,    dz, &
                                .true., .true., .true., &
                                   "d02", "d02", "d02" )
         
        call method4%init(          nx,     ny,    nz, &
                                    dx,     dy,    dz, &
                                .true., .true., .true., &
                                   "d04", "d04", "d04" )
         
        call method5%init(          nx,     ny,    nz, &
                                    dx,     dy,    dz, &
                                .true., .true., .true., &
                                   "d06", "d06", "d06" )
         

        print*, "Initialized all methods"
  
        print*, "FIRST DERIVATIVE TESTS" 
   
        print*, "==========================================="
        print*, "Now trying METHOD 1: CD10"
        print*, "==========================================="
        call method1 % ddx(f,df)
        error = MAXVAL( ABS(df - dfdx_exact))
        print*, "Maximum error = ", error
        cd10_1_error(ind,1) = error

        call method1 % ddy(f,df)
        error = MAXVAL( ABS(df - dfdy_exact))
        print*, "Maximum error = ", error
        cd10_1_error(ind,2) = error

        call method1 % ddz(f,df)
        error = MAXVAL( ABS(df - dfdz_exact))
        print*, "Maximum error = ", error
        cd10_1_error(ind,3) = error


        print*, "==========================================="
        print*, "Now trying METHOD 2: CD06"
        print*, "==========================================="
        call method2 % ddx(f,df)
        error = MAXVAL( ABS(df - dfdx_exact))
        print*, "Maximum error = ", error
        cd06_1_error(ind,1) = error

        call method2 % ddy(f,df)
        error = MAXVAL( ABS(df - dfdy_exact))
        print*, "Maximum error = ", error
        cd06_1_error(ind,2) = error

        call method2 % ddz(f,df)
        error = MAXVAL( ABS(df - dfdz_exact))
        print*, "Maximum error = ", error
        cd06_1_error(ind,3) = error


        print*, "==========================================="
        print*, "Now trying METHOD 3: D02"
        print*, "==========================================="
        call method3 % ddx(f,df)
        error = MAXVAL( ABS(df - dfdx_exact))
        print*, "Maximum error = ", error
        d02_1_error(ind,1) = error

        call method3 % ddy(f,df)
        error = MAXVAL( ABS(df - dfdy_exact))
        print*, "Maximum error = ", error
        d02_1_error(ind,2) = error

        call method3 % ddz(f,df)
        error = MAXVAL( ABS(df - dfdz_exact))
        print*, "Maximum error = ", error
        d02_1_error(ind,3) = error

        print*, "==========================================="
        print*, "Now trying METHOD 4: D04"
        print*, "==========================================="
        call method4 % ddx(f,df)
        error = MAXVAL( ABS(df - dfdx_exact))
        print*, "Maximum error = ", error
        d04_1_error(ind,1) = error

        call method4 % ddy(f,df)
        error = MAXVAL( ABS(df - dfdy_exact))
        print*, "Maximum error = ", error
        d04_1_error(ind,2) = error

        call method4 % ddz(f,df)
        error = MAXVAL( ABS(df - dfdz_exact))
        print*, "Maximum error = ", error
        d04_1_error(ind,3) = error

        print*, "==========================================="
        print*, "Now trying METHOD 5: D06"
        print*, "==========================================="
        call method5 % ddx(f,df)
        error = MAXVAL( ABS(df - dfdx_exact))
        print*, "Maximum error = ", error
        d06_1_error(ind,1) = error

        call method5 % ddy(f,df)
        error = MAXVAL( ABS(df - dfdy_exact))
        print*, "Maximum error = ", error
        d06_1_error(ind,2) = error

        call method5 % ddz(f,df)
        error = MAXVAL( ABS(df - dfdz_exact))
        print*, "Maximum error = ", error
        d06_1_error(ind,3) = error


        print*, "---------------------------------------------"


        print*, "SECOND DERIVATIVE TESTS" 
        deallocate( dfdx_exact )
        deallocate( dfdy_exact )
        deallocate( dfdz_exact )

        print*, "==========================================="
        print*, "Now trying METHOD 1: CD10"
        print*, "==========================================="
        call method1 % d2dx2(f,df)
        error = MAXVAL( ABS(df - d2fdx2_exact))
        print*, "Maximum error = ", error
        cd10_2_error(ind,1) = error

        call method1 % d2dy2(f,df)
        error = MAXVAL( ABS(df - d2fdy2_exact))
        print*, "Maximum error = ", error
        cd10_2_error(ind,2) = error

        call method1 % d2dz2(f,df)
        error = MAXVAL( ABS(df - d2fdz2_exact))
        print*, "Maximum error = ", error
        cd10_2_error(ind,3) = error


        !print*, "==========================================="
        !print*, "Now trying METHOD 2: CD06"
        !print*, "==========================================="
        !call method2 % d2dx2(f,df)
        !error = MAXVAL( ABS(df - d2fdx2_exact))
        !print*, "Maximum error = ", error
        !cd06_2_error(ind,1) = error

        !call method2 % d2dy2(f,df)
        !error = MAXVAL( ABS(df - d2fdy2_exact))
        !print*, "Maximum error = ", error
        !cd06_2_error(ind,2) = error

        !call method2 % d2dz2(f,df)
        !error = MAXVAL( ABS(df - d2fdz2_exact))
        !print*, "Maximum error = ", error
        !cd06_2_error(ind,3) = error


        print*, "==========================================="
        print*, "Now trying METHOD 3: D02"
        print*, "==========================================="
        call method3 % d2dx2(f,df)
        error = MAXVAL( ABS(df - d2fdx2_exact))
        print*, "Maximum error = ", error
        d02_2_error(ind,1) = error

        call method3 % d2dy2(f,df)
        error = MAXVAL( ABS(df - d2fdy2_exact))
        print*, "Maximum error = ", error
        d02_2_error(ind,2) = error

        call method3 % d2dz2(f,df)
        error = MAXVAL( ABS(df - d2fdz2_exact))
        print*, "Maximum error = ", error
        d02_2_error(ind,3) = error

        print*, "==========================================="
        print*, "Now trying METHOD 4: D04"
        print*, "==========================================="
        call method4 % d2dx2(f,df)
        error = MAXVAL( ABS(df - d2fdx2_exact))
        print*, "Maximum error = ", error
        d04_2_error(ind,1) = error

        call method4 % d2dy2(f,df)
        error = MAXVAL( ABS(df - d2fdy2_exact))
        print*, "Maximum error = ", error
        d04_2_error(ind,2) = error

        call method4 % d2dz2(f,df)
        error = MAXVAL( ABS(df - d2fdz2_exact))
        print*, "Maximum error = ", error
        d04_2_error(ind,3) = error

        print*, "==========================================="
        print*, "Now trying METHOD 5: D06"
        print*, "==========================================="
        call method5 % d2dx2(f,df)
        error = MAXVAL( ABS(df - d2fdx2_exact))
        print*, "Maximum error = ", error
        d06_2_error(ind,1) = error

        call method5 % d2dy2(f,df)
        error = MAXVAL( ABS(df - d2fdy2_exact))
        print*, "Maximum error = ", error
        d06_2_error(ind,2) = error

        call method5 % d2dz2(f,df)
        error = MAXVAL( ABS(df - d2fdz2_exact))
        print*, "Maximum error = ", error
        d06_2_error(ind,3) = error


        print*, "---------------------------------------------"
        

        
        print*, "=========================================="
        print*, " Destroying everything"
        print*, "=========================================="

        deallocate( f )
        deallocate( df )
        deallocate( d2fdx2_exact )
        deallocate( d2fdy2_exact )
        deallocate( d2fdz2_exact )
        call method1%destroy
        call method2%destroy
        call method3%destroy
        call method4%destroy
        call method5%destroy

    enddo

    print*, "=========================================="
    print*, "====1st Derivative Results================"
    print*, "=========================================="
    ratio = cd10_1_error(1:3,1)/cd10_1_error(2:4,1) 
    print *, "CD10 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = cd10_1_error(1:3,2)/cd10_1_error(2:4,2) 
    print *, "CD10 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = cd10_1_error(1:3,3)/cd10_1_error(2:4,3) 
    print *, "CD10 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="
    
    print*, "=========================================="
    ratio = cd06_1_error(1:3,1)/cd06_1_error(2:4,1) 
    print *, "CD06 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = cd06_1_error(1:3,2)/cd06_1_error(2:4,2) 
    print *, "CD06 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = cd06_1_error(1:3,3)/cd06_1_error(2:4,3) 
    print *, "CD06 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="

    print*, "=========================================="
    ratio = d02_1_error(1:3,1)/d02_1_error(2:4,1) 
    print *, "D02 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = d02_1_error(1:3,2)/d02_1_error(2:4,2) 
    print *, "D02 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = d02_1_error(1:3,3)/d02_1_error(2:4,3) 
    print *, "D02 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="

    print*, "=========================================="
    ratio = d04_1_error(1:3,1)/d04_1_error(2:4,1) 
    print *, "D04 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = d04_1_error(1:3,2)/d04_1_error(2:4,2) 
    print *, "D04 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = d04_1_error(1:3,3)/d04_1_error(2:4,3) 
    print *, "D04 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="

    print*, "=========================================="
    ratio = d06_1_error(1:3,1)/d06_1_error(2:4,1) 
    print *, "D06 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = d06_1_error(1:3,2)/d06_1_error(2:4,2) 
    print *, "D06 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = d06_1_error(1:3,3)/d06_1_error(2:4,3) 
    print *, "D06 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="

    print*, "=========================================="
    print*, "====2nd Derivative Results================"
    print*, "=========================================="
    ratio = cd10_2_error(1:3,1)/cd10_2_error(2:4,1) 
    print *, "CD10 2nd Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = cd10_2_error(1:3,2)/cd10_2_error(2:4,2) 
    print *, "CD10 2nd Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = cd10_2_error(1:3,3)/cd10_2_error(2:4,3) 
    print *, "CD10 2nd Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="
    
    !print*, "=========================================="
    !ratio = cd06_2_error(1:3,1)/cd06_2_error(2:4,1) 
    !print *, "CD06 2nd Derivative Order of Convergence (x)"
    !print *,  log(ratio) / log(2.0)
    !ratio = cd06_2_error(1:3,2)/cd06_2_error(2:4,2) 
    !print *, "CD06 2nd Derivative Order of Convergence (y)"
    !print *,  log(ratio) / log(2.0)
    !ratio = cd06_2_error(1:3,3)/cd06_2_error(2:4,3) 
    !print *, "CD06 2nd Derivative Order of Convergence (z)"
    !print *,  log(ratio) / log(2.0)
    !print*, "=========================================="

    print*, "=========================================="
    ratio = d02_2_error(1:3,1)/d02_2_error(2:4,1) 
    print *, "D02 2nd Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = d02_2_error(1:3,2)/d02_2_error(2:4,2) 
    print *, "D02 2nd Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = d02_2_error(1:3,3)/d02_2_error(2:4,3) 
    print *, "D02 2nd Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="

    print*, "=========================================="
    ratio = d04_2_error(1:3,1)/d04_2_error(2:4,1) 
    print *, "D04 2nd Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = d04_2_error(1:3,2)/d04_2_error(2:4,2) 
    print *, "D04 2nd Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = d04_2_error(1:3,3)/d04_2_error(2:4,3) 
    print *, "D04 2nd Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="

    print*, "=========================================="
    ratio = d06_2_error(1:3,1)/d06_2_error(2:4,1) 
    print *, "D06 2nd Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = d06_2_error(1:3,2)/d06_2_error(2:4,2) 
    print *, "D06 2nd Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = d06_2_error(1:3,3)/d06_2_error(2:4,3) 
    print *, "D06 2nd Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="
end program
