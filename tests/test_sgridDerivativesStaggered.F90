program test_derivatives_staggered

    use kind_parameters,          only: rkind
    use constants,                only: half,two,pi
    use DerivativesStaggeredMod,  only: derivativesStagg
    implicit none

    integer :: nx, ny, nz
    integer :: n_vec(4) 
    real(rkind) :: error
    real(rkind) :: cd10_N2F_error(4,3), cd10_F2N_error(4,3)
    real(rkind) :: cd06_N2F_error(4,3), cd06_F2N_error(4,3)
    real(rkind) ::  d02_N2F_error(4,3),  d02_F2N_error(4,3)
    real(rkind) ::  d04_N2F_error(4,3),  d04_F2N_error(4,3)
    real(rkind) ::  d06_N2F_error(4,3),  d06_F2N_error(4,3)
    real(rkind) :: ratio(3)

    type( derivativesStagg ) :: method1, method2, method3, method4, method5
    real(rkind), dimension(:,:,:), allocatable :: x, y, z, f, dfdx_exact, dfdy_exact, dfdz_exact, df
    real(rkind), dimension(:,:,:), allocatable :: xF,yF,zF,fF,dfdxF_exact,dfdyF_exact,dfdzF_exact
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
        allocate( xF(nx,ny,nz) )
        allocate( yF(nx,ny,nz) )
        allocate( zF(nx,ny,nz) )
        allocate( f(nx,ny,nz) )
        allocate( fF(nx,ny,nz) )
        allocate( dfdx_exact(nx,ny,nz) )
        allocate( dfdy_exact(nx,ny,nz) )
        allocate( dfdz_exact(nx,ny,nz) )
        allocate( dfdxF_exact(nx,ny,nz) )
        allocate( dfdyF_exact(nx,ny,nz) )
        allocate( dfdzF_exact(nx,ny,nz) )

        dx = two*pi/real(nx,rkind)
        dy = two*pi/real(ny,rkind)
        dz = two*pi/real(nz,rkind)

        do k=1,nz
        do j=1,ny
        do i=1,nx
            x(i,j,k) = real(i-1,rkind)*dx
            y(i,j,k) = real(j-1,rkind)*dy
            z(i,j,k) = real(k-1,rkind)*dz
            xF(i,j,k) = MOD( x(i,j,k) + dx*half , two*pi ) 
            yF(i,j,k) = MOD( y(i,j,k) + dy*half , two*pi )
            zF(i,j,k) = MOD( z(i,j,k) + dz*half , two*pi )
            f(i,j,k)  = sin(omega *  x(i,j,k)) + sin(omega *  y(i,j,k)) + sin(omega *  z(i,j,k))
            fF(i,j,k) = sin(omega * xF(i,j,k)) + sin(omega * yF(i,j,k)) + sin(omega * zF(i,j,k))
            dfdx_exact(i,j,k) = omega * cos( omega * x(i,j,k))
            dfdy_exact(i,j,k) = omega * cos( omega * y(i,j,k)) 
            dfdz_exact(i,j,k) = omega * cos( omega * z(i,j,k)) 
            dfdxF_exact(i,j,k) = omega * cos( omega * xF(i,j,k))
            dfdyF_exact(i,j,k) = omega * cos( omega * yF(i,j,k)) 
            dfdzF_exact(i,j,k) = omega * cos( omega * zF(i,j,k)) 
        end do
        end do 
        end do 

        deallocate( x )
        deallocate( y )
        deallocate( z )
        deallocate( xF )
        deallocate( yF )
        deallocate( zF )
        
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
  
        print*, "FIRST DERIVATIVE TESTS: Node 2 Face" 
   
        print*, "==========================================="
        print*, "Now trying METHOD 1: CD10"
        print*, "==========================================="
        call method1 % ddxN2F(f,df)
        error = MAXVAL( ABS(df - dfdxF_exact))
        print*, "Maximum error = ", error
        cd10_N2F_error(ind,1) = error

        call method1 % ddyN2F(f,df)
        error = MAXVAL( ABS(df - dfdyF_exact))
        print*, "Maximum error = ", error
        cd10_N2F_error(ind,2) = error

        call method1 % ddzN2F(f,df)
        error = MAXVAL( ABS(df - dfdzF_exact))
        print*, "Maximum error = ", error
        cd10_N2F_error(ind,3) = error


        print*, "==========================================="
        print*, "Now trying METHOD 2: CD06"
        print*, "==========================================="
        call method2 % ddxN2F(f,df)
        error = MAXVAL( ABS(df - dfdxF_exact))
        print*, "Maximum error = ", error
        cd06_N2F_error(ind,1) = error

        call method2 % ddyN2F(f,df)
        error = MAXVAL( ABS(df - dfdyF_exact))
        print*, "Maximum error = ", error
        cd06_N2F_error(ind,2) = error

        call method2 % ddzN2F(f,df)
        error = MAXVAL( ABS(df - dfdzF_exact))
        print*, "Maximum error = ", error
        cd06_N2F_error(ind,3) = error


        print*, "==========================================="
        print*, "Now trying METHOD 3: D02"
        print*, "==========================================="
        call method3 % ddxN2F(f,df)
        error = MAXVAL( ABS(df - dfdxF_exact))
        print*, "Maximum error = ", error
        d02_N2F_error(ind,1) = error

        call method3 % ddyN2F(f,df)
        error = MAXVAL( ABS(df - dfdyF_exact))
        print*, "Maximum error = ", error
        d02_N2F_error(ind,2) = error

        call method3 % ddzN2F(f,df)
        error = MAXVAL( ABS(df - dfdzF_exact))
        print*, "Maximum error = ", error
        d02_N2F_error(ind,3) = error

        print*, "==========================================="
        print*, "Now trying METHOD 4: D04"
        print*, "==========================================="
        call method4 % ddxN2F(f,df)
        error = MAXVAL( ABS(df - dfdxF_exact))
        print*, "Maximum error = ", error
        d04_N2F_error(ind,1) = error

        call method4 % ddyN2F(f,df)
        error = MAXVAL( ABS(df - dfdyF_exact))
        print*, "Maximum error = ", error
        d04_N2F_error(ind,2) = error

        call method4 % ddzN2F(f,df)
        error = MAXVAL( ABS(df - dfdzF_exact))
        print*, "Maximum error = ", error
        d04_N2F_error(ind,3) = error

        print*, "==========================================="
        print*, "Now trying METHOD 5: D06"
        print*, "==========================================="
        call method5 % ddxN2F(f,df)
        error = MAXVAL( ABS(df - dfdxF_exact))
        print*, "Maximum error = ", error
        d06_N2F_error(ind,1) = error

        call method5 % ddyN2F(f,df)
        error = MAXVAL( ABS(df - dfdyF_exact))
        print*, "Maximum error = ", error
        d06_N2F_error(ind,2) = error

        call method5 % ddzN2F(f,df)
        error = MAXVAL( ABS(df - dfdzF_exact))
        print*, "Maximum error = ", error
        d06_N2F_error(ind,3) = error


        print*, "---------------------------------------------"
        print*, "---------------------------------------------"

        print*, "FIRST DERIVATIVE TESTS: Face 2 Node" 
   
        print*, "==========================================="
        print*, "Now trying METHOD 1: CD10"
        print*, "==========================================="
        call method1 % ddxF2N(fF,df)
        error = MAXVAL( ABS(df - dfdx_exact))
        print*, "Maximum error = ", error
        cd10_F2N_error(ind,1) = error

        call method1 % ddyF2N(fF,df)
        error = MAXVAL( ABS(df - dfdy_exact))
        print*, "Maximum error = ", error
        cd10_F2N_error(ind,2) = error

        call method1 % ddzF2N(fF,df)
        error = MAXVAL( ABS(df - dfdz_exact))
        print*, "Maximum error = ", error
        cd10_F2N_error(ind,3) = error


        print*, "==========================================="
        print*, "Now trying METHOD 2: CD06"
        print*, "==========================================="
        call method2 % ddxF2N(fF,df)
        error = MAXVAL( ABS(df - dfdx_exact))
        print*, "Maximum error = ", error
        cd06_F2N_error(ind,1) = error

        call method2 % ddyF2N(fF,df)
        error = MAXVAL( ABS(df - dfdy_exact))
        print*, "Maximum error = ", error
        cd06_F2N_error(ind,2) = error

        call method2 % ddzF2N(fF,df)
        error = MAXVAL( ABS(df - dfdz_exact))
        print*, "Maximum error = ", error
        cd06_F2N_error(ind,3) = error


        print*, "==========================================="
        print*, "Now trying METHOD 3: D02"
        print*, "==========================================="
        call method3 % ddxF2N(fF,df)
        error = MAXVAL( ABS(df - dfdx_exact))
        print*, "Maximum error = ", error
        d02_F2N_error(ind,1) = error

        call method3 % ddyF2N(fF,df)
        error = MAXVAL( ABS(df - dfdy_exact))
        print*, "Maximum error = ", error
        d02_F2N_error(ind,2) = error

        call method3 % ddzF2N(fF,df)
        error = MAXVAL( ABS(df - dfdz_exact))
        print*, "Maximum error = ", error
        d02_F2N_error(ind,3) = error

        print*, "==========================================="
        print*, "Now trying METHOD 4: D04"
        print*, "==========================================="
        call method4 % ddxF2N(fF,df)
        error = MAXVAL( ABS(df - dfdx_exact))
        print*, "Maximum error = ", error
        d04_F2N_error(ind,1) = error

        call method4 % ddyF2N(fF,df)
        error = MAXVAL( ABS(df - dfdy_exact))
        print*, "Maximum error = ", error
        d04_F2N_error(ind,2) = error

        call method4 % ddzF2N(fF,df)
        error = MAXVAL( ABS(df - dfdz_exact))
        print*, "Maximum error = ", error
        d04_F2N_error(ind,3) = error

        print*, "==========================================="
        print*, "Now trying METHOD 5: D06"
        print*, "==========================================="
        call method5 % ddxF2N(fF,df)
        error = MAXVAL( ABS(df - dfdx_exact))
        print*, "Maximum error = ", error
        d06_F2N_error(ind,1) = error

        call method5 % ddyF2N(fF,df)
        error = MAXVAL( ABS(df - dfdy_exact))
        print*, "Maximum error = ", error
        d06_F2N_error(ind,2) = error

        call method5 % ddzF2N(fF,df)
        error = MAXVAL( ABS(df - dfdz_exact))
        print*, "Maximum error = ", error
        d06_F2N_error(ind,3) = error


        print*, "---------------------------------------------"
        print*, "---------------------------------------------"
        

        
        print*, "=========================================="
        print*, " Destroying everything"
        print*, "=========================================="

        deallocate( f )
        deallocate( fF )
        deallocate( df )
        deallocate( dfdx_exact )
        deallocate( dfdy_exact )
        deallocate( dfdz_exact )
        deallocate(dfdxF_exact )
        deallocate(dfdyF_exact )
        deallocate(dfdzF_exact )
        call method1%destroy
        call method2%destroy
        call method3%destroy
        call method4%destroy
        call method5%destroy

    enddo

    print*, "=========================================="
    print*, "====1st Derivative Results: N2F =========="
    print*, "=========================================="
    ratio = cd10_N2F_error(1:3,1)/cd10_N2F_error(2:4,1) 
    print *, "CD10 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = cd10_N2F_error(1:3,2)/cd10_N2F_error(2:4,2) 
    print *, "CD10 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = cd10_N2F_error(1:3,3)/cd10_N2F_error(2:4,3) 
    print *, "CD10 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="
    
    print*, "=========================================="
    ratio = cd06_N2F_error(1:3,1)/cd06_N2F_error(2:4,1) 
    print *, "CD06 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = cd06_N2F_error(1:3,2)/cd06_N2F_error(2:4,2) 
    print *, "CD06 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = cd06_N2F_error(1:3,3)/cd06_N2F_error(2:4,3) 
    print *, "CD06 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="

    print*, "=========================================="
    ratio = d02_N2F_error(1:3,1)/d02_N2F_error(2:4,1) 
    print *, "D02 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = d02_N2F_error(1:3,2)/d02_N2F_error(2:4,2) 
    print *, "D02 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = d02_N2F_error(1:3,3)/d02_N2F_error(2:4,3) 
    print *, "D02 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="

    print*, "=========================================="
    ratio = d04_N2F_error(1:3,1)/d04_N2F_error(2:4,1) 
    print *, "D04 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = d04_N2F_error(1:3,2)/d04_N2F_error(2:4,2) 
    print *, "D04 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = d04_N2F_error(1:3,3)/d04_N2F_error(2:4,3) 
    print *, "D04 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="

    print*, "=========================================="
    ratio = d06_N2F_error(1:3,1)/d06_N2F_error(2:4,1) 
    print *, "D06 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = d06_N2F_error(1:3,2)/d06_N2F_error(2:4,2) 
    print *, "D06 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = d06_N2F_error(1:3,3)/d06_N2F_error(2:4,3) 
    print *, "D06 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="

    print*, "=========================================="
    print*, "====1st Derivative Results: F2N =========="
    print*, "=========================================="
    ratio = cd10_F2N_error(1:3,1)/cd10_F2N_error(2:4,1) 
    print *, "CD10 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = cd10_F2N_error(1:3,2)/cd10_F2N_error(2:4,2) 
    print *, "CD10 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = cd10_F2N_error(1:3,3)/cd10_F2N_error(2:4,3) 
    print *, "CD10 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="
    
    print*, "=========================================="
    ratio = cd06_F2N_error(1:3,1)/cd06_F2N_error(2:4,1) 
    print *, "CD06 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = cd06_F2N_error(1:3,2)/cd06_F2N_error(2:4,2) 
    print *, "CD06 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = cd06_F2N_error(1:3,3)/cd06_F2N_error(2:4,3) 
    print *, "CD06 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="

    print*, "=========================================="
    ratio = d02_F2N_error(1:3,1)/d02_F2N_error(2:4,1) 
    print *, "D02 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = d02_F2N_error(1:3,2)/d02_F2N_error(2:4,2) 
    print *, "D02 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = d02_F2N_error(1:3,3)/d02_F2N_error(2:4,3) 
    print *, "D02 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="

    print*, "=========================================="
    ratio = d04_F2N_error(1:3,1)/d04_F2N_error(2:4,1) 
    print *, "D04 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = d04_F2N_error(1:3,2)/d04_F2N_error(2:4,2) 
    print *, "D04 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = d04_F2N_error(1:3,3)/d04_F2N_error(2:4,3) 
    print *, "D04 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="

    print*, "=========================================="
    ratio = d06_F2N_error(1:3,1)/d06_F2N_error(2:4,1) 
    print *, "D06 1st Derivative Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = d06_F2N_error(1:3,2)/d06_F2N_error(2:4,2) 
    print *, "D06 1st Derivative Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = d06_F2N_error(1:3,3)/d06_F2N_error(2:4,3) 
    print *, "D06 1st Derivative Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "=========================================="

end program
