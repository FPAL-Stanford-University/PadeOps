program test_interpolators

    use kind_parameters, only: rkind, clen
    use constants,       only: two, pi, half
    use InterpolatorsMod,  only: interpolators
    implicit none

    integer :: nx, ny, nz
    integer :: n_vec(4) 
    real(rkind) :: error
    real(rkind) :: ei02_N2F_error(4,3), ei02_F2N_error(4,3)
    real(rkind) :: ei06_N2F_error(4,3), ei06_F2N_error(4,3)
    real(rkind) :: ratio(3)

    type( interpolators ) ::method3, method5
    real(rkind), dimension(:,:,:), allocatable :: x,y,z,f,fF, fN
    real(rkind), dimension(:,:,:), allocatable :: xF,yF,zF, fFx_exact, fFy_exact, fFz_exact
    real(rkind) :: dx, dy, dz, Lx, Ly, Lz
    real(rkind), parameter :: omega = 2._rkind, phi = 1.73d0

    integer :: i,j,k, ind
    logical :: periodic_x, periodic_y, periodic_z
    integer :: bc1_x, bcn_x, bc1_y, bcn_y, bc1_z, bcn_z
    character(len=clen) :: inputfile
    integer :: ierr, ioUnit 

    ! Start MPI
    call MPI_Init(ierr)

    ! Get file location 
    call GETARG(1,inputfile)

    namelist /TEST/ bc1_x, bcn_x, bc1_y, bcn_y, bc1_z, bcn_z,  &
                    periodic_x, periodic_y, periodic_z

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=TEST)

    !n_vec = (/16, 32, 64, 128/)
    n_vec = (/32, 64, 128, 256/)
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
        allocate( fFx_exact(nx,ny,nz) )
        allocate( fFy_exact(nx,ny,nz) )
        allocate( fFz_exact(nx,ny,nz) )

        Lx = two * pi
        Ly = two * pi
        Lz = two * pi

        !Grid spacing and placement depends on periodicity
        if (periodic_x) then
            dx = Lx/real(nx,rkind)
        else
            dx = Lx/real(nx-1,rkind)
        endif
        if (periodic_y) then
            dy = Ly/real(ny,rkind)
        else
            dy = Ly/real(ny-1,rkind)
        endif
        if (periodic_z) then
            dz = Lz/real(nz,rkind)
        else
            dz = Lz/real(nz-1,rkind)
        endif

        do k=1,nz
        do j=1,ny
        do i=1,nx
            x(i,j,k) = real(i-1,rkind)*dx
            y(i,j,k) = real(j-1,rkind)*dy
            z(i,j,k) = real(k-1,rkind)*dz
            xF(i,j,k) = MOD( x(i,j,k) + dx*half , two*pi ) 
            yF(i,j,k) = MOD( y(i,j,k) + dy*half , two*pi )
            zF(i,j,k) = MOD( z(i,j,k) + dz*half , two*pi )
            f(i,j,k) = sin(omega * x(i,j,k)) + sin(omega * y(i,j,k)) + sin(omega * z(i,j,k))
            fFx_exact(i,j,k) = sin(omega * xF(i,j,k)) + sin(omega * y(i,j,k))  + sin(omega * z(i,j,k))
            fFy_exact(i,j,k) = sin(omega * x(i,j,k))  + sin(omega * yF(i,j,k)) + sin(omega * z(i,j,k))
            fFz_exact(i,j,k) = sin(omega * x(i,j,k))  + sin(omega * y(i,j,k))  + sin(omega * zF(i,j,k))
            f(i,j,k)         = 1.5d0 + sin(omega *  x(i,j,k)  + 1.1d0*phi) + cos(omega * y(i,j,k)  -1.2d0*phi) - sin(omega * z(i,j,k)  + 1.3d0*phi)
            fFx_exact(i,j,k) = 1.5d0 + sin(omega *  xF(i,j,k) + 1.1d0*phi) + cos(omega * y(i,j,k)  -1.2d0*phi) - sin(omega * z(i,j,k)  + 1.3d0*phi)
            fFy_exact(i,j,k) = 1.5d0 + sin(omega *  x(i,j,k)  + 1.1d0*phi) + cos(omega * yF(i,j,k) -1.2d0*phi) - sin(omega * z(i,j,k)  + 1.3d0*phi)
            fFz_exact(i,j,k) = 1.5d0 + sin(omega *  x(i,j,k)  + 1.1d0*phi) + cos(omega * y(i,j,k)  -1.2d0*phi) - sin(omega * zF(i,j,k) + 1.3d0*phi)
        end do
        end do 
        end do 

        !TODO: delete below comments once done, these are nice for checking that
        !grid is concstructed correctly
        !print *, "x/y/z(1) Nodes"
        !print *, x(1,1,1)
        !print *, y(1,1,1)
        !print *, z(1,1,1)
        !print *, "x/y/z(n) Nodes"
        !print *, x(nx,ny,nz)
        !print *, y(nx,ny,nz)
        !print *, z(nx,ny,nz)

        !print *, "x/y/z(1) Faces"
        !print *, xF(1,1,1)
        !print *, yF(1,1,1)
        !print *, zF(1,1,1)
        !print *, "x/y/z(n) Faces"
        !print *, xF(nx,ny,nz)
        !print *, yF(nx,ny,nz)
        !print *, zF(nx,ny,nz)

        deallocate( x )
        deallocate( y )
        deallocate( z )
        deallocate( xF )
        deallocate( yF )
        deallocate( zF )
        
        allocate( fF(nx,ny,nz) )
        allocate( fN(nx,ny,nz) )

        print*, "Created initial data"


        !Note: periodicity is controlled by true/false statements
        call method3%init(          nx,     ny,    nz, &
                                    dx,     dy,    dz, &
                     periodic_x,periodic_y,periodic_z, &
                                   "ei02", "ei02", "ei02" )
         
        call method5%init(          nx,     ny,    nz, &
                                    dx,     dy,    dz, &
                     periodic_x,periodic_y,periodic_z, &
                                   "ei06", "ei06", "ei06" )
         
        print*, "Initialized all methods"
  
        print*, "==========================================="
        print*, "Now trying METHOD 3: EI02 (N2F)"
        print*, "==========================================="
        call method3 % iN2Fx(f,fF,bc1_x,bcn_x)
        if (periodic_x) then
           error = MAXVAL( ABS(fF - fFx_exact))
        else
           error = MAXVAL( ABS(fF(1:nx-1,:,:) - fFx_exact(1:nx-1,:,:)))
        endif
        print*, "Maximum error = ", error
        ei02_N2F_error(ind,1) = error

        call method3 % iN2Fy(f,fF,bc1_y,bcn_y)
        if (periodic_y) then
           error = MAXVAL( ABS(fF - fFy_exact))
        else
           error = MAXVAL( ABS(fF(:,1:ny-1,:) - fFy_exact(:,1:ny-1,:)))
        endif
        print*, "Maximum error = ", error
        ei02_N2F_error(ind,2) = error

        call method3 % iN2Fz(f,fF,bc1_z,bcn_z)
        if (periodic_z) then
           error = MAXVAL( ABS(fF - fFz_exact))
        else
           error = MAXVAL( ABS(fF(:,:,1:nz-1) - fFz_exact(:,:,1:nz-1)))
        endif
        print*, "Maximum error = ", error
        ei02_N2F_error(ind,3) = error
        
        print*, "==========================================="
        print*, "Now trying METHOD 3: EI02 (F2N)"
        print*, "==========================================="
        call method3 % iF2Nx(fFx_exact,fN,bc1_x,bcn_x)
        error = MAXVAL( ABS(f - fN))
        print*, "Maximum error = ", error
        ei02_F2N_error(ind,1) = error

        call method3 % iF2Ny(fFy_exact,fN,bc1_y,bcn_y)
        error = MAXVAL( ABS(f - fN))
        print*, "Maximum error = ", error
        ei02_F2N_error(ind,2) = error

        call method3 % iF2Nz(fFz_exact,fN,bc1_z,bcn_z)
        error = MAXVAL( ABS(f - fN))
        print*, "Maximum error = ", error
        ei02_F2N_error(ind,3) = error

        print*, "==========================================="
        print*, "Now trying METHOD 5: EI06 (N2F)"
        print*, "==========================================="
        call method5 % iN2Fx(f,fF,bc1_x,bcn_x)
        if (periodic_x) then
           error = MAXVAL( ABS(fF - fFx_exact))
        else
           error = MAXVAL( ABS(fF(1:nx-1,:,:) - fFx_exact(1:nx-1,:,:)))
        endif
        print*, "Maximum error = ", error
        ei06_N2F_error(ind,1) = error

        call method5 % iN2Fy(f,fF,bc1_y,bcn_y)
        if (periodic_y) then
           error = MAXVAL( ABS(fF - fFy_exact))
        else
           error = MAXVAL( ABS(fF(:,1:ny-1,:) - fFy_exact(:,1:ny-1,:)))
        endif
        print*, "Maximum error = ", error
        ei06_N2F_error(ind,2) = error

        call method5 % iN2Fz(f,fF,bc1_z,bcn_z)
        if (periodic_z) then
           error = MAXVAL( ABS(fF - fFz_exact))
        else
           error = MAXVAL( ABS(fF(:,:,1:nz-1) - fFz_exact(:,:,1:nz-1)))
        endif
        print*, "Maximum error = ", error
        ei06_N2F_error(ind,3) = error
        
        print*, "==========================================="
        print*, "Now trying METHOD 5: EI06 (F2N)"
        print*, "==========================================="
        call method5 % iF2Nx(fFx_exact,fN,bc1_x,bcn_x)
        error = MAXVAL( ABS(f - fN))
        print*, "Maximum error = ", error
        ei06_F2N_error(ind,1) = error

        call method5 % iF2Ny(fFy_exact,fN,bc1_y,bcn_y)
        error = MAXVAL( ABS(f - fN))
        print*, "Maximum error = ", error
        ei06_F2N_error(ind,2) = error

        call method5 % iF2Nz(fFz_exact,fN,bc1_z,bcn_z)
        error = MAXVAL( ABS(f - fN))
        print*, "Maximum error = ", error
        ei06_F2N_error(ind,3) = error



        print*, "---------------------------------------------"
        

        
        print*, "=========================================="
        print*, " Destroying everything"
        print*, "=========================================="

        deallocate( f )
        deallocate( fFx_exact )
        deallocate( fFy_exact )
        deallocate( fFz_exact )
        deallocate( fF )
        deallocate( fN )
        call method3%destroy
        call method5%destroy

    enddo

    print*, "==============================================================="
    print*, "====Midpoint interpolation from Nodes to Faces (N2F) Results==="
    print*, "==============================================================="
    print*, "==============================================================="
    print*, "=======================2nd Order Explicit======================"
    print*, "==============================================================="
    ratio = ei02_N2F_error(1:3,1)/ei02_N2F_error(2:4,1) 
    print *, "X Interp Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = ei02_N2F_error(1:3,2)/ei02_N2F_error(2:4,2) 
    print *, "Y Interp Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = ei02_N2F_error(1:3,3)/ei02_N2F_error(2:4,3) 
    print *, "Z Interp Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "==============================================================="
    print*, "=======================6nd Order Explicit======================"
    print*, "==============================================================="
    ratio = ei06_N2F_error(1:3,1)/ei06_N2F_error(2:4,1) 
    print *, "X Interp Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = ei06_N2F_error(1:3,2)/ei06_N2F_error(2:4,2) 
    print *, "Y Interp Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = ei06_N2F_error(1:3,3)/ei06_N2F_error(2:4,3) 
    print *, "Z Interp Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "==============================================================="

    print*, "==============================================================="
    print*, "====Midpoint interpolation from Nodes to Faces (F2N) Results==="
    print*, "==============================================================="
    print*, "==============================================================="
    print*, "=======================2nd Order Explicit======================"
    print*, "==============================================================="
    ratio = ei02_F2N_error(1:3,1)/ei02_F2N_error(2:4,1) 
    print *, "X Interp Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = ei02_F2N_error(1:3,2)/ei02_F2N_error(2:4,2) 
    print *, "Y Interp Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = ei02_F2N_error(1:3,3)/ei02_F2N_error(2:4,3) 
    print *, "Z Interp Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "==============================================================="
    print*, "=======================6nd Order Explicit======================"
    print*, "==============================================================="
    ratio = ei06_F2N_error(1:3,1)/ei06_F2N_error(2:4,1) 
    print *, "X Interp Order of Convergence (x)"
    print *,  log(ratio) / log(2.0)
    ratio = ei06_F2N_error(1:3,2)/ei06_F2N_error(2:4,2) 
    print *, "Y Interp Order of Convergence (y)"
    print *,  log(ratio) / log(2.0)
    ratio = ei06_F2N_error(1:3,3)/ei06_F2N_error(2:4,3) 
    print *, "Z Interp Order of Convergence (z)"
    print *,  log(ratio) / log(2.0)
    print*, "==============================================================="

    ! End the run
    call MPI_Finalize(ierr)

end program
