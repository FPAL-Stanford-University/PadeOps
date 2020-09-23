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
    real(rkind), dimension(:,:,:), allocatable :: x,y,z,fx,fy,fz,fF, fN
    real(rkind), dimension(:,:,:), allocatable :: xF,yF,zF, fFx_exact, fFy_exact, fFz_exact
    real(rkind) :: dx, dy, dz, Lx, Ly, Lz
    real(rkind) :: omega_x, omega_y, omega_z, phi_x, phi_y, phi_z

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

    !set phi and omega to achieve desired BCs    
    if (periodic_x) then
        phi_x = 2.345d0
        omega_x = 4.0d0
    else
        if (bc1_x .eq. 1) then           !1: symmetric
            phi_x = half * pi              
            if (bcn_x .eq. 1) then       !n: symmetric
                omega_x = 2.0d0        
            else if (bcn_x .eq. -1) then !n: anti-symmetric
                omega_x = 2.25d0           
            else                         !n: general
                omega_x = 2.34d0          
            endif
            
        else if (bc1_x .eq. -1) then
            phi_x = 0.0d0                !1: anti-symmetric
            if (bcn_x .eq. 1) then       !n: symmetric
                omega_x = 2.25d0        
            else if (bcn_x .eq. -1) then !n: anti-symmetric
                omega_x = 2.0d0           
            else                         !n: general
                omega_x = 2.113d0          
            endif
        else
            phi_x = 2.374d0          
            if (bcn_x .eq. 1) then       !n: symmetric
                omega_x = (5.0d0*half*pi - phi_x) / (2.0d0 * pi)
            else if (bcn_x .eq. -1) then !n: anti-symmetric
                omega_x = (2.0d0*pi - phi_x) / (2.0d0 * pi)
            else                         !n: general
                omega_x = 2.213d0          
            endif
        endif
    endif
    if (periodic_y) then
        phi_y = 2.345d0
        omega_y = 4.0d0
    else
        if (bc1_y .eq. 1) then           !1: symmetric
            phi_y = half * pi              
            if (bcn_y .eq. 1) then       !n: symmetric
                omega_y = 2.0d0        
            else if (bcn_y .eq. -1) then !n: anti-symmetric
                omega_y = 2.25d0           
            else                         !n: general
                omega_y = 2.34d0          
            endif
            
        else if (bc1_y .eq. -1) then
            phi_y = 0.0d0                !1: anti-symmetric
            if (bcn_y .eq. 1) then       !n: symmetric
                omega_y = 2.25d0        
            else if (bcn_y .eq. -1) then !n: anti-symmetric
                omega_y = 2.0d0           
            else                         !n: general
                omega_y = 2.113d0          
            endif
        else
            phi_y = 2.374d0          
            if (bcn_y .eq. 1) then       !n: symmetric
                omega_y = (5.0d0*half*pi - phi_y) / (2.0d0 * pi)
            else if (bcn_y .eq. -1) then !n: anti-symmetric
                omega_y = (2.0d0*pi - phi_y) / (2.0d0 * pi)
            else                         !n: general
                omega_y = 2.213d0          
            endif
        endif
    endif
    if (periodic_z) then
        phi_z = 2.345d0
        omega_z = 4.0d0
    else
        if (bc1_z .eq. 1) then           !1: symmetric
            phi_z = half * pi              
            if (bcn_z .eq. 1) then       !n: symmetric
                omega_z = 2.0d0        
            else if (bcn_z .eq. -1) then !n: anti-symmetric
                omega_z = 2.25d0           
            else                         !n: general
                omega_z = 2.34d0          
            endif
            
        else if (bc1_z .eq. -1) then
            phi_z = 0.0d0                !1: anti-symmetric
            if (bcn_z .eq. 1) then       !n: symmetric
                omega_z = 2.25d0        
            else if (bcn_z .eq. -1) then !n: anti-symmetric
                omega_z = 2.0d0           
            else                         !n: general
                omega_z = 2.113d0          
            endif
        else
            phi_z = 2.374d0          
            if (bcn_z .eq. 1) then       !n: symmetric
                omega_z = (5.0d0*half*pi - phi_z) / (2.0d0 * pi)
            else if (bcn_z .eq. -1) then !n: anti-symmetric
                omega_z = (2.0d0*pi - phi_z) / (2.0d0 * pi)
            else                         !n: general
                omega_z = 2.213d0          
            endif
        endif
    endif

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
        allocate( fx(nx,ny,nz) )
        allocate( fy(nx,ny,nz) )
        allocate( fz(nx,ny,nz) )
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
            fx(i,j,k)        = sin(omega_x * x(i,j,k)  + phi_x)
            fy(i,j,k)        = sin(omega_y * y(i,j,k)  + phi_y)
            fz(i,j,k)        = sin(omega_z * z(i,j,k)  + phi_z)
            fFx_exact(i,j,k) = sin(omega_x * xF(i,j,k) + phi_x)
            fFy_exact(i,j,k) = sin(omega_y * yF(i,j,k) + phi_y) 
            fFz_exact(i,j,k) = sin(omega_z * zF(i,j,k) + phi_z)
        end do
        end do 
        end do 

        if (ind .eq. 1) then
            print *, y(1,:,1)
            print *, fy(1,:,1)
        endif

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
                                   "ei02", "ei02", "ei02" )
                print *, "Warning: overriding ei06 init until BCs implemented"
!                                   "ei06", "ei06", "ei06" )
         
        print*, "Initialized all methods"
  
        print*, "==========================================="
        print*, "Now trying METHOD 3: EI02 (N2F)"
        print*, "==========================================="
        call method3 % iN2Fx(fx,fF,bc1_x,bcn_x)
        if (periodic_x) then
           error = MAXVAL( ABS(fF - fFx_exact))
        else
           error = MAXVAL( ABS(fF(1:nx-1,:,:) - fFx_exact(1:nx-1,:,:)))
        endif
        print*, "Maximum error = ", error
        ei02_N2F_error(ind,1) = error

        call method3 % iN2Fy(fy,fF,bc1_y,bcn_y)
        if (periodic_y) then
           error = MAXVAL( ABS(fF - fFy_exact))
        else
           error = MAXVAL( ABS(fF(:,1:ny-1,:) - fFy_exact(:,1:ny-1,:)))
        endif
        print*, "Maximum error = ", error
        ei02_N2F_error(ind,2) = error

        call method3 % iN2Fz(fz,fF,bc1_z,bcn_z)
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
        error = MAXVAL( ABS(fx - fN))
        print*, "Maximum error = ", error
        ei02_F2N_error(ind,1) = error

        call method3 % iF2Ny(fFy_exact,fN,bc1_y,bcn_y)
        error = MAXVAL( ABS(fy - fN))
        print*, "Maximum error = ", error
        ei02_F2N_error(ind,2) = error

        call method3 % iF2Nz(fFz_exact,fN,bc1_z,bcn_z)
        error = MAXVAL( ABS(fz - fN))
        print*, "Maximum error = ", error
        ei02_F2N_error(ind,3) = error

        print*, "==========================================="
        print*, "Now trying METHOD 5: EI06 (N2F)"
        print*, "==========================================="
        call method5 % iN2Fx(fx,fF,bc1_x,bcn_x)
        if (periodic_x) then
           error = MAXVAL( ABS(fF - fFx_exact))
        else
           error = MAXVAL( ABS(fF(1:nx-1,:,:) - fFx_exact(1:nx-1,:,:)))
        endif
        print*, "Maximum error = ", error
        ei06_N2F_error(ind,1) = error

        call method5 % iN2Fy(fy,fF,bc1_y,bcn_y)
        if (periodic_y) then
           error = MAXVAL( ABS(fF - fFy_exact))
        else
           error = MAXVAL( ABS(fF(:,1:ny-1,:) - fFy_exact(:,1:ny-1,:)))
        endif
        print*, "Maximum error = ", error
        ei06_N2F_error(ind,2) = error

        call method5 % iN2Fz(fz,fF,bc1_z,bcn_z)
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
        error = MAXVAL( ABS(fx - fN))
        print*, "Maximum error = ", error
        ei06_F2N_error(ind,1) = error

        call method5 % iF2Ny(fFy_exact,fN,bc1_y,bcn_y)
        error = MAXVAL( ABS(fy - fN))
        print*, "Maximum error = ", error
        ei06_F2N_error(ind,2) = error

        call method5 % iF2Nz(fFz_exact,fN,bc1_z,bcn_z)
        error = MAXVAL( ABS(fz - fN))
        print*, "Maximum error = ", error
        ei06_F2N_error(ind,3) = error



        print*, "---------------------------------------------"
        

        
        print*, "=========================================="
        print*, " Destroying everything"
        print*, "=========================================="

        deallocate( fx )
        deallocate( fy )
        deallocate( fz )
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
