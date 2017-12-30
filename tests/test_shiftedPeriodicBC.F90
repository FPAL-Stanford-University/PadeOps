program test_shiftedPeriodicBC
    use kind_parameters, only: rkind, clen
    use constants, only: pi, two, one, imi, zero, half, kappa
    use reductions, only: p_maxval, p_sum, p_minval
    use timer, only: tic, toc
    use decomp_2d
    use decomp_2d_io
    use FringeMethod, only: fringe
    use exits, only: message
    use spectralMod, only: spectral  
    use mpi

    implicit none 

    type(fringe), allocatable :: fringe_x, fringe_x_spbc
    type(spectral), allocatable, target :: spectE, spectC
    integer, parameter :: nx = 192, ny = 96, nz = 128
    character(len=clen) :: inputfile = "/home/nghaisas/runs/PadeOps/tests/spbc/input.in"
    character(len=clen) :: inputfile_spbc = "/home/nghaisas/runs/PadeOps/tests/spbc/input_spbc.in"
    character(len=clen) :: filter_x, filter_y, filter_z          ! What filter to use: "cf90", "gaussian", "lstsq", "spectral"
    real(rkind), dimension(:,:,:), allocatable :: xG, yG, zG
    real(rkind), dimension(:,:,:), allocatable :: u, v, w, rhsu, rhsv, rhsw
    real(rkind), dimension(:,:,:,:), allocatable :: rbuffxC, rbuffyC, rbuffzC
    real(rkind), dimension(:,:,:,:), allocatable :: rbuffxE, rbuffyE, rbuffzE
    complex(rkind), dimension(:,:,:,:), allocatable :: cbuffyC, cbuffyE
    complex(rkind), dimension(:,:,:), pointer :: u_rhs, v_rhs, w_rhs 
    complex(rkind), dimension(:,:,:,:), allocatable, target :: rhsC, rhsE
    real(rkind), parameter :: Lx = pi, Ly = 0.5d0*pi, Lz = 1.0d0
    real(rkind) :: dx, dy, dz, dt = 1.0d0
    type(decomp_info) :: gpC, gpE
    type(decomp_info), pointer :: sp_gpC, sp_gpE 
    integer :: idx, ix1, iy1, iz1, ixn, iyn, izn, i, j, k, ierr, prow = 0, pcol = 0, num_turbines 
    real(rkind) :: xPeriods = 2.d0, yPeriods = 2.d0, zpeak = 0.3d0, epsnd = 5.d0, z0init = 1.d-4 
    real(rkind) :: inst_val(8)
    real(rkind) :: comp1, comp2, maxdiff

    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gpC)
    call decomp_info_init(nx,ny,nz+1,gpE)

    allocate(xG(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3))); allocate(yG(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(zG(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3))); allocate(u (gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(v (gpC%xsz(1),gpC%xsz(2),gpC%xsz(3))); allocate(w (gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
    allocate(rhsu(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3))) 
    allocate(rhsv(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3))) 
    allocate(rhsw(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3))) 

    dx = Lx/real(nx,rkind); dy = Ly/real(ny,rkind); dz = Lz/real(nz,rkind)
    ix1 = gpC%xst(1); iy1 = gpC%xst(2); iz1 = gpC%xst(3)
    ixn = gpC%xen(1); iyn = gpC%xen(2); izn = gpC%xen(3)
    do k=1,size(xG,3)
        do j=1,size(xG,2)
            do i=1,size(xG,1)
                xG(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                yG(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                zG(i,j,k) = real( iz1 + k - 1, rkind ) * dz + dz/two

                u(i,j,k) = real(iy1 + j, rkind);
                v(i,j,k) = real(iy1 + j, rkind)*2.0d0
                w(i,j,k) = real(iy1 + j, rkind)*3.0d0
            end do
        end do
    end do
    xG = xG - dx; yG = yG - dy; zG = zG - dz 

    filter_x = "spectral";     filter_y = "spectral";     filter_z = "spectral"
    allocate(spectC)
    call spectC%init("x",nx,ny,nz,dx,dy,dz,"four",filter_x,2,fixOddball=.false.,exhaustiveFFT=.false., init_periodicInZ=.false.,dealiasF=2.0d0/3.0d0)
    allocate(spectE)
    call spectE%init("x",nx,ny,nz+1,dx,dy,dz,"four",filter_x,2,fixOddball=.false.,exhaustiveFFT=.false.,init_periodicInZ=.false.,dealiasF=2.0d0/3.0d0)

    sp_gpC => spectC%spectdecomp
    sp_gpE => spectE%spectdecomp

    allocate(rbuffxC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3),3))
    allocate(rbuffyC(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3),2))
    allocate(rbuffxE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3),3))
    allocate(rbuffyE(gpE%ysz(1),gpE%ysz(2),gpE%ysz(3),2))
    allocate(cbuffyC(sp_gpC%ysz(1),sp_gpC%ysz(2),sp_gpC%ysz(3),2))
    allocate(cbuffyE(sp_gpE%ysz(1),sp_gpE%ysz(2),sp_gpE%ysz(3),2))

    call spectC%alloc_r2c_out(rhsC,3) 
    call spectE%alloc_r2c_out(rhsE,1)
    u_rhs => rhsC(:,:,:,1); v_rhs => rhsC(:,:,:,2); w_rhs => rhsE(:,:,:,1)

    allocate(fringe_x)
    call fringe_x%init(inputfile, dx, xG(:,1,1), dy, yG(1,:,1), spectC, spectE, gpC, gpE, rbuffxC, rbuffxE, rbuffyC, rbuffyE, cbuffyC, cbuffyE)

    allocate(fringe_x_spbc)
    call fringe_x_spbc%init(inputfile_spbc, dx, xG(:,1,1), dy, yG(1,:,1), spectC, spectE, gpC, gpE, rbuffxC, rbuffxE, rbuffyC, rbuffyE, cbuffyC, cbuffyE)

    call fringe_x%associateFringeTargets(u, v, w) !<-- Link the target velocity array to igp 
    rhsu = 0.d0; rhsv = 0.0d0; rhsw = 0.0d0
    u_rhs = 0.0d0; v_rhs = 0.0d0; w_rhs = 0.0d0
    call mpi_barrier(mpi_comm_world, ierr)
    call tic()
    call fringe_x%addFringeRHS(dt, u_rhs, v_rhs, w_rhs, u, v, w)
    call spectC%ifft(u_rhs, rhsu)
    call spectC%ifft(v_rhs, rhsv)
    call spectE%ifft(w_rhs, rhsw)
    call mpi_barrier(mpi_comm_world, ierr)
    call toc()
    !call decomp_2d_write_one(1,rhs1,"temp_T1.bin", gp)
    call message(0,"No shifted periodic BC")
    call message(1,"Max rhs u:", p_maxval(rhsu))
    call message(1,"Min rhs u:", p_minval(rhsu))
    call message(1,"Max rhs v:", p_maxval(rhsv))
    call message(1,"Min rhs v:", p_minval(rhsv))
    call message(1,"Max rhs w:", p_maxval(rhsw))
    call message(1,"Min rhs w:", p_minval(rhsw))
  
    rhsu = 0.d0; rhsv = 0.0d0; rhsw = 0.0d0
    u_rhs = 0.0d0; v_rhs = 0.0d0; w_rhs = 0.0d0
    call mpi_barrier(mpi_comm_world, ierr)
    call tic()
    call fringe_x_spbc%addFringeRHS(dt, u_rhs, v_rhs, w_rhs, u, v, w)
    call spectC%ifft(u_rhs, rhsu)
    call spectC%ifft(v_rhs, rhsv)
    call spectE%ifft(w_rhs, rhsw)
    call mpi_barrier(mpi_comm_world, ierr)
    call toc()
    !call decomp_2d_write_one(1,rhs1,"temp_T1.bin", gp)
    call message(0,"Shifted periodic BC")
    call message(1,"Max rhs u:", p_maxval(rhsu))
    call message(1,"Min rhs u:", p_minval(rhsu))
    call message(1,"Max rhs v:", p_maxval(rhsv))
    call message(1,"Min rhs v:", p_minval(rhsv))
    call message(1,"Max rhs w:", p_maxval(rhsw))
    call message(1,"Min rhs w:", p_minval(rhsw))


    deallocate(fringe_x_spbc)
    deallocate(fringe_x)
    deallocate(rhsC, rhsE)
    deallocate(cbuffyC, cbuffyE, rbuffyE, rbuffyC, rbuffxC, rbuffxE)
    deallocate(spectE, spectC)
    deallocate(xG, yG, zG, u, v, w, rhsu, rhsv, rhsw)
    call MPI_Finalize(ierr)
end program 
