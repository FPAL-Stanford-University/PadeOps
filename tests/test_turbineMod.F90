program test_actuatorDisk
    use kind_parameters, only: rkind, clen
    use constants, only: pi, two, one, imi, zero, half, kappa
    use reductions, only: p_maxval, p_sum 
    use timer, only: tic, toc
    use decomp_2d
    use decomp_2d_io
    use turbineMod, only: turbineArray
    use spectralMod, only: spectral
    use exits, only: message
    use mpi

    implicit none 

    type(turbineArray), allocatable :: turbArray
    integer, parameter :: nx = 192, ny = 192, nz = 128
    character(len=clen) :: inputDir = "/home/nghaisas/ActuatorDisk/"
    character(len=clen) :: inputFile = "/home/nghaisas/PadeOps/problems/incompressible/pbl_files/input_pbl.dat"
    real(rkind), dimension(:,:,:,:), allocatable :: mesh
    real(rkind), dimension(:,:,:), allocatable :: u, v, w, rhsx, rhsy, rhsz
    complex(rkind), dimension(:,:,:), allocatable :: urhs, vrhs, wrhs
    real(rkind), parameter :: Lx = pi, Ly = pi, Lz = one
    real(rkind) :: dx, dy, dz, dt
    type(decomp_info) :: gpC, gpE
    type(decomp_info), pointer :: sp_gpC, sp_gpE
    type(spectral), allocatable, target :: spectC, spectE 
    integer :: idx, ix1, iy1, iz1, ixn, iyn, izn, i, j, k, ierr, prow = 0, pcol = 0 
    real(rkind) :: xPeriods = 2.d0, yPeriods = 2.d0, zpeak = 0.3d0, epsnd = 5.d0, z0init = 1.d-4 
    real(rkind), dimension(:,:,:,:), allocatable :: rbuffxC
    complex(rkind), dimension(:,:,:,:), allocatable :: cbuffyC, cbuffzC, cbuffyE, cbuffzE
    real(rkind), dimension(:,:,:), allocatable :: rhs_real
    complex(rkind) :: zeroc = 0.d0 + imi*0.0d0

    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(gpC)
    call decomp_info_init(nx,ny,nz+1,gpE)

    allocate(mesh(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3),3))
    allocate(u(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(rhs_real(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(v(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3))); allocate(w(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(rhsx(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3))) 
    allocate(rhsy(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3))) 
    allocate(rhsz(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3))) 
    allocate(urhs(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3))) 
    allocate(vrhs(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3))) 
    allocate(wrhs(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3))) 

    dx = Lx/real(nx,rkind); dy = Ly/real(ny,rkind); dz = Lz/real(nz,rkind)
    ix1 = gpC%xst(1); iy1 = gpC%xst(2); iz1 = gpC%xst(3)
    ixn = gpC%xen(1); iyn = gpC%xen(2); izn = gpC%xen(3)
    do k=1,size(mesh,3)
        do j=1,size(mesh,2)
            do i=1,size(mesh,1)
                mesh(i,j,k,1) = real( ix1 + i - 1, rkind ) * dx
                mesh(i,j,k,2) = real( iy1 + j - 1, rkind ) * dy
                mesh(i,j,k,3) = real( iz1 + k - 1, rkind ) * dz + dz/two
            end do
        end do
    end do
    mesh(:,:,:,1) = mesh(:,:,:,1) - dx; mesh(:,:,:,2) = mesh(:,:,:,2) - dy; mesh(:,:,:,3) = mesh(:,:,:,3) - dz 
    u = one; v = zero; w = zero

    allocate(spectC)
    call spectC%init("x", nx, ny, nz, dx, dy, dz, &
            "four", '2/3rd', 2 , .false.)
    allocate(spectE)
    call spectE%init("x", nx, ny, nz+1, dx, dy, dz, &
            "four", '2/3rd', 2 , .false.)
    sp_gpC => spectC%spectdecomp
    sp_gpE => spectE%spectdecomp
    
    allocate(rbuffxC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3),2))
    allocate(cbuffyC(sp_gpC%ysz(1),sp_gpC%ysz(2),sp_gpC%ysz(3),1))
    allocate(cbuffyE(sp_gpE%ysz(1),sp_gpE%ysz(2),sp_gpE%ysz(3),1))
    allocate(cbuffzC(sp_gpC%zsz(1),sp_gpC%zsz(2),sp_gpC%zsz(3),1))
    allocate(cbuffzE(sp_gpE%zsz(1),sp_gpE%zsz(2),sp_gpE%zsz(3),1))

    allocate(turbArray)
    call turbArray%init(inputFile, gpC, gpE, sp_gpC, sp_gpE, spectC, spectE, rbuffxC, cbuffyC, cbuffyE, cbuffzC, cbuffzE, mesh, dx, dy, dz) 

    dt = 0.1_rkind
    urhs = zeroc; vrhs = zeroc; wrhs = zeroc
    call tic()
    call turbArray%getForceRHS(dt, u, v, w, urhs, vrhs, wrhs)
    call toc()
    call spectC%ifft(urhs,rhs_real) 
    call decomp_2d_write_one(1,rhs_real,"ArrayRHS.bin", gpC)
         
    call turbArray%destroy()
    deallocate(turbArray)
    deallocate(mesh, u, v, w, rhsx, rhsy, rhsz)
    call MPI_Finalize(ierr)
end program 
