program test_actuatorDisk
    use kind_parameters, only: rkind, clen
    use constants, only: pi, two, one, imi, zero, half, kappa
    use reductions, only: p_maxval, p_sum 
    use timer, only: tic, toc
    use decomp_2d
    use decomp_2d_io
    use actuatorDiskMod, only: actuatorDisk
    use exits, only: message
    use mpi

    implicit none 

    type(actuatorDisk), dimension(:), allocatable :: hawts
    integer, parameter :: nx = 192, ny = 192, nz = 128
    character(len=clen) :: inputDir = "/home/aditya90/Codes/PadeOps/data/ActuatorDisk/"
    real(rkind), dimension(:,:,:), allocatable :: xG, yG, zG
    real(rkind), dimension(:,:,:), allocatable :: u, v, w, rhs
    real(rkind), parameter :: Lx = pi, Ly = pi, Lz = one
    real(rkind) :: dx, dy, dz
    type(decomp_info) :: gp 
    integer :: idx, ix1, iy1, iz1, ixn, iyn, izn, i, j, k, ierr, prow = 0, pcol = 0 
    real(rkind) :: xPeriods = 2.d0, yPeriods = 2.d0, zpeak = 0.3d0, epsnd = 5.d0, z0init = 1.d-4 

    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(gp)

    allocate(xG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(yG(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(zG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(u(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(v(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(w(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(rhs(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 

    dx = Lx/real(nx,rkind); dy = Ly/real(ny,rkind); dz = Lz/real(nz,rkind)
    ix1 = gp%xst(1); iy1 = gp%xst(2); iz1 = gp%xst(3)
    ixn = gp%xen(1); iyn = gp%xen(2); izn = gp%xen(3)
    do k=1,size(xG,3)
        do j=1,size(xG,2)
            do i=1,size(xG,1)
                xG(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                yG(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                zG(i,j,k) = real( iz1 + k - 1, rkind ) * dz + dz/two
            end do
        end do
    end do
    xG = xG - dx; yG = yG - dy; zG = zG - dz 

    u = (one/kappa)*log(zG/z0init) + epsnd*cos(Yperiods*two*pi*yG/Ly)*exp(-half*(zG/zpeak/Lz)**2)
    v = epsnd*(zG/Lz)*cos(Xperiods*two*pi*xG/Lx)*exp(-half*(zG/zpeak/Lz)**2)
    w= zero  
    u = one; v = zero; w = zero
    allocate(hawts(6))
    do idx = 1,6
        call hawts(idx)%init(inputDir, idx, xG, yG, zG)
    end do 
    rhs = 0.d0
    call mpi_barrier(mpi_comm_world, ierr)
    call tic()
    do idx = 1,6
        call hawts(idx)%get_RHS(u, v, w, rhs)
    end do 
    call toc()
    call mpi_barrier(mpi_comm_world, ierr)
   
    call decomp_2d_write_one(1,rhs,"temp.bin", gp)
    
    print*, p_sum(sum(rhs))
    do idx = 1,6
    call hawts(idx)%destroy()
    end do 
    deallocate(hawts)
    deallocate(xG, yG, zG, u, v, w, rhs)
    call MPI_Finalize(ierr)
end program 
