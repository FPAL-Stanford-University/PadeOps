program test_actuatorLine
    use kind_parameters, only: rkind, clen
    use constants, only: pi, two, one, imi, zero, half, kappa
    use reductions, only: p_maxval, p_sum 
    use timer, only: tic, toc
    use decomp_2d
    use decomp_2d_io
    use actuatorLineMod, only: actuatorLine
    use exits, only: message
    use mpi

    implicit none 

    type(actuatorLine), dimension(:), allocatable :: hawts
    integer, parameter :: nx = 96, ny = 96, nz = 64
    character(len=clen) :: inputDir = "/home/nghaisas/PadeOps/tests/test_actuatorLine_files"
    real(rkind), dimension(:,:,:), allocatable :: xG, yG, zG, yRightHalo, zRightHalo, zLeftHalo
    real(rkind), dimension(:,:,:), allocatable :: u, v, w, rhs, rhsv, rhsw
    real(rkind), parameter :: Lx = pi, Ly = pi, Lz = one
    real(rkind) :: dx, dy, dz
    type(decomp_info) :: gp 
    integer :: idx, ix1, iy1, iz1, ixn, iyn, izn, i, j, k, ierr, prow = 2, pcol = 1 
    real(rkind) :: xPeriods = 2.d0, yPeriods = 2.d0, zpeak = 0.3d0, epsnd = 5.d0, z0init = 1.d-4 
    real(rkind) :: inst_val(8), dt = 0.1_rkind, xyzPads(6)

    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)

    allocate(xG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(yG(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(zG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(u(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(v(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(w(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(rhs(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(rhsv(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(rhsw(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(yRightHalo(gp%xsz(1),gp%xsz(3),  3))
    allocate(zRightHalo(gp%xsz(1),gp%xsz(2)+1,3))
    allocate(zLeftHalo (gp%xsz(1),gp%xsz(2)+1,3))


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
        if(nrank==0) write(*,*) '-----------Initializing Turbine No.', idx, '------------'
        call hawts(idx)%init(inputDir, idx, xG, yG, zG, xyzPads)
        if(nrank==0) then
            write(*,*) hawts(idx)%Am_I_Active, hawts(idx)%Am_I_Split
        endif
        if(nrank==0) write(*,*) '-----------Initialized Turbine No.', idx, '------------'
    end do 
    rhs = 0.d0
    call mpi_barrier(mpi_comm_world, ierr)
    call tic()
    !call halo_communicate
    do idx = 1,6
        call hawts(idx)%get_RHS(dt, u, v, w, yRightHalo, zLeftHalo, zRightHalo, rhs, rhsv, rhsw, inst_val)
    end do 
    call mpi_barrier(mpi_comm_world, ierr)
    call toc()
 
    call decomp_2d_write_one(1,rhs,"temp.bin", gp)
    
    call message(2,"Computed Source:", p_sum(sum(rhs)) * dx*dy*dz)
    call message(2,"Expected Source:", (6.d0*0.5d0*(pi/4.d0)*(0.08d0**2)*0.65d0))

    call message(2,"error:", 100.d0*abs(p_sum(sum(rhs)) * dx*dy*dz - 6.d0*0.5d0*(pi/4.d0)*(0.08d0**2)*0.65d0) / (6.d0*0.5d0*(pi/4.d0)*(0.08d0**2)*0.65d0))
    do idx = 1,6
    call hawts(idx)%destroy()
    end do 
    deallocate(hawts)
    deallocate(xG, yG, zG, u, v, w, rhs,yRightHalo, zRightHalo, zLeftHalo)
    call MPI_Finalize(ierr)
end program 
