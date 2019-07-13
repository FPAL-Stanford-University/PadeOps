program test_actuatorDiscRot
    use kind_parameters, only: rkind, clen
    use constants, only: pi, two, one, imi, zero, half, kappa
    use reductions, only: p_maxval, p_sum, p_minval 
    use timer, only: tic, toc
    use decomp_2d
    use decomp_2d_io
    use actuatorDisk_RotMod, only: actuatorDisk_Rot
    use exits, only: message
    use mpi

    implicit none 

    type(actuatorDisk_Rot), dimension(:), allocatable :: hawts_Rot
    integer, parameter :: nx = 192, ny = 96, nz = 128
    character(len=clen) :: inputDir = "/fastscratch/nghaisas/runs/PadeOps/wupa/run3/deb/turbInfo"
    real(rkind), dimension(:,:,:), allocatable :: xG, yG, zG
    real(rkind), dimension(:,:,:), allocatable :: u, v, w, rhs1, rhsv, rhsw, rhs2
    real(rkind), parameter :: Lx = 28.8d0, Ly = 4.8d0, Lz = 3.0d0
    real(rkind) :: dx, dy, dz, diam, CT
    type(decomp_info) :: gp 
    integer :: idx, ix1, iy1, iz1, ixn, iyn, izn, i, j, k, ierr, prow = 0, pcol = 0, num_turbines 
    real(rkind) :: xPeriods = 2.d0, yPeriods = 2.d0, zpeak = 0.3d0, epsnd = 5.d0, z0init = 1.d-4 
    real(rkind) :: inst_val(8)
    real(rkind) :: comp1, comp2, maxdiff

    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)

    allocate(xG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(yG(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(zG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(u(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(v(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(w(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(rhs1(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(rhs2(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(rhsv(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(rhsw(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 

    num_turbines = 1

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

    u = one 
    v = zero 
    w = zero  
    u = one; v = zero; w = zero
    allocate(hawts_Rot(num_turbines))
    do idx = 1,num_turbines
        call hawts_Rot(idx)%init(inputDir, idx, xG, yG, zG)
    end do 
call message(1,"Done Init")
    ! reset diam, CT
    diam = hawts_Rot(1)%diam
    !CT   = hawts(1)%CT
 
    rhs1 = 0.d0
    call mpi_barrier(mpi_comm_world, ierr)
    call tic()
    do idx = 1,num_turbines
        call hawts_Rot(idx)%get_RHS(u, v, w, rhs1, rhsv, rhsw, inst_val)
    end do 
    call mpi_barrier(mpi_comm_world, ierr)
    call toc()
    call decomp_2d_write_one(1,rhs1,"temp_T1.bin", gp)
    call message(2,"Computed Source fx (T1):", p_sum(sum(rhs1)) * dx*dy*dz)
    call message(2,"Computed Source fy (T1):", p_sum(sum(rhsv)) * dx*dy*dz)
    call message(2,"Computed Source fz (T1):", p_sum(sum(rhsw)) * dx*dy*dz)
    !call message(3,"absolute error:", p_sum(sum(rhs1)) * dx*dy*dz + (num_turbines*0.5d0*(pi/4.d0)*(diam**2)*CT))
    call message(2,"Relative Error  fx (T1):", (one-p_sum(inst_val(1))/(p_sum(sum(rhs1))*dx*dy*dz)))
    call message(2,"Expected Source fy (T1):", zero)
    call message(2,"Expected Source fz (T1):", zero)
    
    deallocate(hawts_Rot)
    deallocate(xG, yG, zG, u, v, w, rhs1, rhs2)
    call MPI_Finalize(ierr)
end program
