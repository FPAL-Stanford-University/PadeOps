program test_actuatorDisk
    use kind_parameters, only: rkind, clen
    use constants, only: pi, two, one, imi, zero, half, kappa
    use reductions, only: p_maxval, p_sum 
    use timer, only: tic, toc
    use decomp_2d
    use decomp_2d_io
    use actuatorDisk_T2Mod, only: actuatorDisk_T2
    use actuatorDiskMod, only: actuatorDisk
    use actuatorDisk_YawMod, only: actuatorDisk_yaw
    use exits, only: message
    use mpi

    implicit none 

    type(actuatorDisk_T2), dimension(:), allocatable :: hawts_T2
    type(actuatorDisk), dimension(:), allocatable :: hawts
    type(actuatorDisk_yaw), dimension(:), allocatable :: hawts_Tyaw
    integer, parameter :: nx = 320, ny = 240, nz = 192
    !character(len=clen) :: inputDir = "/home/aditya90/Codes/PadeOps/data/AD_Coriolis/"
    !character(len=clen) :: inputDir = "/fastscratch/nghaisas/runs/PadeOps/wupa/run3/deb/turbInfo"
    character(len=clen) :: inputDir = "/home1/05294/mhowland/PadeOps/problems/turbines/neutral_pbl_concurrent_files/turbInfo/6x1array"
    real(rkind), dimension(:,:,:), allocatable :: xG, yG, zG
    real(rkind), dimension(:,:,:), allocatable :: u, v, w, rhs1, rhsv, rhsw, rhs2, rhs3
    !real(rkind), parameter :: Lx = pi, Ly = 0.5d0*pi, Lz = 1.0d0
    real(rkind), parameter :: Lx = 30.d0, Ly = 15.d0, Lz = 6.0d0
    real(rkind) :: dx, dy, dz, diam, CT
    type(decomp_info) :: gp 
    integer :: idx, ix1, iy1, iz1, ixn, iyn, izn, i, j, k, ierr, prow = 0, pcol = 0, num_turbines 
    real(rkind) :: xPeriods = 2.d0, yPeriods = 2.d0, zpeak = 0.3d0, epsnd = 5.d0, z0init = 1.d-4 
    real(rkind) :: inst_val(8)
    real(rkind) :: comp1, comp2, maxdiff
    real(rkind) :: gamma_negative, theta
    real(rkind), dimension(:,:,:), allocatable :: rbuff, blanks, speed
    real(rkind), dimension(:,:,:), allocatable :: scalarSource

    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)

    allocate(xG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(yG(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(zG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(u(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(v(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(w(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(rhs1(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(rhs2(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(rhs3(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(rhsv(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(rhsw(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    ! Yaw stuff
    allocate(rbuff(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(blanks(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(speed(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(scalarSource(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 

    num_turbines = 6

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
    allocate(hawts(num_turbines))
    allocate(hawts_T2(num_turbines))
    allocate(hawts_Tyaw(num_turbines))
    do idx = 1,num_turbines
        call hawts(idx)%init(inputDir, idx, xG, yG, zG, gp)
        call hawts_T2(idx)%init(inputDir, idx, xG, yG, zG)
        call hawts_Tyaw(idx)%init(inputDir, idx, xG, yG, zG)
        call hawts_Tyaw(idx)%link_memory_buffers(rbuff, blanks, speed, scalarSource)
    end do 

    ! reset diam, CT
    diam = hawts(1)%diam
    CT   = hawts(1)%CT
 
    rhs1 = 0.d0
    call mpi_barrier(mpi_comm_world, ierr)
    call tic()
    do idx = 1,num_turbines
        call hawts(idx)%get_RHS(u, v, w, rhs1, rhsv, rhsw, inst_val)
    end do 
    call mpi_barrier(mpi_comm_world, ierr)
    call toc()
    call decomp_2d_write_one(1,rhs1,"temp_T1.bin", gp)
    call message(2,"Computed Source (T1):", p_sum(sum(rhs1)) * dx*dy*dz)
    call message(3,"absolute error:", p_sum(sum(rhs1)) * dx*dy*dz + (num_turbines*0.5d0*(pi/4.d0)*(diam**2)*CT))
    
    rhs2 = 0.d0
    call mpi_barrier(mpi_comm_world, ierr)
    call tic()
    do idx = 1,num_turbines
        call hawts_T2(idx)%get_RHS(u, v, w, rhs2, rhsv, rhsw, inst_val)
    end do 
    call mpi_barrier(mpi_comm_world, ierr)
    call toc()
    
    rhs3 = 0.d0
    call mpi_barrier(mpi_comm_world, ierr)
    call tic()
    gamma_negative = 0.d0 * pi / 180.d0
    theta = 0.d0
    do idx = 1,num_turbines
        call hawts_Tyaw(idx)%get_RHS(u, v, w, rhs3, rhsv, rhsw, gamma_negative, theta)
    end do 
    call mpi_barrier(mpi_comm_world, ierr)
    call toc()
    
    
    call decomp_2d_write_one(1,rhs2,"temp_T2.bin", gp)
    call message(2,"Computed Source (T2):", p_sum(sum(rhs2)) * dx*dy*dz)
    call message(3,"absolute error:", p_sum(sum(rhs2)) * dx*dy*dz + (num_turbines*0.5d0*(pi/4.d0)*(diam**2)*CT))
    call message(2,"Expected Source:", -(num_turbines*0.5d0*(pi/4.d0)*(diam**2)*CT))

    call decomp_2d_write_one(1,rhs3,"temp_T4.bin", gp)
    call message(2,"Computed Source (Tyaw):", p_sum(sum(rhs3)) * dx*dy*dz)
    call message(3,"absolute error:", p_sum(sum(rhs3)) * dx*dy*dz + (num_turbines*0.5d0*(pi/4.d0)*(diam**2)*CT))
    call message(2,"Expected Source:", -(num_turbines*0.5d0*(pi/4.d0)*(diam**2)*CT))
    
    maxDiff = p_maxval(abs(rhs2 - rhs1))*dx*dy*dz
    call message(2,"RHS AD Type2 - Type 1: ", maxDiff)
    maxDiff = p_maxval(abs(rhs3 - rhs1))*dx*dy*dz
    call message(2,"RHS AD Type3 - Type 1: ", maxDiff)
    do idx = 1,num_turbines
    !  call hawts(idx)%destroy()
    end do 
    deallocate(hawts)
    deallocate(xG, yG, zG, u, v, w, rhs1, rhs2, rhs3)
    deallocate(rbuff, blanks, speed, scalarSource)
    call MPI_Finalize(ierr)
end program 
