program test_dynamicYawUncertain
    use kind_parameters, only: rkind, clen
    use constants, only: pi, two, one, imi, zero, half, kappa
    use Gridtools, only: linspace, mytrapz
    use reductions, only: p_maxval, p_sum 
    use timer, only: tic, toc
    use decomp_2d
    use decomp_2d_io
    use dynamicYawMod, only: dynamicYaw
    use actuatorDisk_YawMod, only: actuatorDisk_yaw
    use exits, only: message
    use mpi
    use basic_io, only: read_2d_ascii

    implicit none 

    type(dynamicYaw) :: dyaw
    type(actuatorDisk_yaw), allocatable, dimension(:) :: ad
    character(len=clen) :: inputDir = "/home1/05294/mhowland/dynamicYawFiles/WES_P2_2020_inputs/neutral/dynamicYaw_neutral_debug.inp"
    character(len=clen) :: inputDir_turb = "/home1/05294/mhowland/PadeOps_diurnal/PadeOps/problems/incompressible/diurnal_concurrent_files/turbInfo/4x2array_0"
    character(len=clen) :: input_LUT = "/home1/05294/mhowland/PadeOps_diurnal/PadeOps/problems/incompressible/diurnal_concurrent_files/turbInfo/yawLookupTables/4x2array_0/empirical.txt"

    integer, parameter :: nx = 192, ny = 96, nz = 128
    real(rkind), dimension(:,:,:), allocatable :: xG, yG, zG
    real(rkind), parameter :: Lx = 2.d0, Ly = 2.d0, Lz = 2.0d0
    real(rkind) :: dx, dy, dz, diam, CT, wind_speed, wind_direction, dirStd
    integer :: idx, ix1, iy1, iz1, ixn, iyn, izn, i, j, k, ierr, prow = 0, pcol = 0, num_turbines, dynamicStart, dirType
    type(decomp_info) :: gp 
    real(rkind), dimension(:,:,:), allocatable :: rbuff, blanks, speed, X
    real(rkind), dimension(:,:,:), allocatable :: Y, Z, scalarSource
    real(rkind), dimension(:), allocatable :: xLoc, yLoc, power, Pstd
    logical :: fixedYaw
    logical :: considerAdvection, lookup
    ! Variables for the alpha_check() function
    real(rkind), dimension(4000) :: alpha, t
    real(rkind) :: alpha_m, alpha_std
    integer :: Tf_out
    ! Lookup table
    real(rkind), dimension(:,:), allocatable :: data2read, yaw_LUT
    real(rkind), dimension(:), allocatable :: alpha_LUT, diffv, yaw_setpoints, yaw
    real(rkind) :: alpha_input
    integer :: alpha_index

    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)

    ! number of turbines
    num_turbines = 9

    ! AD stuff
    allocate(xG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(yG(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(zG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); 
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
    allocate(ad(num_turbines))

    !! Yaw stuff
    allocate(rbuff(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(blanks(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(speed(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(X(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(Y(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(Z(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(scalarSource(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    !allocate(yaw(num_turbines)) 
    allocate(xLoc(num_turbines)) 
    allocate(yLoc(num_turbines)) 
    allocate(power(num_turbines)) 
    allocate(Pstd(num_turbines)) 
         
   do i = 1, num_turbines
        call ad(i)%init(inputDir_turb, i, xG, yG, zG)
        call ad(i)%link_memory_buffers(rbuff, blanks, speed, scalarSource)
    end do
    call message(0,"YAWING WIND TURBINE (Type 4) array initialized")


    ! Initialize
    do i=1,num_turbines
        xLoc(i) = ad(i)%xLoc
        yLoc(i) = ad(i)%yLoc
    end do
    call dyaw%init(inputDir, xLoc, yLoc, 0.315d0, num_turbines, fixedYaw, dynamicStart, dirType, considerAdvection, lookup)
    !write(*,*) dyaw%turbCenter(:,1)
    !write(*,*) dyaw%turbCenter(:,2)

    ! Lookup table stuff
    allocate(data2read(17, num_turbines))
    allocate(yaw_LUT(17, num_turbines-1))
    allocate(alpha_LUT(17))
    allocate(diffv(17))
    allocate(yaw_setpoints(num_turbines-1))
    allocate(yaw(num_turbines))

    ! Load lookup table
    call read_2d_ascii(data2read, input_LUT)
    alpha_LUT = data2read(:,1)
    yaw_LUT = data2read(:,2:num_turbines)

    ! Apply lookup table
    alpha_input = 7.5d0
    diffv = abs(alpha_input - alpha_LUT) 
    alpha_index = minloc(diffv, 1)
    yaw_setpoints = yaw_LUT(alpha_index,:)

    !write(*,*) yaw_setpoints

    ! Create yaw setpoints with reference
    yaw(1) = 0.d0
    yaw(2:num_turbines) = yaw_setpoints
    write(*,*) yaw

    call MPI_Finalize(ierr)
end program 
