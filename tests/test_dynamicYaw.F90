program test_dynamicYaw
    use kind_parameters, only: rkind, clen
    use constants, only: pi, two, one, imi, zero, half, kappa
    use reductions, only: p_maxval, p_sum 
    use timer, only: tic, toc
    use decomp_2d
    use decomp_2d_io
    use dynamicYawMod, only: dynamicYaw
    use actuatorDisk_YawMod, only: actuatorDisk_yaw
    use exits, only: message
    use mpi

    implicit none 

    type(dynamicYaw) :: dyaw
    type(actuatorDisk_yaw), allocatable, dimension(:) :: ad
    character(len=clen) :: inputDir = "/home1/05294/mhowland/dynamicYawFiles/dynamicYaw.inp"
    character(len=clen) :: inputDir_turb = "/home1/05294/mhowland/PadeOps/problems/turbines/neutral_pbl_concurrent_files/turbInfo/2x1array"
    integer, parameter :: nx = 192, ny = 96, nz = 128
    real(rkind), dimension(:,:,:), allocatable :: xG, yG, zG
    real(rkind), parameter :: Lx = 2.d0, Ly = 2.d0, Lz = 2.0d0
    real(rkind) :: dx, dy, dz, diam, CT
    integer :: idx, ix1, iy1, iz1, ixn, iyn, izn, i, j, k, ierr, prow = 0, pcol = 0, num_turbines 
    type(decomp_info) :: gp 
    real(rkind), dimension(:,:,:), allocatable :: rbuff, blanks, speed, X
    real(rkind), dimension(:,:,:), allocatable :: Y, Z, scalarSource
    real(rkind), dimension(:), allocatable :: yaw

    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)

    ! number of turbines
    num_turbines = 2

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
    allocate(yaw(num_turbines)) 
         
   do i = 1, num_turbines
        call ad(i)%init(inputDir_turb, i, xG, yG, zG)
        call ad(i)%link_memory_buffers(rbuff, blanks, speed, X, &
                            Y, Z, scalarSource)
    end do
    call message(0,"YAWING WIND TURBINE (Type 4) array initialized")


    ! Initialize
    call dyaw%init(inputDir, ad)

    ! Run the full dynamic yaw state estimation and yaw optimize
    yaw = 10.d0 * pi / 180.d0
    ! Input
    write(*,*) yaw*180.d0/pi
    call dyaw%update_and_yaw(yaw)

    ! Output
    write(*,*) dyaw%yaw*180.d0/pi

    call MPI_Finalize(ierr)
end program 
