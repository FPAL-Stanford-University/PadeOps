program test_averaging
    use mpi
    use decomp_2d,       only: decomp_info, decomp_2d_init, get_decomp_info
    use kind_parameters, only: rkind
    use constants,       only: pi, two
    use AveragingMod,    only: averaging
    use exits,           only: message
    use reductions,      only: P_MAXVAL

    implicit none

    real(rkind), dimension(:,:,:), allocatable :: x, y, z, f, f_avg, f_avg_correct, f_fluct, f_fluct_correct
    type(averaging) :: stats
    type(decomp_info) :: gp

    integer :: nx = 128, ny = 128, nz = 128
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k

    real(rkind), dimension(:), allocatable :: f_local

    real(rkind) :: dx, dy, dz, z_local, f_avg_z, f_avg_z_weighted, z_sum

    call MPI_Init(ierr)

    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)

    ! Initialize input data (y decomposition)
    allocate( x         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( y         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( z         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( f         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )

    dx = two*pi/nx
    dy = two*pi/ny
    dz = two*pi/nz

    ! Generate input data
    do k = 1,gp%ysz(3)
        do j = 1,gp%ysz(2)
            do i = 1,gp%ysz(1)
                x(i,j,k) = real(gp%yst(1) - 1 + i - 1, rkind)*dx
                y(i,j,k) = real(gp%yst(2) - 1 + j - 1, rkind)*dy
                z(i,j,k) = real(gp%yst(3) - 1 + k - 1, rkind)*dz
            end do
        end do
    end do

    ! f = sin(x)*sin(y)*cos(z)
    ! f = (2*x**2 + 3*x + 1) * sin(y) * cos(z)
    ! f = (2*y**2 + 3*y + 1) * sin(x) * cos(z)
    f = (2._rkind*z**2 - 3._rkind*z + 1._rkind) * sin(x) * cos(y)

    allocate(f_local(nz))
    do k = 1,nz
        z_local = real(k - 1, rkind)*dz
        f_local(k) = 2._rkind*z_local**2 - 3._rkind*z_local + 1._rkind
    end do
    f_avg_z = sum(f_local)/real(nz,rkind)

    z_sum = 0._rkind
    do k = 1,nz
        z_local = real(k - 1, rkind)*dz
        f_local(k) = z_local * (2._rkind*z_local**2 - 3._rkind*z_local + 1._rkind)
        z_sum = z_sum + z_local
    end do
    f_avg_z_weighted = sum(f_local)/z_sum
    deallocate(f_local)

    stats =  averaging(gp, 2, [.true., .true., .true.])
    call stats%allocate_average(f_avg)
    call stats%get_average(f, f_avg)
    call message("Max error in average", P_MAXVAL(abs(f_avg)))
    deallocate(f_avg)

    stats =  averaging(gp, 2, [.false., .false., .true.])
    call stats%allocate_average(f_avg)

    call stats%allocate_average(f_avg_correct)

    f_avg_correct(:,:,1) = f_avg_z * sin(x(:,:,1)) * cos(y(:,:,1))
    call stats%get_average(f, f_avg)
    call message("Max error in average", P_MAXVAL(abs(f_avg - f_avg_correct)))

    f_avg_correct(:,:,1) = f_avg_z_weighted * sin(x(:,:,1)) * cos(y(:,:,1))
    call stats%get_weighted_average(z, f, f_avg)
    call message("Max error in average", P_MAXVAL(abs(f_avg - f_avg_correct)))

    allocate( f_fluct        ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( f_fluct_correct( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )

    do k = 1, gp%ysz(3)
        f_fluct_correct(:,:,k) = f(:,:,k) - f_avg_correct(:,:,1)
    end do
    call stats%get_fluctuations(f, f_avg, f_fluct)
    call message("Max error in fluctuations", P_MAXVAL(abs(f_fluct - f_fluct_correct)))

    deallocate( f_fluct        )
    deallocate( f_fluct_correct)


    deallocate(f_avg)
    deallocate(f_avg_correct)

    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(f)

    call mpi_finalize(ierr)
end program
