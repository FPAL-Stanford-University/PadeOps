program test_statistics
    use mpi
    use decomp_2d
    use kind_parameters, only: rkind
    use constants,       only: pi, two
    use StatisticsMod,   only: statistics
    use exits,           only: message
    use reductions,      only: P_MAXVAL

    implicit none

    real(rkind), dimension(:,:,:), allocatable :: x, y, z, f, f_avg
    type(statistics) :: stats
    type(decomp_info) :: gp

    integer :: nx = 128, ny = 128, nz = 128
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k

    real(rkind) :: dx, dy, dz

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

    f = sin(x)*sin(y)*cos(z)

    stats =  statistics(gp, 1, [.true., .true., .true.])
    call stats%allocate_average(f_avg)
    call stats%get_average(f, f_avg)
    call message("Max error in average", P_MAXVAL(abs(f_avg)))
    deallocate(f_avg)

    stats =  statistics(gp, 1, [.false., .false., .true.])
    call stats%allocate_average(f_avg)
    call stats%get_average(f, f_avg)
    call message("Max error in average", P_MAXVAL(abs(f_avg)))
    deallocate(f_avg)

    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(f)

    call mpi_finalize(ierr)
end program
