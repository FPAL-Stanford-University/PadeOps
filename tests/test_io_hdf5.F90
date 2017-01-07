program test_derivatives_parallel
    use mpi
    use kind_parameters, only : rkind
    use decomp_2d
    use constants,       only: zero, two, pi
    use IO_HDF5_stuff,   only: io_hdf5

    implicit none

    real(rkind), dimension(:,:,:,:), allocatable, target :: coords
    real(rkind), dimension(:,:,:),   allocatable :: f,dfdx, dfdy, dfdz
    real(rkind), dimension(:,:,:),   pointer     :: x, y, z
    type(decomp_info) :: gp
    type(io_hdf5)     :: viz

    integer :: nx = 16, ny = 16, nz = 16
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k

    real(rkind) :: dx, dy, dz

    ! double precision :: t0, t1

    call MPI_Init(ierr)

    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)

    ! Initialize input data (y decomposition)
    allocate( coords    ( gp%ysz(1), gp%ysz(2), gp%ysz(3), 3 ) )
    allocate( f         ( gp%ysz(1), gp%ysz(2), gp%ysz(3)    ) )
    allocate( dfdx      ( gp%ysz(1), gp%ysz(2), gp%ysz(3)    ) )
    allocate( dfdy      ( gp%ysz(1), gp%ysz(2), gp%ysz(3)    ) )
    allocate( dfdz      ( gp%ysz(1), gp%ysz(2), gp%ysz(3)    ) )

    x => coords(:,:,:,1); y => coords(:,:,:,2); z => coords(:,:,:,3)

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
    dfdx = cos(x)*sin(y)*cos(z)
    dfdy = sin(x)*cos(y)*cos(z)
    dfdz = -sin(x)*sin(y)*sin(z)

    print *, 1

    ! Initialize everything
    call viz%init(mpi_comm_world, gp, 'y', '.', 'parallel_hdf5_io')

    call mpi_barrier(mpi_comm_world, ierr)

    print *, 2
    ! Write mesh coordinates to file
    call viz%write_coords(coords)

    call mpi_barrier(mpi_comm_world, ierr)

    print *, 3
    ! Start vizualization dump
    call viz%start_viz(zero)

    call mpi_barrier(mpi_comm_world, ierr)

    print *, 4
    ! Add variables to file
    call viz%write_variable(f,    'f'   )
    call viz%write_variable(dfdx, 'dfdx')
    call viz%write_variable(dfdy, 'dfdy')
    call viz%write_variable(dfdz, 'dfdz')

    call mpi_barrier(mpi_comm_world, ierr)

    print *, 5
    ! End vizualization dump
    call viz%end_viz()

    call mpi_barrier(mpi_comm_world, ierr)

    print *, 6
    deallocate(coords, f, dfdx,dfdy, dfdz)
    nullify(x, y, z)

    call viz%destroy()
    call decomp_2d_finalize
    call MPI_Finalize(ierr)

end program
