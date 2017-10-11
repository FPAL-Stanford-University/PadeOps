program test_io_hdf5
    use mpi
    use kind_parameters, only : rkind, single_kind
    use decomp_2d
    use decomp_2d_io
    use constants,       only: eps, two, pi
    use exits,           only: message
    use reductions,      only: P_MAXVAL
    use io_hdf5_stuff,   only: io_hdf5

    implicit none

    real(rkind), dimension(:,:,:,:), allocatable, target :: coords, coords_read
    real(rkind), dimension(:,:,:),   allocatable :: f,dfdx, dfdy, dfdz
    real(rkind), dimension(:,:,:),   pointer     :: x, y, z
    type(decomp_info) :: gp
    type(io_hdf5)     :: viz

    integer :: nx = 128, ny = 192, nz = 129
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k, iter

    real(rkind) :: dx, dy, dz
    real(rkind) :: time
    real(rkind), dimension(1) :: time_arr
    integer, dimension(3) :: subdomain_lo, subdomain_hi

    logical :: reduce_precision = .true.
    real(rkind) :: tolerance = 10._rkind*eps

    double precision :: t0, t1, t2, t3, t4
    real(rkind) :: gbytes
    integer :: elemsize

    call MPI_Init(ierr)

    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)

    ! Initialize input data (y decomposition)
    allocate( coords     ( gp%ysz(1), gp%ysz(2), gp%ysz(3), 3 ) )
    allocate( coords_read( gp%ysz(1), gp%ysz(2), gp%ysz(3), 3 ) )
    allocate( f          ( gp%ysz(1), gp%ysz(2), gp%ysz(3)    ) )
    allocate( dfdx       ( gp%ysz(1), gp%ysz(2), gp%ysz(3)    ) )
    allocate( dfdy       ( gp%ysz(1), gp%ysz(2), gp%ysz(3)    ) )
    allocate( dfdz       ( gp%ysz(1), gp%ysz(2), gp%ysz(3)    ) )

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

    ! subdomain_lo = [      1,      1,      1]
    ! subdomain_hi = [     nx,     ny,     nz]
    subdomain_lo = [   nx/4,   ny/4,   nz/4]
    subdomain_hi = [ 3*nx/4, 3*ny/4, 3*nz/4]

    if (nrank == 0) then
        print *, ""
        print '(A,9(I0,A))', "Writing subdomain [",subdomain_lo(1),":",subdomain_hi(1),",",subdomain_lo(2),":",subdomain_hi(2),",",subdomain_lo(3),":",subdomain_hi(3),"] of ",nx,"x",ny,"x",nz," grid"
        print *, ""
    end if

    if (reduce_precision) then
        elemsize = storage_size(1._single_kind)
        tolerance = 10._rkind*epsilon(real(1.0,single_kind))
    else
        elemsize = storage_size(1._rkind)
        tolerance = 10._rkind*epsilon(real(1.0,rkind))
    end if
    gbytes = ( real(elemsize, rkind) * product(subdomain_hi - subdomain_lo + 1) * 4 ) / ( real(8, rkind) * 1024 * 1024 * 1024 )

    ! Initialize everything
    call viz%init(mpi_comm_world, gp, 'y', '.', 'parallel_hdf5_io', reduce_precision=reduce_precision, read_only=.false., &
                  subdomain_lo=subdomain_lo, subdomain_hi=subdomain_hi)

    ! Write mesh coordinates to file
    call viz%write_coords(coords)

    do iter=1,4
        call message("Starting viz dump", viz%vizcount)

        time = (iter - 1)*0.01D0
        f    =  sin(iter*x)*sin(iter*y)*cos(iter*z)*iter
        dfdx =  cos(iter*x)*sin(iter*y)*cos(iter*z)*iter
        dfdy =  sin(iter*x)*cos(iter*y)*cos(iter*z)*iter
        dfdz = -sin(iter*x)*sin(iter*y)*sin(iter*z)*iter

        call mpi_barrier(mpi_comm_world, ierr)
        t0 = mpi_wtime()

        ! Start vizualization dump
        call viz%start_viz(time)

        call mpi_barrier(mpi_comm_world, ierr)
        t1 = mpi_wtime()

        ! Add variables to file
        call viz%write_variable(f,    'f'   )
        call viz%write_variable(dfdx, 'dfdx')
        call viz%write_variable(dfdy, 'dfdy')
        call viz%write_variable(dfdz, 'dfdz')

        call mpi_barrier(mpi_comm_world, ierr)
        t2 = mpi_wtime()

        ! End vizualization dump
        call viz%end_viz()

        call mpi_barrier(mpi_comm_world, ierr)
        t3 = mpi_wtime()

        call decomp_2d_write_one(2, f,    'f.dat',    gp)
        call decomp_2d_write_one(2, dfdx, 'dfdx.dat', gp)
        call decomp_2d_write_one(2, dfdy, 'dfdy.dat', gp)
        call decomp_2d_write_one(2, dfdz, 'dfdz.dat', gp)

        call mpi_barrier(mpi_comm_world, ierr)
        t4 = mpi_wtime()

        if (nrank == 0) then
            print '(A,ES10.3,A)', "    Time to start viz       = ", t1-t0, " seconds"
            print '(A,ES10.3,A)', "    Time to write variables = ", t2-t1, " seconds"
            print '(A,ES10.3,A)', "    Time to end viz         = ", t3-t2, " seconds"
            print '(A,ES10.3,A)', "    Effective I/O bandwidth = ", gbytes/(t2-t1), " GB/s"
            print '(A,ES10.3,A)', ""
            print '(A,ES10.3,A)', "    Time to write variables using 2decomp = ", t4-t3, " seconds"
            print '(A,ES10.3,A)', "    Effective I/O bandwidth = ", gbytes/(t4-t3), " GB/s"
            print '(A,ES10.3,A)', ""
        end if
    end do

    call viz%destroy()
    
    call message("")
    call message("Now reading in the file written out")

    call viz%init(mpi_comm_world, gp, 'y', '.', 'parallel_hdf5_io', reduce_precision=reduce_precision, read_only=.true., &
                  subdomain_lo=subdomain_lo, subdomain_hi=subdomain_hi)

    iter = viz%find_last_dump()
    if (iter /= 3) then
        call message("ERROR:")
        call message("  Last dump index is incorrect")
        stop
    end if

    call viz%read_coords(coords_read)
    do k = 1,gp%ysz(3)
        if ( ( (gp%yst(3) - 1 + k) < subdomain_lo(3) ) .or. ( (gp%yst(3) - 1 + k) > subdomain_hi(3) ) ) cycle
        do j = 1,gp%ysz(2)
            if ( ( (gp%yst(2) - 1 + j) < subdomain_lo(2) ) .or. ( (gp%yst(2) - 1 + j) > subdomain_hi(2) ) ) cycle
            do i = 1,gp%ysz(1)
                if ( ( (gp%yst(1) - 1 + i) < subdomain_lo(1) ) .or. ( (gp%yst(1) - 1 + i) > subdomain_hi(1) ) ) cycle
                if ( maxval(abs(coords_read(i,j,k,:) - coords(i,j,k,:))) > tolerance ) then
                    call message("ERROR:")
                    call message("  Coordinates read in have incorrect values")
                    stop
                end if
            end do
        end do
    end do


    do iter = 1,4
        f    =  sin(iter*x)*sin(iter*y)*cos(iter*z)*iter
        time = (iter - 1)*0.01D0
        call viz%start_reading(iter-1)
        call viz%read_dataset(dfdx, 'f')
        
        do k = 1,gp%ysz(3)
            if ( ( (gp%yst(3) - 1 + k) < subdomain_lo(3) ) .or. ( (gp%yst(3) - 1 + k) > subdomain_hi(3) ) ) cycle
            do j = 1,gp%ysz(2)
                if ( ( (gp%yst(2) - 1 + j) < subdomain_lo(2) ) .or. ( (gp%yst(2) - 1 + j) > subdomain_hi(2) ) ) cycle
                do i = 1,gp%ysz(1)
                    if ( ( (gp%yst(1) - 1 + i) < subdomain_lo(1) ) .or. ( (gp%yst(1) - 1 + i) > subdomain_hi(1) ) ) cycle
                        if ( abs(f(i,j,k) - dfdx(i,j,k)) > tolerance ) then
                            call message("ERROR:")
                            call message("  Array 'f' read in has incorrect values")
                            stop
                        end if
                    end do
                end do
            end do

        call viz%read_attribute(1, time_arr, 'Time')
        if ( abs(time_arr(1) - time) > 10.D0*eps ) then
            call message("ERROR:")
            call message("  Wrong value of time read in. Expected time",time)
            call message("  Got time",time_arr(1))
        end if
        call viz%end_reading
    end do

    deallocate(coords, coords_read, f, dfdx,dfdy, dfdz)
    nullify(x, y, z)

    call viz%destroy()
    call decomp_2d_finalize
    call MPI_Finalize(ierr)

end program
