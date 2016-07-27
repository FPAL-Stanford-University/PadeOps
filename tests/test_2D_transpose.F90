program test_transpose

    use mpi

    use decomp_2d,       only: decomp_2d_init, decomp_2d_finalize, get_decomp_info, &
                               decomp_info,decomp_info_init,decomp_info_finalize,   &
                               transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y, &
                               nrank, real_type

    use kind_parameters, only: rkind
    use constants,       only: one,two,pi
    use timer,           only: tic, toc
    use reductions,      only: p_maxval
    use DerivativesMod,  only: derivatives
    
    implicit none

    type(decomp_info) :: gp
    type(derivatives) :: der

    integer :: i, j, k, ierr

    integer :: nx = 512, ny = 512, nz = 512
    integer :: prow = 0, pcol = 0

    logical :: periodicx = .TRUE.
    logical :: periodicy = .TRUE.
    logical :: periodicz = .TRUE.

    real(rkind) :: omega = one
    real(rkind) :: dx, dy, dz
    real(rkind), dimension(:,:,:), allocatable :: x,f,dfdx,dfdx_exact, dfdy, dfdz, dummy
    real(rkind), dimension(:,:,:), allocatable :: buffX, buffY, buffZ

    real(rkind), dimension(:,:,:), allocatable :: wkx, dwkx, wkz, dwkz

    real(rkind) :: mymaxerr, maxerr

    double precision :: t0, t1, xtime, time2trans_x, time2compute_x

    call MPI_Init(ierr)

    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)
    
    ! Create X work arrays
    allocate(  wkx( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( dwkx( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )

    allocate(  buffx( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate(  buffy( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate(  buffz( gp%zsz(1), gp%zsz(2), gp%zsz(3) ) )
    
    
    ! Allocate all physical data in Y pencils
    allocate( x         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( f         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( dfdx      ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( dfdy      ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( dfdz      ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( dummy     ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( dfdx_exact( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )

    dx = two*pi/real(nx,rkind)
    dy = dx
    dz = dx

    do k=1,gp%ysz(3)
        do j=1,gp%ysz(2)
            do i=1,gp%ysz(1)
                x(i,j,k) = real( gp%yst(1)-1+i-1,rkind )*dx
                f(i,j,k) = sin(omega*x(i,j,k))
                dfdx_exact(i,j,k) = omega*cos(omega*x(i,j,k))
            end do
        end do
    end do

    call der%init(                         gp, &
                            dx,     dy,    dz, &
                        .TRUE., .TRUE.,.TRUE., &
                        "cd10", "cd10", "cd10" )

    t0 = MPI_WTIME()
    !======== Start the X derivative stuff ========!
    call transpose_y_to_x(f, wkx, gp)
    t1 = MPI_WTIME()
    time2trans_x = t1 - t0
    call der % ddx( wkx, dwkx )
    t0 = MPI_WTIME()
    time2compute_x = t0 - t1
    call transpose_x_to_y(dwkx, dfdx, gp)
    !======== End the X derivative stuff   ========!
    t1 = MPI_WTIME()
    if(nrank == 0) xtime = t1 - t0 + time2trans_x + time2compute_x

    mymaxerr = MAXVAL(ABS(dfdx - dfdx_exact))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    
    
    buffx = 2.5_rkind
    if (nrank == 0) print*, "Now Testing Transposes"
    call mpi_barrier(mpi_comm_world, ierr)
    call tic()
    t0 = MPI_WTIME()
    call transpose_x_to_y(buffx,buffy, gp)
    call toc()
    call tic()
    call transpose_y_to_z(buffy,buffz, gp)
    call toc()
    call tic()
    call transpose_z_to_y(buffz,buffy, gp)
    call toc()
    call tic()
    call transpose_y_to_x(buffy,buffx, gp)
    t1 = MPI_WTIME()
    call toc()
    t1 = t1 - t0
    t0 = p_maxval(t1)
    if (nrank == 0) print*, "Total transpose time:", t0
    deallocate( x          )
    deallocate( f          )
    deallocate( dfdx       )
    deallocate( dfdx_exact )

    call der%destroy

    call decomp_2d_finalize

    call MPI_Finalize(ierr)

end program
