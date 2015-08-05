program test_transpose

    use mpi

    use decomp_2d,       only: decomp_2d_init, decomp_2d_finalize, get_decomp_info, &
                               decomp_info,decomp_info_init,decomp_info_finalize,   &
                               transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y, &
                               nrank, real_type

    use kind_parameters, only: rkind
    use constants,       only: one,two,pi
    use timer,           only: tic, toc

    use cd10stuff,       only: cd10
    
    implicit none

    type(decomp_info) :: gp
    type(cd10)        :: xcd10, ycd10, zcd10

    integer :: i, j, k, ierr

    integer :: nx = 128, ny = 128, nz = 128
    integer :: prow = 0, pcol = 0

    logical :: periodicx = .TRUE.
    logical :: periodicy = .TRUE.
    logical :: periodicz = .TRUE.

    real(rkind) :: omega = one
    real(rkind) :: dx, dy, dz
    real(rkind), dimension(:,:,:), allocatable :: x,f,dfdx,dfdx_exact, dfdy, dfdz, dummy
    
    real(rkind), dimension(:,:,:), allocatable :: wkx, dwkx, wkz, dwkz

    real(rkind) :: mymaxerr, maxerr

    double precision :: t0, t1, xtime, ytime, ztime, dummytime, time2trans_x, time2compute_x

    call MPI_Init(ierr)

    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)
    
    ! Create X work arrays
    allocate(  wkx( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( dwkx( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )

    ! Allocate all physical data in Y pencils
    allocate( x         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( f         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( dfdx      ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( dfdy      ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( dfdz      ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( dummy     ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( dfdx_exact( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )

    ! Create Z work arrays
    allocate(  wkz( gp%zsz(1), gp%zsz(2), gp%zsz(3) ) )
    allocate( dwkz( gp%zsz(1), gp%zsz(2), gp%zsz(3) ) )

    dx = two*pi/real(nx,rkind)
    dy = dx
    dz = dx

    ierr = xcd10%init( nx, dx, periodicx, 0, 0)
    ierr = ycd10%init( ny, dy, periodicy, 0, 0)
    ierr = zcd10%init( nz, dz, periodicz, 0, 0)

    do k=1,gp%ysz(3)
        do j=1,gp%ysz(2)
            do i=1,gp%ysz(1)
                x(i,j,k) = real( gp%yst(1)-1+i-1,rkind )*dx
                f(i,j,k) = sin(omega*x(i,j,k))
                dfdx_exact(i,j,k) = omega*cos(omega*x(i,j,k))
            end do
        end do
    end do

    t0 = MPI_WTIME()
    !======== Start the X derivative stuff ========!
    call transpose_y_to_x(f, wkx, gp)
    t1 = MPI_WTIME()
    time2trans_x = t1 - t0
    call xcd10 % dd1( wkx, dwkx, gp%xsz(2), gp%xsz(3) )
    t0 = MPI_WTIME()
    time2compute_x = t0 - t1
    call transpose_x_to_y(dwkx, dfdx, gp)
    !======== End the X derivative stuff   ========!
    t1 = MPI_WTIME()
    if(nrank == 0) xtime = t1 - t0 + time2compute_x  + time2trans_x
    
    mymaxerr = MAXVAL(ABS(dfdx - dfdx_exact))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    
    t0 = MPI_WTIME()
    !======== Start the Y derivative stuff ========!
    call ycd10 % dd2( f, dfdy, gp%ysz(1), gp%ysz(3) )
    !======== End the Y derivative stuff   ========!
    t1 = MPI_WTIME()
    if(nrank == 0) ytime = t1-t0

    mymaxerr = MAXVAL(ABS(dfdy))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    
    t0 = MPI_WTIME()
    !======== Start the Z derivative stuff ========!
    call transpose_y_to_z(f, wkz, gp)
    call zcd10 % dd3( wkz, dwkz, gp%zsz(1), gp%zsz(2) )
    call transpose_z_to_y(dwkz, dfdz, gp)
    !======== End the Z derivative stuff   ========!
    t1 = MPI_WTIME()
    if(nrank == 0) ztime = t1-t0
   

    mymaxerr = MAXVAL(ABS(dfdz))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    
    t0 = MPI_WTIME()
    !======== Perform some dummy operations========!
    dummy = sin(x)
    
    dummy = dummy + dfdx*dfdx*f + dfdy*dfdx*x + dfdz*dfdx*dfdx_exact + f + x

    t1 = MPI_WTIME()
    if (nrank == 0) dummytime = t1 - t0 

    if (nrank == 0) then
        print*, "Time for X derivative: ", xtime
        print*, "            transpose: ", time2trans_x 
        print*, "            compute  : ", time2compute_x 
        print*, "Time for Y derivative: ", ytime
        print*, "Time for Z derivative: ", ztime
        print*, "Time for dummy calc  : ", dummytime
        print*, "----------------------------------------------------"
        print*, "           Total time: ", xtime+ytime+ztime+dummytime
    end if
    
    call xcd10%destroy()
    call ycd10%destroy()
    call zcd10%destroy()
    
    deallocate( x          )
    deallocate( f          )
    deallocate( dfdx       )
    deallocate( dfdx_exact )

    call decomp_2d_finalize

    call MPI_Finalize(ierr)

end program
