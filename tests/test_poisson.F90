program test_poisson
    use mpi
    use kind_parameters, only : rkind
    use decomp_2d
    use constants, only: pi, two
    use poisson_eq, only: poisson

    real(rkind), dimension(:,:,:), allocatable :: x, y, z, f,d2fdx2, d2fdy2, rhs, fold
    type(poisson) :: myPoisson
    type(decomp_info) :: gp

    integer :: nx = 1024, ny = 512, nz = 256
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k
    real(rkind) :: maxerr, mymaxerr
    
    logical, parameter :: verbose = .false. 
    real(rkind) :: dx, dy, dz

    double precision :: t0, t1

    call MPI_Init(ierr)
    
    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)
   
    ! Initialize input data (y decomposition) 
    allocate( x         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( y         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( z         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( f         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( d2fdx2    ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( d2fdy2    ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( rhs       ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( fold      ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    
   
    dx = two*pi/nx
    dy = two*pi/ny
    dz = two*pi/nz

    ! Initialize Poisson Solver 
    ierr = myPoisson%init(nx, ny, nz,dx, dy, dz,"four", "cd08", 0)   
    if (nrank == 0) print*, "Initialization complete"
    
    
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

    f = sin(x)*sin(y)
    !dfdx = cos(x)*sin(y)
    d2fdx2 = -sin(x)*sin(y)
    d2fdy2 = -sin(x)*sin(y)

    ! Generate manufactured RHS
    rhs = d2fdx2 + d2fdy2
  
    call mpi_barrier(mpi_comm_world, ierr) 
    t0 = MPI_WTIME()
    call myPoisson%solve(rhs,fold)
    t1 = MPI_WTIME()

    if (verbose) then
        call sleep(nrank)
        print*, "-------------------------------------" 
        print*, "My rank:", nrank
        print*, "x edges:", gp%yst(1), gp%yen(1), gp%ysz(1)
        print*, "y edges:", gp%yst(2), gp%yen(2), gp%ysz(2)
        print*, "z edges:", gp%yst(3), gp%yen(3), gp%ysz(3)
        print*, "-------------------------------------" 
        do i = 1,gp%ysz(1)
            write(*,20) (f(i,j,4), j = 1,gp%ysz(2))
        end do 
        print*, "------------------------------------" 
        do i = 1,gp%ysz(1)
            write(*,20) (fold(i,j,4), j = 1,gp%ysz(2))
        end do 
        print*, "------------------------------------" 
    end if  

    mymaxerr = MAXVAL(ABS(f - fold))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    if (nrank == 0) print*, "Time to solve:", t1 - t0
    call myPoisson%destroy

    deallocate(x, y, z, f, fold, d2fdx2,d2fdy2)
    call decomp_2d_finalize
    call MPI_Finalize(ierr)

20 format(1x,16D10.3)
end program 
