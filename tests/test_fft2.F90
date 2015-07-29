program test_fft2d
    use mpi
    use kind_parameters, only : rkind
    use decomp_2d
    use fft_2d_stuff, only: fft_2d 
    use constants, only: pi, two

    real(rkind), dimension(:,:,:), allocatable :: x, y, z, f,d2fdx2, d2fdy2, fold
    type(fft_2d) :: myfft2d
    type(decomp_info) :: gp
    complex(rkind), dimension(:,:,:), allocatable :: fhat

    integer :: nx = 8, ny = 8, nz = 8
    integer :: prow = 2, pcol = 2
    integer :: ierr, i, j, k
    real(rkind) :: maxerr, mymaxerr
    
    logical :: get_exhaustive_plan = .true.
    logical, parameter :: verbose = .false. 
    real(rkind) :: dx, dy

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
    allocate( fold      ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    
   
    dx = two*pi/nx
    dy = two*pi/ny

    ! Initialize fft_2d type  
    ierr = myfft2d%init(nx,ny,nz,"y",dx, dy, get_exhaustive_plan)
    call myfft2d%alloc_output(fhat) 

    if (nrank == 0) print*, "Initialization complete"
    ! Generate input data

    do k = 1,gp%ysz(3)
        do j = 1,gp%ysz(2)
            do i = 1,gp%ysz(1)
                x(i,j,k) = real(gp%yst(1) - 1 + i - 1, rkind)*dx
                y(i,j,k) = real(gp%yst(2) - 1 + j - 1, rkind)*dy 
            end do 
        end do 
    end do 

    f = sin(x)*sin(y)
    !dfdx = cos(x)*sin(y)
    d2fdx2 = -sin(x)*sin(y)
    d2fdy2 = -sin(x)*sin(y)

   
    call MPI_BARRIER(mpi_comm_world,ierr)
    t0 = MPI_WTIME()
    ! Forward transform 
    call myfft2d%fft2(f,fhat)
    t1 = MPI_WTIME()
    
    if (nrank == 0) print*, "Forward transform time:", t1 - t0

    do k = 1,gp%ysz(3)
        fhat(:,:,k) = -myfft2d%k2 * myfft2d%k2 * fhat(:,:,k)
    end do 

    t0 = MPI_WTIME()
    ! Backward transform 
    call myfft2d%ifft2(fhat,fold)
    t1 = MPI_WTIME()

    if (nrank == 0) print*, "Backward transform time:", t1 - t0
    
    if (verbose) then
        call sleep(nrank)
        print*, "-------------------------------------" 
        print*, "My rank:", nrank
        print*, "x edges:", gp%yst(1), gp%yen(1), gp%ysz(1)
        print*, "y edges:", gp%yst(2), gp%yen(2), gp%ysz(2)
        print*, "z edges:", gp%yst(3), gp%yen(3), gp%ysz(3)
        print*, "-------------------------------------" 
            do i = 1,gp%ysz(1)
                write(*,20) (d2fdy2(i,j,4), j = 1,gp%ysz(2))
            end do 
        print*, "------------------------------------" 
        print*, "Recalcuted:"
        print*, "------------------------------------" 
            do i = 1,gp%ysz(1)
                write(*,20) (fold(i,j,4), j = 1,gp%ysz(2))
            end do 
        
    end if 
   
    
    mymaxerr = MAXVAL(ABS(d2fdy2 - fold))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    
    
    
    call myfft2d%destroy


    deallocate(x, y, z, f, fold, fhat, d2fdx2,d2fdy2)
    call decomp_2d_finalize
    call MPI_Finalize(ierr)
    
20 format(1x,16D10.3)
end program 
