program test_fft3d
    use mpi
    use kind_parameters, only : rkind
    use decomp_2d
    use fft_3d_stuff, only: fft_3d 
    use constants, only: pi, two

    implicit none
    real(rkind), dimension(:,:,:), allocatable :: x, y, z, f,d2fdx2, d2fdy2, d2fdz2,fold
    type(fft_3d) :: myfft3d
    type(decomp_info) :: gp
    complex(rkind), dimension(:,:,:), allocatable :: fhat

    integer :: nx = 64, ny = 64, nz =64
    integer :: prow = 2, pcol = 2
    integer :: ierr, i, j, k
    real(rkind) :: maxerr, mymaxerr
    
    logical :: get_exhaustive_plan = .true.
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
    allocate( d2fdz2    ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( fold      ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    
   
    dx = two*pi/nx
    dy = two*pi/ny
    dz = two*pi/nz

    ! Initialize fft_2d type  
    ierr = myfft3d%init(nx,ny,nz,"y",dx, dy,dz, get_exhaustive_plan)
    call myfft3d%alloc_output(fhat) 

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

    f = sin(x)*sin(y)*sin(z)
    d2fdx2 = -sin(x)*sin(y)*sin(z)
    d2fdy2 = -sin(x)*sin(y)*sin(z)
    d2fdz2 = -sin(x)*sin(y)*sin(z)
   
    call MPI_BARRIER(mpi_comm_world,ierr)
    t0 = MPI_WTIME()
    ! Forward transform 
    call myfft3d%fft3(f,fhat)
    t1 = MPI_WTIME()
    
    if (nrank == 0) print*, "Forward transform time:", t1 - t0

    ! Compute the laplacian as a check 
    fhat = -fhat * myfft3d%kabs_sq

    t0 = MPI_WTIME()
    ! Backward transform 
    call myfft3d%ifft3(fhat,fold)
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
                write(*,20) (f(i,j,4), j = 1,gp%ysz(2))
            end do 
        print*, "------------------------------------" 
        print*, "Recalcuted:"
        print*, "------------------------------------" 
            do i = 1,gp%ysz(1)
                write(*,20) (fold(i,j,4), j = 1,gp%ysz(2))
            end do 
        
    end if 
   
    mymaxerr = MAXVAL(ABS(3*d2fdx2 - fold))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    
    
    
    call myfft3d%destroy


    deallocate(x, y, z, f, fold, fhat, d2fdx2,d2fdy2)
    call decomp_2d_finalize
    call MPI_Finalize(ierr)
    
20 format(1x,16D10.3)
end program
