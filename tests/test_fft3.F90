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

    integer :: nx = 128, ny = 128, nz =128
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k
    real(rkind) :: maxerr, mymaxerr
    
    logical :: get_exhaustive_plan = .false.
    logical, parameter :: verbose = .false. 
    real(rkind) :: dx, dy, dz
    character(len=1), parameter :: base_dec = "x"

    double precision :: t0, t1

    call MPI_Init(ierr)
    
    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)
   
    select case (base_dec)
    case ("y")
        ! Initialize input data (y decomposition) 
        allocate( x         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
        allocate( y         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
        allocate( z         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
        allocate( f         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
        allocate( d2fdx2    ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
        allocate( d2fdy2    ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
        allocate( d2fdz2    ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
        allocate( fold      ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    case ("x")
        ! Initialize input data (y decomposition) 
        allocate( x         ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
        allocate( y         ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
        allocate( z         ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
        allocate( f         ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
        allocate( d2fdx2    ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
        allocate( d2fdy2    ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
        allocate( d2fdz2    ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
        allocate( fold      ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    end select 

    dx = two*pi/nx
    dy = two*pi/ny
    dz = two*pi/nz
    ierr = myfft3d%init(nx,ny,nz,base_dec,dx, dy,dz, get_exhaustive_plan)
    call myfft3d%alloc_output(fhat) 
    if (nrank == 0) print*, "Initialization complete"

    select case (base_dec)
    case ("y")
    ! Initialize fft_2d type  
        do k = 1,gp%ysz(3)
            do j = 1,gp%ysz(2)
                do i = 1,gp%ysz(1)
                    x(i,j,k) = real(gp%yst(1) - 1 + i - 1, rkind)*dx
                    y(i,j,k) = real(gp%yst(2) - 1 + j - 1, rkind)*dy 
                    z(i,j,k) = real(gp%yst(3) - 1 + k - 1, rkind)*dz
                end do 
            end do 
        end do 
    case ("x")
    ! Initialize fft_2d type  
        do k = 1,gp%xsz(3)
            do j = 1,gp%xsz(2)
                do i = 1,gp%xsz(1)
                    x(i,j,k) = real(gp%xst(1) - 1 + i - 1, rkind)*dx
                    y(i,j,k) = real(gp%xst(2) - 1 + j - 1, rkind)*dy 
                    z(i,j,k) = real(gp%xst(3) - 1 + k - 1, rkind)*dz
                end do 
            end do 
        end do 
    end select 

    f = sin(x)*sin(y)*sin(z)
    !f = x*y*z
    d2fdx2 = -sin(x)*sin(y)*sin(z)
    d2fdy2 = -sin(x)*sin(y)*sin(z)
    d2fdz2 = -sin(x)*sin(y)*sin(z)
   
    call MPI_BARRIER(mpi_comm_world,ierr)
    t0 = MPI_WTIME()
    ! Forward transform
    select case (base_dec)
    case ("y") 
        call myfft3d%fft3_y2y(f,fhat)
    case ("x")
        call myfft3d%fft3_x2z(f,fhat)
    end select
    t1 = MPI_WTIME()
    
    if (nrank == 0) print*, "Forward transform time:", t1 - t0

    ! Compute the laplacian as a check 
    fhat = -fhat * myfft3d%kabs_sq
   
    t0 = MPI_WTIME()
    ! Backward transform 
    select case (base_dec)
    case ("y")
        call myfft3d%ifft3_y2y(fhat,fold)
    case ("x")
        call myfft3d%ifft3_z2x(fhat,fold)
    end select 
    t1 = MPI_WTIME()

    if (nrank == 0) print*, "Backward transform time:", t1 - t0
   
    if (verbose) then
        call sleep(nrank)
        print*, "-------------------------------------" 
        print*, "My rank:", nrank
        print*, "x edges:", gp%yst(1), gp%yen(1)
        print*, "y edges:", gp%yst(2), gp%yen(2)
        print*, "z edges:", gp%yst(3), gp%yen(3)
        print*, "-------------------------------------" 
            do i = 1,size(f,1)
                write(*,20) (f(i,j,3), j = 1,size(f,2))
            end do 
        print*, "------------------------------------" 
        print*, "Recalcuted:"
        print*, "------------------------------------" 
            do i = 1,size(fold,1)
                write(*,20) (fold(i,j,3), j = 1,size(fold,2))
            end do 
        
    end if 
   
    mymaxerr = MAXVAL(ABS(3._rkind*d2fdx2 - fold))
    !mymaxerr = MAXVAL(ABS(f - fold))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    
    
    
    call myfft3d%destroy


    deallocate(x, y, z, f, fold, fhat, d2fdx2,d2fdy2)
    call decomp_2d_finalize
    call MPI_Finalize(ierr)
    
20 format(1x,16D10.3)
end program
