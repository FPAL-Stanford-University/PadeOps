program dealasingTest3by2rule
    use mpi
    use kind_parameters, only: rkind
    use timer,           only: tic, toc
    use decomp_2d
    use fft_3d_stuff, only: fft_3d 
    use constants, only: pi, two
    use reductions, only: p_maxval
    use exits, only: message

    real(rkind), dimension(:,:,:), allocatable :: x,y,z, f, fsq, fup
    type(fft_3d) :: myfft3d
    type(decomp_info) :: gp

    integer :: nx =128, ny = 128, nz =128
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k
    real(rkind) :: maxerr, mymaxerr
    real(rkind) :: dx, dy, dz
    logical :: get_exhaustive_plan = .true. 

    call MPI_Init(ierr)
    
    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)
   
    allocate( x         ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( y         ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( z         ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( f         ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( fsq       ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    
    dx = two*pi/nx
    dy = two*pi/ny
    dz = two*pi/nz
    ierr = myfft3d%init(nx,ny,nz,"x",dx, dy,dz, get_exhaustive_plan)

     do k = 1,gp%xsz(3)
         do j = 1,gp%xsz(2)
             do i = 1,gp%xsz(1)
                 x(i,j,k) = real(gp%xst(1) - 1 + i - 1, rkind)*dx
                 y(i,j,k) = real(gp%xst(2) - 1 + j - 1, rkind)*dy 
                 z(i,j,k) = real(gp%xst(3) - 1 + k - 1, rkind)*dz
             end do 
         end do 
     end do


    f = sin(8*y)

    call myfft3d%alloc_upsampledArr(fup)
   
    do k = 1,2 
    call myfft3d%upsample(f,fup) 
    fup = fup*fup
    call myfft3d%downsample(fup,fsq)
    fsq = abs(fsq - f*f)
    mymaxerr = maxval(fsq)
    maxerr = p_maxval(mymaxerr)
    call message(0,"Maximum error:",mymaxerr)
    end do 

    call mpi_barrier(mpi_comm_world, ierr)
    call decomp_2d_finalize
    call MPI_Finalize(ierr)


end program 
