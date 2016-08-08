program test_derivatives_parallel
    use mpi
    use kind_parameters, only : rkind
    use decomp_2d
    use constants, only: pi, two
    use DerivativesMod, only: derivatives

    implicit none

    real(rkind), dimension(:,:,:), allocatable :: x, y, z, f,dfdx, dfdy, dfdz, df
    real(rkind), dimension(:,:,:), allocatable :: f_in_x, df_in_x, f_in_z, df_in_z 
    type(derivatives) :: method1, method2, method3
    type(decomp_info) :: gp

    integer :: nx = 128, ny = 128, nz = 128
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k
    real(rkind) :: maxerr, mymaxerr
    
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
    allocate( dfdx      ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( dfdy      ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( dfdz      ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( df       ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    
   
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
 
    ! Initialize everything 

    call method1%init(                         gp, &
                                dx,     dy,    dz, &
                            .TRUE., .TRUE., .TRUE., &
                            "cd10", "cd10", "cd10" )

    !call method1%init(                         gp, &
    !                            dx,     dy,    dz, &
    !                        .FALSE., .FALSE., .FALSE., &
    !                        "cd10", "cd10", "cd10" )

    call method2%init(                         gp, &
                                dx,     dy,    dz, &
                            .TRUE., .TRUE.,  .TRUE., &
                            "cd06", "cd06", "cd06" )
    
    call method3%init(                         gp, &
                                dx,     dy,    dz, &
                            .TRUE., .TRUE.,  .TRUE., &
                            "four", "four", "four" )

    call mpi_barrier(mpi_comm_world, ierr) 
    if (nrank == 0) print*, "Initialized all methods" 

  
    if (nrank == 0) then 
    print*, "==========================================="
    print*, "Now trying METHOD 1: CD10"
    print*, "==========================================="
    end if 

    allocate(df_in_x(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(f_in_x(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    t0 = MPI_WTIME()
    call transpose_y_to_x(f,f_in_x,gp)
    call method1 % ddx(f_in_x,df_in_x)
    call transpose_x_to_y(df_in_x,df,gp)
    t1 = MPI_WTIME()
    mymaxerr = MAXVAL(ABS(df - dfdx))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    if (nrank == 0) print*, "Time to solve in x:", t1 - t0
    deallocate (df_in_x, f_in_x)

    call mpi_barrier(mpi_comm_world, ierr) 
    t0 = MPI_WTIME()
    call method1 % ddy(f,df)
    t1 = MPI_WTIME()
    mymaxerr = MAXVAL(ABS(df - dfdy))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    if (nrank == 0) print*, "Time to solve in y:", t1 - t0

    allocate(df_in_z(gp%zsz(1),gp%zsz(2),gp%zsz(3))) 
    allocate(f_in_z(gp%zsz(1),gp%zsz(2),gp%zsz(3))) 
    call mpi_barrier(mpi_comm_world, ierr) 
    t0 = MPI_WTIME()
    call transpose_y_to_z(f,f_in_z,gp)
    call method1 % ddz(f_in_z,df_in_z)
    call transpose_z_to_y(df_in_z,df,gp)
    t1 = MPI_WTIME()
    mymaxerr = MAXVAL(ABS(df - dfdz))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    if (nrank == 0) print*, "Time to solve in z:", t1 - t0
    deallocate (df_in_z, f_in_z)


    if (nrank == 0) then
    print*, "==========================================="
    print*, "Now trying METHOD 2: CD06"
    print*, "==========================================="
    end if 

    allocate(df_in_x(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(f_in_x(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    call mpi_barrier(mpi_comm_world, ierr) 
    t0 = MPI_WTIME()
    call transpose_y_to_x(f,f_in_x,gp)
    call method2 % ddx(f_in_x,df_in_x)
    call transpose_x_to_y(df_in_x,df,gp)
    t1 = MPI_WTIME()
    mymaxerr = MAXVAL(ABS(df - dfdx))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    if (nrank == 0) print*, "Time to solve in x:", t1 - t0
    deallocate (df_in_x, f_in_x)

    call mpi_barrier(mpi_comm_world, ierr) 
    t0 = MPI_WTIME()
    call method2 % ddy(f,df)
    t1 = MPI_WTIME()
    mymaxerr = MAXVAL(ABS(df - dfdy))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    if (nrank == 0) print*, "Time to solve in y:", t1 - t0

    allocate(df_in_z(gp%zsz(1),gp%zsz(2),gp%zsz(3))) 
    allocate(f_in_z(gp%zsz(1),gp%zsz(2),gp%zsz(3))) 
    call mpi_barrier(mpi_comm_world, ierr) 
    t0 = MPI_WTIME()
    call transpose_y_to_z(f,f_in_z,gp)
    call method2 % ddz(f_in_z,df_in_z)
    call transpose_z_to_y(df_in_z,df,gp)
    t1 = MPI_WTIME()
    mymaxerr = MAXVAL(ABS(df - dfdz))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    if (nrank == 0) print*, "Time to solve in z:", t1 - t0
    deallocate(df_in_z, f_in_z)

    if (nrank == 0) then
    print*, "==========================================="
    print*, "Now trying METHOD 3: FOUR"
    print*, "==========================================="
    end if 

    allocate(df_in_x(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    allocate(f_in_x(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    call mpi_barrier(mpi_comm_world, ierr) 
    t0 = MPI_WTIME()
    call transpose_y_to_x(f,f_in_x,gp)
    call method3 % ddx(f_in_x,df_in_x)
    call transpose_x_to_y(df_in_x,df,gp)
    t1 = MPI_WTIME()
    mymaxerr = MAXVAL(ABS(df - dfdx))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    if (nrank == 0) print*, "Time to solve in x:", t1 - t0
    deallocate (df_in_x, f_in_x)

    call mpi_barrier(mpi_comm_world, ierr) 
    t0 = MPI_WTIME()
    call method3 % ddy(f,df)
    t1 = MPI_WTIME()
    mymaxerr = MAXVAL(ABS(df - dfdy))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    if (nrank == 0) print*, "Time to solve in y:", t1 - t0

    allocate(df_in_z(gp%zsz(1),gp%zsz(2),gp%zsz(3))) 
    allocate(f_in_z(gp%zsz(1),gp%zsz(2),gp%zsz(3))) 
    call mpi_barrier(mpi_comm_world, ierr) 
    t0 = MPI_WTIME()
    call transpose_y_to_z(f,f_in_z,gp)
    call method3 % ddz(f_in_z,df_in_z)
    call transpose_z_to_y(df_in_z,df,gp)
    t1 = MPI_WTIME()
    mymaxerr = MAXVAL(ABS(df - dfdz))
    call MPI_Reduce(mymaxerr, maxerr, 1, real_type, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Maximum error = ", maxerr
    if (nrank == 0) print*, "Time to solve in z:", t1 - t0
    deallocate(df_in_z, f_in_z)


    deallocate(x, y, z, f, dfdx,dfdy, dfdz)
    call decomp_2d_finalize
    call MPI_Finalize(ierr)

end program 
