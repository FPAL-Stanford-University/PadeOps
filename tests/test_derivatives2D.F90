program derivatives2D

    use mpi
    use kind_parameters,  only: rkind
    use constants,        only: two, pi
    use decomp_2d,        only: decomp_info, decomp_2d_init, get_decomp_info, decomp_2d_finalize, nrank
    use mytranspose2DMod, only: mytranspose2D
    use DerivativesMod,   only: derivatives
    use reductions,       only: P_MAXVAL
    use exits,            only: message
    implicit none

    integer, parameter :: nx = 64, ny = 32, nz = 8
    integer, parameter :: prow = 4, pcol = 2

    real(rkind), dimension(:,:,:), allocatable :: x,y

    real(rkind), dimension(:,:), allocatable :: f2D, der, der_exact
    real(rkind), dimension(:,:,:), allocatable :: f2D_trans, der_trans

    real(rkind) :: dx, dy

    type(decomp_info)   :: gp
    type(mytranspose2D) :: gp2D

    type(derivatives)   :: der2D

    integer :: xyproc, xyrank
    integer :: i,j, ierr

    integer :: XY_COMM

    call MPI_Init(ierr)

    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)

    dx = two*pi/real(nx,rkind)
    dy = two*pi/real(ny,rkind)
    call der2D%init(                    gp, &
                         dx,     dy,    dx, &
                     .TRUE., .TRUE.,.TRUE., &
                     "cd10", "cd10", "cd10" )

    call der2D%set_xsz( [gp%xsz(1), gp%xsz(2), 1] )
    call der2D%set_ysz( [gp%ysz(1), gp%ysz(2), 1] )

    call MPI_COMM_SPLIT(MPI_COMM_WORLD, gp%yst(3), nrank, XY_COMM, ierr)
    call MPI_COMM_RANK(XY_COMM, xyrank, ierr)
    call MPI_COMM_SIZE(XY_COMM, xyproc, ierr)
    
    ! Initialize input data (y decomposition) 
    allocate( x         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( y         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( f2D       ( gp%ysz(1), gp%ysz(2)            ) )
    allocate( der       ( gp%ysz(1), gp%ysz(2)            ) )
    allocate( der_exact ( gp%ysz(1), gp%ysz(2)            ) )
    allocate( f2D_trans ( gp%xsz(1), gp%xsz(2), 1         ) )
    allocate( der_trans ( gp%xsz(1), gp%xsz(2), 1         ) )

    ! Initialize gp2D
    call gp2D%init(nx,ny,XY_COMM)
    
    do j = 1,gp%ysz(2)
        do i = 1,gp%ysz(1)
            x(i,j,:) = real((gp%yst(1) - 1 + i),rkind)*dx
            y(i,j,:) = real((gp%yst(2) - 1 + j),rkind)*dy
        end do
    end do

    f2D = sin( x(:,:,1) )*cos( y(:,:,1) ) * real(gp%yst(3),rkind)
    der_exact = cos( x(:,:,1) )*cos( y(:,:,1) ) * real(gp%yst(3),rkind)
   
    !!!! Get 2D derivative !!!!
    call gp2D%transpose_y_to_x(f2D,f2D_trans)
    call der2D%ddx(f2D_trans, der_trans)
    call gp2D%transpose_x_to_y(der_trans,der)
    !!!! ================ !!!!
   
    call message("Maximum error in X derivative", P_MAXVAL(abs(der - der_exact)))

    der_exact = -sin( x(:,:,1) )*sin( y(:,:,1) ) * real(gp%yst(3),rkind)
    call der2D%ddy(f2D,der)
    call message("Maximum error in Y derivative", P_MAXVAL(abs(der - der_exact)))

    deallocate( x         )
    deallocate( y         )
    deallocate( f2D       )
    deallocate( f2D_trans )
    deallocate( der       )
    deallocate( der_exact )
    deallocate( der_trans )

    call der2D%destroy()
    call decomp_2d_finalize()
    call MPI_Finalize(ierr)

end program
