program test_transpose

    use mpi

    use decomp_2d,       only: decomp_2d_init, decomp_2d_finalize, get_decomp_info, &
                               decomp_info,decomp_info_init,decomp_info_finalize,   &
                               transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y

    use kind_parameters, only: rkind
    use constants,       only: zero,two,pi
    use DerivativesMod,  only: derivatives
    use operators,       only: curl
    use reductions,      only: P_MAXVAL, P_AVGZ
    use exits,           only: message
    
    implicit none

    type(decomp_info) :: gp
    type(derivatives) :: der

    integer :: i, j, k, ierr

    integer :: nx = 128, ny = 128, nz = 128
    integer :: prow = 4, pcol = 4

    logical :: periodicx = .FALSE.
    logical :: periodicy = .FALSE.
    logical :: periodicz = .FALSE.

    real(rkind) :: dx, dy, dz
    real(rkind), dimension(:,:,:), allocatable :: x,y,z,u,v,w
    
    real(rkind), dimension(:,:,:,:), allocatable :: vort, vort_exact
    real(rkind), dimension(:,:), allocatable :: avg

    call MPI_Init(ierr)

    call decomp_2d_init(nx, ny, nz, prow, pcol, [periodicx, periodicy, periodicz])
    call get_decomp_info(gp)
    
    ! Allocate all physical data in Y pencils
    allocate( x         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( y         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( z         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( u         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( v         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( w         ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( vort      ( gp%ysz(1), gp%ysz(2), gp%ysz(3), 3 ) )
    allocate( vort_exact( gp%ysz(1), gp%ysz(2), gp%ysz(3), 3 ) )
    allocate( avg       ( gp%ysz(1), gp%ysz(2) ) )

    dx = two*pi/real(nx,rkind)
    dy = dx
    dz = dx

    ! Initialize derivatives 
    call der%init(                                    gp, &
                       dx,            dy,             dz, &
                periodicx,     periodicy,      periodicz, &
                   "cd10",        "cd10",         "cd10", &    ! Derivative type
                  .false.,       .false.,        .false., &    ! Metrics?
                  .false.)                                     ! Curvilinear?

    do k=1,gp%ysz(3)
        do j=1,gp%ysz(2)
            do i=1,gp%ysz(1)
                x(i,j,k) = real( gp%yst(1)-1+i-1,rkind )*dx
                y(i,j,k) = real( gp%yst(2)-1+j-1,rkind )*dy
                z(i,j,k) = real( gp%yst(3)-1+k-1,rkind )*dz
                
                u(i,j,k) =  y(i,j,k)
                v(i,j,k) = -x(i,j,k)
                w(i,j,k) = zero
                
                vort_exact(i,j,k,1) = zero
                vort_exact(i,j,k,2) = zero
                vort_exact(i,j,k,3) = -two
            end do
        end do
    end do

    ! Get the vorticity
    call curl( gp, der, u, v, w, vort)

    vort = vort - vort_exact
    call message("Z vorticity:")
    call message("    X vorticity error",P_MAXVAL(ABS(vort(:,:,:,1))))
    call message("    Y vorticity error",P_MAXVAL(ABS(vort(:,:,:,2))))
    call message("    Z vorticity error",P_MAXVAL(ABS(vort(:,:,:,3))))

    do k=1,gp%ysz(3)
        do j=1,gp%ysz(2)
            do i=1,gp%ysz(1)
                u(i,j,k) = zero
                v(i,j,k) =  z(i,j,k)
                w(i,j,k) = -y(i,j,k)
                
                vort_exact(i,j,k,1) = -two
                vort_exact(i,j,k,2) = zero
                vort_exact(i,j,k,3) = zero
            end do
        end do
    end do

    call P_AVGZ(gp,v,avg)

    call message("Z average error",P_MAXVAL(MAXVAL(ABS(avg - (pi - dx/two)))))

    ! Get the vorticity
    call curl( gp, der, u, v, w, vort)

    vort = vort - vort_exact
    call message("X vorticity:")
    call message("    X vorticity error",P_MAXVAL(ABS(vort(:,:,:,1))))
    call message("    Y vorticity error",P_MAXVAL(ABS(vort(:,:,:,2))))
    call message("    Z vorticity error",P_MAXVAL(ABS(vort(:,:,:,3))))

    do k=1,gp%ysz(3)
        do j=1,gp%ysz(2)
            do i=1,gp%ysz(1)
                u(i,j,k) = -z(i,j,k)
                v(i,j,k) = zero
                w(i,j,k) =  x(i,j,k)
                
                vort_exact(i,j,k,1) = zero
                vort_exact(i,j,k,2) = -two
                vort_exact(i,j,k,3) = zero
            end do
        end do
    end do

    ! Get the vorticity
    call curl( gp, der, u, v, w, vort)

    vort = vort - vort_exact
    call message("Y vorticity:")
    call message("    X vorticity error",P_MAXVAL(ABS(vort(:,:,:,1))))
    call message("    Y vorticity error",P_MAXVAL(ABS(vort(:,:,:,2))))
    call message("    Z vorticity error",P_MAXVAL(ABS(vort(:,:,:,3))))

    call der%destroy()
    
    deallocate( x          )
    deallocate( y          )
    deallocate( z          )
    deallocate( u          )
    deallocate( v          )
    deallocate( w          )
    deallocate( vort       )
    deallocate( vort_exact )

    call decomp_2d_finalize

    call MPI_Finalize(ierr)

end program
