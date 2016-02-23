program test_io_VTK
    use mpi
    use kind_parameters, only : rkind
    use decomp_2d, only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, nrank
    use io_VTK_stuff, only: io_VTK
    use constants, only: pi, two

    implicit none

    type(io_VTK) :: viz

    real(rkind), dimension(:,:,:,:), allocatable, target :: mesh
    real(rkind), dimension(:,:,:), pointer :: x, y, z

    real(rkind), dimension(:,:,:,:), allocatable, target :: primary, secondary
    real(rkind), dimension(:,:,:), pointer :: f, g
    type(decomp_info) :: gp

    logical :: periodicx = .TRUE., periodicy = .TRUE., periodicz = .TRUE.

    integer :: nx = 128, ny = 128, nz = 128
    integer :: nxp, nyp, nzp
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k, step

    real(rkind) :: dx, dy, dz

    call MPI_Init(ierr)
    
    call decomp_2d_init(nx, ny, nz, prow, pcol, [periodicx, periodicy, periodicz])
    call get_decomp_info(gp)
  
    nxp = gp%ysz(1); nyp = gp%ysz(2); nzp = gp%ysz(3)

    ! Initialize input data (y decomposition) 
    allocate( mesh ( gp%ysz(1), gp%ysz(2), gp%ysz(3), 3 ) )
    allocate( primary ( gp%ysz(1), gp%ysz(2), gp%ysz(3), 1 ) )
    allocate( secondary ( gp%ysz(1), gp%ysz(2), gp%ysz(3), 1 ) )
    
    x => mesh(:,:,:,1)
    y => mesh(:,:,:,2)
    z => mesh(:,:,:,3)

    f => primary(:,:,:,1)
    g => secondary(:,:,:,1)

    dx = two*pi/nx
    dy = two*pi/ny
    dz = two*pi/nz

    ! Generate input data
    do k = 1,nzp
        do j = 1,nyp
            do i = 1,nxp
                x(i,j,k) = real(gp%yst(1) - 1 + i - 1, rkind)*dx
                y(i,j,k) = real(gp%yst(2) - 1 + j - 1, rkind)*dy 
                z(i,j,k) = real(gp%yst(3) - 1 + k - 1, rkind)*dz
            end do 
        end do 
    end do
    
    ! Initialize io_VTK object
    call viz%init('.', 'VTK', 1, ['field'])

    ! Create the fields
    f = sin(x)*sin(y)*cos(z)
    g = real(nrank,rkind)
    
    if (nrank == 0) print*, "Created coordinates and fields"

    do step = 1,4
        call viz%WriteViz(gp, mesh, primary, real(0.,rkind), secondary, ['rank'])
    end do
  
    deallocate(mesh, primary, secondary)
    call viz%destroy()
    call decomp_2d_finalize
    call MPI_Finalize(ierr)

end program 
