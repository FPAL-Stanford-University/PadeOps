program test_vectorExponential
    use kind_parameters, only: rkind
    use timer, only: tic, toc
    use mpi
    use decomp_2d

    implicit none
    integer :: idx, ix1, iy1, iz1, ixn, iyn, izn, i, j, k, ierr, prow = 0, pcol = 0 
    real(rkind), dimension(:,:,:), allocatable :: xG, yG, zG
    real(rkind), dimension(:,:,:), allocatable :: u, rhs
    type(decomp_info) :: gp 
    integer, parameter :: nx = 64, ny = 64, nz = 64
    integer :: ntot

    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(gp)

    allocate(u(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    allocate(rhs(gp%xsz(1),gp%xsz(2),gp%xsz(3))) 
    u = 12.d0
    call tic()
    rhs = exp(u)
    call toc()
    ntot = size(u)
    call tic()
    call vsexp (ntot, u, rhs)
    call toc()
    deallocate(u,rhs)

end program 
