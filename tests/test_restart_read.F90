program test_restart_read
    use decomp_2d
    use decomp_2d_io
    use kind_parameters, only: rkind
    use constants, only: one, pi, two, imi

    implicit none 

    integer :: ierr
    real(rkind), dimension(:,:,:), allocatable :: zE
    real(rkind), dimension(:,:,:), allocatable :: fE, dfE
    real(rkind), dimension(:,:,:), allocatable :: fEt, dfEt
    integer :: nx = 500, ny = 500, nz = 60
    real(rkind) :: omega = 2._rkind, dz
    integer :: i, j, k
    type(decomp_info) :: gp

    call mpi_init(ierr)
    call decomp_2d_init(nx,ny,nz,1,4)
    call get_decomp_info(gp)

    dz = one/real(nz,rkind)
    allocate(zE(gp%xsz(1), gp%xsz(2), gp%xsz(3)))
    allocate(fE(gp%xsz(1), gp%xsz(2), gp%xsz(3)))
    allocate(dfE(gp%xsz(1), gp%xsz(2), gp%xsz(3)))
    allocate(fEt(gp%xsz(1), gp%xsz(2), gp%xsz(3)))
    allocate(dfEt(gp%xsz(1), gp%xsz(2), gp%xsz(3)))

    do k = 1,gp%xsz(3)
        do j = 1,gp%xsz(2)
            do i = 1,gp%xsz(1)
                zE(i,j,k) = real((gp%xst(3) + k - 1),rkind)*dz
            end do 
        end do 
    end do 

    fEt = sin(two*pi*omega*zE) 
    dfEt = two*pi*omega*cos(two*pi*omega*zE) 

    call decomp_2d_read_one(1,fE,'fE.dat')
    call decomp_2d_read_one(1,dfE,'dfE.dat')

    call sleep(nrank)
    print*, nrank, maxval(abs(fE-fEt))

    call decomp_2d_finalize
    call mpi_finalize(ierr)
    deallocate(zE, fE, dfE)


end program
