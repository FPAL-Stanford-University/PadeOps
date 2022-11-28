program test_halo
    use mpi
    use decomp_2d
    use kind_parameters 
    use constants

    !integer, parameter :: nx = 128, ny = 65, nz=97
    integer, parameter :: nx = 8, ny = 6, nz=8
    real(rkind), dimension(:,:,:), allocatable :: f, fpad, fpadTrue
    type(decomp_info) :: gp
    integer :: ierr, i, j, k
    real(rkind), parameter :: L = two*pi
    real(rkind), parameter :: dx = L/nx, dy = L/ny, dz = L/nz
    real(rkind) :: x, y, z


    call MPI_Init(ierr)

    call decomp_2d_init(nx, ny, nz, 0, 0, [.TRUE., .TRUE., .TRUE.])
    call get_decomp_info(gp)
    


    allocate(f(gp%xsz(1), gp%xsz(2), gp%xsz(3)))
    allocate(fpadTrue(0:gp%xsz(1)+1, 0:gp%xsz(2)+1, 0:gp%xsz(3)+1))


    do k = 1,size(f,3)
        z = (k + gp%xst(3) - 1)*dz
        do j = 1,size(f,2)
            y = (j + gp%xst(2) - 1)*dy 
            do i = 1,size(f,1)
                x = (i + gp%xst(1) - 1)*dx 
                f(i,j,k) = cos(x)*sin(y)*cos(z)
            end do 
        end do 
    end do 

    do k = 0,size(f,3)+1
        z = (k + gp%xst(3) - 1)*dz
        do j = 0,size(f,2)+1
            y = (j + gp%xst(2) - 1)*dy 
            do i = 0,size(f,1)+1
                x = (i + gp%xst(1) - 1)*dx 
                fpadTrue(i,j,k) = cos(x)*sin(y)*cos(z)
            end do 
        end do 
    end do 

    call update_halo(f, fpad, 1, gp)


    print*, shape(f), shape(fpad)

    print*, f(5,2,:)
    print*, fpadTrue(5,2,:)
    print*, fpad(5,2,:)
    call decomp_2d_finalize

    call MPI_Finalize(ierr)

end program 
