module mytranspose
    use kind_parameters, only: rkind
    implicit none

    integer, parameter :: nx = 16, ny = 8, nz = 1
    integer, parameter :: ax = 4, ay = 8, az = 1

contains

    subroutine print_array(x)
        real(rkind), dimension(:,:,:) :: x
        integer :: i,j,k

        do k = 1,size(x,3)
            do j = 1,size(x,2)
                do i = 1,size(x,1)
                    print *, x(i,j,k)
                end do
            end do
        end do
    end subroutine

end module

program mpi_transpose

    use mpi
    use kind_parameters, only: rkind
    use mytranspose
    implicit none

    integer, dimension(ax,ay,az) :: x
    integer, dimension(nx,ny/4,az) :: x_trans
    integer, dimension(ax,ay) :: x_trans_local

    integer :: i,j,iblock
    integer :: rank, nprocs
    integer :: ierr

    call MPI_Init(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,   rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

    if (nprocs /= 4) then
        if (rank == 0) print *, "Not the right number of processors. Re-run with 4 procs."
        stop
    end if

    do j = 1,ay
        do i = 1,ax
            x(i,j,1) = i + rank*ax + (j-1)*nx
        end do
    end do

    call sleep(rank)
    print *, "x: ", rank
    do j = 1,ay
        write(*,*) (x(I,j,1), I=1,ax)
    end do

    !!!! Y to X transpose !!!!
    call MPI_Alltoall(x, ax*ny/nprocs, MPI_INTEGER, x_trans_local, ax*ny/nprocs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    
    do j = 1,ny/4
        do iblock = 1,4
            do i = 1,ax
                x_trans(i+(iblock-1)*ax,j,1) = x_trans_local(i,j+(iblock-1)*ny/nprocs)
            end do
        end do
    end do
    !!!! ================ !!!!

    if (rank == 0) print *, "Finished Y to X transpose"

    call sleep(rank)
    print *, "x_trans: ", rank
    do j = 1,ny/nprocs
        write(*,*) (x_trans(I,j,1), I=1,nx)
    end do

    !!!! Y to X transpose !!!!
    do j = 1,ny/4
        do iblock = 1,4
            do i = 1,ax
                x_trans_local(i,j+(iblock-1)*ny/nprocs) = x_trans(i+(iblock-1)*ax,j,1)
            end do
        end do
    end do
    
    call MPI_Alltoall(x_trans_local, ax*ny/nprocs, MPI_INTEGER, x, ax*ny/nprocs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    !!!! ================ !!!!
   
    if (rank == 0) print *, "Finished X to Y transpose"

    call sleep(rank)
    print *, "x: ", rank
    do j = 1,ay
        write(*,*) (x(I,j,1), I=1,ax)
    end do

    call MPI_Finalize(ierr)

end program
