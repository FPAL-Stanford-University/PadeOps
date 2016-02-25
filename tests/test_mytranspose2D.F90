program mpi_transpose

    use mpi
    use kind_parameters, only: rkind, clen
    use mytranspose2DMod, only: mytranspose2D
    implicit none

    integer, parameter :: nx = 16, ny = 8
    integer, parameter :: ax = 4, ay = 2

    real(rkind), dimension(ax,ny) :: x
    real(rkind), dimension(nx,ay) :: x_trans

    type(mytranspose2D) :: gp2D

    integer :: rank, nprocs
    integer :: i, j, ierr

    character(len=clen) :: ioformat

    call MPI_Init(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,   rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

    if (nprocs /= 4) then
        if (rank == 0) print *, "Not the right number of processors. Re-run with 4 procs."
        stop
    end if

    ! Initialize gp2D
    call gp2D%init(nx,ny,MPI_COMM_WORLD)

    do j = 1,ny
        do i = 1,ax
            x(i,j) = i + rank*ax + (j-1)*nx
        end do
    end do

    write(ioformat,'(A,I3,A)') '(', ax, 'F6.1)'
    call sleep(rank)
    print *, "x: ", rank
    do j = 1,ny
        write(*,ioformat) (x(i,j), i=1,ax)
    end do

    !!!! Y to X transpose !!!!
    call gp2D%transpose_y_to_x(x,x_trans)
    !!!! ================ !!!!

    if (rank == 0) print *
    if (rank == 0) print *, "Finished Y to X transpose"

    write(ioformat,'(A,I3,A)') '(', nx, 'F6.1)'
    call sleep(rank)
    print *, "x_trans: ", rank
    do j = 1,ay
        write(*,ioformat) (x_trans(i,j), i=1,nx)
    end do

    !!!! Y to X transpose !!!!
    call gp2D%transpose_x_to_y(x_trans,x)
    !!!! ================ !!!!
   
    if (rank == 0) print *
    if (rank == 0) print *, "Finished X to Y transpose"

    write(ioformat,'(A,I3,A)') '(', ax, 'F6.1)'
    call sleep(rank)
    print *, "x: ", rank
    do j = 1,ny
        write(*,ioformat) (x(i,j), i=1,ax)
    end do

    call MPI_Finalize(ierr)

end program
