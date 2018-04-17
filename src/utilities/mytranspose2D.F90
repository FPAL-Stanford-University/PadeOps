module mytranspose2DMod

    use kind_parameters, only: rkind, mpirkind
    use mpi
    implicit none

    private

    public mytranspose2D

    type mytranspose2D

        integer :: comm
        integer :: nx, ny
        integer :: ax, ay
        integer :: nprocs

    contains

        procedure :: init

        procedure :: transpose_x_to_y
        procedure :: transpose_y_to_x

    end type

contains

    subroutine init(this, nx, ny, comm)
        use exits, only: GracefulExit
        class(mytranspose2D), intent(inout) :: this
        integer, intent(in) :: nx, ny
        integer, intent(in) :: comm

        integer :: ierr

        this%nx = nx
        this%ny = ny

        this%comm = comm
        call MPI_COMM_SIZE(this%comm, this%nprocs, ierr)

        this%ax = this%nx / this%nprocs
        this%ay = this%ny / this%nprocs

        if (this%ax*this%nprocs /= this%nx) then
            call GracefulExit("nx not a multiple of nprocs in mytranspose2D",231)
        end if
        if (this%ay*this%nprocs /= this%ny) then
            call GracefulExit("ny not a multiple of nprocs in mytranspose2D",231)
        end if

    end subroutine

    subroutine transpose_x_to_y(this, input, output)
        class(mytranspose2D), intent(in) :: this
        real(rkind), dimension(this%nx,this%ay), intent(in)  :: input
        real(rkind), dimension(this%ax,this%ny), intent(out) :: output
        
        real(rkind), dimension(this%ax,this%ny) :: trans
        integer :: i, j, iblock, ierr

        ! Rearrange local data
        do j = 1,this%ay
            do iblock = 1,this%nprocs
                do i = 1,this%ax
                    trans(i,j+(iblock-1)*this%ay) = input(i+(iblock-1)*this%ax,j)
                end do
            end do
        end do

        ! Do the All-to-All collective
        call MPI_Alltoall(trans, this%ax*this%ay, mpirkind, output, this%ax*this%ay, mpirkind, this%comm, ierr)

    end subroutine

    subroutine transpose_y_to_x(this, input, output)
        class(mytranspose2D), intent(in) :: this
        real(rkind), dimension(this%ax,this%ny), intent(in)  :: input
        real(rkind), dimension(this%nx,this%ay), intent(out) :: output
        
        real(rkind), dimension(this%ax,this%ny) :: trans
        integer :: i, j, iblock, ierr

        ! Do the All-to-All collective
        call MPI_Alltoall(input, this%ax*this%ay, mpirkind, trans, this%ax*this%ay, mpirkind, this%comm, ierr)

        ! Rearrange local data
        do j = 1,this%ay
            do iblock = 1,this%nprocs
                do i = 1,this%ax
                    output(i+(iblock-1)*this%ax,j) = trans(i,j+(iblock-1)*this%ay)
                end do
            end do
        end do

    end subroutine

end module
