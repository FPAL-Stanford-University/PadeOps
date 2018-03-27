module StatisticsMod

    use mpi
    use decomp_2d,       only: decomp_info, nrank
    use kind_parameters, only: rkind, mpirkind
    use constants,       only: half,one
    use exits,           only: GracefulExit

    implicit none

    type :: statistics
        private

        type(decomp_info), pointer :: gp
        logical, dimension(3) :: averaging_directions = [.true., .true., .true.]
        integer, dimension(3) :: buffer_size_x, buffer_size_y, buffer_size_z

        integer :: xy_comm, yz_comm, xz_comm
        integer :: x_comm, y_comm, z_comm

        real(rkind), dimension(:,:,:), allocatable :: buffer_x, buffer_y, buffer_z

        integer, dimension(3), public :: sz, avg_size

    contains

        procedure, private :: avg_x
        procedure, private :: avg_y
        procedure, private :: avg_z
        procedure          :: allocate_average
        procedure          :: get_average
        procedure          :: get_weighted_average
        procedure          :: get_fluctuations
        final              :: destroy

    end type

    interface statistics
        module procedure init
    end interface

contains

    function init(gp, pencil, averaging_directions) result(this)
        type(statistics) :: this
        class(decomp_info), target, intent(in) :: gp
        integer,                    intent(in) :: pencil
        logical, dimension(3),      intent(in) :: averaging_directions

        integer, dimension(3) :: st
        integer :: ierr

        this%gp => gp
        this%averaging_directions = averaging_directions

        if(pencil == 1) then        ! X-pencil
            st = gp%xst
            this%sz = gp%xsz
        else if (pencil == 2) then  ! Y-pencil
            st = gp%yst
            this%sz = gp%ysz
        else if (pencil == 3) then  ! Z-pencil
            st = gp%zst
            this%sz = gp%zsz
        else
            call GracefulExit("Pencil has to be one of 1 (x), 2 (y) or 3(z)",4587)
        end if

        call mpi_comm_split(mpi_comm_world, st(1), nrank, this%yz_comm, ierr)
        call mpi_comm_split(mpi_comm_world, st(2), nrank, this%xz_comm, ierr)
        call mpi_comm_split(mpi_comm_world, st(3), nrank, this%xy_comm, ierr)

        call mpi_comm_split( this%yz_comm,  st(2), nrank,  this%z_comm, ierr)
        call mpi_comm_split( this%yz_comm,  st(3), nrank,  this%y_comm, ierr)
        call mpi_comm_split( this%xz_comm,  st(3), nrank,  this%x_comm, ierr)

        this%buffer_size_x = this%sz
        if ( this%averaging_directions(1) ) then
            this%buffer_size_x(1) = 1

            if (allocated(this%buffer_x)) deallocate(this%buffer_x)
            allocate(this%buffer_x(this%buffer_size_x(1), this%buffer_size_x(2), this%buffer_size_x(3) ))
        end if

        this%buffer_size_y = this%buffer_size_x
        if ( this%averaging_directions(2) ) then
            this%buffer_size_y(2) = 1

            if (allocated(this%buffer_y)) deallocate(this%buffer_y)
            allocate(this%buffer_y(this%buffer_size_y(1), this%buffer_size_y(2), this%buffer_size_y(3) ))
        end if

        this%buffer_size_z = this%buffer_size_y
        if ( this%averaging_directions(3) ) then
            this%buffer_size_z(3) = 1

            if (allocated(this%buffer_z)) deallocate(this%buffer_z)
            allocate(this%buffer_z(this%buffer_size_z(1), this%buffer_size_z(2), this%buffer_size_z(3) ))
        end if

        this%avg_size = this%buffer_size_z

    end function

    impure elemental subroutine destroy(this)
        type(statistics), intent(inout) :: this

        if (allocated(this%buffer_x)) deallocate(this%buffer_x)
        if (allocated(this%buffer_y)) deallocate(this%buffer_y)
        if (allocated(this%buffer_z)) deallocate(this%buffer_z)
    end subroutine

    subroutine allocate_average(this, f_avg)
        class(statistics),                          intent(in)  :: this
        real(rkind), dimension(:,:,:), allocatable, intent(out) :: f_avg

        if (allocated(f_avg)) deallocate(f_avg)
        allocate(f_avg(this%avg_size(1), this%avg_size(2), this%avg_size(3)))
    end subroutine

    subroutine avg_x(this, f, avg)
        class(statistics),                               intent(in)    :: this
        real(rkind), dimension(:,:,:),                   intent(in)    :: f
        real(rkind), dimension(1, size(f,2), size(f,3)), intent(out)   :: avg
        real(rkind), dimension(size(f,2), size(f,3))                   :: p_avg
        integer :: ierr
        
        if (size(f,1) == 1) then
            avg = f
            return
        end if

        if (size(f,1) == this%gp%xsz(1)) then
            avg(1,:,:) = SUM(f,1) / real(this%gp%xsz(1),rkind)
            return
        end if

        p_avg = SUM(f,1)
        call MPI_ALLREDUCE(p_avg, avg(1,:,:), size(f,2)*size(f,3), mpirkind, MPI_SUM, this%x_comm, ierr)
        avg = avg / real(this%gp%xsz(1),rkind)

    end subroutine

    subroutine avg_y(this, f, avg)
        class(statistics),                               intent(in)    :: this
        real(rkind), dimension(:,:,:),                   intent(in)    :: f
        real(rkind), dimension(size(f,1), 1, size(f,3)), intent(out)   :: avg
        real(rkind), dimension(size(f,1), size(f,3))                   :: p_avg
        integer :: ierr
        
        if (size(f,2) == 1) then
            avg = f
            return
        end if

        if (size(f,2) == this%gp%ysz(2)) then
            avg(:,1,:) = SUM(f,2) / real(this%gp%ysz(2),rkind)
            return
        end if

        p_avg = SUM(f,2)
        call MPI_ALLREDUCE(p_avg, avg(:,1,:), size(f,1)*size(f,3), mpirkind, MPI_SUM, this%y_comm, ierr)
        avg = avg / real(this%gp%ysz(2),rkind)

    end subroutine

    subroutine avg_z(this, f, avg)
        class(statistics),                               intent(in)    :: this
        real(rkind), dimension(:,:,:),                   intent(in)    :: f
        real(rkind), dimension(size(f,1), size(f,2), 1), intent(out)   :: avg
        real(rkind), dimension(size(f,1), size(f,2))                   :: p_avg
        integer :: ierr
        
        if (size(f,3) == 1) then
            avg = f
            return
        end if

        if (size(f,3) == this%gp%zsz(3)) then
            avg(:,:,1) = SUM(f,3) / real(this%gp%zsz(3),rkind)
            return
        end if

        p_avg = SUM(f,3)
        call MPI_ALLREDUCE(p_avg, avg(:,:,1), size(f,1)*size(f,2), mpirkind, MPI_SUM, this%z_comm, ierr)
        avg = avg / real(this%gp%zsz(3),rkind)

    end subroutine

    subroutine get_average(this,f,f_avg)
        class(statistics),                                                          target, intent(inout) :: this
        real(rkind), dimension(      this%sz(1),      this%sz(2),      this%sz(3)), target, intent(in)    :: f
        real(rkind), dimension(this%avg_size(1),this%avg_size(2),this%avg_size(3)),         intent(out)   :: f_avg

        real(rkind), dimension(:,:,:), pointer :: tmp

        if ( this%averaging_directions(1) .and. (size(f_avg,1) /= 1) ) then
            call GracefulExit("f_avg needs to have size 1 in the averaging directions", 764)
        end if

        if ( this%averaging_directions(2) .and. (size(f_avg,2) /= 1) ) then
            call GracefulExit("f_avg needs to have size 1 in the averaging directions", 764)
        end if

        if ( this%averaging_directions(3) .and. (size(f_avg,3) /= 1) ) then
            call GracefulExit("f_avg needs to have size 1 in the averaging directions", 764)
        end if

        ! Set tmp to be input array
        tmp => f

        if (this%averaging_directions(1)) then
            call this%avg_x(tmp, this%buffer_x)
            tmp => this%buffer_x
        end if

        if (this%averaging_directions(2)) then
            call this%avg_y(tmp, this%buffer_y)
            tmp => this%buffer_y
        end if

        if (this%averaging_directions(3)) then
            call this%avg_z(tmp, this%buffer_z)
            tmp => this%buffer_z
        end if

        f_avg = tmp

    end subroutine

    subroutine get_weighted_average(this,weights,f,f_avg)
        class(statistics),                                                          intent(inout) :: this
        real(rkind), dimension(      this%sz(1),      this%sz(2),      this%sz(3)), intent(in)    :: f, weights
        real(rkind), dimension(this%avg_size(1),this%avg_size(2),this%avg_size(3)), intent(out)   :: f_avg

        real(rkind), dimension(this%avg_size(1),this%avg_size(2),this%avg_size(3)) :: weights_avg

        call this%get_average(weights*f, f_avg)
        call this%get_average(weights, weights_avg)

        f_avg = f_avg / weights_avg
    end subroutine

    subroutine get_fluctuations(this, f, f_avg, f_fluct)
        class(statistics),                                                          intent(in)  :: this
        real(rkind), dimension(      this%sz(1),      this%sz(2),      this%sz(3)), intent(in)  :: f
        real(rkind), dimension(this%avg_size(1),this%avg_size(2),this%avg_size(3)), intent(in)  :: f_avg
        real(rkind), dimension(      this%sz(1),      this%sz(2),      this%sz(3)), intent(out) :: f_fluct

        integer :: i, j, k, ii, jj, kk

        do k = 1,this%sz(3)
            do j = 1,this%sz(2)
                do i = 1,this%sz(2)
                    ii = min(i, this%avg_size(1))
                    jj = min(j, this%avg_size(2))
                    kk = min(k, this%avg_size(3))
                    f_fluct(i,j,k) = f(i,j,k) - f_avg(ii,jj,kk)
                end do
            end do
        end do
    end subroutine

end module
