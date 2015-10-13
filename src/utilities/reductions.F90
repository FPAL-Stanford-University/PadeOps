module reductions
    use mpi
    use kind_parameters, only: rkind,mpirkind
    use decomp_2d, only: decomp_info, nrank, nproc
    use exits, only: GracefulExit
    
    implicit none

    private
    public :: P_MAXVAL, P_MINVAL, P_SUM, P_MEAN, P_AVGZ

    interface P_MAXVAL
        module procedure P_MAXVAL_arr3, P_MAXVAL_sca
    end interface

    interface P_MINVAL
        module procedure P_MINVAL_double
    end interface

    interface P_SUM
        module procedure P_SUM_double
    end interface

    interface P_MEAN
        module procedure P_MEAN_arr3, P_MEAN_sca
    end interface
contains
    
    function P_MEAN_sca(x) result(mean)
        real(rkind), intent(in) :: x
        real(rkind) :: mean
        real(rkind) :: summation
        integer :: ierr

        call MPI_Allreduce(x, summation, 1, mpirkind, MPI_MAX, MPI_COMM_WORLD, ierr)
        mean = summation/real(nproc,rkind)

    end function
    
    function P_MEAN_arr3(x) result(mean)
        real(rkind), dimension(:,:,:), intent(in) :: x
        real(rkind) :: mean
        real(rkind) :: mysum
        real(rkind) :: summation
        integer :: ierr

        mysum = sum(x)
        call MPI_Allreduce(mysum, summation, 1, mpirkind, MPI_MAX, MPI_COMM_WORLD, ierr)
        mean = summation/real(nproc,rkind)

    end function

    function P_MAXVAL_arr3(x) result(maximum)
        real(rkind), dimension(:,:,:), intent(in) :: x
        real(rkind) :: maximum
        real(rkind) :: mymax
        integer :: ierr

        mymax = MAXVAL(x)
        call MPI_Allreduce(mymax, maximum, 1, mpirkind, MPI_MAX, MPI_COMM_WORLD, ierr)

    end function

    function P_MAXVAL_sca(x) result(maximum)
        real(rkind), intent(in) :: x
        real(rkind) :: maximum
        integer :: ierr

        call MPI_Allreduce(x, maximum, 1, mpirkind, MPI_MAX, MPI_COMM_WORLD, ierr)

    end function
    
    function P_MINVAL_double(x) result(minimum)
        real(rkind), dimension(:,:,:), intent(in) :: x
        real(rkind) :: minimum
        real(rkind) :: mymin
        integer :: ierr

        mymin = MINVAL(x)
        call MPI_Allreduce(mymin, minimum, 1, mpirkind, MPI_MAX, MPI_COMM_WORLD, ierr)

    end function

    function P_SUM_double(x) result(summation)
        real(rkind), dimension(:,:,:), intent(in) :: x
        real(rkind) :: summation
        real(rkind) :: mysum
        integer :: ierr

        mysum = SUM(x)
        call MPI_Allreduce(mysum, summation, 1, mpirkind, MPI_MAX, MPI_COMM_WORLD, ierr)

    end function

    subroutine P_AVGZ(gp, f, avg)
        type(decomp_info), intent(in) :: gp
        real(rkind), dimension(:,:,:), intent(in) :: f
        real(rkind), dimension(size(f,1), size(f,2)), intent(out) :: avg
        real(rkind), dimension(size(f,1), size(f,2))              :: p_avg
        integer :: YZ_COMM, Z_COMM
        integer :: ierr
        
        if(size(f,1) == gp%xsz(1)) then        ! X-pencil
            call MPI_COMM_SPLIT(MPI_COMM_WORLD, gp%xst(1), nrank, YZ_COMM, ierr)
            call MPI_COMM_SPLIT( YZ_COMM,       gp%xst(2), nrank,  Z_COMM, ierr)
        else if (size(f,2) == gp%ysz(2)) then  ! Y-pencil
            call MPI_COMM_SPLIT(MPI_COMM_WORLD, gp%yst(1), nrank, YZ_COMM, ierr)
            call MPI_COMM_SPLIT( YZ_COMM,       gp%yst(2), nrank,  Z_COMM, ierr)
        else if (size(f,3) == gp%zsz(3)) then  ! Z-pencil
            avg = SUM(f,3) / real(gp%zsz(3),rkind)
            return
        else
            call GracefulExit("In P_AVGZ: Input array does not belong to any decomposition",256)
        end if

        p_avg = SUM(f,3)
        call MPI_ALLREDUCE(p_avg, avg, size(f,1)*size(f,2), mpirkind, MPI_SUM, Z_COMM, ierr)
        avg = avg / real(gp%zsz(3),rkind)

    end subroutine


end module
