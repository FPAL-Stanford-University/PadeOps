module reductions
    use kind_parameters, only: rkind,mpirkind
    use decomp_2d, only: nrank, nproc
    use mpi
    
    implicit none

    private
    public :: P_MAXVAL, P_MINVAL, P_SUM, P_MEAN

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

end module
