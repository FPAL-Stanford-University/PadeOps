! Set the floating point precision

module kind_parameters

    use mpi
    implicit none
    
    private
    public :: rkind, mpirkind, mpickind, clen, stdin, stdout, stderr

    integer, parameter :: rkind=kind(0.d0)
    integer, parameter :: mpirkind = MPI_DOUBLE_PRECISION
    integer, parameter :: mpickind = MPI_DOUBLE_COMPLEX
   
    integer, parameter :: clen = 100

    integer, parameter :: stdin  = 5
    integer, parameter :: stdout = 6
    integer, parameter :: stderr = 0

end module
