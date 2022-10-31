! Set the floating point precision

module kind_parameters

    use mpi
    implicit none
    
    private
    public :: single_kind, double_kind, rkind, mpirkind, mpickind, clen, &
      stdin, stdout, stderr, castSingle

    integer, parameter :: single_kind = kind(0.0)
    integer, parameter :: double_kind = kind(0.d0)
    integer, parameter :: rkind = double_kind
    integer, parameter :: mpirkind = MPI_DOUBLE_PRECISION
    integer, parameter :: mpickind = MPI_DOUBLE_COMPLEX
   
    integer, parameter :: clen = 256

    integer, parameter :: stdin  = 5
    integer, parameter :: stdout = 6
    integer, parameter :: stderr = 0

    interface castSingle
      module procedure castSingleDouble, castSingleInt
    end interface
  
    contains
      
        function castSingleDouble(fdouble) result(fsingle)
            real(rkind), intent(in) :: fdouble
            real(single_kind) :: fsingle

            fsingle = real(fdouble,single_kind)
        end function 
        
        function castSingleInt(fInt) result(fsingle)
            integer, intent(in) :: fInt
            real(single_kind) :: fsingle

            fsingle = real(fInt,single_kind)
        end function 
   
end module
