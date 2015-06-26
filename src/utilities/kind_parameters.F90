! Set the floating point precision

module kind_parameters

    use ISO_FORTRAN_ENV, only: real32, real64
    implicit none
    
    private
    public :: rkind

    integer, parameter :: rkind=real64

end module
