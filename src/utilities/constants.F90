module constants

    use kind_parameters, only: rkind
    implicit none

    real(rkind), parameter :: zero = 0._rkind
    real(rkind), parameter :: one = 1._rkind
    real(rkind), parameter :: two = 2._rkind
    real(rkind), parameter :: three = 3._rkind
    real(rkind), parameter :: four = 4._rkind
    real(rkind), parameter :: five = 5._rkind
    real(rkind), parameter :: six = 6._rkind
    real(rkind), parameter :: seven = 7._rkind
    real(rkind), parameter :: eight = 8._rkind
    real(rkind), parameter :: nine = 9._rkind
    
    real(rkind), parameter :: half = 1._rkind/2._rkind
    real(rkind), parameter :: third = 1._rkind/3._rkind
    real(rkind), parameter :: twothird = 2._rkind/3._rkind

    real(rkind), parameter :: pi=4._rkind*atan(1._rkind)
    
    complex(rkind), parameter :: imi=(zero,one)

    real(rkind), parameter :: r_eps = real(1.d-15,rkind)
end module
