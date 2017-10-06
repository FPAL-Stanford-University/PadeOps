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
    real(rkind), parameter :: ten = 10._rkind

    real(rkind), parameter :: half = 1._rkind/2._rkind
    real(rkind), parameter :: third = 1._rkind/3._rkind
    real(rkind), parameter :: fourth = 1._rkind/4._rkind
    real(rkind), parameter :: fifth = 1._rkind/5._rkind
    real(rkind), parameter :: sixth = 1._rkind/6._rkind
    real(rkind), parameter :: seventh = 1._rkind/7._rkind
    real(rkind), parameter :: eighth = 1._rkind/8._rkind

    real(rkind), parameter :: twothird = 2._rkind/3._rkind
    real(rkind), parameter :: fourthird = 4._rkind/3._rkind

    real(rkind), parameter :: pi=4._rkind*atan(1._rkind)
    
    complex(rkind), parameter :: imi=(zero,one)
    complex(rkind), parameter :: im0=(zero,zero)

    real(rkind), parameter :: kappa = 0.40_rkind ! Kappa used in Turbulent BL wall models
    real(rkind), parameter :: eps = epsilon(real(1.0,rkind))
end module
