module fftstuff_f03
    use, intrinsic :: ISO_C_BINDING
    use kind_parameters, only: rkind
    use constants, only: zero,one,half,two,pi,imi

    implicit none
    
    private
    public :: ffts_f03

    include "fftw3.f03"

    type ffts_f03

        type(C_PTR)                            :: plan_fwd
        type(C_PTR)                            :: plan_bwd

        real(rkind), dimension(:), allocatable :: k1d

        real(rkind)                            :: split
        integer                                :: n

        integer                                :: nx, ny, nz
        character(len=1)                       :: dir = "x"

        !contains
        !procedure :: init
        !procedure :: destroy
        !procedure :: dd1
        !procedure :: dd2
        !procedure :: dd2


        !procedure, private :: fft1d
        !procedure, private :: ifft1d
        !procedure, private :: nullifyOddball 
        !procedure, private :: genK
        !procedure, private :: create_Plans

    end type







end module 
