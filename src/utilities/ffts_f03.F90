module fftstuff_f03
    use, intrinsic :: ISO_C_BINDING
    use kind_parameters, only: rkind
    use constants, only: zero,one,half,two,pi,imi

    implicit none
    
    private
    public :: ffts_f03

    include "fftw3.f03"

    type ffts_f03

        integer                                :: fftw_plan_type = FFTW_EXHAUSTIVE
        type(C_PTR)                            :: plan_fwd
        type(C_PTR)                            :: plan_bwd

        real(rkind), dimension(:), allocatable :: k1d

        integer                                :: split
        integer                                :: n

        integer                                :: nx, ny, nz
        character(len=1)                       :: dir = "x"

        contains
        procedure :: init
        !procedure :: destroy
        !procedure :: dd1
        !procedure :: dd2
        !procedure :: dd2


        !procedure, private :: fft1d
        !procedure, private :: ifft1d
        !procedure, private :: 0ifyOddball 
        !procedure, private :: genK
        !procedure, private :: create_Plans

    end type



contains

    function init(this,n_,n1_,n2_,dir_) result (ierr)
        class(ffts_f03) , intent(inout) :: this
        integer         , intent(in   ) :: n_, n1_, n2_
        character(len=1), intent(in   ) :: dir_
        integer                         :: ierr 
        
        real(rkind)    , pointer       :: a1(:,:,:)
        complex(rkind) , pointer       :: a2(:,:,:)
        type(C_PTR)                     :: a1_p, a2_p
        integer(C_SIZE_T)               :: sz   
        
        
        ! Check if even 
        if (mod(n_,2) .NE. 0) then
            ierr = 20
            return
        end if 

        this%n = n_
       
        this%split = this%n/2 + 1 

        this%dir = dir_

        select case (dir_)
        case ("x")
            this%ny = n1_
            this%nz = n2_
            this%nx = 0
            
            ! Input data (real)
            sz = this%n*n1_*n2_
            a1_p = fftw_alloc_real(sz)
            call c_f_pointer(a1_p,a1, [n_,this%ny,this%nz])
        
            ! Output data (complex)
            sz = this%split*n1_*n2_
            a2_p = fftw_alloc_complex(sz)
            call c_f_pointer(a2_p,a2,[this%split,this%ny,this%nz])

            !this%plan_fwd =  fftw_plan_many_dft_r2c(1, this%n, this%ny*this%nz, a1, this%n, 1,this%n, a2, this%split, 1,this%split, this%fftw_plan_type)

        case ("y")
            this%nx = n1_
            this%nz = n2_
            this%ny = 0
        case ("z")
            this%nx = n1_
            this%ny = n2_
            this%nz = 0
        case default 
            ierr = 20
            return 
        end select



    end function




end module 
