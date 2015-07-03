! Module that contains all the derivative subroutines and drivers

module Derivatives

    use kind_parameters, only: rkind
    use fftstuff, only: ffts
    use dctstuff, only: dcts
    implicit none

    private

    real(rkind), parameter :: one=1.0_rkind

    integer          :: nx
    integer          :: ny
    integer          :: nz
    real(rkind)      :: dx=one
    real(rkind)      :: dy=one
    real(rkind)      :: dz=one
   
    real(rkind)      :: onebydx 
    real(rkind)      :: onebydy 
    real(rkind)      :: onebydz
    
    integer          :: bcx1=0
    integer          :: bcxn=0
    integer          :: bcy1=0
    integer          :: bcyn=0
    integer          :: bcz1=0
    integer          :: bczn=0
    
    character(len=4) :: methodx           ! Scheme to use in x. 'CD06', 'CD10', 'FOUR', 'CHEB'
    character(len=4) :: methody           ! Scheme to use in y. 'CD06', 'CD10', 'FOUR', 'CHEB'
    character(len=4) :: methodz           ! Scheme to use in z. 'CD06', 'CD10', 'FOUR', 'CHEB'
    
    logical          :: periodicx=.TRUE.
    logical          :: periodicy=.TRUE.
    logical          :: periodicz=.TRUE.

    real(rkind), allocatable, dimension(:,:) :: LU06X1, LU06Y1, LU06Z1  ! 6th order LU decomp matrices for first derivative
    real(rkind), allocatable, dimension(:,:) :: LU06X2, LU06Y2, LU06Z2  ! 6th order LU decomp matrices for second derivative
    
    real(rkind), allocatable, dimension(:,:) :: LU10X1, LU10Y1, LU10Z1  ! 10th order LU decomp matrices for first derivative
    real(rkind), allocatable, dimension(:,:) :: LU10X2, LU10Y2, LU10Z2  ! 10th order LU decomp matrices for second derivative

    type(ffts) :: xfft, yfft, zfft                                      ! FFT objects for each direction
    real(rkind), allocatable, dimension(:) :: k1, k2, k3                ! Wavenumbers

    type(dcts) :: xdct, ydct, zdct                                      ! DCT objects for each direction
    
    include 'PadeCoeffs.F90'

    public :: InitializeDerivatives

contains

    function InitializeDerivatives(nx_, ny_, nz_,                          &
                                     dx_, dy_, dz_,                        &
                                     periodicx_,periodicy_,periodicz_,     &
                                     methodx_,methody_,methodz_,           &
                                     bcx1_,bcxn_,bcy1_,bcyn_,bcz1_,bczn_ ) &
                                     result ( ierr )
        
        integer         , intent(in)           :: nx_,ny_,nz_
        real(rkind)     , intent(in), optional :: dx_,dy_,dz_
        integer         , intent(in), optional :: bcx1_,bcxn_,bcy1_,bcyn_,bcz1_,bczn_
        character(len=4), intent(in), optional :: methodx_,methody_,methodz_
        logical         , intent(in), optional :: periodicx_,periodicy_,periodicz_
        integer         , intent(out)          :: ierr

        nx = nx_; ny = ny_; nz = nz_

        if (present(dx_)) dx = dx_
        if (present(dy_)) dy = dy_
        if (present(dz_)) dz = dz_
       
        onebydx = one/dx 
        onebydy = one/dy 
        onebydz = one/dz 
        
        if (present(periodicx_)) periodicx = periodicx_
        if (present(periodicy_)) periodicy = periodicy_
        if (present(periodicz_)) periodicz = periodicz_
        
        if (present(methodx_)) methodx = methodx_
        if (present(methody_)) methody = methody_
        if (present(methodz_)) methodz = methodz_
        
        if (present(bcx1_)) bcx1 = bcx1_
        if (present(bcxn_)) bcxn = bcxn_
        if (present(bcy1_)) bcy1 = bcy1_
        if (present(bcyn_)) bcyn = bcyn_
        if (present(bcz1_)) bcz1 = bcz1_
        if (present(bczn_)) bczn = bczn_

        ! Initialize LU matrices for both 6th and 10th order
        ierr = InitLU06()
        if (ierr .NE. 0) return
        ierr = InitLU10()
        if (ierr .NE. 0) return
        
        ! Initialize FFT stuff
        ierr = InitFFT()
        if (ierr .NE. 0) return

        ! Initialize DCT stuff
        ierr = InitDCT()
        if (ierr .NE. 0) return
        
        select case (methodx)
        case ('CD06')
            ierr = InitCD06('x')
            if (ierr .NE. 0) return
        case ('CD10')
            ierr = InitCD10('x')
            if (ierr .NE. 0) return
        case ('FOUR')
            ierr = InitFOUR('x')
            if (ierr .NE. 0) return
        case ('CHEB')
            ierr = InitCHEB('x')
            if (ierr .NE. 0) return
        case default
            ierr = 1
            return
        end select
        
        select case (methody)
        case ('CD06')
            ierr = InitCD06('y')
            if (ierr .NE. 0) return
        case ('CD10')
            ierr = InitCD10('y')
            if (ierr .NE. 0) return
        case ('FOUR')
            ierr = InitFOUR('y')
            if (ierr .NE. 0) return
        case ('CHEB')
            ierr = InitCHEB('y')
            if (ierr .NE. 0) return
        case default
            ierr = 1
            return
        end select
        
        select case (methodz)
        case ('CD06')
            ierr = InitCD06('z')
            if (ierr .NE. 0) return
        case ('CD10')
            ierr = InitCD10('z')
            if (ierr .NE. 0) return
        case ('FOUR')
            ierr = InitFOUR('z')
            if (ierr .NE. 0) return
        case ('CHEB')
            ierr = InitCHEB('z')
            if (ierr .NE. 0) return
        case default
            ierr = 1
            return
        end select
        
        ierr = 0

    end function

end module
