! Module that contains all the derivative subroutines and drivers

module derivativestuff

    use kind_parameters, only: rkind
    use fftstuff, only: ffts
    use dctstuff, only: dcts
    use cd10stuff, only: cd10
    use cd06stuff, only: cd06
    use constants, only: one
    use exits, only: GracefulExit
    implicit none

    private

    public :: derivatives
    
    real(rkind), allocatable, dimension(:) :: k1, k2, k3                ! Wavenumbers
    
    
    type :: derivatives

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

        integer          :: xoprank=1, yoprank=2, zoprank=3


        type(ffts)      :: xfft, yfft, zfft                                      ! FFT objects for each direction
        type(dcts)      :: xdct, ydct, zdct                                      ! DCT objects for each direction
        type(cd10)      :: xcd10, ycd10, zcd10                                   ! CD10 objects for each direction  
        type(cd06)      :: xcd06, ycd06, zcd06                                   ! CD06 objects for each direction  
       
       
        contains

            !! Public Procedures
            procedure :: init
            procedure :: SetXoprank
            procedure :: SetYoprank
            procedure :: SetZoprank
            
            procedure   :: ddx 
            procedure   :: ddy 
            procedure   :: ddz 
            
            procedure   :: d2dx2 
            procedure   :: d2dy2 
            procedure   :: d2dz2 
            
            !! Private Procedures
            ! direction independent procedures
            !procedure, private :: chebder
            !procedure, private :: fourder
            !procedure, private :: cd10der
            !procedure, private :: cd06der

            ! direction dependent procedures
            ! First order derivatives
            procedure, private :: ddx_cheb
            procedure, private :: ddx_four
            procedure, private :: ddx_cd10
            procedure, private :: ddx_cd06
            procedure, private :: ddy_cheb
            procedure, private :: ddy_four
            procedure, private :: ddy_cd10
            procedure, private :: ddy_cd06
            procedure, private :: ddz_cheb
            procedure, private :: ddz_four
            procedure, private :: ddz_cd10
            procedure, private :: ddz_cd06

            ! Second order derivatives
            procedure, private :: d2dx2_cheb
            procedure, private :: d2dx2_four
            procedure, private :: d2dx2_cd10
            procedure, private :: d2dx2_cd06
            procedure, private :: d2dy2_cheb
            procedure, private :: d2dy2_four
            procedure, private :: d2dy2_cd10
            procedure, private :: d2dy2_cd06
            procedure, private :: d2dz2_cheb
            procedure, private :: d2dz2_four
            procedure, private :: d2dz2_cd10
    end type 
    


contains
    
    include "FourierCollocation.F90"
    include "ChebyshevCollocation.F90"
    include "CompactDifference06.F90"
    include "CompactDifference10.F90"
    
    function init(this, nx_, ny_, nz_,                          &
                  dx_, dy_, dz_,                        &
                  periodicx_,periodicy_,periodicz_,     &
                  methodx_,methody_,methodz_,           &
                  bcx1_,bcxn_,bcy1_,bcyn_,bcz1_,bczn_ ) &
                  result ( ierr )
        
        class(derivatives), intent(inout)      :: this
        integer         , intent(in)           :: nx_,ny_,nz_
        real(rkind)     , intent(in), optional :: dx_,dy_,dz_
        integer         , intent(in), optional :: bcx1_,bcxn_,bcy1_,bcyn_,bcz1_,bczn_
        character(len=4), intent(in), optional :: methodx_,methody_,methodz_
        logical         , intent(in), optional :: periodicx_,periodicy_,periodicz_
        integer         , intent(out)          :: ierr

        this%nx = nx_; this%ny = ny_; this%nz = nz_

        if (present(dx_)) this%dx = dx_
        if (present(dy_)) this%dy = dy_
        if (present(dz_)) this%dz = dz_
       
        this%onebydx = one/this%dx 
        this%onebydy = one/this%dy 
        this%onebydz = one/this%dz 
        
        if (present(periodicx_)) this%periodicx = periodicx_
        if (present(periodicy_)) this%periodicy = periodicy_
        if (present(periodicz_)) this%periodicz = periodicz_
        
        if (present(methodx_)) this%methodx = methodx_
        if (present(methody_)) this%methody = methody_
        if (present(methodz_)) this%methodz = methodz_
        
        if (present(bcx1_)) this%bcx1 = bcx1_
        if (present(bcxn_)) this%bcxn = bcxn_
        if (present(bcy1_)) this%bcy1 = bcy1_
        if (present(bcyn_)) this%bcyn = bcyn_
        if (present(bcz1_)) this%bcz1 = bcz1_
        if (present(bczn_)) this%bczn = bczn_

        ! Initialize 10th order CD classes
        ierr = xcd10%init(nx_, dx_, periodicx_, bcx1_, bcx2_)
        if (ierr .NE. 0) return
        ierr = ycd10%init(ny_, dy_, periodicy_, bcy1_, bcy2_)
        if (ierr .NE. 0) return
        ierr = zcd10%init(nz_, dz_, periodicz_, bcz1_, bcz2_)
        if (ierr .NE. 0) return

        ! Initialize  6th order CD classes
        ierr = xcd06%init(nx_, dx_, periodicx_, bcx1_, bcx2_)
        if (ierr .NE. 0) return
        ierr = ycd06%init(ny_, dy_, periodicy_, bcy1_, bcy2_)
        if (ierr .NE. 0) return
        ierr = zcd06%init(nz_, dz_, periodicz_, bcz1_, bcz2_)
        if (ierr .NE. 0) return

        ! Initialize FFT classes
        ierr = xfft % init(nx,dx)
        if (ierr .NE. 0) return
        ierr = yfft % init(ny,dy)
        if (ierr .NE. 0) return
        ierr = zfft % init(nz,dz)
        if (ierr .NE. 0) return

        ! Initialize FFT classes
        ierr = xdct % init(nx)
        if (ierr .NE. 0) return
        ierr = ydct % init(ny)
        if (ierr .NE. 0) return
        ierr = zdct % init(nz)
        if (ierr .NE. 0) return
        
        ierr = 0

    end function

    subroutine SetXoprank(this,xop)
        class(derivatives), intent(inout) :: this
        integer, intent(in)               :: xop

        this%xoprank = xop
    end subroutine

    subroutine SetYoprank(this,yop)
        class(derivatives), intent(inout) :: this
        integer, intent(in) :: yop

        this%yoprank = yop
    end subroutine

    subroutine SetZoprank(this,zop)
        integer, intent(in) :: zop

        this%zoprank = zop
    end subroutine


    function ddx(this,f) result(df)
        class(derivatives), intent(in) :: this
        real(rkind), intent(in), dimension(:,:,:), f
        real(rkind), dimension(size(f,1),size(f,2),size(f,3)) :: df
       
        select case (this%methodx)
        case ("CD06") 
            call this%ddx_cd06(f)
        case ("CD10") 
            call this%ddx_cd10(f)
        case ("FOUR") 
            call this%ddx_four(f)
        case ("CHEB") 
            call this%ddx_cheb(f)
        case default
            call gracefulExit("You intialized the derivative class (X direction) using an &
            incorrect numerical method. Reinitialize using one of the following  &
            methods: 1) CD06, 2) CD10, 3) FOUR, 4) CHEB", 10)
        end select
    end function

    function ddy(this,f) result(df)
        class(derivatives), intent(in) :: this
        real(rkind), intent(in), dimension(:,:,:), f
        real(rkind), dimension(size(f,1),size(f,2),size(f,3)) :: df
       
        select case (this%methody)
        case ("CD06") 
            call this%ddy_cd06(f)
        case ("CD10") 
            call this%ddy_cd10(f)
        case ("FOUR") 
            call this%ddy_four(f)
        case ("CHEB") 
            call this%ddy_cheb(f)
        case default
            call gracefulExit("You intialized the derivative class (Y direction) using an &
            incorrect numerical method. Reinitialize using one of the following  &
            methods: 1) CD06, 2) CD10, 3) FOUR, 4) CHEB", 10)
        end select
    end function

    function ddz(this,f) result(df)
        class(derivatives), intent(in) :: this
        real(rkind), intent(in), dimension(:,:,:), f
        real(rkind), dimension(size(f,1),size(f,2),size(f,3)) :: df
       
        select case (this%methodz)
        case ("CD06") 
            call this%ddz_cd06(f)
        case ("CD10") 
            call this%ddz_cd10(f)
        case ("FOUR") 
            call this%ddz_four(f)
        case ("CHEB") 
            call this%ddz_cheb(f)
        case default
            call gracefulExit("You intialized the derivative class (Z direction) using an &
            incorrect numerical method. Reinitialize using one of the following  &
            methods: 1) CD06, 2) CD10, 3) FOUR, 4) CHEB", 10)
        end select
    end function

    function d2dx2(this,f) result(df)
        class(derivatives), intent(in) :: this
        real(rkind), intent(in), dimension(:,:,:), f
        real(rkind), dimension(size(f,1),size(f,2),size(f,3)) :: df
       
        select case (this%methodx)
        case ("CD06") 
            call this%d2dx2_cd06(f)
        case ("CD10") 
            call this%d2dx2_cd10(f)
        case ("FOUR") 
            call this%d2dx2_four(f)
        case ("CHEB") 
            call this%d2dx2_cheb(f)
        case default
            call gracefulExit("You intialized the derivative class (X direction) using an &
            incorrect numerical method. Reinitialize using one of the following  &
            methods: 1) CD06, 2) CD10, 3) FOUR, 4) CHEB", 10)
        end select
    end function

    function d2dy2(this,f) result(df)
        class(derivatives), intent(in) :: this
        real(rkind), intent(in), dimension(:,:,:), f
        real(rkind), dimension(size(f,1),size(f,2),size(f,3)) :: df
       
        select case (this%methody)
        case ("CD06") 
            call this%d2dy2_cd06(f)
        case ("CD10") 
            call this%d2dy2_cd10(f)
        case ("FOUR") 
            call this%d2dy2_four(f)
        case ("CHEB") 
            call this%d2dy2_cheb(f)
        case default
            call gracefulExit("You intialized the derivative class (Y direction) using an &
            incorrect numerical method. Reinitialize using one of the following  &
            methods: 1) CD06, 2) CD10, 3) FOUR, 4) CHEB", 10)
        end select
    end function

    function d2dz2(this,f) result(df)
        class(derivatives), intent(in) :: this
        real(rkind), intent(in), dimension(:,:,:), f
        real(rkind), dimension(size(f,1),size(f,2),size(f,3)) :: df
       
        select case (this%methodz)
        case ("CD06") 
            call this%d2dz2_cd06(f)
        case ("CD10") 
            call this%d2dz2_cd10(f)
        case ("FOUR") 
            call this%d2dz2_four(f)
        case ("CHEB") 
            call this%d2dz2_cheb(f)
        case default
            call gracefulExit("You intialized the derivative class (Z direction) using an &
            incorrect numerical method. Reinitialize using one of the following  &
            methods: 1) CD06, 2) CD10, 3) FOUR, 4) CHEB", 10)
        end select
    end function


    subroutine check_dimension(this,arr,dir) result(isCorrect)
        class(derivatives), intent(in) :: this
    subroutine check_dimension(this,arr,dir) result(isCorrect)
        class(derivatives), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in) :: arr
        character(len=1), intent(in) :: dir

        select case (dir)
        case ("x")
            if (size(arr,this%xoprank) .ne. this%nx) then
                call gracefulExit("Dimenion of input vector for direction X is
                inconsistent with initialization")
            end if 
        case ("y")
            if (size(arr,this%yoprank) .ne. this%ny) then
                call gracefulExit("Dimenion of input vector for direction Y is
                inconsistent with initialization")
            end if 
            if (size(arr,this%zoprank) .ne. this%nz) then
                call gracefulExit("Dimenion of input vector for direction Z is
                inconsistent with initialization")
            end if 
        case default 
            call GracefulExit ( "Incorrect String entered for SUBOUTINE: &
                        check_dimension", 5)
        end select

    end function 

end module
