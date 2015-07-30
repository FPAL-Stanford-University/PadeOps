module derivatives
    use kind_parameters, only: rkind
    use cd10stuff, only: cd10
    use cd06stuff, only: cd06
    use fftstuff, only: ffts
    use dctstuff, only: dcts
    use exits,    only: gracefulExit, message

    implicit none
    private 
    public :: derivatives

    type derivatives
        private
        character(len=*) :: xmethod, ymethod, zmethod
        type(cd10) :: xcd10, ycd10, zcd10
        type(cd06) :: xcd06, ycd06, zcd06
        type(ffts) :: xfour, yfour, zfour
        type(dcts) :: xcheb, ycheb, zcheb 
        
        logical    :: curvilinear = .false. 
        logical    :: xmetric = .false.
        logical    :: ymetric = .false.
        logical    :: zmetric = .false.

        logical    :: initialized = .false. 
        contains
            procedure :: init
            procedure :: destroy
            procedure :: ddx
            procedure :: ddy
            procedure :: ddz
            procedure :: d2dx2
            procedure :: d2dy2
            procedure :: d2dz2
    end type


contains

    subroutine init(this,       nx,         ny,         nz, &
                                dx,         dy,         dz, &
                        periodic_x, periodic_y, periodic_z, &
                          method_x,   method_y,   method_z, &
                           xmetric,    ymetric,    zmetric, &
                       curvilinear)
        
        class(derivatives), intent(inout)          :: this
        integer           , intent(in)             :: nx, ny, nz  
        real(rkind)       , intent(in)             :: dx, dy, dz
        logical           , intent(in)             :: periodic_x, periodic_y, periodic_z
        character(len=*)  , intent(in),   optional :: method_x
        character(len=*)  , intent(in),   optional :: method_y
        character(len=*)  , intent(in),   optional :: method_z
        logical           , intent(in),   optional :: xmetric
        logical           , intent(in),   optional :: ymetric
        logical           , intent(in),   optional :: zmetric
        logical           , intent(in),   optional :: curvilinear
        
        integer :: ierr

       
        if (this%initialized) then
            call message("WARNING: Reinitializing the DERIVATIVE class!")
        end if  
        
        ! X direction

        if (present(method_x)) then
            select case (method_x)
            case ("cd10")
                ierr = xcd10%init( nx, dx, periodic_x, 0, 0)
                if (ierr .ne. 0) then
                    call GracefulExit("Initializing cd10 failed in X ",11)
                end if 
            case ("cd06")
                ierr = xcd06%init( nx, dx, periodic_x, 0, 0)
                if (ierr .ne. 0) then
                    call GracefulExit("Initializing cd06 failed in X ",12)
                end if 
            case ("four")
                ierr = xfour%init( nx, "x", ny, nz, dx)
                if (ierr .ne. 0) then
                    call GracefulExit("Initializing Fourier Collocation (FOUR) failed in X ",13)
                end if 
            case ("cheb")
                ierr = xcheb%init( nx )
                if (ierr .ne. 0) then
                    call GracefulExit("Initializing Chebyshev Collocation (CHEB) failed in X ",14)
                end if 
            case default 
                call GracefulExit("Invalid method selected in x direction ",01)
            end select
        else
            ierr = xcd10%init( nx, dx, periodic_x, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cd10 (default) failed in X ",11)
            end if 
        end if 
       
        ! Y direction

        if (present(method_y)) then
            select case (method_y)
            case ("cd10")
                ierr = ycd10%init( ny, dy, periodic_y, 0, 0)
                if (ierr .ne. 0) then
                    call GracefulExit("Initializing cd10 failed in Y ",11)
                end if 
            case ("cd06")
                ierr = ycd06%init( ny, dy, periodic_y, 0, 0)
                if (ierr .ne. 0) then
                    call GracefulExit("Initializing cd06 failed in Y",12)
                end if 
            case ("four")
                ierr = yfour%init( ny, "y", nx, nz, dy)
                if (ierr .ne. 0) then
                    call GracefulExit("Initializing Fourier Collocation (FOUR) failed in Y ",13)
                end if 
            case ("cheb")
                ierr = ycheb%init( ny )
                if (ierr .ne. 0) then
                    call GracefulExit("Initializing Chebyshev Collocation (CHEB) failed in Y",14)
                end if 
            case default 
                call GracefulExit("Invalid method selected in y direction",01)
            end select
        else
            ierr = ycd10%init( ny, dy, periodic_y, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cd10 (default) failed in Y ",11)
            end if 
        end if 


        ! Z direction

        if (present(method_z)) then
            select case (method_z)
            case ("cd10")
                ierr = zcd10%init( nz, dz, periodic_z, 0, 0)
                if (ierr .ne. 0) then
                    call GracefulExit("Initializing cd10 failed in Z ",11)
                end if 
            case ("cd06")
                ierr = zcd06%init( nz, dz, periodic_z, 0, 0)
                if (ierr .ne. 0) then
                    call GracefulExit("Initializing cd06 failed in Z",12)
                end if 
            case ("four")
                ierr = zfour%init( nz, "z", nx, ny, dz)
                if (ierr .ne. 0) then
                    call GracefulExit("Initializing Fourier Collocation (FOUR) failed in Z ",13)
                end if 
            case ("cheb")
                ierr = zcheb%init( nz )
                if (ierr .ne. 0) then
                    call GracefulExit("Initializing Chebyshev Collocation (CHEB) failed in Z",14)
                end if 
            case default 
                call GracefulExit("Invalid method selected in z direction",01)
            end select
        else
            ierr = zcd10%init( nz, dz, periodic_z, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cd10 (default) failed in Z ",11)
        end if 

        ! Metrics

        if (present(xMetric)) then
            if (xMetric) then
                call GracefulExit("Code INCOMPLETE! Metric terms are not &
                 &   supported in x direction!",04)
            end if
        end if 

        if (present(yMetric)) then
            if (yMetric) then
                call GracefulExit("Code INCOMPLETE! Metric terms are not &
                 &   supported in y direction!",05)
            end if
        end if 
        
        if (present(zMetric)) then
            if (zMetric) then
                call GracefulExit("Code INCOMPLETE! Metric terms are not &
                 &   supported in z direction!",06)
            end if
        end if 
        
        if (present(Curvilinear)) then
            if (Curvilinear) then
                call GracefulExit("Code INCOMPLETE! Curvilinear Formulation is &
                & not supported!",06)
            end if
        end if 

        this%xmethod = x_method 
        this%ymethod = y_method 
        this%zmethod = z_method 

    end subroutine

    
    
    
end module derivatives
