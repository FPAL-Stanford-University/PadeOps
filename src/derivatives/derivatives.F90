module DerivativesMod
    use kind_parameters, only: rkind, clen
    use cd10stuff, only: cd10
    use cd06stuff, only: cd06
    use fftstuff, only: ffts
    use dctstuff, only: dcts
    use exits,    only: gracefulExit, message
    use decomp_2d, only: decomp_info, nrank
    use constants, only: zero, one, half

    implicit none
    private 
    public :: derivatives

    type derivatives
        private
       
        
        ! Available Methods: Specific to each direction
        ! 1: "CD10", 2: "CD06", 3: "FOUR", 4: "CHEB"
        integer                :: xmethod, ymethod, zmethod 

        type(cd10), allocatable :: xcd10, ycd10, zcd10 
        type(cd06), allocatable :: xcd06, ycd06, zcd06 
        type(ffts), allocatable :: xfour, yfour, zfour
        type(dcts), allocatable :: xcheb, ycheb, zcheb 
       
    
        logical                        :: curvilinear = .false. 
        logical                        :: xmetric = .false.
        logical                        :: ymetric = .false.
        logical                        :: zmetric = .false.

        real(rkind), allocatable, dimension(:,:,:) :: dxidx,   dxidx_sq,   d2xidx2,   dfdxbuf
        real(rkind), allocatable, dimension(:,:,:) :: detady,  detady_sq,  d2etady2,  dfdybuf
        real(rkind), allocatable, dimension(:,:,:) :: dzetadz, dzetadz_sq, d2zetadz2, dfdzbuf

        integer                        :: nxg, nyg, nzg ! Global sizes
        integer, dimension(3)          :: xsz, ysz, zsz ! Local decomposition sizes
        
        logical                        :: initialized = .false. 
        contains
            procedure, private :: init_parallel
            procedure, private :: init_serial
            procedure, private :: init_procedures
            procedure, private :: init_curvilinear
            procedure, private :: init_metric_1d
            generic :: init => init_parallel, init_serial
            procedure :: destroy
            procedure :: ddx
            procedure :: ddy
            procedure :: ddz
            procedure :: d2dx2
            procedure :: d2dy2
            procedure :: d2dz2
            procedure :: getMethodx
            procedure :: getMethody
            procedure :: getMethodz
            procedure :: getmetricx
            procedure :: getmetricy
            procedure :: getmetricz
            procedure :: getcurvilinear
            procedure :: set_xsz
            procedure :: set_ysz
            procedure :: set_zsz
    end type


contains

    subroutine set_xsz(this,xsz)
        class(derivatives), intent(inout) :: this
        integer, dimension(3), intent(in) :: xsz

        if (this%xsz(1) /= xsz(1)) then
            call GracefulExit("Cannot change xsz(1) in set_xsz",453)
        end if
        this%xsz = xsz
    end subroutine

    subroutine set_ysz(this,ysz)
        class(derivatives), intent(inout) :: this
        integer, dimension(3), intent(in) :: ysz

        if (this%ysz(2) /= ysz(2)) then
            call GracefulExit("Cannot change ysz(2) in set_ysz",453)
        end if
        this%ysz = ysz
    end subroutine

    subroutine set_zsz(this,zsz)
        class(derivatives), intent(inout) :: this
        integer, dimension(3), intent(in) :: zsz

        if (this%zsz(3) /= zsz(3)) then
            call GracefulExit("Cannot change zsz(3) in set_zsz",453)
        end if
        this%zsz = zsz
    end subroutine

    function getmetricx(this) result(m)
        class(derivatives), intent(in) :: this
        logical :: m 
        m = this%xmetric
    end function

    function getmetricy(this) result(m)
        class(derivatives), intent(in) :: this
        logical :: m 
        m = this%ymetric
    end function

    function getmetricz(this) result(m)
        class(derivatives), intent(in) :: this
        logical :: m 
        m = this%zmetric
    end function

    function getcurvilinear(this) result(m)
        class(derivatives), intent(in) :: this
        logical :: m 
        m = this%curvilinear
    end function

    function getMethodx(this) result(m)
        class(derivatives), intent(in) :: this
        character(len=clen) :: m 
        select case (this%xmethod)
        case (1)
            m = "cd10"
        case (2)
            m = "cd06"
        case (3)
            m = "four"
        case (4)
            m = "cheb"
        end select 

    end function

    function getMethody(this) result(m)
        class(derivatives), intent(in) :: this
        character(len=clen) :: m 
        select case (this%ymethod)
        case (1)
            m = "cd10"
        case (2)
            m = "cd06"
        case (3)
            m = "four"
        case (4)
            m = "cheb"
        end select 

    end function

    function getMethodz(this) result(m)
        class(derivatives), intent(in) :: this
        character(len=clen) :: m 
        select case (this%zmethod)
        case (1)
            m = "cd10"
        case (2)
            m = "cd06"
        case (3)
            m = "four"
        case (4)
            m = "cheb"
        end select 

    end function


    subroutine init_serial(this,nx,         ny,         nz, &
                                dx,         dy,         dz, &
                        periodic_x, periodic_y, periodic_z, &
                          method_x,   method_y,   method_z, &
                                 x,          y,          z, &
                           xmetric,    ymetric,    zmetric, &
                       curvilinear, inputfile, xi, eta, zeta)
                       
        
        class(derivatives), intent(inout)          :: this
        integer           , intent(in)             :: nx, ny, nz
        real(rkind)       , intent(in)             :: dx, dy, dz
        logical           , intent(in)             :: periodic_x, periodic_y, periodic_z
        character(len=*)  , intent(in)             :: method_x
        character(len=*)  , intent(in)             :: method_y
        character(len=*)  , intent(in)             :: method_z
        real(rkind),  intent(in), dimension(:,:,:) :: x, y, z
        logical           , intent(in),   optional :: xmetric
        logical           , intent(in),   optional :: ymetric
        logical           , intent(in),   optional :: zmetric
        logical           , intent(in),   optional :: curvilinear
        character(len=* ) ,intent(in),   optional :: inputfile
        real(rkind), intent(in), optional, dimension(:,:,:) :: xi, eta, zeta
        
       
        if (this%initialized) then
            call message("WARNING: Reinitializing the DERIVATIVE class!")
            call this%destroy
        end if  
      
        if (present(xmetric)) this%xmetric = xmetric
        if (present(ymetric)) this%ymetric = ymetric
        if (present(zmetric)) this%zmetric = zmetric

        if (present(curvilinear)) this%curvilinear = curvilinear

        this%nxg = nx
        this%nyg = ny
        this%nzg = nz
        this%xsz = [nx, ny, nz]
        this%ysz = this%xsz
        this%zsz = this%xsz

        call this%init_procedures (dx, dy, dz, periodic_x, periodic_y, periodic_z,&
                                            method_x,method_y,method_z)

        if(this%xmetric .or. this%ymetric .or. this%zmetric .or. this%curvilinear) then
            ! confirm that xi, eta, zeta are passed in
            if(.not. (present(inputfile))) call GracefulExit("inputfile missing in derivative", 04)
            if(.not. (present(xi)       )) call GracefulExit("xi        missing in derivative", 04)
            if(.not. (present(eta)      )) call GracefulExit("eta       missing in derivative", 04)
            if(.not. (present(zeta)     )) call GracefulExit("zeta      missing in derivative", 04)

            call this%init_curvilinear(this%xmetric, this%ymetric, this%zmetric, &
                              x, y, z, this%curvilinear, inputfile, xi, eta, zeta)
        endif

    end subroutine

    subroutine init_parallel(this,                      gp, &
                                dx,         dy,         dz, &
                        periodic_x, periodic_y, periodic_z, &
                          method_x,   method_y,   method_z, &
                                 x,          y,          z, &
                           xmetric,    ymetric,    zmetric, &
                       curvilinear, inputfile, xi, eta, zeta)
                       
        
        class(derivatives), intent(inout)          :: this
        class(decomp_info), intent(in)             :: gp 
        real(rkind)       , intent(in)             :: dx, dy, dz
        logical           , intent(in)             :: periodic_x, periodic_y, periodic_z
        character(len=*)  , intent(in)             :: method_x
        character(len=*)  , intent(in)             :: method_y
        character(len=*)  , intent(in)             :: method_z
        real(rkind),  intent(in), dimension(:,:,:) :: x, y, z
        logical           , intent(in),   optional :: xmetric
        logical           , intent(in),   optional :: ymetric
        logical           , intent(in),   optional :: zmetric
        logical           , intent(in),   optional :: curvilinear
        character(len=*)  , intent(in),   optional :: inputfile
        real(rkind), intent(in), optional, dimension(:,:,:) :: xi, eta, zeta

       
        if (this%initialized) then
            call message("WARNING: Reinitializing the DERIVATIVE class!")
            call this%destroy
        end if  
       
        if (present(xmetric)) this%xmetric = xmetric
        if (present(ymetric)) this%ymetric = ymetric
        if (present(zmetric)) this%zmetric = zmetric

        if (present(curvilinear)) this%curvilinear = curvilinear

        this%nxg = gp%xsz(1)
        this%nyg = gp%ysz(2)
        this%nzg = gp%zsz(3)
        this%xsz = gp%xsz
        this%ysz = gp%ysz
        this%zsz = gp%zsz

        call this%init_procedures (dx, dy, dz, periodic_x, periodic_y, periodic_z,&
                                            method_x,method_y,method_z)

        if(this%xmetric .or. this%ymetric .or. this%zmetric .or. this%curvilinear) then
            ! confirm that xi, eta, zeta are passed in
            if(.not. (present(inputfile))) call GracefulExit("inputfile missing in derivative", 04)
            if(.not. (present(xi)       )) call GracefulExit("xi        missing in derivative", 04)
            if(.not. (present(eta)      )) call GracefulExit("eta       missing in derivative", 04)
            if(.not. (present(zeta)     )) call GracefulExit("zeta      missing in derivative", 04)

            call this%init_curvilinear(this%xmetric, this%ymetric, this%zmetric, &
                              x, y, z, this%curvilinear, inputfile, xi, eta, zeta)
        endif

    end subroutine


    subroutine init_curvilinear(this, xmetric, ymetric, zmetric, x, y, z, curvilinear, inputfile, xi, eta, zeta)
        class(derivatives), intent(inout) :: this
        logical    , intent(in)           :: xmetric
        logical    , intent(in)           :: ymetric
        logical    , intent(in)           :: zmetric
        real(rkind), intent(in), dimension(:,:,:) :: x, y, z
        logical    , intent(in)           :: curvilinear
        character(len=*), intent(in)      :: inputfile
        real(rkind),      intent(in), dimension(:,:,:) :: xi, eta, zeta
       
        integer     ::  xmetric_flag,  ymetric_flag, zmetric_flag, ierr
        !real(rkind) :: xmetric_param, ymetric_param, zmetric_param
        real(rkind), allocatable, dimension(:,:) :: metric_params

        !! read input file
        !namelist /METRICS/ xmetric_flag, xmetric_param,    &
        !                   ymetric_flag, ymetric_param,    &
        !                   zmetric_flag, zmetric_param,    &
        !                   curviln_flag, curviln_param
        namelist /METRICS/ xmetric_flag, ymetric_flag, zmetric_flag, metric_params

        allocate(metric_params(3,5))    ! 3 :: (x,y,z); 5 :: max no of parameters
        metric_params = zero

        open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=123, NML=METRICS)
        close(123)

        ! Metrics
        if (xmetric) then
            call GracefulExit("Code INCOMPLETE! Metric terms are not &
             &   supported in x direction!",04)
             this%xmetric = .true. 
             !! allocate(this%dxidx)
        end if

        if (ymetric) then
            !call GracefulExit("Code INCOMPLETE! Metric terms are not &
            ! &   supported in y direction!",05)

             if(allocated(this%dfdxbuf)) deallocate(this%dfdxbuf)
             allocate(this%dfdxbuf(this%ysz(1), this%ysz(2), this%ysz(3)))

             if(allocated(this%detady)) deallocate(this%detady)
             allocate(this%detady(this%ysz(1), this%ysz(2), this%ysz(3)))

             if(allocated(this%detady_sq)) deallocate(this%detady_sq)
             allocate(this%detady_sq(this%ysz(1), this%ysz(2), this%ysz(3)))

             if(allocated(this%d2etady2)) deallocate(this%d2etady2)
             allocate(this%d2etady2(this%ysz(1), this%ysz(2), this%ysz(3)))

             if(allocated(this%dfdybuf)) deallocate(this%dfdybuf)
             allocate(this%dfdybuf(this%ysz(1), this%ysz(2), this%ysz(3)))

             call this%init_metric_1d(y, eta, ymetric_flag, metric_params(2,:), this%detady, this%detady_sq, this%d2etady2)
             !! stretched  :: [x, y, z]  --> uniform :: [xi, eta, zeta]]
             !! [dfdx, dfdy] --> [dfdxi*dxidx + dfdeta*detadx, dfdxi*dxidy + dfdeta*detady]
             !! [ystretch_type, ystretch_param]

             this%ymetric = .true. 
        end if
        
        if (zmetric) then
            call GracefulExit("Code INCOMPLETE! Metric terms are not &
             &   supported in z direction!",06)
             this%zmetric = .true. 
             !! allocate(this%dzetadz)
        end if
        
        if (curvilinear) then
            call GracefulExit("Code INCOMPLETE! Curvilinear Formulation is &
            & not supported!",06)
            this%curvilinear = .true. 
            !! if this is true, compute all 9 components of the Jacobian
        end if

    end subroutine

    subroutine init_procedures(this,dx, dy, dz, periodic_x, periodic_y, periodic_z,&
                                    method_x,method_y,method_z)
                                    
        class(derivatives), intent(inout)          :: this
        real(rkind)       , intent(in)             :: dx, dy, dz
        logical           , intent(in)             :: periodic_x, periodic_y, periodic_z
        character(len=*)  , intent(in)             :: method_x
        character(len=*)  , intent(in)             :: method_y
        character(len=*)  , intent(in)             :: method_z
        integer :: ierr

        select case (method_x)
        case ("cd10")
            allocate(this%xcd10)
            ierr = this % xcd10%init( this%xsz(1), dx, periodic_x, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cd10 failed in X ",11)
            end if
            this%xmethod = 1 
        case ("cd06")
            allocate(this%xcd06)
            ierr = this % xcd06%init( this%xsz(1), dx, periodic_x, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cd06 failed in X ",12)
            end if 
            this%xmethod = 2 
        case ("four")
            allocate(this%xfour)
            ierr = this % xfour%init( this%xsz(1), "x", this%xsz(2), this%xsz(3), dx)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing Fourier Collocation (FOUR) failed in X ",13)
            end if 
            this%xmethod = 3 
        case ("cheb")
            allocate(this%xcheb)
            ierr = this % xcheb%init( this%xsz(1) )
            if (ierr .ne. 0) then
                call GracefulExit("Initializing Chebyshev Collocation (CHEB) failed in X ",14)
            end if 
            this%xmethod = 4 
        case default 
            call GracefulExit("Invalid method selected in x direction ",01)
        end select
       
        ! Y direction
        select case (method_y)
        case ("cd10")
            allocate(this%ycd10)
            ierr = this % ycd10%init( this%ysz(2), dy, periodic_y, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cd10 failed in Y ",11)
            end if 
            this%ymethod = 1 
        case ("cd06")
            allocate(this%ycd06)
            ierr = this % ycd06%init( this%ysz(2), dy, periodic_y, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cd06 failed in Y",12)
            end if 
            this%ymethod = 2 
        case ("four")
            allocate(this%yfour)
            ierr = this % yfour%init( this%ysz(2), "y", this%ysz(1), this%ysz(3), dy)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing Fourier Collocation (FOUR) failed in Y ",13)
            end if 
            this%ymethod = 3 
        case ("cheb")
            allocate(this%ycheb)
            ierr = this % ycheb%init( this%ysz(2) )
            if (ierr .ne. 0) then
                call GracefulExit("Initializing Chebyshev Collocation (CHEB) failed in Y",14)
            end if 
            this%ymethod = 4 
        case default 
            call GracefulExit("Invalid method selected in y direction",01)
        end select


        ! Z direction
        select case (method_z)
        case ("cd10")
            allocate(this%zcd10)
            ierr = this % zcd10%init( this%zsz(3), dz, periodic_z, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cd10 failed in Z ",11)
            end if 
            this%zmethod = 1 
        case ("cd06")
            allocate(this%zcd06)
            ierr = this % zcd06%init( this%zsz(3), dz, periodic_z, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cd06 failed in Z",12)
            end if 
            this%zmethod = 2 
        case ("four")
            allocate(this%zfour)
            ierr = this % zfour%init( this%zsz(3), "z", this%zsz(1), this%zsz(2), dz)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing Fourier Collocation (FOUR) failed in Z ",13)
            end if 
            this%zmethod = 3 
        case ("cheb")
            allocate(this%zcheb)
            ierr = this % zcheb%init( this%zsz(3) )
            if (ierr .ne. 0) then
                call GracefulExit("Initializing Chebyshev Collocation (CHEB) failed in Z",14)
            end if 
            this%zmethod = 4 
        case default 
            call GracefulExit("Invalid method selected in z direction",01)
        end select

    end subroutine


    subroutine destroy(this)
        class(derivatives), intent(inout) :: this

        select case (this%xmethod) 
        case (1)
            call this%xcd10%destroy
            deallocate(this%xcd10)
        case (2)
            call this%xcd06%destroy
            deallocate(this%xcd06)
        case (3)
            call this%xfour%destroy
            deallocate(this%xfour)
        case (4)
            call this%xcheb%destroy
            deallocate(this%xcheb)
        end select 
        
        select case (this%ymethod) 
        case (1)
            call this%ycd10%destroy
            deallocate(this%ycd10)
        case (2)
            call this%ycd06%destroy
            deallocate(this%ycd06)
        case (3)
            call this%yfour%destroy
            deallocate(this%yfour)
        case (4)
            call this%ycheb%destroy
            deallocate(this%ycheb)
        end select 

        select case (this%zmethod) 
        case (1)
            call this%zcd10%destroy
            deallocate(this%zcd10)
        case (2)
            call this%zcd06%destroy
            deallocate(this%zcd06)
        case (3)
            call this%zfour%destroy
            deallocate(this%zfour)
        case (4)
            call this%zcheb%destroy
            deallocate(this%zcheb)
        end select

        this%xmetric = .false. 
        this%ymetric = .false. 
        this%zmetric = .false. 
        this%curvilinear = .false. 

        if(allocated(this%dxidx))      deallocate(this%dxidx)
        if(allocated(this%dxidx_sq))   deallocate(this%dxidx_sq)
        if(allocated(this%d2xidx2))    deallocate(this%d2xidx2)
        if(allocated(this%dfdxbuf))    deallocate(this%dfdxbuf)
        if(allocated(this%detady))     deallocate(this%detady)
        if(allocated(this%detady_sq))  deallocate(this%detady_sq)
        if(allocated(this%d2etady2))   deallocate(this%d2etady2)
        if(allocated(this%dfdybuf))    deallocate(this%dfdybuf)
        if(allocated(this%dzetadz))    deallocate(this%dzetadz)
        if(allocated(this%dzetadz_sq)) deallocate(this%dzetadz_sq)
        if(allocated(this%d2zetadz2))  deallocate(this%d2zetadz2)
        if(allocated(this%dfdzbuf))    deallocate(this%dfdzbuf)

    end subroutine

    subroutine ddx(this,f,dfdx,bc1,bcn)
        class(derivatives), intent(in) :: this
        real(rkind), intent(in), dimension(this%xsz(1),this%xsz(2),this%xsz(3)) :: f
        real(rkind), intent(out),dimension(this%xsz(1),this%xsz(2),this%xsz(3)) :: dfdx
        integer, optional, intent(in) :: bc1, bcn

        select case (this%xmethod)
        case (1)
            if (present(bc1) .AND. present(bcn)) then
                call this%xcd10 % dd1(f,dfdx,this%xsz(2),this%xsz(3),bc1,bcn)
            else
                call this%xcd10 % dd1(f,dfdx,this%xsz(2),this%xsz(3))
            end if
        case (2)
            call this%xcd06 % dd1(f,dfdx,this%xsz(2),this%xsz(3))
        case (3)
            call this%xfour % dd1(f,dfdx)
        case (4)
            call GracefulExit("Chebychev is incomplete right now",21)
        end select 

        if(this%xmetric) then
            if(this%curvilinear) then
                call GracefulExit("Code INCOMPLETE! Curvilinear Formulation is &
                & not supported!",06)
            else
                dfdx = dfdx * this%dxidx
            endif
        endif

    end subroutine 

    subroutine ddy(this,f,dfdx,bc1,bcn)
        class(derivatives), intent(in) :: this
        real(rkind), intent(in), dimension(this%ysz(1),this%ysz(2),this%ysz(3)) :: f
        real(rkind), intent(out),dimension(this%ysz(1),this%ysz(2),this%ysz(3)) :: dfdx
        integer, optional, intent(in) :: bc1, bcn

        select case (this%ymethod)
        case (1)
            if (present(bc1) .AND. present(bcn)) then
                call this%ycd10 % dd2(f,dfdx,this%ysz(1),this%ysz(3),bc1,bcn)
            else
                call this%ycd10 % dd2(f,dfdx,this%ysz(1),this%ysz(3))
            end if
        case (2)
            call this%ycd06 % dd2(f,dfdx,this%ysz(1),this%ysz(3))
        case (3)
            call this%yfour % dd2(f,dfdx)
        case (4)
            call GracefulExit("Chebychev is incomplete right now",21)
        end select 

        if(this%ymetric) then
            if(this%curvilinear) then
                call GracefulExit("Code INCOMPLETE! Curvilinear Formulation is &
                & not supported!",06)
            else
                dfdx = dfdx * this%detady
            endif
        endif

    end subroutine 
    
    subroutine ddz(this,f,dfdx,bc1,bcn)
        class(derivatives), intent(in) :: this
        real(rkind), intent(in), dimension(this%zsz(1),this%zsz(2),this%zsz(3)) :: f
        real(rkind), intent(out),dimension(this%zsz(1),this%zsz(2),this%zsz(3)) :: dfdx
        integer, optional, intent(in) :: bc1, bcn

        select case (this%zmethod)
        case (1)
            if (present(bc1) .AND. present(bcn)) then
                call this%zcd10 % dd3(f,dfdx,this%zsz(1),this%zsz(2),bc1,bcn)
            else
                call this%zcd10 % dd3(f,dfdx,this%zsz(1),this%zsz(2))
            end if
        case (2)
            call this%zcd06 % dd3(f,dfdx,this%zsz(1),this%zsz(2))
        case (3)
            call this%zfour % dd3(f,dfdx)
        case (4)
            call GracefulExit("Chebychev is incomplete right now",21)
        end select 

        if(this%zmetric) then
            if(this%curvilinear) then
                call GracefulExit("Code INCOMPLETE! Curvilinear Formulation is &
                & not supported!",06)
            else
                dfdx = dfdx * this%dzetadz
            endif
        endif

    end subroutine 

    subroutine d2dx2(this,f,d2fdx2,bc1,bcn)
        class(derivatives), intent(inout) :: this
        real(rkind), intent(in), dimension(this%xsz(1),this%xsz(2),this%xsz(3)) :: f
        real(rkind), intent(out),dimension(this%xsz(1),this%xsz(2),this%xsz(3)) :: d2fdx2
        integer, optional, intent(in) :: bc1, bcn

        select case (this%xmethod)
        case (1)
            call this%xcd10 % d2d1(f,d2fdx2,this%xsz(2),this%xsz(3),bc1,bcn)
        case (2)
            call GracefulExit("CD06 is incomplete right now",21)
        case (3)
            call this%xfour % d2d1(f,d2fdx2)
        case (4)
            call GracefulExit("Chebychev is incomplete right now",21)
        end select 

        if(this%xmetric) then
          if(this%curvilinear) then
              call GracefulExit("Code INCOMPLETE! Curvilinear Formulation is &
              & not supported!",06)
          else
              select case (this%xmethod)
              case (1)
                  call this%xcd10 % dd1(f,this%dfdxbuf,this%xsz(2),this%xsz(3),bc1,bcn)
              case (2)
                  call this%xcd06 % dd1(f,this%dfdxbuf,this%xsz(2),this%xsz(3))
              case (3)
                  call this%xfour % dd1(f,this%dfdxbuf)
              case (4)
                  call GracefulExit("Chebychev is incomplete right now",21)
              end select 
              d2fdx2 = d2fdx2 * this%dxidx_sq + this%dfdxbuf * this%d2xidx2
          endif
        endif

    end subroutine 
    
    subroutine d2dy2(this,f,d2fdx2,bc1,bcn)
        class(derivatives), intent(inout) :: this
        real(rkind), intent(in), dimension(this%ysz(1),this%ysz(2),this%ysz(3)) :: f
        real(rkind), intent(out),dimension(this%ysz(1),this%ysz(2),this%ysz(3)) :: d2fdx2
        integer, optional, intent(in) :: bc1, bcn

        select case (this%ymethod)
        case (1)
            call this%ycd10 % d2d2(f,d2fdx2,this%ysz(1),this%ysz(3),bc1,bcn)
        case (2)
            call GracefulExit("CD06 is incomplete right now",21)
        case (3)
            call this%yfour % d2d2(f,d2fdx2)
        case (4)
            call GracefulExit("Chebychev is incomplete right now",21)
        end select 

        if(this%ymetric) then
          if(this%curvilinear) then
              call GracefulExit("Code INCOMPLETE! Curvilinear Formulation is &
              & not supported!",06)
          else
              select case (this%ymethod)
              case (1)
                  call this%ycd10 % dd2(f,this%dfdybuf,this%ysz(1),this%ysz(3),bc1,bcn)
              case (2)
                  call this%ycd06 % dd2(f,this%dfdybuf,this%ysz(1),this%ysz(3))
              case (3)
                  call this%yfour % dd2(f,this%dfdybuf)
              case (4)
                  call GracefulExit("Chebychev is incomplete right now",21)
              end select 
              d2fdx2 = d2fdx2 * this%detady_sq + this%dfdybuf * this%d2etady2
          endif
        endif

    end subroutine

    subroutine d2dz2(this,f,d2fdx2,bc1,bcn)
        class(derivatives), intent(inout) :: this
        real(rkind), intent(in), dimension(this%zsz(1),this%zsz(2),this%zsz(3)) :: f
        real(rkind), intent(out),dimension(this%zsz(1),this%zsz(2),this%zsz(3)) :: d2fdx2
        integer, optional, intent(in) :: bc1, bcn

        select case (this%zmethod)
        case (1)
            call this%zcd10 % d2d3(f,d2fdx2,this%zsz(1),this%zsz(2),bc1,bcn)
        case (2)
            call GracefulExit("CD06 is incomplete right now",21)
        case (3)
            call this%zfour % d2d3(f,d2fdx2)
        case (4)
            call GracefulExit("Chebychev is incomplete right now",21)
        end select 

        if(this%zmetric) then
          if(this%curvilinear) then
              call GracefulExit("Code INCOMPLETE! Curvilinear Formulation is &
              & not supported!",06)
          else
              select case (this%zmethod)
              case (1)
                  call this%zcd10 % dd3(f,this%dfdzbuf,this%zsz(1),this%zsz(2),bc1,bcn)
              case (2)
                  call this%zcd06 % dd3(f,this%dfdzbuf,this%zsz(1),this%zsz(2))
              case (3)
                  call this%zfour % dd3(f,this%dfdzbuf)
              case (4)
                  call GracefulExit("Chebychev is incomplete right now",21)
              end select 
              d2fdx2 = d2fdx2 * this%dzetadz_sq + this%dfdzbuf * this%d2zetadz2
          endif
        endif

    end subroutine 

    subroutine init_metric_1d(this, xstretch, xuniform, flag, params, dxudxs, dxudxs_sq, d2xudxs2)
        class(derivatives),         intent(in) :: this
        integer,                    intent(in) :: flag
        real(rkind), dimension(5), intent(in)  :: params
        real(rkind), dimension(this%ysz(1),this%ysz(2),this%ysz(3)), intent(in)  :: xstretch, xuniform
        real(rkind), dimension(this%ysz(1),this%ysz(2),this%ysz(3)), intent(out) :: dxudxs, dxudxs_sq, d2xudxs2

        real(rkind) :: xfocus, xfocus_adj, xtau, hh, xstart, num, den, BB, xstretch_loc
        real(rkind) :: xbyxfocm1, xuniform_adj
        integer     :: j

        if(flag==1) then
           ! concentrate towards the center -- Pletcher, Tannehill, Anderson
           ! (Section 5.6, Transformation 3, pg. 332) 
           xfocus = params(1);  xtau   = params(2);  xstart = params(3); hh = params(4)
           print '(5(e19.12,1x))', params(:)
           print '(5(e19.12,1x))', xfocus, xtau, xstart, hh
           xfocus_adj = xfocus - xstart
           num = one + (xfocus_adj/hh) * (exp( xtau) - one)
           den = one + (xfocus_adj/hh) * (exp(-xtau) - one)
           BB = half/xtau*log(num/den)
           do j = 1, this%ysz(2)
               ! adjust for starting point
               xuniform_adj = (xuniform(1,j,1) - xstart) !/ hh

               ! stretched location
               num = sinh(xtau*BB)
               xstretch_loc = xfocus_adj * (one + sinh(xtau * (xuniform_adj/hh-BB))/num) + xstart
               xbyxfocm1 = (xstretch_loc-xstart)/xfocus_adj - one

               ! metric for first derivative
               dxudxs(:,j,:) = num * hh / (xtau * xfocus_adj * sqrt(one + (xbyxfocm1 * num)**2))

               ! square of the metric for first derivative
               dxudxs_sq(:,j,:) = dxudxs(:,j,:)**2

               ! metric for second derivative
               d2xudxs2(:,j,:) = -hh * num**3 * xbyxfocm1 / (xtau*xfocus_adj*xfocus_adj * (one + (num*xbyxfocm1)**2)**1.5d0)

               print '(i5,1x,4(e19.12,1x))', j, xstretch(1,j,1), xuniform(1,j,1), dxudxs(1,j,1), d2xudxs2(1,j,1)
               ! compare with xstretch specified in meshgen
               if(abs(xstretch(1,j,1)-xstretch_loc) > 1.0d-12) then
                   print '(i5,1x,2(e19.12,1x))', j, xstretch(1,j,1), xstretch_loc
                   call GracefulExit("flag = 1; metric is not consistent with meshgen. Check details.", 21)
               endif
           enddo
        elseif(flag==2) then
           ! concentrate towards the two ends
            call GracefulExit("flag = 2 (concentrate towards two ends) is incomplete right now",21)
        elseif(flag==3) then
           ! concentrate at arbitrary point
            call GracefulExit("flag = 3 (concentrate at arbitrary point) is incomplete right now",21)
        elseif(flag==10) then
           ! finite-difference evaluation of metrics (reduces order of accuracy)
            call GracefulExit("flag = 4 (finite-difference evaluation of metrics) is incomplete right now",21)
        endif

    end subroutine 


end module
