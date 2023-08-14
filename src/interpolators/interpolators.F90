module InterpolatorsMod
    use kind_parameters, only: rkind, clen
        
    !TODO change names to be interpolators (get rid of all d's)
    use cd10stuff, only: cd10 !-> ci10
    use ci06stuff, only: ci06 !-> ci06
    use ei02stuff, only: ei02 !explicit interpolator, 2nd order
    use d04stuff,  only: d04 !-> ei04
    use ei06stuff, only: ei06

    use exits,    only: gracefulExit, message
    use decomp_2d, only: decomp_info, nrank

    !For purposes of interpolation, (i,j,k) = "Nodes" and (i+1/2,j+1/2,k+1/2) = "Faces"
    !   N2F routine goes from Nodes to Faces (forward )
    !   F2N routine goes from Faces to Nodes (backward)

    !With periodic boundary condtiions, the number of nodes and faces is equal
    !With non-periodic BCs, the last index of the "face" arrays will be garbage

    implicit none
    private 
    public :: interpolators

    type interpolators
        private
       
        
        ! Available Methods: Specific to each direction
        ! 1: "CD10", 2: "CD06", 3: "FOUR", 4: "CHEB", 5: "D02", 6: "D04", 7: "D06"
        integer                :: xmethod, ymethod, zmethod 

        !TODO: change names to be interpolators, list out different ones
        type(cd10), allocatable :: xcd10, ycd10, zcd10 
        type(ci06), allocatable :: xci06, yci06, zci06 
        type(ei02),  allocatable :: xei02, yei02, zei02 
        type(d04),  allocatable :: xd04, yd04, zd04 
        type(ei06),  allocatable :: xei06, yei06, zei06 
    
        logical                        :: curvilinear = .false. 
        logical                        :: xmetric = .false.
        logical                        :: ymetric = .false.
        logical                        :: zmetric = .false.

        integer                        :: nxg, nyg, nzg ! Global sizes
        integer, dimension(3)          :: xsz, ysz, zsz ! Local decomposition sizes
        
        logical                        :: initialized = .false. 
        contains
            procedure, private :: init_parallel
            procedure, private :: init_serial
            procedure, private :: init_procedures
            procedure, private :: init_curvilinear
            generic :: init => init_parallel, init_serial
            procedure :: destroy
            procedure :: iN2Fx
            procedure :: iN2Fy
            procedure :: iN2Fz
            procedure :: iF2Nx
            procedure :: iF2Ny
            procedure :: iF2Nz
            procedure :: getMethodx
            procedure :: getMethody
            procedure :: getMethodz
            procedure :: set_xsz
            procedure :: set_ysz
            procedure :: set_zsz
    end type


contains

    subroutine set_xsz(this,xsz)
        class(interpolators), intent(inout) :: this
        integer, dimension(3), intent(in) :: xsz

        if (this%xsz(1) /= xsz(1)) then
            call GracefulExit("Cannot change xsz(1) in set_xsz",453)
        end if
        this%xsz = xsz
    end subroutine

    subroutine set_ysz(this,ysz)
        class(interpolators), intent(inout) :: this
        integer, dimension(3), intent(in) :: ysz

        if (this%ysz(2) /= ysz(2)) then
            call GracefulExit("Cannot change ysz(2) in set_ysz",453)
        end if
        this%ysz = ysz
    end subroutine

    subroutine set_zsz(this,zsz)
        class(interpolators), intent(inout) :: this
        integer, dimension(3), intent(in) :: zsz

        if (this%zsz(3) /= zsz(3)) then
            call GracefulExit("Cannot change zsz(3) in set_zsz",453)
        end if
        this%zsz = zsz
    end subroutine

    function getMethodx(this) result(m)
        class(interpolators), intent(in) :: this
        character(len=clen) :: m 
        select case (this%xmethod)
        case (1)
            m = "cd10"
        case (2)
            m = "ci06"
        case (3)
            call GracefulExit("Case doesn't exist",453)
        case (4)
            call GracefulExit("Case doesn't exist",453)
        case (5)
            m = "ei02"
        case (6)
            m = "d04"
        case (7)
            m = "ei06"
        end select 

    end function

    function getMethody(this) result(m)
        class(interpolators), intent(in) :: this
        character(len=clen) :: m 
        select case (this%ymethod)
        case (1)
            m = "cd10"
        case (2)
            m = "ci06"
        case (3)
            call GracefulExit("Case doesn't exist",453)
        case (4)
            call GracefulExit("Case doesn't exist",453)
        case (5)
            m = "ei02"
        case (6)
            m = "d04"
        case (7)
            m = "ei06"
        end select 

    end function

    function getMethodz(this) result(m)
        class(interpolators), intent(in) :: this
        character(len=clen) :: m 
        select case (this%zmethod)
        case (1)
            m = "cd10"
        case (2)
            m = "ci06"
        case (3)
            call GracefulExit("Case doesn't exist",453)
        case (4)
            call GracefulExit("Case doesn't exist",453)
        case (5)
            m = "ei02"
        case (6)
            m = "d04"
        case (7)
            m = "ei06"
        end select 

    end function


    subroutine init_serial(this,nx,         ny,         nz, &
                                dx,         dy,         dz, &
                        periodic_x, periodic_y, periodic_z, &
                          method_x,   method_y,   method_z, &
                           xmetric,    ymetric,    zmetric, &
                       curvilinear)
        
        class(interpolators), intent(inout)          :: this
        integer, intent(in)                        :: nx, ny, nz
        real(rkind)       , intent(in)             :: dx, dy, dz
        logical           , intent(in)             :: periodic_x, periodic_y, periodic_z
        character(len=*)  , intent(in)             :: method_x
        character(len=*)  , intent(in)             :: method_y
        character(len=*)  , intent(in)             :: method_z
        logical           , intent(in),   optional :: xmetric
        logical           , intent(in),   optional :: ymetric
        logical           , intent(in),   optional :: zmetric
        logical           , intent(in),   optional :: curvilinear
        


       
        if (this%initialized) then
            call message("WARNING: Reinitializing the INTERPOLATORS class!")
            call this%destroy
        end if  
       
        this%nxg = nx
        this%nyg = ny
        this%nzg = nz
        this%xsz = [nx, ny, nz]
        this%ysz = this%xsz
        this%zsz = this%xsz

        call this%init_procedures (dx, dy, dz, periodic_x, periodic_y, periodic_z,&
                                            method_x,method_y,method_z)

        if (present(xmetric)) then
            if (present(ymetric)) then
                if (present(zmetric)) then
                    if (present(curvilinear)) then
                        call this%init_curvilinear(xMetric, yMetric, zMetric, curvilinear)
                    else
                        call this%init_curvilinear(xMetric, yMetric, zMetric)
                    end if 
                else
                    call this%init_curvilinear(xMetric, yMetric)
                end if 
            else
                call this%init_curvilinear(xMetric)
            end if
        end if 

    end subroutine

    subroutine init_parallel(this,                      gp, &
                                dx,         dy,         dz, &
                        periodic_x, periodic_y, periodic_z, &
                          method_x,   method_y,   method_z, &
                           xmetric,    ymetric,    zmetric, &
                       curvilinear)
        
        class(interpolators), intent(inout)          :: this
        class(decomp_info), intent(in)             :: gp 
        real(rkind)       , intent(in)             :: dx, dy, dz
        logical           , intent(in)             :: periodic_x, periodic_y, periodic_z
        character(len=*)  , intent(in)             :: method_x
        character(len=*)  , intent(in)             :: method_y
        character(len=*)  , intent(in)             :: method_z
        logical           , intent(in),   optional :: xmetric
        logical           , intent(in),   optional :: ymetric
        logical           , intent(in),   optional :: zmetric
        logical           , intent(in),   optional :: curvilinear
        

       
        if (this%initialized) then
            call message("WARNING: Reinitializing the INTERPOLATORS class!")
            call this%destroy
        end if  
       
        this%nxg = gp%xsz(1)
        this%nyg = gp%ysz(2)
        this%nzg = gp%zsz(3)
        this%xsz = gp%xsz
        this%ysz = gp%ysz
        this%zsz = gp%zsz


        call this%init_procedures (dx, dy, dz, periodic_x, periodic_y, periodic_z,&
                                            method_x,method_y,method_z)

        if (present(xmetric)) then
            if (present(ymetric)) then
                if (present(zmetric)) then
                    if (present(curvilinear)) then
                        call this%init_curvilinear(xMetric, yMetric, zMetric, curvilinear)
                    else
                        call this%init_curvilinear(xMetric, yMetric, zMetric)
                    end if 
                else
                    call this%init_curvilinear(xMetric, yMetric)
                end if 
            else
                call this%init_curvilinear(xMetric)
            end if
        end if 

    end subroutine


    subroutine init_curvilinear(this, xMetric, yMetric, zMetric, Curvilinear)
        class(interpolators), intent(inout)        :: this
        logical           , intent(in)           :: xmetric
        logical           , intent(in), optional :: ymetric
        logical           , intent(in), optional :: zmetric
        logical           , intent(in), optional :: curvilinear
        

        ! Metrics
            if (xMetric) then
                call GracefulExit("Code INCOMPLETE! Metric terms are not &
                 &   supported in x direction!",04)
                 this%xMetric = .true. 
            end if

            if (yMetric) then
                call GracefulExit("Code INCOMPLETE! Metric terms are not &
                 &   supported in y direction!",05)
                 this%yMetric = .true. 
            end if
        
            if (zMetric) then
                call GracefulExit("Code INCOMPLETE! Metric terms are not &
                 &   supported in z direction!",06)
                 this%zMetric = .true. 
            end if
        
            if (Curvilinear) then
                call GracefulExit("Code INCOMPLETE! Curvilinear Formulation is &
                & not supported!",06)
                this%Curvilinear = .true. 
            end if

    end subroutine

    subroutine init_procedures(this,dx, dy, dz, periodic_x, periodic_y, periodic_z,&
                                    method_x,method_y,method_z)
                                    
        class(interpolators), intent(inout)          :: this
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
        case ("ci06")
            allocate(this%xci06)
            ierr = this % xci06%init( this%xsz(1), dx, periodic_x, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing ci06 failed in X ",12)
            end if 
            this%xmethod = 2 
        case ("ei02")
            allocate(this%xei02)
            ierr = this % xei02%init( this%xsz(1), dx, periodic_x, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing ei02 failed in X ",12)
            end if 
            this%xmethod = 5 
        case ("d04")
            allocate(this%xd04)
            ierr = this % xd04%init( this%xsz(1), dx, periodic_x, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing d04 failed in X ",12)
            end if 
            this%xmethod = 6 
        case ("ei06")
            allocate(this%xei06)
            ierr = this % xei06%init( this%xsz(1), dx, periodic_x, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing ei06 failed in X ",12)
            end if 
            this%xmethod = 7 
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
        case ("ci06")
            allocate(this%yci06)
            ierr = this % yci06%init( this%ysz(2), dy, periodic_y, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing ci06 failed in Y",12)
            end if 
            this%ymethod = 2 
        case ("ei02")
            allocate(this%yei02)
            ierr = this % yei02%init( this%ysz(2), dy, periodic_y, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing ei02 failed in Y ",12)
            end if 
            this%ymethod = 5 
        case ("d04")
            allocate(this%yd04)
            ierr = this % yd04%init( this%ysz(2), dy, periodic_y, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing d04 failed in Y ",12)
            end if 
            this%ymethod = 6 
        case ("ei06")
            allocate(this%yei06)
            ierr = this % yei06%init( this%ysz(2), dy, periodic_y, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing ei06 failed in Y ",12)
            end if 
            this%ymethod = 7 
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
        case ("ci06")
            allocate(this%zci06)
            ierr = this % zci06%init( this%zsz(3), dz, periodic_z, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing ci06 failed in Z",12)
            end if 
            this%zmethod = 2 
        case ("ei02")
            allocate(this%zei02)
            ierr = this % zei02%init( this%zsz(3), dz, periodic_z, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing ei02 failed in Z ",12)
            end if 
            this%zmethod = 5 
        case ("d04")
            allocate(this%zd04)
            ierr = this % zd04%init( this%zsz(3), dz, periodic_z, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing d04 failed in Z ",12)
            end if 
            this%zmethod = 6 
        case ("ei06")
            allocate(this%zei06)
            ierr = this % zei06%init( this%zsz(3), dz, periodic_z, 0, 0)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing ei06 failed in Z ",12)
            end if 
            this%zmethod = 7 
        case default 
            call GracefulExit("Invalid method selected in z direction",01)
        end select

    end subroutine


    subroutine destroy(this)
        class(interpolators), intent(inout) :: this

        select case (this%xmethod) 
        case (1)
            call this%xcd10%destroy
            deallocate(this%xcd10)
        case (2)
            call this%xci06%destroy
            deallocate(this%xci06)
        case (3)
            call GracefulExit("Case doesn't exist",453)
        case (4)
            call GracefulExit("Case doesn't exist",453)
        case (5)
            call this%xei02%destroy
            deallocate(this%xei02)
        case (6)
            call this%xd04%destroy
            deallocate(this%xd04)
        case (7)
            call this%xei06%destroy
            deallocate(this%xei06)
        end select 
        
        select case (this%ymethod) 
        case (1)
            call this%ycd10%destroy
            deallocate(this%ycd10)
        case (2)
            call this%yci06%destroy
            deallocate(this%yci06)
        case (3)
            call GracefulExit("Case doesn't exist",453)
        case (4)
            call GracefulExit("Case doesn't exist",453)
        case (5)
            call this%yei02%destroy
            deallocate(this%yei02)
        case (6)
            call this%yd04%destroy
            deallocate(this%yd04)
        case (7)
            call this%yei06%destroy
            deallocate(this%yei06)
        end select 

        select case (this%zmethod) 
        case (1)
            call this%zcd10%destroy
            deallocate(this%zcd10)
        case (2)
            call this%zci06%destroy
            deallocate(this%zci06)
        case (3)
            call GracefulExit("Case doesn't exist",453)
        case (4)
            call GracefulExit("Case doesn't exist",453)
        case (5)
            call this%zei02%destroy
            deallocate(this%zei02)
        case (6)
            call this%zd04%destroy
            deallocate(this%zd04)
        case (7)
            call this%zei06%destroy
            deallocate(this%zei06)
        end select

        this%xmetric = .false. 
        this%ymetric = .false. 
        this%zmetric = .false. 
        this%curvilinear = .false. 

    end subroutine

    subroutine iN2Fx(this,fN,fF,bc1,bcn)
        class(interpolators), intent(in) :: this
        real(rkind), intent(in), dimension(this%xsz(1),this%xsz(2),this%xsz(3)) :: fN !fNode
        real(rkind), intent(out),dimension(this%xsz(1),this%xsz(2),this%xsz(3)) :: fF !fFace
        integer, optional, intent(in) :: bc1, bcn

      
       select case (this%xmethod)
        case (1)
            if (present(bc1) .AND. present(bcn)) then
                call this%xcd10 % dd1(fN,fF,this%xsz(2),this%xsz(3),bc1,bcn)
            else
                call this%xcd10 % dd1(fN,fF,this%xsz(2),this%xsz(3))
            end if
        case (2)
            if (present(bc1) .AND. present(bcn)) then
            call this%xci06 % iN2F1(fN,fF,this%xsz(2),this%xsz(3),bc1,bcn)
             else
            call this%xci06 % iN2F1(fN,fF,this%xsz(2),this%xsz(3))
            end if
        case (3)
            call GracefulExit("Case doesn't exist",453)
        case (4)
            call GracefulExit("Chebychev is incomplete right now",21)
        case (5)
            if (present(bc1) .AND. present(bcn)) then
                call this%xei02 % iN2F1(fN,fF,this%xsz(2),this%xsz(3),bc1,bcn)
            else
                call this%xei02 % iN2F1(fN,fF,this%xsz(2),this%xsz(3))
            end if
        case (6)
            if (present(bc1) .AND. present(bcn)) then
                call this%xd04 % dd1(fN,fF,this%xsz(2),this%xsz(3),bc1,bcn)
            else
                call this%xd04 % dd1(fN,fF,this%xsz(2),this%xsz(3))
            end if
        case (7)
            if (present(bc1) .AND. present(bcn)) then
                call this%xei06 % iN2F1(fN,fF,this%xsz(2),this%xsz(3),bc1,bcn)
            else
                call this%xei06 % iN2F1(fN,fF,this%xsz(2),this%xsz(3))

            end if
        end select 

    end subroutine 

    subroutine iN2Fy(this,fN,fF,bc1,bcn)
        class(interpolators), intent(in) :: this
        real(rkind), intent(in), dimension(this%ysz(1),this%ysz(2),this%ysz(3)) :: fN !fNode
        real(rkind), intent(out),dimension(this%ysz(1),this%ysz(2),this%ysz(3)) :: fF !fFace
        integer, optional, intent(in) :: bc1, bcn


        select case (this%ymethod)
        case (1)
            if (present(bc1) .AND. present(bcn)) then
                call this%ycd10 % dd2(fN,fF,this%ysz(1),this%ysz(3),bc1,bcn)
            else
                call this%ycd10 % dd2(fN,fF,this%ysz(1),this%ysz(3))
            end if
        case (2)
            if (present(bc1) .AND. present(bcn)) then
            call this%yci06 % iN2F2(fN,fF,this%ysz(1),this%ysz(3),bc1,bcn)
            else
            call this%yci06 % iN2F2(fN,fF,this%ysz(1),this%ysz(3))
            end if
        case (3)
            call GracefulExit("Case doesn't exist",453)
        case (4)
            call GracefulExit("Chebychev is incomplete right now",21)
        case (5)
            if (present(bc1) .AND. present(bcn)) then
                call this%yei02 % iN2F2(fN,fF,this%ysz(1),this%ysz(3),bc1,bcn)
            else
                call this%yei02 % iN2F2(fN,fF,this%ysz(1),this%ysz(3))
            end if
        case (6)
            if (present(bc1) .AND. present(bcn)) then
                call this%yd04 % dd2(fN,fF,this%ysz(1),this%ysz(3),bc1,bcn)
            else
                call this%yd04 % dd2(fN,fF,this%ysz(1),this%ysz(3))
            end if
        case (7)
            if (present(bc1) .AND. present(bcn)) then
                call this%yei06 % iN2F2(fN,fF,this%ysz(1),this%ysz(3),bc1,bcn)

            else
                call this%yei06 % iN2F2(fN,fF,this%ysz(1),this%ysz(3))

            end if
        end select 

    end subroutine 
    
    subroutine iN2Fz(this,fN,fF,bc1,bcn)
        class(interpolators), intent(in) :: this
        real(rkind), intent(in), dimension(this%zsz(1),this%zsz(2),this%zsz(3)) :: fN !fNode
        real(rkind), intent(out),dimension(this%zsz(1),this%zsz(2),this%zsz(3)) :: fF !fFace
        integer, optional, intent(in) :: bc1, bcn

        select case (this%zmethod)
        case (1)
            if (present(bc1) .AND. present(bcn)) then
                call this%zcd10 % dd3(fN,fF,this%zsz(1),this%zsz(2),bc1,bcn)
            else
                call this%zcd10 % dd3(fN,fF,this%zsz(1),this%zsz(2))
            end if
        case (2)
            call this%zci06 % iN2F3(fN,fF,this%zsz(1),this%zsz(2))
        case (3)
            call GracefulExit("Case doesn't exist",453)
        case (4)
            call GracefulExit("Chebychev is incomplete right now",21)
        case (5)
            if (present(bc1) .AND. present(bcn)) then
                call this%zei02 % iN2F3(fN,fF,this%zsz(1),this%zsz(2),bc1,bcn)
            else
                call this%zei02 % iN2F3(fN,fF,this%zsz(1),this%zsz(2))
            end if
        case (6)
            if (present(bc1) .AND. present(bcn)) then
                call this%zd04 % dd3(fN,fF,this%zsz(1),this%zsz(2),bc1,bcn)
            else
                call this%zd04 % dd3(fN,fF,this%zsz(1),this%zsz(2))
            end if
        case (7)
            if (present(bc1) .AND. present(bcn)) then
                call this%zei06 % iN2F3(fN,fF,this%zsz(1),this%zsz(2),bc1,bcn)

            else
                call this%zei06 % iN2F3(fN,fF,this%zsz(1),this%zsz(2))
            end if
        end select 
    end subroutine 

    subroutine iF2Nx(this,fF,fN,bc1,bcn)
        class(interpolators), intent(in) :: this
        real(rkind), intent(in),dimension(this%xsz(1),this%xsz(2),this%xsz(3)) :: fF !fFace
        real(rkind), intent(out), dimension(this%xsz(1),this%xsz(2),this%xsz(3)) :: fN !fNode
        integer, optional, intent(in) :: bc1, bcn

        select case (this%xmethod)
        case (1)
            if (present(bc1) .AND. present(bcn)) then
                call this%xcd10 % dd1(fF,fN,this%xsz(2),this%xsz(3),bc1,bcn)
            else
                call this%xcd10 % dd1(fF,fN,this%xsz(2),this%xsz(3))
            end if
        case (2)
            call this%xci06 % iF2N1(fF,fN,this%xsz(2),this%xsz(3))
        case (3)
            call GracefulExit("Case doesn't exist",453)
        case (4)
            call GracefulExit("Chebychev is incomplete right now",21)
        case (5)
            if (present(bc1) .AND. present(bcn)) then
                call this%xei02 % iF2N1(fF,fN,this%xsz(2),this%xsz(3),bc1,bcn)
            else
                call this%xei02 % iF2N1(fF,fN,this%xsz(2),this%xsz(3))
            end if
        case (6)
            if (present(bc1) .AND. present(bcn)) then
                call this%xd04 % dd1(fF,fN,this%xsz(2),this%xsz(3),bc1,bcn)
            else
                call this%xd04 % dd1(fF,fN,this%xsz(2),this%xsz(3))
            end if
        case (7)
            if (present(bc1) .AND. present(bcn)) then
                call this%xei06 % iF2N1(fF,fN,this%xsz(2),this%xsz(3),bc1,bcn)

            else
                call this%xei06 % iF2N1(fF,fN,this%xsz(2),this%xsz(3))
            end if
        end select 

    end subroutine 

    subroutine iF2Ny(this,fF,fN,bc1,bcn)
        class(interpolators), intent(in) :: this
        real(rkind), intent(in),dimension(this%ysz(1),this%ysz(2),this%ysz(3)) :: fF !fFace
        real(rkind), intent(out), dimension(this%ysz(1),this%ysz(2),this%ysz(3)) :: fN !fNode
        integer, optional, intent(in) :: bc1, bcn

        select case (this%ymethod)
        case (1)
            if (present(bc1) .AND. present(bcn)) then
                call this%ycd10 % dd2(fF,fN,this%ysz(2),this%ysz(3),bc1,bcn)
            else
                call this%ycd10 % dd2(fF,fN,this%ysz(2),this%ysz(3))
            end if
        case (2)
            call this%yci06 % iF2N2(fF,fN,this%ysz(2),this%ysz(3))
        case (3)
            call GracefulExit("Case doesn't exist",453)
        case (4)
            call GracefulExit("Chebychev is incomplete right now",21)
        case (5)
            if (present(bc1) .AND. present(bcn)) then
                call this%yei02 % iF2N2(fF,fN,this%ysz(2),this%ysz(3),bc1,bcn)
            else
                call this%yei02 % iF2N2(fF,fN,this%ysz(2),this%ysz(3))
            end if
        case (6)
            if (present(bc1) .AND. present(bcn)) then
                call this%yd04 % dd2(fF,fN,this%ysz(2),this%ysz(3),bc1,bcn)
            else
                call this%yd04 % dd2(fF,fN,this%ysz(2),this%ysz(3))
            end if
        case (7)
            if (present(bc1) .AND. present(bcn)) then
                call this%yei06 % iF2N2(fF,fN,this%ysz(2),this%ysz(3),bc1,bcn)

            else
                call this%yei06 % iF2N2(fF,fN,this%ysz(2),this%ysz(3))
            end if
        end select 

    end subroutine 

    subroutine iF2Nz(this,fF,fN,bc1,bcn)
        class(interpolators), intent(in) :: this
        real(rkind), intent(in),dimension(this%zsz(1),this%zsz(2),this%zsz(3)) :: fF !fFace
        real(rkind), intent(out), dimension(this%zsz(1),this%zsz(2),this%zsz(3)) :: fN !fNode
        integer, optional, intent(in) :: bc1, bcn

        select case (this%zmethod)
        case (1)
            if (present(bc1) .AND. present(bcn)) then
                call this%zcd10 % dd3(fF,fN,this%zsz(2),this%zsz(3),bc1,bcn)
            else
                call this%zcd10 % dd3(fF,fN,this%zsz(2),this%zsz(3))
            end if
        case (2)
            call this%zci06 % iF2N3(fF,fN,this%zsz(2),this%zsz(3))
        case (3)
            call GracefulExit("Case doesn't exist",453)
        case (4)
            call GracefulExit("Chebychev is incomplete right now",21)
        case (5)
            if (present(bc1) .AND. present(bcn)) then
                call this%zei02 % iF2N3(fF,fN,this%zsz(2),this%zsz(3),bc1,bcn)
            else
                call this%zei02 % iF2N3(fF,fN,this%zsz(2),this%zsz(3))
            end if
        case (6)
            if (present(bc1) .AND. present(bcn)) then
                call this%zd04 % dd3(fF,fN,this%zsz(2),this%zsz(3),bc1,bcn)
            else
                call this%zd04 % dd3(fF,fN,this%zsz(2),this%zsz(3))
            end if
        case (7)
            if (present(bc1) .AND. present(bcn)) then
                call this%zei06 % iF2N3(fF,fN,this%zsz(2),this%zsz(3),bc1,bcn)
            else
                call this%zei06 % iF2N3(fF,fN,this%zsz(2),this%zsz(3))
            end if
        end select 

    end subroutine 


end module
