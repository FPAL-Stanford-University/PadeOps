module filtersWrapper
    use kind_parameters, only: rkind
    use cf90stuff, only: cf90
    use exits,    only: gracefulExit, message
    use decomp_2d, only: decomp_info, nrank

    implicit none
    private
    public :: filters 
   
    type filters
        private
        
        integer :: xmethod, ymethod, zmethod
        
        type(cf90), allocatable :: xcf90, ycf90, zcf90

        integer, dimension(3)          :: xsz, ysz, zsz ! Local decomposition sizes
        
        logical                        :: initialized = .false. 

        contains
            procedure, private :: init_parallel
            procedure, private :: init_serial
            procedure, private :: init_procedures
            generic :: init => init_parallel, init_serial
            procedure :: destroy
            procedure :: filterx
            procedure :: filtery
            procedure :: filterz 

    end type  

contains 

    subroutine init_procedures(this, nx, ny, nz, &
                                     methodx, methody, methodz, &
                                     periodicx, periodicy, periodicz) 

        class( filters ) , intent(inout) :: this
        integer          , intent(in)    :: nx, ny, nz
        character(len=*) , intent(in)    :: methodx, methody, methodz 
        logical          , intent(in)    :: periodicx, periodicy, periodicz
        
        integer :: ierr

        select case (methodx)
        case ("cf90")
            allocate (this%xcf90)
            ierr = this%xcf90%init( nx, periodicx)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cf90 failed in X ",51)
            end if
            this%xmethod = 1 
        case default
            call GracefulExit("Incorrect method select in direction X", 52)
        end select 
        
        select case (methody)
        case ("cf90")
            allocate (this%ycf90)
            ierr = this%ycf90%init( ny, periodicy)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cf90 failed in Y ",51)
            end if
            this%ymethod = 1 
        case default
            call GracefulExit("Incorrect method select in direction Y", 52)
        end select 

        select case (methodz)
        case ("cf90")
            allocate (this%zcf90)
            ierr = this%zcf90%init( nz, periodicz)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cf90 failed in Z ",51)
            end if
            this%zmethod = 1 
        case default
            call GracefulExit("Incorrect method select in direction Z", 52)
        end select 

    end subroutine 

    subroutine destroy(this)
        class(filters), intent(inout) :: this
        
        select case (this%xmethod)  
        case (1)
            call this%xcf90%destroy
        end select

        select case (this%ymethod)  
        case (1)
            call this%ycf90%destroy
        end select
        
        select case (this%zmethod)  
        case (1)
            call this%zcf90%destroy
        end select

        this%initialized = .false. 
    end subroutine


    subroutine filterx(this, f, ff)
        class (filters), intent(in) :: this
        real(rkind), dimension(this%xsz(1),this%xsz(2), this%xsz(3)), intent(in)  :: f
        real(rkind), dimension(this%xsz(1),this%xsz(2), this%xsz(3)), intent(out) :: ff

        select case (this%xmethod)
        case (1)
            call this%xcf90%filter1( f, ff, this%xsz(2), this%xsz(3))
        end select

    end subroutine

    subroutine filtery(this, f, ff)
        class (filters), intent(in) :: this
        real(rkind), dimension(this%ysz(1),this%ysz(2), this%ysz(3)), intent(in)  :: f
        real(rkind), dimension(this%ysz(1),this%ysz(2), this%ysz(3)), intent(out) :: ff

        select case (this%ymethod)
        case (1)
            call this%ycf90%filter2( f, ff, this%ysz(1), this%ysz(3))
        end select

    end subroutine

    subroutine filterz(this, f, ff)
        class (filters), intent(in) :: this
        real(rkind), dimension(this%zsz(1),this%zsz(2), this%zsz(3)), intent(in)  :: f
        real(rkind), dimension(this%zsz(1),this%zsz(2), this%zsz(3)), intent(out) :: ff

        select case (this%zmethod)
        case (1)
            call this%zcf90%filter3( f, ff, this%zsz(1), this%zsz(2))
        end select

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! INTERNAL SUBROUTINES !!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine init_parallel(this,                              gp  , &
                                     methodx  , methody  , methodz  , &
                                     periodicx, periodicy, periodicz  ) 
        class( filters )   , intent(inout) :: this
        class( decomp_info), intent(in)  :: gp 
        character(len=*)   , intent(in)    :: methodx, methody, methodz 
        logical            , intent(in)    :: periodicx, periodicy, periodicz
        

        if (this%initialized) then
            call message("WARNING: Reinitializing the FILTER class!")
            call this%destroy
        end if  
       
        this%xsz = gp%xsz
        this%ysz = gp%ysz
        this%zsz = gp%zsz

        call this%init_procedures(  this%xsz(1),  this%ysz(2),  this%zsz(3), &
                                     methodx, methody, methodz, &
                                     periodicx, periodicy, periodicz) 

        this%initialized = .true. 
    end subroutine

    subroutine init_serial(this,          nx  ,      ny  ,       nz , &
                                     methodx  , methody  , methodz  , &
                                     periodicx, periodicy, periodicz  )

        class( filters ) , intent(inout) :: this
        integer          , intent(in)    :: nx, ny, nz
        character(len=*) , intent(in)    :: methodx, methody, methodz 
        logical          , intent(in)    :: periodicx, periodicy, periodicz

        if (this%initialized) then
            call message("WARNING: Reinitializing the FILTER class!")
            call this%destroy
        end if  
       
        this%xsz = [nx, ny, nz]
        this%ysz = this%xsz
        this%zsz = this%xsz

        call this%init_procedures(        nx,      ny,      nz, &
                                     methodx, methody, methodz, &
                                     periodicx, periodicy, periodicz) 
    
        this%initialized = .true. 
    end subroutine

end module 
