module FiltersMod
    use kind_parameters, only: rkind, clen
    use cf90stuff,       only: cf90
    use gaussianstuff,   only: gaussian
    use lstsqstuff,      only: lstsq
    use box2stuff,       only: box2
    use exits,           only: gracefulExit, message
    use decomp_2d,       only: decomp_info, nrank
    use MultiBlockTopologyMod, only: multiblocktopol

    implicit none
    private
    public :: filters 
   
    type filters
        private
        
        integer :: xmethod, ymethod, zmethod
        
        type(cf90)    , allocatable, dimension(:) :: xcf90, ycf90, zcf90
        type(gaussian), allocatable, dimension(:) :: xgauf, ygauf, zgauf
        type(lstsq)   , allocatable :: xlsqf, ylsqf, zlsqf
        type(box2)    , allocatable :: xbox2, ybox2, zbox2

        integer, dimension(3)          :: xsz, ysz, zsz ! Local decomposition sizes
        
        class(multiblocktopol), pointer :: mbtopology

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
            procedure :: getMethodx
            procedure :: getMethody
            procedure :: getMethodz

    end type  

contains 

    function getMethodx(this) result(m)
        class(filters), intent(in) :: this
        character(len=clen) :: m 
        select case (this%xmethod)
        case (1)
            m = "cf90"
        case (2)
            m = "gaussian"
        case (3)
            m = "lstsq"
        case (4)
            m = "box2"
        end select 
    end function

    function getMethody(this) result(m)
        class(filters), intent(in) :: this
        character(len=clen) :: m 
        select case (this%ymethod)
        case (1)
            m = "cf90"
        case (2)
            m = "gaussian"
        case (3)
            m = "lstsq"
        case (4)
            m = "box2"
        end select
    end function 
        
    function getMethodz(this) result(m)
        class(filters), intent(in) :: this
        character(len=clen) :: m 
        select case (this%zmethod)
        case (1)
            m = "cf90"
        case (2)
            m = "gaussian"
        case (3)
            m = "lstsq"
        case (4)
            m = "box2"
        end select
    end function 
        
    subroutine init_procedures(this, nx, ny, nz, &
                                     methodx, methody, methodz, &
                                     periodicx, periodicy, periodicz)

        class( filters ) , intent(inout) :: this
        integer          , intent(in)    :: nx, ny, nz
        character(len=*) , intent(in)    :: methodx, methody, methodz 
        logical          , intent(in)    :: periodicx, periodicy, periodicz
        
        integer :: ierr, imb

        select case (methodx)
        
        case ("cf90")
          if(associated(this%mbtopology)) then
            allocate(this%xcf90(this%mbtopology%x_num_blocks))
            do imb = 1, this%mbtopology%x_num_blocks
              ierr = this%xcf90(imb)%init(this%mbtopology%xen(1,imb)-this%mbtopology%xst(1,imb)+1, periodicx)
              if (ierr .ne. 0) then
                  call message("Initializing cf90 block number ", imb)
                  call GracefulExit("Initializing cf90 failed in X ",51)
              end if
            enddo
          else
            allocate (this%xcf90(1))
            ierr = this%xcf90(1)%init( nx, periodicx)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cf90 failed in X ",51)
            end if
          endif
          this%xmethod = 1
        
        case("gaussian") 
          if(associated(this%mbtopology)) then
            allocate(this%xgauf(this%mbtopology%x_num_blocks))
            do imb = 1, this%mbtopology%x_num_blocks
              ierr = this%xgauf(imb)%init(this%mbtopology%xen(1,imb)-this%mbtopology%xst(1,imb)+1, periodicx)
              if (ierr .ne. 0) then
                  call message("Initializing gaussian filter block number ", imb)
                  call GracefulExit("Initializing gaussian filter failed in X ",51)
              end if
            enddo
          else
            allocate (this%xgauf(1))
            ierr = this%xgauf(1)%init( nx, periodicx)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing gaussian filter failed in X ",51)
            end if
          endif
          this%xmethod = 2
        
        case("lstsq") 
            allocate (this%xlsqf)
            ierr = this%xlsqf%init( nx, periodicx)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing least squares filter failed in X ",51)
            end if
            this%xmethod = 3
        
        case("box2") 
            allocate (this%xbox2)
            ierr = this%xbox2%init(nx, periodicx)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing box2 filter failed in X ",51)
            end if
            this%xmethod = 4
        
        case default
            call GracefulExit("Incorrect method select in direction X", 52)
        end select 
        
        
        select case (methody)
        
        case ("cf90")
          if(associated(this%mbtopology)) then
            allocate(this%ycf90(this%mbtopology%y_num_blocks))
            do imb = 1, this%mbtopology%y_num_blocks
              ierr = this%ycf90(imb)%init(this%mbtopology%yen(1,imb)-this%mbtopology%yst(1,imb)+1, periodicy)
              if (ierr .ne. 0) then
                  call message("Initializing cf90 block number ", imb)
                  call GracefulExit("Initializing cf90 failed in Y ",51)
              end if
            enddo
          else
            allocate (this%ycf90(1))
            ierr = this%ycf90(1)%init( ny, periodicy)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cf90 failed in Y ",51)
            end if
          endif
          this%ymethod = 1 
        
        case("gaussian") 
          if(associated(this%mbtopology)) then
            allocate(this%ygauf(this%mbtopology%y_num_blocks))
            do imb = 1, this%mbtopology%y_num_blocks
              ierr = this%ygauf(imb)%init(this%mbtopology%yen(1,imb)-this%mbtopology%yst(1,imb)+1, periodicy)
              if (ierr .ne. 0) then
                  call message("Initializing gaussian filter block number ", imb)
                  call GracefulExit("Initializing gaussian filter failed in Y ",51)
              end if
            enddo
          else
            allocate (this%ygauf(1))
            ierr = this%ygauf(1)%init( ny, periodicy)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing gaussian filter failed in Y ",51)
            end if
          endif
          this%ymethod = 2
        
        case("lstsq") 
            allocate (this%ylsqf)
            ierr = this%ylsqf%init( ny, periodicy)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing least squares filter failed in Y ",51)
            end if
            this%ymethod = 3
        
        case("box2") 
            allocate (this%ybox2)
            ierr = this%ybox2%init(ny, periodicy)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing box2 filter failed in Y ",51)
            end if
            this%ymethod = 4
        
        case default
            call GracefulExit("Incorrect method select in direction Y", 52)
        end select 

        select case (methodz)

        case ("cf90")
          if(associated(this%mbtopology)) then
            allocate(this%zcf90(this%mbtopology%z_num_blocks))
            do imb = 1, this%mbtopology%z_num_blocks
              ierr = this%zcf90(imb)%init(this%mbtopology%zen(1,imb)-this%mbtopology%zst(1,imb)+1, periodicz)
              if (ierr .ne. 0) then
                  call message("Initializing cf90 block number ", imb)
                  call GracefulExit("Initializing cf90 failed in Z ",51)
              end if
            enddo
          else
            allocate (this%zcf90(1))
            ierr = this%zcf90(1)%init( nz, periodicz)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing cf90 failed in Z ",51)
            end if
          endif
          this%zmethod = 1 
        
        case("gaussian") 
          if(associated(this%mbtopology)) then
            allocate(this%zgauf(this%mbtopology%z_num_blocks))
            do imb = 1, this%mbtopology%z_num_blocks
              ierr = this%zgauf(imb)%init(this%mbtopology%zen(1,imb)-this%mbtopology%zst(1,imb)+1, periodicz)
              if (ierr .ne. 0) then
                  call message("Initializing gaussian filter block number ", imb)
                  call GracefulExit("Initializing gaussian filter failed in Z ",51)
              end if
            enddo
          else
            allocate (this%zgauf(1))
            ierr = this%zgauf(1)%init( nz, periodicz)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing gaussian filter failed in Z ",51)
            end if
          endif
          this%zmethod = 2
        
        case("lstsq") 
            allocate (this%zlsqf)
            ierr = this%zlsqf%init( nz, periodicz)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing least squares filter failed in Z ",51)
            end if
            this%zmethod = 3
        
        case("box2") 
            allocate (this%zbox2)
            ierr = this%zbox2%init(nz, periodicz)
            if (ierr .ne. 0) then
                call GracefulExit("Initializing box2 filter failed in Z ",51)
            end if
            this%zmethod = 4
        
        case default
            call GracefulExit("Incorrect method select in direction Z", 52)
        end select 

    end subroutine 

    subroutine destroy(this)
        class(filters), intent(inout) :: this

        integer :: imb
 
        select case (this%xmethod)  
        case (1)
          if(associated(this%mbtopology)) then
            do imb = 1, this%mbtopology%x_num_blocks
              call this%xcf90(imb)%destroy()
            enddo
          else
            call this%xcf90(1)%destroy
          endif
          deallocate(this%xcf90)
        case (2)
          if(associated(this%mbtopology)) then
            do imb = 1, this%mbtopology%x_num_blocks
              call this%xgauf(imb)%destroy()
            enddo
          else
            call this%xgauf(1)%destroy
          endif
          deallocate(this%xgauf)
        case (3)
            call this%xlsqf%destroy
        case (4)
            call this%xbox2%destroy
        end select

        select case (this%ymethod)  
        case (1)
          if(associated(this%mbtopology)) then
            do imb = 1, this%mbtopology%x_num_blocks
              call this%ycf90(imb)%destroy()
            enddo
          else
            call this%ycf90(1)%destroy
          endif
          deallocate(this%ycf90)
        case (2)
          if(associated(this%mbtopology)) then
            do imb = 1, this%mbtopology%x_num_blocks
              call this%ygauf(imb)%destroy()
            enddo
          else
            call this%ygauf(1)%destroy
          endif
          deallocate(this%ygauf)
        case (3)
            call this%ylsqf%destroy
        case (4)
            call this%ybox2%destroy
        end select
        
        select case (this%zmethod)  
        case (1)
          if(associated(this%mbtopology)) then
            do imb = 1, this%mbtopology%x_num_blocks
              call this%zcf90(imb)%destroy()
            enddo
          else
            call this%zcf90(1)%destroy
          endif
          deallocate(this%zcf90)
        case (2)
          if(associated(this%mbtopology)) then
            do imb = 1, this%mbtopology%x_num_blocks
              call this%zgauf(imb)%destroy()
            enddo
          else
            call this%zgauf(1)%destroy
          endif
          deallocate(this%zgauf)
        case (3)
            call this%zlsqf%destroy
        case (4)
            call this%zbox2%destroy
        end select

        if(associated(this%mbtopology)) nullify(this%mbtopology)

        this%initialized = .false. 
    end subroutine


    subroutine filterx(this, f, ff, bc1, bcn)
        class (filters), intent(in) :: this
        real(rkind), dimension(this%xsz(1),this%xsz(2), this%xsz(3)), intent(in)  :: f
        real(rkind), dimension(this%xsz(1),this%xsz(2), this%xsz(3)), intent(out) :: ff
        integer, optional, intent(in) :: bc1, bcn

        integer :: imb, nxst, nyst, nzst, nxen, nyen, nzen

        select case (this%xmethod)
        case (1)
          if(associated(this%mbtopology)) then
            do imb = 1, this%mbtopology%x_num_blocks
              nxst = this%mbtopology%xst(1, imb);   nxen = this%mbtopology%xen(1, imb)
              nyst = this%mbtopology%xst(2, imb);   nyen = this%mbtopology%xen(2, imb)
              nzst = this%mbtopology%xst(3, imb);   nzen = this%mbtopology%xen(3, imb)
              call this%xcf90(imb) % filter1(f(nxst:nxen,nyst:nyen,nzst:nzen), &
                                            ff(nxst:nxen,nyst:nyen,nzst:nzen), &
                                            nyen-nyst+1, nzen-nzst+1, bc1, bcn)
            enddo
          else
            call this%xcf90(1)%filter1( f, ff, this%xsz(2), this%xsz(3), bc1, bcn)
          endif
        case (2)
          if(associated(this%mbtopology)) then
            do imb = 1, this%mbtopology%x_num_blocks
              nxst = this%mbtopology%xst(1, imb);   nxen = this%mbtopology%xen(1, imb)
              nyst = this%mbtopology%xst(2, imb);   nyen = this%mbtopology%xen(2, imb)
              nzst = this%mbtopology%xst(3, imb);   nzen = this%mbtopology%xen(3, imb)
              call this%xgauf(imb) % filter1(f(nxst:nxen,nyst:nyen,nzst:nzen), &
                                            ff(nxst:nxen,nyst:nyen,nzst:nzen), &
                                            nyen-nyst+1, nzen-nzst+1, bc1, bcn)
            enddo
          else
            call this%xgauf(1)%filter1( f, ff, this%xsz(2), this%xsz(3), bc1, bcn)
          endif
        case (3)
            call this%xlsqf%filter1( f, ff, this%xsz(2), this%xsz(3))
        case (4)
            call this%xbox2%filter1( f, ff, this%xsz(2), this%xsz(3), bc1, bcn)
        end select

    end subroutine

    subroutine filtery(this, f, ff, bc1, bcn)
        class (filters), intent(in) :: this
        real(rkind), dimension(this%ysz(1),this%ysz(2), this%ysz(3)), intent(in)  :: f
        real(rkind), dimension(this%ysz(1),this%ysz(2), this%ysz(3)), intent(out) :: ff
        integer, optional, intent(in) :: bc1, bcn

        integer :: imb, nxst, nyst, nzst, nxen, nyen, nzen

        select case (this%ymethod)
        case (1)
          if(associated(this%mbtopology)) then
            do imb = 1, this%mbtopology%y_num_blocks
              nxst = this%mbtopology%yst(1, imb);   nxen = this%mbtopology%yen(1, imb)
              nyst = this%mbtopology%yst(2, imb);   nyen = this%mbtopology%yen(2, imb)
              nzst = this%mbtopology%yst(3, imb);   nzen = this%mbtopology%yen(3, imb)
              call this%ycf90(imb) % filter2(f(nxst:nxen,nyst:nyen,nzst:nzen), &
                                            ff(nxst:nxen,nyst:nyen,nzst:nzen), &
                                            nxen-nxst+1, nzen-nzst+1, bc1, bcn)
            enddo
          else
            call this%ycf90(1)%filter2( f, ff, this%ysz(1), this%ysz(3), bc1, bcn)
          endif
        case (2)
          if(associated(this%mbtopology)) then
            do imb = 1, this%mbtopology%y_num_blocks
              nxst = this%mbtopology%yst(1, imb);   nxen = this%mbtopology%yen(1, imb)
              nyst = this%mbtopology%yst(2, imb);   nyen = this%mbtopology%yen(2, imb)
              nzst = this%mbtopology%yst(3, imb);   nzen = this%mbtopology%yen(3, imb)
              call this%ygauf(imb) % filter2(f(nxst:nxen,nyst:nyen,nzst:nzen), &
                                            ff(nxst:nxen,nyst:nyen,nzst:nzen), &
                                            nxen-nxst+1, nzen-nzst+1, bc1, bcn)
            enddo
          else
            call this%ygauf(1)%filter2( f, ff, this%ysz(1), this%ysz(3), bc1, bcn)
          endif
        case (3)
            call this%ylsqf%filter2( f, ff, this%ysz(1), this%ysz(3))
        case (4)
            call this%ybox2%filter2( f, ff, this%ysz(1), this%ysz(3), bc1, bcn)
        end select

    end subroutine

    subroutine filterz(this, f, ff, bc1, bcn)
        class (filters), intent(in) :: this
        real(rkind), dimension(this%zsz(1),this%zsz(2), this%zsz(3)), intent(in)  :: f
        real(rkind), dimension(this%zsz(1),this%zsz(2), this%zsz(3)), intent(out) :: ff
        integer, optional, intent(in) :: bc1, bcn

        integer :: imb, nxst, nyst, nzst, nxen, nyen, nzen

        select case (this%zmethod)
        case (1)
          if(associated(this%mbtopology)) then
            do imb = 1, this%mbtopology%z_num_blocks
              nxst = this%mbtopology%zst(1, imb);   nxen = this%mbtopology%zen(1, imb)
              nyst = this%mbtopology%zst(2, imb);   nyen = this%mbtopology%zen(2, imb)
              nzst = this%mbtopology%zst(3, imb);   nzen = this%mbtopology%zen(3, imb)
              call this%zcf90(imb) % filter3(f(nxst:nxen,nyst:nyen,nzst:nzen), &
                                            ff(nxst:nxen,nyst:nyen,nzst:nzen), &
                                            nxen-nxst+1, nyen-nyst+1, bc1, bcn)
            enddo
          else
            call this%zcf90(1)%filter3( f, ff, this%zsz(1), this%zsz(2), bc1, bcn)
          endif
        case (2)
          if(associated(this%mbtopology)) then
            do imb = 1, this%mbtopology%z_num_blocks
              nxst = this%mbtopology%zst(1, imb);   nxen = this%mbtopology%zen(1, imb)
              nyst = this%mbtopology%zst(2, imb);   nyen = this%mbtopology%zen(2, imb)
              nzst = this%mbtopology%zst(3, imb);   nzen = this%mbtopology%zen(3, imb)
              call this%zgauf(imb) % filter3(f(nxst:nxen,nyst:nyen,nzst:nzen), &
                                            ff(nxst:nxen,nyst:nyen,nzst:nzen), &
                                            nxen-nxst+1, nyen-nyst+1, bc1, bcn)
            enddo
          else
            call this%zgauf(1)%filter3( f, ff, this%zsz(1), this%zsz(2), bc1, bcn)
          endif
        case (3)
            call this%zlsqf%filter3( f, ff, this%zsz(1), this%zsz(2))
        case (4)
            call this%zbox2%filter3( f, ff, this%zsz(1), this%zsz(2), bc1, bcn)
        end select

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! INTERNAL SUBROUTINES !!!!!!!!!!!!!!!!!!!!!!!!
    ! Don't change things below this point
    
    subroutine init_parallel(this,                              gp  , &
                                     periodicx, periodicy, periodicz, & 
                                     methodx  , methody  , methodz,   &
                                     mbtopology)

        class( filters )   , intent(inout) :: this
        class( decomp_info), intent(in)    :: gp 
        character(len=*)   , intent(in)    :: methodx, methody, methodz 
        logical            , intent(in)    :: periodicx, periodicy, periodicz
        class(multiblocktopol), intent(in), optional, target :: mbtopology

        if (this%initialized) then
            call message("WARNING: Reinitializing the FILTER class!")
            call this%destroy
        end if  
       
        this%xsz = gp%xsz
        this%ysz = gp%ysz
        this%zsz = gp%zsz

        if(present(mbtopology)) then
            this%mbtopology => mbtopology
        else
            this%mbtopology => null()
        endif

        if( (associated(this%mbtopology)) .and. &
            (.not. ( (methodx=='cf90') .or. (methodx=='gaussian') )) ) then
            call GracefulExit("Only cf90 and gaussian are supported in x with multi-block currently", 11)
        endif

        if( (associated(this%mbtopology)) .and. &
            (.not. ( (methody=='cf90') .or. (methody=='gaussian') )) ) then
            call GracefulExit("Only cf90 and gaussian are supported in y with multi-block currently", 11)
        endif

        if( (associated(this%mbtopology)) .and. &
            (.not. ( (methodz=='cf90') .or. (methodz=='gaussian') )) ) then
            call GracefulExit("Only cf90 and gaussian are supported in z with multi-block currently", 11)
        endif

        call this%init_procedures(  this%xsz(1),  this%ysz(2),  this%zsz(3), &
                                     methodx, methody, methodz, &
                                     periodicx, periodicy, periodicz)

        this%initialized = .true. 
    end subroutine

    subroutine init_serial(this,          nx  ,      ny  ,       nz , &
                                     periodicx, periodicy, periodicz, &
                                     methodx  , methody  , methodz)

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
