module GridMod

    use kind_parameters, only: rkind,clen
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit

    type, abstract :: grid

        integer                                              :: nx                ! Number of points in X
        integer                                              :: ny                ! Number of points in Y
        integer                                              :: nz                ! Number of points in Z

        real(rkind)                                          :: dx                ! Grid spacing in X
        real(rkind)                                          :: dy                ! Grid spacing in Y
        real(rkind)                                          :: dz                ! Grid spacing in Z

        logical                                              :: periodicx         ! Periodic in X?
        logical                                              :: periodicy         ! Periodic in Y?
        logical                                              :: periodicz         ! Periodic in Z?

        type( derivatives )                                  :: der               ! Derivative object
        character(len=clen)                                  :: derivative_x      ! What derivative to use in X: "cd10", "cd06", "four", "cheb"
        character(len=clen)                                  :: derivative_y      ! What derivative to use in Y: "cd10", "cd06", "four", "cheb"
        character(len=clen)                                  :: derivative_z      ! What derivative to use in Z: "cd10", "cd06", "four", "cheb"

        type( filters )                                      :: fil
        character(len=clen)                                  :: filter_x          ! What filter to use in X: "cf90", "gaussian", "lstsq", "spectral"
        character(len=clen)                                  :: filter_y          ! What filter to use in X: "cf90", "gaussian", "lstsq", "spectral" 
        character(len=clen)                                  :: filter_z          ! What filter to use in X: "cf90", "gaussian", "lstsq", "spectral" 

        character(len=clen)                                  :: outputdir         ! Directory for output files

        type( decomp_info )                                  :: decomp

        real(rkind), dimension(:,:,:,:), allocatable         :: mesh
        real(rkind), dimension(:,:,:,:), allocatable         :: fields

        integer                                              :: nx_proc, ny_proc, nz_proc
    contains

        procedure(init_interface),    deferred :: init
        procedure(destroy_interface), deferred :: destroy
        procedure(laplacian_interface), deferred :: laplacian
        procedure(gradient_interface), deferred :: gradient

    end type

    abstract interface

        subroutine init_interface(this, inputfile)
            import :: grid
            import :: clen
            class(grid), intent(inout) :: this
            character(len=clen), intent(in) :: inputfile
        end subroutine

        subroutine destroy_interface(this)
            import :: grid
            class(grid), intent(inout) :: this
        end subroutine
    
        subroutine laplacian_interface(this, f, lapf)
            import :: grid
            import :: rkind
            class(grid), intent(in) :: this
            real(rkind), dimension(this%nx_proc, this%ny_proc, this%nz_proc), intent(in):: f   
            real(rkind), dimension(this%nx_proc, this%ny_proc, this%nz_proc), intent(out):: lapf
        end subroutine

        subroutine gradient_interface(this, f, dfdx, dfdy, dfdz)
            import :: grid
            import :: rkind
            class(grid), intent(in) :: this
            real(rkind), dimension(this%nx_proc, this%ny_proc, this%nz_proc), intent(in):: f   
            real(rkind), dimension(this%nx_proc, this%ny_proc, this%nz_proc), intent(out):: dfdx
            real(rkind), dimension(this%nx_proc, this%ny_proc, this%nz_proc), intent(out):: dfdy
            real(rkind), dimension(this%nx_proc, this%ny_proc, this%nz_proc), intent(out):: dfdz

        end subroutine

    end interface

    interface destroy_buffs
        module procedure destroy_buffs_real, destroy_buffs_complex 
    end interface
    
    interface alloc_buffs
        module procedure alloc_buffs_real, alloc_buffs_complex 
    end interface
contains

    subroutine destroy_buffs_real(buff)
        real(rkind), dimension(:,:,:,:),allocatable, intent(inout) :: buff

        if (allocated(buff)) deallocate(buff)
    end subroutine

    subroutine destroy_buffs_complex(buff)
        complex(rkind), dimension(:,:,:,:),allocatable, intent(inout) :: buff

        if (allocated(buff)) deallocate(buff)
    end subroutine

    subroutine alloc_buffs_real(buff,vars,dir,decomp)
        character(len=1), intent(in) :: dir
        class(decomp_info), intent(in) :: decomp
        integer, intent(in) :: vars
        real(rkind), dimension(:,:,:,:), allocatable, intent(out) :: buff

        if (allocated(buff)) deallocate(buff)

        select case (dir)
        case ("x")
            allocate(buff(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3),vars))
        case("y")
            allocate(buff(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),vars))
        case("z")
            allocate(buff(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3),vars))
        case default
            call GracefulExit("Incorrect direction selected in ALLOC_BUFFS_REAL subroutine", 13123)
        end select

    end subroutine

    subroutine alloc_buffs_complex(buff,vars,dir,decomp)
        character(len=1), intent(in) :: dir
        class(decomp_info), intent(in) :: decomp
        integer, intent(in) :: vars
        complex(rkind), dimension(:,:,:,:), allocatable, intent(out) :: buff

        if (allocated(buff)) deallocate(buff)

        select case (dir)
        case ("x")
            allocate(buff(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3),vars))
        case("y")
            allocate(buff(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),vars))
        case("z")
            allocate(buff(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3),vars))
        case default
            call GracefulExit("Incorrect direction selected in ALLOC_BUFFS_REAL subroutine", 13123)
        end select

    end subroutine

end module
