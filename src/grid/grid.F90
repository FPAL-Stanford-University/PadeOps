module GridMod

    use kind_parameters, only: rkind,clen
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters

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

    contains

        procedure(init_interface),    deferred :: init
        procedure(destroy_interface), deferred :: destroy

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
    
    end interface

end module
