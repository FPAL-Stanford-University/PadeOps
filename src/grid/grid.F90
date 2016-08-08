module GridMod

    use kind_parameters, only: rkind,clen
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use io_VTK_stuff,    only: io_VTK
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

        type( derivatives ), allocatable                     :: der               ! Derivative object
        character(len=clen)                                  :: derivative_x      ! What derivative to use in X: "cd10", "cd06", "four", "cheb"
        character(len=clen)                                  :: derivative_y      ! What derivative to use in Y: "cd10", "cd06", "four", "cheb"
        character(len=clen)                                  :: derivative_z      ! What derivative to use in Z: "cd10", "cd06", "four", "cheb"

        type( filters )    , allocatable                     :: fil
        character(len=clen)                                  :: filter_x          ! What filter to use in X: "cf90", "gaussian", "lstsq", "spectral"
        character(len=clen)                                  :: filter_y          ! What filter to use in X: "cf90", "gaussian", "lstsq", "spectral" 
        character(len=clen)                                  :: filter_z          ! What filter to use in X: "cf90", "gaussian", "lstsq", "spectral" 

        character(len=clen)                                  :: outputdir         ! Directory for output files

        type( decomp_info ), allocatable                     :: decomp

        real(rkind), dimension(:,:,:,:), allocatable         :: mesh
        real(rkind), dimension(:,:,:,:), allocatable         :: fields

        real(rkind)                                          :: tstop, dt, dtfixed, tsim, CFL
        integer                                              :: step, nsteps
        integer                                              :: nxp, nyp, nzp

        integer                                              :: t_dataDump, t_restartDump
        
        logical                                              :: SkewSymm 
        logical                                              :: ViscConsrv         ! Is the viscous term being computed using the conservative formulation? 
        
        type( io_VTK ), allocatable                          :: viz
        real(rkind)                                          :: tviz

    contains

        procedure(init_interface),    deferred :: init
        procedure(destroy_interface), deferred :: destroy_grid
        ! procedure(laplacian_interface), deferred :: laplacian
        ! procedure(gradient_interface), deferred :: gradient

    end type

    abstract interface

        subroutine init_interface(this, inputfile)
            import :: grid
            import :: clen
            class(grid),target, intent(inout) :: this
            character(len=clen), intent(in) :: inputfile
        end subroutine

        subroutine destroy_interface(this)
            import :: grid
            class(grid), intent(inout) :: this
        end subroutine
    
        ! subroutine laplacian_interface(this, f, lapf, x_bc, y_bc, z_bc)
        !     import :: grid
        !     import :: rkind
        !     class(grid), target, intent(inout) :: this
        !     real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in):: f   
        !     real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(out):: lapf
        !     integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc
        ! end subroutine

        ! subroutine gradient_interface(this, f, dfdx, dfdy, dfdz, x_bc, y_bc, z_bc)
        !     import :: grid
        !     import :: rkind
        !     class(grid), target, intent(inout) :: this
        !     real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(in):: f   
        !     real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(out):: dfdx
        !     real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(out):: dfdy
        !     real(rkind), dimension(this%nxp, this%nyp, this%nzp), intent(out):: dfdz
        !     integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc

        ! end subroutine

    end interface

contains

end module
