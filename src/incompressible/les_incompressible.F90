module les_incompressible
    use kind_parameters, only: rkind
    use standard_smagorinsky, only: smag
    use dynamic_smagorinsky, only: dsmag
    use exits, only: GracefulExit
    use decomp_2d, only: decomp_info
    use constants, only: zero
    use spectralMod, only: spectral  

    implicit none 
    private

    public :: sgs

    type :: sgs
        private
        type(smag), allocatable :: my_smag
        type(dsmag), allocatable :: my_dsmag
        
        integer :: model = 1                    ! 1: Standard Smagorinsky, 2: Dynamic Smagorinsky
        integer :: test_filter = 1              ! 1: least squares, 2: gaussian, 3: sharp spectral
        
        integer :: nxhat, nyhat, nzhat
        integer :: nx, ny, nz
        
        logical :: isInitialized

        contains
            procedure :: init
            procedure :: destroy

            procedure :: getVelSource    

    end type 

contains
    
    subroutine init(this,method,isZperiodic, decomp, spect)
        class(sgs), intent(inout) :: this
        character(len=*), intent(in) :: method
        logical, intent(in) :: isZperiodic
        type(decomp_info), intent(in) :: decomp
        class(spectral), intent(in) :: spect 

        if (.not. isZperiodic) then
            call GracefulExit("CODE INCOMPLETE: non-periodic z is not supported yet", 312)
        end if 

        select case (method)
        case ("NONE")
            this%model = 0
        case ("SMAG")
            allocate(this%my_smag)
            this%model = 1
        case ("DSMAG")
            allocate(this%my_dsmag)
            this%model = 2
        case default
            call GracefulExit("Incorrect Subgrid model selected. Check the &
                    & available options is les_incompressible.F90",101)
        end select
        
        this%nx = decomp%xsz(1)
        this%ny = decomp%xsz(2)
        this%nz = decomp%xsz(3)

        this%nxhat = size(spect%kabs_sq,1)
        this%nyhat = size(spect%kabs_sq,2)
        this%nzhat = size(spect%kabs_sq,3)

        this%isInitialized = .true. 

    end subroutine 

    subroutine destroy(this)
        class(sgs), intent(inout) :: this

        select case (this%model) 
        case (1)
            call this%my_smag%destroy()
            deallocate(this%my_smag)
        case (2)
            call this%my_dsmag%destroy()
            deallocate(this%my_dsmag)
        end select 
        
        this%isInitialized = .false. 

    end subroutine 

    subroutine getVelSource(this,Sfields,SourceVal)
        class(sgs), intent(in) :: this 
        complex(rkind), intent(in), dimension(this%nxhat, this%nyhat, this%nzhat,3) :: Sfields
        complex(rkind), intent(out), dimension(this%nxhat, this%nyhat, this%nzhat,3) :: SourceVal

        select case (this%model) 
        case (0)
            SourceVal = zero
        case (1)
            SourceVal = Sfields - Sfields 
        case (2)
            SourceVal = Sfields - Sfields 
        end select 

    end subroutine



end module 
