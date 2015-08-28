module dynamic_smagorinsky
    use kind_parameters, only: rkind
    
    implicit none
    private

    public :: dsmag

    type :: dsmag
        private
        real(rkind), dimension(:,:,:,:), allocatable :: Sij 
        
        contains
            procedure :: init
            procedure :: destroy 
    end type 

contains

    subroutine init(this)
        class(dsmag), intent(inout) :: this

    end subroutine

    subroutine destroy(this)
        class(dsmag), intent(inout) :: this

    end subroutine 

end module 
