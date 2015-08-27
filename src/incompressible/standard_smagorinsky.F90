module standard_smagorinsky
    use kind_parameters, only: rkind
    
    implicit none
    private

    public :: smag

    type :: smag
        private
        real(rkind), dimension(:,:,:,:), allocatable :: Sij 
        
        contains
            procedure :: init
            procedure :: destroy

    end type 

contains

    subroutine init(this)
        class(smag), intent(inout) :: this

    end subroutine


    subroutine destroy(this)
        class(smag), intent(inout) :: this

    end subroutine 

end module 
