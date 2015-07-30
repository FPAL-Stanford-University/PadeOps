module exits

    use decomp_2d, only: decomp_2d_abort
    
    implicit none
    private
    public :: GracefulExit
     
contains

    subroutine GracefulExit(message, errcode)
        
        character(len=*), intent(in) :: message
        integer, intent(in) :: errcode

        call decomp_2d_abort(errcode, message)
    
    end subroutine
    
end module 
