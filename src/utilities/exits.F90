module exits

    use kind_parameters, only: rkind,stdout
    use decomp_2d, only: nrank, decomp_2d_abort
    
    implicit none
    private
    public :: GracefulExit, message
        
    interface message
        module procedure message_char, message_char_double
    end interface
     
contains

    subroutine GracefulExit(message, errcode)
        
        character(len=*), intent(in) :: message
        integer, intent(in) :: errcode

        call decomp_2d_abort(errcode, message)
    
    end subroutine

    subroutine message_char(mess)
        character(len=*), intent(in) :: mess
        if (nrank == 0) write(stdout,*) mess
    end subroutine
    
    subroutine message_char_double(mess,val)
        character(len=*), intent(in) :: mess
        real(rkind), intent(in) :: val
        if (nrank == 0) write(stdout,*) mess, " = ", val
    end subroutine
    
end module 
