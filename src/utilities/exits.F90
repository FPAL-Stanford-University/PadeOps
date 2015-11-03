module exits

    use kind_parameters, only: rkind,clen,stdout,stderr
    use decomp_2d, only: nrank, decomp_2d_abort
    
    implicit none
    private
    public :: GracefulExit, message, warning, newline
        
    interface message
        module procedure message_char, message_char_double, message_level_char, message_level_char_double, message_level_char_int
    end interface

    interface warning
        module procedure warning_char
    end interface
     
contains

    subroutine GracefulExit(message, errcode)
        
        character(len=*), intent(in) :: message
        integer, intent(in) :: errcode

        call decomp_2d_abort(errcode, message)
    
    end subroutine

    subroutine newline()
        if (nrank == 0) write(stdout,*)
    end subroutine

    subroutine message_char(mess)
        character(len=*), intent(in) :: mess
        if (nrank == 0) write(stdout,*) mess
    end subroutine
    
    subroutine warning_char(mess)
        character(len=*), intent(in) :: mess
        if (nrank == 0) write(stderr,*) mess
    end subroutine
    
    subroutine message_char_double(mess,val)
        character(len=*), intent(in) :: mess
        real(rkind), intent(in) :: val
        if (nrank == 0) write(stdout,*) mess, " = ", val
    end subroutine
    
    subroutine message_level_char(level,mess)
        integer, intent(in) :: level
        character(len=*), intent(in) :: mess
        character(:), allocatable :: full_message
        integer :: i

        full_message = ""
        do i=1,level
            full_message = full_message // "    "
        end do

        full_message = full_message // "> " // mess
        if (nrank == 0) write(stdout,*) full_message

    end subroutine
    
    subroutine message_level_char_double(level,mess,val)
        integer, intent(in) :: level
        character(len=*), intent(in) :: mess
        real(rkind), intent(in) :: val
        character(:), allocatable :: full_message
        integer :: i
        
        full_message = ""
        do i=1,level
            full_message = full_message // "    "
        end do

        full_message = full_message // "> " // mess
        if (nrank == 0) write(stdout,*) full_message, " = ", val
    end subroutine
    
    subroutine message_level_char_int(level,mess,val)
        integer, intent(in) :: level
        character(len=*), intent(in) :: mess
        integer, intent(in) :: val
        character(:), allocatable :: full_message
        integer :: i
        
        full_message = ""
        do i=1,level
            full_message = full_message // "    "
        end do

        full_message = full_message // "> " // mess
        if (nrank == 0) write(stdout,*) full_message, " = ", val
    end subroutine
    
end module 
