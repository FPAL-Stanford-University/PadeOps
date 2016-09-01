module exits

    use kind_parameters, only: rkind,clen,stdout,stderr
    use constants, only: one
    use decomp_2d, only: nrank, decomp_2d_abort
    
    implicit none
    private
    public :: GracefulExit, message, warning, newline, nancheck
        
    interface message
        module procedure message_char, message_char_double, message_level_char, message_level_char_double, message_level_char_int
    end interface

    interface warning
        module procedure warning_char
    end interface

    interface nancheck
        module procedure nancheck_std, nancheck_ind, nancheck_std_arr3, nancheck_ind_arr3
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

    logical function nancheck_std(f) result(nancheck)
        real(rkind), dimension(:,:,:,:), intent(in) :: f
        integer :: i,j,k,l
        
        nancheck = .FALSE.
        do l = 1,size(f,4)
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    do i = 1,size(f,1)
                        if ( isnan(f(i,j,k,l)) .OR. ( f(i,j,k,l) + one == f(i,j,k,l) ) ) then
                            nancheck = .TRUE.
                            exit
                        end if
                    end do
                    if (nancheck) exit
                end do
                if (nancheck) exit
            end do
            if (nancheck) exit
        end do

    end function

    logical function nancheck_std_arr3(f) result(nancheck)
        real(rkind), dimension(:,:,:), intent(in) :: f
        integer :: i,j,k
        
        nancheck = .FALSE.
        do k = 1,size(f,3)
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    if ( isnan(f(i,j,k)) .OR. ( f(i,j,k) + one == f(i,j,k) ) ) then
                        nancheck = .TRUE.
                        exit
                    end if
                end do
                if (nancheck) exit
            end do
            if (nancheck) exit
        end do

    end function

    logical function nancheck_ind(f,i,j,k,l) result(nancheck)
        real(rkind), dimension(:,:,:,:), intent(in) :: f
        integer, intent(out) :: i,j,k,l
        
        nancheck = .FALSE.
        do l = 1,size(f,4)
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    do i = 1,size(f,1)
                        if ( isnan(f(i,j,k,l)) .OR. ( f(i,j,k,l) + one == f(i,j,k,l) ) ) then
                            nancheck = .TRUE.
                            exit
                        end if
                    end do
                    if (nancheck) exit
                end do
                if (nancheck) exit
            end do
            if (nancheck) exit
        end do

    end function
    
    logical function nancheck_ind_arr3(f,i,j,k) result(nancheck)
        real(rkind), dimension(:,:,:), intent(in) :: f
        integer, intent(out) :: i,j,k
       
        nancheck = .FALSE.
        do k = 1,size(f,3)
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    if ( isnan(f(i,j,k)) .OR. ( f(i,j,k) + one == f(i,j,k) ) ) then
                        nancheck = .TRUE.
                        exit
                    end if
                end do
                if (nancheck) exit
            end do
            if (nancheck) exit
        end do

    end function
    
end module 
