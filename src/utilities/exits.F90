module exits

    use kind_parameters, only: rkind,clen,stdout,stderr
    use constants, only: one
    use decomp_2d, only: nrank, decomp_2d_abort
    
    implicit none
    private
    public :: GracefulExit, message, message_min_max, warning, newline, nancheck, check_exit
        
    interface message
        module procedure message_char, message_char_double, message_char_int, &
          message_level_char, message_level_char_double, message_level_char_int, &
          message_char_int_char_int 
    end interface

    interface message_min_max
        module procedure message_minmax_int, message_minmax_real
    end interface
    interface warning
        module procedure warning_char
    end interface

    interface nancheck
        module procedure nancheck_std, nancheck_ind, nancheck_std_arr3
    end interface
     
contains

    subroutine GracefulExit(message, errcode)
        use mpi, only: MPI_COMM_WORLD, MPI_Barrier, MPI_Abort 
        character(len=*), intent(in) :: message
        integer, intent(in) :: errcode
        integer :: ierr
        if (nrank == 0) then
          print*, trim(message)
        end if
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
        !call decomp_2d_abort(errcode, message)
    end subroutine

    subroutine newline()
        if (nrank == 0) write(stdout,*)
    end subroutine

    subroutine message_char(mess)
        character(len=*), intent(in) :: mess
        if (nrank == 0) write(stdout,'(A)') mess
    end subroutine
    
    subroutine warning_char(mess)
        character(len=*), intent(in) :: mess
        if (nrank == 0) write(stderr,'(A)') mess
    end subroutine
    
    subroutine message_char_double(mess,val)
        character(len=*), intent(in) :: mess
        real(rkind), intent(in) :: val
        if (nrank == 0) write(stdout,'(A,A,ES26.16)') mess, " = ", val
    end subroutine
    
    subroutine message_char_int(mess,val)
        character(len=*), intent(in) :: mess
        integer, intent(in) :: val
        if (nrank == 0) write(stdout,'(A,A,I0)') mess, " = ", val
    end subroutine
    
    subroutine message_char_int_char_int(mess1,val1,mess2,val2)
        character(len=*), intent(in) :: mess1, mess2
        integer, intent(in) :: val1, val2
        if (nrank == 0) write(stdout,'(A,I0,A,I0)') mess1, val1, mess2, val2
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
        if (nrank == 0) write(stdout,'(A)') full_message

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
        if (nrank == 0) write(stdout,'(A,A,ES26.16)') full_message, " = ", val
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
        if (nrank == 0) write(stdout,'(A,A,I0)') full_message, " = ", val
    end subroutine

    subroutine message_minmax_int(level,mess,valmin, valmax)
        integer, intent(in) :: level
        character(len=*), intent(in) :: mess
        integer, intent(in) :: valmin, valmax
        character(:), allocatable :: full_message
        integer :: i
        
        full_message = ""
        do i=1,level
            full_message = full_message // "    "
        end do

        full_message = full_message // "> " // mess
        if (nrank == 0) write(stdout,'(A,A,2(I0,A))') full_message, " = (", valmin,",",valmax,")"
    end subroutine


    subroutine message_minmax_real(level,mess,valmin, valmax)
        integer, intent(in) :: level
        character(len=*), intent(in) :: mess
        real(rkind), intent(in) :: valmin, valmax
        character(:), allocatable :: full_message
        integer :: i
        
        full_message = ""
        do i=1,level
            full_message = full_message // "    "
        end do

        full_message = full_message // "> " // mess
        if (nrank == 0) write(stdout,'(A,A,2(ES26.16,A))') full_message, " = (", valmin,",",valmax,")"
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
    
    logical function check_exit(dir)
         character(len=*), intent(in) :: dir

         integer :: ierr

         open(777,file=trim(dir)//"/exitpdo",status='old',iostat=ierr)
         if(ierr==0) then
             check_exit = .TRUE.
         else
             check_exit = .FALSE.
         endif
         close(777)

    end function

end module 
