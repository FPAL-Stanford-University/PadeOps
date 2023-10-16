module fortran_assert
  use exits, only: gracefulExit
  use kind_parameters, only: clen, rkind
  implicit none
contains

  subroutine assert(expression, description, val)
    logical, intent(in) :: expression
    character(len=*), intent(in), optional :: description
    real(rkind), intent(in), optional :: val
    integer :: err_code = 1
    character(len=clen) :: root_stmt, mssg

    root_stmt = 'assert statement failed:'
    if (.not. expression) then
      if (present(description)) then
        if (present(val)) then
          write(mssg,'(A,F4.9)')trim(root_stmt)//" "//trim(description),val
        else
          mssg = trim(root_stmt)//" "//trim(description)
        end if
      else
        mssg = trim(root_stmt)
      end if
      call gracefulExit(trim(mssg), err_code)
    end if
  end subroutine assert
  
end module fortran_assert
