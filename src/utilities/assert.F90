module fortran_assert
  use exits, only: gracefulExit
  use kind_parameters, only: clen
  implicit none
contains

  subroutine assert(expression, description)
    logical, intent(in) :: expression
    character(len=*), intent(in), optional :: description
    integer :: err_code
    character(len=clen) :: root_stmt

    root_stmt = 'assert statement failed:'
    if (.not. expression) then
      if (present(description)) then
        call gracefulExit(trim(root_stmt)//" "//description, err_code)
      else
        call gracefulExit(trim(root_stmt), err_code)
      end if
    end if
  end subroutine assert
  
end module fortran_assert
