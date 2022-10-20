module fortran_assert
  use exits, only: gracefulExit
  use kind_parameters, only: clen
  implicit none
contains

  subroutine assert(expression, description, pe)
    logical, intent(in) :: expression
    character(len=*), intent(in), optional :: description
    integer, intent(in), optional :: pe
    integer :: err_code
    character(len=clen) :: root_stmt

    if (present(pe)) then
      write(root_stmt,'(A,I2.2,A)')'PE', pe, ' -- assert statement failed:'
    else
      root_stmt = 'Assert statement FAILED:'
    end if
    
    if (.not. expression) then
      if (present(description)) then
        call gracefulExit(trim(root_stmt)//" "//description, err_code)
      else
        call gracefulExit(trim(root_stmt), err_code)
      end if
    end if
  end subroutine assert
  
  subroutine check(expression, description, pass, pe)
    logical, intent(in) :: expression
    character(len=*), intent(in), optional :: description
    logical, intent(inout) :: pass
    integer, intent(in), optional :: pe
    integer :: err_code
    character(len=clen) :: root_stmt

    if (present(pe)) then
      write(root_stmt,'(A,I2.2,A)')'PE', pe, ' -- WARNING:'
    else
      root_stmt = 'WARNING:'
    end if
    if (.not. expression) then
      pass = .false.
      if (present(description)) then
        print*, trim(root_stmt)//" "//description
      else
        print*, trim(root_stmt)
      end if
    end if
  end subroutine

end module fortran_assert
