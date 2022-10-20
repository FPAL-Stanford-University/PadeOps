program test_find
  use findMod, only: find
  use fortran_assert, only: assert
  implicit none

  integer, dimension(10) :: x
  integer :: xcompare = 3
  integer, dimension(:), allocatable :: idx

  x = [1,2,3,1,2,3,4,5,6,3]
  call find(x,xcompare,idx)
  call assert(size(idx) == 3)
  print*, "Test PASSED! idx = ", idx

end program
