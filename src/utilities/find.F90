module findMod
  use kind_parameters, only: rkind
  implicit none
  contains
    pure subroutine find(xin,xcompare,xout)
      ! Returns array with the index of every item in xin that is equal to
      ! xcompare
      integer, dimension(:), intent(in) :: xin
      integer, intent(in) :: xcompare
      integer, dimension(:), allocatable, intent(out) :: xout
      integer :: i, counter, idx

      counter = 0
      do i = 1,size(xin)
        if (xin(i) == xcompare) counter = counter + 1
      end do

      allocate(xout(counter))
      idx = 1

      do i = 1,size(xin)
        if (xin(i) == xcompare) then
          xout(idx) = i
          idx = idx + 1
        end if
      end do
    end subroutine 
end module
