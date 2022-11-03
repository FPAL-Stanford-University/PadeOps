module sorting_mod 
use kind_parameters, only: rkind 
implicit none

! QUICK SORT: Code taken from rosettacode.org
! url: https://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran

type sortgroup
    real(rkind) :: value  = 0._rkind  ! values to be sorted by
    integer     :: zpos   = 1         ! z coordinate 
end type sortgroup

 interface binary_sort
    !--------------------------------------------------------------!
    ! Binary sort compares elementwise neighbors and sorts them.
    ! It recursively passes over the array until the entire array
    ! is sorted. The 2D implementations take a 2D array and sort 
    ! the rows according to the value of the specified column (i.e. 
    ! "arg" in the input parameters).
    !--------------------------------------------------------------!
    module procedure integer_binary_sort1D, integer_binary_sort2D,&
                        real_binary_sort1D,    real_binary_sort2D
  end interface
 
contains
 
recursive subroutine QSort(a,na)
 
! DUMMY ARGUMENTS
integer, intent(in) :: nA
type (sortgroup), dimension(nA), intent(in out) :: A
 
! LOCAL VARIABLES
integer :: left, right
real(rkind) :: random
real(rkind) :: pivot
type (sortgroup) :: temp
integer :: marker
 
    if (nA > 1) then
 
        call random_number(random)
        pivot = A(int(random*real(nA-1))+1)%value   ! random pivor (not best performance, but avoids worst-case)
        left = 0
        right = nA + 1
 
        do while (left < right)
            right = right - 1
            do while (A(right)%value > pivot)
                right = right - 1
            end do
            left = left + 1
            do while (A(left)%value < pivot)
                left = left + 1
            end do
            if (left < right) then
                temp = A(left)
                A(left) = A(right)
                A(right) = temp
            end if
        end do
 
        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if
 
        call QSort(A(:marker-1),marker-1)
        call QSort(A(marker:),nA-marker+1)
 
    end if
 
end subroutine QSort

    subroutine integer_binary_sort1D(data)
      integer, dimension(:), intent(inout) :: data
      logical :: cont
      integer :: i, a, hold

      cont = .true.
      do while (cont)
        a = 0
        do i = 1,size(data)-1
          hold = data(i)
          if (data(i) > data(i+1)) then
            data(i)   = data(i+1)
            data(i+1) = hold
            a = 1 
          end if
        end do
        if (a == 0) cont = .false.
      end do
    end subroutine

    subroutine integer_binary_sort2D(data,arg)
      integer, dimension(:,:), intent(inout) :: data
      logical :: cont
      integer :: i, a
      integer, intent(in) :: arg 
      integer, dimension(2) :: hold
    
      cont = .true.
      do while (cont)
        a = 0 
        do i = 1,size(data,1)-1
          hold = data(i,:)
          if (data(i,arg) > data(i+1,arg)) then
            data(i,:)   = data(i+1,:)
            data(i+1,:) = hold
            a = 1 
          end if
        end do
        if (a == 0) cont = .false.
      end do
    end subroutine

    subroutine real_binary_sort1D(data)
      real(rkind), dimension(:), intent(inout) :: data
      logical :: cont
      integer :: i, a
      real(rkind) :: hold

      cont = .true.
      do while (cont)
        a = 0 
        do i = 1,size(data)-1
          hold = data(i)
          if (data(i) > data(i+1)) then
            data(i)   = data(i+1)
            data(i+1) = hold
            a = 1
          end if
        end do
        if (a == 0) cont = .false.
      end do
    end subroutine

    subroutine real_binary_sort2D(data,arg)
      real(rkind), dimension(:,:), intent(inout) :: data
      logical :: cont
      integer :: i, a
      integer, intent(in) :: arg
      real(rkind), dimension(2) :: hold

      cont = .true.
      do while (cont)
        a = 0
        do i = 1,size(data,1)-1
          hold = data(i,:)
          if (data(i,arg) > data(i+1,arg)) then
            data(i,:)   = data(i+1,:)
            data(i+1,:) = hold
            a = 1
          end if
        end do
        if (a == 0) cont = .false.
      end do
    end subroutine
 
end module sorting_mod 
