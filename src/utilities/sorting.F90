module sorting_mod 
use kind_parameters, only: rkind 
implicit none

! QUICK SORT: Code taken from rosettacode.org
! url: https://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran

type sortgroup
    real(rkind) :: value  = 0._rkind  ! values to be sorted by
    integer     :: zpos   = 1         ! z coordinate 
end type sortgroup
 
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
 
end module sorting_mod 
