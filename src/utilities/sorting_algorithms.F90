module sorting_algorithms
  use kind_parameters, only: rkind
  implicit none

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

    pure subroutine integer_binary_sort1D(data)
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

    pure subroutine integer_binary_sort2D(data,arg)
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

    pure subroutine real_binary_sort1D(data)
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

    pure subroutine real_binary_sort2D(data,arg)
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
end module sorting_algorithms
