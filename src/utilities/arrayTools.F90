module arrayTools
  use kind_parameters, only: rkind
  implicit none
  
  interface findEQ
    module procedure findEQ1
  end interface
  interface findGL
    module procedure findGL3
  end interface
  interface findLE
    module procedure findLE1first
  end interface
  interface findGE
    module procedure findGE1first
  end interface
 
contains
  
  pure subroutine findGE1first(xin,xcompare,idx)
    ! Returns array with the index of every item (or the first N)
    ! in xin that is greater than or equal to xcompare
    real(rkind), dimension(:), intent(in) :: xin 
    real(rkind), intent(in) :: xcompare
    integer,  intent(out) :: idx
    integer :: i
    
    do i = 1,size(xin)
      if (xin(i) >= xcompare) then
        idx = i 
        exit
      end if
    end do
  end subroutine

  pure subroutine findLE1first(xin,xcompare,idx)
    ! Returns array with the index of every item (or the first N)
    ! in xin that is less than or equal to xcompare
    real(rkind), dimension(:), intent(in) :: xin 
    real(rkind), intent(in) :: xcompare
    integer,  intent(out) :: idx
    integer :: i
    
    do i = 1,size(xin)
      if (xin(i) <= xcompare) then
        idx = i 
        exit
      end if
    end do
  end subroutine

  pure subroutine findEQ1(xin,xcompare,xout)
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
  
  pure subroutine findGL3(xin,xmin,xmax,iarr,jarr,karr)
    ! Returns array with the index of every item in xin that is equal to
    ! xcompare
    real(rkind), dimension(:,:,:), intent(in) :: xin 
    real(rkind), intent(in) :: xmin, xmax
    integer, dimension(:), allocatable, intent(out) :: iarr, jarr, karr
    integer :: i, j, k, counter, idx, nx, ny, nz

    nx = size(xin,1)
    ny = size(xin,2)
    nz = size(xin,3)
 
    counter = 0 
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (xin(i,j,k) >= xmin .and. xin(i,j,k) <= xmax) then
            counter = counter + 1
          end if
        end do
      end do
    end do
 
    allocate(iarr(counter),jarr(counter),karr(counter))
    idx = 1 

    if (counter > 0) then 
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            if (xin(i,j,k) >= xmin .and. xin(i,j,k) <= xmax) then
              iarr(idx) = i
              jarr(idx) = j
              karr(idx) = k
              idx = idx + 1
            end if 
          end do
        end do
      end do
    end if
  end subroutine

  subroutine swap(A,i,j)
    ! Swap A(i) with A(j)
    integer, dimension(:), intent(inout) :: A
    integer, intent(in) :: i, j
    integer, dimension(size(A)) :: B
    
    B = A
    B(i) = A(j)
    B(j) = A(i)
    A = B
  end subroutine 
end module
