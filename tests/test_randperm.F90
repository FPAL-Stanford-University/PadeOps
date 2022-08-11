program test_randperm
  use kind_parameters, only: rkind
  use random, only: randperm
  implicit none

  real(rkind), dimension(100) :: A
  integer, dimension(100) :: B
  integer :: n

  do n = 1,100
    A(n) = real(n,rkind)
    B(n) = n
  end do

  print*, "A:", A
  print*, "B:", B
  call randperm(A,1)
  call randperm(B,1)
  print*, "A:", A
  print*, "B:", B
  call randperm(A,45)
  call randperm(B,45)
  print*, "A:", A
  print*, "B:", B

  print*, "huge(0):", huge(0)
end program test_randperm
