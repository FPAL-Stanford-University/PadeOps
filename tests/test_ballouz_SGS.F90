program test_ballouz_SGS
  use mpi
  use decomp_2d
  use kind_parameters, only: rkind, clen
  use sgsmod_igrid,    only: getRotationTensor, getInverseRotationTensor, &
                             contract_Rik_Skl_invRlj
  use constants,       only: pi
  use random,          only: uniform_random
  use fortran_assert,  only: assert
  implicit none
  
  real(rkind), dimension(3,3) :: A, B
  real(rkind), dimension(1,1,1,3,3) :: A2, C
  real(rkind), dimension(1,1,1,6) :: B2
  real(rkind), dimension(1,1,1) :: D
  real(rkind), dimension(1,1,1,3,3) :: Rij, invRij
  real(rkind), dimension(1,1,1) :: phi, theta, psi
  integer :: i, j, k

  A(1,1) = 1.d0
  A(1,2) = 0.d0
  A(1,3) = 0.d0

  A(2,1) = 0.d0
  A(2,2) = 1.d0
  A(2,3) = 0.d0
  
  A(3,1) = 0.d0
  A(3,2) = 0.d0
  A(3,3) = 1.d0

  phi = pi/2.d0
  theta = 0.d0
  psi = 0.d0

  call getRotationTensor(Rij,phi,theta,psi)

  do j = 1,3
    do  i = 1,3
      B(i,j) = 0.d0
      do k = 1,3
        B(i,j) = B(i,j) + Rij(1,1,1,i,k)*A(k,j)
      end do
    end do
  end do

  print*, "For phi = pi/2, B = "
  print*, nint(B(1,:))
  print*, nint(B(2,:))
  print*, nint(B(3,:))
      
  phi = 0.d0
  theta = pi/2.d0
  psi = 0.d0

  call getRotationTensor(Rij,phi,theta,psi)

  do j = 1,3
    do  i = 1,3
      B(i,j) = 0.d0
      do k = 1,3
        B(i,j) = B(i,j) + Rij(1,1,1,i,k)*A(k,j)
      end do
    end do
  end do

  print*, "For theta = pi/2, B = "
  print*, nint(B(1,:))
  print*, nint(B(2,:))
  print*, nint(B(3,:))
      
  phi = 0.d0
  theta = 0.d0
  psi = pi/2.d0

  call getRotationTensor(Rij,phi,theta,psi)

  do j = 1,3
    do  i = 1,3
      B(i,j) = 0.d0
      do k = 1,3
        B(i,j) = B(i,j) + Rij(1,1,1,i,k)*A(k,j)
      end do
    end do
  end do

  print*, "For psi = pi/2, B = "
  print*, nint(B(1,:))
  print*, nint(B(2,:))
  print*, nint(B(3,:))
      
  phi = pi/2.d0
  theta = pi/2.d0
  psi = 0.d0

  call getRotationTensor(Rij,phi,theta,psi)

  do j = 1,3
    do  i = 1,3
      B(i,j) = 0.d0
      do k = 1,3
        B(i,j) = B(i,j) + Rij(1,1,1,i,k)*A(k,j)
      end do
    end do
  end do

  print*, "For (pi/2, pi/2, 0) B = "
  print*, nint(B(1,:))
  print*, nint(B(2,:))
  print*, nint(B(3,:))
      
  phi = pi/2.d0
  theta = 0.d0
  psi = pi/2.d0

  call getRotationTensor(Rij,phi,theta,psi)

  do j = 1,3
    do  i = 1,3
      B(i,j) = 0.d0
      do k = 1,3
        B(i,j) = B(i,j) + Rij(1,1,1,i,k)*A(k,j)
      end do
    end do
  end do

  print*, "For (pi/2, 0, pi/2) B = "
  print*, nint(B(1,:))
  print*, nint(B(2,:))
  print*, nint(B(3,:))
      
  phi = 0.d0
  theta = pi/2.d0
  psi = pi/2.d0

  call getRotationTensor(Rij,phi,theta,psi)

  do j = 1,3
    do  i = 1,3
      B(i,j) = 0.d0
      do k = 1,3
        B(i,j) = B(i,j) + Rij(1,1,1,i,k)*A(k,j)
      end do
    end do
  end do

  print*, "For (0, pi/2, pi/2) B = "
  print*, nint(B(1,:))
  print*, nint(B(2,:))
  print*, nint(B(3,:))
      
  phi = pi/2.d0
  theta = pi/2.d0
  psi = pi/2.d0

  call getRotationTensor(Rij,phi,theta,psi)

  do j = 1,3
    do  i = 1,3
      B(i,j) = 0.d0
      do k = 1,3
        B(i,j) = B(i,j) + Rij(1,1,1,i,k)*A(k,j)
      end do
    end do
  end do

  print*, "For (pi/2, pi/2, pi/2) B = "
  print*, nint(B(1,:))
  print*, nint(B(2,:))
  print*, nint(B(3,:))

  call uniform_random(phi,0.d0,2.d0*pi)
  call uniform_random(theta,0.d0,2.d0*pi)
  call uniform_random(psi,0.d0,2.d0*pi)
  call getRotationTensor(Rij,phi,theta,psi)
  call getInverseRotationTensor(invRij,phi,theta,psi)

  do j = 1,3
    do  i = 1,3
      B(i,j) = 0.d0
      do k = 1,3
        B(i,j) = B(i,j) + Rij(1,1,1,i,k)*invRij(1,1,1,k,j)
      end do
    end do
  end do

  print*, "Rij*invRij = "
  print*, nint(B(1,:))
  print*, nint(B(2,:))
  print*, nint(B(3,:))
      
  do j = 1,3
    do  i = 1,3
      B(i,j) = 0.d0
      do k = 1,3
        B(i,j) = B(i,j) + invRij(1,1,1,i,k)*Rij(1,1,1,k,j)
      end do
    end do
  end do

  print*, "invRij*Rij = "
  print*, nint(B(1,:))
  print*, nint(B(2,:))
  print*, nint(B(3,:))

  print*, " "
  print*, "Testing RSR^-1 contraction"
  print*, "--------------------------"
  print*, " "

  A2(1,1,1,1,1) = 1.d0
  A2(1,1,1,1,2) = 4.d0
  A2(1,1,1,1,3) = 9.d0
  A2(1,1,1,2,1) = 3.d0
  A2(1,1,1,2,2) = 9.d0
  A2(1,1,1,2,3) = 7.d0
  A2(1,1,1,3,1) = 4.d0
  A2(1,1,1,3,2) = 2.d0
  A2(1,1,1,3,3) = 1.d0
! A = [1 4 9;
!      3 9 7;
!      4 2 1];
  B2(1,1,1,1) = 5.d0
  B2(1,1,1,2) = 8.d0
  B2(1,1,1,3) = 0.d0
  B2(1,1,1,4) = 3.d0
  B2(1,1,1,5) = 1.d0
  B2(1,1,1,6) = 1.d0
! B = [5 8 0;
!      8 3 1;
!      0 1 1];
  C(1,1,1,1,1) = 8.d0
  C(1,1,1,1,2) = 2.d0
  C(1,1,1,1,3) = 0.d0
  C(1,1,1,2,1) = 1.d0
  C(1,1,1,2,2) = 2.d0
  C(1,1,1,2,3) = 3.d0
  C(1,1,1,3,1) = 0.d0
  C(1,1,1,3,2) = 4.d0
  C(1,1,1,3,3) = 5.d0
!C = [8 2 0;
!     1 2 3;
!     0 4 5];
  call contract_Rik_Skl_invRlj(D,A2,C,B2,1,1)
  call assert(nint(D(1,1,1)) == 325,'D = 325')
  call contract_Rik_Skl_invRlj(D,A2,C,B2,1,2)
  call assert(nint(D(1,1,1)) == 184,'D = 184')
  call contract_Rik_Skl_invRlj(D,A2,C,B2,1,3)
  call assert(nint(D(1,1,1)) == 152,'D = 152')
  call contract_Rik_Skl_invRlj(D,A2,C,B2,2,1)
  call assert(nint(D(1,1,1)) == 754,'D = 754')
  call contract_Rik_Skl_invRlj(D,A2,C,B2,2,2)
  call assert(nint(D(1,1,1)) == 354,'D = 354')
  call contract_Rik_Skl_invRlj(D,A2,C,B2,2,3)
  call assert(nint(D(1,1,1)) == 254,'D = 254')
  call contract_Rik_Skl_invRlj(D,A2,C,B2,3,1)
  call assert(nint(D(1,1,1)) == 327,'D = 327')
  call contract_Rik_Skl_invRlj(D,A2,C,B2,3,2)
  call assert(nint(D(1,1,1)) == 162,'D = 162')
  call contract_Rik_Skl_invRlj(D,A2,C,B2,3,3)
  call assert(nint(D(1,1,1)) == 132,'D = 132')

  print*, "Test PASSED!"
! D =
!    325   184   152
!    754   354   254
!    327   162   132 
end program
