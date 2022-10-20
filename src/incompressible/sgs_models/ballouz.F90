subroutine init_ballouz(this,alpha)
   use fortran_assert, only: assert
   class(sgs_igrid), intent(inout) :: this
   real(rkind), intent(in) :: alpha

   call assert(this%isPeriodic,'Ballouz SGS model should only be used'//&
     'for triply periodic flows')
   this%useCglobal = .true. 
   this%cmodel_global = alpha
   this%isEddyViscosityModel = .false.
   call message(1,"Ballouz model initialized")
end subroutine

subroutine destroy_ballouz(this)
   class(sgs_igrid), intent(inout) :: this

   this%isEddyViscosityModel = .false.
end subroutine

subroutine compute_tauij_ballouz(this)
   use random,         only: uniform_random 
   use fortran_assert, only: assert
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)) :: &
     phiE, thetaE, psiE
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)) :: &
     phiC, thetaC, psiC
   !real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)) :: &
   !  R11E, R12E, R13E, R21E, R22E, R23E
   !real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)) :: &
   !  R11C, R12C, R13C, R21C, R22C, R23C, R31C, R32C, R33C
   real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),3,3)&
     :: RijE, invRijE
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),3,3)&
     :: RijC, invRijC
   real(rkind) :: phi, theta, psi
   
   ! Randomly sample Euler angles
   call uniform_random(phi,  0.d0,2.d0*pi)
   call uniform_random(theta,0.d0,2.d0*pi)
   call uniform_random(psi,  0.d0,2.d0*pi)

   phiC   = phi;   phiE   = phi
   thetaC = theta; thetaE = theta
   psiC   = psi;   psiE   = psi
   
   ! Get the rotation tensor corresponding to the Euler angles
   call getRotationTensor(RijE,phiE,thetaE,psiE)
   call getRotationTensor(RijC,phiC,thetaC,psiC)
   call getInverseRotationTensor(invRijE,phiE,thetaE,psiE)
   call getInverseRotationTensor(invRijC,phiC,thetaC,psiC)

   ! Do matrix multiplication
   call contract_Rik_Skl_invRlj(this%tau_11,RijC,invRijC,this%S_ij_C,1,1)
   call contract_Rik_Skl_invRlj(this%tau_12,RijC,invRijC,this%S_ij_C,1,2)
   call contract_Rik_Skl_invRlj(this%tau_13,RijE,invRijE,this%S_ij_E,1,3)
   call contract_Rik_Skl_invRlj(this%tau_22,RijC,invRijC,this%S_ij_C,2,2)
   call contract_Rik_Skl_invRlj(this%tau_23,RijE,invRijE,this%S_ij_E,2,3)
   call contract_Rik_Skl_invRlj(this%tau_33,RijC,invRijC,this%S_ij_C,3,3)

   ! Remove trace
   call removeTrace(this%tau_11,this%tau_22,this%tau_33)

   ! Multiply by model constant
   call assert(this%useCglobal,"Ballouz model only supports a global"//&
    " model constant -- ballouz.F90")
   this%tau_11 = -this%tau_11*this%cmodel_global
   this%tau_12 = -this%tau_12*this%cmodel_global
   this%tau_13 = -this%tau_13*this%cmodel_global
   this%tau_22 = -this%tau_22*this%cmodel_global
   this%tau_23 = -this%tau_23*this%cmodel_global
   this%tau_33 = -this%tau_33*this%cmodel_global

end subroutine

subroutine removeTrace(tau11,tau22,tau33)
  real(rkind), dimension(:,:,:), intent(inout) :: tau11, tau22, tau33
  real(rkind), dimension(size(tau11,1),size(tau11,2),size(tau11,3)) :: taukk

  taukk = tau11 + tau22 + tau33

  tau11 = tau11 - taukk/3.d0
  tau22 = tau22 - taukk/3.d0
  tau33 = tau33 - taukk/3.d0
end subroutine
  
subroutine contract_Rik_Skl_invRlj(tauij,Rij,invRij,Sij,i,j)
  ! This computes the following triple product:
  ! tau_ij = (R_ik)(S_kl)(invR_lj)
  ! The key to interpreting the code below is noting that 
  ! Sij is symmetric and size (nx,ny,nz,6)
  real(rkind), dimension(:,:,:), intent(inout) :: tauij
  real(rkind), dimension(:,:,:,:,:), intent(in) :: Rij, invRij
  real(rkind), dimension(:,:,:,:), intent(in) :: Sij
  integer, intent(in) :: i, j

  tauij = Rij(:,:,:,i,1)*Sij(:,:,:,1)*invRij(:,:,:,1,j) + &
          Rij(:,:,:,i,2)*Sij(:,:,:,2)*invRij(:,:,:,1,j) + &
          Rij(:,:,:,i,3)*Sij(:,:,:,3)*invRij(:,:,:,1,j) + &
          Rij(:,:,:,i,1)*Sij(:,:,:,2)*invRij(:,:,:,2,j) + &
          Rij(:,:,:,i,2)*Sij(:,:,:,4)*invRij(:,:,:,2,j) + &
          Rij(:,:,:,i,3)*Sij(:,:,:,5)*invRij(:,:,:,2,j) + &
          Rij(:,:,:,i,1)*Sij(:,:,:,3)*invRij(:,:,:,3,j) + &
          Rij(:,:,:,i,2)*Sij(:,:,:,5)*invRij(:,:,:,3,j) + &
          Rij(:,:,:,i,3)*Sij(:,:,:,6)*invRij(:,:,:,3,j)
end subroutine

subroutine getRotationTensor(Rij,phi,theta,psi)
  use fortran_assert, only: assert
  real(rkind), dimension(:,:,:,:,:), intent(out) :: Rij
  real(rkind), dimension(:,:,:), intent(in) :: phi, theta, psi
  integer, dimension(3,3) :: idx
  integer :: i, j

  idx(1,1) = 11; idx(1,2) = 12; idx(1,3) = 13
  idx(2,1) = 21; idx(2,2) = 22; idx(2,3) = 23
  idx(3,1) = 31; idx(3,2) = 32; idx(3,3) = 33

  call assert(size(Rij,4) == 3)
  call assert(size(Rij,5) == 3)
  do j = 1,3
    do i = 1,3
      call getRotationTensorComponent(Rij(:,:,:,i,j),phi,theta,psi,idx(i,j))
    end do
  end do
end subroutine

subroutine getInverseRotationTensor(invRij,phi,theta,psi)
  use fortran_assert, only: assert
  real(rkind), dimension(:,:,:,:,:), intent(out) :: invRij
  real(rkind), dimension(:,:,:), intent(in) :: phi, theta, psi
  integer, dimension(3,3) :: idx
  integer :: i, j

  idx(1,1) = 11; idx(1,2) = 12; idx(1,3) = 13
  idx(2,1) = 21; idx(2,2) = 22; idx(2,3) = 23
  idx(3,1) = 31; idx(3,2) = 32; idx(3,3) = 33

  call assert(size(invRij,4) == 3)
  call assert(size(invRij,5) == 3)
  do j = 1,3
    do i = 1,3
      call getInverseRotationTensorComponent(invRij(:,:,:,i,j),phi,theta,psi,idx(i,j))
    end do
  end do
end subroutine

subroutine getRotationTensorComponent(Rij,phi,theta,psi,comp)
  real(rkind), dimension(:,:,:), intent(out) :: Rij
  real(rkind), dimension(:,:,:), intent(in) :: phi, theta, psi
  integer, intent(in) :: comp

  ! See Ch.13 of "Classical Dynamics of Particles and Systems" by 
  ! Stephen Thornton & Jerry Marion
  select case(comp)
  case(11)
    Rij = cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi)
  case(21)
    Rij = -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi) 
  case(31)
    Rij = sin(theta)*sin(phi)
  case(12)
    Rij = cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi)
  case(22)
    Rij = -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi)
  case(32)
    Rij = -sin(theta)*cos(phi)
  case(13)
    Rij = sin(psi)*sin(theta)
  case(23)
    Rij = cos(psi)*sin(theta)
  case(33)
    Rij = cos(theta)
  end select
end subroutine

subroutine getInverseRotationTensorComponent(invRij,phi,theta,psi,comp)
  real(rkind), dimension(:,:,:), intent(out) :: invRij
  real(rkind), dimension(:,:,:), intent(in) :: phi, theta, psi
  integer, intent(in) :: comp
  select case(comp)
  case(11)
    invRij = (cos(phi)*cos(psi)*cos(theta)**2.d0 + cos(phi)*cos(psi)*sin(theta)**2.d0 - &
    cos(theta)*sin(phi)*sin(psi))/(cos(phi)**2.d0*cos(psi)**2.d0*sin(theta)**2.d0 + &
    cos(phi)**2.d0*cos(theta)**2.d0*sin(psi)**2.d0 + &
    cos(psi)**2.d0*cos(theta)**2.d0*sin(phi)**2.d0 + &
    cos(phi)**2.d0*sin(psi)**2.d0*sin(theta)**2.d0 + &
    cos(psi)**2.d0*sin(phi)**2.d0*sin(theta)**2.d0 + &
    cos(theta)**2.d0*sin(phi)**2.d0*sin(psi)**2.d0 + &
    sin(phi)**2.d0*sin(psi)**2.d0*sin(theta)**2.d0 + &
    cos(phi)**2.d0*cos(psi)**2.d0*cos(theta)**2.d0)
  case(21)
    invRij = (cos(psi)*cos(theta)**2.d0*sin(phi) + cos(psi)*sin(phi)*sin(theta)**2.d0 + &
    cos(phi)*cos(theta)*sin(psi))/(cos(phi)**2.d0*cos(psi)**2.d0*sin(theta)**2.d0 + &
    cos(phi)**2.d0*cos(theta)**2.d0*sin(psi)**2.d0 + &
    cos(psi)**2.d0*cos(theta)**2.d0*sin(phi)**2.d0 + &
    cos(phi)**2.d0*sin(psi)**2.d0*sin(theta)**2.d0 + &
    cos(psi)**2.d0*sin(phi)**2.d0*sin(theta)**2.d0 + &
    cos(theta)**2.d0*sin(phi)**2.d0*sin(psi)**2.d0 + &
    sin(phi)**2.d0*sin(psi)**2.d0*sin(theta)**2.d0 + &
    cos(phi)**2.d0*cos(psi)**2.d0*cos(theta)**2.d0)
  case(31)
    invRij = (sin(psi)*sin(theta))/(cos(psi)**2.d0*cos(theta)**2.d0 + &
    cos(psi)**2.d0*sin(theta)**2.d0 + &
    cos(theta)**2.d0*sin(psi)**2.d0 + sin(psi)**2.d0*sin(theta)**2.d0)
  case(12)
    invRij = -(cos(phi)*cos(theta)**2.d0*sin(psi) + cos(phi)*sin(psi)*sin(theta)**2.d0 + &
    cos(psi)*cos(theta)*sin(phi))/(cos(phi)**2.d0*cos(psi)**2.d0*sin(theta)**2.d0 + &
    cos(phi)**2.d0*cos(theta)**2.d0*sin(psi)**2.d0 + &
    cos(psi)**2.d0*cos(theta)**2.d0*sin(phi)**2.d0 + &
    cos(phi)**2.d0*sin(psi)**2.d0*sin(theta)**2.d0 + &
    cos(psi)**2.d0*sin(phi)**2.d0*sin(theta)**2.d0 + &
    cos(theta)**2.d0*sin(phi)**2.d0*sin(psi)**2.d0 + &
    sin(phi)**2.d0*sin(psi)**2.d0*sin(theta)**2.d0 + &
    cos(phi)**2.d0*cos(psi)**2.d0*cos(theta)**2.d0)
  case(22)
    invRij = -(cos(theta)**2.d0*sin(phi)*sin(psi) + sin(phi)*sin(psi)*sin(theta)**2.d0 - &
    cos(phi)*cos(psi)*cos(theta))/(cos(phi)**2.d0*cos(psi)**2.d0*sin(theta)**2.d0 + &
    cos(phi)**2.d0*cos(theta)**2.d0*sin(psi)**2.d0 + &
    cos(psi)**2.d0*cos(theta)**2.d0*sin(phi)**2.d0 + &
    cos(phi)**2.d0*sin(psi)**2.d0*sin(theta)**2.d0 + &
    cos(psi)**2.d0*sin(phi)**2.d0*sin(theta)**2.d0 + &
    cos(theta)**2.d0*sin(phi)**2.d0*sin(psi)**2.d0 + &
    sin(phi)**2.d0*sin(psi)**2.d0*sin(theta)**2.d0 + &
    cos(phi)**2.d0*cos(psi)**2.d0*cos(theta)**2.d0)
  case(32)
    invRij = (cos(psi)*sin(theta))/(cos(psi)**2.d0*cos(theta)**2.d0 + &
    cos(psi)**2.d0*sin(theta)**2.d0 + &
    cos(theta)**2.d0*sin(psi)**2.d0 + sin(psi)**2.d0*sin(theta)**2.d0)
  case(13)
    invRij = (sin(phi)*sin(theta))/(cos(phi)**2.d0*cos(theta)**2.d0 + &
    cos(phi)**2.d0*sin(theta)**2.d0 + &
    cos(theta)**2.d0*sin(phi)**2.d0 + sin(phi)**2.d0*sin(theta)**2.d0)
  case(23)
    invRij = -(cos(phi)*sin(theta))/(cos(phi)**2.d0*cos(theta)**2.d0 + &
    cos(phi)**2.d0*sin(theta)**2.d0 + &
    cos(theta)**2.d0*sin(phi)**2.d0 + sin(phi)**2.d0*sin(theta)**2.d0)
  case(33)
    invRij = cos(theta)/(cos(theta)**2.d0 + sin(theta)**2.d0)
  end select
end subroutine
