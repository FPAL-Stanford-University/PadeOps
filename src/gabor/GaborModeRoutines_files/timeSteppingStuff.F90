subroutine rk4Step(uR,uI,k,x,dt,nuMod,KE,L,nu,dudx,U)
  real(rkind), dimension(3), intent(inout) :: uR, uI, k, x
  real(rkind), intent(in) :: dt, nuMod, KE, L, nu
  real(rkind), dimension(3,3), intent(in) :: dudx
  real(rkind), dimension(3), intent(in) :: U
  real(rkind), dimension(4) :: a, b
  real(rkind), dimension(3) :: duR, duI, dk, rhsR, rhsI, rhsK, uRstar, uIStar, &
    kstar
  integer :: rk

  ! RK4 coefficients
  a = [0.d0, 0.5d0, 0.5d0, 1.d0]
  b = [1.d0, 2.d0, 2.d0, 1.d0]/6.d0

  duR = 0.d0
  duI = 0.d0
  dk = 0.d0
  rhsR = 0.d0
  rhsI = 0.d0
  rhsK = 0.d0

  do rk = 1,4
    uRstar = uR + dt*a(rk)*rhsR
    uIstar = uI + dt*a(rk)*rhsI
    kstar  = k  + dt*a(rk)*rhsK

    call getUrhs(rhsR,uRstar,dudx,kstar,nuMod,KE,L,nu)
    call getUrhs(rhsI,uIstar,dudx,kstar,nuMod,KE,L,nu)
    call getKrhs(rhsK,kstar,dudx)

    duR = duR + dt*b(rk)*rhsR
    duI = duI + dt*b(rk)*rhsI
    dk  = dk  + dt*b(rk)*rhsK
  end do
  uR = uR + duR
  uI = uI + duI
  k  = k  + dk
  x = x + dt*U
end subroutine

subroutine getUrhs(rhs,u,dudx,k,nuMod,KE,L,nu)
  real(rkind), dimension(3), intent(out) :: rhs
  real(rkind), dimension(3), intent(in) :: u
  real(rkind), dimension(3,3), intent(in) :: dudx
  real(rkind), dimension(3), intent(in) :: k
  real(rkind), intent(in) :: nuMod, KE, L, nu
  real(rkind), dimension(3,3) :: delta
  real(rkind) :: ksq, nuk
  integer :: jj, kk, ll

  rhs = 0.d0
  delta = 0.d0
  delta(1,1) = 1.d0; delta(2,2) = 1.d0; delta(3,3) = 1.d0
  ksq = sum(k*k)

  do jj = 1,3
    do kk = 1,3
      do ll = 1,3 ! Indices correspond to Pope eqn (11.83)
        rhs(jj) = rhs(jj) - u(kk)*dudx(ll,kk)*(delta(jj,ll) - 2.d0*k(jj)*k(ll)/ksq)
      end do
    end do
  end do
  
  ! Spectral eddy viscosity. Eqn (3.22) in Aditya's thesis
  nuk = sqrt(nu*nu + nuMod*KE*L*3.d0/8.d0*ksq**(-4.d0/3.d0)) - nu
 
  rhs = rhs - nuk*ksq*u
end subroutine

subroutine getKrhs(rhs,k,dudx)
  real(rkind), dimension(3), intent(out) :: rhs
  real(rkind), dimension(3), intent(in) :: k
  real(rkind), dimension(3,3), intent(in) :: dudx

  rhs(1) = -(k(1)*dudx(1,1) + k(2)*dudx(2,1) + k(3)*dudx(3,1))
  rhs(2) = -(k(1)*dudx(1,2) + k(2)*dudx(2,2) + k(3)*dudx(3,2))
  rhs(3) = -(k(1)*dudx(1,3) + k(2)*dudx(2,3) + k(3)*dudx(3,3))
end subroutine
    
subroutine getDtMax(dt,uhatR,uhatI,S,kmag)
  ! Compute the maximum stable time step size
  ! Inputs:
  !     uhatR, uhatI --> vector of complex velocity amplitudes
  !     S --> Frobenius norm of large scale velocity gradient at GM location
  ! Outputs:
  !     dt --> Maximum stable time step based on time-scales of the problem

  real(rkind), dimension(3) :: uhatR, uhatI
  real(rkind), intent(in) :: S, kmag
  real(rkind), intent(out) :: dt
  real(rkind) :: umax, tau1, tau2

  umax = sqrt(uhatR(1)*uhatR(1) + uhatI(1)*uhatI(1) + &
              uhatR(2)*uhatR(2) + uhatI(2)*uhatI(2) + &
              uhatR(3)*uhatR(3) + uhatI(3)*uhatI(3))

  call assert(kmag > small,'kmag > small -- getDtMax()') 
  if (umax > small) then
    tau1 = 1.d0/(kmag*umax)
  else
    tau1 = big
  end if

  if (S > small) then 
    tau2 = 1.d0/S
  else
    tau2 = big
  end if

  dt = minval([tau1,tau2])/3.d0
end subroutine
