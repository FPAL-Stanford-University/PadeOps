  subroutine RK4(this, AnuRNG, nu, epsKE, dudx, dudy, dudz, dvdx, dvdy, dvdz, &
      dwdx, dwdy, dwdz, kmin, dt)
      class(gaborMode), intent(inout), target :: this
      real(rkind), intent(in) :: AnuRNG, nu, epsKE
      real(rkind), intent(in) :: dudx, dudy, dudz
      real(rkind), intent(in) :: dvdx, dvdy, dvdz
      real(rkind), intent(in) :: dwdx, dwdy, dwdz
      real(rkind), intent(in) :: kmin, dt
      real(rkind), dimension(3) :: uR, uRin, uRrhs1, uRrhs2, uRrhs3, uRrhs4
      real(rkind), dimension(3) :: uI, uIin, uIrhs1, uIrhs2, uIrhs3, uIrhs4
      real(rkind), dimension(3) :: k, kin, krhs1, krhs2, krhs3, krhs4
      real(rkind), dimension(3,3) :: dUidXj
            
      k(1)        = this%kx;    k(2)        = this%ky;    k(3)        = this%kz
      uR(1)       = this%uhatR; uR(2)       = this%vhatR; uR(3)       = this%whatR
      uI(1)       = this%uhatI; uI(2)       = this%vhatI; uI(3)       = this%whatI
      dUidXj(1,1) = dudx;       dUidXj(1,2) = dudy;       dUidXj(1,3) = dudz
      dUidXj(2,1) = dvdx;       dUidXj(2,2) = dvdy;       dUidXj(2,3) = dvdz
      dUidXj(3,1) = dwdx;       dUidXj(3,2) = dwdy;       dUidXj(3,3) = dwdz
      
      kin = k
      uRin = uR
      uIin = uI
      call gaborEvoRHS(epsKE, AnuRNG, nu, kin, uRin, uIin, dUidXj, &
           uRrhs1, uIrhs1, krhs1, kmin) 
      
      kin = krhs1*dt*0.5D0 + k
      uRin = uRrhs1*dt*0.5D0 + uR
      uIin = uIrhs1*dt*0.5D0 + uI
      call gaborEvoRHS(epsKE, AnuRNG, nu, kin, uRin, uIin, dUidXj, &
           uRrhs2, uIrhs2, krhs2, kmin) 
      
      kin = krhs2*dt*0.5D0 + k
      uRin = uRrhs2*dt*0.5D0 + uR
      uIin = uIrhs2*dt*0.5D0 + uI
      call gaborEvoRHS(epsKE, AnuRNG, nu, kin, uRin, uIin, dUidXj, &
           uRrhs3, uIrhs3, krhs3, kmin) 
      
      kin = krhs3*dt + k
      uRin = uRrhs3*dt + uR
      uIin = uIrhs3*dt + uI
      call gaborEvoRHS(epsKE, AnuRNG, nu, kin, uRin, uIin, dUidXj, &
           uRrhs4, uIrhs4, krhs4, kmin) 
  
      this%uhatR = this%uhatR + (dt/6.D0)*(uRrhs1(1) + 2.D0*(uRrhs2(1) + uRrhs3(1)) + uRrhs4(1))
      this%uhatI = this%uhatI + (dt/6.D0)*(uIrhs1(1) + 2.D0*(uIrhs2(1) + uIrhs3(1)) + uIrhs4(1))
      this%vhatR = this%vhatR + (dt/6.D0)*(uRrhs1(2) + 2.D0*(uRrhs2(2) + uRrhs3(2)) + uRrhs4(2))
      this%vhatI = this%vhatI + (dt/6.D0)*(uIrhs1(2) + 2.D0*(uIrhs2(2) + uIrhs3(2)) + uIrhs4(2))
      this%whatR = this%whatR + (dt/6.D0)*(uRrhs1(3) + 2.D0*(uRrhs2(3) + uRrhs3(3)) + uRrhs4(3))
      this%whatI = this%whatI + (dt/6.D0)*(uIrhs1(3) + 2.D0*(uIrhs2(3) + uIrhs3(3)) + uIrhs4(3))

      this%kx = this%kx + (dt/6.D0)*(krhs1(1) + 2.D0*(krhs2(1) + krhs3(1)) + krhs4(1))
      this%ky = this%ky + (dt/6.D0)*(krhs1(2) + 2.D0*(krhs2(2) + krhs3(2)) + krhs4(2))
      this%kz = this%kz + (dt/6.D0)*(krhs1(3) + 2.D0*(krhs2(3) + krhs3(3)) + krhs4(3))
  end subroutine

  pure subroutine gaborEvoRHS(epsKE, AnuRNG, nuMolec, k, uhatR, uhatI, dudx, uRrhs, uIrhs, krhs, kmin)
      real(rkind), dimension(3), intent(in) :: k
      real(rkind), dimension(3), intent(in) :: uhatR, uhatI
      real(rkind), dimension(3), intent(out) :: uRrhs, uIrhs
      real(rkind), dimension(3), intent(out) :: krhs
      real(rkind), dimension(3,3), intent(in) :: dudx
      real(rkind), dimension(3,3) :: kronD
      real(rkind) :: Ksq, Onebyksq, Kmag
      real(rkind), intent(in)  :: epsKE, kmin, AnuRNG, nuMolec
      integer :: jj, kk, ll
      real(rkind) :: nuK, nuFact
      
      kronD = 0.d0 
      kronD(1,1) = 1.d0; kronD(2,2) = 1.d0; kronD(3,3) = 1.d0

      uRrhs = 0.d0 
      uIrhs = 0.d0 
      krhs = 0.d0 
      Ksq = (k(1)*k(1) + k(2)*k(2) + k(3)*k(3) + 1.d-14)
      Onebyksq = 1.0d0/Ksq
      Kmag = sqrt(Ksq)
      nuFact = 1.d0
      if (Kmag < kmin) nuFact = 10.d0

      do jj = 1,3
          do kk = 1,3
              do ll = 1,3
                  uRrhs(jj) = uRrhs(jj) + (2.d0*k(jj)*k(ll)*OneByksq - &
                              kronD(jj,ll)) * uhatR(kk)*dudx(ll,kk)
                  uIrhs(jj) = uIrhs(jj) + (2.d0*k(jj)*k(ll)*OneByksq - &
                              kronD(jj,ll)) * uhatI(kk)*dudx(ll,kk)
              end do 
          end do 
      end do 
      
      nuK = AnuRNG*epsKE**(1.0/3.0)*(Ksq**(-2.0/3.0))*nuFact
      
      uRrhs(1) = uRrhs(1) - (nuK + nuMolec)*Ksq*uhatR(1) 
      uIrhs(1) = uIrhs(1) - (nuK + nuMolec)*Ksq*uhatI(1)
      uRrhs(2) = uRrhs(2) - (nuK + nuMolec)*Ksq*uhatR(2)
      uIrhs(2) = uIrhs(2) - (nuK + nuMolec)*Ksq*uhatI(2)
      uRrhs(3) = uRrhs(3) - (nuK + nuMolec)*Ksq*uhatR(3)
      uIrhs(3) = uIrhs(3) - (nuK + nuMolec)*Ksq*uhatI(3)

      do ll = 1,3
          do jj = 1,3
              krhs(ll) = krhs(ll) + (-k(jj)*dudx(jj,ll))
          end do 
      end do 
  end subroutine

  subroutine periodicBClocation(this)
    use domainSetup, only: xDom, yDom, zDom
    class(GaborMode), intent(inout), target :: this
    
    if (this%x < xDom(1)) then 
      this%x = xDom(2) - (xDom(1) - this%x)
    elseif (this%x > xDom(2)) then
      this%x = xDom(1) + (this%x - xDom(2))
    end if
    
    if (this%y < yDom(1)) then 
      this%y = yDom(2) - (yDom(1) - this%y)
    elseif (this%y > yDom(2)) then
      this%y = yDom(1) + (this%y - yDom(2))
    end if
    
    if (this%z < zDom(1)) then 
      this%z = zDom(2) - (zDom(1) - this%z)
    elseif (this%z > zDom(2)) then
      this%z = zDom(1) + (this%z - zDom(2))
    end if
  end subroutine 
