subroutine strainIsotropicModes()
  real(rkind), dimension(nmodes) :: kabs, input
  real(rkind) :: output, tauEddy, Sgm, Lgm, KEgm, UpGM, VpGM, WpGM
  real(rkind), dimension(3,3) :: dudx
  integer :: n, tid
  real(rkind) :: dt
  real(rkind), dimension(3) :: k, uRtmp, uItmp
  character(len=clen) :: mssg

  kabs = sqrt(kx*kx + ky*ky + kz*kz)
  input = -(kabs)**(-2.0_rkind)

  !!$OMP PARALLEL SHARED(input,kx,ky,kz,uhatR,uhatI,vhatR,vhatI,whatR,whatI) &
  !!$OMP SHARED(Anu,numolec,S,ctau,kabs)
  !!$OMP PRIVATE(output,n,k,uRtmp,uItmp,tid,tauEddy,dt,KEgm,Lgm,dudx,Sgm)
  !!$OMP DO
  do n = 1,nmodes
    CALL PFQ(input(n),output)
    call getLargeScaleDataAtModeLocation(n,Sgm,dudx,Lgm,KEgm,UpGM,VpGM,WpGM)
    call getDtMax(dt,n,Sgm,kabs(n),UpGM,VpGM,WpGM)
    
    k = [kx(n), ky(n), kz(n)]
    uRtmp = [uhatR(n), vhatR(n), whatR(n)]
    uItmp = [uhatI(n), vhatI(n), whatI(n)]

    tauEddy = ctau/Sgm*((kabs(n))**(-2.0_rkind/3.0_rkind))/sqrt(output) 
    
    do tid = 1,nint(tauEddy/dt)
      call rk4Step(uRtmp,uItmp,k,dt,Anu,KEgm,Lgm,numolec,dudx)
    end do
    kx(n) = k(1)
    ky(n) = k(2)
    kz(n) = k(3)
    uhatR(n) = uRtmp(1)
    uhatI(n) = uItmp(1)
    vhatR(n) = uRtmp(2)
    vhatI(n) = uItmp(2)
    whatR(n) = uRtmp(3)
    whatI(n) = uItmp(3)
    
    if (mod(n,10000) == 0 .and. nrank == 0) then
      write(mssg,'(F8.5,A)') real(n,rkind)/real(nmodes,rkind)*100.d0,'% complete'
      print*, trim(mssg)
    end if
  end do
  !!$OMP END DO
  !!$OMP END PARALLEL
end subroutine

subroutine getLargeScaleDataAtModeLocation(gmID,Sgm,dudx,Lgm,KEgm,&
    UpGM,VpGM,WpGM)
  ! Get large scale data interpolated to Gabor mode location
  ! Inputs:
  !     S --> Frobenius norm of velocity gradient tensor (i.e. dUidxj*dUidxj)
  !     gmID --> Gabor mode ID
  ! Outputs:
  !     Sgm --> S at GM location
  !     dudx --> Velocity gradient at GM location
  !     Lgm --> Integral length scale computed at QH region center
  !     KEgm --> Large scale kinetic energy computed at QH region center
  !     UpGM, etc. --> Large scale fluctuating velocity at GM location

  integer, intent(in) :: gmID
  real(rkind), intent(out) :: Sgm, Lgm, KEgm, UpGM, VpGM, WpGM
  real(rkind), dimension(3,3), intent(out) :: dudx
  integer :: i, j, QHx, QHy, QHz 

  ! Interpolate large scale velocity/velocity-gradient to mode location 
  call interpToMode(Uph,UpGM,gmID,gpLESe,dxLES,dyLES,dzLES,xLESe,yLESe,zLESe)
  call interpToMode(Vph,VpGM,gmID,gpLESe,dxLES,dyLES,dzLES,xLESe,yLESe,zLESe)
  call interpToMode(Wph,WpGM,gmID,gpLESe,dxLES,dyLES,dzLES,xLESe,yLESe,zLESe)

  call interpToMode(Sh,Sgm,gmID,gpLESe,dxLES,dyLES,dzLES,xLESe,yLESe,zLESe)
  do j = 1,3
    do i = 1,3
      call interpToMode(gradUh(i,j,:,:,:),dudx(i,j),gmID,gpLESe,dxLES,dyLES,&
        dzLES,xLESe,yLESe,zLESe)
    end do
  end do
  
  ! Find the index of the mode's QH region
  QHx = findMeshIdx(gmxloc(gmID),xQHedge)
  QHy = findMeshIdx(gmyloc(gmID),yQHedge)
  QHz = findMeshIdx(gmzloc(gmID),zQHedge)
  
  ! Use L and KE for the QH region that the mode resides in
  Lgm = L(QHx,QHy,QHz)
  KEgm = KE(QHx,QHy,QHz)
 
end subroutine

integer function findMeshIdx(xloc,x)
  ! Returns the index of cell x where xloc is found

  real(rkind), intent(in) :: xloc
  real(rkind), dimension(:), intent(in) :: x
  integer :: i, N

  N = size(x)-1
  findMeshIdx = 0
  do i = 1,N
    if (xloc > x(i) .and. xloc <= x(i+1)) then
      findMeshIdx = i
    end if
  end do
end function

subroutine getDtMax(dt,gmID,S,kmag,UpGM,VpGM,WpGM)
  ! Compute the maximum stable time step size
  ! Inputs:
  !     gmID --> Gabor mode (GM) ID
  !     S --> Frobenius norm of large scale velocity gradient at GM location
  !     UpGM, etc. --> Fluctuating large scale velocity at GM location
  ! Outputs:
  !     dt --> Maximum stable time step based on time-scales of the problem

  real(rkind), intent(in) :: S, kmag, UpGM, VpGM, WpGM
  real(rkind), intent(out) :: dt
  integer, intent(in) :: gmID
  real(rkind) :: umax, tau1, tau2, tau3, UlrgMax, tol = 1.d-13

  umax = sqrt(uhatR(gmID)*uhatR(gmID) + uhatI(gmID)*uhatI(gmID) + &
              vhatR(gmID)*vhatR(gmID) + vhatI(gmID)*vhatI(gmID) + &
              whatR(gmID)*whatR(gmID) + whatI(gmID)*whatI(gmID))
  
  UlrgMax =  sqrt(UpGM*UpGM + VpGM*VpGM + WpGM*WpGM)
  tau1 = 1.d0/(kmag*umax)
  if (UlrgMax < tol) then
    tau2 = 100.d0
  else
    tau2 = 1.d0/(kmag*UlrgMax)
  end if
  tau3 = 1.d0/S

  dt = minval([tau1,tau2,tau3])/3.d0
end subroutine

SUBROUTINE PFQ(X_ACT,Z_ACT)
  use hygfx_mod, only: HYGFX
  REAL(kind=8) :: A, B, C, X, HF
  REAL(kind=8) :: A_ACT, B_ACT, C_ACT
  REAL(rkind), intent(in) :: X_ACT
  real(rkind), intent(out) :: Z_ACT
  A_ACT = 0.333333d0
  B_ACT = 2.833333d0
  C_ACT = 1.333333d0
  A = A_ACT
  B = B_ACT
  C = C_ACT
  X = real(X_ACT,kind=8)
  
  CALL HYGFX(A,B,C,X,HF)
  Z_ACT = real(HF,rkind)
  RETURN
END SUBROUTINE PFQ 
