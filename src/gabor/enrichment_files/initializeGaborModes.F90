subroutine generateIsotropicModes(this)
  use gridtools, only: logspace
  use random,    only: uniform_random
  use GaborModeRoutines, only: normalizeVec, isOrthogonal, getModelSpectrum, &
                               doWarning, small
  use mpi
  class(enrichmentOperator), intent(inout) :: this
  character(len=clen) :: mssg1, mssg2, mssg3, mssg4
  real(rkind) :: xmin, ymin, zmin
  integer :: i, j, k, kid, thetaID, nk, ntheta, nmodes
  integer :: seed1=1, seed2=2, seed3=3, seed4=4, seed5=5, seed6=6, seed7=7
  integer :: isz, jsz, ksz
  real(rkind), dimension(:),         allocatable :: kedge, kmag, dk, rand1, &
                                                    theta, kztemp, r, E
  real(rkind), dimension(:,:), allocatable :: orientationX, orientationY, &
                                              orientationZ, p1x, p1y, p1z, &
                                              p2x, p2y, p2z, rand2, uRmag, uImag, &
                                              umag, thetaVel
  
  ! Physical location of Gabor modes
  real(rkind), dimension(:,:,:,:), allocatable :: gmx, gmy, gmz

  ! Wave-vector components
  real(rkind), dimension(:,:,:,:,:), allocatable :: k1, k2, k3

  ! Velocity-vector amplitudes
  real(rkind), dimension(:,:,:,:,:), allocatable :: uR, uI, vR, vI, wR, wI
  
  ! QHmesh indices of each mode
  real(rkind), dimension(:,:,:,:,:), allocatable :: QHi, QHj, QHk

  ! Misc
  integer :: ierr

  call message("                                                                ")
  call message("================================================================")
  call message("============== Initializing isotropic Gabor modes ==============")
  write(mssg1,'(A)') 'Total modes to initialize:'
  write(mssg2,'(A,I2,A,I2,A)') '(', this%nk,' shells per QH region) X (',&
    this%ntheta, ' modes per shell) X ...'
  write(mssg3,'(A,I2,A,I2,A,I2,A)') &
    '(', this%QHgrid%nx,' QH regions in x) X (', &
    this%QHgrid%ny,' QH regions in y) X (', this%QHgrid%nz, ' QH regions in z) = ...'
  write(mssg4,'(I8,A)') this%nmodesGlobal, ' modes'
  call message(                          trim(mssg1)                             )
  call message(                          trim(mssg2)                             )
  call message(                          trim(mssg3)                             )
  call message(                          trim(mssg4)                             )
  call message("                                                                ")
  
  ! Allocate memory
  nk = this%nk
  ntheta = this%ntheta 
  isz = this%QHgrid%gpC%xsz(1) 
  jsz = this%QHgrid%gpC%xsz(2)
  ksz = this%QHgrid%gpC%xsz(3)

  allocate(kedge(nk+1))
  allocate(kmag(nk),dk(nk),E(nk))
  allocate(rand1(nk*ntheta))
  allocate(theta(ntheta),kztemp(ntheta),r(ntheta))
  allocate(umag(nk,ntheta), thetaVel(nk,ntheta), &
           uRmag(nk,ntheta), uImag(nk,ntheta),&
           rand2(nk,ntheta), p1x(nk,ntheta), &
           p1y(nk,ntheta), p1z(nk,ntheta), p2x(nk,ntheta), p2y(nk,ntheta), &
           p2z(nk,ntheta), orientationX(nk,ntheta), &
           orientationY(nk,ntheta), orientationZ(nk,ntheta))

  ! Allocate memory for ND attribute arrays
  allocate(gmx(isz,jsz,ksz,nk*ntheta), gmy(isz,jsz,ksz,nk*ntheta), &
           gmz(isz,jsz,ksz,nk*ntheta))
  allocate(uR(isz,jsz,ksz,nk,ntheta),uI(isz,jsz,ksz,nk,ntheta),&
           vR(isz,jsz,ksz,nk,ntheta),vI(isz,jsz,ksz,nk,ntheta),&
           wR(isz,jsz,ksz,nk,ntheta),wI(isz,jsz,ksz,nk,ntheta))
  allocate(k1(isz,jsz,ksz,nk,ntheta),k2(isz,jsz,ksz,nk,ntheta),&
           k3(isz,jsz,ksz,nk,ntheta))
  allocate(QHi(isz,jsz,ksz,nk,ntheta), QHj(isz,jsz,ksz,nk,ntheta), &
           QHk(isz,jsz,ksz,nk,ntheta))

  
    ! Assign wave-vector magnitudes based on logarithmically spaced shells
  ! (non-dimensionalized with L)
  !kedge = logspace(kmin,kmax,nk+1)
  kedge = logspace(log10(this%kmin),log10(this%kmax),nk+1)
  kmag = 0.5d0*(kedge(1:nk) + kedge(2:nk+1))
  dk = kedge(2:nk+1) - kedge(1:nk)
  
  do k = 1,this%QHgrid%gpC%xsz(3)
    zmin = this%QHgrid%zE(k)
    do j = 1,this%QHgrid%gpC%xsz(2)
      ymin = this%QHgrid%yE(j)
      do i = 1,this%QHgrid%gpC%xsz(1)
        xmin = this%QHgrid%xE(i)

        call this%updateSeeds(seed1,seed2,seed3,seed4,seed5,seed6,seed7,&
          this%QHgrid%gpC%xst(1)+i-1,&
          this%QHgrid%gpC%xst(2)+j-1,&
          this%QHgrid%gpC%xst(3)+k-1)

        ! Uniformily distribute modes in QH region
        call uniform_random(rand1,0.d0,1.d0,seed1)
        gmx(i,j,k,:) = xmin + this%QHgrid%dx*rand1
        
        call uniform_random(rand1,0.d0,1.d0,seed2)
        gmy(i,j,k,:) = ymin + this%QHgrid%dy*rand1

        call uniform_random(rand1,0.d0,1.d0,seed3)
        gmz(i,j,k,:) = zmin + this%QHgrid%dz*rand1
        
        do kid = 1,nk
          ! Isotropically sample wave-vector components
          call uniform_random(theta,0.d0,2.d0*pi,seed4)
          call uniform_random(kztemp,-1.d0,1.d0,seed5)

          ! Assign wave-vector components
          k3(i,j,k,kid,:) = kmag(kid)*kztemp
          r = sqrt(1.d0 - kztemp*kztemp)
          k1(i,j,k,kid,:) = kmag(kid)*r*cos(theta)
          k2(i,j,k,kid,:) = kmag(kid)*r*sin(theta)
        end do

        ! Non-dimensional model energy spectrum
        E = getModelSpectrum(kmag,this%QHgrid%KE(i,j,k),&
          this%QHgrid%L(i,j,k),nk)

        ! Amplitude of each mode such that the sum of mode energies gives
        ! the correct kinetic energy
        umag(:,1) = sqrt(2*E*dk/ntheta)
        do thetaID = 2,ntheta
          umag(:,thetaID) = umag(:,1)
        end do
  
        ! Get velocity vector orientations
        ! Generate basis of tangent plane to wavenumber shell
        p1x = 0.d0
        p1y = -k3(i,j,k,:,:)
        p1z = k2(i,j,k,:,:)

        p2x = k2(i,j,k,:,:)**2 + k3(i,j,k,:,:)**2
        p2y = -k1(i,j,k,:,:)*k2(i,j,k,:,:)
        p2z = -k1(i,j,k,:,:)*k3(i,j,k,:,:)

        call normalizeVec(p1x,p1y,p1z)
        call normalizeVec(p2x,p2y,p2z)

        ! Assign orientation of the velocity vectors
        call uniform_random(thetaVel,0.d0,2.d0*pi,seed6)
        orientationX = cos(thetaVel)*p1x + sin(thetaVel)*p2x
        orientationY = cos(thetaVel)*p1y + sin(thetaVel)*p2y
        orientationZ = cos(thetaVel)*p1z + sin(thetaVel)*p2z

        ! We have the modulus and orientation of the complex-valued amplitudes.
        ! Now assigng the real and imaginary components
        call uniform_random(rand2,0.d0,1.d0/sqrt(3.d0),seed7)
        uRmag = rand2*umag
        uImag = sqrt(umag*umag/3.d0 - uRmag*uRmag)

        uR(i,j,k,:,:) = uRmag*orientationX
        uI(i,j,k,:,:) = uImag*orientationX

        vR(i,j,k,:,:) = uRmag*orientationY
        vI(i,j,k,:,:) = uImag*orientationY

        wR(i,j,k,:,:) = uRmag*orientationZ
        wI(i,j,k,:,:) = uImag*orientationZ

        ! Mark the QH region in which these modes were initialized
        QHi(i,j,k,:,:) = i
        QHj(i,j,k,:,:) = j
        QHk(i,j,k,:,:) = k
      end do
    end do
  end do

  ! Replicate array to match array dimensions for multiplication below
  !do thetaID = 2,ntheta
  !  umag(:,:,:,:,thetaID) = umag(:,:,:,:,1)
  !end do

  ! Collapse the Gabor mode data into a 1D vector
  nmodes = this%nmodes 
  this%uhatR = this%scalefact*reshape(uR,(/nmodes/))
  this%uhatI = this%scalefact*reshape(uI,(/nmodes/))
  
  this%vhatR = this%scalefact*reshape(vR,(/nmodes/))
  this%vhatI = this%scalefact*reshape(vI,(/nmodes/))
  
  this%whatR = this%scalefact*reshape(wR,(/nmodes/))
  this%whatI = this%scalefact*reshape(wI,(/nmodes/))

  this%kx = reshape(k1,(/nmodes/))
  this%ky = reshape(k2,(/nmodes/))
  this%kz = reshape(k3,(/nmodes/))

  this%x = reshape(gmx,(/nmodes/))
  this%y = reshape(gmy,(/nmodes/))
  this%z = reshape(gmz,(/nmodes/))

  this%QHi = reshape(QHi,(/nmodes/))
  this%QHj = reshape(QHj,(/nmodes/))
  this%QHk = reshape(QHk,(/nmodes/))
  
  ! Rescale amplitudes using interpolated values of KE and L
  call this%rescaleAmplitudesUsingLocalParameters()

  ! Confirm modes are divergence free
  call assert(isOrthogonal(this%uhatR,this%vhatR,this%whatR,&
    this%kx,this%ky,this%kz),'Velocity not divergence free (R)')
  call assert(isOrthogonal(this%uhatI,this%vhatI,this%whatI,&
    this%kx,this%ky,this%kz),'Velocity not divergence free (I)')

  ! Debug checks
  if (this%debugChecks) then
    call this%doDebugChecks()
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
  end if
  
  ! Reset warning trigger
  doWarning = .true.

  ! Deallocate temporary arrays
  deallocate(kedge,kmag,dk,E,rand1,theta,kztemp,r,umag,thetaVel,uRmag,uImag,&
             rand2, p1x, p1y, p2x, p2y, p2z, orientationX, orientationY, &
             orientationZ, uR, uI, vR, vI, wR, wI, k1, k2, k3, gmx, gmy, gmz, &
             QHi, QHj, QHk)
  
  ! Acknowledge that isotropic modes are initialized
  call message(1,'Finished generating isotropic modes')

end subroutine

subroutine rescaleAmplitudesUsingLocalParameters(this)
  use GaborModeRoutines, only: interpToLocation
  class(enrichmentOperator), intent(inout), target :: this
  real(rkind), dimension(:,:,:), allocatable :: KE_h, L_h
  real(rkind), dimension(:,:,:), pointer :: KE_local, L_local
  real(rkind) :: KE_mode, L_mode, kmag, alpha
  integer :: n, QHi, QHj, QHk

  ! Get halo-padded arrays for KE and L
  call this%largeScales%pg%alloc_array(KE_h)
  call this%largeScales%pg%alloc_array(L_h)
 
  ! Compute local large scale kinetic energy
  KE_local => this%largeScales%rbuffxC(:,:,:,1)
  L_local  => this%largeScales%rbuffxC(:,:,:,2)

  KE_local = 1.d0!0.5d0*(this%largeScales%u**2 + this%largeScales%v**2 + this%largeScales%wC**2)
  L_local = 1.d0

  call this%largeScales%haloUpdateField(KE_local, KE_h)
  call this%largeScales%haloUpdateField(L_local, L_h)

  do n = 1,this%nmodes
    QHi = this%QHi(n)
    QHj = this%QHj(n)
    QHk = this%QHk(n)
    
    ! Interoplate KE and L to mode location
    call interpToLocation(KE_h,KE_mode,this%largeScales%dx, this%largeScales%dy, &
      this%largeScales%dz, this%largeScales%mesh(1,1,1,1), &
      this%largeScales%mesh(1,1,1,2), this%largeScales%mesh(1,1,1,3), &
      this%x(n), this%y(n), this%z(n))
    call interpToLocation(L_h,L_mode,this%largeScales%dx, this%largeScales%dy, &
      this%largeScales%dz, this%largeScales%mesh(1,1,1,1), &
      this%largeScales%mesh(1,1,1,2), this%largeScales%mesh(1,1,1,3), &
      this%x(n), this%y(n), this%z(n))

    kmag = sqrt(this%kx(n)**2 + this%ky(n)**2 + this%kz(n)**2)

    ! Compute the local scaling parameter based on QH-region quantities and local
    ! LES quantites
    alpha = 2*this%QHgrid%KE(QHi,QHj,QHk)*this%QHgrid%L(QHi,QHj,QHk)**5/&
      (2*KE_mode*L_mode**5)*&
      ((1.d0 + (kmag*this%QHgrid%L(QHi,QHj,QHk))**2)/&
      (1.d0 + (kmag*L_mode)**2))**(17.d0/6.d0)
    
    this%uhatR(n) = alpha*this%uhatR(n)
    this%uhatI(n) = alpha*this%uhatI(n)
    this%vhatR(n) = alpha*this%vhatR(n)
    this%vhatI(n) = alpha*this%vhatI(n)
    this%whatR(n) = alpha*this%whatR(n)
    this%whatI(n) = alpha*this%whatI(n)
  end do

  deallocate(KE_h, L_h)
  nullify(KE_local, L_local)
end subroutine

subroutine strainModes(this)
  use GaborModeRoutines, only: PFQ, getDtMax, rk4Step, small
  use decomp_2D,         only: nrank
  class(enrichmentOperator), intent(inout) :: this 
  real(rkind), dimension(this%nmodes) :: kabs, input
  real(rkind) :: output, tauEddy, S, L, KE, U, V, W
  real(rkind), dimension(3,3) :: dudx
  integer :: n, tid
  real(rkind) :: dt
  real(rkind), dimension(3) :: k, uRtmp, uItmp
  character(len=clen) :: mssg

  kabs = sqrt(this%kx*this%kx + this%ky*this%ky + this%kz*this%kz)
  input = -(kabs)**(-2.0_rkind)
 
  !$OMP PARALLEL SHARED(input) &
  !$OMP PRIVATE(n, output, dudx, L, KE, U, V, W, S, k, uRtmp, uItmp, tauEddy) &
  !$OMP PRIVATE(dt)
  !$OMP DO
  do n = 1,this%nmodes
    CALL PFQ(input(n),output)
    call this%getLargeScaleDataAtModeLocation(n,dudx,L,KE,U,V,W)
    S = sqrt(dudx(1,1)*dudx(1,1) + dudx(1,2)*dudx(1,2) + dudx(1,3)*dudx(1,3) + &
             dudx(2,1)*dudx(2,1) + dudx(2,2)*dudx(2,2) + dudx(2,3)*dudx(2,3) + &
             dudx(3,1)*dudx(3,1) + dudx(3,2)*dudx(3,2) + dudx(3,3)*dudx(3,3))
    call getDtMax(dt,[this%uhatR(n), this%vhatR(n), this%whatR(n)],&
                     [this%uhatI(n), this%vhatI(n), this%whatI(n)],&
                     S,kabs(n))
    
    k     = [this%kx(n),    this%ky(n),    this%kz(n)]
    uRtmp = [this%uhatR(n), this%vhatR(n), this%whatR(n)]
    uItmp = [this%uhatI(n), this%vhatI(n), this%whatI(n)]
    
    if (S < small) then
      tauEddy = 0.d0
    else
      tauEddy = this%ctauGlobal/S*((kabs(n))**(-2.0_rkind/3.0_rkind))/sqrt(output) 
    end if
    
    do tid = 1,nint(tauEddy/dt)
      call rk4Step(uRtmp,uItmp,k,dt,this%Anu,KE,L,this%numolec,dudx)
    end do
    this%kx(n) = k(1)
    this%ky(n) = k(2)
    this%kz(n) = k(3)
    this%uhatR(n) = uRtmp(1)
    this%uhatI(n) = uItmp(1)
    this%vhatR(n) = uRtmp(2)
    this%vhatI(n) = uItmp(2)
    this%whatR(n) = uRtmp(3)
    this%whatI(n) = uItmp(3)
    
    if (mod(n,10000) == 0 .and. nrank == 0) then
      write(mssg,'(F8.5,A)') real(n,rkind)/real(this%nmodes,rkind)*100.d0,'% complete'
      print*, trim(mssg)
    end if
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  call message(1,'Finished straining isotropic modes')
end subroutine

subroutine getLargeScaleDataAtModeLocation(this,gmID,dudx,L,KE,U,V,W)
  use GaborModeRoutines, only: interpToLocation, getNearestNeighborValue!, findMeshIdx

  class(enrichmentOperator), intent(inout) :: this
  integer, intent(in) :: gmID
  real(rkind), dimension(3,3), intent(out) :: dudx
  real(rkind), intent(out) :: L, KE, U, V, W
  integer :: i, j, idx, QHx, QHy, QHz

  call interpToLocation(this%uh, U,&
    this%largeScales%dx, this%largeScales%dy, this%largeScales%dz,&
    this%largeScales%mesh(1,1,1,1), this%largeScales%mesh(1,1,1,2), &
    this%largeScales%mesh(1,1,1,3), &
    this%x(gmID), this%y(gmID), this%z(gmID))
  call interpToLocation(this%vh, V,&
    this%largeScales%dx, this%largeScales%dy, this%largeScales%dz,&
    this%largeScales%mesh(1,1,1,1), this%largeScales%mesh(1,1,1,2), &
    this%largeScales%mesh(1,1,1,3), &
    this%x(gmID), this%y(gmID), this%z(gmID))
  call interpToLocation(this%wh, W,&
    this%largeScales%dx, this%largeScales%dy, this%largeScales%dz,&
    this%largeScales%mesh(1,1,1,1), this%largeScales%mesh(1,1,1,2), &
    this%largeScales%mesh(1,1,1,3), &
    this%x(gmID), this%y(gmID), this%z(gmID))
  
  idx = 1
  do i = 1,3
    do j = 1,3
      call interpToLocation(this%duidxj_h(:,:,:,idx), dudx(i,j),&
        this%largeScales%dx, this%largeScales%dy, this%largeScales%dz,&
        this%largeScales%mesh(1,1,1,1), this%largeScales%mesh(1,1,1,2), &
        this%largeScales%mesh(1,1,1,3), &
        this%x(gmID), this%y(gmID), this%z(gmID))
      idx = idx + 1
      ! dudx(1,1) = dudx = duidxj_h(:,:,:,1) see compute_duidxj() in igrid
      ! dudx(1,2) = dudy = duidxj_h(:,:,:,2) see compute_duidxj() in igrid
      ! dudx(1,3) = dudz = duidxj_h(:,:,:,3) see compute_duidxj() in igrid
      
      ! dudx(2,1) = dvdx = duidxj_h(:,:,:,4) see compute_duidxj() in igrid
      ! dudx(2,2) = dvdy = duidxj_h(:,:,:,5) see compute_duidxj() in igrid
      ! dudx(2,3) = dvdz = duidxj_h(:,:,:,6) see compute_duidxj() in igrid
      
      ! dudx(3,1) = dwdx = duidxj_h(:,:,:,7) see compute_duidxj() in igrid
      ! dudx(3,2) = dwdy = duidxj_h(:,:,:,8) see compute_duidxj() in igrid
      ! dudx(3,3) = dwdz = duidxj_h(:,:,:,9) see compute_duidxj() in igrid
    end do
  end do

  ! Find the index of the mode's QH region
  QHx = this%QHi(gmID)!findMeshIdx(this%x(gmID),this%QHgrid%xE)
  QHy = this%QHj(gmID)!findMeshIdx(this%y(gmID),this%QHgrid%yE)
  QHz = this%QHk(gmID)!findMeshIdx(this%z(gmID),this%QHgrid%zE)
  
  ! Use L and KE for the QH region that the mode resides in
  L = this%QHgrid%L(QHx,QHy,QHz)
  KE = this%QHgrid%KE(QHx,QHy,QHz)
end subroutine

subroutine updateSeeds(this,seed1,seed2,seed3,seed4,seed5,seed6,seed7,&
    QHi,QHj,QHk)
  class(enrichmentOperator), intent(inout) :: this
  integer, intent(inout) :: seed1, seed2, seed3, seed4, seed5, seed6, seed7
  integer, intent(in) :: QHi, QHj, QHk
  integer :: nxQH, nyQH, QHidx

  nxQH = this%QHgrid%gpC%xsz(1)
  nyQH = this%QHgrid%gpC%ysz(2)
  QHidx = (QHk-1)*nxQH*nyQH + (QHj-1)*nxQH + QHi
  
  seed1 = 1000*QHidx
  seed2 = 2000*QHidx
  seed3 = 3000*QHidx
  seed4 = 4000*QHidx
  seed5 = 5000*QHidx
  seed6 = 6000*QHidx
  seed7 = 7000*QHidx
end subroutine
