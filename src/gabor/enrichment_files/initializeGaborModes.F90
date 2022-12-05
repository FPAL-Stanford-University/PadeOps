subroutine generateIsotropicModes(this)
  use constants, only: rmaxInt
  use gridtools, only: logspace
  use arrayTools, only: sub2idx
  use random,    only: uniform_random
  use GaborModeRoutines, only: normalizeVec, isOrthogonal, getModelSpectrum, &
                               doWarning, interpToLocation, initializeSeeds, &
                               updateSeeds
  use mpi
  class(enrichmentOperator), intent(inout) :: this
  character(len=clen) :: mssg1, mssg2, mssg3, mssg4
  real(rkind) :: xmin, ymin, zmin
  integer :: i, j, k, kid, thetaID, nk, ntheta, nmodes
  real(rkind), dimension(7) :: seedFact
  integer :: isz, jsz, ksz
  real(rkind), dimension(:), allocatable :: kedge, kmag, dk, dkmodes, rand1, &
                                            theta, kztemp, r, rand2, thetaVel
  real(rkind) :: orientationX, orientationY, &
                 orientationZ, p1x, p1y, p1z, &
                 p2x, p2y, p2z, uRmag, uImag, &
                 umag, E
  
  ! Physical location of Gabor modes
  real(rkind), dimension(:,:,:,:), allocatable :: gmx, gmy, gmz

  ! Wave-vector components
  real(rkind), dimension(:,:,:,:), allocatable :: k1, k2, k3
  real(rkind) :: myKmag

  ! Velocity-vector amplitudes
  real(rkind), dimension(:,:,:,:), allocatable :: uR, uI, vR, vI, wR, wI

  ! Large scale quantities
  real(rkind) :: KE_loc, L_loc

  ! Misc
  integer :: ierr, kst, ken, n, thetaSt, thetaEn
  integer :: nxQH, nyQH, istQH, jstQH, kstQH, idxOld, idxCurrent

  call message("                                                                ")
  call message("================================================================")
  call message("============== Initializing isotropic Gabor modes ==============")
  write(mssg1,'(A)') 'Total modes to initialize:'
  write(mssg2,'(A,I2,A,I2,A)') '(', this%nk,' shells per QH region) X (',&
    this%ntheta, ' modes per shell) X ...'
  write(mssg3,'(A,I3,A,I3,A,I3,A)') &
    '(', this%QHgrid%nx,' QH regions in x) X (', &
    this%QHgrid%ny,' QH regions in y) X (', this%QHgrid%nz, ' QH regions in z) = ...'
  write(mssg4,'(I12,A)') this%nmodesGlobal, ' modes'
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
  allocate(kmag(nk),dk(nk),dkmodes(nk*ntheta))
  allocate(rand1(nk*ntheta))
  allocate(theta(ntheta),kztemp(ntheta),r(ntheta))
  allocate(thetaVel(nk*ntheta))
  allocate(rand2(nk*ntheta))

  ! Allocate memory for ND attribute arrays
  allocate(gmx(nk*ntheta,isz,jsz,ksz), gmy(nk*ntheta,isz,jsz,ksz), &
           gmz(nk*ntheta,isz,jsz,ksz))
  allocate(uR(nk*ntheta,isz,jsz,ksz),uI(nk*ntheta,isz,jsz,ksz),&
           vR(nk*ntheta,isz,jsz,ksz),vI(nk*ntheta,isz,jsz,ksz),&
           wR(nk*ntheta,isz,jsz,ksz),wI(nk*ntheta,isz,jsz,ksz))
  allocate(k1(nk*ntheta,isz,jsz,ksz),k2(nk*ntheta,isz,jsz,ksz),&
           k3(nk*ntheta,isz,jsz,ksz))

  nxQH = this%QHgrid%gpC%xsz(1)
  nyQH = this%QHgrid%gpC%ysz(2)

  istQH = this%QHgrid%gpC%xst(1)
  jstQH = this%QHgrid%gpC%xst(2)
  kstQH = this%QHgrid%gpC%xst(3)
  
  call initializeSeeds(seedFact, istQH, jstQH, kstQH, nxQH, nyQH)

  idxOld = sub2idx(istQH, jstQH, kstQH, nxQH, nyQH) - 1
  do k = 1,this%QHgrid%gpC%xsz(3)
    zmin = this%QHgrid%zE(k)
    do j = 1,this%QHgrid%gpC%xsz(2)
      ymin = this%QHgrid%yE(j)
      do i = 1,this%QHgrid%gpC%xsz(1)
        xmin = this%QHgrid%xE(i)
        
        idxCurrent = sub2idx(istQH + i - 1, jstQH + j - 1, kstQH + k - 1, nxQH, nyQH) 
        call updateSeeds(seedFact,idxOld,idxCurrent)
        idxOld = idxCurrent
!print*, (istQH + i - 1), (jstQH + j - 1), (kstQH + k - 1), nint(seedFact(1)*rmaxInt)
         
        ! Uniformily distribute modes in QH region
        call uniform_random(rand1,0.d0,1.d0,nint(seedFact(1)*rmaxInt))
        gmx(:,i,j,k) = xmin + this%QHgrid%dx*rand1
        
        call uniform_random(rand1,0.d0,1.d0,nint(seedFact(2)*rmaxInt))
        gmy(:,i,j,k) = ymin + this%QHgrid%dy*rand1

        call uniform_random(rand1,0.d0,1.d0,nint(seedFact(3)*rmaxInt))
        gmz(:,i,j,k) = zmin + this%QHgrid%dz*rand1
          
        call uniform_random(thetaVel,0.d0,2.d0*pi,nint(seedFact(6)*rmaxInt))
          
        call uniform_random(rand2,0.d0,1.d0/sqrt(3.d0),nint(seedFact(7)*rmaxInt))
        
        ! TODO: kmin and kmax should be functions of space
        ! call this%getKminKmax(i,j,k)
        kedge = logspace(log10(this%kmin),log10(this%kmax),nk+1)
        kmag = 0.5d0*(kedge(1:nk) + kedge(2:nk+1))
        dk = kedge(2:nk+1) - kedge(1:nk)
          
        do kid = 1,nk
          kst = (kid - 1)*ntheta + 1
          ken = kst + ntheta - 1

          ! Isotropically sample wave-vector components
          call uniform_random(theta,0.d0,2.d0*pi,nint(seedFact(4)*rmaxInt))
          call uniform_random(kztemp,-1.d0,1.d0,nint(seedFact(5)))

          ! Assign wave-vector components
          k3(kst:ken,i,j,k) = kmag(kid)*kztemp
          r = sqrt(1.d0 - kztemp*kztemp)
          k1(kst:ken,i,j,k) = kmag(kid)*r*cos(theta)
          k2(kst:ken,i,j,k) = kmag(kid)*r*sin(theta)
          dkmodes(kst:ken) = dk(kid)
        end do

        do n = 1,this%nk*this%ntheta

          ! get KE_loc and L_loc
          call interpToLocation(this%KEh,KE_loc,this%largeScales%dx,this%largeScales%dy,&
            this%largeScales%dz,this%largeScales%mesh(1,1,1,1),&
            this%largeScales%mesh(1,1,1,2),this%largeScales%mesh(1,1,1,3),&
            gmx(n,i,j,k),gmy(n,i,j,k),gmz(n,i,j,k))
          call interpToLocation(this%Lh,L_loc,this%largeScales%dx,this%largeScales%dy,&
            this%largeScales%dz,this%largeScales%mesh(1,1,1,1),&
            this%largeScales%mesh(1,1,1,2),this%largeScales%mesh(1,1,1,3),&
            gmx(n,i,j,k),gmy(n,i,j,k),gmz(n,i,j,k))
          
          myKmag = sqrt(k1(n,i,j,k)**2 + k2(n,i,j,k)**2 + k3(n,i,j,k)**2)
          
          E = getModelSpectrum(myKmag,KE_loc,L_loc)

          ! Amplitude of each mode such that the sum of mode energies gives
          ! the correct kinetic energy
          umag = sqrt(2*E*dkmodes(n)/ntheta)
  
          ! Get velocity vector orientations
          ! Generate basis of tangent plane to wavenumber shell
          p1x = 0.d0
          p1y = -k3(n,i,j,k)
          p1z =  k2(n,i,j,k)

          p2x =  k2(n,i,j,k)**2 + k3(n,i,j,k)**2
          p2y = -k1(n,i,j,k)*k2(n,i,j,k)
          p2z = -k1(n,i,j,k)*k3(n,i,j,k)

          call normalizeVec(p1x,p1y,p1z)
          call normalizeVec(p2x,p2y,p2z)

          ! Assign orientation of the velocity vectors
          orientationX = cos(thetaVel(n))*p1x + sin(thetaVel(n))*p2x
          orientationY = cos(thetaVel(n))*p1y + sin(thetaVel(n))*p2y
          orientationZ = cos(thetaVel(n))*p1z + sin(thetaVel(n))*p2z
          ! We have the modulus and orientation of the complex-valued amplitudes.
          ! Now assigning the real and imaginary components
          uRmag = rand2(n)*umag
          uImag = sqrt(umag*umag/3.d0 - uRmag*uRmag)

          uR(n,i,j,k) = uRmag*orientationX
          uI(n,i,j,k) = uImag*orientationX
              
          vR(n,i,j,k) = uRmag*orientationY
          vI(n,i,j,k) = uImag*orientationY
              
          wR(n,i,j,k) = uRmag*orientationZ
          wI(n,i,j,k) = uImag*orientationZ
        end do
        
      end do
    end do
  end do

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
  deallocate(kedge,kmag,dk,rand1,theta,kztemp,r,thetaVel,rand2,&
             uR, uI, vR, vI, wR, wI, k1, k2, k3, gmx, gmy, gmz)
  
  ! Acknowledge that isotropic modes are initialized
  call message(1,'Finished generating isotropic modes')

end subroutine

subroutine strainModes(this)
  use GaborModeRoutines, only: PFQ, getDtMax, rk4Step, small
  use decomp_2D,         only: nrank
  class(enrichmentOperator), intent(inout) :: this 
  real(rkind), dimension(this%nmodes) :: kabs, input
  real(rkind) :: output, tauEddy, S, L, KE
  real(rkind), dimension(3) :: Ui
  real(rkind), dimension(3,3) :: dudx
  integer :: n, tid
  real(rkind) :: dt
  real(rkind), dimension(3) :: k, uRtmp, uItmp, x
  character(len=clen) :: mssg

  kabs = sqrt(this%kx*this%kx + this%ky*this%ky + this%kz*this%kz)
  input = -(kabs)**(-2.0_rkind)
  x     = 0.d0
 
  !$OMP PARALLEL SHARED(input) &
  !$OMP PRIVATE(n, output, dudx, L, KE, Ui, S, k, uRtmp, uItmp, tauEddy) &
  !$OMP PRIVATE(dt)
  !$OMP DO
  do n = 1,this%nmodes
    CALL PFQ(input(n),output)
    call this%getLargeScaleDataAtModeLocation(n,dudx,L,KE,Ui)
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
  
    if (nint(tauEddy/dt) < 50) then
      dt = tauEddy/50.d0
    end if  

    do tid = 1,nint(tauEddy/dt)
      call rk4Step(uRtmp,uItmp,k,x,dt,this%Anu,KE,L,this%numolec,dudx,Ui)
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

subroutine getLargeScaleDataAtModeLocation(this,gmID,dudx,L,KE,Ui)
  use GaborModeRoutines, only: interpToLocation, getNearestNeighborValue!, findMeshIdx

  class(enrichmentOperator), intent(inout) :: this
  integer, intent(in) :: gmID
  real(rkind), dimension(3,3), intent(out) :: dudx
  real(rkind), intent(out) :: L, KE
  real(rkind), dimension(3), intent(out) :: Ui
  integer :: i, j, idx

  call interpToLocation(this%uh, Ui(1),&
    this%largeScales%dx, this%largeScales%dy, this%largeScales%dz,&
    this%largeScales%mesh(1,1,1,1), this%largeScales%mesh(1,1,1,2), &
    this%largeScales%mesh(1,1,1,3), &
    this%x(gmID), this%y(gmID), this%z(gmID))
  call interpToLocation(this%vh, Ui(2),&
    this%largeScales%dx, this%largeScales%dy, this%largeScales%dz,&
    this%largeScales%mesh(1,1,1,1), this%largeScales%mesh(1,1,1,2), &
    this%largeScales%mesh(1,1,1,3), &
    this%x(gmID), this%y(gmID), this%z(gmID))
  call interpToLocation(this%wh, Ui(3),&
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

  ! Use L and KE for the QH region that the mode resides in
  call interpToLocation(this%KEh,KE,this%largeScales%dx,this%largeScales%dy,&
    this%largeScales%dz,this%largeScales%mesh(1,1,1,1),&
    this%largeScales%mesh(1,1,1,2),this%largeScales%mesh(1,1,1,3),&
    this%x(gmID),this%y(gmID),this%z(gmID))
  call interpToLocation(this%Lh,L,this%largeScales%dx,this%largeScales%dy,&
    this%largeScales%dz,this%largeScales%mesh(1,1,1,1),&
    this%largeScales%mesh(1,1,1,2),this%largeScales%mesh(1,1,1,3),&
    this%x(gmID),this%y(gmID),this%z(gmID))
end subroutine
