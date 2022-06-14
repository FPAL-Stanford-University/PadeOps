subroutine generateIsotropicModes(this)
  use gridtools, only: logspace
  use random,    only: uniform_random
  use GaborModeRoutines, only: normalizeVec, isOrthogonal, getModelSpectrum, &
                               doWarning
  class(enrichmentOperator), intent(inout) :: this
  real(rkind) :: xmin, ymin, zmin
  integer :: i, j, k, kid, thetaID, nk, ntheta, nmodes, nx, ny, nz
  integer :: seed1=1, seed2=2, seed3=3, seed4=4, seed5=5, seed6=6, seed7=7
  real(rkind), dimension(:),         allocatable :: kedge, kmag, dk, rand1, &
                                                    theta, kztemp, r, E
  real(rkind), dimension(:,:,:,:,:), allocatable :: umag, thetaVel, &
                                                    uRmag, uImag, rand2, &
                                                    p1x, p1y, p1z, &
                                                    p2x, p2y, p2z, &
                                                    orientationX, orientationY, orientationZ
  
  ! Physical location of Gabor modes
  real(rkind), dimension(:,:,:,:), allocatable :: gmx, gmy, gmz

  ! Wave-vector components
  real(rkind), dimension(:,:,:,:,:), allocatable :: k1, k2, k3

  ! Velocity-vector amplitudes
  real(rkind), dimension(:,:,:,:,:), allocatable :: uR, uI, vR, vI, wR, wI


  call message("Initializing isotropic modes")
  
  ! Allocate memory
  nk = this%nk
  ntheta = this%ntheta 
  nx = this%QHgrid%nx 
  ny = this%QHgrid%ny 
  nz = this%QHgrid%nz

  allocate(kedge(nk+1))
  allocate(kmag(nk),dk(nk),E(nk))
  allocate(rand1(nk*ntheta))
  allocate(theta(ntheta),kztemp(ntheta),r(ntheta))
  allocate(umag(nx,ny,nz,nk,ntheta), thetaVel(nx,ny,nz,nk,ntheta), &
           uRmag(nx,ny,nz,nk,ntheta), uImag(nx,ny,nz,nk,ntheta),&
           rand2(nx,ny,nz,nk,ntheta), p1x(nx,ny,nz,nk,ntheta), &
           p1y(nx,ny,nz,nk,ntheta), p1z(nx,ny,nz,nk,ntheta), &
           p2x(nx,ny,nz,nk,ntheta), p2y(nx,ny,nz,nk,ntheta), &
           p2z(nx,ny,nz,nk,ntheta), orientationX(nx,ny,nz,nk,ntheta), &
           orientationY(nx,ny,nz,nk,ntheta), orientationZ(nx,ny,nz,nk,ntheta))

  ! Allocate memory for ND attribute arrays
  allocate(gmx(nx,ny,nz,nk*ntheta), gmy(nx,ny,nz,nk*ntheta), &
           gmz(nx,ny,nz,nk*ntheta))
  allocate(uR(nx,ny,nz,nk,ntheta),uI(nx,ny,nz,nk,ntheta),&
           vR(nx,ny,nz,nk,ntheta),vI(nx,ny,nz,nk,ntheta),&
           wR(nx,ny,nz,nk,ntheta),wI(nx,ny,nz,nk,ntheta))
  allocate(k1(nx,ny,nz,nk,ntheta),k2(nx,ny,nz,nk,ntheta),&
           k3(nx,ny,nz,nk,ntheta))

  ! Assign wave-vector magnitudes based on logarithmically spaced shells
  ! (non-dimensionalized with L)
  !kedge = logspace(kmin,kmax,nk+1)
  kedge = logspace(log10(this%kmin),log10(this%kmax),nk+1)
  kmag = 0.5d0*(kedge(1:nk) + kedge(2:nk+1))
  dk = kedge(2:nk+1) - kedge(1:nk)
  
  do k = 1,this%QHgrid%nx
    zmin = this%QHgrid%zE(k)
    do j = 1,this%QHgrid%ny
      ymin = this%QHgrid%yE(j)
      do i = 1,this%QHgrid%nz
        xmin = this%QHgrid%xE(i)

        !TODO: Need better seed selection 
        seed1 = 1 + i + j + k
        seed2 = 2 + i + j + k
        seed3 = 3 + i + j + k
        seed4 = 4 + i + j + k
        seed5 = 5 + i + j + k

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

        ! Amplitude of each mode such that the sum of mode amplitudes gives
        ! the correct kinetic energy
        umag(i,j,k,:,1) = sqrt(2*E*dk/ntheta)
  
      end do
    end do
  end do

  ! Replicate array to match array dimensions for multiplication below
  do thetaID = 2,ntheta
    umag(:,:,:,:,thetaID) = umag(:,:,:,:,1)
  end do

  ! Get velocity vector orientations
    ! Generate basis of tangent plane to wavenumber shell
    p1x = 0.d0
    p1y = -k3
    p1z = k2

    p2x = k2*k2 + k3*k3
    p2y = -k1*k2
    p2z = -k1*k3

    call normalizeVec(p1x,p1y,p1z)
    call normalizeVec(p2x,p2y,p2z)

    ! Assign orientation of the velocity vectors
    call uniform_random(thetaVel,0.d0,2.d0*pi,seed6)
    orientationX = cos(thetaVel)*p1x + sin(thetaVel)*p2x
    orientationY = cos(thetaVel)*p1y + sin(thetaVel)*p2y
    orientationZ = cos(thetaVel)*p1z + sin(thetaVel)*p2z

  ! We have the modulus and orientation of the complex-valued amplitudes.
  ! Now assing the real and imaginary components
    call uniform_random(rand2,0.d0,1.d0,seed7)
    uRmag = rand2*umag
    uImag = sqrt(umag*umag - uRmag*uRmag)

    uR = uRmag*orientationX
    uI = uImag*orientationX

    vR = uRmag*orientationY
    vI = uImag*orientationY

    wR = uRmag*orientationZ
    wI = uImag*orientationZ
  
  ! Collapse the Gabor mode data into a 1D vector
  nmodes = nk*ntheta
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
    this%kz,this%ky,this%kz),'Velocity not divergence free (R)')
  call assert(isOrthogonal(this%uhatI,this%vhatI,this%whatI,&
    this%kx,this%ky,this%kz),'Velocity not divergence free (I)')
  
  ! Reset warning trigger
  doWarning = .true.

  ! Deallocate temporary arrays
  deallocate(kedge,kmag,dk,E,rand1,theta,kztemp,r,umag,thetaVel,uRmag,uImag,&
             rand2, p1x, p1y, p2x, p2y, p2z, orientationX, orientationY, &
             orientationZ, uR, uI, vR, vI, wR, wI, k1, k2, k3, gmx, gmy, gmz)
  
  ! Acknowledge that isotropic modes are initialized
  call message('Finished generating isotropic modes')

end subroutine 

