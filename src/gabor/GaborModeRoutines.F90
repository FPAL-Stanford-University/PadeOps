module GaborModeRoutines
  ! This module is used to program the initialization of Gabor mode attributes
  ! (e.g. uhat, kx, etc.) without committing to an OOP paradigm. Will port this
  ! as need-be to the GaborMode and GaborModePop classes
  use kind_parameters, only: rkind, clen
  use domainSetup, only: getStartAndEndIndices, gpQHcent, kmin, kmax, dxQH, &
    dyQH, dzQH, xQHedge, yQHedge, zQHedge, finishDomainSetup, nprocX, nprocY, &
    nprocZ, nrankX, nrankY, nrankZ, nzF, nxF, nyF, nxsupp, nysupp, nzsupp, &
    periodic, gpFE, gpFC, dxF, dyF, dzF, xFh, yFh, zFh, decomp2Dpencil
  use constants, only: pi
  use decomp_2D
  use mpi
  use fortran_assert, only: assert
  use exits, only: gracefulExit
  use PadeDerOps, only: Pade6stagg
  implicit none

  ! Complex vector amplitudes
  real(rkind), dimension(:), allocatable :: uhatR, uhatI, vhatR, vhatI, whatR, whatI
  real(rkind), dimension(:,:,:,:,:), allocatable, private :: uR, uI, vR, vI, wR, wI

  ! Wave-vector components
  real(rkind), dimension(:), allocatable :: kx, ky, kz
  real(rkind), dimension(:,:,:,:,:), allocatable, private :: k1, k2, k3

  ! Physical location of Gabor modes
  real(rkind), dimension(:), allocatable :: gmxloc, gmyloc, gmzloc
  real(rkind), dimension(:,:,:,:), allocatable, private :: gmx, gmy, gmz

  ! Misc variables
  integer :: nmodes
  integer, dimension(:), allocatable :: nmodesALL

  ! Problem parameters
  integer :: nk, ntheta
  real(rkind) :: scalefact = 1.d0, ctau, Anu = 1.d-4

  ! Warning handling
  logical, private :: doWarning = .true.
  logical :: modeMemoryInitialized = .false., isotropicModesInitialized = .false.
  !logical, private :: modeMemoryInitialized = .false., isotropicModesInitialized = .false.
  logical, private :: finishModeInitialization = .false. ! <-- After straining procedure

  ! Array dimensions for the particular MPI partition. This is based on gpQHcent
  integer, private :: istQH, ienQH, jstQH, jenQH, kstQH, kenQH
  integer, private :: iszQH, jszQH, kszQH

  ! Array dimensions for the particular MPI partition base on gpFE
  integer, private :: istFE, ienFE, jstFE, jenFE, kstFE, kenFE
  integer, private :: iszFE, jszFE, kszFE

  ! Array indices for halo-padded arrays
  integer, private :: istFh, ienFh, jstFh, jenFh, kstFh, kenFh
  integer, private :: iszFh, jszFh, kszFh

  ! Grid partition for halo-padded arrays
  !type(decomp_info), allocatable :: gpFh
  
  ! start and end indices for all MPI ranks based on gpFE
  integer, dimension(:), allocatable :: istAll, ienAll, jstAll, &
    jenAll, kstAll, kenAll

  ! Send and receive requests for MPI_Isend and MPI_Irecv
  integer, dimension(:), allocatable :: sendReqHi, sendReqLo, &
    recvReqHi, recvReqLo

  ! Halo padded velocity fields
  real(rkind), dimension(:,:,:), allocatable :: uFh, vFh, wFh
  real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, &
    buff3, buff4

  ! Rendered velocity fields
  real(rkind), dimension(:,:,:), allocatable :: uG, vG, wG
  real(rkind), dimension(:,:,:), allocatable :: uGout, vGout, wGout

  ! Interpolator for output velocity field
  type(Pade6stagg) :: interp

  ! Number of OpenMP threads
  integer :: nthreads = 1

  interface mpiIsendIrecv
    module procedure mpiIsendIrecv2D, mpiIsendIrecv3D
  end interface

  contains
    include "GaborMode_files/renderVelocity.F90"
    include "GaborMode_files/GaborMode_MPI.F90"

    subroutine initializeModes(inputfile)
      use domainSetup
      character(len=*), intent(in) :: inputfile
      integer :: ierr, ioUnit, scheme
      namelist /GABOR/ nk, ntheta, scalefact, ctau, Anu

      ! Verify that the domain has been setup
      call assert(finishDomainSetup,'Attempting to initialize Gabor modes '//&
        'without setting up domain first.')
      
      ! Read inputfile
      ioUnit = 1
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=ioUnit, NML=GABOR)
      close(ioUnit)

      call allocateMemory()

      if (.not. periodic(3)) then
        scheme = 1 ! 0: fd02, 1: cd06, 2: fourierColl
        call interp%init(gpFC,gpFC,gpFE,gpFE,dzF,scheme,periodic(3))
      end if

      modeMemoryInitialized = .true.
    end subroutine

    subroutine allocateMemory()
      integer :: ierr
      integer, dimension(:), allocatable :: myIst, myIen, myJst, myJen, myKst, myKen
      integer :: ist, ien, jst, jen, kst, ken, isz, jsz, ksz

      call getStartAndEndIndices(gpQHcent,istQH,ienQH,jstQH,jenQH,kstQH,kenQH,iszQH,jszQH,kszQH)
      
      ! ND arrays
      allocate(uR(iszQH,jszQH,kszQH,nk,ntheta), uI(iszQH,jszQH,kszQH,nk,ntheta))
      allocate(vR(iszQH,jszQH,kszQH,nk,ntheta), vI(iszQH,jszQH,kszQH,nk,ntheta))
      allocate(wR(iszQH,jszQH,kszQH,nk,ntheta), wI(iszQH,jszQH,kszQH,nk,ntheta))

      allocate(k1(iszQH,jszQH,kszQH,nk,ntheta))
      allocate(k2(iszQH,jszQH,kszQH,nk,ntheta))
      allocate(k3(iszQH,jszQH,kszQH,nk,ntheta))

      allocate(gmx(iszQH,jszQH,kszQH,nk*ntheta))
      allocate(gmy(iszQH,jszQH,kszQH,nk*ntheta))
      allocate(gmz(iszQH,jszQH,kszQH,nk*ntheta))
      
      ! Allocate memory for 1D vectors
      nmodes = iszQH*jszQH*kszQH*nk*ntheta
      allocate(uhatR(nmodes), uhatI(nmodes))
      allocate(vhatR(nmodes), vhatI(nmodes))
      allocate(whatR(nmodes), whatI(nmodes))

      allocate(kx(nmodes), ky(nmodes), kz(nmodes))

      allocate(gmxloc(nmodes), gmyloc(nmodes), gmzloc(nmodes))
      
      ! Send the number of modes on each partition
      allocate(nmodesALL(nproc))
      call MPI_Allgather(nmodes,1,MPI_INTEGER,nmodesALL,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      ! Allocate memory for the velocity field to be written
      if (.not. periodic(3)) then
        call getStartAndEndIndices(gpFC,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
        allocate(uGout(isz,jsz,ksz),&
                 vGout(isz,jsz,ksz),&
                 wGout(isz,jsz,ksz))
      end if
      call getStartAndEndIndices(gpFE,istFE,ienFE,jstFE,jenFE,kstFE,kenFE,iszFE,jszFE,kszFE)
      allocate(uG(iszFE,jszFE,kszFE),&
               vG(iszFE,jszFE,kszFE),&
               wG(iszFE,jszFE,kszFE))
      uG = 0.d0; vG = 0.d0; wG = 0.d0

! I need things below this line if implementing my own update_halo
!-----------------------------------------------------------------
      istFh = istFE - nxsupp/2; ienFh = ienFE + nxsupp/2
      jstFh = jstFE - nysupp/2; jenFh = jenFE + nysupp/2
      kstFh = kstFE - nzsupp/2; kenFh = kenFE + nzsupp/2
      if (.not. periodic(1)) then
        istFh = max(1,istFh)
        ienFh = min(nxF+1,ienFh)
      end if
      if (.not. periodic(2)) then
        jstFh = max(1,jstFh)
        jenFh = min(nyF+1,jenFh)
      end if
      if (.not. periodic(3)) then
        kstFh = max(1,kstFh)
        kenFh = min(nzF+1,kenFh)
      end if
      iszFh = ienFh - istFh + 1
      jszFh = jenFh - jstFh + 1
      kszFh = kenFh - kstFh + 1
      allocate(uFh(istFh:ienFh,jstFh:jenFh,kstFh:kenFh))
      allocate(vFh(istFh:ienFh,jstFh:jenFh,kstFh:kenFh))
      allocate(wFh(istFh:ienFh,jstFh:jenFh,kstFh:kenFh))

      select case (decomp2Dpencil)
      case ('x')
        allocate(buff1(iszFh,jszFh,nzsupp/2))
        allocate(buff2(iszFh,jszFh,nzsupp/2))
        allocate(buff3(iszFh,nysupp/2,kszFh))
        allocate(buff4(iszFh,nysupp/2,kszFh))

        allocate(sendReqHi(kszFh),sendReqLo(kszFh))
        allocate(recvReqHi(kszFh),recvReqLo(kszFh))
      case ('y')
        call assert(.false.,'TODO: y pencil stuff')
      case ('z')
        call assert(.false.,'TODO: z pencil stuff')
      case default 
        call gracefulExit("decomp2Dpencil must be 'x', 'y', or 'z'",ierr)
      end select

      ! We need start and end indices of each MPI rank for halo exchange 
      ! during velocity rendering
      allocate(myIst(nproc),myIen(nproc),myJst(nproc),myJen(nproc),&
        myKst(nproc),myKen(nproc))
      allocate(istAll(nproc),ienAll(nproc),jstAll(nproc),jenAll(nproc),&
        kstAll(nproc),kenAll(nproc))

      myIst = istFE; myIen = ienFE
      myJst = jstFE; myJen = jenFE
      myKst = kstFE; myKen = kenFE

      call MPI_Alltoall(myIst,1,MPI_INTEGER,istAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Alltoall(myIen,1,MPI_INTEGER,ienAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Alltoall(myJst,1,MPI_INTEGER,jstAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Alltoall(myJen,1,MPI_INTEGER,jenAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Alltoall(myKst,1,MPI_INTEGER,kstAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Alltoall(myKen,1,MPI_INTEGER,kenAll,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      deallocate(myIst,myIen,myJst,myJen,myKst,myKen)
    end subroutine

    subroutine generateIsotropicModes()
      use exits, only: message
      use gridtools, only: logspace
      use largeScalesMod, only: KE, L, computeLrgSclQOIs
      use random, only: uniform_random
      integer :: i, j, k, kid, thetaID
      real(rkind), dimension(nk+1) :: kedge
      real(rkind), dimension(nk) :: kmag, dk
      real(rkind), dimension(nk*ntheta) :: rand1
      real(rkind), dimension(ntheta) :: theta, kztemp, r
      real(rkind), dimension(nk) :: E
      real(rkind), dimension(iszQH,jszQH,kszQH,nk,ntheta) :: umag, thetaVel
      real(rkind), dimension(iszQH,jszQH,kszQH,nk,ntheta) :: uRmag, uImag, rand2
      real(rkind), dimension(iszQH,jszQH,kszQH,nk,ntheta) :: p1x, p1y, p1z
      real(rkind), dimension(iszQH,jszQH,kszQH,nk,ntheta) :: p2x, p2y, p2z
      real(rkind), dimension(iszQH,jszQH,kszQH,nk,ntheta) :: orientationX, &
                                                             orientationY, &
                                                             orientationZ
      integer :: seed1=1, seed2=2, seed3=3, seed4=4, seed5=5, seed6=6, seed7=7
      real(rkind) :: xmin, ymin, zmin


      ! First make sure we have the large scale info we need
      call assert(computeLrgSclQOIs,"Can't initialize Gabor modes w/o "//&
        "large scale quantities")

      ! Confirm we have allocated the required memory
      call assert(modeMemoryInitialized,"Gabor mode memory not allocated")

      call message("Initializing isotropic modes")
      
      ! Assign wave-vector magnitudes based on logarithmically spaced shells
      ! (non-dimensionalized with L)
      !kedge = logspace(kmin,kmax,nk+1)
      kedge = logspace(log10(kmin),log10(kmax),nk+1)
      kmag = 0.5d0*(kedge(1:nk) + kedge(2:nk+1))
      dk = kedge(2:nk+1) - kedge(1:nk)
      
      do k = 1,kszQH
        zmin = zQHedge(k)
        do j = 1,jszQH
          ymin = yQHedge(j)
          do i = 1,iszQH
            xmin = xQHedge(i)
    
            seed1 = 1 + (istQH + i - 1) + (jstQH + j - 1) + (kstQH + k - 1)
            seed2 = 2 + (istQH + i - 1) + (jstQH + j - 1) + (kstQH + k - 1)
            seed3 = 3 + (istQH + i - 1) + (jstQH + j - 1) + (kstQH + k - 1)
            seed4 = 4 + (istQH + i - 1) + (jstQH + j - 1) + (kstQH + k - 1)
            seed5 = 5 + (istQH + i - 1) + (jstQH + j - 1) + (kstQH + k - 1)

            ! Uniformily distribute modes in QH region
            call uniform_random(rand1,0.d0,1.d0,seed1)
            gmx(i,j,k,:) = xmin + dxQH*rand1
            
            call uniform_random(rand1,0.d0,1.d0,seed2)
            gmy(i,j,k,:) = ymin + dyQH*rand1

            call uniform_random(rand1,0.d0,1.d0,seed3)
            gmz(i,j,k,:) = zmin + dzQH*rand1
            
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
            E = getModelSpectrum(kmag,KE(i,j,k),L(i,j,k),nk)

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
      uhatR = scalefact*reshape(uR,(/nmodes/))
      uhatI = scalefact*reshape(uI,(/nmodes/))
      
      vhatR = scalefact*reshape(vR,(/nmodes/))
      vhatI = scalefact*reshape(vI,(/nmodes/))
      
      whatR = scalefact*reshape(wR,(/nmodes/))
      whatI = scalefact*reshape(wI,(/nmodes/))

      kx = reshape(k1,(/nmodes/))
      ky = reshape(k2,(/nmodes/))
      kz = reshape(k3,(/nmodes/))

      gmxloc = reshape(gmx,(/nmodes/))
      gmyloc = reshape(gmy,(/nmodes/))
      gmzloc = reshape(gmz,(/nmodes/))

      ! Confirm modes are divergence free
      call assert(isOrthogonal(uhatR,vhatR,whatR,kx,ky,kz),'Velocity not '//&
        'divergence free (R)')
      call assert(isOrthogonal(uhatI,vhatI,whatI,kx,ky,kz),'Velocity not '//&
        'divergence free (I)')
      
      ! Reset warning trigger
      doWarning = .true.

      ! Acknowledge that isotropic modes are initialized
      isotropicModesInitialized = .true.
      call message('Finished generating isotropic modes')
    end subroutine

    subroutine finalizeGaborModes()
      if (allocated(uhatR)) deallocate(uhatR)
      if (allocated(uhatI)) deallocate(uhatI)
      if (allocated(vhatR)) deallocate(vhatR)
      if (allocated(vhatI)) deallocate(vhatI)
      if (allocated(whatR)) deallocate(whatR)
      if (allocated(whatI)) deallocate(whatI)

      if (allocated(uR)) deallocate(uR)
      if (allocated(uI)) deallocate(uI)
      if (allocated(vR)) deallocate(vR)
      if (allocated(vI)) deallocate(vI)
      if (allocated(wR)) deallocate(wR)
      if (allocated(wI)) deallocate(wI)

      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kz)) deallocate(kz)

      if (allocated(k1)) deallocate(k1)
      if (allocated(k2)) deallocate(k2)
      if (allocated(k3)) deallocate(k3)

      if (allocated(gmxloc)) deallocate(gmxloc)
      if (allocated(gmyloc)) deallocate(gmyloc)
      if (allocated(gmzloc)) deallocate(gmzloc)

      if (allocated(gmx)) deallocate(gmx)
      if (allocated(gmy)) deallocate(gmy)
      if (allocated(gmz)) deallocate(gmz)

      if (allocated(nmodesALL)) deallocate(nmodesALL)
        
      !if (allocated(gpFh)) then
      !  call decomp_info_finalize(gpFh)
      !  deallocate(gpFh)
      !end if

      if (allocated(uFh)) deallocate(uFh)
      if (allocated(vFh)) deallocate(vFh)
      if (allocated(wFh)) deallocate(wFh)

      if (allocated(uG)) deallocate(uG)
      if (allocated(vG)) deallocate(vG)
      if (allocated(wG)) deallocate(wG)

      if (allocated(uGout)) deallocate(uGout)
      if (allocated(vGout)) deallocate(vGout)
      if (allocated(wGout)) deallocate(wGout)

      if (allocated(istAll)) deallocate(istAll)
      if (allocated(ienAll)) deallocate(ienAll)
      if (allocated(jstAll)) deallocate(jstAll)
      if (allocated(jenAll)) deallocate(jenAll)
      if (allocated(kstAll)) deallocate(kstAll)
      if (allocated(kenAll)) deallocate(kenAll)

      if (allocated(sendReqLo)) deallocate(sendReqLo)
      if (allocated(sendReqHi)) deallocate(sendReqHi)
      if (allocated(recvReqLo)) deallocate(recvReqLo)
      if (allocated(recvReqHi)) deallocate(recvReqHi)

      if (.not. periodic(3)) call interp%destroy()
    end subroutine

    function getModelSpectrum(k,KE,L,nk) result(E) 
      use exits, only: warning
      real(rkind), dimension(:), intent(in) :: k
      integer, intent(in) :: nk
      real(rkind), dimension(nk) :: kL, fL, feta
      real(rkind), intent(in) :: KE, L
      real(rkind) :: C
      real(rkind), dimension(nk) :: E

      ! Non-dimensional wavenumber
      kL = k*L

      ! Model coefficient chosen such that the integrated spectrum equals total
      ! kinetic energy
      C = 0.4843d0

      ! Energetic scalse
      fL = kL**4.d0/(1.d0 + kL*kL)**(17.d0/6.d0)

      ! TODO: Dissipative scales
      feta = 1.d0
      if (doWarning) then
        call warning("WARNING: Finite Re model spectrum is not implemented")
        doWarning = .false.
      end if

      ! Model spectrum
      E = C*2.d0*KE*L*fL*feta
    end function

    pure subroutine normalizeVec(x,y,z)
      real(rkind), dimension(:,:,:,:,:), intent(inout) :: x, y, z
      real(rkind), dimension(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5)) :: mag
    
      mag = sqrt(x*x + y*y + z*z)
      x = x/mag
      y = y/mag
      z = z/mag
    end subroutine

    pure function isOrthogonal(a1,a2,a3,b1,b2,b3) result(TF)
      real(rkind), dimension(:), intent(in) :: a1, a2, a3, b1, b2, b3
      logical :: TF
      real(rkind) :: maxDiv, small

      small = 1.d-14
    
      maxDiv = maxval(a1*b1 + a2*b2 + a3*b3)
      TF = .false.
      if (maxDiv < small) TF = .true.
    end function
end module
