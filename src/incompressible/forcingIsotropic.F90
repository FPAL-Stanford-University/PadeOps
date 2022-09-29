module forcingmod
   use kind_parameters, only: rkind, mpirkind, clen
   use decomp_2d
   use decomp_2d_io
   use constants,       only: im0, one, zero, two, pi, ieps
   use spectralMod,     only: spectral 
   use exits,           only: GracefulExit, message, warning
   use mpi
   use fortran_assert,  only: assert 
   use random,          only: uniform_random, randperm
   use arrayTools,      only: findGL 
   use PadePoissonMod,  only: padePoisson 
   use reductions,      only: p_sum, p_maxval

   implicit none
   private
   public :: HIT_shell_forcing

   integer :: ierr

   type :: HIT_shell_forcing
      private
      type(decomp_info), pointer :: sp_gpC, sp_gpE, gpC
      real(rkind) :: kmin, kmax
      integer :: Nwaves, Nmodes
      real(rkind), public :: EpsAmplitude
      class(spectral), pointer :: spectC, spectE

      integer, dimension(:), allocatable :: wave_x, wave_y, wave_z
      complex(rkind), dimension(:,:,:), allocatable, public :: uhat, vhat, what, &
        fxhat_old, fyhat_old, fzhat_old, fxhatwrite, fyhatwrite, fzhatwrite
      complex(rkind), dimension(:,:,:), pointer :: fxhat, fyhat, fzhat, cbuffzE, cbuffyE, cbuffyC

      real(rkind), dimension(:), allocatable :: kabs_sample, zeta_sample, theta_sample
      integer, public :: seed0, seed1, seed2, seed3
      real(rkind) :: Nwaves_rkind
      real(rkind), dimension(:), allocatable :: tmpModes
      real(rkind) :: alpha_t = 1.d0 
      real(rkind) :: normfact = 1.d0, A_force = 1.d0 
      integer     :: DomAspectRatioZ
      logical :: useLinearForcing, firstCall
      
      ! Ryan's additions/changes:
      real(rkind), dimension(:,:,:), allocatable :: k1, k2, k3, kmag
      integer, dimension(:), allocatable :: k1_1d, k2_1d, k3_1d
      integer, dimension(:), allocatable :: i, j, k
      integer :: version
      real(rkind), dimension(:,:,:), pointer :: rbuffxC
      complex(rkind), dimension(:,:,:), allocatable :: cbuffzC
      logical :: storeForce = .false.
      logical :: confirmEnergyInjectionRate

      ! Domain info
      real(rkind), pointer :: dx, dy, dz, Lx, Ly, Lz
    
      contains
      procedure          :: init
      procedure          :: destroy
      procedure, private :: update_seeds
      procedure, private :: pick_random_wavenumbers
      procedure, private :: pick_random_wavenumbersV2
      procedure, private :: compute_forcing
      procedure, private :: embed_forcing_mode
      procedure          :: getRHS_HITforcing
      procedure, private :: scrubConjPairs
      !procedure, private :: scrubConjPairsV2
      procedure, private :: replaceConjPartners
      procedure          :: dumpForcing
      procedure, private :: prepAndDumpField
      procedure          :: computeEnergyInjectionRate
      procedure          :: takeFFTz
      procedure, private :: takeIFFT3
   end type 

contains

subroutine init(this, inputfile, gpC, sp_gpC, sp_gpE, spectC, spectE, cbuffyE, cbuffyC, &
    cbuffzE, cbuffzC, rbuffxC, tidStart)
   class(HIT_shell_forcing), intent(inout) :: this
   character(len=*), intent(in) :: inputfile
   type(decomp_info), intent(in), target :: sp_gpC, sp_gpE, gpC
   integer, intent(in) :: tidStart
   complex(rkind), dimension(:,:,:  ), intent(inout), target :: cbuffzE, cbuffyE, cbuffyC
   complex(rkind), dimension(:,:,:,:), intent(inout), target :: cbuffzC
   real(rkind), dimension(:,:,:), intent(inout), target :: rbuffxC
   class(spectral), intent(in), target :: spectC, spectE
   integer :: RandSeedToAdd = 0, DomAspectRatioZ = 1
   real(rkind) :: alpha_t = 1.d0 
   integer :: Nwaves = 20
   real(rkind) :: kmin = 2.d0, kmax = 10.d0, EpsAmplitude = 0.1d0, A_force = 1.d0 
   logical :: useLinearForcing = .false. 
   real(rkind) :: filtfact_linForcing = 0.5d0
   integer :: nforce, version = 1
   logical :: storeForce = .false.
   logical :: confirmEnergyInjectionRate = .false.

   namelist /HIT_Forcing/ kmin, kmax, Nwaves, EpsAmplitude, RandSeedToAdd, &
     DomAspectRatioZ, alpha_t, useLinearForcing, filtfact_linForcing, &
     version, storeForce, confirmEnergyInjectionRate, A_force

   open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=123, NML=HIT_Forcing)
   close(123)

   if(DomAspectRatioZ < 1) then
       call GracefulExit("Aspect ratio in z must be greater than or equal to 1", 111)
   endif

   this%A_force = A_force
   this%kmin = kmin
   this%kmax = kmax
   this%EpsAmplitude = EpsAmplitude
   this%Nwaves = Nwaves
   this%DomAspectRatioZ = DomAspectRatioZ
   this%sp_gpC => sp_gpC
   this%sp_gpE => sp_gpE
   this%gpC    => gpC
   this%spectC => spectC
   this%spectE => spectE
   this%useLinearForcing = useLinearForcing
   this%version = version
   this%storeForce = storeForce
   this%confirmEnergyInjectionRate = confirmEnergyInjectionRate

   allocate(this%uhat (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%vhat (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%what (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%fxhat_old (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%fyhat_old (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%fzhat_old (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   !if (storeForce) then
     allocate(this%fxhatwrite (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
     allocate(this%fyhatwrite (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
     allocate(this%fzhatwrite (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   !end if
   allocate(this%cbuffzC(sp_gpC%zsz(1),sp_gpC%zsz(2),sp_gpC%zsz(3)))

   this%seed0 = tidStart + RandSeedToAdd
   call this%update_seeds()
  
   ! Define the wave-vector grid
   nforce = ceiling(this%kmax)*2
   if (mod(ceiling(this%kmax),2) .ne. 0) nforce = ceiling(this%kmax + 1.d0)*2

   select case (this%version)
     case (1)
       allocate(this%kabs_sample(Nwaves))
       allocate(this%theta_sample(Nwaves))
       allocate(this%zeta_sample(Nwaves))
       allocate(this%tmpModes(Nwaves))

       allocate(this%wave_x(Nwaves))
       allocate(this%wave_y(Nwaves))
       allocate(this%wave_z(Nwaves))

       call message("Initializating VERSION 1 of HIT_shell_forcing")
     case (2)
       call message("Initializating VERSION 2 of HIT_shell_forcing")
       allocate(this%k1  (nforce/2,nforce,nforce))
       allocate(this%k2  (nforce/2,nforce,nforce))
       allocate(this%k3  (nforce/2,nforce,nforce))
       allocate(this%kmag(nforce/2,nforce,nforce))
       call getIsotropicWaveVectorComponents(this%k1,this%k2,this%k3,nforce) 
       this%kmag = sqrt(this%k1*this%k1 + this%k2*this%k2 + this%k3*this%k3)

       call assert(nforce/2 < gpC%xsz(1)/2,&
         'Need to allocated nforce/2+1 -- forcingIsotropic.F90')
       allocate(this%k1_1d(nforce/2), this%k2_1d(nforce), this%k3_1d(nforce))
       this%k1_1d = nint(this%k1(:,1,1))
       this%k2_1d = nint(this%k2(1,:,1))
       this%k3_1d = nint(this%k3(1,1,:))

       ! Find the indices of each mode that satisfies kmin <= k <= kmax
       call findGL(this%kmag,this%kmin,this%kmax,this%i,this%j,this%k)
       call assert(size(this%i) >= this%Nwaves,"There are fewer than Nwaves modes "//&
         "between kmin and kmax -- forcingIsotropic.F90")

       ! Shuffle the modes before scrubbing conjugate pairs so we don't get biased statistics
       ! (do it a few times to make sure it is well shuffled)
       !call randomShuffle(this%i, this%j, this%k, this%seed0)
       !call randomShuffle(this%i, this%j, this%k, this%seed1)
       !call randomShuffle(this%i, this%j, this%k, this%seed2)
       !call randomShuffle(this%i, this%j, this%k, this%seed3)
  
       ! Remove complex conjugate pairs
       !call this%scrubConjPairsV2(this%i,this%j,this%k,&
       !  nint(this%k1(:,1,1)),nint(this%k2(1,:,1)),nint(this%k3(1,1,:)))
      
       ! Allocate Nwaves + a buffer of Nwaves 
       allocate(this%wave_x(Nwaves+50))
       allocate(this%wave_y(Nwaves+50))
       allocate(this%wave_z(Nwaves+50))
     case default
       call assert(.false.,'Invalid version number -- forcingIsotropic.F90')
   end select

   this%firstCall = .true.
   this%fxhat    => cbuffzC(:,:,:,1)
   this%fyhat    => cbuffzC(:,:,:,2)
   this%fzhat    => cbuffzC(:,:,:,3)
   this%cbuffzE  => cbuffzE 
   this%cbuffyE  => cbuffyE
   this%cbuffyC  => cbuffyC
   this%rbuffxC  => rbuffxC

   this%normfact = (real(spectC%physdecomp%xsz(1),rkind)*&
                  & real(spectC%physdecomp%ysz(2),rkind)*&
                  & real(spectC%physdecomp%zsz(3),rkind))**2.d0
  
   this%Nwaves_rkind = real(this%Nwaves, rkind)

   if (this%useLinearForcing) then
      call this%spectC%init_HIT_linearForcing(filtfact_linForcing)
   end if 
end subroutine 

subroutine update_seeds(this)
   class(HIT_shell_forcing), intent(inout) :: this
   integer :: big

   big = huge(0)
   if (big - this%seed0 < 10000000 .or. &
       big - this%seed1 < 10000000 .or. &
       big - this%seed2 < 10000000 .or. &
       big - this%seed3 < 10000000 ) then
     this%seed0 = 1
   end if
   this%seed0 = abs(this%seed0 + 2223345) 
   this%seed1 = abs(this%seed0 + 1423246)
   this%seed2 = abs(this%seed0 + 8723446)
   this%seed3 = abs(this%seed0 + 3423444)
end subroutine


subroutine pick_random_wavenumbers(this)
   class(HIT_shell_forcing), intent(inout) :: this

   call uniform_random(this%kabs_sample, this%kmin, this%kmax, this%seed1)
   call uniform_random(this%zeta_sample, -one, one, this%seed2)
   call uniform_random(this%theta_sample, zero, two*pi, this%seed3)
 
   call convertCylindricalToSphericalWaveVectors(&
     this%wave_x, this%wave_y, this%wave_z, this%tmpModes, &
     this%kabs_sample, this%zeta_sample, this%theta_sample)

   call this%scrubConjPairs(this%wave_x, this%wave_y, this%wave_z)

end subroutine

subroutine convertCylindricalToSphericalWaveVectors(kx, ky, kz, tmp, kabs, &
    zeta, theta)
  integer, dimension(:), intent(inout) :: kx, ky, kz
  real(rkind), dimension(:), intent(inout) :: tmp
  real(rkind), dimension(:), intent(in) :: kabs, zeta, theta
   
  tmp = kabs*sqrt(1 - zeta**2)*cos(theta)
  where(tmp < 0) tmp = -tmp
  kx = nint(tmp)
  
  tmp = kabs*sqrt(1 - zeta**2)*sin(theta)
  ky = nint(tmp)
  
  tmp = kabs*zeta
  kz = nint(tmp)

end subroutine

subroutine pick_random_wavenumbersV2(this,Nmodes)
  class(HIT_shell_forcing), intent(inout) :: this
  integer, intent(out) :: Nmodes
  integer :: n, Nadd1, Nadd2
  !integer, dimension(this%Nwaves) :: i, j, k
  
  ! Randomly sort the arrays
  call randomShuffle(this%i,this%j,this%k,this%seed3)

  ! Choose the first N modes
  Nmodes = this%Nwaves
  Nadd1 = 1
  Nadd2 = 0
  do while ((Nadd1 - Nadd2) /= 0)
    Nadd1 = NmodesToAdd(this%i, this%j, this%k, this%k1_1d, this%k2_1d, &
     this%k3_1d, Nmodes + Nadd2)
    Nadd2 = NmodesToAdd(this%i, this%j, this%k, this%k1_1d, this%k2_1d, &
     this%k3_1d, Nmodes + Nadd1)
  end do
  Nmodes = Nmodes + Nadd2
  !i = this%i(1:this%Nwaves)
  !j = this%j(1:this%Nwaves)
  !k = this%k(1:this%Nwaves)
  
  this%wave_x = 0
  this%wave_y = 0
  this%wave_z = 0
  call assert(Nmodes <= this%Nwaves+50,&
    'Did not allocated enough memory for wave_x, wave_y, and wave_z '//&
    '-- forcingIsotropic.F90')

  do n = 1,Nmodes
    this%wave_x(n) = this%k1(this%i(n),this%j(n),this%k(n)) 
    this%wave_y(n) = this%k2(this%i(n),this%j(n),this%k(n)) 
    this%wave_z(n) = this%k3(this%i(n),this%j(n),this%k(n)) 
  end do

end subroutine

function NmodesToAdd(i,j,k,k1,k2,k3,Ninit) result(Nadd)
  integer, dimension(:), intent(in) :: i, j, k, k1, k2, k3
  integer, intent(in) :: Ninit
  integer :: Nadd, n, m

  Nadd = 0
  do n = 1,Ninit-1
    if (k1(i(n)) == 0) then
      do m = n+1,Ninit
        if(k1(i(m)) == 0         .and. &
           k2(j(n)) == -k2(j(m)) .and. &
           k3(k(n)) == -k3(k(m))) then
          Nadd = Nadd + 1
        end if
      end do
    else
      continue
    end if
  end do
end function

subroutine randomShuffle(i,j,k,seed)
  integer, dimension(:), intent(inout) :: i, j, k
  integer, intent(in) :: seed
  integer, dimension(size(i)) :: idx

  call randperm(size(i),seed,idx)

  i = i(idx)
  j = j(idx)
  k = k(idx)
  
end subroutine

!subroutine scrubConjPairsV2(this,i,j,k,k1,k2,k3)
!  use sorting_mod, only: binary_sort
!  ! Ensure none of the wavevectors indexed by i, j, and k correspond to
!  ! complex conjugate (C.C.) modes
!
!  class(HIT_shell_forcing), intent(inout) :: this
!  integer, dimension(:), allocatable, intent(inout) :: i, j, k
!  integer, dimension(:), intent(in) :: k1, k2, k3
!  integer, dimension(:), allocatable :: icpy, jcpy, kcpy
!  integer :: counter, m, n, id, Nmodes, big
!  integer, dimension(:), allocatable :: idCC
!
!  counter = 0
!  Nmodes = size(i)
!
!  do n = 1,Nmodes-1
!    do m = n+1,Nmodes
!      if(k1(i(n)) == 0 .and. k2(j(n)) == -k2(j(m)) .and. k3(k(n)) == -k3(k(m))) then
!        counter = counter + 1
!      end if
!    end do
!  end do
!
!  if (counter > 0) then
!    call assert(Nmodes - counter > this%Nwaves)
!    allocate(idCC(counter))
!
!    id = 1
!    do n = 1,size(i)-1
!      do m = n+1,size(i)
!      if(k1(i(n)) == 0 .and. k2(j(n)) == -k2(j(m)) .and. k3(k(n)) == -k3(k(m))) then
!          idCC(id) = m
!          id = id + 1
!        end if
!      end do
!    end do
!    
!    allocate(icpy(Nmodes),jcpy(Nmodes),kcpy(Nmodes))
!    icpy = i
!    jcpy = j
!    kcpy = k
!    deallocate(i,j,k)
!    allocate(i(Nmodes-counter),j(Nmodes-counter),k(Nmodes-counter))
!   
!    call binary_sort(idCC) 
!    call removeItems(icpy,jcpy,kcpy,idCC,i,j,k)
!
!    deallocate(icpy,jcpy,kcpy)
!    deallocate(idCC)
!  end if
!end subroutine

subroutine removeItems(ii,ji,ki,idrmv,io,jo,ko)
  integer, dimension(:), intent(in) :: ii, ji, ki, idrmv
  integer, dimension(:), intent(out) :: io, jo, ko
  integer :: m, n, id

  call assert(size(io) == size(ii) - size(idrmv))
  
  m = 1
  id = 1
  do n = 1,size(ii)
    if (n == idrmv(id)) then
      id = min(id + 1,size(idrmv))
    else
      io(m) = ii(n)
      jo(m) = ji(n)
      ko(m) = ki(n)
      m = m + 1
    end if
  end do
end subroutine

subroutine scrubConjPairs(this,kx,ky,kz)
  class(HIT_shell_forcing), intent(inout) :: this
  integer, dimension(:), intent(inout) :: kx, ky, kz
  integer :: id, m, n, counter
  integer, dimension(:), allocatable :: idCC

  counter = 0
  do n = 1,size(ky) - 1
    do m = n+1,size(ky)
      if (kx(n) == 0 .and. ky(n) == -ky(m) .and. kz(n) == -kz(m)) then
        counter = counter + 1
      end if 
    end do
  end do

  if (counter > 0) then
    allocate(idCC(counter))
    id = 1
    do n = 1,size(ky) - 1
      do m = n+1,size(ky)
        if (kx(n) == 0 .and. ky(n) == -ky(m) .and. kz(n) == -kz(m)) then
          idCC(id) = m
          id = id + 1
        end if 
      end do
    end do
    call this%replaceConjPartners(kx,ky,kz,idCC)
    deallocate(idCC)
  end if
end subroutine

subroutine replaceConjPartners(this,kx,ky,kz,idCC)
  class(HIT_shell_forcing), intent(inout) :: this
  integer, dimension(:), intent(inout) :: kx, ky, kz
  integer, dimension(:), intent(in) :: idCC
  integer :: id
  real(rkind), dimension(size(idCC)) :: tmp, kabs, zeta, theta
  integer, dimension(size(idCC)) :: kxtmp, kytmp, kztmp

  call this%update_seeds()
  call uniform_random(kabs, this%kmin, this%kmax, this%seed1)
  call uniform_random(zeta, -one, one, this%seed2)
  call uniform_random(theta, zero, two*pi, this%seed3)
  
  call convertCylindricalToSphericalWaveVectors(&
    kxtmp, kytmp, kztmp, tmp, kabs, zeta, theta)

  do id = 1,size(idCC)
    kx(idCC(id)) = kxtmp(id)
    ky(idCC(id)) = kytmp(id)
    kz(idCC(id)) = kztmp(id)
  end do
  
  ! Now check that we didn't inadvertantly introduce a new conjugate pair
  call this%scrubConjPairs(kx,ky,kz)
end subroutine

subroutine destroy(this)
   class(HIT_shell_forcing), intent(inout) :: this

   deallocate(this%uhat, this%vhat, this%what)
   nullify(this%fxhat, this%fyhat, this%fzhat, this%cbuffzE)
   if (allocated(this%kabs_sample)) deallocate(this%kabs_sample)
   if (allocated(this%theta_sample)) deallocate(this%theta_sample)
   if (allocated(this%zeta_sample)) deallocate(this%zeta_sample)
   deallocate(this%wave_x, this%wave_y, this%wave_z)
   nullify(this%sp_gpC, this%spectC)
   if (allocated(this%k1)) deallocate(this%k1)
   if (allocated(this%k2)) deallocate(this%k2)
   if (allocated(this%k3)) deallocate(this%k3)
   if (allocated(this%kmag)) deallocate(this%kmag)

   if (allocated(this%i)) deallocate(this%i)
   if (allocated(this%j)) deallocate(this%j)
   if (allocated(this%k)) deallocate(this%k)

   if (allocated(this%fxhatwrite)) deallocate(this%fxhatwrite)
   if (allocated(this%fyhatwrite)) deallocate(this%fyhatwrite)
   if (allocated(this%fzhatwrite)) deallocate(this%fzhatwrite)
   if (allocated(this%cbuffzC)) deallocate(this%cbuffzC)
end subroutine 

subroutine compute_forcing(this,Nmodes)
   class(HIT_shell_forcing), intent(inout) :: this
   integer, intent(in) :: Nmodes
   integer :: ik

   this%fxhat = im0
   this%fyhat = im0
   this%fzhat = im0
   do ik = 1,Nmodes
      call this%embed_forcing_mode(this%wave_x(ik),  this%wave_y(ik),  this%wave_z(ik))
      if (this%wave_x(ik) == 0) then
        call this%embed_forcing_mode(this%wave_x(ik),  -this%wave_y(ik),  -this%wave_z(ik))
      end if
   end do
end subroutine 

subroutine embed_forcing_mode(this, kx, ky, kz)
   class(HIT_shell_forcing), intent(inout) :: this
   integer, intent(in) :: kx, ky, kz
   
   real(rkind) :: den, fac
   integer :: gid_x, gid_y, gid_z
   integer :: lid_x, lid_y
   integer :: ierr

   ! Get global ID of the mode and conjugate
   if (kx >=0 ) then
     gid_x  = kx + 1
   else
     gid_x = this%sp_gpC%xsz(1) ! Since the only negative wavenumber is the oddball
   end if
   if (ky >= 0) then
     gid_y  = ky + 1
   else
     gid_y  = this%sp_gpC%ysz(2) + ky + 1 
   end if
   if (kz >= 0) then
     gid_z  = this%DomAspectRatioZ*kz + 1
   else
     gid_z  = this%sp_gpC%zsz(3) - this%DomAspectRatioZ*(-kz) + 1
   end if 

   ! Get local ID of the mode and conjugate
   lid_x  = gid_x  - this%sp_gpC%zst(1) + 1
   lid_y  = gid_y  - this%sp_gpC%zst(2) + 1
  
   select case (this%version)
   case (1) 
     if ((lid_x >= 1).and.(lid_x <= this%sp_gpC%zsz(1))) then
        if ((lid_y >= 1).and.(lid_y <= this%sp_gpC%zsz(2))) then
           den = abs(this%uhat(lid_x,lid_y,gid_z))**2.d0 + &
                 abs(this%vhat(lid_x,lid_y,gid_z))**2.d0 + &
                 abs(this%what(lid_x,lid_y,gid_z))**2.d0 + 1.d-14 
           
           fac = 0.5d0*this%normfact*this%EpsAmplitude/den/this%Nwaves_rkind
           this%fxhat(lid_x, lid_y, gid_z ) = this%fxhat(lid_x, lid_y, gid_z )&
             + fac*(this%uhat(lid_x, lid_y, gid_z ))
           this%fyhat(lid_x, lid_y, gid_z ) = this%fyhat(lid_x, lid_y, gid_z )&
             + fac*(this%vhat(lid_x, lid_y, gid_z ))
           this%fzhat(lid_x, lid_y, gid_z ) = this%fzhat(lid_x, lid_y, gid_z )&
             + fac*(this%what(lid_x, lid_y, gid_z ))
        end if 
     end if 
   case (2)
     if ((lid_x >= 1).and.(lid_x <= this%sp_gpC%zsz(1))) then
        if ((lid_y >= 1).and.(lid_y <= this%sp_gpC%zsz(2))) then
           den = abs(this%uhat(lid_x,lid_y,gid_z))**2.d0 + &
                 abs(this%vhat(lid_x,lid_y,gid_z))**2.d0 + &
                 abs(this%what(lid_x,lid_y,gid_z))**2.d0 + 1.d-14 
           
           fac = 0.5d0*this%normfact*this%EpsAmplitude/den/this%Nwaves_rkind
           this%fxhat(lid_x, lid_y, gid_z ) = fac*(this%uhat(lid_x, lid_y, gid_z ))
           this%fyhat(lid_x, lid_y, gid_z ) = fac*(this%vhat(lid_x, lid_y, gid_z ))
           this%fzhat(lid_x, lid_y, gid_z ) = fac*(this%what(lid_x, lid_y, gid_z ))
        end if 
     end if 
   case default
     call gracefulExit('Must select version 1 or 2 -- forcingIsotropic.F90',ierr)
   end select

end subroutine 

subroutine getRHS_HITforcing(this, urhs_xy, vrhs_xy, wrhs_xy, uhat_xy, vhat_xy,&
    what_xy, newTimestep, tid)

   class(HIT_shell_forcing), intent(inout) :: this
   
   complex(rkind), dimension(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), &
     this%sp_gpC%ysz(3)), intent(in) :: uhat_xy, vhat_xy
   complex(rkind), dimension(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), &
     this%sp_gpE%ysz(3)), intent(in) :: what_xy
   complex(rkind), dimension(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), &
     this%sp_gpC%ysz(3)), intent(inout) :: urhs_xy, vrhs_xy
   complex(rkind), dimension(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), &
     this%sp_gpE%ysz(3)), intent(inout) :: wrhs_xy
   
   logical, intent(in) :: newTimestep
   integer, intent(in), optional :: tid
   real(rkind) :: alpha_t

!   integer, dimension(this%Nwaves) :: i, j, k ! Indices of the N random modes
   integer :: Nmodes

   ! Run-time checks
   real(rkind) :: EnergyInjection

   if (this%useLinearForcing) then
        this%cbuffyC = uhat_xy
        call this%spectC%hitForceFilter(this%cbuffyC)
        urhs_xy = urhs_xy + (1.d0/this%A_force)*uhat_xy!this%cbuffyC 
       
        this%cbuffyC = vhat_xy
        call this%spectC%hitForceFilter(this%cbuffyC)
        vrhs_xy = vrhs_xy + (1.d0/this%A_force)*vhat_xy!this%cbuffyC
        
        call transpose_y_to_z(what_xy, this%cbuffzE, this%sp_gpE)
        call this%spectC%hitForceFilter_edgeField(this%cbuffzE)
        call transpose_z_to_y(this%cbuffzE, this%cbuffyE, this%sp_gpE)
        wrhs_xy = wrhs_xy + (1.d0/this%A_force)*what_xy!this%cbuffyE
   else


        ! STEP 1: Populate wave_x, wave_y, wave_z

! TODO: Should we pick random wavenumber every sub-step or not?
        if (newTimestep) then
           this%fxhat_old = this%fxhat
           this%fyhat_old = this%fyhat
           this%fzhat_old = this%fzhat
           select case (this%version)
             case (1) 
               call this%pick_random_wavenumbers()
               this%Nmodes = this%Nwaves
             case (2) 
               call this%pick_random_wavenumbersV2(this%Nmodes)
           end select
           call this%update_seeds()
        end if
        
        ! STEP 2: Take FFTz
        call this%takeFFTz(uhat_xy,vhat_xy,what_xy)


        ! STEP 3a: embed into fhat
        call this%compute_forcing(this%Nmodes)
        
        if (newTimeStep .and. this%firstCall) then
            this%fxhat_old = this%fxhat
            this%fyhat_old = this%fyhat
            this%fzhat_old = this%fzhat
        end if 

        ! STEP 3b: Time filter
        this%fxhat = this%alpha_t*this%fxhat + (1.d0 - this%alpha_t)*this%fxhat_old
        this%fyhat = this%alpha_t*this%fyhat + (1.d0 - this%alpha_t)*this%fyhat_old
        this%fzhat = this%alpha_t*this%fzhat + (1.d0 - this%alpha_t)*this%fzhat_old

        ! Store a publicly-available copy of the force at the time of computation
        this%fxhatwrite = this%fxhat
        this%fyhatwrite = this%fyhat
        this%fzhatwrite = this%fzhat
        
        ! On-the-fly checks
        if (this%confirmEnergyInjectionRate) then
          call this%computeEnergyInjectionRate(this%uhat,this%vhat,this%what,&
            this%fxhat, this%fyhat, this%fzhat, EnergyInjection)
          if (abs(this%EpsAmplitude - EnergyInjection) > 1.d-11) then
            call message(1,'Energy injection in forcing routine: ', EnergyInjection)
            call assert(abs(this%EpsAmplitude - EnergyInjection) < 1.d-11,&
              'Energy injection rate is not consistent with input file')
          end if
        end if

        ! STEP 4: Take ifft of fx, fy, fz and add to RHS
        call this%spectC%take_ifft1d_z2z_ip(this%fxhat)
        call transpose_z_to_y(this%fxhat, this%cbuffyC, this%sp_gpC)
        call this%spectC%fix_kx0(this%cbuffyC)
        urhs_xy = urhs_xy + this%cbuffyC

        call this%spectC%take_ifft1d_z2z_ip(this%fyhat)
        call transpose_z_to_y(this%fyhat, this%cbuffyC, this%sp_gpC)
        call this%spectC%fix_kx0(this%cbuffyC)
        vrhs_xy = vrhs_xy + this%cbuffyC

        call this%spectC%shiftz_C2E(this%fzhat)
        call this%spectC%take_ifft1d_z2z_ip(this%fzhat)

        this%cbuffzE(:,:,1:this%sp_gpC%zsz(3)) = this%fzhat
        this%cbuffzE(:,:,this%sp_gpC%zsz(3)+1) = this%cbuffzE(:,:,1)
        call transpose_z_to_y(this%cbuffzE, this%cbuffyE, this%sp_gpE)
        call this%spectE%fix_kx0(this%cbuffyE)
        wrhs_xy = wrhs_xy + this%cbuffyE
    
   end if 
end subroutine

subroutine getIsotropicWaveVectorComponents(k1,k2,k3,n)
  real(rkind), dimension(:,:,:), intent(inout) :: k1, k2, k3
  integer, intent(in) :: n
  integer :: i, j, k, id 

  do i = 1,n/2
    k1(i,:,:) = real(i - 1,rkind)
  end do

  do j = 1,n/2
    k2(:,j,:) = real(j - 1,rkind)
  end do

  id = n/2
  do j = n/2+1,n
    k2(:,j,:) = -k2(:,id,:) - 1.d0
    id = id - 1
  end do

  do k = 1,n/2
    k3(:,:,k) = real(k - 1,rkind)
  end do
  
  id = n/2
  do k = n/2+1,n
    k3(:,:,k) = -k3(:,:,id) - 1.d0
    id = id - 1
  end do

end subroutine

subroutine dumpForcing(this,dumpdir,runID,tid,dumpSpec)
  class(HIT_shell_forcing), intent(inout) :: this
  character(len=*), intent(in) :: dumpdir
  integer, intent(in) :: runID, tid
  logical, intent(in), optional :: dumpSpec
  logical :: dumpSpectral

  dumpSpectral = .false.
  if (present(dumpSpec)) dumpSpectral = dumpSpec

  if (this%storeForce) then
    call this%prepAndDumpField(this%fxhatwrite,dumpdir,runID,tid,'ufrc',.false.)
    call this%prepAndDumpField(this%fyhatwrite,dumpdir,runID,tid,'vfrc',.false.)
    call this%prepAndDumpField(this%fzhatwrite,dumpdir,runID,tid,'wfrc',.false.)
    if (dumpSpectral) then
      call this%prepAndDumpField(this%fxhatwrite,dumpdir,runID,tid,'ufcI',.true.,comp='imag')
      call this%prepAndDumpField(this%fyhatwrite,dumpdir,runID,tid,'vfcI',.true.,comp='imag')
      call this%prepAndDumpField(this%fzhatwrite,dumpdir,runID,tid,'wfcI',.true.,comp='imag')
      
      call this%prepAndDumpField(this%fxhatwrite,dumpdir,runID,tid,'ufcR',.true.,comp='real')
      call this%prepAndDumpField(this%fyhatwrite,dumpdir,runID,tid,'vfcR',.true.,comp='real')
      call this%prepAndDumpField(this%fzhatwrite,dumpdir,runID,tid,'wfcR',.true.,comp='real')
    end if

  end if
end subroutine

subroutine prepAndDumpField(this,arr,dir,runID,tid,desc,dumpSpectral,comp)
  class(HIT_shell_forcing), intent(inout) :: this
  complex(rkind), dimension(this%sp_gpC%zsz(1),this%sp_gpC%zsz(2),&
    this%sp_gpC%zsz(3)), intent(in) :: arr
  character(len=*), intent(in) :: dir
  character(len=4), intent(in) :: desc
  integer, intent(in) :: runID, tid
  logical, intent(in) :: dumpSpectral
  character(len=4), intent(in), optional :: comp
  character(len=clen) :: fname
  character(len=4) :: component
    
  component = 'real'
  if (present(comp)) component = comp

  this%cbuffzC = arr

  if (dumpSpectral) then
    call assert(present(comp),'Must specify which component (real or imag)'//&
     ' to write to disk -- forcingIsotropic.F90')

    write(fname,'(A,I2.2,A1,A4,A2,I6.6,A4)')&
      trim(dir)//'/'//'Run',runID,'_',desc,'_t',tid,'.out'
    
    select case(component)
    case('imag')
      call decomp_2d_write_one(3,dimag(this%cbuffzC),trim(fname),this%sp_gpC)
    case('real')
      call decomp_2d_write_one(3,real(this%cbuffzC,rkind),trim(fname),this%sp_gpC)
    end select
  else
    ! Take 3D ifft
    call this%takeIFFT3(this%cbuffzC,this%rbuffxC)

    ! Write to disk
    write(fname,'(A,I2.2,A1,A4,A2,I6.6,A4)')&
      trim(dir)//'/'//'Run',runID,'_',desc,'_t',tid,'.out'
    call decomp_2d_write_one(1,this%rbuffxC,trim(fname),this%gpC)
  end if
end subroutine

subroutine computeEnergyInjectionRate(this,uhat,vhat,what,&
    fxhat,fyhat,fzhat,EnergyInjection)
  ! uhat, vhat, what, fxhat, fyhat, & fzhat are all 3D FFT cell-centered quantities
  class(HIT_Shell_forcing), intent(inout) :: this
  complex(rkind), dimension(this%sp_gpC%zsz(1),this%sp_gpC%zsz(2), &
    this%sp_gpC%zsz(3)), intent(in) :: uhat, vhat, what, fxhat, fyhat, fzhat
  real(rkind), intent(out) :: EnergyInjection
  complex(rkind) :: EnergyInjectionI
  complex(rkind), dimension(this%sp_gpC%zsz(1),this%sp_gpC%zsz(2), &
    this%sp_gpC%zsz(3)) :: udotf

  udotf = (conjg(uhat)*fxhat + conjg(vhat)*fyhat + &
    conjg(what)*fzhat)/this%normfact

  EnergyInjectionI = p_sum(sum(udotf))

  ! Add complex conjugate components 
  udotf = (uhat*conjg(fxhat) + vhat*conjg(fyhat) + &
    what*conjg(fzhat))/this%normfact
  
  ! Remove the kx=0 contribution soas not to double-count
  if (this%sp_gpC%zst(1) == 1) then ! Corresponds to kx = 0
    udotf(1,:,:) = im0
  end if 

  EnergyInjectionI = EnergyInjectionI + p_sum(sum(udotf))

  ! Verify that EnergyInjection is real valued
  !if (abs(aimag(EnergyInjectionI)) > 1.d-7) then
  !  call warning("The imaginary component of the computed Energy "//&
  !    "injection rate exceeds the threshold of 1.d-7.")
  !  call assert(.false.)
  !end if

  EnergyInjection = real(EnergyInjectionI,rkind)
end subroutine

subroutine takeFFTz(this,uhat_xy,vhat_xy,what_xy)
  class(HIT_Shell_forcing), intent(inout) :: this
  complex(rkind), dimension(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), &
    this%sp_gpC%ysz(3)), intent(in) :: uhat_xy, vhat_xy
  complex(rkind), dimension(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), &
    this%sp_gpE%ysz(3)), intent(in) :: what_xy
  
  call transpose_y_to_z(uhat_xy, this%uhat, this%sp_gpC)
  call this%spectC%take_fft1d_z2z_ip(this%uhat)

  call transpose_y_to_z(vhat_xy, this%vhat, this%sp_gpC)
  call this%spectC%take_fft1d_z2z_ip(this%vhat)

  call transpose_y_to_z(what_xy, this%cbuffzE, this%sp_gpE)
  this%what = this%cbuffzE(:,:,1:this%sp_gpC%zsz(3))
  call this%spectC%take_fft1d_z2z_ip(this%what)
  call this%spectC%shiftz_E2C(this%what)
end subroutine

subroutine takeIFFT3(this,fhat,f)
  class(HIT_shell_forcing), intent(inout) :: this
  complex(rkind), dimension(:,:,:), intent(inout) :: fhat
  real(rkind), dimension(:,:,:), intent(out) :: f
  call this%spectC%take_ifft1d_z2z_ip(fhat)
  call transpose_z_to_y(fhat, this%cbuffyC, this%sp_gpC)
  call this%spectC%ifft(this%cbuffyC, f)
end subroutine

end module 
