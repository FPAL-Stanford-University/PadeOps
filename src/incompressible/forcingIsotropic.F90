module forcingmod
   use kind_parameters, only: rkind, mpirkind
   use decomp_2d
   use constants,       only: im0, one, zero, two, pi
   use spectralMod,     only: spectral 
   use exits,           only: GracefulExit, message
   use mpi
   use fortran_assert,  only: assert 
   use random,          only: uniform_random, randperm
   use arrayTools,         only: findGL 
   use gridtools,       only: loc2glob, glob2loc

   implicit none
   private
   public :: HIT_shell_forcing

   type :: HIT_shell_forcing
      private
      type(decomp_info), pointer :: sp_gpC, sp_gpE
      real(rkind) :: kmin, kmax
      integer(rkind) :: Nwaves
      real(rkind) :: EpsAmplitude
      class(spectral), pointer :: spectC

      integer, dimension(:), allocatable :: wave_x, wave_y, wave_z
      complex(rkind), dimension(:,:,:), allocatable :: uhat, vhat, what, fxhat_old, fyhat_old, fzhat_old
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
      integer, dimension(:), allocatable :: i, j, k
    
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
      procedure, private :: scrubConjPairsV2
      procedure, private :: replaceConjPartners
   end type 

contains

subroutine init(this, inputfile, sp_gpC, sp_gpE, spectC, cbuffyE, cbuffyC, cbuffzE, cbuffzC, tidStart)
   class(HIT_shell_forcing), intent(inout) :: this
   character(len=*), intent(in) :: inputfile
   type(decomp_info), intent(in), target :: sp_gpC, sp_gpE
   integer, intent(in) :: tidStart
   complex(rkind), dimension(:,:,:  ), intent(inout), target :: cbuffzE, cbuffyE, cbuffyC
   complex(rkind), dimension(:,:,:,:), intent(inout), target :: cbuffzC
   class(spectral), intent(in), target :: spectC
   integer :: RandSeedToAdd = 0, ierr, DomAspectRatioZ = 1
   real(rkind) :: alpha_t = 1.d0 
   integer :: Nwaves = 20
   real(rkind) :: kmin = 2.d0, kmax = 10.d0, EpsAmplitude = 0.1d0, A_force = 1.d0 
   logical :: useLinearForcing = .false. 
   real(rkind) :: filtfact_linForcing = 0.5d0
   integer :: nforce
   namelist /HIT_Forcing/ kmin, kmax, Nwaves, EpsAmplitude, RandSeedToAdd, DomAspectRatioZ, alpha_t, useLinearForcing, filtfact_linForcing  

   open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=123, NML=HIT_Forcing)
   close(123)

   if(DomAspectRatioZ < 1.0d0) then
       call GracefulExit("Aspect ratio in z must be greater than 1", 111)
   endif

   this%A_force = A_force
   this%kmin = kmin
   this%kmax = kmax
   this%EpsAmplitude = EpsAmplitude
   this%Nwaves = Nwaves
   this%DomAspectRatioZ = DomAspectRatioZ
   this%sp_gpC => sp_gpC
   this%sp_gpE => sp_gpE
   this%spectC => spectC
   this%useLinearForcing = useLinearForcing
   allocate(this%uhat (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%vhat (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%what (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%fxhat_old (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%fyhat_old (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%fzhat_old (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
  
   ! Define the wave-vector grid
   nforce = ceiling(this%kmax)*2
   if (mod(ceiling(this%kmax),2) .ne. 0) nforce = ceiling(this%kmax + 1.d0)*2

   allocate(this%k1  (nforce/2+1,nforce,nforce))
   allocate(this%k2  (nforce/2+1,nforce,nforce))
   allocate(this%k3  (nforce/2+1,nforce,nforce))
   allocate(this%kmag(nforce/2+1,nforce,nforce))
   call getIsotropicWaveVectorComponents(this%k1,this%k2,this%k3,nforce) 
   this%kmag = sqrt(this%k1*this%k1 + this%k2*this%k2 + this%k3*this%k3)

   ! Find the indices of each mode that satisfies kmin <= k <= kmax
   call findGL(this%kmag,this%kmin,this%kmax,this%i,this%j,this%k)
  
   ! Remove complex conjugate pairs
   call this%scrubConjPairsV2(this%i,this%j,this%k,&
     this%k1(:,1,1),this%k2(1,:,1),this%k3(1,1,:))
  
   this%firstCall = .true.
   this%fxhat    => cbuffzC(:,:,:,1)
   this%fyhat    => cbuffzC(:,:,:,2)
   this%fzhat    => cbuffzC(:,:,:,3)
   this%cbuffzE  => cbuffzE 
   this%cbuffyE  => cbuffyE
   this%cbuffyC  => cbuffyC

   this%seed0 = tidStart + RandSeedToAdd
   call this%update_seeds()
  
   this%normfact = (real(spectC%physdecomp%xsz(1),rkind)*real(spectC%physdecomp%ysz(2),rkind)*real(spectC%physdecomp%zsz(3),rkind))**2
  
   allocate(this%kabs_sample(Nwaves))
   allocate(this%theta_sample(Nwaves))
   allocate(this%zeta_sample(Nwaves))
   allocate(this%tmpModes(Nwaves))

   allocate(this%wave_x(Nwaves))
   allocate(this%wave_y(Nwaves))
   allocate(this%wave_z(Nwaves))

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

subroutine pick_random_wavenumbersV2(this)
  class(HIT_shell_forcing), intent(inout) :: this
  integer :: n, ierr
  integer, dimension(this%Nwaves) :: i, j, k
  
  ! Randomly sort the arrays
  call randomShuffle(this%i,this%j,this%k,this%seed3)

  ! Choose the first N modes
  i = this%i(1:this%Nwaves)
  j = this%j(1:this%Nwaves)
  k = this%k(1:this%Nwaves)

  do n = 1,this%Nwaves
    this%wave_x(n) = this%k1(i(n),j(n),k(n)) 
    this%wave_y(n) = this%k2(i(n),j(n),k(n)) 
    this%wave_z(n) = this%k3(i(n),j(n),k(n)) 
  end do

end subroutine

subroutine randomShuffle(i,j,k,seed)
  integer, dimension(:), intent(inout) :: i, j, k
  integer, intent(in) :: seed
  integer, dimension(size(i)) :: idx

  call randperm(size(i),seed,idx)

  i = i(idx)
  j = j(idx)
  k = k(idx)
  
end subroutine

subroutine scrubConjPairsV2(this,i,j,k,k1,k2,k3)
  use sorting_mod, only: binary_sort
  ! Ensure none of the wavevectors indexed by i, j, and k correspond to
  ! complex conjugate (C.C.) modes

  class(HIT_shell_forcing), intent(inout) :: this
  integer, dimension(:), allocatable, intent(inout) :: i, j, k
  real(rkind), dimension(:), intent(in) :: k1, k2, k3
  integer, dimension(:), allocatable :: icpy, jcpy, kcpy
  integer :: counter, m, n, id, Nmodes, big
  integer, dimension(:), allocatable :: idCC

  counter = 0
  Nmodes = size(i)

  do n = 1,Nmodes-1
    do m = n+1,Nmodes
      if(k1(i(n)) == 0 .and. k2(j(n)) == -k2(j(m)) .and. k3(k(n)) == -k3(k(m))) then
        counter = counter + 1
      end if
    end do
  end do

  if (counter > 0) then
    allocate(idCC(counter))

    id = 1
    do n = 1,size(i)-1
      do m = n+1,size(i)
      if(k1(i(n)) == 0 .and. k2(j(n)) == -k2(j(m)) .and. k3(k(n)) == -k3(k(m))) then
          idCC(id) = m
          id = id + 1
        end if
      end do
    end do
    
    allocate(icpy(Nmodes),jcpy(Nmodes),kcpy(Nmodes))
    icpy = i
    jcpy = j
    kcpy = k
    deallocate(i,j,k)
    allocate(i(Nmodes-counter),j(Nmodes-counter),k(Nmodes-counter))
   
    call binary_sort(idCC) 
    call removeItems(icpy,jcpy,kcpy,idCC,i,j,k)

    deallocate(icpy,jcpy,kcpy)
    deallocate(idCC)
  end if
end subroutine

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
   if (nrank == 0) then
      deallocate(this%kabs_sample, this%theta_sample, this%zeta_sample)
   end if
   deallocate(this%wave_x, this%wave_y, this%wave_z)
   nullify(this%sp_gpC, this%spectC)
   if (allocated(this%k1)) deallocate(this%k1)
   if (allocated(this%k2)) deallocate(this%k2)
   if (allocated(this%k3)) deallocate(this%k3)
   if (allocated(this%kmag)) deallocate(this%kmag)

   if (allocated(this%i)) deallocate(this%i)
   if (allocated(this%j)) deallocate(this%j)
   if (allocated(this%k)) deallocate(this%k)
end subroutine 

subroutine compute_forcing(this)
   class(HIT_shell_forcing), intent(inout) :: this
   integer :: ik

   this%fxhat = im0
   this%fyhat = im0
   this%fzhat = im0
   do ik = 1,this%Nwaves
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
   integer :: gid_x, gid_y, gid_yC, gid_z, gid_zC
   integer :: lid_x, lid_y, lid_yC

   ! Get global ID of the mode and conjugate
   gid_x  = kx + 1
   if (ky >= 0) then
     gid_y  = ky + 1
     gid_yC = this%sp_gpC%ysz(2) - ky + 1 
   else
     gid_yC = -ky + 1
     gid_y  = this%sp_gpC%ysz(2) + ky + 1 
   end if
   if (kz >= 0) then
     gid_z  = this%DomAspectRatioZ*kz + 1
     gid_zC = this%sp_gpC%zsz(3) - this%DomAspectRatioZ*kz + 1
   else
     gid_zC = this%DomAspectRatioZ*(-kz) + 1
     gid_z  = this%sp_gpC%zsz(3) - this%DomAspectRatioZ*(-kz) + 1
   end if 

   ! Get local ID of the mode and conjugate
   lid_x  = gid_x  - this%sp_gpC%zst(1) + 1
   lid_y  = gid_y  - this%sp_gpC%zst(2) + 1
   lid_yC = gid_yC - this%sp_gpC%zst(2) + 1
   

   if ((lid_x >= 1).and.(lid_x <= this%sp_gpC%zsz(1))) then
      if ((lid_y >= 1).and.(lid_y <= this%sp_gpC%zsz(2))) then
         den = abs(this%uhat(lid_x,lid_y,gid_z))**2 + &
               abs(this%vhat(lid_x,lid_y,gid_z))**2 + &
               abs(this%what(lid_x,lid_y,gid_z))**2 + 1.d-14
         
         fac = 0.5d0*this%normfact*this%EpsAmplitude/den/this%Nwaves_rkind
         this%fxhat(lid_x, lid_y, gid_z ) = this%fxhat(lid_x, lid_y, gid_z )&
           + fac*this%uhat(lid_x, lid_y, gid_z )
         this%fyhat(lid_x, lid_y, gid_z ) = this%fyhat(lid_x, lid_y, gid_z )&
           + fac*this%vhat(lid_x, lid_y, gid_z )
         this%fzhat(lid_x, lid_y, gid_z ) = this%fzhat(lid_x, lid_y, gid_z )&
           + fac*this%what(lid_x, lid_y, gid_z )

      end if 
   end if 

end subroutine 

subroutine getRHS_HITforcing(this, urhs_xy, vrhs_xy, wrhs_xy, uhat_xy, vhat_xy, what_xy, newTimestep)
   class(HIT_shell_forcing), intent(inout) :: this
   
   complex(rkind), dimension(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), this%sp_gpC%ysz(3)), intent(in)    :: uhat_xy, vhat_xy
   complex(rkind), dimension(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), this%sp_gpE%ysz(3)), intent(in)    :: what_xy
   complex(rkind), dimension(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), this%sp_gpC%ysz(3)), intent(inout) :: urhs_xy, vrhs_xy
   complex(rkind), dimension(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), this%sp_gpE%ysz(3)), intent(inout) :: wrhs_xy
   
   logical, intent(in) :: newTimestep
   real(rkind) :: alpha_t

   integer, dimension(this%Nwaves) :: i, j, k ! Indices of the N random modes

   if (this%useLinearForcing) then
        this%cbuffyC = uhat_xy
        call this%spectC%hitForceFilter(this%cbuffyC)
        urhs_xy = urhs_xy + (1.d0/this%A_force)*this%cbuffyC 
       
        this%cbuffyC = vhat_xy
        call this%spectC%hitForceFilter(this%cbuffyC)
        vrhs_xy = vrhs_xy + (1.d0/this%A_force)*this%cbuffyC
        
        call transpose_y_to_z(what_xy, this%cbuffzE, this%sp_gpE)
        call this%spectC%hitForceFilter_edgeField(this%cbuffzE)
        call transpose_z_to_y(this%cbuffzE, this%cbuffyE, this%sp_gpE)
        wrhs_xy = wrhs_xy + (1.d0/this%A_force)*this%cbuffyE
   else


        ! STEP 1: Populate wave_x, wave_y, wave_z
        
        if (newTimestep) then
           this%fxhat_old = this%fxhat
           this%fyhat_old = this%fyhat
           this%fzhat_old = this%fzhat
           !call this%pick_random_wavenumbers()
           call this%pick_random_wavenumbersV2()
           call this%update_seeds()
        end if


        ! STEP 2: Take FFTz
        call transpose_y_to_z(uhat_xy, this%uhat, this%sp_gpC)
        call this%spectC%take_fft1d_z2z_ip(this%uhat)

        call transpose_y_to_z(vhat_xy, this%vhat, this%sp_gpC)
        call this%spectC%take_fft1d_z2z_ip(this%vhat)

        call transpose_y_to_z(what_xy, this%cbuffzE, this%sp_gpE)
        this%what = this%cbuffzE(:,:,1:this%sp_gpC%zsz(3))
        call this%spectC%take_fft1d_z2z_ip(this%what)
        call this%spectC%shiftz_E2C(this%what)


        ! STEP 3a: embed into fhat
        call this%compute_forcing()
        
        if (newTimeStep .and. this%firstCall) then
            this%fxhat_old = this%fxhat
            this%fyhat_old = this%fyhat
            this%fzhat_old = this%fzhat
        end if 

        ! STEP 3b: Time filter
        this%fxhat = this%alpha_t*this%fxhat + (1.d0 - this%alpha_t)*this%fxhat_old
        this%fyhat = this%alpha_t*this%fyhat + (1.d0 - this%alpha_t)*this%fyhat_old
        this%fzhat = this%alpha_t*this%fzhat + (1.d0 - this%alpha_t)*this%fzhat_old

        ! STEP 4: Take ifft of fx, fy, fz and add to RHS
        call this%spectC%take_ifft1d_z2z_ip(this%fxhat)
        call transpose_z_to_y(this%fxhat, this%cbuffyC, this%sp_gpC)
        urhs_xy = urhs_xy + this%cbuffyC

        call this%spectC%take_ifft1d_z2z_ip(this%fyhat)
        call transpose_z_to_y(this%fyhat, this%cbuffyC, this%sp_gpC)
        vrhs_xy = vrhs_xy + this%cbuffyC

        call this%spectC%shiftz_C2E(this%fzhat)
        call this%spectC%take_ifft1d_z2z_ip(this%fzhat)
        this%cbuffzE(:,:,1:this%sp_gpC%zsz(3)) = this%fzhat
        this%cbuffzE(:,:,this%sp_gpC%zsz(3)+1) = this%cbuffzE(:,:,1)
        call transpose_z_to_y(this%cbuffzE, this%cbuffyE, this%sp_gpE)
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
  k1(n/2+1,:,:) = real(-n/2,rkind)

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

end module 
