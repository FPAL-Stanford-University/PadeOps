module forcingmod
   use kind_parameters, only: rkind
   use decomp_2d
   use constants, only: im0, one, zero, two, pi
   use spectralMod, only: spectral 
   use mpi 

   implicit none
   private
   public :: HIT_shell_forcing

   type :: HIT_shell_forcing
      private
      type(decomp_info), pointer :: sp_gpC, sp_gpE
      real(rkind) :: kmin, kmax
      real(rkind) :: Nwaves
      real(rkind) :: EpsAmplitude
      class(spectral), pointer :: spectC

      integer, dimension(:), allocatable :: wave_x, wave_y, wave_z
      complex(rkind), dimension(:,:,:), allocatable :: uhat, vhat, what
      complex(rkind), dimension(:,:,:), pointer :: fxhat, fyhat, fzhat, cbuffzE, cbuffyE, cbuffyC
      integer, dimension(:), allocatable :: k1inZ, k2inZ, k3inZ

      real(rkind), dimension(:), allocatable :: kabs_sample, zeta_sample, theta_sample
      integer :: seed0, seed1, seed2, seed3
      integer :: myk1min, myk1max, myk2min, myk2max
      real(rkind) :: Nwaves_rkind
      real(rkind), dimension(:), allocatable :: tmpModes

      real(rkind) :: normfact = 1.d0

      contains
      procedure          :: init
      procedure          :: destroy
      procedure, private :: update_seeds
      procedure, private :: pick_random_wavenumbers
      procedure, private :: compute_forcing
      procedure, private :: embed_forcing_mode
      procedure          :: getRHS_HITforcing
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
   real(rkind), dimension(:,:,:), allocatable :: rbuffzC
   integer :: RandSeedToAdd = 0

   integer :: nWaves = 20, ierr
   real(rkind) :: kmin = 2.d0, kmax = 10.d0, EpsAmplitude = 0.1d0
   namelist /HIT_Forcing/ kmin, kmax, Nwaves, EpsAmplitude, RandSeedToAdd 

   open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=123, NML=HIT_Forcing)
   close(123)

   this%kmin = kmin
   this%kmax = kmax
   this%EpsAmplitude = EpsAmplitude
   this%Nwaves = Nwaves
   this%sp_gpC => sp_gpC
   this%sp_gpE => sp_gpE
   this%spectC => spectC

   allocate(this%uhat (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%vhat (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%what (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   
   this%fxhat    => cbuffzC(:,:,:,1)
   this%fyhat    => cbuffzC(:,:,:,2)
   this%fzhat    => cbuffzC(:,:,:,3)
   this%cbuffzE  => cbuffzE 
   this%cbuffyE  => cbuffyE
   this%cbuffyC  => cbuffyC

   this%seed0 = tidStart + RandSeedToAdd
   call this%update_seeds()
  
   this%normfact = real(spectC%physdecomp%xsz(1),rkind)*real(spectC%physdecomp%ysz(2),rkind)*real(spectC%physdecomp%zsz(3),rkind)
   
   allocate(this%kabs_sample(Nwaves))
   allocate(this%theta_sample(Nwaves))
   allocate(this%zeta_sample(Nwaves))
   allocate(this%tmpModes(Nwaves))

   allocate(this%wave_x(Nwaves))
   allocate(this%wave_y(Nwaves))
   allocate(this%wave_z(Nwaves))

   allocate(rbuffzC(sp_gpC%zsz(1), sp_gpC%zsz(2), sp_gpC%zsz(3)))
   allocate(this%k1inZ(sp_gpC%zsz(1)))
   allocate(this%k2inZ(sp_gpC%zsz(2)))
   allocate(this%k3inZ(sp_gpC%zsz(3)))

   call transpose_y_to_z(spectC%k1, rbuffzC, sp_gpC)
   this%k1inZ = nint(rbuffzC(:,1,1)) 
   
   call transpose_y_to_z(spectC%k2, rbuffzC, sp_gpC)
   this%k2inZ = nint(rbuffzC(1,:,1)) 
  
   call transpose_y_to_z(spectC%k3, rbuffzC, sp_gpC)
   this%k3inZ = nint(rbuffzC(1,1,:))

   this%myk1min = minval(this%k1inZ)
   this%myk2min = minval(this%k2inZ)
   
   this%myk1max = maxval(this%k1inZ)
   this%myk2max = maxval(this%k2inZ)

   deallocate(rbuffzC)

   this%Nwaves_rkind = real(this%Nwaves, rkind)

end subroutine 

subroutine update_seeds(this)
   class(HIT_shell_forcing), intent(inout) :: this
   this%seed0 = this%seed0 + 1
   this%seed1 = this%seed0*12 + mod(this%seed0,83)*123 + mod(this%seed0,1196)*1234 + mod(14326455,this%seed0)
   this%seed2 = this%seed0*21 + mod(this%seed0,41)*321 + mod(this%seed0,2434)*4321 + mod(87646314,this%seed1)
   this%seed3 = this%seed0*16 + mod(this%seed0,17)*256 + mod(this%seed0,1353)*4096 + mod(34844875,this%seed2)
end subroutine


subroutine pick_random_wavenumbers(this)
   use random, only: uniform_random 
   class(HIT_shell_forcing), intent(inout) :: this
   integer :: ierr

   call uniform_random(this%kabs_sample, this%kmin, this%kmax, this%seed1)
   call uniform_random(this%zeta_sample, -one, one, this%seed2)
   call uniform_random(this%theta_sample, zero, two*pi, this%seed3)
   
   this%tmpModes = this%kabs_sample*sqrt(1 - this%zeta_sample**2)*cos(this%theta_sample)
   where(this%tmpModes < 0) this%tmpModes = -this%tmpModes
   this%wave_x = ceiling(this%tmpModes)
   
   this%tmpModes = this%kabs_sample*sqrt(1 - this%zeta_sample**2)*sin(this%theta_sample)
   where(this%tmpModes < 0) this%tmpModes = -this%tmpModes
   this%wave_y = ceiling(this%tmpModes)
   
   this%tmpModes = this%kabs_sample*this%zeta_sample
   where(this%tmpModes < 0) this%tmpModes = -this%tmpModes
   this%wave_z = ceiling(this%tmpModes)


end subroutine 


subroutine destroy(this)
   class(HIT_shell_forcing), intent(inout) :: this

   deallocate(this%uhat, this%vhat, this%what)
   nullify(this%fxhat, this%fyhat, this%fzhat, this%cbuffzE)
   deallocate(this%k1inZ, this%k2inZ, this%k3inZ)
   if (nrank == 0) then
      deallocate(this%kabs_sample, this%theta_sample, this%zeta_sample)
   end if
   deallocate(this%wave_x, this%wave_y, this%wave_z)
   nullify(this%sp_gpC, this%spectC)
end subroutine 


!subroutine embed_forcing(this)
!   class(HIT_shell_forcing), intent(inout) :: this
!
!!!!!NOTE:::::!!!! 
!!   I am assuming wave_x, wave_y and wave_z are in z decomp, and are
!!   integers from (1:N) irrespective of domain size.
!
!   integer :: ik, indx, indy, indz
!   real(rkind) :: Nwaves_rkind, den, fac
!
!   Nwaves_rkind = real(this%Nwaves, rkind)
!   this%fxhat = im0; this%fyhat = im0; this%fzhat = im0
!   do ik = 1, size(this%wave_x)
!       ! check if kx is on this processor
!       indx = this%wave_x(ik) - this%sp_gpC%zst(1) + 1
!       if(indx > this%sp_gpC%zsz(1)) then
!         ! not on this processor
!         cycle
!       endif
!
!       ! check if ky is on this processor
!       indy = this%wave_y(ik) - this%sp_gpC%zst(2) + 1
!       if(indy > this%sp_gpC%zsz(2)) then
!         ! not on this processor
!         cycle
!       endif
!
!       ! kz must be on this processor because we are in zdecomp
!       indz = this%wave_z(ik) - this%sp_gpC%zst(3) + 1
!
!      ! now indx, indy, indz are indices for this wavenumber
!      den = abs(this%uhat(indx, indy, indz))**2 + &
!            abs(this%vhat(indx, indy, indz))**2 + &
!            abs(this%what(indx, indy, indz))**2
!      fac = this%EpsAmplitude/den/Nwaves_rkind
!      this%fxhat(indx, indy, indz) = this%fxhat(indx, indy, indz) + fac*conjg(this%uhat(indx,indy,indz))
!      this%fyhat(indx, indy, indz) = this%fyhat(indx, indy, indz) + fac*conjg(this%vhat(indx,indy,indz))
!      this%fzhat(indx, indy, indz) = this%fzhat(indx, indy, indz) + fac*conjg(this%what(indx,indy,indz))
!   enddo
!
!end subroutine 


subroutine compute_forcing(this)
   class(HIT_shell_forcing), intent(inout) :: this
   integer :: ik, ierr

   this%fxhat = im0
   this%fyhat = im0
   this%fzhat = im0
   do ik = 1,this%Nwaves
      call this%embed_forcing_mode(this%wave_x(ik),  this%wave_y(ik),  this%wave_z(ik))
   end do 
end subroutine 

subroutine embed_forcing_mode(this, kx, ky, kz)
   class(HIT_shell_forcing), intent(inout) :: this
   integer, intent(in) :: kx, ky, kz
   
   real(rkind) :: den, fac
   integer :: gid_x, gid_y, gid_yC, gid_z, gid_zC
   integer :: lid_x, lid_y, lid_yC, ierr

   ! Get global ID of the mode and conjugate
   gid_x  = kx + 1
   gid_y  = ky + 1
   gid_yC = this%sp_gpC%ysz(2) - ky + 1 
   gid_z  = kz + 1
   gid_zC = this%sp_gpC%zsz(3) - kz + 1 

   ! Get local ID of the mode and conjugate
   lid_x  = gid_x  - this%sp_gpC%zst(1) + 1
   lid_y  = gid_y  - this%sp_gpC%zst(2) + 1
   lid_yC = gid_yC - this%sp_gpC%zst(2) + 1
   

   if ((lid_x >= 1).and.(lid_x <= this%sp_gpC%zsz(1))) then
      if ((lid_y >= 1).and.(lid_y <= this%sp_gpC%zsz(2))) then
         den = abs(this%uhat(lid_x,lid_y,gid_z))**2 + &
               abs(this%vhat(lid_x,lid_y,gid_z))**2 + &
               abs(this%what(lid_x,lid_y,gid_z))**2 + 1.d-14
         
         fac = this%normfact*this%EpsAmplitude/den/this%Nwaves_rkind
         this%fxhat(lid_x, lid_y, gid_z ) = this%fxhat(lid_x, lid_y, gid_z ) + fac*conjg(this%uhat(lid_x, lid_y, gid_z ))
         this%fxhat(lid_x, lid_y, gid_zC) = this%fxhat(lid_x, lid_y, gid_zC) + fac*conjg(this%uhat(lid_x, lid_y, gid_zC))
                                            
         this%fyhat(lid_x, lid_y, gid_z ) = this%fyhat(lid_x, lid_y, gid_z ) + fac*conjg(this%vhat(lid_x, lid_y, gid_z ))
         this%fyhat(lid_x, lid_y, gid_zC) = this%fyhat(lid_x, lid_y, gid_zC) + fac*conjg(this%vhat(lid_x, lid_y, gid_zC))
                                            
         this%fzhat(lid_x, lid_y, gid_z ) = this%fzhat(lid_x, lid_y, gid_z ) + fac*conjg(this%what(lid_x, lid_y, gid_z ))
         this%fzhat(lid_x, lid_y, gid_zC) = this%fzhat(lid_x, lid_y, gid_zC) + fac*conjg(this%what(lid_x, lid_y, gid_zC))
         
      end if

      if ((lid_yC >= 1).and.(lid_yC <= this%sp_gpC%zsz(2))) then
         den = abs(this%uhat(lid_x,lid_yC,gid_z))**2 + &
               abs(this%vhat(lid_x,lid_yC,gid_z))**2 + &
               abs(this%what(lid_x,lid_yC,gid_z))**2 + 1.d-14
         
         fac = this%normfact*this%EpsAmplitude/den/this%Nwaves_rkind
         this%fxhat(lid_x, lid_yC, gid_z ) = this%fxhat(lid_x, lid_yC, gid_z ) + fac*conjg(this%uhat(lid_x, lid_yC, gid_z ))
         this%fxhat(lid_x, lid_yC, gid_zC) = this%fxhat(lid_x, lid_yC, gid_zC) + fac*conjg(this%uhat(lid_x, lid_yC, gid_zC))
                                             
         this%fyhat(lid_x, lid_yC, gid_z ) = this%fyhat(lid_x, lid_yC, gid_z ) + fac*conjg(this%vhat(lid_x, lid_yC, gid_z ))
         this%fyhat(lid_x, lid_yC, gid_zC) = this%fyhat(lid_x, lid_yC, gid_zC) + fac*conjg(this%vhat(lid_x, lid_yC, gid_zC))
                                             
         this%fzhat(lid_x, lid_yC, gid_z ) = this%fzhat(lid_x, lid_yC, gid_z ) + fac*conjg(this%what(lid_x, lid_yC, gid_z ))
         this%fzhat(lid_x, lid_yC, gid_zC) = this%fzhat(lid_x, lid_yC, gid_zC) + fac*conjg(this%what(lid_x, lid_yC, gid_zC))
         
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


   ! STEP 1: Populate wave_x, wave_y, wave_z
   if (newTimestep) then
      call this%pick_random_wavenumbers()
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


   ! STEP 3: embed into fhat
   call this%compute_forcing()

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


end subroutine




end module 
