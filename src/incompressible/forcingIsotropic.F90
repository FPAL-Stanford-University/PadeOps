module forcingmod
   use kind_parameters, only: rkind
   use decomp_2d, only: decomp_info
   use constants, only: im0

   implicit none
   private
   public :: HIT_shell_forcing

   type :: HIT_shell_forcing
      private
      type(decomp_info), pointer :: sp_gpC
      real(rkind), dimension(:,:,:), allocatable :: k1, k2, k3, kshell
      real(rkind) :: kmin, kmax
      real(rkind) :: Nwaves
      real(rkind) :: EpsAmplitude

      real(rkind), dimension(:), allocatable :: wave_x, wave_y, wave_z
      complex(rkind), dimension(:,:,:), allocatable :: uhat, vhat, what, fxhat, fyhat, fzhat
      contains
      procedure          :: init
      procedure          :: destroy
      procedure, private :: pick_random_wavenumbers
      procedure, private :: embed_forcing
      procedure          :: getRHS_HITforcing
   end type 

contains

subroutine init(this, inputfile, sp_gpC)
   class(HIT_shell_forcing), intent(inout) :: this
   character(len=*), intent(in) :: inputfile
   type(decomp_info), intent(in), target :: sp_gpC
   integer :: nWaves = 20, ierr
   real(rkind) :: kmin = 2.d0, kmax = 10.d0, EpsAmplitude = 0.1d0
   namelist /HIT_Forcing/ kmin, kmax, Nwaves, EpsAmplitude 

   open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=123, NML=HIT_Forcing)
   close(123)

   this%kmin = kmin
   this%kmax = kmax
   this%EpsAmplitude = EpsAmplitude
   this%Nwaves = Nwaves
   this%sp_gpC => sp_gpC

   allocate(this%uhat (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%vhat (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%what (this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%fxhat(this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%fyhat(this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))
   allocate(this%fzhat(this%sp_gpC%zsz(1), this%sp_gpC%zsz(2), this%sp_gpC%zsz(3)))

end subroutine 


subroutine pick_random_wavenumbers(this)
   class(HIT_shell_forcing), intent(inout) :: this

end subroutine 


subroutine destroy(this)
   class(HIT_shell_forcing), intent(inout) :: this

   deallocate(this%uhat, this%vhat, this%what, this%fxhat, this%fyhat, this%fzhat)

   nullify(this%sp_gpC)
end subroutine 


subroutine embed_forcing(this)
   class(HIT_shell_forcing), intent(inout) :: this

!!!!NOTE:::::!!!! 
!   I am assuming wave_x, wave_y and wave_z are in z decomp, and are
!   integers from (1:N) irrespective of domain size.

   integer :: ik, indx, indy, indz
   real(rkind) :: Nwaves_rkind, den, fac

   Nwaves_rkind = real(this%Nwaves, rkind)
   this%fxhat = im0; this%fyhat = im0; this%fzhat = im0
   do ik = 1, size(this%wave_x)
       ! check if kx is on this processor
       indx = this%wave_x(ik) - this%sp_gpC%zst(1) + 1
       if(indx > this%sp_gpC%zsz(1)) then
         ! not on this processor
         cycle
       endif

       ! check if ky is on this processor
       indy = this%wave_y(ik) - this%sp_gpC%zst(2) + 1
       if(indy > this%sp_gpC%zsz(2)) then
         ! not on this processor
         cycle
       endif

       ! kz must be on this processor because we are in zdecomp
       indz = this%wave_z(ik) - this%sp_gpC%zst(3) + 1

      ! now indx, indy, indz are indices for this wavenumber
      den = abs(this%uhat(indx, indy, indz))**2 + &
            abs(this%vhat(indx, indy, indz))**2 + &
            abs(this%what(indx, indy, indz))**2
      fac = this%EpsAmplitude/den/Nwaves_rkind
      this%fxhat(indx, indy, indz) = fac*conjg(this%uhat(indx,indy,indz))
      this%fyhat(indx, indy, indz) = fac*conjg(this%vhat(indx,indy,indz))
      this%fzhat(indx, indy, indz) = fac*conjg(this%what(indx,indy,indz))
   enddo

end subroutine 


subroutine getRHS_HITforcing(this, urhs_xy, vrhs_xy, wrhs_xy, uhat_xy, vhat_xy, what_xy)
   class(HIT_shell_forcing), intent(inout) :: this
   complex(rkind), dimension(:,:,:), intent(in) :: uhat_xy, vhat_xy, what_xy
   complex(rkind), dimension(:,:,:), intent(inout) :: urhs_xy, vrhs_xy, wrhs_xy
   
   integer :: idx



   ! STEP 1: Populate wave_x, wave_y, wave_z


   ! STEP 2: embed into fhat
   call this%embed_forcing()
   ! transform from fhat -> fhat_xy
   ! then add to urhs



end subroutine




end module 
