module forcingmod
   use kind_parameters, only: rkind
   use decomp_2d, only: decomp_info

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
      complex(rkind), dimension(:,:,:), allocatable :: uhat, vhat, what, fhat
      contains
      procedure          :: init
      procedure          :: destroy
      procedure, private :: pick_random_wavenumbers
      procedure, private :: embed_forcing
      procedure          :: getRHS_HITforcing
   end type 

contains

subroutine init(this, inputfile, sp_gp)
   class(HIT_shell_forcing), intent(inout) :: this
   character(len=*), intent(in) :: inputfile
   type(decomp_info), intent(in), target :: sp_gp
   integer :: nWaves = 20
   real(rkind) :: kmin = 2.d0, kmax = 10.d0, EpsAmplitude = 0.1d0
   namelist /HIT_Forcing/ kmin, kmax, Nwaves, EpsAmplitude 

   open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=123, NML=HIT_Forcing)
   close(123)

   this%kmin = kmin
   this%kmax = kmax
   this%EpsAmplitude = EpsAmplitude
   this%Nwaves = Nwaves
   this%sp_gp => sp_gp

end subroutine 


subroutine destroy(this)
   class(HIT_shell_forcing), intent(inout) :: this

   nullify(this%sp_gpC)
end subroutine 

end module



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
