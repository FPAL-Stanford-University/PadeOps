module PadeDerOps
   use kind_parameters, only: rkind
   use cd06staggstuff, only: cd06stagg
   use decomp_2d
   use staggOpsMod, only: staggOps
   use constants, only: zero, one, two, three, five
   use exits, only: gracefulExit, message
   use spectralMod, only: spectral

   implicit none
   
   private

   public :: Pade6stagg, fd02, cd06, fourierColl

   integer, parameter :: fd02 = 0
   integer, parameter :: cd06 = 1
   integer, parameter :: fourierColl = 2
   complex(rkind), parameter :: zeroC = dcmplx(zero,zero)
   !logical :: isPeriodic = .false. 

   type Pade6stagg
      type(cd06stagg), allocatable :: derOO, derEE, derOE, derEO, derOS, derSO, derSE, derES, derSS
      type(cd06stagg), allocatable :: derPeriodic
      type(staggOps) :: fd02_ss, fd02_ns, fd02_sn, fd02_nn, fd02_periodic
      type(decomp_info), pointer :: gp, sp_gp
      type(spectral), pointer :: spectC
      integer :: scheme = 1
      real(rkind) :: dz
      logical, public :: isPeriodic = .false. 
      contains
      procedure          :: init
      procedure          :: destroy
      procedure          :: getModifiedWavenumbers
      procedure          :: getApproxPoincareConstant 
      procedure, private :: ddz_C2E_real
      procedure, private :: ddz_C2E_cmplx
      procedure, private :: ddz_E2C_real
      procedure, private :: ddz_E2C_cmplx
      procedure, private :: interpz_C2E_real
      procedure, private :: interpz_C2E_cmplx
      procedure, private :: interpz_E2C_real
      procedure, private :: interpz_E2C_cmplx
      procedure, private :: d2dz2_C2C_real
      procedure, private :: d2dz2_C2C_cmplx
      procedure, private :: d2dz2_E2E_real
      procedure, private :: d2dz2_E2E_cmplx
      generic            :: ddz_C2E => ddz_C2E_real, ddz_C2E_cmplx
      generic            :: ddz_E2C => ddz_E2C_real, ddz_E2C_cmplx
      procedure          :: ddz_C2C
      generic            :: d2dz2_C2C => d2dz2_C2C_real, d2dz2_C2C_cmplx
      generic            :: d2dz2_E2E => d2dz2_E2E_real, d2dz2_E2E_cmplx
      generic            :: interpz_C2E => interpz_C2E_real, interpz_C2E_cmplx
      generic            :: interpz_E2C => interpz_E2C_real, interpz_E2C_cmplx
      procedure          :: ddz_1d_C2C
      procedure          :: ddz_1d_C2E
      procedure          :: interp_1d_E2C
  end type

contains

subroutine init(this, gpC, sp_gpC, gpE, sp_gpE, dz, scheme, isPeriodic, spectC)
   class(Pade6stagg), intent(out) :: this
   type(decomp_info), intent(in), target :: gpC, sp_gpC, gpE, sp_gpE
   integer, intent(in) :: scheme 
   real(rkind), intent(in) :: dz
   logical, intent(in) :: isPeriodic 
   type(spectral), intent(in), target, optional :: spectC

   this%gp => gpC
   this%sp_gp => sp_gpC
   this%scheme = scheme
   this%dz = dz

   this%isPeriodic = isPeriodic

   if (this%isPeriodic) then
      select case(this%scheme) 
      case(fourierColl)
         if (present(spectC)) then
            this%spectC => spectC
            call message(0,"Fourier collocation successfully initialized in Z")
         else
            call GracefulExit("You need to pass in a spectral derived type if you want to use Fourier differentiation in z", 43)
         end if
      case(cd06)
         allocate(this%derPeriodic)
         call this%derPeriodic%init(this%gp%zsz(3),dz)
      case(fd02)
         call this%fd02_periodic%init(gpC, gpE, 1, one, one, dz, sp_gpC, sp_gpE,.true.,.true.,isPeriodic)
      case default
         call GracefulExit("Invalid choice of numerical scheme in vertical",434)
      end select
   else
      select case(this%scheme)
      case(cd06)
         allocate(this%derES, this%derSE, this%derOS, this%derSO, this%derSS, this%derEE, this%derOO, this%derEO, this%derOE)
         call this%derES%init(this%gp%zsz(3), dz, isTopEven = .true., isBotEven = .true., &
                                   isTopSided = .false., isBotSided = .true.)
         call this%derSE%init(this%gp%zsz(3), dz, isTopEven = .true., isBotEven = .true., &
                                   isTopSided = .true., isBotSided = .false.)
         call this%derSS%init(this%gp%zsz(3), dz, isTopEven = .true., isBotEven = .true., &
                                   isTopSided = .true., isBotSided = .true.)
         call this%derSO%init(this%gp%zsz(3), dz, isTopEven = .false., isBotEven = .false., &
                                   isTopSided = .true., isBotSided = .false.)
         call this%derOS%init(this%gp%zsz(3), dz, isTopEven = .false., isBotEven = .false., &
                                   isTopSided = .false., isBotSided = .true.)
         call this%derOO%init(this%gp%zsz(3), dz, isTopEven = .false., isBotEven = .false., &
                                   isTopSided = .false., isBotSided = .false.)
         call this%derEE%init(this%gp%zsz(3), dz, isTopEven = .true., isBotEven = .true., &
                                   isTopSided = .false., isBotSided = .false.)
         call this%derOE%init(this%gp%zsz(3), dz, isTopEven = .false., isBotEven = .true., &
                                   isTopSided = .false., isBotSided = .false.)
         call this%derEO%init(this%gp%zsz(3), dz, isTopEven = .true., isBotEven = .false., &
                                   isTopSided = .false., isBotSided = .false.)
         call message(0,"6th order stagerred Compact Finite difference schemes&
                         & initialized in the Z direction")
      case(fd02)
         call this%fd02_ss%init(gpC, gpE, 1, one, one, dz, sp_gpC, sp_gpE,.true.,.true.)
         call this%fd02_sn%init(gpC, gpE, 1, one, one, dz, sp_gpC, sp_gpE,.true.,.false.)
         call this%fd02_ns%init(gpC, gpE, 1, one, one, dz, sp_gpC, sp_gpE,.false.,.true.)
         call this%fd02_nn%init(gpC, gpE, 1, one, one, dz, sp_gpC, sp_gpE,.false.,.false.)
         call message(0,"2nd order stagerred Explicit Finite difference schemes&
                         & initialized in the Z direction")
      case default
         call gracefulExit("Invalid choice for numerical scheme in vertical direction", 323)
      end select
   end if 
end subroutine


subroutine destroy(this)
   class(Pade6stagg), intent(inout) :: this
  
   if (this%isPeriodic) then
      if (this%scheme == cd06 .OR. this%scheme == fd02) then
         Call this%derPeriodic%destroy()
         Deallocate(this%derPeriodic)
      else if (this%scheme == fourierColl) then
         nullify(this%spectC)
      end if
   else
      if (this%scheme == cd06) then
         deallocate(this%derES, this%derSE, this%derOS, this%derSO, this%derSS, this%derEE, this%derOO, this%derEO, this%derOE)
      end if
   end if 
   nullify(this%gp, this%sp_gp)

end subroutine

subroutine d2dz2_C2C_real(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)  ), intent(in)  :: input
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)  ), intent(out) :: output
   integer, intent(in) :: bot, top

   if (this%isPeriodic) then
      select case(this%scheme) 
      case (fourierColl) 
         call this%spectC%d2dz2_C2C_spect(input, output)
      case (cd06)
         call this%derPeriodic%d2dz2_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
      case (fd02)
         call this%fd02_periodic%d2dz2_C2C(input, output,.false.,.false.)
      end select
   else
      select case (this%scheme) 
      case(fd02)
         select case (bot) 
            case(-1)
               if (top == -1) then
                  call this%fd02_nn%d2dz2_C2C(input,output,.false.,.false.)
               else if (top == 1) then
                  call this%fd02_nn%d2dz2_C2C(input,output,.true.,.false.)
               else
                  output = 0.d0
               end if 
            case(1)
               if (top == -1) then
                  call this%fd02_nn%d2dz2_C2C(input,output,.false.,.true.)
               else if (top == 1) then
                  call this%fd02_nn%d2dz2_C2C(input,output,.true.,.true.)
               else
                  output = 0.d0
               end if 
            case default 
               output = 0.d0
         end select
      case(cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%d2dz2_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derEO%d2dz2_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%d2dz2_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derEE%d2dz2_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select
   end if 

end subroutine 


subroutine d2dz2_C2C_cmplx(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)  ), intent(in)  :: input
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)  ), intent(out) :: output
   integer, intent(in) :: bot, top

   if (this%isPeriodic) then
      select case(this%scheme) 
      case (fourierColl) 
         call this%spectC%d2dz2_C2C_spect(input, output)
      case (cd06)
         call this%derPeriodic%d2dz2_C2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
      case (fd02)
         call this%fd02_periodic%d2dz2_C2C(input, output,.false.,.false.)
      end select
   else
      select case (this%scheme) 
      case(fd02)
         select case (bot) 
            case(-1)
               if (top == -1) then
                  call this%fd02_nn%d2dz2_C2C(input,output,.false.,.false.)
               else if (top == 1) then
                  call this%fd02_nn%d2dz2_C2C(input,output,.true.,.false.)
               else
                  output = 0.d0
               end if 
            case(1)
               if (top == -1) then
                  call this%fd02_nn%d2dz2_C2C(input,output,.false.,.true.)
               else if (top == 1) then
                  call this%fd02_nn%d2dz2_C2C(input,output,.true.,.true.)
               else
                  output = 0.d0
               end if 
            case default 
               output = 0.d0
         end select
      case(cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%d2dz2_C2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derEO%d2dz2_C2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%d2dz2_C2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derEE%d2dz2_C2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select
   end if 

end subroutine 

subroutine d2dz2_E2E_real(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3) + 1), intent(in)  :: input
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3) + 1), intent(out) :: output
   integer, intent(in) :: bot, top

   if (this%isPeriodic) then
      select case(this%scheme) 
      case (fourierColl) 
         call this%spectC%d2dz2_E2E_spect(input, output)
      case (cd06)
         call this%derPeriodic%d2dz2_E2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
      case (fd02)
         call this%fd02_periodic%d2dz2_E2E(input, output,.false.,.false.)
      end select 
   else
      select case (this%scheme)
      case(fd02)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%fd02_nn%d2dz2_E2E(input, output, .false.,.false.)
            elseif (top ==  1) then
               call this%fd02_nn%d2dz2_E2E(input, output, .true.,.false.)
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%fd02_nn%d2dz2_E2E(input, output, .false.,.true.)
            elseif (top ==  1) then
               call this%fd02_nn%d2dz2_E2E(input, output, .true.,.true.)
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      case(cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%d2dz2_E2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derEO%d2dz2_E2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%d2dz2_E2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derEE%d2dz2_E2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select
   end if  

end subroutine 

subroutine d2dz2_E2E_cmplx(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)+1), intent(in)  :: input
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)+1), intent(out) :: output
   integer, intent(in) :: bot, top

   if (this%isPeriodic) then
      select case(this%scheme) 
      case (fourierColl) 
         call this%spectC%d2dz2_E2E_spect(input, output)
      case (cd06)
         call this%derPeriodic%d2dz2_E2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
      case (fd02)
         call this%fd02_periodic%d2dz2_E2E(input, output,.false.,.false.)
      end select 
   else
      select case (this%scheme)
      case(fd02)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%fd02_nn%d2dz2_E2E(input, output, .false.,.false.)
            elseif (top ==  1) then
               call this%fd02_nn%d2dz2_E2E(input, output, .true.,.false.)
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%fd02_nn%d2dz2_E2E(input, output, .false.,.true.)
            elseif (top ==  1) then
               call this%fd02_nn%d2dz2_E2E(input, output, .true.,.true.)
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      case(cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%d2dz2_E2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derEO%d2dz2_E2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%d2dz2_E2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derEE%d2dz2_E2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select
   end if 

end subroutine 


subroutine ddz_C2C(this, input, output, bot, top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)), intent(in)  :: input
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)), intent(out) :: output
   integer, intent(in) :: bot, top

   if (this%isPeriodic) then
      select case(this%scheme) 
      case (fourierColl)
         call this%spectC%ddz_C2C_spect(input, output)
      case (cd06)
         call this%derPeriodic%ddz_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
      case (fd02)    
         call this%fd02_periodic%ddz_C2C(input, output,.false.,.false.)
      end select
   else 
      select case (this%scheme) 
      case(fd02)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%fd02_nn%ddz_C2C(input,output,.false.,.false.)
            elseif (top ==  0) then
               call this%fd02_sn%ddz_C2C(input,output,.true.,.false.)
            elseif (top ==  1) then
               call this%fd02_nn%ddz_C2C(input,output,.true.,.false.)
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%fd02_ns%ddz_C2C(input,output,.false.,.true.)
            elseif (top ==  0) then
               call this%fd02_ss%ddz_C2C(input,output,.false.,.true.)
            elseif (top ==  1) then
               call this%fd02_ns%ddz_C2C(input,output,.true.,.true.)
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%fd02_nn%ddz_C2C(input,output,.false.,.true.)
            elseif (top ==  0) then
               call this%fd02_sn%ddz_C2C(input,output,.false.,.true.)
            elseif (top ==  1) then
               call this%fd02_nn%ddz_C2C(input,output,.true.,.true.)
            end if
         end select
      case(cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%ddz_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSO%ddz_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derEO%ddz_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%derOS%ddz_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSS%ddz_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derES%ddz_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%ddz_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSE%ddz_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derEE%ddz_C2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select 
   end if 


end subroutine 





subroutine ddz_C2E_real(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)  ), intent(in)  :: input
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)+1), intent(out) :: output
   integer, intent(in) :: bot, top

   if (this%isPeriodic) then
      select case(this%scheme) 
      case (fourierColl)
         call this%spectC%ddz_C2E_spect(input, output)
      case (cd06)
         call this%derPeriodic%ddz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
      case (fd02)    
         call this%fd02_periodic%ddz_C2E(input, output,.false.,.false.)
      end select
   else
      select case (this%scheme) 
      case(fd02)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%fd02_nn%ddz_C2E(input,output,.false.,.false.)
            elseif (top ==  0) then
               call this%fd02_sn%ddz_C2E(input,output,.true.,.false.)
            elseif (top ==  1) then
               call this%fd02_nn%ddz_C2E(input,output,.true.,.false.)
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%fd02_ns%ddz_C2E(input,output,.false.,.true.)
            elseif (top ==  0) then
               call this%fd02_ss%ddz_C2E(input,output,.false.,.true.)
            elseif (top ==  1) then
               call this%fd02_ns%ddz_C2E(input,output,.true.,.true.)
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%fd02_nn%ddz_C2E(input,output,.false.,.true.)
            elseif (top ==  0) then
               call this%fd02_sn%ddz_C2E(input,output,.false.,.true.)
            elseif (top ==  1) then
               call this%fd02_nn%ddz_C2E(input,output,.true.,.true.)
            end if
         end select
      case(cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%ddz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSO%ddz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derEO%ddz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%derOS%ddz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSS%ddz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derES%ddz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%ddz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSE%ddz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derEE%ddz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select
   end if 

end subroutine 

subroutine ddz_C2E_cmplx(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)  ), intent(in)  :: input
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)+1), intent(out) :: output
   integer, intent(in) :: bot, top

   if (this%isPeriodic) then
      select case(this%scheme) 
      case (fourierColl)
         call this%spectC%ddz_C2E_spect(input, output)
      case (cd06)
         call this%derPeriodic%ddz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
      case (fd02)    
         call this%fd02_periodic%ddz_C2E(input, output,.false.,.false.)
      end select 
   else
      select case (this%scheme) 
      case(fd02)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%fd02_nn%ddz_C2E(input,output,.false.,.false.)
            elseif (top ==  0) then
               call this%fd02_sn%ddz_C2E(input,output,.true.,.false.)
            elseif (top ==  1) then
               call this%fd02_nn%ddz_C2E(input,output,.true.,.false.)
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%fd02_ns%ddz_C2E(input,output,.false.,.true.)
            elseif (top ==  0) then
               call this%fd02_ss%ddz_C2E(input,output,.false.,.true.)
            elseif (top ==  1) then
               call this%fd02_ns%ddz_C2E(input,output,.true.,.true.)
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%fd02_nn%ddz_C2E(input,output,.false.,.true.)
            elseif (top ==  0) then
               call this%fd02_sn%ddz_C2E(input,output,.false.,.true.)
            elseif (top ==  1) then
               call this%fd02_nn%ddz_C2E(input,output,.true.,.true.)
            end if
         end select
      case(cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%ddz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  0) then
               call this%derSO%ddz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derEO%ddz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%derOS%ddz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  0) then
               call this%derSS%ddz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derES%ddz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%ddz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  0) then
               call this%derSE%ddz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derEE%ddz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select
   end if 

end subroutine 

subroutine ddz_E2C_real(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)+1), intent(in)  :: input
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)  ), intent(out) :: output
   integer, intent(in) :: bot, top

   if (this%isPeriodic) then
      select case(this%scheme) 
      case (fourierColl) 
         call this%spectC%ddz_E2C_spect(input, output)
      case (cd06)
         call this%derPeriodic%ddz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
      case (fd02)    
         call this%fd02_periodic%ddz_E2C(input, output)
      end select 
   else
      select case(this%scheme)
      case(fd02)
         call this%fd02_nn%ddz_E2C(input,output)
      case(cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%ddz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSO%ddz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derEO%ddz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%derOS%ddz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSS%ddz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derES%ddz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%ddz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSE%ddz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derEE%ddz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select
   end if 

end subroutine 

subroutine ddz_E2C_cmplx(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)+1), intent(in)  :: input
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)  ), intent(out) :: output
   integer, intent(in) :: bot, top

   if (this%isPeriodic) then
      select case(this%scheme) 
      case (fourierColl) 
         call this%spectC%ddz_E2C_spect(input, output)
      case (cd06)
         call this%derPeriodic%ddz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
      case (fd02)    
         call this%fd02_periodic%ddz_E2C(input, output)
      end select 
   else
      select case(this%scheme)
      case(fd02)
         call this%fd02_nn%ddz_E2C(input,output)
      case(cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%ddz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  0) then
               call this%derSO%ddz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derEO%ddz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%derOS%ddz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  0) then
               call this%derSS%ddz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derES%ddz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%ddz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  0) then
               call this%derSE%ddz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derEE%ddz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select
   end if 
end subroutine 

subroutine interpz_C2E_real(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)  ), intent(in)  :: input
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)+1), intent(out) :: output
   integer, intent(in) :: bot, top

   if (this%isPeriodic) then
      select case(this%scheme) 
      case (fourierColl) 
         call this%spectC%interp_C2E_spect(input, output)
      case (cd06)
         call this%derPeriodic%interpz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
      case (fd02)
         call this%fd02_periodic%InterpZ_Cell2Edge(input, output,0.d0,0.d0)
      end select 
   else
      select case(this%scheme)
      case (fd02)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zero,zero)
            elseif (top ==  0) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zero,zero)
               output(:,:,this%gp%zsz(3)+1) = two*input(:,:,this%gp%zsz(3)) - output(:,:,this%gp%zsz(3))
            elseif (top ==  1) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zero,zero)
               output(:,:,this%gp%zsz(3)+1) = input(:,:,this%gp%zsz(3)) 
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zero,zero)
               output(:,:,1) = two*input(:,:,1) - output(:,:,2)
            elseif (top ==  0) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zero,zero)
               output(:,:,1) = two*input(:,:,1) - output(:,:,2)
               output(:,:,this%gp%zsz(3)+1) = two*input(:,:,this%gp%zsz(3)) - output(:,:,this%gp%zsz(3))
            elseif (top ==  1) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zero,zero)
               output(:,:,1) = two*input(:,:,1) - output(:,:,2)
               output(:,:,this%gp%zsz(3)+1) = input(:,:,this%gp%zsz(3)) 
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zero,zero)
               output(:,:,1) = input(:,:,1) 
            elseif (top ==  0) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zero,zero)
               output(:,:,1) = input(:,:,1) 
               output(:,:,this%gp%zsz(3)+1) = two*input(:,:,this%gp%zsz(3)) - output(:,:,this%gp%zsz(3))
            elseif (top ==  1) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zero,zero)
               output(:,:,1) = input(:,:,1) 
               output(:,:,this%gp%zsz(3)+1) = input(:,:,this%gp%zsz(3)) 
            end if
         end select
      case (cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%interpz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSO%interpz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derEO%interpz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%derOS%interpz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSS%interpz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derES%interpz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%interpz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSE%interpz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derEE%interpz_C2E(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select
   end if 
end subroutine
 
subroutine interpz_C2E_cmplx(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)  ), intent(in)  :: input
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)+1), intent(out) :: output
   integer, intent(in) :: bot, top

   if (this%isPeriodic) then
      select case(this%scheme) 
      case (fourierColl) 
         call this%spectC%interp_C2E_spect(input, output)
      case (cd06)
         call this%derPeriodic%interpz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
      case (fd02)
         call this%fd02_periodic%InterpZ_Cell2Edge(input, output, (0.d0,0.d0), (0.d0,0.d0))
      end select 
   else
      select case(this%scheme)
      case (fd02)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zeroC,zeroC)
            elseif (top ==  0) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zeroC,zeroC)
               output(:,:,this%gp%zsz(3)+1) = two*input(:,:,this%gp%zsz(3)) - output(:,:,this%gp%zsz(3))
            elseif (top ==  1) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zeroC,zeroC)
               output(:,:,this%gp%zsz(3)+1) = input(:,:,this%gp%zsz(3)) 
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zeroC,zeroC)
               output(:,:,1) = two*input(:,:,1) - output(:,:,2)
            elseif (top ==  0) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zeroC,zeroC)
               output(:,:,1) = two*input(:,:,1) - output(:,:,2)
               output(:,:,this%gp%zsz(3)+1) = two*input(:,:,this%gp%zsz(3)) - output(:,:,this%gp%zsz(3))
            elseif (top ==  1) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zeroC,zeroC)
               output(:,:,1) = two*input(:,:,1) - output(:,:,2)
               output(:,:,this%gp%zsz(3)+1) = input(:,:,this%gp%zsz(3)) 
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zeroC,zeroC)
               output(:,:,1) = input(:,:,1) 
            elseif (top ==  0) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zeroC,zeroC)
               output(:,:,1) = input(:,:,1) 
               output(:,:,this%gp%zsz(3)+1) = two*input(:,:,this%gp%zsz(3)) - output(:,:,this%gp%zsz(3))
            elseif (top ==  1) then
               call this%fd02_nn%interpZ_Cell2Edge(input,output,zeroC,zeroC)
               output(:,:,1) = input(:,:,1) 
               output(:,:,this%gp%zsz(3)+1) = input(:,:,this%gp%zsz(3)) 
            end if
         end select
      case (cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%interpz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  0) then
               call this%derSO%interpz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derEO%interpz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%derOS%interpz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  0) then
               call this%derSS%interpz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derES%interpz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%interpz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  0) then
               call this%derSE%interpz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derEE%interpz_C2E(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select
   end if 
end subroutine 

subroutine interpz_E2C_real(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)+1), intent(in)  :: input
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)  ), intent(out) :: output
   integer, intent(in) :: bot, top


   if (this%isPeriodic) then
      select case(this%scheme) 
      case (fourierColl) 
         call this%spectC%interp_E2C_spect(input, output)
      case (cd06)
         call this%derPeriodic%interpz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
      case (fd02)
         call this%fd02_periodic%InterpZ_Edge2Cell(input, output)
      end select 
   else
      select case(this%scheme)
      case(fd02)
         call this%fd02_nn%interpZ_Edge2Cell(input,output)
      case(cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%interpz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSO%interpz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derEO%interpz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%derOS%interpz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSS%interpz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derES%interpz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%interpz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  0) then
               call this%derSE%interpz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            elseif (top ==  1) then
               call this%derEE%interpz_E2C(input,output,this%gp%zsz(1),this%gp%zsz(2))
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select 
   end if 

end subroutine

subroutine interpz_E2C_cmplx(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)+1), intent(in)  :: input
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)  ), intent(out) :: output
   integer, intent(in) :: bot, top

   if (this%isPeriodic) then
      select case(this%scheme) 
      case (fourierColl) 
         call this%spectC%interp_E2C_spect(input, output)
      case (cd06)
         call this%derPeriodic%interpz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
      case (fd02)
         call this%fd02_periodic%InterpZ_Edge2Cell(input, output)
      end select 
   else
      select case(this%scheme)
      case(fd02)
         call this%fd02_nn%interpZ_Edge2Cell(input,output)
      case(cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%interpz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  0) then
               call this%derSO%interpz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derEO%interpz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%derOS%interpz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  0) then
               call this%derSS%interpz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derES%interpz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%interpz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  0) then
               call this%derSE%interpz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            elseif (top ==  1) then
               call this%derEE%interpz_E2C(input,output,this%sp_gp%zsz(1),this%sp_gp%zsz(2))
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select
   end if 
end subroutine

subroutine getModifiedWavenumbers(this, k, kp)
    class(Pade6stagg), intent(in) :: this
    real(rkind), dimension(:), intent(in)  :: k
    real(rkind), dimension(:), intent(out) :: kp


    select case(this%scheme)
    case (fourierColl) 
       kp = k
    case (cd06)
       call getmodCD06stagg(k,this%dz,kp)
    case (fd02)
       call getmodFD02stagg(k,this%dz,kp)
    case default
       call GracefulExit("You have not provided the modified wavenumber function &
                        & for the chosen stagerred scheme",23)
    end select


end subroutine

pure function getApproxPoincareConstant(this) result(const)
   use constants, only: pi
   class(Pade6stagg), intent(in) :: this
   real(rkind) :: const   

   select case(this%scheme)
   case (fourierColl)
      const = 1.d0/sqrt(pi**2)
   case (cd06) ! For wall bounded flows the code gives second order near wall
      const = 1.d0/sqrt(3.d0)
   case (fd02)
      const = 1.d0/sqrt(3.d0)
   end select 

end function


pure subroutine getmodCD06stagg(k,dx,kp)
    real(rkind), dimension(:), intent(in)  :: k
    real(rkind), intent(in) :: dx
    real(rkind), dimension(:), intent(out) :: kp
    real(rkind), dimension(:), allocatable :: omega
    real(rkind), parameter :: alpha = 9._rkind/62._rkind
    real(rkind), parameter :: beta = 0._rkind
    real(rkind), parameter :: a = 63._rkind/62._rkind
    real(rkind), parameter :: b = 17._rkind/62._rkind
    real(rkind), parameter :: c = 0._rkind

    allocate(omega(size(k)))
    omega = k*dx
    kp = (two*a*sin(omega/two) + (two/three)*b*sin(three*omega/two) + &
         (two/five)*c*sin(five*omega/two))/(one + two*alpha*cos(omega) +&
          two*beta*cos(two*omega))
    kp = kp/dx
    deallocate(omega)

end subroutine

pure subroutine getmodFD02stagg(k,dx,kp)
    real(rkind), dimension(:), intent(in)  :: k
    real(rkind), intent(in) :: dx
    real(rkind), dimension(:), intent(out) :: kp
    real(rkind), dimension(:), allocatable :: omega
    real(rkind), parameter :: alpha = 0._rkind
    real(rkind), parameter :: beta = 0._rkind
    real(rkind), parameter :: a = 1._rkind
    real(rkind), parameter :: b = 0._rkind 
    real(rkind), parameter :: c = 0._rkind

    allocate(omega(size(k)))
    omega = k*dx
    kp = (two*a*sin(omega/two) + (two/three)*b*sin(three*omega/two) + &
         (two/five)*c*sin(five*omega/two))/(one + two*alpha*cos(omega) +&
          two*beta*cos(two*omega))
    kp = kp/dx
    deallocate(omega)

end subroutine

subroutine ddz_1d_C2C(this, input, output, bot, top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(1,1,this%gp%zsz(3)), intent(in)  :: input
   real(rkind), dimension(1,1,this%gp%zsz(3)), intent(out) :: output
   integer, intent(in) :: bot, top

   if (this%isPeriodic) then
        output = 0.d0 
   else 
      select case (this%scheme) 
      case(fd02)
        call message(0,"WARNING: Second order FD cannot be used for 1D derivative evaluations.")
        output = 0.d0
      case(cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%ddz_C2C(input,output,1,1)
            elseif (top ==  0) then
               call this%derSO%ddz_C2C(input,output,1,1)
            elseif (top ==  1) then
               call this%derEO%ddz_C2C(input,output,1,1)
            else 
               output = 0.d0
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%derOS%ddz_C2C(input,output,1,1)
            elseif (top ==  0) then
               call this%derSS%ddz_C2C(input,output,1,1)
            elseif (top ==  1) then
               call this%derES%ddz_C2C(input,output,1,1)
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%ddz_C2C(input,output,1,1)
            elseif (top ==  0) then
               call this%derSE%ddz_C2C(input,output,1,1)
            elseif (top ==  1) then
               call this%derEE%ddz_C2C(input,output,1,1)
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select 
   end if 

end subroutine 

subroutine ddz_1d_C2E(this, input, output, bot, top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(1,1,this%gp%zsz(3)), intent(in)  :: input
   real(rkind), dimension(1,1,this%gp%zsz(3)+1), intent(out) :: output
   integer, intent(in) :: bot, top

   if (this%isPeriodic) then
        output = 0.d0 
   else 
      select case (this%scheme) 
      case(fd02)
        call message(0,"WARNING: Second order FD cannot be used for 1D derivative evaluations.")
        output = 0.d0
      case(cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%ddz_C2E(input,output,1,1)
            elseif (top ==  0) then
               call this%derSO%ddz_C2E(input,output,1,1)
            elseif (top ==  1) then
               call this%derEO%ddz_C2E(input,output,1,1)
            else 
               output = 0.d0
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%derOS%ddz_C2E(input,output,1,1)
            elseif (top ==  0) then
               call this%derSS%ddz_C2E(input,output,1,1)
            elseif (top ==  1) then
               call this%derES%ddz_C2E(input,output,1,1)
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%ddz_C2E(input,output,1,1)
            elseif (top ==  0) then
               call this%derSE%ddz_C2E(input,output,1,1)
            elseif (top ==  1) then
               call this%derEE%ddz_C2E(input,output,1,1)
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select 
   end if 

end subroutine 

subroutine interp_1d_E2C(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(1,1,this%gp%zsz(3)+1), intent(in)  :: input
   real(rkind), dimension(1,1,this%gp%zsz(3)  ), intent(out) :: output
   integer, intent(in) :: bot, top


   if (this%isPeriodic) then
        output = 0.d0 
   else
      select case(this%scheme)
      case(fd02)
            call message(0,"WARNING: Second order FD cannot be used for 1D derivative evaluations.")
      case(cd06)
         select case (bot)
         case(-1) 
            if     (top == -1) then
               call this%derOO%interpz_E2C(input,output,1,1)
            elseif (top ==  0) then
               call this%derSO%interpz_E2C(input,output,1,1)
            elseif (top ==  1) then
               call this%derEO%interpz_E2C(input,output,1,1)
            else 
               output = 0.d0
            end if
         case(0)  ! bottom = sided
            if     (top == -1) then
               call this%derOS%interpz_E2C(input,output,1,1)
            elseif (top ==  0) then
               call this%derSS%interpz_E2C(input,output,1,1)
            elseif (top ==  1) then
               call this%derES%interpz_E2C(input,output,1,1)
            else 
               output = 0.d0
            end if
         case(1)  ! bottom = even
            if     (top == -1) then
               call this%derOE%interpz_E2C(input,output,1,1)
            elseif (top ==  0) then
               call this%derSE%interpz_E2C(input,output,1,1)
            elseif (top ==  1) then
               call this%derEE%interpz_E2C(input,output,1,1)
            else 
               output = 0.d0
            end if
         case default
            output = 0.d0
         end select
      end select 
   end if 

end subroutine
end module 
