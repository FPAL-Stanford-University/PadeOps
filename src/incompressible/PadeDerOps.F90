module PadeDerOps
   use kind_parameters, only: rkind
   use cd06staggstuff, only: cd06stagg
   use decomp_2d

   implicit none
   
   private

   public :: Pade6stagg

   type Pade6stagg
      type(cd06stagg), allocatable :: derOO, derEE, derOE, derEO, derOS, derSO, derSE, derES, derSS
      type(decomp_info), pointer :: gp, sp_gp
      contains
      procedure          :: init
      procedure          :: destroy
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
      generic            :: d2dz2_C2C => d2dz2_C2C_real, d2dz2_C2C_cmplx
      generic            :: d2dz2_E2E => d2dz2_E2E_real, d2dz2_E2E_cmplx
      generic            :: interpz_C2E => interpz_C2E_real, interpz_C2E_cmplx
      generic            :: interpz_E2C => interpz_E2C_real, interpz_E2C_cmplx
   end type
contains

subroutine init(this, gp, sp_gp, dz)
   class(Pade6stagg), intent(out) :: this
   type(decomp_info), intent(in), target :: gp, sp_gp
   real(rkind), intent(in) :: dz

   this%gp => gp
   this%sp_gp => sp_gp
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

end subroutine


subroutine destroy(this)
   class(Pade6stagg), intent(inout) :: this
   
   deallocate(this%derES, this%derSE, this%derOS, this%derSO, this%derSS, this%derEE, this%derOO, this%derEO, this%derOE)
   nullify(this%gp, this%sp_gp)

end subroutine

subroutine d2dz2_C2C_real(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)  ), intent(in)  :: input
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)  ), intent(out) :: output
   integer, intent(in) :: bot, top

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

end subroutine 


subroutine d2dz2_C2C_cmplx(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)  ), intent(in)  :: input
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)  ), intent(out) :: output
   integer, intent(in) :: bot, top

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

end subroutine 


subroutine d2dz2_E2E_real(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3) + 1), intent(in)  :: input
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3) + 1), intent(out) :: output
   integer, intent(in) :: bot, top

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

end subroutine 

subroutine d2dz2_E2E_cmplx(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)+1), intent(in)  :: input
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)+1), intent(out) :: output
   integer, intent(in) :: bot, top

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

end subroutine 


subroutine ddz_C2E_real(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)  ), intent(in)  :: input
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)+1), intent(out) :: output
   integer, intent(in) :: bot, top

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

end subroutine 

subroutine ddz_C2E_cmplx(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)  ), intent(in)  :: input
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)+1), intent(out) :: output
   integer, intent(in) :: bot, top

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

end subroutine 

subroutine ddz_E2C_real(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)+1), intent(in)  :: input
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)  ), intent(out) :: output
   integer, intent(in) :: bot, top

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

end subroutine 

subroutine ddz_E2C_cmplx(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)+1), intent(in)  :: input
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)  ), intent(out) :: output
   integer, intent(in) :: bot, top

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

end subroutine 

subroutine interpz_C2E_real(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)  ), intent(in)  :: input
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)+1), intent(out) :: output
   integer, intent(in) :: bot, top

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

end subroutine
 
subroutine interpz_C2E_cmplx(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)  ), intent(in)  :: input
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)+1), intent(out) :: output
   integer, intent(in) :: bot, top

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

end subroutine 

subroutine interpz_E2C_real(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)+1), intent(in)  :: input
   real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)  ), intent(out) :: output
   integer, intent(in) :: bot, top

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

end subroutine

subroutine interpz_E2C_cmplx(this,input,output,bot,top)
   class(Pade6stagg), intent(in) :: this
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)+1), intent(in)  :: input
   complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)  ), intent(out) :: output
   integer, intent(in) :: bot, top

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

end subroutine
end module 
