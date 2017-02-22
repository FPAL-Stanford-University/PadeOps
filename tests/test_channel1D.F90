program test_channel1D
   use kind_parameters, only: rkind
   use PadeDerOps, only: Pade6stagg
   use constants, only: pi
   use cd06staggstuff, only: cd06stagg
   use timer, only: tic, toc

   implicit none

   integer, parameter :: nz = 32
   real(rkind), dimension(1,1,nz) :: u, uexact,z, urhs, u1, u2
   real(rkind), parameter :: Re = 100.d0, F = 0.02d0, tstop = 1500.d0, CFL = 0.2d0
   real(rkind) :: dt
   integer :: k
   real(rkind) :: dz, time
   type(cd06stagg) :: der1

   dz = 2.d0/real(nz,rkind)
   do k = 1,nz
      z(1,1,k) = real(k,rkind)*dz
   end do 
   z = z - dz/2.d0

   u = sin(pi*z/2.d0)
   uexact = z*(2.d0 - z)
   !u = uexact

   dt = CFL*Re*(dz**2)

   print*, "DT:", dt
   call der1%init(nz, dz, isTopEven = .false., isBotEven = .false., &
                       isTopSided = .false., isBotSided = .false.)
   time = 0.d0

   call tic()
   do while(time < tstop)
      call get_rhs(u,urhs)
      u1 = u + dt*urhs
      call get_rhs(u1,urhs)
      u2 = (3.d0/4.d0)*u + (1.d0/4.d0)*u1 + (1.d0/4.d0)*dt*urhs
      call get_rhs(u2,urhs)
      u = (1.d0/3.d0)*u + (2.d0/3.d0)*u2 + (2.d0/3.d0)*dt*urhs
      time = time + dt
   end do
   call toc()
   print*, "Error:", maxval(abs(u - uexact))
   call der1%destroy()

contains
   subroutine get_rhs(uin,urhsin)
      real(rkind), dimension(1,1,nz), intent(in) :: uin
      real(rkind), dimension(1,1,nz), intent(out) :: urhsin
      
      call der1%d2dz2_C2C(uin,urhs,size(u,1),size(u,2))
      urhsin = (1.d0/Re)*urhsin + F
   

   end subroutine


end program
