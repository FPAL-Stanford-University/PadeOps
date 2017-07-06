program test_PressureProj_PeriodicCD06
   use kind_parameters, only: rkind
   use mpi
   use decomp_2d
   use decomp_2d_io
   use exits, only: message 
   use spectralMod, only: spectral
   use constants, only: pi, two
   use PadePoissonMod, only: padepoisson
   use PadeDerOps, only: Pade6stagg

   implicit none
   real   (rkind), dimension(:,:,:), allocatable :: u, v, w, uTrue, vTrue, wTrue
   complex(rkind), dimension(:,:,:), allocatable :: uhat, vhat, what
   logical, parameter :: isPeriodic = .true. 
   integer, parameter :: scheme = 1

   integer, parameter :: nx = 128, ny = 128, nz = 128
   real(rkind), parameter :: dx = two*pi/real(nx,rkind), dy = two*pi/real(ny,rkind), dz = two*pi/real(nz,rkind) 
   type(decomp_info) :: gp, gpE
   integer :: ierr 
   type(spectral), allocatable :: spectC, spectE
   type(padepoisson), allocatable :: poiss
   type(Pade6stagg), allocatable :: der


   call mpi_init(ierr)
   call decomp_2d_init(nx,ny,nz,0,0)
   call get_decomp_info(gp)
   call decomp_info_init(nx,ny,nz+1, gpE)

   allocate(u(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(v(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(w(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))
   allocate(uTrue(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(vTrue(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(wTrue(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))


   call decomp_2d_read_one(1,uTrue,'/home/aditya90/Codes/PadeOps/data/HIT_u_PostProj.dat',gp )
   call decomp_2d_read_one(1,vTrue,'/home/aditya90/Codes/PadeOps/data/HIT_v_PostProj.dat',gp )
   call decomp_2d_read_one(1,wTrue,'/home/aditya90/Codes/PadeOps/data/HIT_w_PostProj.dat',gpE)
   
   call decomp_2d_read_one(1,u,'/home/aditya90/Codes/PadeOps/data/HIT_u_PreProj.dat',gp )
   call decomp_2d_read_one(1,v,'/home/aditya90/Codes/PadeOps/data/HIT_v_PreProj.dat',gp )
   call decomp_2d_read_one(1,w,'/home/aditya90/Codes/PadeOps/data/HIT_w_PreProj.dat',gpE)
   
   call message(0, "Completed data read")

   allocate(spectC, spectE)
   call spectC%init("x", nx, ny, nz  , dx, dy, dz, "four", "2/3rd", 2,.false.)
   call spectE%init("x", nx, ny, nz+1, dx, dy, dz, "four", "2/3rd", 2,.false.)

   allocate(der)
   call der%init( gp, spectC%spectdecomp, gpE, spectE%spectdecomp, dz, scheme, isPeriodic)
   allocate(poiss)
   call poiss%init(dx, dy, dz, spectC, spectE, .false., two*pi, .true., gp, der, isPeriodic)

   call spectC%alloc_r2c_out(uhat)
   call spectC%alloc_r2c_out(vhat)
   call spectE%alloc_r2c_out(what)

   call spectC%fft(u, uhat)
   call spectC%fft(v, vhat)
   call spectE%fft(w, what)


   print*, shape(uhat)
   print*, shape(u)

   call poiss%PressureProjection(uhat, vhat, what)

   print*, u(12,12,12)
   print*, v(12,12,12)
   print*, w(12,12,12)


   call spectC%ifft(uhat, u)
   call spectC%ifft(vhat, v)
   call spectE%ifft(what, w)

   print*, u(12,12,12), uTrue(12,12,12)
   print*, v(12,12,12), vTrue(12,12,12)
   print*, w(12,12,12), wTrue(12,12,12)



   if ((maxval(abs(u - uTrue)) > 1.d-12) .or. (maxval(abs(v - vTrue)) > 1.d-12) .or. (maxval(abs(w - wTrue)) > 1.d-12)) then
      call message(0,"TEST FAILED")
   else 
      call message(0,"TEST PASSED")
   end if

   deallocate(der)
   deallocate(poiss)
   deallocate(uhat, vhat, what)
   deallocate(u, v, w, uTrue, vTrue, wTrue)
   
   call MPI_Finalize(ierr)

end program 
