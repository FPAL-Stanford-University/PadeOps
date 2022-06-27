program test_HITforcing
   use kind_parameters, only: rkind, clen
   use mpi
   use decomp_2d
   use decomp_2d_io
   use exits,           only: message 
   use spectralMod,     only: spectral
   use constants,       only: pi, im0, two
   use PadePoissonMod,  only: padepoisson
   use PadeDerOps,      only: Pade6stagg
   use forcingMod,      only: HIT_shell_forcing 
   use reductions,      only: p_maxval
   use fortran_assert,  only: assert

   implicit none
   real   (rkind), dimension(:,:,:), allocatable :: Pressure, u, v, w, wC 
   complex(rkind), dimension(:,:,:), allocatable :: uhat, vhat, what, wChat, wcHatz, whatz, urhs, vrhs, wrhs
   real(rkind),    dimension(:,:,:), allocatable :: divergence
   logical, parameter :: isPeriodic = .true. 
   integer, parameter :: scheme = 2
   complex(rkind), dimension(:,:,:,:), allocatable :: cbuffzC
   complex(rkind), dimension(:,:,:), allocatable :: cbuffzE1, cbuffyE, cbuffyC

   integer, parameter :: nx = 640, ny = 64, nz = 64
   real(rkind), parameter :: dx = 20.d0*pi/real(nx,rkind), &
                             dy = two*pi/real(ny,rkind), &
                             dz = two*pi/real(nz,rkind) 
   type(decomp_info) :: gp, gpE
   integer :: ierr 
   type(spectral), allocatable :: spectC, spectE
   type(padepoisson), allocatable :: poiss
   type(Pade6stagg), allocatable :: der
   type(HIT_shell_forcing) :: HITForce
   integer :: tidStart = 0 ! Time index to initialize seeding in HITForce
   character(len=clen) :: datadir
   
   call mpi_init(ierr)
   call decomp_2d_init(nx,ny,nz,0,0)
   call get_decomp_info(gp)
   call decomp_info_init(nx,ny,nz+1, gpE)

   allocate(u (gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(v (gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(wC(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(w (gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))
   allocate(Pressure(gp%xsz(1), gp%xsz(2), gp%xsz(3)))
   allocate(divergence(gp %xsz(1), gp %xsz(2), gp %xsz(3)))

   divergence = 0.d0

   datadir = '/work2/06632/ryanhass/stampede2/LES/SpatialDecayDEBUG'
   call decomp_2d_read_one(1,u ,trim(datadir)//'/Run01_uVel_t012500.out',gp )
   call decomp_2d_read_one(1,v ,trim(datadir)//'/Run01_vVel_t012500.out',gp )
   call decomp_2d_read_one(1,wC,trim(datadir)//'/Run01_wVel_t012500.out',gp )
   
   call message(0, "Completed data read")

   allocate(spectC, spectE)
   call spectC%init("x", nx, ny, nz  , dx, dy, dz, "four", "2/3rd", 2,  fixOddball=.false., init_periodicInZ=.true. )
   call spectE%init("x", nx, ny, nz+1, dx, dy, dz, "four", "2/3rd", 2,  fixOddball=.false., init_periodicInZ=.false.)

   allocate(der)
   call der%init( gp, spectC%spectdecomp, gpE, spectE%spectdecomp, dz, scheme, isPeriodic, spectC)
   allocate(poiss)
   call poiss%init(dx, dy, dz, spectC, spectE, .false., two*pi, .true., gp, der, isPeriodic)
  

   call spectC%alloc_r2c_out(uhat)
   call spectC%alloc_r2c_out(vhat)
   call spectC%alloc_r2c_out(wChat)
   call spectE%alloc_r2c_out(what)
   
   call spectC%alloc_r2c_out(urhs)
   call spectC%alloc_r2c_out(vrhs)
   call spectE%alloc_r2c_out(wrhs)

   allocate(wChatz(spectC%spectdecomp%zsz(1),spectC%spectdecomp%zsz(2), &
     spectC%spectdecomp%zsz(3)))
   allocate(whatz (spectE%spectdecomp%zsz(1),spectE%spectdecomp%zsz(2), &
     spectE%spectdecomp%zsz(3)))
   
   allocate(cbuffzC (spectC%spectdecomp%zsz(1),spectC%spectdecomp%zsz(2), &
     spectC%spectdecomp%zsz(3),3))
   allocate(cbuffzE1(spectE%spectdecomp%zsz(1),spectE%spectdecomp%zsz(2), &
     spectE%spectdecomp%zsz(3)  ))
   allocate(cbuffyE(spectE%spectdecomp%ysz(1),spectE%spectdecomp%ysz(2), &
     spectE%spectdecomp%ysz(3)  ))
   allocate(cbuffyC(spectC%spectdecomp%ysz(1),spectC%spectdecomp%ysz(2), &
     spectC%spectdecomp%ysz(3)  ))

   call HITForce%init(trim(datadir)//"/HIT_forcing_input_test.dat", &
     spectC%spectdecomp, spectE%spectdecomp, spectC, cbuffyE, cbuffyC, &
     cbuffzE1, cbuffzC, tidStart)  
   call message(0, "Initialized HIT Forcing")
   call assert(HITForce%DomAspectRatioX == 10)
   call assert(HITForce%DomAspectRatioZ == 1)

   urhs = im0; vrhs = im0; wrhs = im0

   call spectC%fft(u, uhat)
   call spectC%fft(v, vhat)
   call spectC%fft(wC, wChat)
   call transpose_y_to_z(wChat, wChatz, spectC%spectdecomp)
   call der%interpz_C2E(wChatz, whatz, 0, 0)
   call transpose_z_to_y(whatz, what, spectE%spectdecomp)
   call spectE%ifft(what, w)
   call message(0, "Finished transforming velocities")

   !call poiss%getPressure(uhat, vhat, what, pressure)
   !call poiss%PressureProjection(uhat, vhat, what)
   !print*, uhat(3:4,34,43)
   
   call poiss%getPressureAndUpdateRHS(uhat, vhat, what, pressure)
   call poiss%DivergenceCheck(uhat, vhat, what, divergence)
   print*, maxval(abs(divergence))

   call spectC%dealias(uhat)
   call spectC%dealias(vhat)
   call transpose_y_to_z(what, whatz, spectE%spectdecomp)
   call spectC%dealias_edgeField(whatz)
   call transpose_z_to_y(whatz, what, spectE%spectdecomp)
   call poiss%DivergenceCheck(uhat, vhat, what, divergence)
 
   call message(0, "Finished Projection + Dealiasing")
   call HITForce%getRHS_HITforcing(urhs, vrhs, wrhs, uhat, vhat, what, .true.) 
   call message(0, "Finished computing HIT Force")
   
   call poiss%DivergenceCheck(urhs, vrhs, wrhs, divergence)
   print*, maxval(abs(divergence))

   call transpose_y_to_z(wrhs, whatz, spectE%spectdecomp)
   call der%interpz_E2C(whatz, wChatz, 0, 0)
   call transpose_z_to_y(wChatz, wChat, spectC%spectdecomp)
   call spectC%ifft(urhs, u)
   call spectC%ifft(vrhs, v)
   call spectC%ifft(wChat, wC)

   call decomp_2d_write_one(1,u ,trim(datadir)//"/u_HITforcing2.dat", gp)
   call decomp_2d_write_one(1,v ,trim(datadir)//"/v_HITforcing2.dat", gp)
   call decomp_2d_write_one(1,wC,trim(datadir)//"/w_HITforcing2.dat", gp)

   print*, "Max urhs:", p_maxval(abs(u))
   print*, "Max vrhs:", p_maxval(abs(v))
   print*, "Max wrhs:", p_maxval(abs(wC))

   deallocate(der)
   deallocate(poiss)
   deallocate(uhat, vhat, what)
   deallocate(u, v, w, wChat, wChatz, whatz)
   
   call MPI_Finalize(ierr)

end program 
