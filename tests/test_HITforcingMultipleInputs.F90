program test_HITforcing
   use kind_parameters, only: rkind, clen
   use mpi
   use decomp_2d
   use decomp_2d_io
   use exits, only: message 
   use spectralMod, only: spectral
   use constants, only: pi, im0, two
   use PadePoissonMod, only: padepoisson
   use PadeDerOps, only: Pade6stagg
   use forcingMod, only: HIT_shell_forcing 
   use reductions, only: p_maxval, p_mean, p_sum
   use fortran_assert, only: assert

   implicit none
   real   (rkind), dimension(:,:,:), allocatable :: Pressure, u, v, w, wC, &
                                                    fx, fy, fz 
   complex(rkind), dimension(:,:,:), allocatable :: uhat, vhat, what, wChat, wcHatz, whatz, urhs, vrhs, wrhs
   real(rkind),    dimension(:,:,:), allocatable :: divergence
   logical, parameter :: isPeriodic = .true. 
   integer, parameter :: scheme = 2
   complex(rkind), dimension(:,:,:,:), allocatable :: cbuffzC
   complex(rkind), dimension(:,:,:), allocatable :: cbuffzE1, cbuffyE, cbuffyC
   real(rkind), dimension(:,:,:), allocatable :: rbuffxC

   integer, parameter :: nx = 128, ny = 128, nz = 128
   real(rkind), parameter :: dx = two*pi/real(nx,rkind), dy = two*pi/real(ny,rkind), dz = two*pi/real(nz,rkind) 
   type(decomp_info) :: gp, gpE
   integer :: ierr 
   type(spectral), allocatable :: spectC, spectE
   type(padepoisson), allocatable :: poiss
   type(Pade6stagg), allocatable :: der
   type(HIT_shell_forcing) :: HITForce
   integer :: tidStart = 0 ! Time index to initialize seeding in HITForce
   integer :: tid = 0
   character(len=clen) :: mssg
   integer :: i, n, Ninputs = 100, Nrealizations = 100
   integer :: dumpFreq = 1
   character(len=clen) :: fname, outputdir, inputdir
   real(rkind) :: energyInjectionR

   call mpi_init(ierr)
   call decomp_2d_init(nx,ny,nz,0,0)
   call get_decomp_info(gp)
   call decomp_info_init(nx,ny,nz+1, gpE)

   allocate(u (gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(v (gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(wC(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(w (gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))
   allocate(fx(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(fy(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(fz(gp %xsz(1), gp %xsz(2), gp %xsz(3)))
   allocate(Pressure(gp%xsz(1), gp%xsz(2), gp%xsz(3)))
   allocate(divergence(gp %xsz(1), gp %xsz(2), gp %xsz(3)))

   divergence = 0.d0

   inputdir = '/anvil/scratch/x-ryanhass/HITforced/HIT128_DNS2/output'
   outputdir = '/anvil/projects/x-atm170028/ryanhass/HITforced/HITforcingDebug/HIT128_DNS2_InputData/forcingData' 
   

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

   allocate(wChatz(spectC%spectdecomp%zsz(1),spectC%spectdecomp%zsz(2), spectC%spectdecomp%zsz(3)))
   allocate(whatz (spectE%spectdecomp%zsz(1),spectE%spectdecomp%zsz(2), spectE%spectdecomp%zsz(3)))
   
   allocate(cbuffzC (spectC%spectdecomp%zsz(1),spectC%spectdecomp%zsz(2), spectC%spectdecomp%zsz(3),3))
   allocate(cbuffzE1(spectE%spectdecomp%zsz(1),spectE%spectdecomp%zsz(2), spectE%spectdecomp%zsz(3)  ))
   allocate(cbuffyE(spectE%spectdecomp%ysz(1),spectE%spectdecomp%ysz(2), spectE%spectdecomp%ysz(3)  ))
   allocate(cbuffyC(spectC%spectdecomp%ysz(1),spectC%spectdecomp%ysz(2), spectC%spectdecomp%ysz(3)  ))
   
   allocate(rbuffxC(gp%xsz(1),gp%xsz(2),gp%xsz(3)))

   call HITForce%init("/anvil/projects/x-atm170028/ryanhass/HITforced"//&
     "/HITforcingDebug/HIT128_DNS2_InputData/HIT_forcing_input_test.dat", gp, &
     spectC%spectdecomp, spectE%spectdecomp, spectC, cbuffyE, cbuffyC, &
     cbuffzE1, cbuffzC, rbuffxC, tidStart, poiss)  
   call message(0, "Initialized HIT Forcing")

   urhs = im0; vrhs = im0; wrhs = im0
  
   do n = 1,Ninputs 
     tid = tid + 100
     write(fname,'(A,I6.6,A)')trim(inputdir)//'/Run02_uVel_t',tid,'.out'
     call decomp_2d_read_one(1,u ,trim(fname),gp )
     write(fname,'(A,I6.6,A)')trim(inputdir)//'/Run02_vVel_t',tid,'.out'
     call decomp_2d_read_one(1,v ,trim(fname),gp )
     write(fname,'(A,I6.6,A)')trim(inputdir)//'/Run02_wVel_t',tid,'.out'
     call decomp_2d_read_one(1,wC,trim(fname),gp )
     call message(0, "Completed data read")
     call message(0, "max u:", p_maxval(maxval(u)))

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

     call spectC%dealias(uhat)
     call spectC%dealias(vhat)
     call transpose_y_to_z(what, whatz, spectE%spectdecomp)
     call spectC%dealias_edgeField(whatz)
     call transpose_z_to_y(whatz, what, spectE%spectdecomp)
     call poiss%DivergenceCheck(uhat, vhat, what, divergence)
 
     call message(0, "Finished Projection + Dealiasing")
     do i = 1,Nrealizations
       urhs = im0; vrhs = im0; wrhs = im0
       call HITForce%getRHS_HITforcing(urhs, vrhs, wrhs, uhat, vhat, what, .true.)
       call message(0, "seed3:", HITForce%seed3)
       
       write(mssg,'(A,I5,A,I5,A,I4,A,I4)')'Realization ',i,' of ',Nrealizations, &
         '; Input ',n,' of ',Ninputs
       call message(0, trim(mssg))
   
       call poiss%DivergenceCheck(urhs, vrhs, wrhs, divergence)
       call message(0, "Max divergence: ", p_maxval(abs(divergence)))
  
       call transpose_y_to_z(wrhs, whatz, spectE%spectdecomp)
       call der%interpz_E2C(whatz, wChatz, 0, 0)
       call transpose_z_to_y(wChatz, wChat, spectC%spectdecomp)
       call spectC%ifft(urhs, fx)
       call spectC%ifft(vrhs, fy)
       call spectC%ifft(wChat, fz)

       energyInjectionR = p_sum(sum(u*fx + v*fy + wC*fz)*dx*dy*dz/((two*pi)**3.d0))
       call message(0, "Physical space energy injection rate:", energyInjectionR)
       
       call message(0, "Max urhs: ",p_maxval(abs(fx)))
       call message(0, "Max vrhs: ",p_maxval(abs(fy)))
       call message(0, "Max wrhs: ",p_maxval(abs(fz)))
       call message(0, "fxfy: ", p_mean(fx*fy))
       call message(0, "fxfz: ", p_mean(fx*fz))
       call message(0, "fyfz: ", p_mean(fy*fz))
       call message(0, "fxfx: ", p_mean(fx*fx))
       call message(0, "fyfy: ", p_mean(fy*fy))
       call message(0, "fzfz: ", p_mean(fz*fz))

       if (mod(n,dumpFreq) == 0 .and. i == 1) then
         call HITForce%dumpForcing(outputdir,0,tid)
       end if
       if (i == Nrealizations .and. n == Ninputs) then
         print*, "Max urhs:", p_maxval(abs(fx))
         print*, "Max vrhs:", p_maxval(abs(fy))
         print*, "Max wrhs:", p_maxval(abs(fz))
       end if
       fx = 0.d0; fy = 0.d0; fz = 0.d0
     end do
   end do

   deallocate(der)
   deallocate(poiss)
   deallocate(uhat, vhat, what)
   deallocate(fx, fy, fz, u, v, w, wChat, wChatz, whatz)
   call HITForce%destroy()   
   call MPI_Finalize(ierr)

end program 
