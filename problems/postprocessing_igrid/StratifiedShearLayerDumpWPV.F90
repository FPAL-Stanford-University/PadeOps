#include "StratifiedShearLayerWPV_files/PVroutines.F90"

program StratifiedShearLayerWPV
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use decomp_2d, only: nrank
   use timer, only: tic, toc
   use exits, only: message
   use PVroutines
   use PadePoissonMod, only: padepoisson
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, buff3
   real(rkind), dimension(:,:,:), allocatable :: buff4
   real(rkind), dimension(:,:,:), allocatable :: Fx, Fy, Fz
   real(rkind), dimension(:,:,:), allocatable :: ufluct, vfluct, w, T
   real(rkind) :: dx, dy, dz, Re = 3000.d0, Rib = 0.05d0, Tref = 100.d0
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   type(padepoisson) :: poiss 
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0
   integer :: idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
   logical :: isZPeriodic = .false.
   integer :: VizDump_Schedule = 0
   logical :: computeStokesPressure = .true.
   logical :: UseTrueWavenumbers = .true.
   integer, dimension(:), allocatable :: timesteps
   real(rkind), dimension(:), allocatable :: times
   integer :: nt

   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep, nx, ny, nz, Rib, NumericalSchemeVert, VizDump_Schedule

   call MPI_Init(ierr)
   call GETARG(1,inputfile)
   open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=99, NML=INPUT)
   close(unit=99)

   dx =     Lx/real(nx,rkind)
   dy =     Ly/real(ny,rkind)
   dz = two*Lz/real(nz,rkind)

 
   call ops%init(nx, ny, nz, dx, dy, dz, InputDir, OutputDir, RunID, isZPeriodic, NumericalSchemeVert)
   call poiss%init(dx, dy, dz, ops%spect, ops%spectE, computeStokesPressure, two*Lz, .false., ops%gp, ops%derZ, .false., UseTrueWavenumbers)


   ! Allocate all the needed memory
   call ops%allocate3DField(ufluct)
   call ops%allocate3DField(vfluct)
   call ops%allocate3DField(w)
   call ops%allocate3DField(T)
   call ops%allocate3DField(Fx)
   call ops%allocate3DField(Fy)
   call ops%allocate3DField(Fz)

   call ops%allocate3DField(buff1)
   call ops%allocate3DField(buff2)
   call ops%allocate3DField(buff3)   
   call ops%allocate3DField(buff4)

   if (VizDump_Schedule == 1) then
      call ops%Read_VizSummary(times, timesteps)
      nt = size(timesteps)
   else
      nt = (tstop - tstart)/tstep
      allocate(times(nt))
   end if
   

   call message(0,"Number of snapshots to read:", nt)

   idx = 1

   do while(idx <= nt)

      if (VizDump_Schedule == 1) then
         tidx = timesteps(idx)
      else
         tidx = tstart + tstep * (idx - 1)
         times(idx) = ops%getSimTime(tidx)
      end if
      

      call message(0, "Reading fields for tid:", TIDX)
      call tic()
      call ops%ReadField3D(buff1,"uVel",TIDX)
      call ops%ReadField3D(buff2,"vVel",TIDX)
      call ops%ReadField3D(w,"wVel",TIDX)
      call ops%ReadField3D(T,"potT",TIDX)
      call message(0, "Read simulation data at time:", times(idx))

      ! STEP 0: Compute Fluctuations

      call ops%getFluct_from_MeanZ(buff1,ufluct)
      call ops%getFluct_from_MeanZ(buff2,vfluct)


      ! STEP 1: Compute Wave-PV Decomposition of Pi and then Xi Velocities

      call ComputeF(ops,.TRUE.,ufluct,vfluct,w,T,Fx,Fy,Fz,buff1,buff2,buff3,buff4)
      
      !!!! DIAGNOSTIC
      call ops%WriteField3D(buff1, "oXPi", tidx) ! Omega_Pi before Decomposition
      call ops%WriteField3D(buff2, "oYPi", tidx)
      call ops%WriteField3D(buff3, "oZPi", tidx)
      !!!

      call poiss%PoissonSolver_HomogeneousNeumannBCz(Fx,buff1)
      call poiss%PoissonSolver_HomogeneousNeumannBCz(Fy,buff2)
      call poiss%PoissonSolver_HomogeneousNeumannBCz(Fz,buff3)

      call ops%WriteField3D(buff1, "uVPi", tidx)
      call ops%WriteField3D(buff2, "vVPi", tidx)
      call ops%WriteField3D(buff3, "wVPi", tidx)
      
      !!!! DIAGNOSTIC
      call ops%getCurl(buff1,buff2,buff3, Fx,Fy,Fz,1,1,1,1)
      call ops%WriteField3D(Fx, "vXPi", tidx) ! Omega_Pi after Decomposition
      call ops%WriteField3D(Fy, "vYPi", tidx)
      call ops%WriteField3D(Fz, "vZPi", tidx)
      !!!!

      
      call ComputeF(ops,.FALSE.,ufluct,vfluct,w,T,Fx,Fy,Fz,buff1,buff2,buff3,buff4)
      
      !!!! DIAGNOSTIC
      call ops%WriteField3D(buff1, "oXXi", tidx) ! Omega_Xi before Decomposition
      call ops%WriteField3D(buff2, "oYXi", tidx)
      call ops%WriteField3D(buff3, "oZXi", tidx)
      !!!

      call poiss%PoissonSolver_HomogeneousNeumannBCz(Fx,buff1)
      call poiss%PoissonSolver_HomogeneousNeumannBCz(Fy,buff2)
      call poiss%PoissonSolver_HomogeneousNeumannBCz(Fz,buff3)

      call ops%WriteField3D(buff1, "uVXi", tidx)
      call ops%WriteField3D(buff2, "vVXi", tidx)
      call ops%WriteField3D(buff3, "wVXi", tidx)
      
      !!!! DIAGNOSTIC
      call ops%getCurl(buff1,buff2,buff3, Fx,Fy,Fz,1,1,1,1)
      call ops%WriteField3D(Fx, "vXXi", tidx) ! Omega_Xi after Decomposition
      call ops%WriteField3D(Fy, "vYXi", tidx)
      call ops%WriteField3D(Fz, "vZXi", tidx)
      !!!!
      


      idx = idx + 1
      call toc()
   end do
   
   call MPI_Finalize(ierr)


end program
