program StratifiedShearLayerDomainIntegrals
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use decomp_2d, only: nrank 
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: buff1
   real(rkind), dimension(:,:,:), allocatable :: u, vfluct3D, wfluct, umean, ufluct2D, wfluct2D, ufluct3D, wfluct3D
   real(rkind) :: dx, dy, dz, Re = 3000.d0, Rib = 0.05d0, Pr = 1.d0, Tref = 100.d0 
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0
   real(rkind), dimension(10000) :: time, MKE, TKE2D, TKE3D 
   integer :: idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
   logical :: isZPeriodic = .false. 
   real(rkind), dimension(:,:), allocatable :: data2write
   integer :: VizDump_Schedule = 0
   integer, dimension(:), allocatable :: timesteps
   real(rkind), dimension(:), allocatable :: times
   integer :: nt 

   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep, nx, ny, nz, Re, Rib, Pr, Tref, NumericalSchemeVert, VizDump_Schedule 
   
   call MPI_Init(ierr)               
   call GETARG(1,inputfile)          
   open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=99, NML=INPUT)
   close(unit=99)

   dx =     Lx/real(nx,rkind) 
   dy =     Ly/real(ny,rkind) 
   dz = two*Lz/real(nz,rkind)

   ! Initialize the operator class
   call ops%init(nx, ny, nz, dx, dy, dz, InputDir, OutputDir, RunID, isZPeriodic, NUmericalSchemeVert)

   ! Allocate all the needed memory 
   call ops%allocate3DField(u)
   call ops%allocate3DField(vfluct3D)
   call ops%allocate3DField(wfluct)
   call ops%allocate3DField(umean)
   call ops%allocate3DField(ufluct2D)
   call ops%allocate3DField(wfluct2D)
   call ops%allocate3DField(ufluct3D)
   call ops%allocate3DField(wfluct3D)
   call ops%allocate3DField(buff1)

   idx = 1
   
   if (VizDump_Schedule == 1) then
      call ops%Read_VizSummary(times, timesteps)
      nt = size(timesteps)
   else
      nt = (tstop - tstart)/tstep
   end if
   
   do while(idx <= nt)
      
      if (VizDump_Schedule == 1) then
         tidx = timesteps(idx)
      else
         tidx = tstart + tstep * (idx - 1)
      end if
      
      call message(0, "Reading fields for tid:", TIDX)
      call tic()
      call ops%ReadField3D(u,"uVel",TIDX)
      call ops%ReadField3D(vfluct3D,"vVel",TIDX)
      call ops%ReadField3D(wfluct,"wVel",TIDX)


      time(idx) = ops%getSimTime(tidx)
      call message(0, "Read simulation data at time:", time(idx))

      ! Mean Fields
      call ops%getFluct_from_MeanZ(u,buff1)
      umean = u - buff1 ! mean
      buff1 = 0.5d0*umean*umean
      MKE(idx) = ops%getVolumeIntegral(buff1)
      
      ! 2D Fluctuations
      buff1 = u - umean
      call ops%TakeMean_y(buff1,ufluct2D)
      call ops%TakeMean_y(wfluct,wfluct2D)

      ! 3D Fluctuations
      ufluct3D = u - umean - ufluct2D
      wfluct3D = wfluct - wfluct2D

      
      ! TKE
      buff1 = (ufluct2D*ufluct2D + wfluct2D*wfluct2D)
      TKE2D(idx) = 0.5d0*ops%getVolumeIntegral(buff1)
      buff1 = (ufluct3D*ufluct3D + vfluct3D*vfluct3D + wfluct3D*wfluct3D)
      TKE3D(idx) = 0.5d0*ops%getVolumeIntegral(buff1)


      if (nrank == 0) then
         print*, time(idx), TKE2D(idx), TKE3D(idx)
      end if 
     
      idx = idx + 1
   end do 

   idx = idx - 1
   allocate(data2write(idx,9))
   data2write(:,1) = time(1:idx)
   data2write(:,2) = MKE(1:idx)
   data2write(:,3) = TKE2D(1:idx)
   data2write(:,4) = TKE3D(1:idx)

   if (nrank == 0) then
      call ops%WriteASCII_2D(data2write, "TKE")
   end if 

   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

