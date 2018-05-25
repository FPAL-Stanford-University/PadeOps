program StratifiedShearLayerDomainIntegrals
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use decomp_2d, only: nrank 
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: buff2, buff3, buff4
   real(rkind), dimension(:,:,:), allocatable :: u, v, w, ufluct, vfluct, T, Tfluct
   real(rkind) :: dx, dy, dz, Re, Rib, Pr = 1.d0, Tref = 100.d0 
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   real(rkind) :: Lx, Ly, Lz
   real(rkind), dimension(10000) :: P, B, D, Dv, Dsgs, time, IEL, MKE, TKE 
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
   call ops%allocate3DField(ufluct)
   call ops%allocate3DField(v)
   call ops%allocate3DField(vfluct)
   call ops%allocate3DField(w)
   call ops%allocate3DField(buff2)
   call ops%allocate3DField(buff3)
   call ops%allocate3DField(buff4)

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
      call ops%ReadField3D(v,"vVel",TIDX)
      call ops%ReadField3D(w,"wVel",TIDX)

      time(idx) = ops%getSimTime(tidx)
      call message(0, "Read simulation data at time:", time(idx))
  
      ! STEP 4: Compute the production
      call ops%getFluct_from_MeanZ(u,ufluct)
      buff3 = u - ufluct ! mean
      call ops%ddz(buff3,buff2,1,1) ! dUdz (no stress BC)
      buff3 = 0.5d0*buff3*buff3
      MKE(idx) = ops%getVolumeIntegral(buff3)
      buff3 = -ufluct*w      ! u'w' (since wmean = 0)
      buff2 = buff2*buff3 
      P(idx) = ops%getVolumeIntegral(buff2)
      
      call ops%getFluct_from_MeanZ(v,vfluct)
      buff3 = v - vfluct ! mean
      call ops%ddz(buff3,buff2,1,1) ! dVdz (no stress BC)
      buff3 = 0.5d0*buff3*buff3
      MKE(idx) = MKE(idx) + ops%getVolumeIntegral(buff3)
      buff3 = -vfluct*w      ! u'w' (since wmean = 0)
      buff2 = buff2*buff3 
      P(idx) = P(idx) + ops%getVolumeIntegral(buff2)

      ! STEP 5a: Compute TKE
      buff2 = (ufluct*ufluct + vfluct*vfluct + w*w)
      TKE(idx) = 0.5d0*ops%getVolumeIntegral(buff2)

      call toc()

      if (nrank == 0) then
         print*, "ddt_TKE:", P(idx) + B(idx) - D(idx)
      end if 
     
      idx = idx + 1
   end do 

   idx = idx - 1
   allocate(data2write(idx,2))
   data2write(:,1) = time(1:idx)
   data2write(:,2) = TKE(1:idx)

   if (nrank == 0) then
      call ops%WriteASCII_2D(data2write, "TKE")
   end if 

   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

