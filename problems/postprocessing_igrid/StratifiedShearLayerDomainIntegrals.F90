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
   real(rkind), dimension(:,:,:), allocatable :: u, v, w, ufluct, vfluct, T, nuSGS, Tfluct
   real(rkind) :: dx, dy, dz, Re = 3000.d0, Rib = 0.05d0, Pr = 1.d0, Tref = 100.d0 
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 26.d0, Ly = 26.d0, Lz =12.d0
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
   call ops%allocate3DField(T)
   call ops%allocate3DField(Tfluct)
   call ops%allocate3DField(nuSGS)
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
      !call ops%ReadField3D(T,"potT",TIDX)
      !call ops%ReadField3D(nuSGS,"nSGS",TIDX)

      !T = Rib*(T - Tref)  ! Rescale Potential temperature to buoyancy variable: b 

      time(idx) = ops%getSimTime(tidx)
      call message(0, "Read simulation data at time:", time(idx))
     
      ! Compute TKE 
      call ops%getFluct_from_MeanZ(u,ufluct)
      call ops%getFluct_from_MeanZ(v,vfluct)
      buff2 = (ufluct*ufluct + vfluct*vfluct + w*w)
      TKE(idx) = 0.5d0*ops%getVolumeIntegral(buff2)
     
      idx = idx + 1
   end do 

   idx = idx - 1
   allocate(data2write(idx,2))
   data2write(:,1) = time(1:idx)
   data2write(:,2) = TKE(1:idx)

   if (nrank == 0) then
      call ops%WriteASCII_2D(data2write, "avgV")
   end if 

   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

