program StratifiedShearLayerDumpQcrit
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: dudx, dudy, dudz 
   real(rkind), dimension(:,:,:), allocatable :: dwdx, dwdy, dwdz 
   real(rkind), dimension(:,:,:), allocatable :: dvdx, dvdy, dvdz 
   real(rkind), dimension(:,:,:), allocatable :: u, v, w, buff1 
<<<<<<< HEAD
   real(rkind) :: dx, dy, dz, Re = 3000.d0
=======
   real(rkind) :: dx, dy, dz, Re = 3000.d0, Rib = 0.05d0
>>>>>>> origin/igridSGS
   integer :: nx, ny, nz, nt, RunID, TIDX, tstart=0, tstop=0, tstep=0, VizDump_Schedule=0
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0
   integer :: ierr, tsnapshot = 0, NumericalSchemeVert = 1, idx 
   logical :: isZPeriodic = .false. 
   integer, dimension(:), allocatable :: timesteps
   real(rkind), dimension(:), allocatable :: times

<<<<<<< HEAD
   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep, nx, ny, nz, Re, NumericalSchemeVert, VizDump_Schedule 
=======
   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep, nx, ny, nz, Re, Rib, NumericalSchemeVert, VizDump_Schedule 
>>>>>>> origin/igridSGS
   
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
   call ops%allocate3DField(v)
   call ops%allocate3DField(w)
   
   call ops%allocate3DField(dudx)
   call ops%allocate3DField(dudy)
   call ops%allocate3DField(dudz)
   call ops%allocate3DField(dvdx)
   call ops%allocate3DField(dvdy)
   call ops%allocate3DField(dvdz)
   call ops%allocate3DField(dwdx)
   call ops%allocate3DField(dwdy)
   call ops%allocate3DField(dwdz)
   call ops%allocate3DField(buff1)
   
   if (VizDump_Schedule == 1) then
      call ops%Read_VizSummary(times, timesteps)
      nt = size(timesteps)
   else
      nt = (tstop - tstart)/tstep
   end if

   idx = 1

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
      
      call ops%getFluct_from_MeanZ(v,buff1)
      v = buff1
      
      call ops%getFluct_from_MeanZ(u,buff1)
      u = buff1

      ! Get gradients and add norms
      call ops%GetGradient(u, dudx, dudy, dudz, 1, 1)
      call ops%GetGradient(v, dvdx, dvdy, dvdz, 1, 1)
      call ops%GetGradient(w, dwdx, dwdy, dwdz, -1, -1)

      buff1 = 0.5d0*(dudx*dudx  + dvdy*dvdy + dwdz*dwdz)
      buff1 = buff1 + dudy*dvdx + dudz*dwdx + dvdz*dwdy
      buff1 = -buff1

      call ops%WriteField3D(buff1, "qcri", tidx)

      idx = idx + 1
  end do 

  call ops%destroy()
  call MPI_Finalize(ierr)           


end program 

