program HIT_AD_get3Dmeans
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use decomp_2d, only: nrank 
   use timer, only: tic, toc
   use exits, only: message
   implicit none
   real(rkind), dimension(:,:,:), allocatable :: usum, vsum, wsum 
   real(rkind), dimension(:,:,:), allocatable :: u, v, w
   real(rkind) :: dx, dy, dz
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 10.d0*pi, Ly = 2.d0*pi, Lz = 2.d0
   integer :: idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 2
   logical :: isZPeriodic = .true. 
   integer :: nt 

   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep, nx, ny, nz
   
   call MPI_Init(ierr)               
   call GETARG(1,inputfile)          
   open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=99, NML=INPUT)
   close(unit=99)

   dx = Lx/real(nx,rkind) 
   dy = Ly/real(ny,rkind) 
   dz = Lz/real(nz,rkind)

   ! Initialize the operator class
   call ops%init(nx, ny, nz, dx, dy, dz, InputDir, OutputDir, RunID, isZPeriodic, NUmericalSchemeVert)

   ! Allocate all the needed memory 
   call ops%allocate3DField(u)
   call ops%allocate3DField(v)
   call ops%allocate3DField(w)
   
   call ops%allocate3DField(usum)
   call ops%allocate3DField(vsum)
   call ops%allocate3DField(wsum)
  
   usum = 0.d0 
   vsum = 0.d0 
   wsum = 0.d0 

   nt = (tstop - tstart)/tstep + 1

   call message(0,"Number of snapshots to read:", nt)

   tidx = tstart
   idx = 0
   do while(tidx <= tstop)
      call message(0, "Reading fields for tid:", TIDX)
      call tic()
      call ops%ReadField3D(u,"uVel",TIDX)
      call ops%ReadField3D(v,"vVel",TIDX)
      call ops%ReadField3D(w,"wVel",TIDX)

      usum = usum + u
      vsum = vsum + v
      wsum = wsum + w
      
      tidx = tidx + tstep
      idx = idx + 1
      call toc()
   end do 

   usum = usum/real(idx,rkind)
   vsum = vsum/real(idx,rkind)
   wsum = wsum/real(idx,rkind)

   call ops%WriteField3D(usum,"uMmn",tidx)
   call ops%WriteField3D(vsum,"vMmn",tidx)
   call ops%WriteField3D(wsum,"wMmn",tidx)
   
   
   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

