program StratifiedShearLayerDomainIntegrals
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use decomp_2d, only: nrank 
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: u, T
   real(rkind) :: dx, dy, dz 
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0
   integer :: idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
   logical :: isZPeriodic = .false.
   real(rkind), dimension(:,:), allocatable :: umn_set, Tmn_set
   real(rkind), dimension(:,:), allocatable :: times
   integer :: nt 

   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep, nx, ny, nz, NumericalSchemeVert 
   
   call MPI_Init(ierr)               
   call GETARG(1,inputfile)          
   open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=99, NML=INPUT)
   close(unit=99)

   dx =     Lx/real(nx,rkind) 
   dy =     Ly/real(ny,rkind) 
   dz = two*Lz/real(nz,rkind)

   ! Initialize the operator class
   call ops%init(nx, ny, nz, dx, dy, dz, InputDir, OutputDir, RunID, isZPeriodic, NumericalSchemeVert)

   ! Allocate all the needed memory 
   call ops%allocate3DField(u)
   call ops%allocate3DField(T)
   

   nt = (tstop - tstart)/tstep + 1

   call message(0,"Number of snapshots to read:", nt)

   !if (nrank == 0) then
      allocate(umn_set(nz,nt))  
      allocate(Tmn_set(nz,nt))  

      allocate(times(nt,1))
   !end if 

   tidx = tstart
   idx = 1
   do while(tidx <= tstop)
      call message(0, "Reading fields for tid:", TIDX)
      call tic()
      call ops%ReadField3D(u,"uVel",TIDX)
      call ops%ReadField3D(T,"potT",TIDX)
      times(idx,1) = ops%getSimTime(tidx)
      call message(0, "Read simulation data at time:", times(idx,1))

      
      ! STEP 0: Compute means
      call ops%TakeMean_xy(u,umn_set(:,idx))
      call ops%TakeMean_xy(T,Tmn_set(:,idx))
      

      tidx = tidx + tstep
      idx = idx + 1
      call toc()
   end do 

   if (nrank == 0) then
      call ops%WriteASCII_2D(umn_set, "umnZ")
      call ops%WriteASCII_2D(Tmn_set, "TmnZ")
      
      call ops%WriteASCII_2D(times, "time")
   end if 

   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

