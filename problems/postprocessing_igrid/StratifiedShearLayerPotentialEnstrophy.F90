program StratifiedShearLayerDomainIntegrals
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, buff3, buff4, buff5, buff6, buff7
   real(rkind), dimension(:,:,:), allocatable :: u, v, w, potT 
   real(rkind) :: dx, dy, dz, Re = 3000.d0, Rib = 0.2d0
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0
   integer :: idx, ierr, tstart = 0, tstop, tstep, NumericalSchemeVert = 1 
   logical :: isZPeriodic = .false. 
   integer :: nt

   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep, nx, ny, nz, Re, Rib, NumericalSchemeVert, tidx  
   
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
   call ops%allocate3DField(potT)

   call ops%allocate3DField(buff1)
   call ops%allocate3DField(buff2)
   call ops%allocate3DField(buff3)
   call ops%allocate3DField(buff4) 
   call ops%allocate3DField(buff5)
   call ops%allocate3DField(buff6)
   call ops%allocate3DField(buff7)
   
   nt = (tstop - tstart)/tstep + 1

   call message(0,"Number of snapshots to read:", nt)

   tidx = tstart
   idx = 1
   do while(tidx <= tstop)
      
      call message(0, "Reading fields for tid:", TIDX)
      call tic()
      call ops%ReadField3D(u,"uVel",TIDX)
      call ops%ReadField3D(v,"vVel",TIDX)
      call ops%ReadField3D(w,"wVel",TIDX) 
      call ops%ReadField3D(potT,"potT",TIDX)
      
      call ops%getCurl(u,v,w,   buff1,buff2,buff3,1,1,1,1) 
      call ops%getGradient(potT,buff4,buff5,buff6,1,1)
      
      ! Potential Enstrophy (Total)
      buff7 = buff1 * buff4 + buff2 * buff5 + buff3 * buff6
      buff7 = buff7 * buff7
      !call ops%writeField3D(buff7,"pEns",TIDX)
      call ops%dump_plane(buff7,1,nx/2,TIDX,"pEns")
      call ops%dump_plane(buff7,2,ny/2,TIDX,"pEns")
      call ops%dump_plane(buff7,3,nz/2,TIDX,"pEns")

      ! Potential Enstrophy with background P
      !call ops%getFluct_from_MeanZ(buff6,buff7)
      !buff7 = buff7 - buff6
      !buff7 = buff7 * buff3
      !buff7 = buff7 * buff7
      !call ops%writeField3D(buff7,"ptE1",TIDX)
      
      ! Potential Enstrophy with background Vorticity
      ! NOTE: Ignore cross term between 1 and 2
      !call ops%getFluct_from_MeanZ(buff2,buff7)
      !buff7 = buff7 - buff2
      !buff7 = buff7 * buff5
      !buff7 = buff7 * buff7
      !call ops%writeField3D(buff7,"ptE2",TIDX)
     
      ! Perturbation Enstrophy
      call ops%getFluct_from_MeanZ(buff2,buff7)
      buff7 = buff1*buff1 + buff7*buff7 + buff3*buff3
      !call ops%writeField3D(buff7,"Ens1",TIDX)
      call ops%dump_plane(buff7,1,nx/2,TIDX,"Ens1")
      call ops%dump_plane(buff7,2,ny/2,TIDX,"Ens1")
      call ops%dump_plane(buff7,3,nz/2,TIDX,"Ens1")
      
      ! Enstrophy (Total)
      buff7 = buff1*buff1 + buff2*buff2 + buff3*buff3
      !call ops%writeField3D(buff7,"Enst",TIDX)
      call ops%dump_plane(buff7,1,nx/2,TIDX,"Enst")
      call ops%dump_plane(buff7,2,ny/2,TIDX,"Enst")
      call ops%dump_plane(buff7,3,nz/2,TIDX,"Enst")
     
      ! (GradT)^2
      call ops%getFluct_from_MeanZ(buff6,buff7)
      buff7 = buff4*buff4 + buff5*buff5 + buff7*buff7
      !call ops%writeField3D(buff7,"Enst",TIDX)
      call ops%dump_plane(buff7,1,nx/2,TIDX,"dpT2")
      call ops%dump_plane(buff7,2,ny/2,TIDX,"dpT2")
      call ops%dump_plane(buff7,3,nz/2,TIDX,"dpT2")

      tidx = tidx + tstep

      call toc()
   
   end do


   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

