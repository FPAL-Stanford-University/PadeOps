program StratifiedShearLayerDomainIntegrals
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, buff3, buff4
   real(rkind), dimension(:,:,:), allocatable :: u, v, w 
   real(rkind) :: dx, dy, dz, Re = 3000.d0
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0
   integer :: ierr, tsnapshot = 0, NumericalSchemeVert = 1
   logical :: isZPeriodic = .false. 

   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tsnapshot, nx, ny, nz, Re, NumericalSchemeVert 
   
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
   
   call ops%allocate3DField(buff1)
   call ops%allocate3DField(buff2)
   call ops%allocate3DField(buff3)
   call ops%allocate3DField(buff4)
   


   tidx = tsnapshot
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
   call ops%GetGradient(u, buff1, buff2, buff3, 1, 1)
   buff4 = buff1*buff1 + buff2*buff2 + buff3*buff3

   call ops%GetGradient(v, buff1, buff2, buff3, 1, 1)
   buff4 = buff4 + buff1*buff1 + buff2*buff2 + buff3*buff3
   
   call ops%GetGradient(w, buff1, buff2, buff3, -1, -1)
   buff4 = buff4 + buff1*buff1 + buff2*buff2 + buff3*buff3

   buff4 = (1.d0/Re)*buff4


   call message(0,"Done computing dissipation rate")
   call ops%WriteField3D(buff4, "diss", tidx)
   call message(0, "Dissipation field written.")
   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

