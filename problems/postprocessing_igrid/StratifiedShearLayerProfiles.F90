program StratifiedShearLayerProfiles
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use decomp_2d, only: nrank
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, buff3
   real(rkind), dimension(:,:,:), allocatable :: buff4, buff5, buff6
   real(rkind), dimension(:,:,:), allocatable :: buff7, buff8
   real(rkind), dimension(:,:,:), allocatable :: u, v, w, ufluct, vfluct, T,  nuSGS
   real(rkind) :: dx, dy, dz, Re = 3000.d0, Rib = 0.05d0, Tref = 100.d0
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0
   integer :: idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
   logical :: isZPeriodic = .false.
   real(rkind), dimension(:,:), allocatable :: Prod, Buoy, Diss
   real(rkind), dimension(:,:), allocatable :: OmZ, PiZ, XiZ
   real(rkind), dimension(:,:), allocatable :: uXiZ, vXiZ, wXiZ, uPiZ, vPiZ, wPiZ
   real(rkind), dimension(:,:), allocatable :: OmXXiZ, OmYXiZ, OmZXiZ
   real(rkind), dimension(:,:), allocatable :: KEXiZ, KEPiZ
   real(rkind), dimension(:,:), allocatable :: OmXPiZ, OmYPiZ, OmZPiZ
   integer :: VizDump_Schedule = 0
   integer, dimension(:), allocatable :: timesteps
   real(rkind), dimension(:), allocatable :: times
   real(rkind), dimension(:,:), allocatable :: timewrite
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

   ! Initialize the operator class
   call ops%init(nx, ny, nz, dx, dy, dz, InputDir, OutputDir, RunID, isZPeriodic, NUmericalSchemeVert)

   ! Allocate all the needed memory
   call ops%allocate3DField(u)
   call ops%allocate3DField(ufluct)
   call ops%allocate3DField(v)
   call ops%allocate3DField(vfluct)
   call ops%allocate3DField(w)
   call ops%allocate3DField(T)
   call ops%allocate3DField(nuSGS)

   call ops%allocate3DField(buff1)
   call ops%allocate3DField(buff2)
   call ops%allocate3DField(buff3)
   call ops%allocate3DField(buff4)
   call ops%allocate3DField(buff5)
   call ops%allocate3DField(buff6)
   call ops%allocate3DField(buff7)
   call ops%allocate3DField(buff8)



   if (VizDump_Schedule == 1) then
      call ops%Read_VizSummary(times, timesteps)
      nt = size(timesteps)
   else
      nt = (tstop - tstart)/tstep
   end if



   call message(0,"Number of snapshots to read:", nt)

   !if (nrank == 0) then
      allocate(Prod(nz,nt))
      allocate(Buoy(nz,nt))
      allocate(Diss(nz,nt))
      allocate(OmZ(nz,nt))
      allocate(PiZ(nz,nt))
      allocate(XiZ(nz,nt))
      allocate(uXiZ(nz,nt))
      allocate(vXiZ(nz,nt))
      allocate(wXiZ(nz,nt))
      allocate(uPiZ(nz,nt))
      allocate(vPiZ(nz,nt))
      allocate(wPiZ(nz,nt))
      allocate(OmXXiZ(nz,nt))
      allocate(OmYXiZ(nz,nt))
      allocate(OmZXiZ(nz,nt))
      allocate(KEXiZ(nz,nt))
      allocate(KEPiZ(nz,nt))
      allocate(OmXPiZ(nz,nt))
      allocate(OmYPiZ(nz,nt))
      allocate(OmZPiZ(nz,nt))
      allocate(timewrite(nt,1))
   !end if

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
      call ops%ReadField3D(T,"potT",TIDX)
      call message(0, "Read simulation data at time:", times(idx))
      if (tidx == 0) then
        nuSGS = 0
      else
        call ops%ReadField3D(nuSGS,"nSGS",TIDX)
      end if

      T = Rib*(T - Tref)  ! Rescale Potential temperature to buoyancy variable: b

      ! STEP 0: Compute Fluctuations
      call ops%getFluct_from_MeanZ(v,vfluct)
      call ops%getFluct_from_MeanZ(u,ufluct)


      ! STEP 1: Compute Production, Buoyancy
      buff3 = u - ufluct
      call ops%ddz(buff3,buff2, 1, 1)
      buff3 = -ufluct*w
      buff2 = buff2*buff3
      call ops%TakeMean_xy(buff2,Prod(:,idx))

      call ops%getFluct_from_MeanZ(T,buff1)
      buff2 = buff1*w
      call ops%TakeMean_xy(buff2,Buoy(:,idx))


      ! STEP 2: Compute Dissipation
      buff3 = 0.d0

      call ops%ddx(ufluct,buff2)
      buff3 = buff2*buff2
      call ops%ddy(ufluct,buff2)
      buff3 = buff3 + buff2*buff2
      call ops%ddz(ufluct,buff2, 1, 1)
      buff3 = buff3 + buff2*buff2

      call ops%ddx(vfluct,buff2)
      buff3 = buff3 + buff2*buff2
      call ops%ddy(vfluct,buff2)
      buff3 = buff3 + buff2*buff2
      call ops%ddz(vfluct,buff2, 1, 1)
      buff3 = buff3 + buff2*buff2

      call ops%ddx(w,buff2)
      buff3 = buff3 + buff2*buff2
      call ops%ddy(w,buff2)
      buff3 = buff3 + buff2*buff2
      call ops%ddz(w,buff2, -1, -1)  ! no-penetration BCs
      buff3 = buff3 + buff2*buff2
      
      buff8 = (1.d0/Re)*buff3

      ! SGS sink term
      ! s11*s11
      call ops%ddx(ufluct,buff3)
      buff2 = buff3*buff3
      
      ! 2*s12*s12
      call ops%ddy(ufluct,buff3)
      call ops%ddx(vfluct,buff4)
      buff3 = 0.5d0*(buff3 + buff4)
      buff2 = buff2 + 2.d0*buff3*buff3 

      ! 2*s13*s13 
      call ops%ddz(ufluct,buff3, 1, 1)
      call ops%ddx(w     ,buff4)
      buff3 = 0.5d0*(buff3 + buff4)
      buff2 = buff2 + 2.d0*buff3*buff3 

      ! 2*s23*s23 
      call ops%ddz(vfluct,buff3, 1, 1)
      call ops%ddy(w     ,buff4)
      buff3 = 0.5d0*(buff3 + buff4)
      buff2 = buff2 + 2.d0*buff3*buff3 

      ! s22*s22
      call ops%ddy(vfluct,buff3)
      buff2 = buff2 + buff3*buff3

      ! s33*s33
      call ops%ddz(w,buff3, -1, -1)
      buff2 = buff2 + buff3*buff3
      buff2 = 2.d0*nuSGS*buff2 + buff8
      call ops%TakeMean_xy(buff2,Diss(:,idx))



      ! STEP 3: Compute Wave-PV Decomposition

      call ops%getCurl(ufluct,vfluct,w, buff1,buff2,buff3,1,1,1,1)
      call ops%getGradient(T,           buff4,buff5,buff6,1,1)

      call ops%TakeMean_xy(buff3,OmZ(:,idx))

      buff7 = buff1 * buff4 + buff2 * buff5 + buff3 * buff6 ! Pi
      call ops%TakeMean_xy(buff7,PiZ(:,idx))

      buff8 = buff4*buff4 + buff5*buff5 + buff6*buff6
      buff8 = buff8**(1.d0/2.d0) ! |gradT|

      buff1 = buff1 - buff7*buff4/buff8**2
      buff2 = buff2 - buff7*buff5/buff8**2
      buff3 = buff3 - buff7*buff6/buff8**2

      buff1 = buff1*buff1 + buff2*buff2 + buff3*buff3
      buff1 = buff1**(1.d0/2.d0) * buff8 ! Xi
      call ops%TakeMean_xy(buff1,XiZ(:,idx))


      buff7 = ufluct*buff4 + vfluct*buff5 + w*buff6
      buff7 = buff7 / buff8**2
      buff1 = buff7*buff4
      call ops%TakeMean_xy(buff1,uXiZ(:,idx))
      buff2 = buff7*buff5
      call ops%TakeMean_xy(buff2,vXiZ(:,idx))
      buff3 = buff7*buff6
      call ops%TakeMean_xy(buff3,wXiZ(:,idx))

      buff8 = 0.5d0*(buff1*buff1+buff2*buff2+buff3*buff3)
      call ops%TakeMean_xy(buff8,KEXiZ(:,idx))

      call ops%getCurl(buff1,buff2,buff3, buff4,buff5,buff6,1,1,1,1)
      call ops%TakeMean_xy(buff4,OmXXiZ(:,idx))
      call ops%TakeMean_xy(buff5,OmYXiZ(:,idx))
      call ops%TakeMean_xy(buff6,OmZXiZ(:,idx))

      buff4 = ufluct - buff1
      call ops%TakeMean_xy(buff4,uPiZ(:,idx))
      buff5 = vfluct - buff2
      call ops%TakeMean_xy(buff5,vPiZ(:,idx))
      buff6 = w - buff3
      call ops%TakeMean_xy(buff6,wPiZ(:,idx))

      buff8 = 0.5d0*(buff4*buff4+buff5*buff5+buff6*buff6)
      call ops%TakeMean_xy(buff8,KEPiZ(:,idx))

      call ops%getCurl(buff4,buff5,buff6, buff1,buff2,buff3,1,1,1,1)
      call ops%TakeMean_xy(buff1,OmXPiZ(:,idx))
      call ops%TakeMean_xy(buff2,OmYPiZ(:,idx))
      call ops%TakeMean_xy(buff3,OmZPiZ(:,idx))


      idx = idx + 1
      call toc()
   end do

   if (nrank == 0) then
      call ops%WriteASCII_2D(Prod, "prod")
      call ops%WriteASCII_2D(Buoy, "buoy")
      call ops%WriteASCII_2D(Diss, "diss")

      call ops%WriteASCII_2D(OmZ, "OmeZ")
      call ops%WriteASCII_2D(PiZ, "PiiZ")
      call ops%WriteASCII_2D(XiZ, "XiiZ")

      call ops%WriteASCII_2D(uXiZ, "uXiZ")
      call ops%WriteASCII_2D(vXiZ, "vXiZ")
      call ops%WriteASCII_2D(wXiZ, "wXiZ")

      call ops%WriteASCII_2D(uPiZ, "uPiZ")
      call ops%WriteASCII_2D(vPiZ, "vPiZ")
      call ops%WriteASCII_2D(wPiZ, "wPiZ")

      call ops%WriteASCII_2D(OmXXiZ, "OmXXiZ")
      call ops%WriteASCII_2D(OmYXiZ, "OmYXiZ")
      call ops%WriteASCII_2D(OmZXiZ, "OmZXiZ")

      call ops%WriteASCII_2D(KEXiZ, "KEXiZ")
      call ops%WriteASCII_2D(KEPiZ, "KEPiZ")

      call ops%WriteASCII_2D(OmXPiZ, "OmXPiZ")
      call ops%WriteASCII_2D(OmYPiZ, "OmYPiZ")
      call ops%WriteASCII_2D(OmZPiZ, "OmZPiZ")


      timewrite(:,1) = times
      call ops%WriteASCII_2D(timewrite, "time")
   end if

   call ops%destroy()
   call MPI_Finalize(ierr)


end program
