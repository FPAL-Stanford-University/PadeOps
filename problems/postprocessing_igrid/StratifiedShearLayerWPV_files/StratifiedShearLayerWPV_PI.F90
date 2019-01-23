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
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, buff3
   real(rkind), dimension(:,:,:), allocatable :: buff4, buff5, buff6
   real(rkind), dimension(:,:,:), allocatable :: buff7, buff8, buff9
   real(rkind), dimension(:,:,:), allocatable :: ufluct, vfluct, w, T
   real(rkind) :: dx, dy, dz, Re = 3000.d0, Rib = 0.05d0, Tref = 100.d0
   real(rkind) :: err
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0
   integer :: idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
   logical :: isZPeriodic = .false.
   integer :: VizDump_Schedule = 0
   logical :: computeStokesPressure = .true.
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

   call ops%allocate3DField(buff1)
   call ops%allocate3DField(buff2)
   call ops%allocate3DField(buff3)
   call ops%allocate3DField(buff4)
   call ops%allocate3DField(buff5)
   call ops%allocate3DField(buff6)
   call ops%allocate3DField(buff7)
   call ops%allocate3DField(buff8)
   call ops%allocate3DField(buff9)
   call ops%allocate3DField(ufluct)
   call ops%allocate3DField(vfluct)
   call ops%allocate3DField(w)
   call ops%allocate3DField(T)

   
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
      call ops%ReadField3D(buff4,"uVPi",TIDX)
      call ops%ReadField3D(buff5,"vVPi",TIDX)
      call ops%ReadField3D(buff6,"wVPi",TIDX)
      call message(0, "Read simulation data at time:", times(idx))

      ! STEP 0: Compute Fluctuations

      call ops%getFluct_from_MeanZ(buff1,ufluct)
      call ops%getFluct_from_MeanZ(buff2,vfluct)

      
      ! Check Om_Pi dot Om_Xi is Point-wise 0

      call ops%getCurl(buff4,buff5,buff6, buff1,buff2,buff3,1,1,1,1)
      buff7 = ufluct - buff4
      buff8 = vfluct - buff5
      buff9 = w - buff6
      call ops%getCurl(buff7,buff8,buff9, buff4,buff5,buff6,1,1,1,1)
      call ops%WriteField3D(buff1, "vXPi", tidx)
      call ops%WriteField3D(buff2, "vYPi", tidx)
      call ops%WriteField3D(buff3, "vZPi", tidx)
      call ops%WriteField3D(buff4, "vXXi", tidx)
      call ops%WriteField3D(buff5, "vYXi", tidx)
      call ops%WriteField3D(buff6, "vZXi", tidx)

      !buff7 = buff1*buff4 + buff2*buff5 + buff3*buff6
      !err = buff7(512,512,512)!maxval(abs(buff7))
      !call message(0, "Error in orthogonality:", buff7(512,512,512))





      ! Check Om_Pi here is same as decomposition Om_Pi



      call ops%getCurl(ufluct,vfluct,w, buff1,buff2,buff3,1,1,1,1)
      call ops%getGradient(T,           buff4,buff5,buff6,1,1)
      buff7 = buff1 * buff4 + buff2 * buff5 + buff3 * buff6 ! Pi
      buff8 = buff4*buff4 + buff5*buff5 + buff6*buff6
      buff8 = sqrt(buff8) ! |gradT|
      buff7 = buff7/buff8 ! |omega_Pi|
      
      buff1 = buff7*buff4/buff8
      buff2 = buff7*buff5/buff8
      buff3 = buff7*buff6/buff8

      call ops%WriteField3D(buff1, "oXPi", tidx)
      call ops%WriteField3D(buff2, "oYPi", tidx)
      call ops%WriteField3D(buff3, "oZPi", tidx)


      !call ops%ReadField3D(buff4,"uVPi",TIDX)
      !call ops%ReadField3D(buff5,"vVPi",TIDX)
      !call ops%ReadField3D(buff6,"wVPi",TIDX)
      !call ops%getCurl(buff4,buff5,buff6, buff7,buff8,buff9,1,1,1,1)


      !buff1 = abs(buff7 - buff1) + abs(buff8 - buff2) + abs(buff9 - buff3)
      !err = buff1(512,512,512)!maxval(abs(buff1))
      !call message(0, "Error in comparing Omega_Pi:", buff1(512,512,512))






      call ops%getCurl(ufluct,vfluct,w, buff4,buff5,buff6,1,1,1,1)

      buff1 = buff4 - buff1
      buff2 = buff5 - buff2
      buff3 = buff6 - buff3 ! omega decomp

      call ops%WriteField3D(buff1, "oXXi", tidx)
      call ops%WriteField3D(buff2, "oYXi", tidx)
      call ops%WriteField3D(buff3, "oZXi", tidx)

      !call ops%ReadField3D(buff4,"uVPi",TIDX)
      !call ops%ReadField3D(buff5,"vVPi",TIDX)
      !call ops%ReadField3D(buff6,"wVPi",TIDX)

      !buff7 = ufluct - buff4
      !buff8 = vfluct - buff5
      !buff9 = w - buff6
      !call ops%getCurl(buff7,buff8,buff9, buff4,buff5,buff6,1,1,1,1)

      !buff1 = abs(buff4 - buff1) + abs(buff5 - buff2) + abs(buff6 - buff3)
      !err = buff1(512,512,512)!maxval(abs(buff1))
      !call message(0, "Error in comparing Omega_Xi:", buff1(512,512,512))

      idx = idx + 1
      call toc()
   end do
   
   
   call MPI_Finalize(ierr)


end program
