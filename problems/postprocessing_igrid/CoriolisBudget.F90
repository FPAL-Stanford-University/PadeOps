program CoriolisBudget
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use decomp_2d, only: nrank 
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, buff3, buff4
   real(rkind), dimension(:,:,:), allocatable :: dudx, dudy, dudz
   real(rkind), dimension(:,:,:), allocatable :: dvdx, dvdy, dvdz
   real(rkind), dimension(:,:,:), allocatable :: dwdx, dwdy, dwdz
   real(rkind), dimension(:,:,:), allocatable :: u, v, w, ufluct, vfluct, T,  Tfluct, P, kSGS, nSGS, Pfluct
   real(rkind) :: dx, dy, dz, Rib = 0.05d0, Pr = 1.d0, Tref = 100.d0 
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0
   integer :: idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
   logical :: isZPeriodic = .false. 
   real(rkind), dimension(:,:), allocatable :: R11, R12, R13, R22, R23, R33, TT
   real(rkind), dimension(:,:), allocatable :: d11, d12, d13, d22, d23, d33
   real(rkind), dimension(:,:), allocatable :: DD11, DD12, DD13, DD22, DD23, DD33
   real(rkind), dimension(:,:), allocatable :: umn_set, Tmn_set, wT, uT, vT
   integer :: VizDump_Schedule = 0
   integer, dimension(:), allocatable :: timesteps
   real(rkind), dimension(:), allocatable :: times, buff1d_1, buff1d_2, buff1d_3
   real(rkind), dimension(:,:), allocatable :: timewrite
   integer :: nt 

   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep, nx, ny, nz, NumericalSchemeVert, VizDump_Schedule 
   
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
   call ops%allocate3DField(wfluct)
   call ops%allocate3DField(w)
   call ops%allocate3DField(T)
   call ops%allocate3DField(Tfluct)
   call ops%allocate3DField(P)
   call ops%allocate3DField(kSGS)       
   call ops%allocate3DField(nSGS)
   call ops%allocate3Dfield(dudx)
   call ops%allocate3Dfield(dudy)
   call ops%allocate3Dfield(dudz)
   call ops%allocate3Dfield(dvdx)
   call ops%allocate3Dfield(dvdy)
   call ops%allocate3Dfield(dvdz)
   call ops%allocate3Dfield(dwdx)
   call ops%allocate3Dfield(dwdy)
   call ops%allocate3Dfield(dwdz)
   
   call ops%allocate3DField(buff1)
   call ops%allocate3DField(buff2)
   call ops%allocate3DField(buff3)
   call ops%allocate3DField(buff4)
   call ops%allocate1DField(buff1d_1)
   call ops%allocate1DField(buff1d_2)
   call ops%allocate1DField(buff1d_3) 
   call ops%allocate1DField(buff1d_4) 


   if (VizDump_Schedule == 1) then
      call ops%Read_VizSummary(times, timesteps)
      nt = size(timesteps)
   else
      nt = (tstop - tstart)/tstep
   end if

   

   call message(0,"Number of snapshots to read:", nt)

   !if (nrank == 0) then
      ! Means
      allocate(umn_set(nz,nt))  
      allocate(Tmn_set(nz,nt))  
      !allocate(Nsq(nz,nt))
      !allocate(Ssq(nz,nt))
      ! Reynolds stresses
      allocate(R11(nz,nt))
      allocate(R12(nz,nt))
      allocate(R13(nz,nt))
      allocate(R22(nz,nt))
      allocate(R23(nz,nt))
      allocate(R33(nz,nt))
      ! Velocity temperature corr
      allocate(uT(nz,nt))
      allocate(vT(nz,nt))
      allocate(wT(nz,nt))
      allocate(timewrite(nt,1))
      ! Pressure and vel grad corr
      allocate(dudzP(nz,nt))
      allocate(dwdxP(nz,nt))
      allocate(dvdzP(nz,nt))
      allocate(dwdyP(nz,nt))
      allocate(dwPdx(nz,nt))
      allocate(duPdz(nz,nt))
      allocate(dwPdy(nz,nt))
      allocate(dvPdz(nz,nt))
      ! Diffusion
      allocate(d13(nz,nt))
      allocate(d23(nz,nt))
      ! Dissipation
      allocate(DD13(nz,nt))
      allocate(DD23(nz,nt))
      ! Triple corr
      allocate(triple13(nz,nt))
      allocate(triple23(nz,nt))
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
      call ops%ReadField3D(kSGS,'kSGS',TIDX)
      call ops%ReadField3D(nSGS,'nSGS',TIDX)
      call ops%ReadField3D(prss,'P', TIDX)
      call message(0, "Read simulation data at time:", times(idx))

      !T = Rib*(T - Tref)  ! Rescale Potential temperature to buoyancy variable: b 
      
      ! STEP 0: Compute means, mean gradients and fluctuations
      call ops%TakeMean_xy(u,umn_set(:,idx))
      call ops%TakeMean_xy(T,Tmn_set(:,idx))
      
      !call ops%ddz_1d(umn_set(:,idx),Ssq(:,idx))
      !call ops%ddz_1d(Tmn_set(:,idx),Nsq(:,idx))

      call ops%getFluct_from_MeanZ(v,vfluct)
      call ops%getFluct_from_MeanZ(u,ufluct)
      call ops%getFluct_from_MeanZ(w,wfluct)
      call ops%getFluct_from_MeanZ(T,Tfluct)
      call ops%getFluct_from_MeanZ(P,Pfluct)

      ! STEP 1: Compute correlations
      buff2 = ufluct*ufluct
      call ops%TakeMean_xy(buff2,R11(:,idx))
      buff2 = ufluct*vfluct
      call ops%TakeMean_xy(buff2,R12(:,idx))
      buff2 = ufluct*w
      call ops%TakeMean_xy(buff2,R13(:,idx))
      buff2 = vfluct*vfluct
      call ops%TakeMean_xy(buff2,R22(:,idx))
      buff2 = vfluct*w
      call ops%TakeMean_xy(buff2,R23(:,idx))
      buff2 = w*w
      call ops%TakeMean_xy(buff2,R33(:,idx))
      buff2 = Tfluct*Tfluct
      call ops%TakeMean_xy(buff2,TT(:,idx))
      buff2 = ufluct*Tfluct
      call ops%TakeMean_xy(buff2,uT(:,idx))
      buff2 = vfluct*Tfluct
      call ops%TakeMean_xy(buff2,vT(:,idx))
      buff2 = w*Tfluct
      call ops%TakeMean_xy(buff2,wT(:,idx))

      ! STEP 2: Compute velocity gradients
      call ops%GetGradient(ufluct, dudx, dudy, dudz, -1, 0)
      call ops%GetGradient(vfluct, dvdx, dvdy, dvdz, -1, 0)
      call ops%GetGradient(wfluct, dwdx, dwdy, dwdz, -1, 0)

      ! STEP 3: Compute the pressure-gradient correlations
      buff4 = -Pfluct * wfluct
      call ops%TakeMean_xy(buff4, buff4)
      call ops%GetGradient(buff4, dwPdx(:,idx), dwPdy(:,idx), buff1d_3, 0, 0) 
      buff4 = -Pfluct * ufluct
      call ops%TakeMean_xy(buff4, buff4)
      call ops%ddz_1d(buff4, duPdz(:,idx), 0, 0) 
      buff4 = -Pfluct * vfluct
      call ops%TakeMean_xy(buff4, buff4)
      call ops%ddz_1d(buff4, dvPdz(:,idx), 0, 0)
       
      buff1 = Pfluct * dudz
      call ops%TakeMean_xy(buff1, dudzP(:,idx))
      buff1 = Pfluct * dwdx
      call ops%TakeMean_xy(buff1, dwdxP(:,idx))
      buff1 = Pfluct * dvdz
      call ops%TakeMean_xy(buff1, dvdzP(:,idx))
      buff1 = Pfluct * dwdy
      call ops%TakeMean_xy(buff1, dwdyP(:,idx))

      ! STEP 2: Compute Temperature gradient correlations 
      call ops%GetGradient(T, buff1, buff2, buff3, 1, 1)
      buff4 = buff1*buff1
      call ops%TakeMean_xy(buff4,G11(:,idx))
      buff4 = buff1*buff2
      call ops%TakeMean_xy(buff4,G12(:,idx))
      buff4 = buff1*buff3
      call ops%TakeMean_xy(buff4,G13(:,idx))
      buff4 = buff2*buff2
      call ops%TakeMean_xy(buff4,G22(:,idx))
      buff4 = buff2*buff3
      call ops%TakeMean_xy(buff4,G23(:,idx))
      buff4 = buff3*buff3
      call ops%TakeMean_xy(buff4,G33(:,idx))

      ! Step 3: Production
      call ops%TakeMean_xy(w, buff1d_1)
      call ops%GetGradient(buff1d_1, buff1d_2, buff1d_3, buff1d_4, 0, 0)
      P13 = -R11(:,idx) * buff1d_1 - R12(:,idx) * buff1d_2 - R13(:,idx) * buff1d_3 
      P23 = -R12(:,idx) * buff1d_1 - R22(:,idx) * buff1d_2 - R23(:,idx) * buff1d_3
      call ops%TakeMean_xy(u, buff1d_1)
      call ops%GetGradient(buff1d_1, buff1d_2, buff1d_3, buff1d_4, 0, 0)
      P31 = -R11(:,idx) * buff1d_1 - R12(:,idx) * buff1d_2 - R13(:,idx) * buff1d_3
      call ops%TakeMean_xy(v, buff1d_1)
      call ops%GetGradient(buff1d_1, buff1d_2, buff1d_3, buff1d_4, 0, 0)
      P32 = -R13(:,idx) * buff1d_1 - R23(:,idx) * buff1d_2 - R33(:,idx) * buff1d_3

      ! Dissipation
      buff4 = -2.*nSGS*(dudx*dudx  + dudy*dudy + dudz*dudz)
      call ops%TakeMean_xy(buff4,DD11(:,idx))
      buff4 = -2.*nSGS*(dudx*dvdx  + dudy*dvdy + dudz*dvdz)
      call ops%TakeMean_xy(buff4,DD12(:,idx))
      buff4 = -2.*nSGS*(dudx*dwdx  + dudy*dwdy + dudz*dwdz)
      call ops%TakeMean_xy(buff4,DD13(:,idx))
      buff4 = -2.*nSGS*(dvdx*dvdx  + dvdy*dvdy + dvdz*dvdz)
      call ops%TakeMean_xy(buff4,DD22(:,idx))
      buff4 = -2.*nSGS*(dvdx*dwdx  + dvdy*dwdy + dvdz*dwdz)
      call ops%TakeMean_xy(buff4,DD23(:,idx))
      buff4 = -2.*nSGS*(dwdx*dwdx  + dwdy*dwdy + dwdz*dwdz)
      call ops%TakeMean_xy(buff4,DD33(:,idx))

      ! Diffusion
      !call ops%GetGradient(R11, buff1d_1, buff1d_2, buff1d_3, 0, 0)
      !call ops%GetGradient(buff1d_1, 




      idx = idx + 1
      call toc()
   end do 

   if (nrank == 0) then
      call ops%WriteASCII_2D(umn_set, "umnZ")
      call ops%WriteASCII_2D(Tmn_set, "TmnZ")
      call ops%WriteASCII_2D(Ssq, "SsqZ")
      call ops%WriteASCII_2D(Nsq, "NsqZ")
      call ops%WriteASCII_2D(R11, "R11Z")
      call ops%WriteASCII_2D(R12, "R12Z")
      call ops%WriteASCII_2D(R13, "R13Z")
      call ops%WriteASCII_2D(R22, "R22Z")
      call ops%WriteASCII_2D(R23, "R23Z")
      call ops%WriteASCII_2D(R33, "R33Z")
      
      call ops%WriteASCII_2D(TT, "TTmZ")
      call ops%WriteASCII_2D(uT, "uTmZ")
      call ops%WriteASCII_2D(vT, "vTmZ")
      call ops%WriteASCII_2D(wT, "wTmZ")

      call ops%WriteASCII_2D(G11, "G11Z")
      call ops%WriteASCII_2D(G12, "G12Z")
      call ops%WriteASCII_2D(G13, "G13Z")
      call ops%WriteASCII_2D(G22, "G22Z")
      call ops%WriteASCII_2D(G23, "G23Z")
      call ops%WriteASCII_2D(G33, "G33Z")

      call ops%WriteASCII_2D(d11, "d11Z")
      call ops%WriteASCII_2D(d12, "d12Z")
      call ops%WriteASCII_2D(d13, "d13Z")
      call ops%WriteASCII_2D(d22, "d22Z")
      call ops%WriteASCII_2D(d23, "d23Z")
      call ops%WriteASCII_2D(d33, "d33Z")
      
      call ops%WriteASCII_2D(v11, "v11Z")
      call ops%WriteASCII_2D(v12, "v12Z")
      call ops%WriteASCII_2D(v13, "v13Z")
      call ops%WriteASCII_2D(v22, "v22Z")
      call ops%WriteASCII_2D(v23, "v23Z")
      call ops%WriteASCII_2D(v33, "v33Z")
      

      timewrite(:,1) = times
      call ops%WriteASCII_2D(timewrite, "time")
   end if 

   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

