program StratifiedShearLayerDomainIntegrals
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
   real(rkind), dimension(:,:,:), allocatable :: u, v, w, ufluct, vfluct, T,  Tfluct
   real(rkind) :: dx, dy, dz, Rib = 0.05d0, Pr = 1.d0, Tref = 100.d0 
   integer :: nx, ny, nz, RunID, TIDX
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0
   integer :: idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 1
   logical :: isZPeriodic = .false. 
   real(rkind), dimension(:,:), allocatable :: Nsq, Ssq, R11, R12, R13, R22, R23, R33, TT, G11, G12, G13, G22, G23, G33 
   real(rkind), dimension(:,:), allocatable :: v11, v12, v13, v22, v23, v33, d11, d12, d13, d22, d23, d33, umn_set, Tmn_set, wT, uT, vT
   real(rkind), dimension(:,:), allocatable :: times
   integer :: nt 

   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep, nx, ny, nz, Rib, NumericalSchemeVert 
   
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
   

   nt = (tstop - tstart)/tstep + 1

   call message(0,"Number of snapshots to read:", nt)

   !if (nrank == 0) then
      allocate(umn_set(nz,nt))  
      allocate(Tmn_set(nz,nt))  
      allocate(Nsq(nz,nt))
      allocate(Ssq(nz,nt))
      allocate(R11(nz,nt))
      allocate(R12(nz,nt))
      allocate(R13(nz,nt))
      allocate(R22(nz,nt))
      allocate(R23(nz,nt))
      allocate(R33(nz,nt))
      allocate(G11(nz,nt))
      allocate(G12(nz,nt))
      allocate(G13(nz,nt))
      allocate(G22(nz,nt))
      allocate(G23(nz,nt))
      allocate(G33(nz,nt))
      allocate(d11(nz,nt))
      allocate(d12(nz,nt))
      allocate(d13(nz,nt))
      allocate(d22(nz,nt))
      allocate(d23(nz,nt))
      allocate(d33(nz,nt))
      allocate(v11(nz,nt))
      allocate(v12(nz,nt))
      allocate(v13(nz,nt))
      allocate(v22(nz,nt))
      allocate(v23(nz,nt))
      allocate(v33(nz,nt))
      allocate(TT(nz,nt))
      allocate(uT(nz,nt))
      allocate(vT(nz,nt))
      allocate(wT(nz,nt))

      allocate(times(nt,1))
   !end if 

   tidx = tstart
   idx = 1
   do while(tidx <= tstop)
      call message(0, "Reading fields for tid:", TIDX)
      call tic()
      call ops%ReadField3D(u,"uVel",TIDX)
      call ops%ReadField3D(v,"vVel",TIDX)
      call ops%ReadField3D(w,"wVel",TIDX)
      call ops%ReadField3D(T,"potT",TIDX)
      times(idx,1) = ops%getSimTime(tidx)
      call message(0, "Read simulation data at time:", times(idx,1))

      T = Rib*(T - Tref)  ! Rescale Potential temperature to buoyancy variable: b 
      
      ! STEP 0: Compute means, mean gradients and fluctuations
      call ops%TakeMean_xy(u,umn_set(:,idx))
      call ops%TakeMean_xy(T,Tmn_set(:,idx))
      
      call ops%ddz_1d(umn_set(:,idx),Ssq(:,idx))
      call ops%ddz_1d(Tmn_set(:,idx),Nsq(:,idx))

      call ops%getFluct_from_MeanZ(v,vfluct)
      call ops%getFluct_from_MeanZ(u,ufluct)
      call ops%getFluct_from_MeanZ(T,Tfluct)

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

      ! STEP 3: Compute velocity gradients
      call ops%GetGradient(ufluct, dudx, dudy, dudz, 1, 1)
      call ops%GetGradient(vfluct, dvdx, dvdy, dvdz, 1, 1)
      call ops%GetGradient(w     , dwdx, dwdy, dwdz, -1, -1)

      buff4 = dudx*dudx  + dudy*dudy + dudz*dudz
      call ops%TakeMean_xy(buff4,d11(:,idx))
      buff4 = dudx*dvdx  + dudy*dvdy + dudz*dvdz
      call ops%TakeMean_xy(buff4,d12(:,idx))
      buff4 = dudx*dwdx  + dudy*dwdy + dudz*dwdz
      call ops%TakeMean_xy(buff4,d13(:,idx))
      buff4 = dvdx*dvdx  + dvdy*dvdy + dvdz*dvdz
      call ops%TakeMean_xy(buff4,d22(:,idx))
      buff4 = dvdx*dwdx  + dvdy*dwdy + dvdz*dwdz
      call ops%TakeMean_xy(buff4,d23(:,idx))
      buff4 = dwdx*dwdx  + dwdy*dwdy + dwdz*dwdz
      call ops%TakeMean_xy(buff4,d33(:,idx))

      ! STEP 4: Compute the vorticity anisotropy
      buff1 = dwdy - dvdz
      buff2 = dudz - dwdx
      buff3 = dvdx - dudy

      buff4 = buff1*buff1
      call ops%TakeMean_xy(buff4,v11(:,idx))
      buff4 = buff1*buff2
      call ops%TakeMean_xy(buff4,v12(:,idx))
      buff4 = buff1*buff3
      call ops%TakeMean_xy(buff4,v13(:,idx))
      buff4 = buff2*buff2
      call ops%TakeMean_xy(buff4,v22(:,idx))
      buff4 = buff2*buff3
      call ops%TakeMean_xy(buff4,v23(:,idx))
      buff4 = buff3*buff3
      call ops%TakeMean_xy(buff4,v33(:,idx))

      tidx = tidx + tstep
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
      
      call ops%WriteASCII_2D(times, "time")
   end if 

   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

