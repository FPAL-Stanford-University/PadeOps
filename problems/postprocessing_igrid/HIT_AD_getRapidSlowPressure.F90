program HIT_AD_getRapidSlowPressure
   use kind_parameters, only: rkind, clen
   use igrid_Operators, only: igrid_ops
   use constants, only: pi, two
   use mpi
   use timer, only: tic, toc
   use exits, only: message, gracefulExit
   implicit none
   real(rkind), dimension(:,:,:), allocatable :: uM, vM, wM 
   real(rkind), dimension(:,:,:), allocatable :: uF, vF, wF
   real(rkind), dimension(:,:,:), allocatable :: u, v, w
   real(rkind), dimension(:,:,:), allocatable :: dudx, dudy, dudz 
   real(rkind), dimension(:,:,:), allocatable :: dvdx, dvdy, dvdz 
   real(rkind), dimension(:,:,:), allocatable :: dwdx, dwdy, dwdz 
   
   real(rkind), dimension(:,:,:), allocatable :: dudxM, dudyM, dudzM 
   real(rkind), dimension(:,:,:), allocatable :: dvdxM, dvdyM, dvdzM 
   real(rkind), dimension(:,:,:), allocatable :: dwdxM, dwdyM, dwdzM 
   
   real(rkind), dimension(:,:,:), allocatable :: R11, R12, R13, R22, R23, R33
   real(rkind), dimension(:,:,:), allocatable :: rbuff1, rbuff2
   real(rkind), dimension(:,:,:), allocatable :: pslow, prapid 
   
   real(rkind), dimension(:,:,:), allocatable :: psps, pspr, prpr
   real(rkind), dimension(:,:,:), allocatable :: T11r, T12r, T13r, T22r, T23r, T33r
   real(rkind), dimension(:,:,:), allocatable :: T11s, T12s, T13s, T22s, T23s, T33s

   real(rkind) :: dx, dy, dz
   integer :: nx, ny, nz, RunID, TIDX, tDumpRestart = 1000
   type(igrid_ops) :: ops
   character(len=clen) ::  inputdir, outputdir, RestartDir
   character(len=clen) :: inputfile
   real(rkind) :: Lx = 10.d0*pi, Ly = 2.d0*pi, Lz = 2.d0*pi
   integer :: idx, ierr, tstart, tstop, tstep, NumericalSchemeVert = 2, RestartIDX
   logical :: isZPeriodic = .true. 
   integer :: nt, MeanTIDX 
   logical :: file_found_u, file_found_v, file_found_w, useRestart = .false. 
   integer :: topBC = 0, botBC = 0

   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RestartDir, RunID, tstart, tstop, tstep, nx, ny, nz, MeanTIDX, useRestart, RestartIDX, tDumpRestart
   
   call MPI_Init(ierr)               
   call GETARG(1,inputfile)          
   open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=99, NML=INPUT)
   close(unit=99)

   dx = Lx/real(nx,rkind) 
   dy = Ly/real(ny,rkind) 
   dz = Lz/real(nz,rkind)

   ! Initialize the operator class
   call ops%init(nx, ny, nz, dx, dy, dz, InputDir, OutputDir, RunID, isZPeriodic, NUmericalSchemeVert, RestartDir)


   ! Allocate all the needed memory 
   call ops%allocate3DField(u)
   call ops%allocate3DField(v)
   call ops%allocate3DField(w)
   
   call ops%allocate3DField(uM)
   call ops%allocate3DField(vM)
   call ops%allocate3DField(wM)
   
   call ops%allocate3DField(uF)
   call ops%allocate3DField(vF)
   call ops%allocate3DField(wF)
   
   call ops%allocate3DField(R11)
   call ops%allocate3DField(R12)
   call ops%allocate3DField(R13)
   call ops%allocate3DField(R22)
   call ops%allocate3DField(R23)
   call ops%allocate3DField(R33)
  
   call ops%allocate3DField(dudx)
   call ops%allocate3DField(dudy)
   call ops%allocate3DField(dudz)
   call ops%allocate3DField(dvdx)
   call ops%allocate3DField(dvdy)
   call ops%allocate3DField(dvdz)
   call ops%allocate3DField(dwdx)
   call ops%allocate3DField(dwdy)
   call ops%allocate3DField(dwdz)

   call ops%allocate3DField(dudxM)
   call ops%allocate3DField(dudyM)
   call ops%allocate3DField(dudzM)
   call ops%allocate3DField(dvdxM)
   call ops%allocate3DField(dvdyM)
   call ops%allocate3DField(dvdzM)
   call ops%allocate3DField(dwdxM)
   call ops%allocate3DField(dwdyM)
   call ops%allocate3DField(dwdzM)

   call ops%allocate3DField(rbuff1)
   call ops%allocate3DField(rbuff2)
   call ops%allocate3DField(pslow)
   call ops%allocate3DField(prapid)
   
   call ops%allocate3DField(psps)
   call ops%allocate3DField(pspr)
   call ops%allocate3DField(prpr)
   
   call ops%allocate3DField(T11r) 
   call ops%allocate3DField(T12r)
   call ops%allocate3DField(T13r)
   call ops%allocate3DField(T22r)
   call ops%allocate3DField(T23r)
   call ops%allocate3DField(T33r)
   
   call ops%allocate3DField(T11s) 
   call ops%allocate3DField(T12s)
   call ops%allocate3DField(T13s)
   call ops%allocate3DField(T22s)
   call ops%allocate3DField(T23s)
   call ops%allocate3DField(T33s)
   
   call ops%ReadField3D(uM,"uMmn",MeanTIDX)
   call ops%ReadField3D(vM,"vMmn",MeanTIDX)
   call ops%ReadField3D(wM,"wMmn",MeanTIDX)

   call ops%ReadField3D(R11,"R11m",MeanTIDX)
   call ops%ReadField3D(R12,"R12m",MeanTIDX)
   call ops%ReadField3D(R13,"R13m",MeanTIDX)
   call ops%ReadField3D(R22,"R22m",MeanTIDX)
   call ops%ReadField3D(R23,"R23m",MeanTIDX)
   call ops%ReadField3D(R33,"R33m",MeanTIDX)

   ! Get mean flow gradients 
   call ops%GetGradient(uM, dudxM, dudyM, dudzM, botBC, topBC)
   call ops%GetGradient(vM, dvdxM, dvdyM, dvdzM, botBC, topBC)
   call ops%GetGradient(wM, dwdxM, dwdyM, dwdzM, botBC, topBC)
   
  
   call ops%initPoissonSolver(dx, dy, dz)

   if (useRestart) tstart = RestartIDX

   nt = (tstop - tstart)/tstep + 1
   call message(0,"Number of snapshots to read:", nt)

   tidx = tstart
   call message(0,"Now checking for existence of all files")
   do while(tidx <= tstop) 
     file_found_u = ops%check_dump_existence("uVel",TIDX)
     file_found_v = ops%check_dump_existence("uVel",TIDX)
     file_found_w = ops%check_dump_existence("uVel",TIDX)
      
     if (file_found_u .and. file_found_v .and. file_found_w) then
         !call message(0, "File succcessfully found for:", TIDX)
     else
         call message(0, "File missing for tid:", TIDX)
         call mpi_barrier(mpi_comm_world, ierr) 
         call gracefulExit("Some file was missing.", 324)
     end if 

     tidx = tidx + tstep
   end do 
   call message(1,"All data files exist.")
   if (mod(tDumpRestart,tstep) == 0) then
      call message(0,"tDumpRestart input is Legal. Can dump restart files.")
   else
      call message(0,"tDumpRestart input is illegal.")
      call gracefulExit("tDumpRestart input is illegal.",145)
   end if
  

   if (useRestart) then
      call ops%ReadSummingRestartInfo(RestartIDX, idx)

      call ops%ReadSummingRestart(T11r, "T11r", RestartIDX)    
      call ops%ReadSummingRestart(T22r, "T22r", RestartIDX)   
      call ops%ReadSummingRestart(T33r, "T33r", RestartIDX)   
      call ops%ReadSummingRestart(T12r, "T12r", RestartIDX)   
      call ops%ReadSummingRestart(T13r, "T13r", RestartIDX)   
      call ops%ReadSummingRestart(T23r, "T23r", RestartIDX)   

      call ops%ReadSummingRestart(T11s, "T11s", RestartIDX)  
      call ops%ReadSummingRestart(T22s, "T22s", RestartIDX)  
      call ops%ReadSummingRestart(T33s, "T33s", RestartIDX)  
      call ops%ReadSummingRestart(T12s, "T12s", RestartIDX)  
      call ops%ReadSummingRestart(T13s, "T13s", RestartIDX)  
      call ops%ReadSummingRestart(T23s, "T23s", RestartIDX)  
     
      call ops%ReadSummingRestart(prpr, "prpr", RestartIDX)  
      call ops%ReadSummingRestart(pspr, "pspr", RestartIDX)  
      call ops%ReadSummingRestart(psps, "psps", RestartIDX)  
      
      call message(0,"Restart files read successfully.")
      call message(1,"Number of visualizations summed in restart", idx)
      call message(1,"New starting index:", tstart)
   else
      T11r = 0.d0 
      T22r = 0.d0
      T33r = 0.d0
      T12r = 0.d0
      T13r = 0.d0
      T23r = 0.d0

      T11s = 0.d0
      T22s = 0.d0
      T33s = 0.d0
      T12s = 0.d0
      T13s = 0.d0
      T23s = 0.d0
     
      prpr = 0.d0 
      pspr = 0.d0 
      psps = 0.d0 

      idx = 0
   end if

   tidx = tstart
   do while(tidx <= tstop)
      call message(0, "Reading fields for tid:", TIDX)
      call tic()
      call ops%ReadField3D(u,"uVel",TIDX)
      call ops%ReadField3D(v,"vVel",TIDX)
      call ops%ReadField3D(w,"wVel",TIDX)
    
      uF = u - uM
      vF = v - vM
      wF = w - wM
   
      call ops%GetGradient(uF, dudx, dudy, dudz, botBC, topBC)
      call ops%GetGradient(vF, dvdx, dvdy, dvdz, botBC, topBC)
      call ops%GetGradient(wF, dwdx, dwdy, dwdz, botBC, topBC)

      prapid = -2.d0*(dudxM*dudx + dudyM*dvdx + dudzM*dwdx &
                 &  + dvdxM*dudy + dvdyM*dvdy + dvdzM*dwdy & 
                 &  + dwdxM*dudz + dwdyM*dvdz + dwdzM*dwdz )  

      call ops%PoissonSolvePeriodic_inplace(prapid,dealiasRHS=.true.)

      rbuff1 = R11 - uF*uF
      call ops%ddx(rbuff1,rbuff2)
      call ops%ddx(rbuff2,rbuff1)
      pslow = rbuff1

      rbuff1 = 2.d0*(R12 - uF*vF)
      call ops%ddx(rbuff1,rbuff2)
      call ops%ddy(rbuff2,rbuff1)
      pslow = pslow + rbuff1

      rbuff1 = 2.d0*(R13 - uF*wF)
      call ops%ddx(rbuff1,rbuff2)
      call ops%ddz(rbuff2,rbuff1, botBC, topBC)
      pslow = pslow + rbuff1

      rbuff1 = R22 - vF*vF
      call ops%ddy(rbuff1,rbuff2)
      call ops%ddy(rbuff2,rbuff1)
      pslow = pslow + rbuff1

      rbuff1 = 2.d0*(R23 - vF*wF)
      call ops%ddy(rbuff1,rbuff2)
      call ops%ddz(rbuff2,rbuff1, botBC, topBC)
      pslow = pslow + rbuff1

      rbuff1 = R33 - wF*wF
      call ops%ddz(rbuff1,rbuff2, botBC, topBC)
      call ops%ddz(rbuff2,rbuff1, botBC, topBC)
      pslow = pslow + rbuff1

      call ops%PoissonSolvePeriodic_inplace(pslow,dealiasRHS=.true.)

      T11r = T11r + prapid*(dudx + dudx) 
      T22r = T22r + prapid*(dvdy + dvdy)
      T33r = T33r + prapid*(dwdz + dwdz)
      T12r = T12r + prapid*(dudy + dvdx)
      T13r = T13r + prapid*(dudz + dwdx)
      T23r = T23r + prapid*(dvdz + dwdy)

      T11s = T11s + pslow*(dudx + dudx) 
      T22s = T22s + pslow*(dvdy + dvdy)
      T33s = T33s + pslow*(dwdz + dwdz)
      T12s = T12s + pslow*(dudy + dvdx)
      T13s = T13s + pslow*(dudz + dwdx)
      T23s = T23s + pslow*(dvdz + dwdy)
     
      prpr = prpr + prapid*prapid
      pspr = pspr + pslow*prapid
      psps = psps + pslow*pslow
   
      tidx = tidx + tstep
      idx = idx + 1
    
      if (mod(tidx,tDumpRestart) == 0) then
        call ops%WriteSummingRestartInfo(tidx,idx) 
        
        call ops%WriteSummingRestart(T11r, "T11r", tidx)    
        call ops%WriteSummingRestart(T22r, "T22r", tidx)   
        call ops%WriteSummingRestart(T33r, "T33r", tidx)   
        call ops%WriteSummingRestart(T12r, "T12r", tidx)   
        call ops%WriteSummingRestart(T13r, "T13r", tidx)   
        call ops%WriteSummingRestart(T23r, "T23r", tidx)   

        call ops%WriteSummingRestart(T11s, "T11s", tidx)  
        call ops%WriteSummingRestart(T22s, "T22s", tidx)  
        call ops%WriteSummingRestart(T33s, "T33s", tidx)  
        call ops%WriteSummingRestart(T12s, "T12s", tidx)  
        call ops%WriteSummingRestart(T13s, "T13s", tidx)  
        call ops%WriteSummingRestart(T23s, "T23s", tidx)  
     
        call ops%WriteSummingRestart(prpr, "prpr", tidx)  
        call ops%WriteSummingRestart(pspr, "pspr", tidx)  
        call ops%WriteSummingRestart(psps, "psps", tidx) 
        call message(0,"Summing restart files dumped.")
      end if 

      call toc()
   end do 

   T11r = T11r/real(idx,rkind)
   T22r = T22r/real(idx,rkind)
   T33r = T33r/real(idx,rkind)
   T12r = T12r/real(idx,rkind)
   T13r = T13r/real(idx,rkind)
   T23r = T23r/real(idx,rkind)


   T11s = T11s/real(idx,rkind)
   T22s = T22s/real(idx,rkind)
   T33s = T33s/real(idx,rkind)
   T12s = T12s/real(idx,rkind)
   T13s = T13s/real(idx,rkind)
   T23s = T23s/real(idx,rkind)

   prpr = prpr/real(idx,rkind)
   pspr = pspr/real(idx,rkind)
   psps = psps/real(idx,rkind)

   call ops%writeField3D(T11r,"T11r",tidx)
   call ops%writeField3D(T22r,"T22r",tidx)
   call ops%writeField3D(T33r,"T33r",tidx)
   call ops%writeField3D(T12r,"T12r",tidx)
   call ops%writeField3D(T13r,"T13r",tidx)
   call ops%writeField3D(T23r,"T23r",tidx)
   
   call ops%writeField3D(T11s,"T11s",tidx)
   call ops%writeField3D(T22s,"T22s",tidx)
   call ops%writeField3D(T33s,"T33s",tidx)
   call ops%writeField3D(T12s,"T12s",tidx)
   call ops%writeField3D(T13s,"T13s",tidx)
   call ops%writeField3D(T23s,"T23s",tidx)
  
   call ops%writeField3D(prpr,"prpr",tidx)
   call ops%writeField3D(pspr,"pspr",tidx)
   call ops%writeField3D(psps,"psps",tidx)

   call ops%destroy()
   call MPI_Finalize(ierr)           

end program 

