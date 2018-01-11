program TimeAvgVort_Periodic
   use kind_parameters, only: rkind, clen
   use igrid_Operators_Periodic, only: Ops_Periodic
   use constants, only: pi, two
   use mpi
   use decomp_2d
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   !real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, buff3, buff4, buff5
   real(rkind), dimension(:,:,:), allocatable :: u, v, w, usum, vsum, wsum, uusum, vvsum, wwsum
   !real(rkind), dimension(:),     allocatable :: Prod, Transp_Conv, Transp_Press, Transp_Visc, Dissp, DisspSGS, DisspTheta, DisspThetaSGS, Buoy

   real(rkind) :: Lx = 10.d0*pi, Ly = 2.d0*pi, Lz = 2.d0*pi
   real(rkind) :: dx, dy, dz
   integer :: nx, ny, nz, RunID, tidx, tstart1, tstop1, tstep1, tstart2, tstop2,tstep2
   type(decomp_info) :: gp
   type(Ops_Periodic) :: ops
   logical :: periodicbcs(3)
   integer :: ierr, tstep, idx 

   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile

   namelist /INPUT/ InputDir, OutputDir, RunID, nx, ny, nz, Lx, Ly, Lz,tstart1, tstart2, tstep1, tstep2, tstop1, tstop2
    
   call MPI_Init(ierr)               
   call GETARG(1,inputfile)          
   open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=99, NML=INPUT)
   close(unit=99)
   periodicbcs(1) = .true.; periodicbcs(2) = .true.; periodicbcs(3) = .true.  
   call decomp_2d_init(nx, ny, nz, 0, 0, periodicbcs)
   call get_decomp_info(gp)

   dx = Lx/real(nx,rkind) 
   dy = Ly/real(ny,rkind) 
   dz = Lz/real(nz,rkind)

   call ops%init(nx, ny, nz, dx, dy, dz, gp, InputDir, OutputDir)

   call ops%allocate3DField(u)
   call ops%allocate3DField(v)
   call ops%allocate3DField(w)
   
   call ops%allocate3DField(usum)
   call ops%allocate3DField(vsum)
   call ops%allocate3DField(wsum)
   
   call ops%allocate3DField(uusum)
   call ops%allocate3DField(vvsum)
   call ops%allocate3DField(wwsum)


   usum = 0.d0
   vsum = 0.d0
   wsum = 0.d0

   uusum = 0.d0
   vvsum = 0.d0
   wwsum = 0.d0

   idx = 0

   tidx = tstart1
   tstep = tstep1
   do while ( tidx <= tstop1) 
      call message("TID:", tidx)
      call tic()
      call ops%ReadField3D(u,"omgX",tidx, RunID)
      call ops%ReadField3D(v,"omgY",tidx, RunID)
      call ops%ReadField3D(w,"omgZ",tidx, RunID)
      usum = usum + u
      vsum = vsum + v
      wsum = wsum + w
      
      uusum = uusum + u*u 
      vvsum = vvsum + v*v
      wwsum = wwsum + w*w 

      idx = idx + 1
      tidx = tidx + tstep
      call toc()
   end do 

   tidx = tstart2
   tstep = tstep2
   do while ( tidx <= tstop2) 
      call message("TID:", tidx)
      call tic()
      call ops%ReadField3D(u,"omgX",tidx, RunID)
      call ops%ReadField3D(v,"omgY",tidx, RunID)
      call ops%ReadField3D(w,"omgZ",tidx, RunID)
      usum = usum + u
      vsum = vsum + v
      wsum = wsum + w
      
      uusum = uusum + u*u 
      vvsum = vvsum + v*v
      wwsum = wwsum + w*w 
   
      idx = idx + 1
      tidx = tidx + tstep
      call toc()
   end do 

   usum = usum/real(idx,rkind)
   vsum = vsum/real(idx,rkind)
   wsum = wsum/real(idx,rkind)

   uusum = uusum/real(idx,rkind)
   vvsum = vvsum/real(idx,rkind)
   wwsum = wwsum/real(idx,rkind)

   call ops%WriteField3D(usum, "oxmn", tidx, RunID) 
   call ops%WriteField3D(vsum, "oymn", tidx, RunID) 
   call ops%WriteField3D(wsum, "ozmn", tidx, RunID) 
   
   uusum = uusum - usum*usum
   vvsum = vvsum - vsum*vsum
   wwsum = wwsum - wsum*wsum

   call ops%WriteField3D(uusum, "oxox", tidx, RunID) 
   call ops%WriteField3D(vvsum, "oyoy", tidx, RunID) 
   call ops%WriteField3D(wwsum, "ozoz", tidx, RunID) 
   
   
   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

