program TimeAvgFields_Periodic
   use kind_parameters, only: rkind, clen
   use igrid_Operators_Periodic, only: Ops_Periodic
   use constants, only: pi, two
   use mpi
   use decomp_2d
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   !real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, buff3, buff4, buff5
   real(rkind), dimension(:,:,:), allocatable :: u, v, w, usum, vsum, wsum, uusum, vvsum, wwsum, uvsum, uwsum, vwsum !
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

   call ops%allocate3DField(uvsum)
   call ops%allocate3DField(uwsum)
   call ops%allocate3DField(vwsum)

   usum = 0.d0
   vsum = 0.d0
   wsum = 0.d0

   uusum = 0.d0
   vvsum = 0.d0
   wwsum = 0.d0

   uvsum = 0.d0
   uwsum = 0.d0
   vwsum = 0.d0

   idx = 0

   tidx = tstart1
   tstep = tstep1
   do while ( tidx <= tstop1) 
      call message("TID:", tidx)
      call tic()
      call ops%ReadField3D(u,"uVel",tidx, RunID)
      call ops%ReadField3D(v,"vVel",tidx, RunID)
      call ops%ReadField3D(w,"wVel",tidx, RunID)
      usum = usum + u
      vsum = vsum + v
      wsum = wsum + w
      
      uusum = uusum + u*u 
      vvsum = vvsum + v*v
      wwsum = wwsum + w*w 
   
      uvsum = uvsum + u*v 
      uwsum = uwsum + u*w
      vwsum = vwsum + v*w 

      idx = idx + 1
      tidx = tidx + tstep
      call toc()
   end do 

   tidx = tstart2
   tstep = tstep2
   do while ( tidx <= tstop2) 
      call message("TID:", tidx)
      call tic()
      call ops%ReadField3D(u,"uVel",tidx, RunID)
      call ops%ReadField3D(v,"vVel",tidx, RunID)
      call ops%ReadField3D(w,"wVel",tidx, RunID)
      usum = usum + u
      vsum = vsum + v
      wsum = wsum + w
      
      uusum = uusum + u*u 
      vvsum = vvsum + v*v
      wwsum = wwsum + w*w 
   
      uvsum = uvsum + u*v 
      uwsum = uwsum + u*w
      vwsum = vwsum + v*w 

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

   uvsum = uvsum/real(idx,rkind)
   uwsum = uwsum/real(idx,rkind)
   vwsum = vwsum/real(idx,rkind)


   call ops%WriteField3D(usum, "u_mn", tidx, RunID) 
   call ops%WriteField3D(vsum, "v_mn", tidx, RunID) 
   call ops%WriteField3D(wsum, "w_mn", tidx, RunID) 
   
   uusum = uusum - usum*usum
   vvsum = vvsum - vsum*vsum
   wwsum = wwsum - wsum*wsum

   uvsum = uvsum - usum*vsum
   uwsum = uwsum - usum*wsum
   vwsum = vwsum - vsum*wsum
   
   call ops%WriteField3D(uusum, "uumn", tidx, RunID) 
   call ops%WriteField3D(vvsum, "vvmn", tidx, RunID) 
   call ops%WriteField3D(wwsum, "wwmn", tidx, RunID) 
   
   call ops%WriteField3D(uvsum, "uvmn", tidx, RunID) 
   call ops%WriteField3D(uwsum, "uwmn", tidx, RunID) 
   call ops%WriteField3D(vwsum, "vwmn", tidx, RunID) 
   
   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

