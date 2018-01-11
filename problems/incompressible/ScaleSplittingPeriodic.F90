program ScaleSplittingPeriodic
   use kind_parameters, only: rkind, clen
   use igrid_Operators_Periodic, only: Ops_Periodic
   use constants, only: pi, two
   use mpi
   use decomp_2d
   use timer, only: tic, toc
   use exits, only: message
   use spectralMod, only: spectral
   use fof_mod, only: fof
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: u, v, w, u1, v1, w1, u2, v2, w2, u3, v3, w3
   real(rkind), dimension(:,:,:), allocatable :: R11_1, R12_1, R13_1, R22_1, R23_1, R33_1
   real(rkind), dimension(:,:,:), allocatable :: R11_2, R12_2, R13_2, R22_2, R23_2, R33_2
   real(rkind), dimension(:,:,:), allocatable :: R11_3, R12_3, R13_3, R22_3, R23_3, R33_3
   real(rkind), dimension(:,:,:), allocatable :: umean, vmean, wmean
   complex(rkind), dimension(:,:,:), allocatable :: uhat, vhat, what, cbuffy, cbuffz
   real(rkind) :: Lx = 10.d0*pi, Ly = 2.d0*pi, Lz = 2.d0*pi
   real(rkind) :: dx, dy, dz
   integer :: nx, ny, nz, RunID, tidx, tstart1, tstop1, tstep1, tstart2, tstop2,tstep2
   type(decomp_info) :: gp
   type(Ops_Periodic) :: ops
   logical :: periodicbcs(3)
   integer :: ierr, tstep, idx, tidx_mean 
   type(spectral), pointer :: spect

   type(fof) :: filt1, filt2
   character(len=clen) ::  inputdir, outputdir, MeanDir
   character(len=clen) :: inputfile, Fof_dir

   namelist /INPUT/ InputDir, OutputDir, MeanDir, fof_dir,  RunID, tidx_mean, nx, ny, nz, Lx, Ly, Lz,tstart1, tstart2, tstep1, tstep2, tstop1, tstop2
    
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
   call ops%link_spect(spect)

   call ops%allocate3DField(u)
   call ops%allocate3DField(v)
   call ops%allocate3DField(w)
   call ops%allocate3DField(u1)
   call ops%allocate3DField(v1)
   call ops%allocate3DField(w1)
   call ops%allocate3DField(u2)
   call ops%allocate3DField(v2)
   call ops%allocate3DField(w2)
   call ops%allocate3DField(u3)
   call ops%allocate3DField(v3)
   call ops%allocate3DField(w3)
   
   call ops%allocate3DField(umean)
   call ops%allocate3DField(vmean)
   call ops%allocate3DField(wmean)

   call ops%allocate3DField(R11_1)
   call ops%allocate3DField(R12_1)
   call ops%allocate3DField(R13_1)
   call ops%allocate3DField(R22_1)
   call ops%allocate3DField(R23_1)
   call ops%allocate3DField(R33_1)
   
   call ops%allocate3DField(R11_2)
   call ops%allocate3DField(R12_2)
   call ops%allocate3DField(R13_2)
   call ops%allocate3DField(R22_2)
   call ops%allocate3DField(R23_2)
   call ops%allocate3DField(R33_2)
   
   call ops%allocate3DField(R11_3)
   call ops%allocate3DField(R12_3)
   call ops%allocate3DField(R13_3)
   call ops%allocate3DField(R22_3)
   call ops%allocate3DField(R23_3)
   call ops%allocate3DField(R33_3)
   
   call ops%ReadField3D(umean,"u_mn",tidx_mean, RunID, meandir)
   call ops%ReadField3D(vmean,"v_mn",tidx_mean, RunID, meandir)
   call ops%ReadField3D(wmean,"w_mn",tidx_mean, RunID, meandir)

   allocate(cbuffy(spect%spectdecomp%ysz(1),spect%spectdecomp%ysz(2),spect%spectdecomp%ysz(3)))
   allocate(cbuffz(spect%spectdecomp%zsz(1),spect%spectdecomp%zsz(2),spect%spectdecomp%zsz(3)))
   
   call filt1%init(runID, fof_dir, meandir, 1, spect, cbuffy, cbuffz, .true., gp)
   call filt2%init(runID, fof_dir, meandir, 2, spect, cbuffy, cbuffz, .true., gp)

   allocate(uhat(spect%spectdecomp%ysz(1),spect%spectdecomp%ysz(2),spect%spectdecomp%ysz(3)))
   allocate(vhat(spect%spectdecomp%ysz(1),spect%spectdecomp%ysz(2),spect%spectdecomp%ysz(3)))
   allocate(what(spect%spectdecomp%ysz(1),spect%spectdecomp%ysz(2),spect%spectdecomp%ysz(3)))
   
   R11_1 = 0.d0
   R12_1 = 0.d0
   R13_1 = 0.d0
   R22_1 = 0.d0
   R23_1 = 0.d0
   R33_1 = 0.d0

   R11_2 = 0.d0
   R12_2 = 0.d0
   R13_2 = 0.d0
   R22_2 = 0.d0
   R23_2 = 0.d0
   R33_2 = 0.d0

   R11_3 = 0.d0
   R12_3 = 0.d0
   R13_3 = 0.d0
   R22_3 = 0.d0
   R23_3 = 0.d0
   R33_3 = 0.d0
   
   idx = 0
   tidx = tstart1
   tstep = tstep1
   do while ( tidx <= tstop1) 
      call message("TID:", tidx)
      call tic()
      call ops%ReadField3D(u,"uVel",tidx, RunID)
      call ops%ReadField3D(v,"vVel",tidx, RunID)
      call ops%ReadField3D(w,"wVel",tidx, RunID)

      call toc()
      call tic()

      u = u - umean
      v = v - vmean
      w = w - wmean

      call spect%fft(u,uhat)
      call spect%fft(v,vhat)
      call spect%fft(w,what)

      ! 1st filter 
      call filt1%filter_Complex2Real(uhat, u1)
      call filt1%filter_Complex2Real(vhat, v1)
      call filt1%filter_Complex2Real(what, w1)


      ! 2nd filter 
      call filt2%filter_Complex2Real(uhat, u2)
      call filt2%filter_Complex2Real(vhat, v2)
      call filt2%filter_Complex2Real(what, w2)
  
      u3 = u - u2
      v3 = v - v2
      w3 = w - w2

      u2 = u2 - u1
      v2 = v2 - v1
      w2 = w2 - w1


      call ops%WriteField3D(u, "uflt", 0, RunID) 
      call ops%WriteField3D(u1, "usc1", 0, RunID) 
      call ops%WriteField3D(u2, "usc2", 0, RunID) 
      call ops%WriteField3D(u3, "usc3", 0, RunID) 
      call mpi_barrier(mpi_comm_world, ierr)

      print*, "Done"
      stop 

      R11_1 = R11_1 + u1*u1 
      R12_1 = R12_1 + u1*v1 
      R13_1 = R13_1 + u1*w1 
      R22_1 = R22_1 + v1*v1 
      R23_1 = R23_1 + v1*w1 
      R33_1 = R33_1 + w1*w1 

      R11_2 = R11_2 + u2*u2 
      R12_2 = R12_2 + u2*v2 
      R13_2 = R13_2 + u2*w2 
      R22_2 = R22_2 + v2*v2 
      R23_2 = R23_2 + v2*w2 
      R33_2 = R33_2 + w2*w2 

      R11_3 = R11_3 + u3*u3 
      R12_3 = R12_3 + u3*v3 
      R13_3 = R13_3 + u3*w3 
      R22_3 = R22_3 + v3*v3 
      R23_3 = R23_3 + v3*w3 
      R33_3 = R33_3 + w3*w3 

      

      !call ops%WriteField3D(u1,"uSc1",tidx, RunID)
      !call ops%WriteField3D(v1,"vSc1",tidx, RunID)
      !call ops%WriteField3D(w1,"wSc1",tidx, RunID)

      !call ops%WriteField3D(u2,"uSc2",tidx, RunID)
      !call ops%WriteField3D(v2,"vSc2",tidx, RunID)
      !call ops%WriteField3D(w2,"wSc2",tidx, RunID)

      !call ops%WriteField3D(u3,"uSc3",tidx, RunID)
      !call ops%WriteField3D(v3,"vSc3",tidx, RunID)
      !call ops%WriteField3D(w3,"wSc3",tidx, RunID)

      !call mpi_barrier(mpi_comm_world, ierr)
      !stop 

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

      call toc()
      call tic()

      u = u - umean
      v = v - vmean
      w = w - wmean

      call spect%fft(u,uhat)
      call spect%fft(v,vhat)
      call spect%fft(w,what)

      ! 1st filter 
      call filt1%filter_Complex2Real(uhat, u1)
      call filt1%filter_Complex2Real(vhat, v1)
      call filt1%filter_Complex2Real(what, w1)


      ! 2nd filter 
      call filt2%filter_Complex2Real(uhat, u2)
      call filt2%filter_Complex2Real(vhat, v2)
      call filt2%filter_Complex2Real(what, w2)
  
      u3 = u - u2
      v3 = v - v2
      w3 = w - w2

      u2 = u2 - u1
      v2 = v2 - v1
      w2 = w2 - w1

      R11_1 = R11_1 + u1*u1 
      R12_1 = R12_1 + u1*v1 
      R13_1 = R13_1 + u1*w1 
      R22_1 = R22_1 + v1*v1 
      R23_1 = R23_1 + v1*w1 
      R33_1 = R33_1 + w1*w1 

      R11_2 = R11_2 + u2*u2 
      R12_2 = R12_2 + u2*v2 
      R13_2 = R13_2 + u2*w2 
      R22_2 = R22_2 + v2*v2 
      R23_2 = R23_2 + v2*w2 
      R33_2 = R33_2 + w2*w2 

      R11_3 = R11_3 + u3*u3 
      R12_3 = R12_3 + u3*v3 
      R13_3 = R13_3 + u3*w3 
      R22_3 = R22_3 + v3*v3 
      R23_3 = R23_3 + v3*w3 
      R33_3 = R33_3 + w3*w3 

      idx = idx + 1
      tidx = tidx + tstep
      call toc()
   end do 

   R11_1 = R11_1/real(idx) 
   R12_1 = R12_1/real(idx) 
   R13_1 = R13_1/real(idx) 
   R22_1 = R22_1/real(idx) 
   R23_1 = R23_1/real(idx) 
   R33_1 = R33_1/real(idx) 
   
   R11_2 = R11_2/real(idx) 
   R12_2 = R12_2/real(idx) 
   R13_2 = R13_2/real(idx) 
   R22_2 = R22_2/real(idx) 
   R23_2 = R23_2/real(idx) 
   R33_2 = R33_2/real(idx) 
   
   R11_3 = R11_3/real(idx) 
   R12_3 = R12_3/real(idx) 
   R13_3 = R13_3/real(idx) 
   R22_3 = R22_3/real(idx) 
   R23_3 = R23_3/real(idx) 
   R33_3 = R33_3/real(idx) 
   
   call ops%WriteField3D(R11_1, "R11a", tidx, RunID) 
   call ops%WriteField3D(R12_1, "R12a", tidx, RunID) 
   call ops%WriteField3D(R13_1, "R13a", tidx, RunID) 
   call ops%WriteField3D(R22_1, "R22a", tidx, RunID) 
   call ops%WriteField3D(R23_1, "R23a", tidx, RunID) 
   call ops%WriteField3D(R33_1, "R33a", tidx, RunID) 
   
   call ops%WriteField3D(R11_2, "R11b", tidx, RunID) 
   call ops%WriteField3D(R12_2, "R12b", tidx, RunID) 
   call ops%WriteField3D(R13_2, "R13b", tidx, RunID) 
   call ops%WriteField3D(R22_2, "R22b", tidx, RunID) 
   call ops%WriteField3D(R23_2, "R23b", tidx, RunID) 
   call ops%WriteField3D(R33_2, "R33b", tidx, RunID) 
   
   call ops%WriteField3D(R11_3, "R11c", tidx, RunID) 
   call ops%WriteField3D(R12_3, "R12c", tidx, RunID) 
   call ops%WriteField3D(R13_3, "R13c", tidx, RunID) 
   call ops%WriteField3D(R22_3, "R22c", tidx, RunID) 
   call ops%WriteField3D(R23_3, "R23c", tidx, RunID) 
   call ops%WriteField3D(R33_3, "R33c", tidx, RunID) 
   
   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

