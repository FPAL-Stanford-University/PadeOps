program TimeAvgFields_Periodic
   use kind_parameters, only: rkind, clen
   use igrid_Operators_Periodic, only: Ops_Periodic
   use constants, only: pi, two
   use mpi
   use decomp_2d
   use timer, only: tic, toc
   use exits, only: message
   use spectralMod, only: spectral
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: u, v, w, umean, vmean, wmean, buff1, buff2, buff3
   real(rkind), dimension(:,:,:), allocatable :: R11, R22, R33, frapid, fslow, prapid, pslow 
   real(rkind), dimension(:,:,:), allocatable :: rapidp_sum, slowp_sum, R11mn, R12mn, R13mn, R22mn, R23mn, R33mn
   
   real(rkind), dimension(:,:,:), allocatable, target :: R12, R13, R23
   real(rkind), dimension(:,:,:), pointer :: R21, R31, R32
   complex(rkind), dimension(:,:,:), allocatable :: cbuff1, cbuff2, cbuff3

   real(rkind), dimension(:,:,:,:,:), allocatable :: duidxj, duidxj_mean
   real(rkind) :: Lx = 10.d0*pi, Ly = 2.d0*pi, Lz = 2.d0*pi
   real(rkind) :: dx, dy, dz
   integer :: nx, ny, nz, RunID, tidx, tstart1, tstop1, tstep1, tstart2, tstop2,tstep2
   type(decomp_info) :: gp
   type(Ops_Periodic) :: ops
   type(spectral), pointer :: spect
   logical :: periodicbcs(3)
   integer :: ierr, tstep, idx, i, j, tidx_mean

   character(len=clen) ::  inputdir, outputdir, meandir
   character(len=clen) :: inputfile

   namelist /INPUT/ InputDir, OutputDir, meandir, RunID, tidx_mean, nx, ny, nz, Lx, Ly, Lz,tstart1, tstart2, tstep1, tstep2, tstop1, tstop2
    
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
   
   call ops%allocate3DField(umean)
   call ops%allocate3DField(vmean)
   call ops%allocate3DField(wmean)
   
   call ops%allocate3DField(R11)
   call ops%allocate3DField(R12)
   call ops%allocate3DField(R13)
   call ops%allocate3DField(R22)
   call ops%allocate3DField(R23)
   call ops%allocate3DField(R33)
   
   call ops%allocate3DField(R11mn)
   call ops%allocate3DField(R12mn)
   call ops%allocate3DField(R13mn)
   call ops%allocate3DField(R22mn)
   call ops%allocate3DField(R23mn)
   call ops%allocate3DField(R33mn)
   
   call ops%allocate3DField(buff1)
   call ops%allocate3DField(buff2)
   call ops%allocate3DField(buff3)
  
   call ops%allocate3DField(frapid)
   call ops%allocate3DField(fslow)
   
   call ops%allocate3DField(prapid)
   call ops%allocate3DField(pslow)
   
   
   R21 => R12
   R31 => R13
   R32 => R23

   allocate(duidxj(gp%xsz(1),gp%xsz(2),gp%xsz(3),3,3))
   allocate(duidxj_mean (gp%xsz(1),gp%xsz(2),gp%xsz(3),3,3))
   
   allocate(cbuff1(spect%spectdecomp%ysz(1),spect%spectdecomp%ysz(2),spect%spectdecomp%ysz(3)))
   allocate(cbuff2(spect%spectdecomp%ysz(1),spect%spectdecomp%ysz(2),spect%spectdecomp%ysz(3)))
   allocate(cbuff3(spect%spectdecomp%ysz(1),spect%spectdecomp%ysz(2),spect%spectdecomp%ysz(3)))


   allocate(rapidp_sum(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
   allocate(slowp_sum(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
   rapidp_sum = 0.d0
   slowp_sum = 0.d0

   call ops%ReadField3D(umean,"u_mn",tidx_mean, RunID, meandir)
   call ops%ReadField3D(vmean,"v_mn",tidx_mean, RunID, meandir)
   call ops%ReadField3D(wmean,"w_mn",tidx_mean, RunID, meandir)
   


   call ops%ddx(umean, duidxj_mean(:,:,:,1,1))
   call ops%ddx(vmean, duidxj_mean(:,:,:,2,1))
   call ops%ddx(wmean, duidxj_mean(:,:,:,3,1))

   call ops%ddy(umean, duidxj_mean(:,:,:,1,2))
   call ops%ddy(vmean, duidxj_mean(:,:,:,2,2))
   call ops%ddy(wmean, duidxj_mean(:,:,:,3,2))

   call ops%ddz(umean, duidxj_mean(:,:,:,1,3))
   call ops%ddz(vmean, duidxj_mean(:,:,:,2,3))
   call ops%ddz(wmean, duidxj_mean(:,:,:,3,3))
   
   call ops%ReadField3D(R11mn,"uumn",tidx_mean, RunID, meandir)
   call ops%ReadField3D(R22mn,"vvmn",tidx_mean, RunID, meandir)
   call ops%ReadField3D(R33mn,"wwmn",tidx_mean, RunID, meandir)
   call ops%ReadField3D(R12mn,"uvmn",tidx_mean, RunID, meandir)
   call ops%ReadField3D(R13mn,"uwmn",tidx_mean, RunID, meandir)
   call ops%ReadField3D(R23mn,"vwmn",tidx_mean, RunID, meandir)
 
   R11mn = 0.d0
   R12mn = 0.d0
   R13mn = 0.d0
   R22mn = 0.d0
   R23mn = 0.d0
   R33mn = 0.d0

   
   idx = 0
   tidx = tstart1
   tstep = tstep1
   do while ( tidx <= tstop1) 
      call message("TID:", tidx)
      call tic()
      call ops%ReadField3D(u,"uVel",tidx, RunID)
      call ops%ReadField3D(v,"vVel",tidx, RunID)
      call ops%ReadField3D(w,"wVel",tidx, RunID)

      

      u = u - umean
      v = v - vmean
      w = w - wmean

      !call ops%ddx(u, duidxj(:,:,:,1,1))
      !call ops%ddx(v, duidxj(:,:,:,2,1))
      !call ops%ddx(w, duidxj(:,:,:,3,1))
      
      call spect%fft(u, cbuff1)
      call spect%mtimes_ik1_oop(cbuff1,cbuff2)
      call spect%mtimes_ik2_oop(cbuff1,cbuff3)
      call ops%ddz_cmplx2cmplx(cbuff1)
      call spect%ifft(cbuff2,duidxj(:,:,:,1,1))
      call spect%ifft(cbuff3,duidxj(:,:,:,1,2))
      call spect%ifft(cbuff1,duidxj(:,:,:,1,3))

      !call ops%ddy(u, duidxj(:,:,:,1,2))
      !call ops%ddy(v, duidxj(:,:,:,2,2))
      !call ops%ddy(w, duidxj(:,:,:,3,2))
      
      call spect%fft(v, cbuff1)
      call spect%mtimes_ik1_oop(cbuff1,cbuff2)
      call spect%mtimes_ik2_oop(cbuff1,cbuff3)
      call ops%ddz_cmplx2cmplx(cbuff1)
      call spect%ifft(cbuff2,duidxj(:,:,:,2,1))
      call spect%ifft(cbuff3,duidxj(:,:,:,2,2))
      call spect%ifft(cbuff1,duidxj(:,:,:,2,3))

      !call ops%ddz(u, duidxj(:,:,:,1,3))
      !call ops%ddz(v, duidxj(:,:,:,2,3))
      !call ops%ddz(w, duidxj(:,:,:,3,3))
      
      call spect%fft(w, cbuff1)
      call spect%mtimes_ik1_oop(cbuff1,cbuff2)
      call spect%mtimes_ik2_oop(cbuff1,cbuff3)
      call ops%ddz_cmplx2cmplx(cbuff1)
      call spect%ifft(cbuff2,duidxj(:,:,:,3,1))
      call spect%ifft(cbuff3,duidxj(:,:,:,3,2))
      call spect%ifft(cbuff1,duidxj(:,:,:,3,3))

      frapid = 0.d0
      do i = 1,3
         do j = 1,3
            frapid = frapid + duidxj(:,:,:,i,j)*duidxj_mean(:,:,:,j,i)
         end do 
      end do 
      frapid = -2.d0*frapid
      call ops%dealiasField(frapid)
      call ops%SolvePoisson_oop(frapid, prapid)
     
      R11 = R11mn - u*u
      R12 = R12mn - u*v
      R13 = R13mn - u*w
      R22 = R22mn - v*v
      R23 = R23mn - v*w
      R33 = R33mn - w*w

      fslow = 0.d0
    
      ! i = 1
      call spect%fft(R11,cbuff1)
      call spect%mtimes_ik1_ip(cbuff1)
      cbuff2 = cbuff1
      call spect%fft(R12,cbuff1)
      call spect%mtimes_ik2_ip(cbuff1)
      cbuff2 = cbuff2 +cbuff1
      call spect%fft(R13,cbuff1)
      call ops%ddz_cmplx2cmplx(cbuff1)
      cbuff2 = cbuff2 +cbuff1
      call spect%mtimes_ik1_oop(cbuff2, cbuff3)

      ! i = 2
      call spect%fft(R21,cbuff1)
      call spect%mtimes_ik1_ip(cbuff1)
      cbuff2 = cbuff1
      call spect%fft(R22,cbuff1)
      call spect%mtimes_ik2_ip(cbuff1)
      cbuff2 = cbuff2 +cbuff1
      call spect%fft(R23,cbuff1)
      call ops%ddz_cmplx2cmplx(cbuff1)
      cbuff2 = cbuff2 +cbuff1
      call spect%mtimes_ik2_ip(cbuff2)
      cbuff3 = cbuff3 + cbuff2

      ! i = 3
      call spect%fft(R31,cbuff1)
      call spect%mtimes_ik1_ip(cbuff1)
      cbuff2 = cbuff1
      call spect%fft(R32,cbuff1)
      call spect%mtimes_ik2_ip(cbuff1)
      cbuff2 = cbuff2 +cbuff1
      call spect%fft(R33,cbuff1)
      call ops%ddz_cmplx2cmplx(cbuff1)
      cbuff2 = cbuff2 +cbuff1
      call ops%ddz_cmplx2cmplx(cbuff2)
      cbuff3 = cbuff3 + cbuff2

      call spect%dealias(cbuff3)
      call spect%ifft(cbuff3,fslow)
      

      call ops%SolvePoisson_oop(fslow, pslow)

      !slowp_sum = slowp_sum + pslow*pslow
      !rapidp_sum = rapidp_sum + prapid*prapid
      call ops%WriteField3D(pslow, "PSLO", tidx, RunID) 
      call ops%WriteField3D(prapid, "PRAP", tidx, RunID) 
      call mpi_barrier(mpi_comm_world, ierr)
      print*, "Just wrote fields"
      stop  

      idx = idx + 1
      tidx = tidx + tstep
      call toc()
   end do 



   slowp_sum = slowp_sum/real(idx,rkind)
   rapidp_sum = rapidp_sum/real(idx,rkind)
   
   call ops%WriteField3D(slowp_sum, "pslo", tidx, RunID) 
   call ops%WriteField3D(rapidp_sum, "prap", tidx, RunID) 
  
   
   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

