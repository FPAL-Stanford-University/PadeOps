program JHU_channel_postprocess
   use kind_parameters, only: rkind, clen
   use mpi
   use decomp_2d
   use timer, only: tic, toc
   use exits, only: message, gracefulexit
   use decomp_2d_io
   use spectralmod, only: spectral
   use constants, only: pi
   use basic_io 

   implicit none
   character(len=clen) :: inputfile, inputdir, outputdir, tempname, fname, zLocsFile
   real(rkind), parameter :: Lx = 8.d0*pi, Ly = 3.d0*pi, Lz = 2.d0
   integer, parameter :: nx = 2048, ny = 1536, nz = 512
   real(rkind), parameter :: dx = Lx/real(nx,rkind), dy = Ly/real(ny,rkind), dz = 1.d0
   type(spectral), allocatable, target :: spect
   integer :: tid, ierr, ioUnit
   type(decomp_info) :: gp
   type(decomp_info), pointer :: sp_gp
   real(rkind), dimension(:,:,:), allocatable :: u, v, w, ufil, vfil, wfil
   complex(rkind), dimension(:,:,:), allocatable :: cbuffy
   real(rkind), dimension(:,:,:), allocatable :: rbuffx, rbuffy, rbuffz
   real(rkind) :: filterfact_xy = 0.1d0
   real(rkind), dimension(nz,12) :: meanprofiles
   real(rkind), dimension(:,:), allocatable :: zLocs

   namelist /INPUT/ inputdir, outputdir, zLocsFile,tid, filterfact_xy

   call MPI_init(ierr)

   call GETARG(1,inputfile)
   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=ioUnit, NML=INPUT)
   close(ioUnit)
   
   
   call decomp_2d_init(nx,ny,nz,0,0,[.TRUE.,.TRUE.,.FALSE.])
   call get_decomp_info(gp)

   ! Get the location of collocation points
   call read_2d_ascii(zLocs,zLocsFile)
   meanprofiles(:,1) = zLocs(:,1)
   deallocate(zLocs)

   !! Initialize the Spectral derived type
   allocate(spect)
   call spect%init("x", nx, ny, nz, dx, dy, dz, &
                "four", "2/3rd", 2 , .false.)
   sp_gp => spect%spectdecomp
   call spect%initPP(filterfact_xy, filterfact_xy, dx, dy)

   call message(0,"Spectral Derived type initialized.")
   
   allocate(u(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
   allocate(v(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
   allocate(w(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
   allocate(rbuffx(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
   allocate(rbuffy(gp%ysz(1),gp%ysz(2),gp%ysz(3)))
   allocate(rbuffz(gp%zsz(1),gp%zsz(2),gp%zsz(3)))
   allocate(ufil(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
   allocate(vfil(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
   allocate(wfil(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
   call spect%alloc_r2c_out(cbuffy)
   call message(0,"Allocated the data arrays.")
   

   write(tempname,"(A1,A2,I3.3,A15)") "u","_t",tid,"_JHUchannel.bin"
   fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
   call decomp_2d_read_one(1,u,fname, gp)
   
   write(tempname,"(A1,A2,I3.3,A15)") "v","_t",tid,"_JHUchannel.bin"
   fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
   call decomp_2d_read_one(1,v,fname, gp)

   write(tempname,"(A1,A2,I3.3,A15)") "w","_t",tid,"_JHUchannel.bin"
   fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
   call decomp_2d_read_one(1,w,fname, gp)
   
   call message(0,"Data reading complete.")
   
   call message("============================================")
   call message(0,"Now starting analysis")
   call message("============================================")

   call message(0,"Step 1: Get the 3 mean profiles")
   call meanZ(u,meanprofiles(:,2))
   call meanZ(v,meanprofiles(:,3))
   call meanZ(w,meanprofiles(:,4))
   
   call message(0,"Step 2: FFT(), filter(), IFFT()")
   call filter(u,ufil)
   call filter(v,vfil)
   call filter(w,wfil)

   call message(0,"Step 3: get correlations for baseline")
   call getcorrelation(u,u,meanprofiles(:,5))
   call getcorrelation(v,v,meanprofiles(:,6))
   call getcorrelation(w,w,meanprofiles(:,7))
   call getcorrelation(u,w,meanprofiles(:,8))

   call message(0,"Step 4: get correlations for filtered")
   call getcorrelation(ufil,ufil,meanprofiles(:,9))
   call getcorrelation(vfil,vfil,meanprofiles(:,10))
   call getcorrelation(wfil,wfil,meanprofiles(:,11))
   call getcorrelation(ufil,wfil,meanprofiles(:,12))

   call message(0,"Step 5: Write statistics to an ascii file")
   write(tempname,"(A7,I3.3,A15)") "stats_t",tid,"_JHUchannel.txt"
   fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
   if (nrank == 0) call write_2d_ascii(meanprofiles,fname)

   call message(0,"Verification - write filtered fields.")
   call message(1,"u - velocity.")
   write(tempname,"(A2,A2,I3.3,A15)") "uF","_t",tid,"_JHUchannel.bin"
   fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
   call decomp_2d_write_one(1,ufil,fname, gp)
  
   call message(1,"v - velocity.")
   write(tempname,"(A2,A2,I3.3,A15)") "vF","_t",tid,"_JHUchannel.bin"
   fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
   call decomp_2d_write_one(1,vfil,fname, gp)

   call message(1,"w - velocity.")
   write(tempname,"(A2,A2,I3.3,A15)") "wF","_t",tid,"_JHUchannel.bin"
   fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
   call decomp_2d_write_one(1,wfil,fname, gp)

   deallocate(u,v,w, ufil, vfil, wfil, cbuffy, rbuffy, rbuffx, rbuffz)
   call spect%destroy()
   deallocate(spect)


   call MPI_Barrier(mpi_comm_world, ierr)
   call MPI_Finalize(ierr)   

contains
  
   subroutine getcorrelation(f1,f2,correlation)
      real(rkind), intent(in), dimension(gp%xsz(1),gp%xsz(2),gp%xsz(3)) :: f1, f2
      real(rkind), dimension(nz) :: mn1, mn2, mncross
      real(rkind), intent(out), dimension(nz) :: correlation
      call meanZ(f1,mn1)
      call meanZ(f2,mn2)
      
      rbuffx = f1*f2
      call meanZ(rbuffx,mncross)
      correlation = mncross - mn1*mn2

   end subroutine

   subroutine filter(f,ffil)
      real(rkind), intent(in), dimension(gp%xsz(1),gp%xsz(2),gp%xsz(3)) :: f
      real(rkind), intent(out), dimension(gp%xsz(1),gp%xsz(2),gp%xsz(3)) :: ffil

      call spect%fft(f,cbuffy)
      call spect%spectralfilter_ip(cbuffy)
      call spect%ifft(cbuffy,ffil)
      

   end subroutine

   subroutine meanZ(f,meanf)
      use reductions, only: p_sum
      real(rkind), intent(in), dimension(gp%xsz(1),gp%xsz(2),gp%xsz(3)) :: f
      real(rkind), intent(out), dimension(nz) :: meanf
      integer :: k
      real(rkind) :: xysum

      call transpose_x_to_y(f,rbuffy,gp)
      call transpose_y_to_z(rbuffy,rbuffz,gp)
      
      do k = 1,nz
         xysum = p_sum(rbuffz(:,:,k))
         meanf(k) = xysum/real(real(nx)*real(ny),rkind)
      end do 

   end subroutine

end program
