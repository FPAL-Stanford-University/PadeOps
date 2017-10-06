program filter_fields_periodic
   use kind_parameters, only: rkind, clen
   use mpi 
   use decomp_2d
   use decomp_2d_io
   use fof_mod, only: fof
   use spectralmod, only: spectral
   use exits, only: message

   implicit none

   character(len=clen) :: fof_dir, datadir, inputfile, tempname, fname
   real(rkind), dimension(:,:,:), allocatable  :: u, v, w, p, rbuffx
   complex(rkind), dimension(:,:,:), allocatable :: cbuffy, cbuffz, uhat, vhat, what
   type(decomp_info) :: gpC
   type(spectral) :: spectC
  
   type(fof), dimension(:), allocatable :: filt
   real(rkind) :: Lx, Ly, Lz, dx, dy, dz
   integer :: nfilters, nx, ny, nz, runID, ierr, fid, tid, ioUnit, nplanesy
   integer, dimension(:), allocatable :: pid
   namelist /FILTER_INFO/ nfilters, fof_dir
   namelist /INPUT/ nx, ny, nz, Lx, Ly, Lz,  tid, runID, datadir


   call MPI_Init(ierr)               !<-- Begin MPI
   call GETARG(1,inputfile)          !<-- Get the location of the input file

   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
   read(unit=ioUnit, NML=INPUT)
   read(unit=ioUnit, NML=FILTER_INFO)
   close(ioUnit)


   call decomp_2d_init(nx, ny, nz, 0, 0)
   call get_decomp_info(gpC)

   dx = Lx/real(nx,rkind)
   dy = Ly/real(ny,rkind)
   dz = Lz/real(nz,rkind)

   call spectC%init("x",nx,ny,nz,dx,dy,dz,"four","2/3rd",2,.false.,.false.,.false.,.true.,.true.,.true.)

   allocate(u(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
   allocate(v(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
   allocate(w(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
   allocate(p(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
   allocate(rbuffx(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
   call spectC%alloc_r2c_out(cbuffy)
   call spectC%alloc_r2c_out(uhat)
   call spectC%alloc_r2c_out(vhat)
   call spectC%alloc_r2c_out(what)
   
   allocate(cbuffz(spectC%spectdecomp%zsz(1),spectC%spectdecomp%zsz(2),spectC%spectdecomp%zsz(3)))
   allocate(filt(nfilters))

   nplanesy = 1
   allocate(pid(nplanesy))
   pid(1) = 128

   ! Read in the fields here; transform to their fourier tranforms
   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",runID, "_","uVel","_t",tid,".out"
   fname = DataDir(:len_trim(DataDir))//"/"//trim(tempname)
   call decomp_2d_read_one(1,u,fname, gpC)
   
   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",runID, "_","vVel","_t",tid,".out"
   fname = DataDir(:len_trim(DataDir))//"/"//trim(tempname)
   call decomp_2d_read_one(1,v,fname, gpC)

   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",runID, "_","wVel","_t",tid,".out"
   fname = DataDir(:len_trim(DataDir))//"/"//trim(tempname)
   call decomp_2d_read_one(1,w,fname, gpC)

   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",runID, "_","prss","_t",tid,".out"
   fname = DataDir(:len_trim(DataDir))//"/"//trim(tempname)
   call decomp_2d_read_one(1,p,fname, gpC)


   call spectC%fft(u,uhat)
   call spectC%fft(v,vhat)
   call spectC%fft(w,what)

   do fid = 1,nfilters
      call filt(fid)%init(runID, fof_dir, datadir, fid, spectC, cbuffy, cbuffz, .true., gpC)
   end do 


   ! Do the filtering here, along with plane and full field writing
   do fid = 1,nfilters
      call message(1,"Now applying filter number:",fid)
      call filt(fid)%filter_Complex2Real(uhat, rbuffx)
      call filt(fid)%dumpFullField(rbuffx, "uVel",tid)
      call filt(fid)%dumpYplanes(rbuffx, "uVel",pid, tid)

      call filt(fid)%filter_Complex2Real(vhat, rbuffx)
      call filt(fid)%dumpFullField(rbuffx, "vVel",tid)
      call filt(fid)%dumpYplanes(rbuffx, "vVel",pid, tid)

      call filt(fid)%filter_Complex2Real(what, rbuffx)
      call filt(fid)%dumpFullField(rbuffx, "wVel",tid)
      call filt(fid)%dumpYplanes(rbuffx, "wVel",pid, tid)

      call filt(fid)%filter_Real2Real(p, rbuffx)
      call filt(fid)%dumpFullField(rbuffx, "press",tid)
      call filt(fid)%dumpYplanes(rbuffx, "press",pid, tid)
   end do 


   do fid = 1,nfilters
      call filt(fid)%destroy()
   end do 


   deallocate(u, v, w, p, uhat, vhat, what, rbuffx, cbuffy, cbuffz)


   call spectC%destroy()
   call MPI_Finalize(ierr)    

end program
