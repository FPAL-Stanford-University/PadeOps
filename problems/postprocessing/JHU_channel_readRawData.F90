module binary_io_singleprec
   implicit none
contains
   subroutine readBinarySinglePrec(filename,dat)
      character(len=*), intent(in) :: filename
      real(kind=4), dimension(:,:,:), intent(out) :: dat
      integer :: nx, ny, nz, ioUnit

      nx = size(dat,1); ny = size(dat,2); nz = size(dat,3);
      ioUnit = 25
      open(unit=ioUnit,file=trim(filename),form='unformatted')
      read(ioUnit) dat(:,:,:)
      close(ioUnit)

   end subroutine
end module

program JHU_channel
   use kind_parameters, only: rkind, clen
   use mpi
   use decomp_2d
   use timer, only: tic, toc
   use exits, only: message, gracefulexit
   use binary_io_singleprec, only: readBinarySinglePrec
   use decomp_2d_io

   implicit none
   character(len=clen) :: inputfile, inputdir, outputdir, tempname, fname
   integer :: tid, ierr, nprocsTot, ioUnit, idx
   integer :: prow, pcol
   character(len=1) :: fieldname
   type(decomp_info) :: gp
   integer, parameter :: nxChunk = 2048, nyChunk = 512, nzChunk = 6 
   real(kind=4), dimension(nxChunk,nyChunk,nzChunk) :: dataRead
   integer, parameter :: nxIn = 2048, nyIn = 512, nzIn = 1536, nchunks = nzIn/nzChunk
   integer, parameter :: nxOut = 2048, nyOut = 1536, nzOut = 512
   
   integer :: mychunksize, my_yidx_start, yidx, j, k
   real(rkind), dimension(:,:,:), allocatable :: rawData

   namelist /INPUT/ inputdir, outputdir, tid, fieldname

   call MPI_init(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocsTot, ierr)

   call GETARG(1,inputfile)
   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=ioUnit, NML=INPUT)
   close(ioUnit)
   
   
   pcol = 1; prow = nprocsTot;
   call decomp_2d_init(nxOut,nyOut,nzOut,prow,pcol,[.TRUE.,.TRUE.,.FALSE.])
   call get_decomp_info(gp)

   if (mod(gp%xsz(2),nzChunk) .ne. 0) then
      call message(0,"y-size of array:",gp%xsz(2))
      call gracefulexit("Check the number of processors.",0)
   end if
   
   call message(0, "Decomp initialized successfully.")
   mychunksize = gp%xsz(2)/nzChunk
   call message(1, "Number of chunks on each processor:",mychunksize)
   call message(1, "nx (local):",gp%xsz(1))
   call message(1, "ny (local):",gp%xsz(2))
   call message(1, "nz (local):",gp%xsz(3))
   my_yidx_start = nrank*gp%xsz(2) + 1

   allocate(rawData(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
   call message(0,"Successfully allocated the double precision data storage.")
   call tic;
   yidx = 1
   do idx = 1, mychunksize 
      write(tempname,"(A6,I3.3,A1,A1,A2,I3.3,A4,I4.4,A4,I4.4,A4)") "data_t",tid,"/",fieldname,"_t",tid,"_zst", &
                                    & (idx-1)*nzChunk+my_yidx_start,"_zen",idx*nzChunk+my_yidx_start-1,".bin"
      fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
      call readBinarySinglePrec(fname,dataRead)
      do j = 1,nzChunk
         do k = 1,gp%xsz(3)
            rawData(:,yidx,k) = real(dataRead(:,k,j),rkind)
         end do
         yidx = yidx + 1
      end do
   end do 
   call toc
   call message(0,"Data reading complete.")

   !! Dump planes
   call decomp_2d_write_plane(1,rawData, 1, nxOut/2, "TestPlane_x.dat")
   call decomp_2d_write_plane(1,rawData, 2, nyOut/2, "TestPlane_y.dat")
   call decomp_2d_write_plane(1,rawData, 3, nzOut/4, "TestPlane_z.dat")


   call MPI_Barrier(mpi_comm_world, ierr)
   call MPI_Finalize(ierr)   
end program 
