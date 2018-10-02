program createHITrestart 
   use kind_parameters, only: rkind, clen
   use mpi 
   use decomp_2d
   use decomp_2d_io
   use exits, only: message, gracefulExit
   use PadeDerOps, only: Pade6stagg
   use basic_io, only: read_2d_ascii 
   implicit none

   character(len=clen) :: inputfile, outputdir, matlabfile
   real(rkind), dimension(:,:,:), allocatable  :: fE, fC
   type(decomp_info) :: gpC, gpE
   type(Pade6stagg) :: der
   real(rkind) :: Lz, dz
   integer :: nfilters, nx, ny, nz, RID, ierr, tid, ioUnit
   integer, dimension(:), allocatable :: pid
   character(len=clen) :: tempname, fname1, fname2
   real(rkind), dimension(:,:), allocatable :: data2read
   namelist /INPUT/ nx, ny, nz, Lz, tid, RID,  outputdir, matlabfile


   call MPI_Init(ierr)               !<-- Begin MPI
   call GETARG(1,inputfile)          !<-- Get the location of the input file

   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
   read(unit=ioUnit, NML=INPUT)
   close(ioUnit)


   call decomp_2d_init(nx, ny, nz, 1, 1)
   call get_decomp_info(gpC)
   call decomp_info_init(nx,ny,nz+1,gpE)
   dz = Lz/real(nz,rkind)

   call read_2d_ascii(data2read,matlabfile)

   allocate(fE(gpE%zsz(1), gpE%zsz(2), gpE%zsz(3)))
   allocate(fC(gpC%zsz(1), gpC%zsz(2), gpC%zsz(3)))
   
   ! u fields 
   fC = reshape(data2read(:,1),[nx,ny,nz])
   write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_u.",tid
   fname1 = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
   call decomp_2d_write_one(1,fC, fname1, gpC)

   ! v fields 
   fC = reshape(data2read(:,2),[nx,ny,nz])
   write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_v.",tid
   fname1 = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
   call decomp_2d_write_one(1,fC, fname1, gpC)

   ! w fields
   fC = reshape(data2read(:,3),[nx,ny,nz])
   write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_w.",tid
   call der%init(gpC, gpC, gpE, gpE, dz, 2, .true.)
   call der%interpz_C2E(fC, fE, 0, 0)
   call decomp_2d_write_one(1,fE, fname1, gpE)

   call message(0, "Initialization RESTART files generated")
   call MPI_Finalize(ierr)    
end program
