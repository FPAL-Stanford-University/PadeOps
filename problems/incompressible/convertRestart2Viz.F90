program convert_restart_to_viz
   use kind_parameters, only: rkind, clen
   use mpi 
   use decomp_2d
   use decomp_2d_io
   use exits, only: message, gracefulExit
   use PadeDerOps, only: Pade6stagg 
   implicit none

   character(len=clen) :: inputfile, outputdir, inputdir
   real(rkind), dimension(:,:,:), allocatable  :: fE, fC
   type(decomp_info) :: gpC, gpE
   type(Pade6stagg) :: der
   logical :: isPeriodicInZ = .false., isStratified = .false.  
   real(rkind) :: Lz, dz
   integer :: nfilters, nx, ny, nz, RID, ierr, tid, ioUnit, interpolationScheme
   integer, dimension(:), allocatable :: pid
   character(len=clen) :: tempname, fname1, fname2
   namelist /INPUT/ nx, ny, nz, Lz, tid, RID, inputdir, outputdir, isPeriodicInZ, isStratified, interpolationScheme


   call MPI_Init(ierr)               !<-- Begin MPI
   call GETARG(1,inputfile)          !<-- Get the location of the input file

   ioUnit = 11
   open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
   read(unit=ioUnit, NML=INPUT)
   close(ioUnit)


   call decomp_2d_init(nx, ny, nz, 0, 0)
   call get_decomp_info(gpC)
   call decomp_info_init(nx,ny,nz+1,gpE)
   dz = Lz/real(nz,rkind)

   if (isPeriodicInZ) then
        call gracefulExit("Currently this program only supports non-periodic in z problems",314)
   end if 

   
   write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_u.",tid
   fname1 = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",rid, "_","uVel","_t",tid,".out"
   fname2 = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
   if (nrank == 0) then 
        call system("cp " // trim(fname1) //" "// trim(fname2))
   end if 
 
   write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_v.",tid
   fname1 = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",rid, "_","vVel","_t",tid,".out"
   fname2 = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
   if (nproc > 1) then
      if (nrank == 1) then 
           call system("cp " // trim(fname1) //" "// trim(fname2))
      end if 
   else
      if (nrank == 0) then 
           call system("cp " // trim(fname1) //" "// trim(fname2))
      end if 
   end if  


   if (isStratified) then
      write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_T.",tid
      fname1 = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
      write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",rid, "_","potT","_t",tid,".out"
      fname2 = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
      if (nproc > 2) then
         if (nrank == 2) then 
              call system("cp " // trim(fname1) //" "// trim(fname2))
         end if 
      else
         if (nrank == 0) then 
              call system("cp " // trim(fname1) //" "// trim(fname2))
         end if 
      end if  
   end if 
   call mpi_barrier(mpi_comm_world, ierr) 
   call message(0, "Now interpolating w from edge to cells")
   write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",rid, "_w.",tid
   fname1 = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",rid, "_","wVel","_t",tid,".out"
   fname2 = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)

   allocate(fE(gpE%zsz(1), gpE%zsz(2), gpE%zsz(3)))
   allocate(fC(gpC%zsz(1), gpC%zsz(2), gpC%zsz(3)))
   call decomp_2d_read_one(3,fE,fname1, gpE)
  
   call der%init(gpC, gpC, gpE, gpE, dz, interpolationScheme, isPeriodicinZ)
   call der%interpz_E2C(fE, fC, -1, -1)
   call message(1, "Done interpolating. Now writing to disk.")
   call decomp_2d_write_one(3,fC,fname2, gpC)

   call message(0, "Conversion complete.")
   call MPI_Finalize(ierr)    

end program
