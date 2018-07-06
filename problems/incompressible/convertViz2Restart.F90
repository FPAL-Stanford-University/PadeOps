program convertViz2Restart 
    use kind_parameters, only: rkind, clen
    use gridtools, only: alloc_buffs, destroy_buffs
    use decomp_2d
    use exits, only: GracefulExit, message
    use mpi 
    use timer, only: tic, toc
    use decomp_2d_io
    use constants, only: pi 
    use PadeDerOps, only: Pade6stagg
    use spectralMod, only: spectral  
    implicit none

    character(len=clen) :: inputfile 
    character(len=clen) :: outputdir, inputdir
    integer :: ioUnit, nx, ny, nz, ierr, inputFile_TID, inputFile_RID, restartfile_TID, restartfile_RID
    logical :: isStratified = .false. 
    type(decomp_info) :: gpC, gpE
    type(decomp_info), pointer :: sp_gpC, sp_gpE
    type(spectral), target :: spectC, spectE
    character(len=clen) :: tempname, fname
    real(rkind), dimension(:,:,:), allocatable :: rbuffxC, rbuffyC, rbuffzC, rbuffzE
    real(rkind) :: tsim 
    logical :: periodicInZ = .false. 
    real(rkind) :: Lz = 2.d0*pi, dz
    type(Pade6stagg) :: DerOps
    integer :: NumericalSchemeVert = 2

    namelist /INPUT/ nx, ny, nz, inputdir, outputdir, inputFile_TID, inputFile_RID, NumericalSchemeVert, &
      & restartfile_TID, restartfile_RID, isStratified, PeriodicInZ, Lz

    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=INPUT)
    close(ioUnit)

    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(gpC)
    call decomp_info_init(nx,ny,nz+1,gpE)

    allocate(rbuffxC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(rbuffyC(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3)))
    allocate(rbuffzC(gpC%zsz(1),gpC%zsz(2),gpC%zsz(3)))
    allocate(rbuffzE(gpE%zsz(1),gpE%zsz(2),gpE%zsz(3)))
    dz = Lz/real(nz,rkind)

    call spectC%init("x", nx, ny, nz  , 1.d0, 1.d0, dz, &
           & "FOUR", "2/3rd", 2 , fixOddball=.false., exhaustiveFFT=.true., init_periodicInZ=periodicinZ, dealiasF=0.667d0)
    call spectE%init("x", nx, ny, nz+1, 1.d0, 1.d0, dz, &
           & "FOUR", "2/3rd", 2 , fixOddball=.false., exhaustiveFFT=.true., init_periodicInZ=.false., dealiasF=0.667d0)

    sp_gpC => spectC%spectdecomp
    sp_gpE => spectE%spectdecomp
    call DerOps%init(gpC, sp_gpC, gpE, sp_gpE, dz, NumericalSchemeVert, PeriodicInZ, spectC)

    !!!!!! READ HEADER !!!!!!!!!!!
    if (nrank == 0) then
      !write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",inputFile_RID, "_info.",inputFile_TID
      write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",inputFile_RID, "_","info","_t",inputFile_TID,".out"
      fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
      open(unit=10,file=fname,access='sequential',form='formatted')
      read (10, *)  tsim
      close(10)
      call message(0, "Reading visualization data dumped at tSIM=", tsim)
    end if 

    call mpi_barrier(mpi_comm_world, ierr)
    !!!!! CELL FIELDS !!!!!!!!!!!!

    ! Read and Write u - field
    write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",inputFile_RID, "_","uVel","_t",inputFile_TID,".out"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,rbuffxC,fname, gpC)
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",restartfile_RID, "_u.",restartfile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
    call decomp_2d_write_one(1,rbuffxC,fname,gpC)

    write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",inputFile_RID, "_","vVel","_t",inputFile_TID,".out"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,rbuffxC,fname, gpC)
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",restartfile_RID, "_v.",restartfile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
    call decomp_2d_write_one(1,rbuffxC,fname,gpC)
    
    if (isStratified) then
      write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",inputFile_RID, "_","potT","_t",inputFile_TID,".out"
      fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
      call decomp_2d_read_one(1,rbuffxC,fname, gpC)
      write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",restartfile_RID, "_T.",restartfile_TID
      fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
      call decomp_2d_write_one(1,rbuffxC,fname,gpC)
    end if
   
    call message(0,"Now converting w to an edge field.")
    call mpi_barrier(mpi_comm_world, ierr)

    !!!!!!!!!!!!! EDGE FIELDS !!!!!!!!!!!!!!!
    ! Read and Write w - field
    write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",inputFile_RID, "_","wVel","_t",inputFile_TID,".out"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,rbuffxC,fname, gpC)
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",restartfile_RID, "_w.",restartfile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
   
    call transpose_x_to_y(rbuffxC, rbuffyC, gpC)
    call transpose_y_to_z(rbuffyC, rbuffzC, gpC)
    call DerOps%interpz_C2E(rbuffzC, rbuffzE, -1, -1)
    call decomp_2d_write_one(3,rbuffzE,fname,gpE)
    
    call mpi_barrier(mpi_comm_world, ierr)

    !!!!!! WRITE HEADER !!!!!!!!!!!
    if (nrank == 0) then
        write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",restartfile_RID, "_info.",restartfile_TID
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
        OPEN(UNIT=10, FILE=trim(fname))
        write(10,"(100g15.5)") tsim
        close(10)
    end if 
    call message(0,"Restart files generated.")

    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program 
