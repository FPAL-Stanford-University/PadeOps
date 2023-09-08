module convertViz2RestartMod
  use kind_parameters, only: clen, rkind
  use decomp_2d
  use decomp_2d_io
  use PadeDerOps, only: Pade6stagg
  implicit none
  character(len=clen) :: outputdir, inputdir
  integer :: inputFile_RID, inputFile_TID, restartfile_TID, restartfile_RID
  type(Pade6stagg) :: DerOps

  contains
    subroutine convertCellField(vizName,restartName,rbuffx,gp)
        character(len=4), intent(in) :: vizName
        character(len=*), intent(in) :: restartName
        real(rkind), dimension(:,:,:), intent(inout) :: rbuffx
        type(decomp_info), intent(inout) :: gp
        character(len=clen) :: tempname, fname

        write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",inputFile_RID, "_",vizNAme,"_t",inputFile_TID,".out"
        fname = inputdir(:len_trim(inputdir))//"/"//trim(tempname)
        call decomp_2d_read_one(1,rbuffx,fname, gp)
        write(tempname,"(A7,A4,I2.2,A1,A,A1,I6.6)") "RESTART", "_Run",restartfile_RID, "_",restartName,".",restartfile_TID
        fname = outputdir(:len_trim(outputdir))//"/"//trim(tempname)
        call decomp_2d_write_one(1,rbuffx,fname,gp)
    end subroutine

    subroutine convertEdgeField(vizName,restartName,rbuffxC,rbuffyC,rbuffzC,&
        rbuffzE,gpC,gpE)
        character(len=4), intent(in) :: vizName
        character(len=*), intent(in) :: restartName
        real(rkind), dimension(:,:,:), intent(inout) :: rbuffxC, rbuffyC, &
          rbuffzC, rbuffzE
        type(decomp_info), intent(inout) :: gpC, gpE
        character(len=clen) :: tempname, fname
      
        write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",inputFile_RID, "_",&
          vizName,"_t",inputFile_TID,".out"
        fname = inputdir(:len_trim(inputdir))//"/"//trim(tempname)
        call decomp_2d_read_one(1,rbuffxC,fname, gpC)
        write(tempname,"(A7,A4,I2.2,A1,A,A1,I6.6)") "RESTART", "_Run",restartfile_RID, "_",&
          restartName,".",restartfile_TID
        fname = outputdir(:len_trim(outputdir))//"/"//trim(tempname)
   
        call transpose_x_to_y(rbuffxC, rbuffyC, gpC)
        call transpose_y_to_z(rbuffyC, rbuffzC, gpC)
        call DerOps%interpz_C2E(rbuffzC, rbuffzE, -1, -1)
        call decomp_2d_write_one(3,rbuffzE,fname,gpE)
    end subroutine

end module convertViz2RestartMod

program convertViz2Restart 
    use kind_parameters, only: rkind, clen
    use gridtools, only: alloc_buffs, destroy_buffs
    use decomp_2d
    use exits, only: GracefulExit, message
    use mpi 
    use timer, only: tic, toc
    use decomp_2d_io
    use constants, only: pi 
    use spectralMod, only: spectral  
    use convertViz2RestartMod, only: inputdir, outputdir, convertCellField, &
      inputFile_TID, inputFile_RID, restartfile_TID, restartfile_RID, &
      convertEdgeField, DerOps
    implicit none

    character(len=clen) :: inputfile 
    integer :: ioUnit, nx, ny, nz, ierr, n 
    logical :: isStratified = .false. 
    type(decomp_info) :: gpC, gpE
    type(decomp_info), pointer :: sp_gpC, sp_gpE
    type(spectral), target :: spectC, spectE
    character(len=clen) :: tempname, fname
    real(rkind), dimension(:,:,:), allocatable :: rbuffxC, rbuffyC, rbuffzC, rbuffzE
    real(rkind) :: tsim 
    character(len=4) :: vizname
    character(len=8) :: rstname
    logical :: periodicInZ = .false. 
    real(rkind) :: Lz = 2.d0*pi, dz
    integer :: NumericalSchemeVert = 2, num_scalars = 0
    logical :: usingLocalizedForceLayer = .false.

    namelist /INPUT/ nx, ny, nz, inputdir, outputdir, inputFile_TID, inputFile_RID, NumericalSchemeVert, &
      & restartfile_TID, restartfile_RID, isStratified, PeriodicInZ, Lz, usingLocalizedForceLayer, &
      num_scalars

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
      fname = inputdir(:len_trim(inputdir))//"/"//trim(tempname)
      open(unit=10,file=fname,access='sequential',form='formatted')
      read (10, *)  tsim
      close(10)
      call message(0, "Reading visualization data dumped at tSIM=", tsim)
    end if 

    call mpi_barrier(mpi_comm_world, ierr)
    
    !!!!! CELL FIELDS !!!!!!!!!!!
    call convertCellField('uVel','u',rbuffxC,gpC)
    call convertCellField('vVel','v',rbuffxC,gpC)
    if (isStratified) call convertCellField('potT','T',rbuffxC,gpC)
    if (num_scalars > 0) then
        do n = 1,num_scalars
            write(vizname,'(A2,I2.2)')'sc',n
            write(rstname,'(A6,I2.2)')'SCALAR',n
            call convertCellField(vizname,rstname,rbuffxC,gpC)
        end do
    end if

    if (usingLocalizedForceLayer) then
      call convertCellField('frcx','fx',rbuffxC,gpC)
      call convertCellField('frcy','fy',rbuffyC,gpC)
    end if
   
    call message(0,"Now converting w to an edge field.")
    call mpi_barrier(mpi_comm_world, ierr)

    !!!!!!!!!!!!! EDGE FIELDS !!!!!!!!!!!!!!!
    call convertEdgeField('wVel','w',rbuffxC,rbuffyC,rbuffzC,rbuffzE,gpC,gpE)
    if (usingLocalizedForceLayer) call convertEdgeField('frcz','fz',rbuffxC,&
      rbuffyC,rbuffzC,rbuffzE,gpC,gpE)
    
    call mpi_barrier(mpi_comm_world, ierr)

    !!!!!! WRITE HEADER !!!!!!!!!!!
    if (nrank == 0) then
        write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",restartfile_RID, "_info.",restartfile_TID
        fname = outputdir(:len_trim(outputdir))//"/"//trim(tempname)
        OPEN(UNIT=10, FILE=trim(fname))
        write(10,"(100g15.5)") tsim
        close(10)
    end if 
    call message(0,"Restart files generated.")

    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program 
