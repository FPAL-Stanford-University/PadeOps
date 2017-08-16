program upsampleFields
    use kind_parameters, only: rkind, clen
    use gridtools, only: alloc_buffs, destroy_buffs
    use decomp_2d
    use exits, only: GracefulExit, message
    use mpi 
    use timer, only: tic, toc
    use decomp_2d_io
    use reductions, only: p_sum
    implicit none

    character(len=clen) :: inputfile 
    character(len=clen) :: outputdir, inputdir
    integer :: ioUnit, n, ierr, inputFile_TID, inputFile_RID, outputFile_TID, outputFile_RID
    type(decomp_info) :: gpC, gpE
    real(rkind), dimension(:,:,:), allocatable :: f
    character(len=clen) :: tempname, fname
    real(rkind) :: setSimTime, tsim, mean_f, newFrameSpeed , ScalingFact = 1.d0
    namelist /INPUT/ n, inputdir, outputdir, inputFile_TID, inputFile_RID, &
    outputFile_TID, outputFile_RID, newFrameSpeed, setSimTime, scalingFact

    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=INPUT)
    close(ioUnit)


    call decomp_2d_init(n, n, n, 0, 0)
    call get_decomp_info(gpC)
    call decomp_info_init(n,n,n+1,gpE)
    
    !!!!!! READ HEADER !!!!!!!!!!!
    write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",inputFile_RID, "_info.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    open(unit=10,file=fname,access='sequential',form='formatted')
    read (10, *)  tsim
    close(10)
    call message(0, "Input File data dumped at tSIM=", tsim)

    !!!!!!!!!!!!! CELL FIELDS !!!!!!!!!!!!!!!
    allocate(f(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
   
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_u.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,f,fname, gpC)
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_u.",outputFile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
    mean_f = p_sum(sum(f))/(real(n)**3)
    call message(0, "Original frame speed", mean_f)
    !f = f + (newFrameSpeed - mean_f)
    f = f - mean_f
    f = f*scalingfact
    f = f + newFrameSpeed
    mean_f = p_sum(sum(f))/(real(n)**3)
    call message(0, "New frame speed", mean_f)
    call decomp_2d_write_one(1,f,fname, gpC)
    call message(0, "Done writing u")
   
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_v.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,f,fname, gpC)
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_v.",outputFile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
    f = f*scalingfact
    mean_f = p_sum(sum(f))/(real(n)**3)
    f = f - mean_f
    call decomp_2d_write_one(1,f,fname, gpC)
    
    deallocate(f)
    call message(0, "Done writing v")
    
    allocate(f(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
    
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_w.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,f,fname, gpE)
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_w.",outputFile_TID
    f = f*scalingfact
    mean_f = p_sum(sum(f))/(real(n)**3)
    f = f - mean_f
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
    call decomp_2d_write_one(1,f,fname, gpE)
    
    deallocate(f)

    call message(0, "Done writing w")
    !!!!!! WRITE HEADER !!!!!!!!!!!
    if (nrank == 0) then
        write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",outputFile_RID, "_info.",outputFile_TID
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
        OPEN(UNIT=10, FILE=trim(fname))
        write(10,"(100g15.5)") setSimTime
        close(10)
    end if 

    call message(0, "Done writing header")
    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program 
