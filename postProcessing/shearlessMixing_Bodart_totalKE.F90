program totalKE
    use mpi
    use kind_parameters,    only: clen, rkind
    use IncompressibleGrid, only: igrid
    use exits,              only: message

    implicit none

    type(igrid), allocatable, target :: SM
    character(len=clen) :: inputfile
    integer :: ierr, ioUnit, tidst, tiden, tstep
   
    ! Initialize MPI 
    call MPI_Init(ierr)                                                
    namelist/postProcess/ tidst, tiden, tstep    

    ! Read input file to get parameters for post processing
    call GETARG(1,inputfile)                                            
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=postProcess)
    close(ioUnit)    

    allocate(SM)                                                     
    
    call SM%init(inputfile, .true.)                              
    call SM%start_io(.true.)                                          
    call SM%printDivergence()

    do tid = tidst,tiden,tstep


    end do

    call SM%finalize_io()
    call SM%destroy()

    call MPI_Finalize(ierr)
end program totalKE
