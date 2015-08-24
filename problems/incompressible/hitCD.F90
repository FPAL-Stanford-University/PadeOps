#include "hitCD_files/meshgen.F90"
#include "hitCD_files/hitCD_IO.F90"
#include "hitCD_files/initfields.F90"
#include "hitCD_files/debug.F90"
program hitcd
    use mpi
    use decomp_2d,          only: nrank 
    use kind_parameters,  only: rkind,clen,stdout,stderr
    use exits,              only: message
    use reductions,         only: P_MAXVAL, P_SUM, P_MEAN  
    use GridMod,            only: alloc_buffs
    use IncompressibleGrid, only: igrid
    use timer, only: tic, toc
    use hitCD_IO, only: dumpData4Matlab, write_matlab_header, closeMatlabFile
    use debugMod, only: debug
    use constants, only: half 
    implicit none

    type(igrid), allocatable, target :: igp
    character(len=clen) :: inputfile
    integer :: ierr
    integer :: TSTEP_DUMP = 600, TDIV_CHECK = 100

    ! Start MPI
    call MPI_Init(ierr)

    call GETARG(1,inputfile)

    allocate(igp)

    ! Initialize the grid object
    call igp%init(inputfile)
    call message("Kinetic Energy:",half*P_MAXVAL(igp%u**2 + igp%v**2 + igp%w**2))
    call igp%printDivergence()
    call write_matlab_header(igp%outputDir)
  
    ! Perform the initial 
    call dumpData4Matlab(igp%step,igp%outputDir,igp%fields(:,:,:,1:3),igp%decomp)
    call tic()
    do while (igp%tsim < igp%tstop) 
        call igp%AdamsBashforth()
        igp%tsim = igp%tsim + igp%dt
        igp%step = igp%step + 1
        
        call message(0,"Time",igp%tsim)
        call message(1,"Kinetic Energy:",half*P_MAXVAL(igp%u**2 + igp%v**2 + igp%w**2))
        if (mod(igp%step,TSTEP_DUMP) == 0) then
            call dumpData4Matlab(igp%step,igp%outputDir,igp%fields(:,:,:,1:3),igp%decomp)
        end if
        if (mod(igp%step,TDIV_CHECK) == 0) then
            call igp%printDivergence()
        end if
    end do 
    call toc()
    call igp%printDivergence()
    call closeMatlabFile

    call mpi_barrier( mpi_comm_world, ierr)  
    
    ! Destroy everything before ending
    call igp%destroy()
    deallocate(igp)

    ! End the run
    call MPI_Finalize(ierr)

end program
