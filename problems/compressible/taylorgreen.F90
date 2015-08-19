#include "taylorgreen_files/meshgen.F90"
#include "taylorgreen_files/initfields.F90"

program taylorgreen

    use kind_parameters,  only: rkind,clen,stdout
    use CompressibleGrid, only: cgrid,u_index
    use GridMod,          only: alloc_buffs, destroy_buffs
    use reductions,       only: P_MAXVAL
    use exits,            only: message
    use timer,            only: tic, toc
    implicit none

    type(cgrid) :: cgp
    character(len=clen) :: inputfile, iters
    integer :: ierr 


    ! Start MPI
    call MPI_Init(ierr)

    ! Get file location 
    call GETARG(1,inputfile)
    
    ! Initialize the grid object
    call cgp%init(inputfile)

    ! Time advance  
    do while (cgp%tsim < cgp%tstop) 
        call cgp%advance_RK45()
        cgp%tsim = cgp%tsim + cgp%dt
        cgp%step = cgp%step + 1
    end do 
        
    ! Destroy everythin before ending
    call cgp%destroy()

    ! End the run
    call MPI_Finalize(ierr)

end program
