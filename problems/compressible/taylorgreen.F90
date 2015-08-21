#include "taylorgreen_files/meshgen.F90"
#include "taylorgreen_files/initfields.F90"

program taylorgreen

    use kind_parameters,  only: rkind,clen,stdout
    use constants,        only: one
    use decomp_2d,        only: nrank
    use CompressibleGrid, only: cgrid,u_index
    use GridMod,          only: alloc_buffs, destroy_buffs
    use reductions,       only: P_MAXVAL,P_MINVAL
    use exits,            only: message,GracefulExit
    use timer,            only: tic, toc
    implicit none

    type(cgrid) :: cgp
    character(len=clen) :: inputfile
    integer :: ierr 


    ! Start MPI
    call MPI_Init(ierr)

    ! Get file location 
    call GETARG(1,inputfile)
    
    ! Initialize the grid object
    call cgp%init(inputfile)

    ! call message(1,"Maximum u velocity",P_MAXVAL(cgp%u))
    ! call cgp%filter(cgp%u,cgp%gfil,2)
    ! call message(1,"Maximum u velocity",P_MAXVAL(cgp%u))

    call message(0,"Time",cgp%tsim)
    call message(1,"Minimum density",P_MINVAL(cgp%rho))
    call message(1,"Maximum u velocity",P_MAXVAL(cgp%u))

    ! Time advance  
    do while (cgp%tsim < cgp%tstop) 
        call cgp%advance_RK45()

        ! call message(0,"Time",cgp%tsim)
        ! call message(1,"Minimum density",P_MINVAL(cgp%rho))
        ! call message(1,"Maximum u velocity",P_MAXVAL(cgp%u))
    end do 
        
    ! Destroy everythin before ending
    call cgp%destroy()

    ! End the run
    call MPI_Finalize(ierr)

end program
