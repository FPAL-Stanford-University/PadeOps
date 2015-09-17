#include "taylorgreen_files/meshgen.F90"
#include "taylorgreen_files/initfields.F90"

program taylorgreen

    use kind_parameters,  only: rkind,clen
    use CompressibleGrid, only: cgrid
    use gridtools,        only: alloc_buffs, destroy_buffs
    use reductions,       only: P_MAXVAL,P_MINVAL
    use exits,            only: message,GracefulExit
    use timer,            only: tic, toc
    implicit none

    type(cgrid) :: cgp
    character(len=clen) :: inputfile
    integer :: ierr 

    real(rkind) :: tke0,tke


    ! Start MPI
    call MPI_Init(ierr)

    ! Get file location 
    call GETARG(1,inputfile)
    
    ! Initialize the grid object
    call cgp%init(inputfile)

    tke0 = P_MAXVAL( cgp%rho * (cgp%u*cgp%u + cgp%v*cgp%v + cgp%w*cgp%w) )

    ! Time advance  
    do while (cgp%tsim < cgp%tstop) 
        call cgp%advance_RK45()
        
        tke = P_MAXVAL( cgp%rho * (cgp%u*cgp%u + cgp%v*cgp%v + cgp%w*cgp%w) )

        call message(2,"TKE",tke/tke0)
    end do 
        
    ! Destroy everythin before ending
    call cgp%destroy()

    ! End the run
    call MPI_Finalize(ierr)

end program
