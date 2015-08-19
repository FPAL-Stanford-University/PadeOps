#include "hitCD_files/meshgen.F90"
#include "hitCD_files/hitCD_IO.F90"
#include "hitCD_files/initfields.F90"
#include "hitCD_files/poisson.F90"

program hitcd

    use kind_parameters,  only: rkind,clen,stdout,stderr
    use exits,              only: message
    use reductions,         only: P_MAXVAL, P_SUM 
    use GridMod,            only: alloc_buffs
    use IncompressibleGrid, only: igrid
    use poissonMod, only: poisson
    use timer, only: tic, toc
    implicit none

    type(igrid), target :: igp
    character(len=clen) :: inputfile
    integer :: ierr
    logical :: isBoxFilter = .FALSE.
    type(poisson) :: POIS

    ! Start MPI
    call MPI_Init(ierr)

    call GETARG(1,inputfile)

    print*, inputfile

    ! Initialize the grid object
    call igp%init(inputfile)
    call message("Kinetic Energy:",p_sum(igp%u*igp%u + igp%v*igp%v + igp%w*igp%w)*igp%dx*igp%dy*igp%dz)
    call POIS%init( igp%decomp, igp%dx, igp%dy, igp%dz, "cd10", isBoxFilter) 

    ! Perform the initial 
    call POIS%pressureProjection(igp%u,igp%v,igp%w)
    call igp%printDivergence()

    call tic()
    do while (igp%tsim < igp%tstop) 
        call igp%AdamsBashforth()
        call POIS%pressureProjection(igp%u,igp%v,igp%w)
        igp%tsim = igp%tsim + igp%dt
        igp%step = igp%step + 1
        call message("Kinetic Energy:",p_sum(igp%u*igp%u + igp%v*igp%v + igp%w*igp%w)*igp%dx*igp%dy*igp%dz)
    end do 
    call toc()

    call igp%printDivergence()
    
    ! Destroy everything before ending
    call POIS%destroy()
    call igp%destroy()

    ! End the run
    call MPI_Finalize(ierr)

end program
