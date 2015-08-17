#include "hitCD_files/meshgen.F90"
#include "hitCD_files/initfields.F90"

program hitcd

    use kind_parameters,  only: clen,stdout
    use IncompressibleGrid, only: igrid
    implicit none

    type(igrid) :: igp
    character(len=clen) :: inputfile
    integer :: ierr

    ! Start MPI
    call MPI_Init(ierr)

    call GETARG(1,inputfile)

    ! Initialize the grid object
    call igp%init(inputfile)

    associate( x => igp%mesh(:,:,:,1) )
        write(stdout,*) x(:,1,1)
    end associate

    ! Destroy everything before ending
    call igp%destroy()

    ! End the run
    call MPI_Finalize(ierr)

end program
