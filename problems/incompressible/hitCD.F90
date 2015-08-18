#include "hitCD_files/meshgen.F90"
#include "hitCD_files/hitCD_IO.F90"
#include "hitCD_files/initfields.F90"
#include "hitCD_files/poisson.F90"

program hitcd

    use kind_parameters,  only: rkind,clen,stdout,stderr
    use IncompressibleGrid, only: igrid
    use poissonMod, only: poisson
    use decomp_2d, only: nrank
    implicit none

    type(igrid), target :: igp
    character(len=clen) :: inputfile
    integer :: ierr
    logical :: isBoxFilter = .FALSE.
    type(poisson) :: POIS
    real(rkind), dimension(:,:,:), pointer :: u,v,w

    ! Start MPI
    call MPI_Init(ierr)

    call GETARG(1,inputfile)

    print*, inputfile

    ! Initialize the grid object
    call igp%init(inputfile)
    call POIS%init( igp%decomp, igp%dx, igp%dy, igp%dz, "cd10", isBoxFilter) 

    ! Associate the pointers for ease of use
    u => igp%fields(:,:,:,1); v => igp%fields(:,:,:,2); w => igp%fields(:,:,:,3)

    ! Destroy everything before ending
    call POIS%destroy()
    call igp%destroy()

    ! End the run
    call MPI_Finalize(ierr)

end program
