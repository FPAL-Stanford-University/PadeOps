#include "hitCD_files/meshgen.F90"
#include "hitCD_files/hitCD_IO.F90"
#include "hitCD_files/initfields.F90"
#include "hitCD_files/poisson.F90"

program hitcd

    use kind_parameters,  only: rkind,clen,stdout,stderr
    use exits,              only: message
    use reductions,         only: P_MAXVAL
    use GridMod,            only: alloc_buffs
    use IncompressibleGrid, only: igrid
    use poissonMod, only: poisson
    use decomp_2d, only: nrank
    implicit none

    type(igrid), target :: igp
    character(len=clen) :: inputfile
    integer :: ierr
    logical :: isBoxFilter = .FALSE.
    type(poisson) :: POIS
    real(rkind), dimension(:,:,:), pointer :: u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    real(rkind), dimension(:,:,:,:), allocatable, target :: duidxj

    ! Start MPI
    call MPI_Init(ierr)

    call GETARG(1,inputfile)

    print*, inputfile

    ! Initialize the grid object
    call igp%init(inputfile)
    call POIS%init( igp%decomp, igp%dx, igp%dy, igp%dz, "cd10", isBoxFilter) 

    ! Allocate duidxj
    call alloc_buffs(duidxj,9,'y',igp%decomp)

    ! Associate the pointers for ease of use
    u => igp%fields(:,:,:,1); v => igp%fields(:,:,:,2); w => igp%fields(:,:,:,3)

    dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
    dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
    dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);

    call POIS%pressureProjection(igp%fields(:,:,:,1:3))

    call igp%gradient(u,dudx,dudy,dudz)
    call igp%gradient(v,dvdx,dvdy,dvdz)
    call igp%gradient(w,dwdx,dwdy,dwdz)
    call message("Maximum divergence",P_MAXVAL(ABS(dudx+dvdy+dwdz)))
    
    ! Destroy everything before ending
    call POIS%destroy()
    call igp%destroy()
    deallocate(duidxj)

    ! End the run
    call MPI_Finalize(ierr)

end program
