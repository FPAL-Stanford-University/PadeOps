#include "ShearLayer_files/hooks.F90"

program ShearLayer_seed

    use hdf5
    use kind_parameters,  only: rkind,clen
    use CompressibleGrid, only: cgrid
    use exits,            only: GracefulExit, message
    implicit none

    type(cgrid) :: cgp,seed
    character(len=clen) :: inputfile,seedfile
    integer :: ierr, error 

    ! Start MPI
    call MPI_Init(ierr)

    ! Initialize the HDF5 library and Fortran interfaces
    call h5open_f(error)
    if (error /= 0) call GracefulExit("Could not initialize HDF5 library and Fortran interfaces.",7356)

    ! Get file location 
    call message("Initializing grid objects...")
    call GETARG(1,inputfile)
    call cgp%init(inputfile)
    
    ! Initialize the grid objects
    call message("Initializing seed objects...")
    call GETARG(2,seedfile)
    call seed%init(seedfile)

    ! Seed
    call message("Seeding from file...")
    call cgp%seed_turb(seed)
    
    ! Destroy everythin before ending
    call message("Destroying grid objects...")
    call cgp%destroy_grid()
    call seed%destroy_grid()

    ! Close Fortran interfaces and HDF5 library
    call h5close_f(error)

    ! End the run
    call MPI_Finalize(ierr)

end program
