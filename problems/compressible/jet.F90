#include "jet_files/hooks.F90"

program jet

    use hdf5
    use kind_parameters,  only: rkind,clen
    use CompressibleGrid, only: cgrid
    use exits,            only: GracefulExit
    implicit none

    type(cgrid) :: cgp
    character(len=clen) :: inputfile
    integer :: ierr, error

    ! Start MPI
    call MPI_Init(ierr)

    ! Initialize the HDF5 library and Fortran interfaces
    call h5open_f(error)
    if (error /= 0) call GracefulExit("Could not initialize HDF5 library and Fortran interfaces.",7356)

    ! Get file location 
    call GETARG(1,inputfile)
    
    ! Initialize the grid object
    call cgp%init(inputfile)

    ! Time advance
    call cgp%simulate()
        
    ! Destroy everythin before ending
    call cgp%destroy_grid()

    ! Close Fortran interfaces and HDF5 library
    call h5close_f(error)

    ! End the run
    call MPI_Finalize(ierr)

end program
