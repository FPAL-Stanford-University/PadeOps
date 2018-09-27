#include "inclinedRM_TKEBudget_files/hooks.F90"

program inclinedRM_TKEBudget

    use mpi
    use hdf5
    use kind_parameters,  only: rkind,clen
    use CompressibleGrid, only: cgrid
    use exits,            only: message, GracefulExit
    implicit none

    type(cgrid) :: cgp
    character(len=clen) :: inputfile, stability, arg
    integer :: ierr, error, nargs, i
    integer :: nrestarts, step
    integer, dimension(:), allocatable :: restart_steps
    double precision :: start_time, elapsed_time

    ! Start MPI
    call MPI_Init(ierr)

    ! Initialize the HDF5 library and Fortran interfaces
    call h5open_f(error)
    if (error /= 0) call GracefulExit("Could not initialize HDF5 library and Fortran interfaces.",7356)

    nargs = command_argument_count()
    if (nargs < 1) call GracefulExit("Usage: mpiexec -n <nprocs> ./inclinedRM_TKEBudget <inputfile> [restart steps to process (ints)]",46)

    ! Get file location 
    call get_command_argument(1,inputfile)
    
    ! Initialize the grid object
    call cgp%init(inputfile)

    ! Setup postprocessing to be able to read restart files
    call cgp%setup_postprocessing(nrestarts)

    if (nargs > 1) then
        allocate(restart_steps(nargs-1))
        do i = 1,nargs-1
            call get_command_argument(i+1, arg)
            read(arg,*) restart_steps(i)
        end do
    else
        allocate(restart_steps(nrestarts+1))
        do i = 0,nrestarts
            restart_steps(i+1) = i
        end do
    end if

    do i = 1, size(restart_steps, 1)
        step = restart_steps(i)

        call message("Reading restart step", step)

        start_time = MPI_WTIME()

        ! Read restart file
        call cgp%read_restart(step)

        ! Time advance for one step
        call cgp%get_dt(stability)
        call cgp%advance_RK45(.true.)

        elapsed_time = MPI_WTIME() - start_time
        call message("Time to process step in seconds", elapsed_time)
    end do

    deallocate(restart_steps)
        
    ! Destroy everythin before ending
    call cgp%destroy_grid()

    ! Close Fortran interfaces and HDF5 library
    call h5close_f(error)

    ! End the run
    call MPI_Finalize(ierr)

end program
