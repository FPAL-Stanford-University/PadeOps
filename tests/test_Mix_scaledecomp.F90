program test_Mix_scaledecomp
    use mpi
    use hdf5
    use decomp_2d,             only: nrank, decomp_info
    use kind_parameters,       only: rkind, clen
    use constants,             only: zero
    use miranda_reader_mod,    only: miranda_reader
    use DerivativesMod,        only: derivatives
    use FiltersMod,            only: filters
    use exits,                 only: message, GracefulExit
    use RKCoeffs,              only: RK45_steps
    use ScaleDecompositionMod, only: scaleDecomposition
    use reductions,            only: P_MAXVAL, P_MINVAL

    implicit none

    type(miranda_reader)     :: mir
    type(scaleDecomposition) :: scaledecomp
    type(derivatives)        :: der
    type(filters)            :: gfil

    character(len=clen) :: jobdir, inputfile, arg
    integer, dimension(:), allocatable :: steps
    integer :: prow = 0, pcol = 0
    integer :: nargs

    integer, dimension(2) :: x_bc = [0, 1]
    integer, dimension(2) :: y_bc = [0, 0]
    integer, dimension(2) :: z_bc = [0, 0]

    logical :: periodicx = .false., periodicy = .true., periodicz = .true.

    real(rkind) :: dt = real(1.0D-4, rkind)

    real(rkind), dimension(:,:,:),     allocatable :: rhoPsi_SD_old        ! Temporary variable for large scale entropic mixing rate
    real(rkind), dimension(:,:,:,:),   allocatable :: rhoPsi_SD_prefilter  ! Temporary variable for large scale entropic mixing generation
    real(rkind), dimension(:,:,:,:),   allocatable :: rhoPsi_SD_postfilter ! Temporary variable for large scale entropic mixing generation

    integer :: ierr, error, step, i

    call MPI_Init(ierr)

    nargs = command_argument_count()

#ifndef __bgq__
    if( nargs < 2 ) then
        call GracefulExit("Usage: "//NEW_LINE('A')//"    mpiexec -n 8 ./test_Mix_scaledecomp <jobdir> <inputfile> [restart steps to process (ints)]", 1729)
    end if
#endif


    ! Initialize the HDF5 library and Fortran interfaces
    call h5open_f(error)
    if (error /= 0) call GracefulExit("Could not initialize HDF5 library and Fortran interfaces.",7356)

    call get_command_argument(1,jobdir)
    call get_command_argument(2,inputfile)

    ! Initialize miranda_reader object
    call mir%init(jobdir, prow, pcol, periodicx, periodicy, periodicz)

    if (nrank == 0) then
        print *, "prow = ", prow, ", pcol = ", pcol
        print *, "Jobdir is "//trim(jobdir)
    end if

    if (nargs > 2) then
        allocate(steps(nargs-2))
        do i = 1,nargs-2
            call get_command_argument(i+2, arg)
            read(arg,*) steps(i)
        end do
    else
        allocate(steps(mir%nsteps))
        do i = 0,mir%nsteps-1
            steps(i+1) = i
        end do
    end if

    ! Read in the grid
    call mir%read_grid()

    ! Initialize the derivative object
    call der%init(                            mir%gp, &
                        mir%dx,    mir%dy,    mir%dz, &
                     periodicx, periodicy, periodicz, &
                        "cd10",    "cd10",    "cd10" )

    ! Initialize the filter object
    call gfil%init(                               mir%gp, &
                periodicx,     periodicy,      periodicz, &
               "gaussian",    "gaussian",     "gaussian"  )      

    ! Initialize budget object
    call scaledecomp%init(mir%gp, der, gfil, mir%mesh, mir%dx, mir%dy, mir%dz, mir%ns, x_bc, y_bc, z_bc, inputfile)

    allocate( rhoPsi_SD_old       (mir%gp%ysz(1),mir%gp%ysz(2),mir%gp%ysz(3)) )
    allocate( rhoPsi_SD_prefilter (mir%gp%ysz(1),mir%gp%ysz(2),mir%gp%ysz(3),RK45_steps) )
    allocate( rhoPsi_SD_postfilter(mir%gp%ysz(1),mir%gp%ysz(2),mir%gp%ysz(3),RK45_steps) )

    rhoPsi_SD_old = zero
    rhoPsi_SD_prefilter = zero
    rhoPsi_SD_postfilter = zero

    ! Read in data for time step 0
    do i = 1, size(steps, 1)
        step = steps(i)

        call message("Processing step ", step)
        call mir%read_data(step)

        call message("Max Ys_3", P_MAXVAL(mir%Ys(:,:,:,3)))
        call message("Min Ys_3", P_MINVAL(mir%Ys(:,:,:,3)))
        call scaledecomp%mix_budget(mir%rho, mir%u, mir%v, mir%w, mir%Ys, mir%Diff, &
                                    rhoPsi_SD_old, rhoPsi_SD_prefilter, rhoPsi_SD_postfilter, step*dt, dt)
    end do

    if ( allocated( rhoPsi_SD_old        ) ) deallocate( rhoPsi_SD_old        )
    if ( allocated( rhoPsi_SD_prefilter  ) ) deallocate( rhoPsi_SD_prefilter  )
    if ( allocated( rhoPsi_SD_postfilter ) ) deallocate( rhoPsi_SD_postfilter )

    call mir%destroy()
    call gfil%destroy()
    call der%destroy()

    deallocate(steps)
        
    ! Close Fortran interfaces and HDF5 library
    call h5close_f(error)

    call MPI_Finalize(ierr)

end program 
