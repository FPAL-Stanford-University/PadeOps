program IRM_vorticity
    use mpi
    use kind_parameters, only: rkind, clen
    use miranda_tools,   only: miranda_reader
    use io_VTK_stuff,    only: io_VTK
    use DerivativesMod,  only: derivatives
    use gridtools,       only: alloc_buffs
    use operators,       only: curl
    use exits,           only: message, GracefulExit

    implicit none

    type(miranda_reader) :: mir
    type(io_VTK)         :: viz
    type(derivatives)    :: der

    character(len=clen) :: inputfile
    integer             :: iounit = 67

    character(len=clen) :: inputdir, outputdir
    integer :: prow, pcol
    logical :: periodicx = .FALSE., periodicy = .FALSE., periodicz = .FALSE.
    character(len=4) :: derivative_x = 'cd10', derivative_y = 'cd10', derivative_z = 'cd10'

    character(len=clen), dimension(:), allocatable :: varnames
    logical :: writeviz

    real(rkind), dimension(:,:,:,:), allocatable, target :: buffer
    real(rkind), dimension(:,:,:,:), pointer :: vort

    integer :: ierr, step

    namelist /INPUT/ inputdir, outputdir, &
                     periodicx, periodicy, periodicz, &
                     derivative_x, derivative_y, derivative_z, &
                     prow, pcol, writeviz

    call MPI_Init(ierr)

    if( command_argument_count() .LT. 1 ) then
        call GracefulExit("Usage: "//NEW_LINE('A')//"    mpiexec -n 8 ./test_Miranda_reader <input file>", 1729)
    end if

    call get_command_argument(1,inputfile)
    open(unit=iounit, file=trim(inputfile), form='FORMATTED')
    read(unit=iounit, NML=INPUT)
    close(iounit)

    call message("Jobdir is " // adjustl(trim(inputdir)) )
    
    ! Initialize miranda_reader object
    call mir%init(inputdir, prow, pcol, periodicx, periodicy, periodicz)

    ! Read in the grid
    call mir%read_grid()
    
    ! Initialize the derivative routine
    call der%init(                                      mir%gp, &
                          mir%dx,        mir%dy,        mir%dz, &
                   mir%periodicx, mir%periodicy, mir%periodicz, &
                    derivative_x,  derivative_y,  derivative_z  )

    ! Allocate buffer to store vorticity in Y decomposition
    call alloc_buffs(buffer, 3, 'y', mir%gp)
    vort => buffer(:,:,:,1:3)

    ! Initialize visualization stuff
    if ( writeviz ) then
        allocate( varnames(3) )
        varnames = ['X-vorticity','Y-vorticity','Z-vorticity']
        call viz%init(outputdir, 'vorticity', 3, varnames)
    end if

    ! Read in data for time step 0
    do step = 0,mir%nsteps-1
        call mir%read_data(step)
        call curl( mir%gp, der, mir%u, mir%v, mir%w, vort)

        if ( writeviz ) then
            ! Set vizcount to be same as Miranda step
            call viz%SetVizcount(step)
            call viz%WriteViz(mir%gp, mir%mesh, vort)
        end if
    end do

    ! Destroy all variables and exit cleanly
    if (allocated(buffer)) deallocate( buffer )
    call der%destroy()
    call mir%destroy()
    if ( writeviz ) then
        deallocate( varnames )
        call viz%destroy()
    end if
    call MPI_Finalize(ierr)

end program 
