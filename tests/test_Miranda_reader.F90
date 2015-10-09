program test_Miranda_reader
    use mpi
    use kind_parameters, only : clen
    use miranda_tools, only: miranda_reader
    use io_VTK_stuff, only: io_VTK

    implicit none

    type(miranda_reader) :: mir
    type(io_VTK)         :: viz

    character(len=clen) :: jobdir = '/home/akshays/Data/miranda/BW/IRM_30MODES_128_CORRECTED'
    integer :: prow = 8, pcol = 1

    character(len=clen), dimension(:), allocatable :: varnames
    logical :: writeviz = .TRUE.

    integer :: ierr, step

    call MPI_Init(ierr)
    
    ! Initialize miranda_reader object
    call mir%init(jobdir, prow, pcol)

    ! Read in the grid
    call mir%read_grid()

    if ( writeviz ) then
        allocate( varnames(mir%nvars+mir%ns) )
        varnames = ['u   ','v   ','w   ','rho ','e   ','p   ','T   ','sos ','mu  ','bulk','ktc ','Diff','N2  ','CO2 ']
        call viz%init('.', 'MIR', mir%nvars+mir%ns, varnames)
    end if

    ! Read in data for time step 0
    step = 0
    call mir%read_data(step)

    if ( writeviz ) then
        call viz%WriteViz(mir%gp, mir%mesh, mir%fields)
    end if

    call mir%destroy()
    if ( writeviz ) then
        deallocate( varnames )
        call viz%destroy()
    end if
    call MPI_Finalize(ierr)

end program 
