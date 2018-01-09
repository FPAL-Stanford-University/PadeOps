program test_Miranda_restart
    use mpi
    use kind_parameters,     only: rkind, clen
    use decomp_2d,           only: decomp_info, decomp_2d_init, get_decomp_info
    use miranda_restart_mod, only: miranda_restart
    use io_VTK_stuff,        only: io_VTK
    use exits,               only: message, GracefulExit

    implicit none

    type(decomp_info)     :: gp
    type(miranda_restart) :: mir
    type(io_VTK)          :: viz

    character(len=clen) :: jobdir, resfile, dummy
    integer :: prow = 0, pcol = 0
    integer :: nx, ny, nz

    character(len=clen), dimension(:), allocatable :: varnames
    logical :: writeviz = .TRUE.
    logical :: periodicx = .false., periodicy = .false., periodicz = .false.

    real(rkind), dimension(:,:,:,:), allocatable :: mesh, resdata
    real(rkind) :: tsim, dt

    integer :: ierr, step, i

    call MPI_Init(ierr)

    if( iargc() .LT. 6 ) then
        call GracefulExit("Usage: "//NEW_LINE('A')//"    mpiexec -n 8 ./test_Miranda_restart <jobdir> <resfile> <nx> <ny> <nz> <step>", 1729)
    end if

    ! Get command line arguments
    call getarg(1,jobdir)
    call getarg(2,resfile)
    call getarg(3,dummy); read(dummy, '(I)') nx
    call getarg(4,dummy); read(dummy, '(I)') ny
    call getarg(5,dummy); read(dummy, '(I)') nz
    call getarg(6,dummy); read(dummy, '(I)') step

    ! Initialize the grid partition object
    call decomp_2d_init(nx, ny, nz, prow, pcol, [periodicx, periodicy, periodicz])
    call get_decomp_info(gp)

    call message("Jobdir is "//trim(jobdir))
    call message("Restart file is "//trim(resfile))
    call message("nx", nx)
    call message("ny", ny)
    call message("nz", nz)

    ! Initialize miranda_restart object
    call mir%init(gp, jobdir, resfile, prow, pcol, periodicx, periodicy, periodicz)

    allocate( mesh   (gp%ysz(1), gp%ysz(2), gp%ysz(3), 3       ) )
    allocate( resdata(gp%ysz(1), gp%ysz(2), gp%ysz(3), mir%nres) )

    ! Read in the grid
    call mir%read_grid(mesh)

    if ( writeviz ) then
        allocate( varnames(mir%nres) )
        varnames(1:5) = ['u', 'v', 'w', 'rho', 'e']
        do i = 1,mir%ns
            write(dummy,'(A,I2.2)') "Y_", i
            varnames(5+i) = trim(dummy)
        end do
        varnames(5+mir%ns+1:5+mir%ns+2) = ['p', 'T']
        call viz%init('.', 'mir_restart', mir%nres, varnames)
    end if

    ! Read in data for time step 0
    call mir%read_data(step, resdata, tsim, dt)
    
    call message("Step", step)
    call message("Simulation time at restart", tsim)
    call message("Last time step", dt)

    if ( writeviz ) then
        call viz%WriteViz(gp, mesh, resdata)
    end if

    deallocate( mesh    )
    deallocate( resdata )

    call mir%destroy()
    if ( writeviz ) then
        deallocate( varnames )
        call viz%destroy()
    end if
    call MPI_Finalize(ierr)

end program 
