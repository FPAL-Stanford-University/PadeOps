#include "hitEnrich_files/initialize.F90"

program hitEnrich
  use kind_parameters,         only: clen, rkind
  use enrichmentMod,           only: enrichmentOperator, nthreads
  use IncompressibleGrid,      only: igrid
  use HIT_Periodic_parameters, only: Lx, Ly, Lz
  use auxiliary_openmp_subs,   only: GetArguments
  use mpi
  use decomp_2d,               only: nrank
  use fortran_assert,          only: assert
  use exits,                   only: message
  implicit none

  character(len=clen) :: inputfileLS, inputfileSS, inputfileGM 
  character(len=clen) :: datadir, fname, outputdir
  integer :: ioUnit, ierr, provided
  type(igrid) :: largeScales, smallScales
  type(enrichmentOperator) :: enrich
  real(rkind) :: tol = 1.d-13
  logical :: doQHtests = .false.
      
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  
  call GETARG(1,inputfileLS)
  call GETARG(2,inputfileSS)
  call GETARG(3,inputfileGM)
  call GetArguments(nthreads)
  
  call largeScales%init(inputfileLS, .true.) 
  call smallScales%init(inputfileSS, .false.)

  call enrich%init(smallScales,largeScales,inputfileGM,Lx,Ly,Lz)

  call message(" ")
  call message( "Running tests ...")
  call message( "                 ")
  if (doQHtests) then
    call message( "QH mesh tests ...")
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    ! QH-region checks
    call assert(abs(enrich%QHgrid%dx - 1.570796326794897d0) < tol)
    call assert(abs(enrich%QHgrid%dy - 1.570796326794897d0) < tol)
    call assert(abs(enrich%QHgrid%dz - 1.570796326794897d0) < tol)

    ! Visually confirm the QHmesh is set up correctly
    call message( " ")
    call message( "QH edge mesh. Visually verify this is what you expect")
    call message( " ")
    call message( "           nrank |    enrich%QHgrid%xE: ")
    call message( "________________________________________")
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    print*, nrank, "xE:", enrich%QHgrid%xE
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call message( " ")
    call message( "           nrank |    enrich%QHgrid%yE: ")
    call message( "________________________________________")
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    print*, nrank, "yE:", enrich%QHgrid%yE
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call message( " ")
    call message( "           nrank |    enrich%QHgrid%zE: ")
    call message( "________________________________________")
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    print*, nrank, "zE:", enrich%QHgrid%zE
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call message( " ")
    call message( "QH center mesh. Visually verify this is what you expect")
    call message( " ")
    call message( "           nrank |    enrich%QHgrid%xC: ")
    call message( "________________________________________")
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    print*, nrank, "xC:", enrich%QHgrid%xC
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call message( " ")
    call message( "           nrank |    enrich%QHgrid%yC: ")
    call message( "________________________________________")
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    print*, nrank, "yC:", enrich%QHgrid%yC
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call message( " ")
    call message( "           nrank |    enrich%QHgrid%zC: ")
    call message( "________________________________________")
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    print*, nrank, "zC:", enrich%QHgrid%zC
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    ! Check support window size
    call assert(enrich%nxsupp == 32)
    call assert(enrich%nysupp == 32)
    call assert(enrich%nzsupp == 32)
    
    ! Check kmin and kmax
    call assert(enrich%kmin == 4)
    call assert(enrich%kmax == 32)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call message(" ")
    call message("Verify this is the total number of modes you expect:")
    print*, nrank, "modes: ", enrich%nmodes
  end if

  do while (enrich%continueSimulation())
    call enrich%updateLargeScales()
    call enrich%advanceTime()
    call enrich%wrapupTimeStep()
  end do

  call enrich%destroy()
  call largeScales%destroy()
  call smallScales%destroy()

  call MPI_Finalize(ierr)
end program hitEnrich
