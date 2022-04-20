module turbHalfChanMod
  use fortran_assert, only: assert
  use kind_parameters, only: rkind, clen
  implicit none
  contains
    subroutine initializeProblem(inputfile)
      use domainSetup, only: setupDomainXYperiodic
      use largeScalesMod, only: initLargeScales
      use GaborModeRoutines, only: initializeModes
      implicit none
      character(len=clen), intent(inout) :: inputfile
      call initializeExternalLibraries()
      call GETARG(1,inputfile)
      call setupDomainXYperiodic(inputfile)
      call initLargeScales(inputfile)
      call initializeModes(inputfile)
    end subroutine
    
    subroutine finalizeProblem()
      use domainSetup, only: finalizeDomainSetup
      use largeScalesMod, only: finalizeLargeScales
      use GaborModeRoutines, only: finalizeGaborModes
      implicit none
      call finalizeGaborModes()
      call finalizeLargeScales()
      call finalizeDomainSetup()
      call finalizeExternalLibraries()
    end subroutine
    
    subroutine initializeExternalLibraries()
      use mpi
      use hdf5
      implicit none
      integer :: ierr
    
      call MPI_Init(ierr)
      call assert(ierr == MPI_SUCCESS,'MPI initialization failed.')
      call H5open_f(ierr)
      call assert(ierr == 0,'HDF5 initialization failed.')
    end subroutine
    
    subroutine finalizeExternalLibraries()
      use mpi
      use hdf5
      implicit none
      integer :: ierr
    
      call H5close_f(ierr)
      call assert(ierr == 0,'HDF5 finalization failed.')
      call MPI_Finalize(ierr)
      call assert(ierr == MPI_SUCCESS,'MPI finalization failed.')
    end subroutine
   
    subroutine computeLargeScaleParams(fname,KEname,Lname)
      use largeScalesMod, only: KE, L, readUVW, computeLrgSclQOIs, cL
      use domainSetup, only: gpQHcent
      use gaborIO_mod, only: readFields
      use exits, only: gracefulExit
      implicit none
      character(len=*), intent(in), optional :: fname, KEname, Lname
      integer :: ierr

      call assert(readUVW,'Trying to compute large scale parameters without velocity data')
    
      if (present(fname)) then
        call assert(present(KEname),'Must specify dataset name for KE')
        call assert(present(Lname), 'Must specify dataset name for L')
        call readFields(fname,KE,L,KEname,Lname,gpQHcent)
        L = cL*L
      else
        call gracefulExit("Have not yet implemented routine to compute large "//&
          "scale quantities (i.e. kinetic energy and length scale). Must provide "//&
          "these in an HDF5 file.",ierr)
      end if
    
      computeLrgSclQOIs = .true.
    end subroutine
end module
