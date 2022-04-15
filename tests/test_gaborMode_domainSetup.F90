program test_gaborMode_domainSetup
  use mpi
  use decomp_2d
  use kind_parameters, only: rkind, clen
  use domainSetup, only: setupDomainXYperiodic, finalizeDomainSetup, xLESb, yLESb, zLESb, gpLESb, &
    xQHedge, yQHedge, zQHedge, gpQHcent, nxLES, nyLES, nzLES, nxF, nyF, nzF, nxQH, nyQH, nzQH, &
    xQHcent, yQHcent, zQHcent, xF, yF, zF, gpF, xFh, yFh, zFh, nxsupp, nysupp, nzsupp, &
    Lx, Ly, Lz, xLES, yLES, zLES, gpLES
  use fortran_assert, only: assert
  use basic_io, only: read_1d_ascii

  implicit none
  character(len=clen) :: inputfile, datadir, fname, mssg
  integer :: i, ist, ien, jst, jen, kst, ken
  integer :: istF, ienF, jstF, jenF, kstF, kenF
  integer :: istM, ienM, jstM, jenM, kstM, kenM
  integer :: ierr, ioUnit
  real(rkind), dimension(:), allocatable :: xLESbTrue, yLESbTrue, zLESbTrue
  real(rkind), dimension(:), allocatable :: xLEStrue, yLEStrue, zLEStrue
  real(rkind), dimension(:), allocatable :: xQHedgeTrue, yQHedgeTrue, zQHedgeTrue
  real(rkind), dimension(:), allocatable :: xQHcentTrue, yQHcentTrue, zQHcentTrue
  real(rkind), dimension(:), allocatable :: xFtrue, yFtrue, zFtrue
  real(rkind), dimension(:), allocatable :: xFpTrue, yFpTrue, zFbTrue
  real(rkind) :: small = 1.d-7
  integer :: nxLESperQH, nyLESperQH, nzLESperQH
  integer :: pcol, prow
  namelist /IO/ datadir
  
  ! Initialize MPI
  call MPI_Init(ierr)

  ! Get the input file path and file name
  call GETARG(1,inputfile)

  ! Setup the domain
  call setupDomainXYperiodic(inputfile)

  ! Get ground truth from MATLAB
    ! Read inputfile
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
!      read(unit=ioUnit, NML=DOMAIN)
      read(unit=ioUnit, NML=IO)
      close(ioUnit)
    
    ! Read LES mesh data
      allocate(xLESbTrue(nxLES+1),yLESbTrue(nyLES+1),zLESbTrue(nzLES+2))
      write(fname,'(A9)') 'xLESb.dat'
      call read_1d_ascii(xLESbTrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A9)') 'yLESb.dat'
      call read_1d_ascii(yLESbTrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A9)') 'zLESb.dat'
      call read_1d_ascii(zLESbTrue,trim(datadir)//'/'//trim(fname))
  
      allocate(xLEStrue(nxLES),yLEStrue(nyLES),zLEStrue(nzLES))
      write(fname,'(A8)') 'xLES.dat'
      call read_1d_ascii(xLEStrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A8)') 'yLES.dat'
      call read_1d_ascii(yLEStrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A8)') 'zLES.dat'
      call read_1d_ascii(zLEStrue,trim(datadir)//'/'//trim(fname))
  
    ! Read QHmesh data
      allocate(xQHcentTrue(nxQH),yQHcentTrue(nyQH),zQHcentTrue(nzQH))
      write(fname,'(A11)') 'xQHcent.dat'
      call read_1d_ascii(xQHcentTrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A11)') 'yQHcent.dat'
      call read_1d_ascii(yQHcentTrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A11)') 'zQHcent.dat'
      call read_1d_ascii(zQHcentTrue,trim(datadir)//'/'//trim(fname))

      allocate(xQHedgeTrue(nxQH+1),yQHedgeTrue(nyQH+1),zQHedgeTrue(nzQH+1))
      write(fname,'(A11)') 'xQHedge.dat'
      call read_1d_ascii(xQHedgeTrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A11)') 'yQHedge.dat'
      call read_1d_ascii(yQHedgeTrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A11)') 'zQHedge.dat'
      call read_1d_ascii(zQHedgeTrue,trim(datadir)//'/'//trim(fname))
       
    ! Read high resolution mesh data
      allocate(xFtrue(nxF),yFtrue(nyF),zFtrue(nzF))
      write(fname,'(A6)') 'xF.dat'
      call read_1d_ascii(xFtrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A6)') 'yF.dat'
      call read_1d_ascii(yFtrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A6)') 'zF.dat'
      call read_1d_ascii(zFtrue,trim(datadir)//'/'//trim(fname))
 
      allocate(xFpTrue(nxF+nxsupp),yFpTrue(nyF+nysupp),zFbTrue(nzF+1))
      write(fname,'(A7)') 'xFp.dat'
      call read_1d_ascii(xFpTrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A7)') 'yFp.dat'
      call read_1d_ascii(yFpTrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A7)') 'zFb.dat'
      call read_1d_ascii(zFbTrue,trim(datadir)//'/'//trim(fname))
 
  ! Ensure the values match the expected values
  do i = 1,nproc
    if (nrank == i - 1) then
      ! LES mesh:
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", xLESb discrepency. Max difference = ",&
        maxval(xLESb - xLESbTrue(gpLESb%xst(1):gpLESb%xen(1)))
      call assert(maxval(xLESb - xLESbTrue(gpLESb%xst(1):gpLESb%xen(1))) < small,trim(mssg))
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", yLESb discrepency. Max difference = ",&
        maxval(yLESb - yLESbTrue(gpLESb%xst(2):gpLESb%xen(2)))
      call assert(maxval(yLESb - yLESbTrue(gpLESb%xst(2):gpLESb%xen(2))) < small,trim(mssg))
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", zLESb discrepency. Max difference = ",&
        maxval(zLESb - zLESbTrue(gpLESb%xst(3):gpLESb%xen(3)))
      call assert(maxval(zLESb - zLESbTrue(gpLESb%xst(3):gpLESb%xen(3))) < small,trim(mssg))

      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", xLES discrepency. Max difference = ",&
        maxval(xLES - xLEStrue(gpLES%xst(1):gpLES%xen(1)))
      call assert(maxval(xLES - xLEStrue(gpLES%xst(1):gpLES%xen(1))) < small,trim(mssg))
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", yLES discrepency. Max difference = ",&
        maxval(yLES - yLEStrue(gpLES%xst(2):gpLES%xen(2)))
      call assert(maxval(yLES - yLEStrue(gpLES%xst(2):gpLES%xen(2))) < small,trim(mssg))
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", zLES discrepency. Max difference = ",&
        maxval(zLES - zLEStrue(gpLES%xst(3):gpLES%xen(3)))
      call assert(maxval(zLES - zLEStrue(gpLES%xst(3):gpLES%xen(3))) < small,trim(mssg))

      ! QH edges:
      ist = gpQHcent%xst(1); ien = gpQHcent%xen(1)+1
      jst = gpQHcent%xst(2); jen = gpQHcent%xen(2)+1
      kst = gpQHcent%xst(3); ken = gpQHcent%xen(3)+1
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", xQHedge discrepency. Max difference = ",&
        maxval(xQHedge - xQHedgeTrue(ist:ien))
      call assert(maxval(xQHedge - xQHedgeTrue(ist:ien)) < small,trim(mssg))
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", yQHedge discrepency. Max difference = ",&
        maxval(yQHedge - yQHedgeTrue(jst:jen))
      call assert(maxval(yQHedge - xQHedgeTrue(jst:jen)) < small,trim(mssg))
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", zQHedge discrepency. Max difference = ",&
        maxval(zQHedge - zQHedgeTrue(kst:ken))
      call assert(maxval(zQHedge - xQHedgeTrue(kst:ken)) < small,trim(mssg))
      
      ! QH centers:
      ist = gpQHcent%xst(1); ien = gpQHcent%xen(1)
      jst = gpQHcent%xst(2); jen = gpQHcent%xen(2)
      kst = gpQHcent%xst(3); ken = gpQHcent%xen(3)
      call assert(size(xQHcent) == size(xQHcentTrue(ist:ien)),'xQHcent size comparison',nrank)
      call assert(size(yQHcent) == size(yQHcentTrue(jst:jen)),'yQHcent size comparison',nrank)
      call assert(size(zQHcent) == size(zQHcentTrue(kst:ken)),'zQHcent size comparison',nrank)
      write(mssg,'(A,F10.8)') "xQHcent discrepency. Max difference = ",&
        maxval(xQHcent - xQHcentTrue(ist:ien))
      call assert(maxval(xQHcent - xQHcentTrue(ist:ien)) < small,trim(mssg),nrank)
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", yQHcent discrepency. Max difference = ",&
        maxval(yQHcent - yQHcentTrue(jst:jen))
      call assert(maxval(yQHcent - yQHcentTrue(jst:jen)) < small,trim(mssg))
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", zQHcent discrepency. Max difference = ",&
        maxval(zQHcent - zQHcentTrue(kst:ken))
      call assert(maxval(zQHcent - zQHcentTrue(kst:ken)) < small,trim(mssg))

      ! High resolution mesh:
      ist = gpF%xst(1); ien = gpF%xen(1)
      jst = gpF%xst(2); jen = gpF%xen(2)
      kst = gpF%xst(3); ken = gpF%xen(3)
      write(mssg,'(A,F10.8)') "xF discrepency. Max difference = ",&
        maxval(xF - xFtrue(ist:ien))
      call assert(maxval(xF - xFtrue(ist:ien)) < small,trim(mssg),nrank)
      write(mssg,'(A,F10.8)') "yF discrepency. Max difference = ",&
        maxval(yF - yFtrue(jst:jen))
      call assert(maxval(yF - yFtrue(jst:jen)) < small,trim(mssg),nrank)
      write(mssg,'(A,F10.8)') "zF discrepency. Max difference = ",&
        maxval(zF - zFtrue(kst:ken))
      call assert(maxval(zF - zFtrue(kst:ken)) < small,trim(mssg),nrank)
      
      istF = ist-nxsupp/2; ienF = ien+nzsupp/2 
      jstF = jst-nysupp/2; jenF = jen+nysupp/2 
      kst = max(1,kst-nzsupp/2); ken = min(nzF+1,ken+nzsupp/2)

      istM = istF + nxsupp/2; ienM = ienF + nxsupp/2 
      jstM = jstF + nysupp/2; jenM = jenF + nysupp/2 
      
      write(mssg,'(A,F10.8)') "xFh discrepency. Max difference = ",&
        maxval(xFh - xFpTrue(istM:ienM))
      call assert(maxval(xFh - xFpTrue(istM:ienM)) < small,trim(mssg),nrank)
      write(mssg,'(A,F10.8)') "yFh discrepency. Max difference = ",&
        maxval(yFh - yFpTrue(jstM:jenM))
      call assert(maxval(yFh - yFpTrue(jstM:jenM)) < small,trim(mssg),nrank)
      write(mssg,'(A,F10.8)') "zFh discrepency. Max difference = ",&
        maxval(zFh - zFbTrue(kst:ken))
      call assert(maxval(zFh - zFbTrue(kst:ken)) < small,trim(mssg),nrank)

    end if
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
  end do
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (nrank == 0) then
    print*, "TEST PASSED!"
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! deallocate memory used by domain setup
  deallocate(xQHedgeTrue,yQHedgeTrue,zQHedgeTrue)
  deallocate(xQHcentTrue,yQHcentTrue,zQHcentTrue)
  deallocate(xLESbTrue,yLESbTrue,zLESbTrue)
  call finalizeDomainSetup()

  ! Finalize MPI
  call MPI_Finalize(ierr)  
end program
