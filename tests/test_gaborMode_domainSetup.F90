program test_gaborMode_domainSetup
  use mpi
  use decomp_2d
  use kind_parameters, only: rkind, clen
  use domainSetup, only: setupDomainXYperiodic, finalizeDomainSetup, xLESe, yLESe, zLESe, gpLESe, &
    xQHedge, yQHedge, zQHedge, gpQHcent, nxLES, nyLES, nzLES, nxF, nyF, nzF, nxQH, nyQH, nzQH, &
    xQHcent, yQHcent, zQHcent, xF, yF, zFC, gpFC, xFh, yFh, zFh, nxsupp, nysupp, nzsupp, &
    Lx, Ly, Lz, xLESc, yLESc, zLESc, gpLESc, kmin, kmax, nprocX, nprocY, nprocZ, &
    getStartAndEndIndices, nrankX, nrankY, nrankZ
  use fortran_assert, only: assert
  use basic_io, only: read_1d_ascii

  implicit none
  character(len=clen) :: inputfile, datadir, fname, mssg
  integer :: i, ist, ien, jst, jen, kst, ken
  integer :: isz, jsz, ksz
  integer :: istF, ienF, jstF, jenF
  integer :: istM, ienM, jstM, jenM
  integer :: ierr, ioUnit = 1
  real(rkind), dimension(:), allocatable :: xLESeTrue, yLESeTrue, zLESeTrue
  real(rkind), dimension(:), allocatable :: xLESctrue, yLESctrue, zLESctrue
  real(rkind), dimension(:), allocatable :: xQHedgeTrue, yQHedgeTrue, zQHedgeTrue
  real(rkind), dimension(:), allocatable :: xQHcentTrue, yQHcentTrue, zQHcentTrue
  real(rkind), dimension(:), allocatable :: xFtrue, yFtrue, zFtrue
  real(rkind), dimension(:), allocatable :: xFpTrue, yFpTrue, zFbTrue
  real(rkind) :: small = 1.d-7
  real(rkind) :: kminTrue = 8.d0, kmaxTrue = 32.d0
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
      allocate(xLESeTrue(nxLES+1),yLESeTrue(nyLES+1),zLESeTrue(nzLES+1))
      write(fname,'(A9)') 'xLESe.dat'
      call read_1d_ascii(xLESeTrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A9)') 'yLESe.dat'
      call read_1d_ascii(yLESeTrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A9)') 'zLESe.dat'
      call read_1d_ascii(zLESeTrue,trim(datadir)//'/'//trim(fname))
  
      allocate(xLESctrue(nxLES),yLESctrue(nyLES),zLESctrue(nzLES))
      write(fname,'(A8)') 'xLES.dat'
      call read_1d_ascii(xLESctrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A8)') 'yLES.dat'
      call read_1d_ascii(yLESctrue,trim(datadir)//'/'//trim(fname))
      write(fname,'(A8)') 'zLES.dat'
      call read_1d_ascii(zLESctrue,trim(datadir)//'/'//trim(fname))
  
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
      call getStartAndEndIndices(gpLESe,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", xLESe discrepency. Max difference = ",&
        maxval(xLESe - xLESeTrue(ist:ien))
      call assert(maxval(xLESe - xLESeTrue(ist:ien)) < small,trim(mssg))
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", yLESe discrepency. Max difference = ",&
        maxval(yLESe - yLESeTrue(jst:jen))
      call assert(maxval(yLESe - yLESeTrue(jst:jen)) < small,trim(mssg))
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", zLESe discrepency. Max difference = ",&
        maxval(zLESe - zLESeTrue(kst:ken))
      call assert(maxval(zLESe - zLESeTrue(kst:ken)) < small,trim(mssg))

      call getStartAndEndIndices(gpLESc,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", xLESc discrepency. Max difference = ",&
        maxval(xLESc - xLESctrue(ist:ien))
      call assert(maxval(xLESc - xLESctrue(ist:ien)) < small,trim(mssg))
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", yLESc discrepency. Max difference = ",&
        maxval(yLESc - yLESctrue(jst:jen))
      call assert(maxval(yLESc - yLESctrue(jst:jen)) < small,trim(mssg))
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", zLESc discrepency. Max difference = ",&
        maxval(zLESc - zLESctrue(kst:ken))
      call assert(maxval(zLESc - zLESctrue(kst:ken)) < small,trim(mssg))
      
      ! QH centers:
      call getStartAndEndIndices(gpQHcent,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
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

      ! QH edges:
      ien = ien + 1
      jen = jen + 1
      ken = ken + 1
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", xQHedge discrepency. Max difference = ",&
        maxval(xQHedge - xQHedgeTrue(ist:ien))
      call assert(maxval(xQHedge - xQHedgeTrue(ist:ien)) < small,trim(mssg))
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", yQHedge discrepency. Max difference = ",&
        maxval(yQHedge - yQHedgeTrue(jst:jen))
      call assert(maxval(yQHedge - xQHedgeTrue(jst:jen)) < small,trim(mssg))
      write(mssg,'(A,I2,A,F10.8)') "rank ", nrank, ", zQHedge discrepency. Max difference = ",&
        maxval(zQHedge - zQHedgeTrue(kst:ken))
      call assert(maxval(zQHedge - xQHedgeTrue(kst:ken)) < small,trim(mssg))

      ! High resolution mesh:
      call getStartAndEndIndices(gpFC,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
      write(mssg,'(A,F10.8)') "xF discrepency. Max difference = ",&
        maxval(xF - xFtrue(ist:ien))
      call assert(maxval(xF - xFtrue(ist:ien)) < small,trim(mssg),nrank)
      write(mssg,'(A,F10.8)') "yF discrepency. Max difference = ",&
        maxval(yF - yFtrue(jst:jen))
      call assert(maxval(yF - yFtrue(jst:jen)) < small,trim(mssg),nrank)
      write(mssg,'(A,F10.8)') "zF discrepency. Max difference = ",&
        maxval(zFC - zFtrue(kst:ken))
      call assert(maxval(zFC - zFtrue(kst:ken)) < small,trim(mssg),nrank)
      
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

      ! Confirm kmin and kmax
      call assert(abs(kmin - kminTrue) < small,'kmin discrepency')
      call assert(abs(kmax - kmaxTrue) < small,'kmax discrepency')
    end if
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
  end do
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (nrank == 0) then
    print*, "Number of partitions in x:", nprocX
    print*, "Number of partitions in y:", nprocY
    print*, "Number of partitions in z:", nprocZ
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  do i = 1,nproc
    if (nrank == i-1) then
      call getStartAndEndIndices(gpLESc,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
      print*, nrank, "ist: ", ist, "nrankX: ", nrankX
      print*, nrank, "jst: ", jst, "nrankY: ", nrankY
      print*, nrank, "kst: ", kst, "nrankZ: ", nrankZ
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
  deallocate(xLESeTrue,yLESeTrue,zLESeTrue)
  deallocate(xLEScTrue,yLEScTrue,zLEScTrue)
  call finalizeDomainSetup()

  ! Finalize MPI
  call MPI_Finalize(ierr)  
end program
