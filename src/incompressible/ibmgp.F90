module ibmgpmod
    use kind_parameters, only: rkind, clen

    implicit none

    private
    public :: ibmgp

    type :: ibmgp
        private 
        !class(decomp_info), pointer :: gpC, gpE
        !class(spectral), pointer :: spectC, spectE

        integer :: num_surfelem
        real(rkind), allocatable, dimension(:,:,:) :: surfelem

        contains 
            !! ALL INIT PROCEDURES
            procedure          :: init
            procedure          :: destroy
            procedure, private :: dosomething
    end type 

contains

subroutine init(this, inputDir, inputFile)
  class(ibmgp),     intent(inout) :: this
  character(len=*), intent(in)    :: inputFile, inputDir

  character(len=clen)    :: surfaceMeshFile, fname, dumstr
  integer :: io, ii, jj, ioUnit, nlines
  real(rkind) :: rnum1, rnum2, rnum3

  namelist /IBMGP/ surfaceMeshFile

  ioUnit = 11
  open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
  read(unit=ioUnit, NML=IBMGP)
  close(ioUnit)

  ! -------Read surfaceMeshFile (x, y, z) for triangles
  ! -- First count the number of elements ---
  fname = inputDir(:len_trim(inputDir))//"/"//trim(surfaceMeshFile)
  ioUnit = 11
  open(unit=ioUnit, file=fname, form='FORMATTED')
  nlines = 0
  do 
      read(ioUnit, *, iostat=io)
      if(io .ne. 0) exit
      nlines = nlines+1
  enddo
  nlines = nlines-2
  this%num_surfelem = nlines/7
  close(ioUnit)

  allocate(this%surfelem(this%num_surfelem, 3, 3))

  ! -- Now read in the vertex data ---
  ioUnit = 11
  open(unit=ioUnit, file=fname, form='FORMATTED')
  read(ioUnit, *)
  do ii = 1, this%num_surfelem 
      ! ignore two lines
      read(ioUnit, *);       read(ioUnit, *)

      ! read first vertex of element ii
      jj = 1
      read(ioUnit, *) dumstr, rnum1, rnum2, rnum3
      this%surfelem(ii,jj,1) = rnum1; this%surfelem(ii,jj,2) = rnum2; this%surfelem(ii,jj,3) = rnum3

      ! read second vertex of element ii
      jj = 2
      read(ioUnit, *) dumstr, rnum1, rnum2, rnum3
      this%surfelem(ii,jj,1) = rnum1; this%surfelem(ii,jj,2) = rnum2; this%surfelem(ii,jj,3) = rnum3

      ! read third vertex of element ii
      jj = 3
      read(ioUnit, *) dumstr, rnum1, rnum2, rnum3
      this%surfelem(ii,jj,1) = rnum1; this%surfelem(ii,jj,2) = rnum2; this%surfelem(ii,jj,3) = rnum3

      ! ignore two lines
      read(ioUnit, *);       read(ioUnit, *)
  enddo
  close(ioUnit)

end subroutine


subroutine destroy(this)
  class(ibmgp), intent(inout) :: this
end subroutine

subroutine dosomething(this)
  class(ibmgp), intent(inout) :: this
end subroutine

end module
