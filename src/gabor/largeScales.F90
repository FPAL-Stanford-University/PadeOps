module largeScalesMod
  use kind_parameters, only: rkind
  use domainSetup, only: gpLES, gpLESb, getStartAndEndIndices, decomp2Dpencil, &
    isBCperiodic, dxLES, dyLES, dzLES, gpQHcent
  use DerivativesMod, only: derivatives
  use decomp_2d
  use exits, only: gracefulExit, message
  use fortran_assert, only: assert
  use hdf5
  use hdf5_fcns, only: read_h5_chunk_data, write_h5_chunk_data, &
    createAndOpenFile, closeFileResources
  use mpi, only: MPI_COMM_WORLD

  implicit none

  real(rkind), dimension(:,:,:), allocatable :: U, V, W
  real(rkind), dimension(:,:,:,:,:), allocatable :: gradU
  real(rkind), dimension(:,:,:), allocatable :: Ubuff1, Ubuff2
  real(rkind), dimension(:,:,:), allocatable :: dUbuff1, dUbuff2
  type(derivatives) :: grad
  logical :: gradientPossible = .true.

  ! Large scale integral quantities
  real(rkind), dimension(:,:,:), allocatable :: KE, L
  real(rkind) :: cL = 1.d0
  
  ! Error checks
  logical :: initialized = .false., readUVW = .false., gradUcomputed = .false.
  logical :: computeLrgSclQOIs = .false. ! computed large scale (Q)uantities
                                         ! (O)f (I)nterest? See
                                         ! computeLargeScaleParams.F90 in
                                         ! problems/gabor/<current problem files>/
  
  interface readFields
    module procedure read3Fields3D, read2Fields3D, read1Field5D
  end interface
  interface writeFields
    module procedure write3Fields3D, write1Field3D, write1Field5D
  end interface

  contains
    include "largeScales_files/generalIO.F90"

    subroutine initLargeScales(inputfile,gradPossible)
      ! Initialize variables and allocate memory for the large scale field
      ! Inputs:
      !     gradPossible --> If we decide to read the gradient from disk
      !                      (instead of compute on the fly) this will 
      !                      be false. Default is true
      logical, intent(in), optional :: gradPossible
      character(len=*), intent(in) :: inputfile
      character(len=4) :: ddxMethod, ddyMethod, ddzMethod
      integer :: ist, ien, jst, jen, kst, ken
      integer :: isz, jsz, ksz
      integer :: ierr, ioUnit
      logical, dimension(3) :: periodic
      
      namelist /LARGESCALES/ cL

      ! Read cL from inputfile
      ioUnit = 1
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=ioUnit, NML=LARGESCALES)
      close(ioUnit)

      if (present(gradPossible)) gradientPossible = gradPossible

      if (gradientPossible) then
        call getStartAndEndIndices(gpLES,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
      else
        call getStartAndEndIndices(gpLESb,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
      end if
      allocate(U(isz,jsz,ksz),V(isz,jsz,ksz),W(isz,jsz,ksz))
      allocate(gradU(3,3,isz,jsz,ksz))

      if (gradientPossible) then
        select case (decomp2Dpencil)
          case('x')
            allocate(Ubuff1 (gpLES%ysz(1),gpLES%ysz(2),gpLES%ysz(3)))
            allocate(dUbuff1(gpLES%ysz(1),gpLES%ysz(2),gpLES%ysz(3)))
            
            allocate(Ubuff2 (gpLES%zsz(1),gpLES%zsz(2),gpLES%zsz(3)))
            allocate(dUbuff2(gpLES%zsz(1),gpLES%zsz(2),gpLES%zsz(3)))
          case('y')
            allocate(Ubuff1 (gpLES%xsz(1),gpLES%xsz(2),gpLES%xsz(3)))
            allocate(dUbuff1(gpLES%xsz(1),gpLES%xsz(2),gpLES%xsz(3)))
            
            allocate(Ubuff2 (gpLES%zsz(1),gpLES%zsz(2),gpLES%zsz(3)))
            allocate(dUbuff2(gpLES%zsz(1),gpLES%zsz(2),gpLES%zsz(3)))
          case('z')
            allocate(Ubuff2 (gpLES%xsz(1),gpLES%xsz(2),gpLES%xsz(3)))
            allocate(dUbuff2(gpLES%xsz(1),gpLES%xsz(2),gpLES%xsz(3)))
            
            allocate(Ubuff1 (gpLES%ysz(1),gpLES%ysz(2),gpLES%ysz(3)))
            allocate(dUbuff1(gpLES%ysz(1),gpLES%ysz(2),gpLES%ysz(3)))
          case default
            call gracefulExit("Must select 'x', 'y', or 'z' domain decomposition",ierr)
        end select

        ! Initialize the derivative class for computing velocity gradients
        ddxMethod = 'cd06'
        ddyMethod = 'cd06'
        ddzMethod = 'cd06'

        call isBCperiodic(periodic)
        if (periodic(1)) ddxMethod = 'four'
        if (periodic(2)) ddyMethod = 'four'
        if (periodic(3)) ddzMethod = 'four'

        call grad%init(gpLES,dxLES,dyLES,dzLES,periodic(1),periodic(2),periodic(3),&
          ddxMethod,ddyMethod,ddzMethod)

        ! Initialize memory for the large scale parameters
        call getStartAndEndIndices(gpQHcent,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
        allocate(KE(isz,jsz,ksz), L(isz,jsz,ksz))
        initialized = .true.
      end if 
    end subroutine

    subroutine finalizeLargeScales()
      if (allocated(U)) deallocate(U)
      if (allocated(V)) deallocate(V)
      if (allocated(W)) deallocate(W)
      if (allocated(gradU)) deallocate(gradU)
      if (allocated(Ubuff1)) deallocate(Ubuff1)
      if (allocated(Ubuff2)) deallocate(Ubuff2)
      if (allocated(dUbuff1)) deallocate(dUbuff1)
      if (allocated(dUbuff2)) deallocate(dUbuff2)
      if (gradientPossible) call grad%destroy()
      if (allocated(KE)) deallocate(KE)
      if (allocated(L)) deallocate(L)
    end subroutine

    subroutine getLargeScaleData(fname,useBoundaryData)
      ! This routine reads in the large scale velocity field
      ! Inputs:
      !     fname --> file name (including full path)
      !     useBoundaryData --> Ultimately we need large scale data at the
      !                         domain boundaries. Ideally, we use PadeOps
      !                         routines to interpolate/extrapolate this, but in
      !                         the interest of time we include the option to
      !                         simply read the data in. If this is true, we
      !                         cannot allow computLargeScaleGradient() below to
      !                         operate on this velocity since the grid is
      !                         nonuniform in z and includes the periodic plane
      !                         in x and y
      character(len=*), intent(in) :: fname
      logical, intent(in), optional :: useBoundaryData
      logical :: bdryData = .false.

      call assert(initialized,'Trying to read data without initializing large scales')
      if (present(useBoundaryData)) bdryData = useBoundaryData
      if (bdryData) then
        gradientPossible = .false.
        call readFields(fname,U,V,W,'/u','/v','/w',gpLESb)
      else
        call readFields(fname,U,V,W,'/u','/v','/w',gpLES)
      end if

      readUVW = .true.
    end subroutine

    subroutine computeLargeScaleGradient(fname)
      character(len=*), intent(in), optional :: fname
      integer :: ierr
      call assert(initialized,'Trying to read data without initializing large scales')
      call assert(readUVW,'Trying to compute gradient without velocity data')
      if (gradientPossible) then
        select case (decomp2Dpencil) 
          case('x')
            ! x derivatives
            call grad%ddx(U,gradU(1,1,:,:,:))
            call grad%ddx(V,gradU(2,1,:,:,:))
            call grad%ddx(W,gradU(3,1,:,:,:))
  
            ! y derivatives
            call transpose_x_to_y(U,Ubuff1,gpLES)
            call grad%ddy(Ubuff1,dUbuff1)
            call transpose_y_to_x(dUbuff1,gradU(1,2,:,:,:),gpLES)
            
            call transpose_x_to_y(V,Ubuff1,gpLES)
            call grad%ddy(Ubuff1,dUbuff1)
            call transpose_y_to_x(dUbuff1,gradU(2,2,:,:,:),gpLES)
            
            call transpose_x_to_y(W,Ubuff1,gpLES)
            call grad%ddy(Ubuff1,dUbuff1)
            call transpose_y_to_x(dUbuff1,gradU(3,2,:,:,:),gpLES)
  
            ! z derivatives
            call transpose_x_to_y(U,Ubuff1,gpLES)
            call transpose_y_to_z(Ubuff1,Ubuff2,gpLES)
            call grad%ddz(Ubuff2,dUbuff2)
            call transpose_z_to_y(dUbuff2,dUbuff1,gpLES)
            call transpose_y_to_x(dUbuff1,gradU(1,3,:,:,:),gpLES)

            call transpose_x_to_y(V,Ubuff1,gpLES)
            call transpose_y_to_z(Ubuff1,Ubuff2,gpLES)
            call grad%ddz(Ubuff2,dUbuff2)
            call transpose_z_to_y(dUbuff2,dUbuff1,gpLES)
            call transpose_y_to_x(dUbuff1,gradU(2,3,:,:,:),gpLES)

            call transpose_x_to_y(W,Ubuff1,gpLES)
            call transpose_y_to_z(Ubuff1,Ubuff2,gpLES)
            call grad%ddz(Ubuff2,dUbuff2)
            call transpose_z_to_y(dUbuff2,dUbuff1,gpLES)
            call transpose_y_to_x(dUbuff1,gradU(3,3,:,:,:),gpLES)
          case('y')
            ! x derivatives
            call transpose_y_to_x(U,Ubuff1,gpLES)
            call grad%ddx(Ubuff1,dUbuff1)
            call transpose_x_to_y(dUbuff1,gradU(1,1,:,:,:),gpLES)

            call transpose_y_to_x(V,Ubuff1,gpLES)
            call grad%ddx(Ubuff1,dUbuff1)
            call transpose_x_to_y(dUbuff1,gradU(2,1,:,:,:),gpLES)

            call transpose_y_to_x(W,Ubuff1,gpLES)
            call grad%ddx(Ubuff1,dUbuff1)
            call transpose_x_to_y(dUbuff1,gradU(3,1,:,:,:),gpLES)

            ! y derivatives
            call grad%ddy(U,gradU(1,2,:,:,:))
            call grad%ddy(V,gradU(2,2,:,:,:))
            call grad%ddy(W,gradU(3,2,:,:,:))
           
            ! z derivatives
            call transpose_y_to_z(U,Ubuff2,gpLES)
            call grad%ddz(Ubuff2,dUbuff2)
            call transpose_z_to_y(dUbuff2,gradU(1,3,:,:,:),gpLES)

            call transpose_y_to_z(V,Ubuff2,gpLES)
            call grad%ddz(Ubuff2,dUbuff2)
            call transpose_z_to_y(dUbuff2,gradU(2,3,:,:,:),gpLES)

            call transpose_y_to_z(W,Ubuff2,gpLES)
            call grad%ddz(Ubuff2,dUbuff2)
            call transpose_z_to_y(dUbuff2,gradU(3,3,:,:,:),gpLES)
          case('z')
            ! x derivatives
            call transpose_z_to_y(U,Ubuff1,gpLES)
            call transpose_y_to_x(Ubuff1,Ubuff2,gpLES)
            call grad%ddx(Ubuff2,dUbuff2)
            call transpose_x_to_y(dUbuff2,dUbuff1,gpLES)
            call transpose_y_to_z(dUbuff1,gradU(1,1,:,:,:),gpLES)

            call transpose_z_to_y(V,Ubuff1,gpLES)
            call transpose_y_to_x(Ubuff1,Ubuff2,gpLES)
            call grad%ddx(Ubuff2,dUbuff2)
            call transpose_x_to_y(dUbuff2,dUbuff1,gpLES)
            call transpose_y_to_z(dUbuff1,gradU(2,1,:,:,:),gpLES)

            call transpose_z_to_y(W,Ubuff1,gpLES)
            call transpose_y_to_x(Ubuff1,Ubuff2,gpLES)
            call grad%ddx(Ubuff2,dUbuff2)
            call transpose_x_to_y(dUbuff2,dUbuff1,gpLES)
            call transpose_y_to_z(dUbuff1,gradU(3,1,:,:,:),gpLES)
  
            ! y derivatives
            call transpose_z_to_y(U,Ubuff1,gpLES)
            call grad%ddy(Ubuff1,dUbuff1)
            call transpose_y_to_z(dUbuff1,gradU(1,2,:,:,:),gpLES)
            
            call transpose_z_to_y(V,Ubuff1,gpLES)
            call grad%ddy(Ubuff1,dUbuff1)
            call transpose_y_to_z(dUbuff1,gradU(2,2,:,:,:),gpLES)
            
            call transpose_z_to_y(W,Ubuff1,gpLES)
            call grad%ddy(Ubuff1,dUbuff1)
            call transpose_y_to_z(dUbuff1,gradU(3,2,:,:,:),gpLES)
  
            ! z derivatives
            call grad%ddz(U,gradU(1,3,:,:,:))
            call grad%ddz(V,gradU(2,3,:,:,:))
            call grad%ddz(W,gradU(3,3,:,:,:))
          case default
            call gracefulExit("ERROR: must set decomp2Dpencil to 'x', 'y', or 'z'", ierr)
        end select
      else
        call readFields(fname,gradU,'/gradU',gpLESb)
      end if
      gradUcomputed = .true.
    end subroutine
end module
