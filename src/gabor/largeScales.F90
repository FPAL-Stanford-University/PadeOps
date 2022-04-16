module largeScalesMod
  use kind_parameters, only: rkind
  use domainSetup, only: gpLES, gpLESb, getStartAndEndIndices, decomp2Dpencil, &
    periodic, dxLES, dyLES, dzLES
  use DerivativesMod, only: derivatives
  use decomp_2d
  use exits, only: gracefulExit
  use gaborIO_mod, only: readFields, writeFields
  use fortran_assert, only: assert

  implicit none

  real(rkind), dimension(:,:,:), allocatable :: U, V, W
  real(rkind), dimension(:,:,:,:,:), allocatable :: gradU
  real(rkind), dimension(:,:,:), allocatable :: Ubuff1, Ubuff2
  real(rkind), dimension(:,:,:), allocatable :: dUbuff1, dUbuff2
  type(derivatives) :: grad
  logical :: gradientPossible = .true.
  contains
    subroutine initLargeScales(gradPossible)
      ! Initialize variables and allocate memory for the large scale field
      ! Inputs:
      !     gradPossible --> If we decide to read the gradient from disk
      !                      (instead of compute on the fly) this will 
      !                      be false. Default is true
      logical, intent(in), optional :: gradPossible
      character(len=4) :: ddxMethod, ddyMethod, ddzMethod
      integer :: ist, ien, jst, jen, kst, ken
      integer :: isz, jsz, ksz
      integer :: ierr

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
            call gracefulExit('y-decomposition not implemented',ierr)
          case('z')
            call gracefulExit('z-decomposition not implemented',ierr)
          case default
            call gracefulExit("Must select 'x', 'y', or 'z' domain decomposition",ierr)
        end select

        ! Initialize the derivative class for computing velocity gradients
        ddxMethod = 'cd06'
        ddyMethod = 'cd06'
        ddzMethod = 'cd06'
        if (periodic(1)) ddxMethod = 'four'
        if (periodic(2)) ddyMethod = 'four'
        if (periodic(3)) ddzMethod = 'four'

        call grad%init(gpLES,dxLES,dyLES,dzLES,periodic(1),periodic(2),periodic(3),&
          ddxMethod,ddyMethod,ddzMethod)
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

      if (present(useBoundaryData)) bdryData = useBoundaryData
      if (bdryData) then
        gradientPossible = .false.
        call readFields(fname,U,V,W,'/u','/v','/w',gpLESb)
      else
        call readFields(fname,U,V,W,'/u','/v','/w',gpLES)
      end if
    end subroutine

    subroutine computeLargeScaleGradient(fname)
      character(len=*), intent(in), optional :: fname
      integer :: ierr
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
          case default
            call gracefulExit("Code does not support y or z "//&
              "domain decompositions", ierr)
        end select
      else
        call readFields(fname,gradU,'/gradU',gpLESb)
      end if
    end subroutine
end module
