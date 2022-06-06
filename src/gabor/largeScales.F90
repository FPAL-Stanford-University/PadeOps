module largeScalesMod
  use kind_parameters, only: rkind
  use domainSetup, only: gpLESc, gpLESe, getStartAndEndIndices, decomp2Dpencil, &
    isBCperiodic, dxLES, dyLES, dzLES, gpQHcent
  use DerivativesMod, only: derivatives
  use decomp_2d
  use exits, only: gracefulExit, message
  use fortran_assert, only: assert
  use hdf5
  use hdf5_fcns, only: read_h5_chunk_data, write_h5_chunk_data, &
    createAndOpenFile, closeFileResources
  use mpi, only: MPI_COMM_WORLD
  use PadeDerOps, only: Pade6stagg

  implicit none

  real(rkind), dimension(:,:,:), allocatable :: U, V, W, Utmp, Vtmp, Wtmp
  real(rkind), dimension(:,:,:), allocatable :: Up, Vp, Wp ! Fluctuating quantities
  real(rkind), dimension(:,:,:,:,:), allocatable :: gradU
  real(rkind), dimension(:,:,:), allocatable :: S ! Frobenius norm of gradU

  ! Each of the above large scale quantities requires a halo-padded counterpart 
  ! so that when interpolating to mode locations we have access to neighboring 
  ! processes' large scale data
  real(rkind), dimension(:,:,:), allocatable :: Uh, Vh, Wh
  real(rkind), dimension(:,:,:), allocatable :: Uph, Vph, Wph
  real(rkind), dimension(:,:,:,:,:), allocatable :: gradUh
  real(rkind), dimension(:,:,:), allocatable :: gradUbuff
  real(rkind), dimension(:,:,:), allocatable :: Sh

!  real(rkind), dimension(:,:,:), allocatable :: Ubuff1, Ubuff2, Ubuff3, Ubuff4
  real(rkind), dimension(:,:,:), allocatable :: Ubuff2Yc, Ubuff2Ye, Ubuff2Zc, &
    Ubuff2Ze, Ubuff2Xe
  real(rkind), dimension(:,:,:), allocatable :: dUbuff1, dUbuff2
  type(derivatives) :: grad
  logical :: gradientPossible = .true.

  ! Large scale integral quantities
  real(rkind), dimension(:,:,:), allocatable :: KE, L
  real(rkind) :: cL = 1.d0
  
  ! Interpolator for LES field velocity field
  type(Pade6stagg) :: interpLES
  
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
      integer :: ierr, ioUnit, scheme
      logical, dimension(3) :: periodic
      
      namelist /LARGESCALES/ cL

      ! Read cL from inputfile
      ioUnit = 1
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=ioUnit, NML=LARGESCALES)
      close(ioUnit)

      if (present(gradPossible)) gradientPossible = gradPossible
      call assert(.not. gradientPossible,'Interpolation from cell to edge '//&
        'needs to be debugged -- largeScales.F90')

      if (gradientPossible) then
        call getStartAndEndIndices(gpLESc,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
        allocate(Utmp(isz,jsz,ksz),Vtmp(isz,jsz,ksz),Wtmp(isz,jsz,ksz))
      end if
      call getStartAndEndIndices(gpLESe,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
      allocate(U(isz,jsz,ksz),V(isz,jsz,ksz),W(isz,jsz,ksz))
      allocate(Up(isz,jsz,ksz),Vp(isz,jsz,ksz),Wp(isz,jsz,ksz))
      allocate(gradU(3,3,isz,jsz,ksz))
      allocate(S(isz,jsz,ksz))

      ! Buffer arrays for data transpositions
      allocate(Ubuff2Yc(gpLESc%ysz(1),gpLESc%ysz(2),gpLESc%ysz(3)))
      allocate(Ubuff2Ye(gpLESe%ysz(1),gpLESe%ysz(2),gpLESe%ysz(3)))

      allocate(Ubuff2Zc(gpLESc%zsz(1),gpLESc%zsz(2),gpLESc%zsz(3)))
      allocate(Ubuff2Ze(gpLESe%zsz(1),gpLESe%zsz(2),gpLESe%zsz(3)))

      allocate(Ubuff2Xe(gpLESe%xsz(1),gpLESe%xsz(2),gpLESe%xsz(3)))
      ! Additional buffer arrays for differentiations
      if (gradientPossible) then
        select case (decomp2Dpencil)
          case('x')
            !allocate(Ubuff1 (gpLESe%ysz(1),gpLESe%ysz(2),gpLESe%ysz(3)))
            allocate(dUbuff1(gpLESe%ysz(1),gpLESe%ysz(2),gpLESe%ysz(3)))
            !
            !allocate(Ubuff2 (gpLESe%zsz(1),gpLESe%zsz(2),gpLESe%zsz(3)))
            allocate(dUbuff2(gpLESe%zsz(1),gpLESe%zsz(2),gpLESe%zsz(3)))
          case('y')
            !allocate(Ubuff1 (gpLESe%xsz(1),gpLESe%xsz(2),gpLESe%xsz(3)))
            allocate(dUbuff1(gpLESe%xsz(1),gpLESe%xsz(2),gpLESe%xsz(3)))
            
            !allocate(Ubuff2 (gpLESe%zsz(1),gpLESe%zsz(2),gpLESe%zsz(3)))
            allocate(dUbuff2(gpLESe%zsz(1),gpLESe%zsz(2),gpLESe%zsz(3)))
          case('z')
            !allocate(Ubuff2 (gpLESe%xsz(1),gpLESe%xsz(2),gpLESe%xsz(3)))
            allocate(dUbuff2(gpLESe%xsz(1),gpLESe%xsz(2),gpLESe%xsz(3)))
            
            !allocate(Ubuff1 (gpLESe%ysz(1),gpLESe%ysz(2),gpLESe%ysz(3)))
            allocate(dUbuff1(gpLESe%ysz(1),gpLESe%ysz(2),gpLESe%ysz(3)))
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

        call grad%init(gpLESc,dxLES,dyLES,dzLES,periodic(1),periodic(2),periodic(3),&
          ddxMethod,ddyMethod,ddzMethod)
      end if 
      
      ! Initialize memory for the large scale parameters
      call getStartAndEndIndices(gpQHcent,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
      allocate(KE(isz,jsz,ksz), L(isz,jsz,ksz))
      
      ! Initialize inerpolator
      if (.not. periodic(3)) then
        scheme = 1 ! 0: fd02, 1: cd06, 2: fourierColl
      else
        scheme = 2
      end if
      call interpLES%init(gpLESc,gpLESc,gpLESe,gpLESe,dzLES,scheme,periodic(3))
      initialized = .true.
    end subroutine

    subroutine finalizeLargeScales()
      if (allocated(U)) deallocate(U)
      if (allocated(V)) deallocate(V)
      if (allocated(W)) deallocate(W)
      if (allocated(Utmp)) deallocate(Utmp)
      if (allocated(Vtmp)) deallocate(Vtmp)
      if (allocated(Wtmp)) deallocate(Wtmp)
      if (allocated(Up)) deallocate(Up)
      if (allocated(Vp)) deallocate(Vp)
      if (allocated(Wp)) deallocate(Wp)
      if (allocated(gradU)) deallocate(gradU)
      if (allocated(S)) deallocate(S)

      ! Halo-padded large scale arrays
      if (allocated(Uh)) deallocate(Uh)
      if (allocated(Vh)) deallocate(Vh)
      if (allocated(Wh)) deallocate(Wh)
      if (allocated(Uph)) deallocate(Uph)
      if (allocated(Vph)) deallocate(Vph)
      if (allocated(Wph)) deallocate(Wph)
      if (allocated(gradUh)) deallocate(gradUh)
      if (allocated(gradUbuff)) deallocate(gradUbuff)
      if (allocated(Sh)) deallocate(Sh)

      ! Memory buffers for data transpositions
!      if (allocated(Ubuff1)) deallocate(Ubuff1)
!      if (allocated(Ubuff2)) deallocate(Ubuff2)
      if (allocated(Ubuff2Ye)) deallocate(Ubuff2Ye)
      if (allocated(Ubuff2Yc)) deallocate(Ubuff2Yc)
      if (allocated(Ubuff2Ze)) deallocate(Ubuff2Ze)
      if (allocated(Ubuff2Zc)) deallocate(Ubuff2Zc)
      if (allocated(Ubuff2Xe)) deallocate(Ubuff2Xe)
      if (allocated(dUbuff1)) deallocate(dUbuff1)
      if (allocated(dUbuff2)) deallocate(dUbuff2)

      ! Gradient operator
      if (gradientPossible) call grad%destroy()

      ! Large scale computed quantities
      if (allocated(KE)) deallocate(KE)
      if (allocated(L)) deallocate(L)
    end subroutine

    subroutine getLargeScaleData(fname)
      ! This routine reads in the large scale velocity field
      ! Inputs:
      !     fname --> file name (including full path)
      
      character(len=*), intent(in) :: fname
      real(rkind), dimension(:), allocatable :: Ubar, Vbar, Wbar
      integer :: ist, ien, jst, jen, kst, ken, isz, jsz, ksz, k

      call assert(initialized,'Trying to read data without initializing large scales')
      if (gradientPossible) then
  
        call assert(allocated(Utmp),'allocated(Utmp)')
        call assert(allocated(Vtmp),'allocated(Vtmp)')
        call assert(allocated(Wtmp),'allocated(Wtmp)')
 
        call readFields(fname,Utmp,Vtmp,Wtmp,'/u','/v','/w',gpLESc)
        readUVW = .true.
        call computeLargeScaleGradient()

        ! Interpolate fields to 0:dzLES:Lz mesh
        select case(decomp2Dpencil)
        case ('x')
          call transpose_x_to_y(Utmp,Ubuff2Yc,gpLESc)
          call transpose_y_to_z(Ubuff2Yc,Ubuff2Zc,gpLESc)
          call interpLES%interpz_C2E(Ubuff2Zc,Ubuff2Ze,-1,1)
          call transpose_z_to_y(Ubuff2Ze,Ubuff2Ye,gpLESe)
          call transpose_y_to_x(Ubuff2Ye,U,gpLESe)

          call transpose_x_to_y(Vtmp,Ubuff2Yc,gpLESc)
          call transpose_y_to_z(Ubuff2Yc,Ubuff2Zc,gpLESc)
          call interpLES%interpz_C2E(Ubuff2Zc,Ubuff2Ze,-1,1)
          call transpose_z_to_y(Ubuff2Ze,Ubuff2Ye,gpLESe)
          call transpose_y_to_x(Ubuff2Ye,V,gpLESe)
          
          call transpose_x_to_y(Wtmp,Ubuff2Yc,gpLESc)
          call transpose_y_to_z(Ubuff2Yc,Ubuff2Zc,gpLESc)
          call interpLES%interpz_C2E(Ubuff2Zc,Ubuff2Ze,-1,-1)
          call transpose_z_to_y(Ubuff2Ze,Ubuff2Ye,gpLESe)
          call transpose_y_to_x(Ubuff2Ye,W,gpLESe)
        case ('y')
          call transpose_y_to_z(Utmp,Ubuff2Zc,gpLESc)
          call interpLES%interpz_C2E(Ubuff2Zc,Ubuff2Ze,-1,1)
          call transpose_z_to_y(Ubuff2Ze,U,gpLESe)

          call transpose_y_to_z(Vtmp,Ubuff2Zc,gpLESc)
          call interpLES%interpz_C2E(Ubuff2Zc,Ubuff2Ze,-1,1)
          call transpose_z_to_y(Ubuff2Ze,V,gpLESe)

          call transpose_y_to_z(Wtmp,Ubuff2Zc,gpLESc)
          call interpLES%interpz_C2E(Ubuff2Zc,Ubuff2Ze,-1,1)
          call transpose_z_to_y(Ubuff2Ze,W,gpLESe)
        case ('z')
          call interpLES%interpz_C2E(Utmp,U,-1,1)
          call interpLES%interpz_C2E(Vtmp,V,-1,1)
          call interpLES%interpz_C2E(Wtmp,W,-1,-1)
        end select
        call message("Finished interpolating LES velocity "//&
          "to cell edges") 
      else
        call readFields(fname,U,V,W,'/u','/v','/w',gpLESe)
        readUVW = .true.
        call computeLargeScaleGradient(fname)
      end if

      call message('Warning: Assuming x and y are homogeneous directions')
      call getStartAndEndIndices(gpLESe,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
      allocate(Ubar(ksz),Vbar(ksz),Wbar(ksz))
      call meanXY(U,Ubar)
      call meanXY(V,Vbar)
      call meanXY(W,Wbar)
      do k = 1,ksz
        Up(:,:,k) = U(:,:,k) - Ubar(k)
        Vp(:,:,k) = V(:,:,k) - Vbar(k)
        Wp(:,:,k) = W(:,:,k) - Wbar(k)
      end do
      deallocate(Ubar,Vbar,Wbar)

      call update_halo(U,Uh,1)
!print*, minval(abs(Uh(:,jsz+1,:))) !DEBUG
!call MPI_Barrier(MPI_COMM_WORLD,ierr) !DEBUG
!call assert(.false.) !DEBUG
      call update_halo(V,Vh,1)
      call update_halo(W,Wh,1)

      call update_halo(Up,Uph,1)
      call update_halo(Vp,Vph,1)
      call update_halo(Wp,Wph,1)

    end subroutine

    subroutine meanXY(f,fbar)
      real(rkind), dimension(:,:,:), intent(in) :: f
      real(rkind), dimension(:), intent(out) :: fbar

      fbar = 0.d0
      fbar = sum(sum(f,1),1)/real(size(f,1),rkind)/real(size(f,2),rkind)
    end subroutine

    subroutine computeLargeScaleGradient(fname)
      character(len=*), intent(in), optional :: fname
      integer :: ierr, i, j
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
            call transpose_x_to_y(U,Ubuff2Ye,gpLESe)
            call grad%ddy(Ubuff2Ye,dUbuff1)
            call transpose_y_to_x(dUbuff1,gradU(1,2,:,:,:),gpLESe)
            
            call transpose_x_to_y(V,Ubuff2Ye,gpLESe)
            call grad%ddy(Ubuff2Ye,dUbuff1)
            call transpose_y_to_x(dUbuff1,gradU(2,2,:,:,:),gpLESe)
            
            call transpose_x_to_y(W,Ubuff2Ye,gpLESe)
            call grad%ddy(Ubuff2Ye,dUbuff1)
            call transpose_y_to_x(dUbuff1,gradU(3,2,:,:,:),gpLESe)
  
            ! z derivatives
            call transpose_x_to_y(U,Ubuff2Ye,gpLESe)
            call transpose_y_to_z(Ubuff2Ye,Ubuff2Ze,gpLESe)
            call grad%ddz(Ubuff2Ze,dUbuff2)
            call transpose_z_to_y(dUbuff2,dUbuff1,gpLESe)
            call transpose_y_to_x(dUbuff1,gradU(1,3,:,:,:),gpLESe)

            call transpose_x_to_y(V,Ubuff2Ye,gpLESe)
            call transpose_y_to_z(Ubuff2Ye,Ubuff2Ze,gpLESe)
            call grad%ddz(Ubuff2Ze,dUbuff2)
            call transpose_z_to_y(dUbuff2,dUbuff1,gpLESe)
            call transpose_y_to_x(dUbuff1,gradU(2,3,:,:,:),gpLESe)

            call transpose_x_to_y(W,Ubuff2Ye,gpLESe)
            call transpose_y_to_z(Ubuff2Ye,Ubuff2Ze,gpLESe)
            call grad%ddz(Ubuff2Ze,dUbuff2)
            call transpose_z_to_y(dUbuff2,dUbuff1,gpLESe)
            call transpose_y_to_x(dUbuff1,gradU(3,3,:,:,:),gpLESe)
          case('y')
            ! x derivatives
            call transpose_y_to_x(U,Ubuff2Xe,gpLESe)
            call grad%ddx(Ubuff2Xe,dUbuff1)
            call transpose_x_to_y(dUbuff1,gradU(1,1,:,:,:),gpLESe)

            call transpose_y_to_x(V,Ubuff2Xe,gpLESe)
            call grad%ddx(Ubuff2Xe,dUbuff1)
            call transpose_x_to_y(dUbuff1,gradU(2,1,:,:,:),gpLESe)

            call transpose_y_to_x(W,Ubuff2Xe,gpLESe)
            call grad%ddx(Ubuff2Xe,dUbuff1)
            call transpose_x_to_y(dUbuff1,gradU(3,1,:,:,:),gpLESe)

            ! y derivatives
            call grad%ddy(U,gradU(1,2,:,:,:))
            call grad%ddy(V,gradU(2,2,:,:,:))
            call grad%ddy(W,gradU(3,2,:,:,:))
           
            ! z derivatives
            call transpose_y_to_z(U,Ubuff2Ze,gpLESe)
            call grad%ddz(Ubuff2Ze,dUbuff2)
            call transpose_z_to_y(dUbuff2,gradU(1,3,:,:,:),gpLESe)

            call transpose_y_to_z(V,Ubuff2Ze,gpLESe)
            call grad%ddz(Ubuff2Ze,dUbuff2)
            call transpose_z_to_y(dUbuff2,gradU(2,3,:,:,:),gpLESe)

            call transpose_y_to_z(W,Ubuff2Ze,gpLESe)
            call grad%ddz(Ubuff2Ze,dUbuff2)
            call transpose_z_to_y(dUbuff2,gradU(3,3,:,:,:),gpLESe)
          case('z')
            ! x derivatives
            call transpose_z_to_y(U,Ubuff2Ye,gpLESe)
            call transpose_y_to_x(Ubuff2Ye,Ubuff2Xe,gpLESe)
            call grad%ddx(Ubuff2Xe,dUbuff2)
            call transpose_x_to_y(dUbuff2,dUbuff1,gpLESe)
            call transpose_y_to_z(dUbuff1,gradU(1,1,:,:,:),gpLESe)

            call transpose_z_to_y(V,Ubuff2Ze,gpLESe)
            call transpose_y_to_x(Ubuff2Ze,Ubuff2Xe,gpLESe)
            call grad%ddx(Ubuff2Xe,dUbuff2)
            call transpose_x_to_y(dUbuff2,dUbuff1,gpLESe)
            call transpose_y_to_z(dUbuff1,gradU(2,1,:,:,:),gpLESe)

            call transpose_z_to_y(W,Ubuff2Ze,gpLESe)
            call transpose_y_to_x(Ubuff2Ze,Ubuff2Xe,gpLESe)
            call grad%ddx(Ubuff2Xe,dUbuff2)
            call transpose_x_to_y(dUbuff2,dUbuff1,gpLESe)
            call transpose_y_to_z(dUbuff1,gradU(3,1,:,:,:),gpLESe)
  
            ! y derivatives
            call transpose_z_to_y(U,Ubuff2Ye,gpLESe)
            call grad%ddy(Ubuff2Ye,dUbuff1)
            call transpose_y_to_z(dUbuff1,gradU(1,2,:,:,:),gpLESe)
            
            call transpose_z_to_y(V,Ubuff2Ye,gpLESe)
            call grad%ddy(Ubuff2Ye,dUbuff1)
            call transpose_y_to_z(dUbuff1,gradU(2,2,:,:,:),gpLESe)
            
            call transpose_z_to_y(W,Ubuff2Ye,gpLESe)
            call grad%ddy(Ubuff2Ye,dUbuff1)
            call transpose_y_to_z(dUbuff1,gradU(3,2,:,:,:),gpLESe)
  
            ! z derivatives
            call grad%ddz(U,gradU(1,3,:,:,:))
            call grad%ddz(V,gradU(2,3,:,:,:))
            call grad%ddz(W,gradU(3,3,:,:,:))
          case default
            call gracefulExit("ERROR: must set decomp2Dpencil to 'x', 'y', or 'z'", ierr)
        end select
      else
        call readFields(fname,gradU,'/gradU',gpLESe)
      end if
      S = sqrt(sum(sum(gradU*gradU,1),1))

      ! Get halo-padded array      
      do j = 1,3
        do i = 1,3
          call update_halo(gradU(i,j,:,:,:),gradUbuff,1)
          if (i == 1 .and. j == 1) allocate(gradUh(3,3,size(gradUbuff,1),&
            size(gradUbuff,2),size(gradUbuff,3)))
          gradUh(i,j,:,:,:) = gradUbuff
        end do
      end do
      deallocate(gradUbuff)
      call update_halo(S,Sh,1)
      
      gradUcomputed = .true.
    end subroutine
end module
