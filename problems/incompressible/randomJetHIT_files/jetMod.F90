module jetMod
  use kind_parameters, only: rkind, clen
  use random,          only: gaussian_random, uniform_random
  use exits,           only: message
  implicit none
  private

  public :: jet, sourceFraction, muOn, sigmaFact, outputdir, wJet

  real(rkind) :: sourceFraction = 0.125d0
  real(rkind) :: muOn = 3.d0
  real(rkind) :: sigmaFact = 0.3333333333333333d0
  real(rkind) :: muOff, sigmaOff, sigmaOn
  logical :: moduleParametersAreSet = .false.
  logical :: dumpJetInfo = .false.
  real(rkind) :: wJet = 1.d0, uHjet = 0.d0

  ! Scaling fact for tanh(t) profile for smooth on/off of the jet
  real(rkind) :: rampOnOffFact = 1.e9
  
  ! 99.9% thickness of tanh profile due to rampOnOffFact
  real(rkind) :: tanhDelta

  character(len=clen) :: outputdir

  type jet
    private
    real(rkind), public :: xsz, ysz, zsz
    real(rkind), public :: xloc, yloc ! Center of the jet
    real(rkind), public :: zst ! The start of the jet
    real(rkind) :: tstart, duration, lastWeight
    real(rkind), public :: weight
    real(rkind), public :: u, v, w
    logical :: IamOn
    logical :: initialized = .false.
    integer :: gID

    ! Random number generation
    integer :: seed
    real(rkind) :: logMap

    contains
      procedure          :: init
      procedure          :: destroy
      procedure          :: isOn
      procedure, private :: getOnOffDuration
      procedure, private :: setJetVelocity
      procedure, private :: updateSeed
      procedure          :: getStatus
      procedure          :: getTstart
      procedure          :: getDuration
      procedure          :: getID
      procedure          :: getLocation
      procedure          :: print
      procedure          :: testRandomSampling
      procedure          :: testUpdateWeight
      procedure, private :: updateWeight
  end type

  contains

    subroutine setModuleParameters(inputfile)
      character(len=*), intent(in) :: inputfile
      integer :: ioUnit
    
      namelist /JET_INPUT/ muOn, sourceFraction, sigmaFact, dumpJetInfo, &
        rampOnOffFact, uHjet, wJet, outputdir

      ioUnit = 11
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
      read(unit=ioUnit, NML=JET_INPUT)
      close(ioUnit)    
   
      muOff = (1.d0/sourceFraction - 1.d0)*muOn
      sigmaOff = sigmaFact*muOff
      sigmaOn = sigmaFact*muOn
      moduleParametersAreSet = .true.
      tanhDelta = max(getDelta(rampOnOffFact), 0.d0)
    end subroutine

    subroutine init(this,xloc,yloc,zst,xsz,ysz,zsz,gID,inputfile,NjetsTotal,tInit)
      use seedGen, only: initializeLogMap
      class(jet), intent(inout) :: this
      real(rkind), intent(in) :: xloc, yloc, zst, xsz, ysz, zsz, tInit
      integer, intent(in) :: gID ! unique identifier for the given jet
      character(len=*), intent(in) :: inputfile
      integer, intent(in) :: NjetsTotal
      real(rkind) :: logMapInit

      if (.not. moduleParametersAreSet) call setModuleParameters(inputfile)

      this%xsz = xsz
      this%ysz = ysz
      this%zsz = zsz

      this%xloc = xloc
      this%yloc = yloc
      this%zst  = zst
      
      this%gID = gID

      ! Initialize the seed machinery (i.e. the logistic map sequence)
      logMapInit = real(gID,rkind)/(real(NjetsTotal,rkind) + 0.5d0)
      this%logMap = initializeLogMap(gID,x0=logMapInit)

      ! Initialize the jet state
      this%tstart = 0.d0
      this%duration = 0.d0
      this%w = wJet
      call this%getOnOffDuration(tInit)
      call this%setJetVelocity()
      call this%print()

      this%initialized = .true.
    end subroutine

    subroutine destroy(this)
      class(jet), intent(inout) :: this

      ! Don't really do anything
      this%xloc = 0.d0
    end subroutine

    subroutine getOnOffDuration(this,tnow)
      class(jet), intent(inout) :: this
      real(rkind), intent(in) :: tnow
      real(rkind), dimension(2) :: rand

      if (.not. this%initialized) then ! Choose on or off
        call this%updateSeed()
        call uniform_random(rand(1),0.d0,1.d0,this%seed)

        if (rand(1) > sourceFraction) then
          this%IamOn = .true.  ! This gets turned to false below
        else
          this%IamOn = .false. ! This gets turned to true below
        end if
        this%lastWeight = 0.d0
      end if
      
      if (this%tstart + this%duration > tnow) then
        continue
      else
        call this%updateSeed()
        if (this%IamOn) then
          this%IamOn = .false.
          call gaussian_random(rand,muOff,sigmaOff,this%seed)
        else
          this%IamOn = .true.
          call gaussian_random(rand,muOn,sigmaOn,this%seed)
        end if
        ! Note, this does allow for the rare event that duration < 0.d0. All
        ! this means is that the jet will turn on/off for a single time step
        ! before being switched on and given a new duration (since tstart +
        ! duration < tnow the next time this routine is called). The same result
        ! would happen if we "clipped" the value to 0 which is what is suggested
        ! in the paper
        this%lastWeight = this%weight
        this%duration = rand(1)
        this%tstart = tnow
      end if

    end subroutine
    
    subroutine setJetVelocity(this)    
      class(jet), intent(inout) :: this
      integer :: sgn

      if (uHjet > 1.d-14 .and. this%IamOn) then
        call uniform_random(this%u,-uHjet,uHjet,this%seed)
        sgn = nint(abs(this%u)/uHjet)
        if (sgn == 0) sgn = -1
        this%v = sgn*sqrt(this%w*this%w - this%u*this%u)
      else
        this%u = 0.d0
        this%v = 0.d0
      end if
    end subroutine

    subroutine print(this)
      class(jet), intent(inout) :: this

      if (dumpJetInfo) then
        call message(0,'Attributes for jet ID ',this%gID)
        if (this%Iamon) then
          call message(1,'On duration', this%duration)
          call message(2,'uJet', this%u)
          call message(2,'vJet', this%v)
          call message(2,'wJet', this%w)
        else
          call message(1,'Off duration', this%duration)
        end if
        call message(2,'(x,y) location', [this%xloc, this%yloc])
        call message(2,'(zst,zen)', [this%zst, this%zst + this%zsz])

      end if
    end subroutine

    function isOn(this,time) result(IamOn)
      class(jet), intent(inout) :: this
      real(rkind), intent(in) :: time
      logical :: IamOn

      call this%getOnOffDuration(time)
      call this%setJetVelocity()
      call this%updateWeight(time)
      IamOn = this%IamOn
    end function

    subroutine updateWeight(this,tnow)
      class(jet), intent(inout) :: this
      real(rkind), intent(in) :: tnow

      if (this%IamOn) then
        this%weight = rampOn(rampOnOffFact,tnow,this%tstart,tanhDelta,this%lastWeight)
      else
        this%weight = rampOff(rampOnOffFact,tnow,this%tstart,tanhDelta,this%lastWeight)
      end if

    end subroutine


    subroutine testUpdateWeight(this)
      class(jet), intent(inout) :: this
      real(rkind) :: dt, time
      real(rkind), dimension(3) :: avec
      integer :: n, tid
      logical :: thisStatus

      avec = [50.d0,100.d0,200.d0]
      dt = 1.d-3
      do n = 1,size(avec)
        this%IamOn = .true.
        rampOnOffFact = avec(n)
        tanhDelta = getDelta(rampOnOffFact)
        this%tstart = 0.d0
        this%duration = 1.d0
        this%initialized = .true.
        this%lastWeight = 0.d0
        call message('a',rampOnOffFact)
        call message('------------------------------')
        do tid = 1,nint(3.d0/dt)
          time = real(tid-1,rkind)*dt
          thisStatus = this%isOn(time)
          print*, this%weight
        end do
        call message(' ')
      end do
    end subroutine

    subroutine testRandomSampling(this,Nsamples)
      use fortran_assert, only: assert
      class(jet), intent(inout) :: this
      integer, intent(in) :: Nsamples
      integer :: n

      call assert(dumpJetInfo,'Must set "dumpJetInfo" to true to run this test')
      do n = 1,Nsamples
        call this%getOnOffDuration(this%tstart + this%duration + 0.1d0)
        call this%print()
      end do
    end subroutine

    subroutine updateSeed(this)
      use constants, only: rmaxInt
      use seedGen, only: incrementLogisticMap
      class(jet), intent(inout) :: this
      real(rkind) :: old
      
      old = this%logMap
      call incrementLogisticMap(old,this%logMap)

      this%seed = nint(this%logMap*rmaxInt)
    end subroutine

    function getStatus(this) result(myStatus)
      class(jet), intent(inout) :: this
      logical :: myStatus
      myStatus = this%IamOn
    end function
    
    function getDuration(this) result(myDuration)
      class(jet), intent(inout) :: this
      real(rkind) :: myDuration
      myDuration = this%duration
    end function
    
    function getID(this) result(gID)
      class(jet), intent(inout) :: this
      integer :: gID
      gID = this%gID
    end function
    
    function getLocation(this) result(myLocation)
      class(jet), intent(inout) :: this
      real(rkind), dimension(2) :: myLocation
      myLocation(1) = this%xloc
      myLocation(2) = this%yloc
    end function
    
    function getTstart(this) result(myTstart)
      class(jet), intent(inout) :: this
      real(rkind) :: myTstart
      myTstart = this%tstart
    end function
    
    function rampOn(a,x,tOn,delta,flast) result(f)
      real(rkind), intent(in) :: a, x, tOn, delta, flast
      real(rkind) :: f, scaleFact
   
      f = 0.5d0*(1.d0 + tanh(a*(x-(tOn + 0.5d0*delta))))
      scaleFact = (1.d0 - flast)
      f = flast + scaleFact*f
    end function

    function rampOff(a,x,tOff,delta,flast) result(f)
      real(rkind), intent(in) :: a, x, tOff, delta, flast
      real(rkind) :: f

      f = 1-rampOn(a,x,tOff,delta,1-flast)
    end function

    function getDelta(a) result(delta)
      use constants, only: pi
      real(rkind), intent(in) :: a
      real(rkind) :: delta
      real(rkind), parameter :: topCutOff = 0.999d0, botCutOff = 0.001d0
      real(rkind) :: dx, x1, x2, f

      dx = 1e-3
      x1 = -pi
      f = rampOn(a,x1,0.d0,0.d0,0.d0)
      if (f > botCutOff + 1e-5) then
        do while (f > botCutOff)
          x1 = x1 - dx
          f = rampOn(a,x1,0.d0,0.d0,0.d0)
        end do
      else if (f < botCutOff - 1e-5) then
        do while (f < botCutOff)
          x1 = x1 + dx
          f = rampOn(a,x1,0.d0,0.d0,0.d0)
        end do
      end if
      
      x2 = pi
      f = rampOn(a,x2,0.d0,0.d0,0.d0)
      if (f < topCutOff - 1e-5) then
        do while (f < topCutOff)
          x2 = x2 + dx
          f = rampOn(a,x2,0.d0,0.d0,0.d0)
        end do
      else if (f > topCutOff + 1e-5) then
        do while (f > topCutOff)
          x2 = x2 - dx
          f = rampOn(a,x2,0.d0,0.d0,0.d0)
        end do
      end if
      
      delta = x2 - x1;
    end function


end module
    
