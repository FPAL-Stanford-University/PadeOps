module actuatorDisk_YawMod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa
    use decomp_2d
    use exits, only: GracefulExit, message
    use mpi 
    use reductions, only: p_maxval, p_sum
    use timer, only: tic, toc
    use Gridtools, only: linspace

    implicit none

    private
    public :: actuatorDisk_yaw
    
    real(rkind), parameter :: alpha_Smooth = 1.d0 ! 0.9d0 ! Exonential smoothing constant
    integer, parameter :: xReg = 8, yReg = 8, zReg = 8

    type :: actuatorDisk_yaw
        ! Actuator Disk_T2 Info
        integer :: xLoc_idx, ActutorDisk_T2ID
        integer, dimension(:,:), allocatable :: tag_face 
        real(rkind) :: yaw, tilt, ut
        real(rkind) :: xLoc, yLoc, zLoc, dx, dy, dz
        real(rkind) :: diam, cT, pfactor, normfactor, OneBydelSq
        real(rkind) :: uface = 0.d0, vface = 0.d0, wface = 0.d0
        integer :: totPointsOnFace
        real(rkind), dimension(:,:,:), allocatable :: eta_delta, dsq
        real(rkind), dimension(:,:), allocatable :: xp, yp, zp
        real(rkind), dimension(:), allocatable :: xs, ys, zs
        integer, dimension(:,:), allocatable :: startEnds

        ! Grid Info
        integer :: nxLoc, nyLoc, nzLoc 
        real(rkind) :: delta ! Smearing size
        real(rkind) :: alpha_tau = 1.d0! Smoothing parameter (set to 1 for initialization) 
        real(rkind), dimension(:,:,:), allocatable :: rbuff
        real(rkind), dimension(:), allocatable :: dline, xline, yline, zline
        real(rkind), dimension(:,:,:), pointer :: xG, yG, zG
        real(rkind), dimension(:,:,:), allocatable :: smearing_base

        ! MPI communicator stuff
        logical :: Am_I_Active, Am_I_Split
        integer :: color, myComm, myComm_nproc, myComm_nrank

    contains
        procedure :: init
        procedure :: destroy
        procedure :: get_RHS
        procedure, private :: AD_force
        procedure, private :: Sfunc 
        procedure :: get_power
    end type


contains

subroutine init(this, inputDir, ActuatorDisk_T2ID, xG, yG, zG)
    class(actuatorDisk_yaw), intent(inout) :: this
    real(rkind), intent(in), dimension(:,:,:), target :: xG, yG, zG
    integer, intent(in) :: ActuatorDisk_T2ID
    character(len=*), intent(in) :: inputDir
    character(len=clen) :: tempname, fname
    integer :: ioUnit, tmpSum, totSum
    real(rkind) :: xLoc=1.d0, yLoc=1.d0, zLoc=0.1d0, diam=0.08d0, cT=0.65d0, yaw=0.d0, tilt=0.d0
      
    ! Read input file for this turbine    
    namelist /ACTUATOR_DISK/ xLoc, yLoc, zLoc, diam, cT, yaw, tilt
    write(tempname,"(A13,I4.4,A10)") "ActuatorDisk_", ActuatorDisk_T2ID, "_input.inp"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)

    ioUnit = 55
    open(unit=ioUnit, file=trim(fname), form='FORMATTED', action="read")
    read(unit=ioUnit, NML=ACTUATOR_DISK)
    close(ioUnit)
   
    this%dx=xG(2,1,1)-xG(1,1,1); this%dy=yG(1,2,1)-yG(1,1,1); this%dz=zG(1,1,2)-zG(1,1,1)
    this%xLoc = xLoc; this%yLoc = yLoc; this%zLoc = zLoc
    this%cT = cT; this%diam = diam; this%yaw = this%yaw
    this%nxLoc = size(xG,1); this%nyLoc = size(xG,2); this%nzLoc = size(xG,3)

    !Allocate stuff
    allocate(this%tag_face(size(xG,2),size(xG,3)))
    allocate(this%xLine(size(xG,1)))
    allocate(this%yLine(size(xG,2)))
    allocate(this%zLine(size(xG,3)))
    allocate(this%dsq(2*xReg+1,2*yReg+1,2*zReg+1))
    allocate(this%dline(2*xReg+1))
    allocate(this%rbuff(size(xG,1), size(xG,2), size(xG,3)))


    ! Declarations and targets
    this%dsq = 0.d0
    this%xG => xG; this%yG => yG; this%zG => zG
   
 
end subroutine 

subroutine destroy(this)
    class(actuatordisk_yaw), intent(inout) :: this

    if (Allocated(this%rbuff))  deallocate(this%rbuff)
    if (Allocated(this%tag_face))  deallocate(this%tag_face)
    
    nullify(this%xG, this%yG, this%zG)
end subroutine 

subroutine get_RHS(this, u, v, w, rhsxvals, rhsyvals, rhszvals, gamma_negative, theta)
    class(actuatordisk_yaw), intent(inout) :: this
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsxvals, rhsyvals, rhszvals
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(in)    :: u, v, w
    real(rkind) :: usp_sq, force, gamma
    real(rkind), dimension(3,3) :: R, T
    real(rkind), dimension(3,1) :: xn, Ft, n
    real(rkind), intent(in) :: gamma_negative, theta
    integer :: tmpSum, totSum
    real(rkind) :: xLoc=1.d0, yLoc=1.d0, zLoc=0.1d0, diam=0.08d0, cT=0.65d0
    real(rkind) :: tilt=0.d0, epsFact = 1.5d0
    real(rkind), dimension(:,:,:), allocatable :: tmp
    integer, dimension(:,:,:), allocatable :: tmp_tag
    integer :: j, locator(1)
    integer :: xLc(1), yLc(1), zLc(1), xst, xen, yst, yen, zst, zen, ierr, xlen
    integer  :: ntry = 100
    real(rkind) :: time2initialize = 0, correction_factor = 1.0d0, normfact_p, numPoints
    real(rkind), dimension(:,:,:), allocatable :: blanks, speed, X, Y, Z, scalarSource
    real(rkind), dimension(:,:,:), allocatable :: Xnew, Ynew, Znew
     
    ! allocate stuff
    allocate(blanks(size(this%xG,1), size(this%xG,2), size(this%xG,3)))
    allocate(speed(size(this%xG,1), size(this%xG,2), size(this%xG,3)))
    allocate(X(size(this%xG,1), size(this%xG,2), size(this%xG,3)))
    allocate(Y(size(this%xG,1), size(this%xG,2), size(this%xG,3)))
    allocate(Z(size(this%xG,1), size(this%xG,2), size(this%xG,3)))
    allocate(Xnew(size(this%xG,1), size(this%xG,2), size(this%xG,3)))
    allocate(Ynew(size(this%xG,1), size(this%xG,2), size(this%xG,3)))
    allocate(Znew(size(this%xG,1), size(this%xG,2), size(this%xG,3)))
    allocate(scalarSource(size(this%xG,1), size(this%xG,2), size(this%xG,3)))
    allocate(tmp(size(this%xG,1), size(this%xG,2), size(this%xG,3)))

    ! Normal vector
    ! Adjust the yaw misalignment angle as per the convention of
    ! Howland et al. 2019 and Shapiro et al. 2019
    gamma = -gamma_negative
    xn = reshape((/1.0d0, 0.d0, 0.d0 /), shape(xn))
    R = reshape((/cos(gamma), sin(gamma), 0.d0, -sin(gamma), cos(gamma), 0.d0, 0.d0, 0.d0, 1.d0/), shape(R))
    T = reshape((/1.d0, 0.d0, 0.d0, 0.d0, cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta)/), shape(T)) 
    n = matmul(matmul(transpose(R), transpose(T)), xn)
    ! Translate and rotate domain
    X = this%xG-this%xLoc 
    Y = this%yG-this%yLoc
    Z = this%zG-this%zLoc
    Xnew = X*cos(gamma) - Y*sin(gamma)
    Ynew = X*sin(gamma) + Y*cos(gamma)
    Znew = Z
    X = Xnew
    Y = Ynew*cos(theta) - Znew*sin(theta)
    Z = Ynew*sin(theta) + Znew*cos(theta)
    call this%AD_force(X, Y, Z, scalarSource)
    ! Get the mean velocities at the turbine face
    speed = u*n(1,1) + v*n(2,1) + w*n(3,1)
    blanks = 1.
    where(abs(scalarSource)<1D-10)
        blanks = 0.
    end where
    numPoints = p_sum(blanks)
    this%ut = p_sum(blanks*speed)/numPoints    
    ! Mean speed at the turbine
    usp_sq = (this%ut)**2
    force = -this%pfactor*this%normfactor*0.5d0*this%cT*(pi*(this%diam**2)/4.d0)*usp_sq
    Ft = reshape((/ force, 0.d0, 0.d0 /), shape(Ft))
    Ft = matmul(matmul(transpose(R), transpose(T)), Ft)
    rhsxvals = rhsxvals + Ft(1,1) * scalarSource 
    rhsyvals = rhsyvals + Ft(2,1) * scalarSource
    rhszvals = rhszvals + Ft(3,1) * scalarSource

    ! Deallocate
    deallocate(blanks)
    deallocate(speed)
    deallocate(X)
    deallocate(Y)
    deallocate(Z)
    deallocate(Xnew)
    deallocate(Ynew)
    deallocate(Znew)
    deallocate(scalarSource)
    deallocate(tmp)

end subroutine

subroutine get_power(this)
    class(actuatordisk_yaw), intent(inout) :: this

    this%power = -this%pfactor*this%normfactor*0.5d0*this%cT*(pi*(this%diam**2)/4.d0)*(this%ut)**3

end subroutine

subroutine AD_force(this, X, Y, Z, scalarSource)
    class(actuatordisk_yaw), intent(inout) :: this
    real(rkind), intent(inout), dimension(:,:,:) :: scalarSource
    real(rkind), intent(in), dimension(:,:,:) :: X, Y, Z
    real(rkind) :: delta_r = 0.22d0, smear_x = 1.5d0, delta, R, sumVal

    ! X,Y,Z are shifted to xc, yc, zx zero center as per AD location
    R = this%diam / 2.
    this%rbuff = sqrt((Y/R)**2 + (Z/R)**2)
    this%rbuff = (this%rbuff-1.)/delta_r + 1.
    call this%Sfunc(this%rbuff, this%rbuff)
    this%rbuff = 1.-this%rbuff
    delta = (this%dx)*smear_x
    scalarSource = this%rbuff * (1./(delta*sqrt(2.*pi))) * exp(-0.5*(X**2)/(delta**2))
    sumVal = p_sum(scalarSource) * this%dx**3
    scalarSource = scalarSource / sumVal

end subroutine

subroutine Sfunc(this, x, val)
    class(actuatordisk_yaw), intent(inout) :: this
    real(rkind), dimension(:,:,:), intent(in)  :: x
    real(rkind), dimension(:,:,:), intent(out) :: val
    
    val = 1. / (1. + exp(min(1./(x-1.+1D-18) + 1./(x + 1D-18), 50.d0)))
    where(x>1.)
        val = 1.
    end where
    where(x<0.)
        val = 0.
    end where

end subroutine

end module 
