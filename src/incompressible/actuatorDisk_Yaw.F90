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
        real(rkind) :: yaw, tilt, ut, power
        real(rkind) :: xLoc, yLoc, zLoc, dx, dy, dz
        real(rkind) :: diam, cT, pfactor, normfactor, OneBydelSq, Cp
        real(rkind) :: uface = 0.d0, vface = 0.d0, wface = 0.d0
        integer :: totPointsOnFace
        real(rkind), dimension(:,:), allocatable :: xp, yp, zp
        real(rkind), dimension(:), allocatable :: xs, ys, zs
        integer, dimension(:,:), allocatable :: startEnds

        ! Grid Info
        integer :: nxLoc, nyLoc, nzLoc 
        real(rkind) :: delta ! Smearing size
        real(rkind) :: alpha_tau = 1.d0! Smoothing parameter (set to 1 for initialization) 
        real(rkind), dimension(:), allocatable :: xline, yline, zline
        real(rkind), dimension(:,:,:), pointer :: xG, yG, zG
        
        ! Pointers to memory buffers 
        logical :: memory_buffers_linked = .false.
        real(rkind), dimension(:,:,:), pointer :: rbuff, blanks, speed, scalarSource

        ! MPI communicator stuff
        logical :: Am_I_Active, Am_I_Split
        integer :: color, myComm, myComm_nproc, myComm_nrank

    contains
        procedure :: init
        procedure :: destroy
        procedure :: get_RHS
        procedure :: link_memory_buffers
        procedure, private :: AD_force
        procedure, private :: AD_force_point
        !procedure, private :: Sfunc 
        procedure :: get_power
        procedure :: dumpPower
        procedure :: dumpPowerUpdate
    end type


contains

subroutine init(this, inputDir, ActuatorDisk_T2ID, xG, yG, zG)
    class(actuatorDisk_yaw), intent(inout) :: this
    real(rkind), intent(in), dimension(:,:,:), target :: xG, yG, zG
    integer, intent(in) :: ActuatorDisk_T2ID
    character(len=*), intent(in) :: inputDir
    character(len=clen) :: tempname, fname
    integer :: ioUnit
    real(rkind) :: xLoc=1.d0, yLoc=1.d0, zLoc=0.1d0
    real(rkind) :: diam=0.08d0, cT=0.65d0, yaw=0.d0, tilt=0.d0, Cp = 0.3
 
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
    this%cT = cT; this%diam = diam; this%yaw = yaw; this%ut = 1.d0; this%Cp = Cp;
    this%nxLoc = size(xG,1); this%nyLoc = size(xG,2); this%nzLoc = size(xG,3)

    !Allocate stuff
    allocate(this%xLine(size(xG,1)))
    allocate(this%yLine(size(xG,2)))
    allocate(this%zLine(size(xG,3)))
    allocate(this%rbuff(size(xG,1), size(xG,2), size(xG,3)))


    this%xG => xG; this%yG => yG; this%zG => zG
    this%memory_buffers_linked = .false.
   
 
end subroutine 

! Need to create pointers instead of allocating fresh memory for scaling (in
! termsof numturbines) and performance (allocating/deallocating is expensive)
subroutine link_memory_buffers(this, rbuff, blanks, speed, scalarSource)
    class(actuatordisk_yaw), intent(inout) :: this
    real(rkind), dimension(:,:,:), intent(in), target :: rbuff, blanks, speed, scalarSource
    this%rbuff => rbuff
    this%speed => speed
    this%scalarSource => scalarSource
    this%blanks => blanks
    this%memory_buffers_linked = .true. 
end subroutine 

subroutine destroy(this)
    class(actuatordisk_yaw), intent(inout) :: this

    this%memory_buffers_linked = .false.
    
    nullify(this%rbuff, this%blanks, this%speed, this%scalarSource)
    nullify(this%xG, this%yG, this%zG)
end subroutine 

subroutine get_RHS(this, u, v, w, rhsxvals, rhsyvals, rhszvals, gamma_negative, theta)
    class(actuatordisk_yaw), intent(inout) :: this
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsxvals, rhsyvals, rhszvals
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(in)    :: u, v, w
    real(rkind), intent(in) :: gamma_negative, theta
    real(rkind) :: usp_sq, force, gamma
    real(rkind), dimension(3,3) :: R, T
    real(rkind), dimension(3,1) :: xn, Ft, n
    real(rkind) ::  numPoints, x, y, z, scalarSource, sumVal
    real(rkind) :: xnew, ynew, znew, cgamma, sgamma, ctheta, stheta, x2, y2, z2
    integer :: i,j,k
     
    if (this%memory_buffers_linked) then
        ! Normal vector
        ! Adjust the yaw misalignment angle as per the convention of
        ! Howland et al. 2019 and Shapiro et al. 2019
        gamma = -gamma_negative
        xn = reshape([1.0d0, 0.d0, 0.d0], shape(xn))
        R = reshape([cos(gamma), sin(gamma), 0.d0, -sin(gamma), cos(gamma), 0.d0, 0.d0, 0.d0, 1.d0], shape(R))
        T = reshape([1.d0, 0.d0, 0.d0, 0.d0, cos(theta), sin(theta), 0.d0, -sin(theta), cos(theta)], shape(T)) 
        n = matmul(matmul(transpose(R), transpose(T)), xn)
        
        ! Above this is fast
        ! Translate and rotate domain
        ctheta = cos(theta); stheta = sin(theta)
        cgamma = cos(gamma); sgamma = sin(gamma)
        this%scalarSource = 0.d0

        do k = 1,size(this%xG,3)
            z = this%zG(1,1,k) - this%zLoc
            do j = 1,size(this%xG,2)
                y = this%yG(1,j,1) - this%yLoc
                !$omp simd
                do i = 1,size(this%xG,1)  
                  x = this%xG(i,1,1) - this%xLoc
                  xnew = x * cgamma - y * sgamma
                  ynew = x * sgamma + y * cgamma
                  znew = z
                  x2 = xnew
                  y2 = ynew*ctheta - znew*stheta
                  z2 = ynew*stheta + znew*ctheta
                  call this%AD_force_point(x2, y2, z2, scalarSource)
                  this%scalarSource(i,j,k) = scalarSource
                end do
            end do
        end do
        ! this part needs to be done at the end of the loops
        sumVal = p_sum(this%scalarSource) * this%dx*this%dy*this%dz
        this%scalarSource = this%scalarSource / sumVal
        ! Get the mean velocities at the turbine face
        this%speed = u*n(1,1) + v*n(2,1) + w*n(3,1)
        this%blanks = 1.d0
        do k = 1,size(this%xG,3)
            do j = 1, size(this%xG,2)
                do i = 1,size(this%xG,1)
                    if (abs(this%scalarSource(i,j,k))<1D-10) then
                        this%blanks(i,j,k) = 0.d0
                    end if
                end do
            end do
        end do
        numPoints = p_sum(this%blanks)
        this%rbuff = this%blanks*this%speed
        this%ut = p_sum(this%rbuff)/numPoints    
        ! Mean speed at the turbine
        usp_sq = (this%ut)**2
        force = -0.5d0*this%cT*(pi*(this%diam**2)/4.d0)*usp_sq
        Ft = reshape([force, 0.d0, 0.d0], shape(Ft))
        Ft = matmul(matmul(transpose(R), transpose(T)), Ft)
        ! Can we avoid the above matmul? Just write it out
        rhsxvals = rhsxvals + Ft(1,1) * this%scalarSource 
        rhsyvals = rhsyvals + Ft(2,1) * this%scalarSource
        rhszvals = rhszvals + Ft(3,1) * this%scalarSource 

    end if 

end subroutine

subroutine get_power(this)
    class(actuatordisk_yaw), intent(inout) :: this

    this%power = 0.5d0*this%Cp*(pi*(this%diam**2)/4.d0)*(this%ut)**3

end subroutine

subroutine AD_force_point(this, X, Y, Z, scalarSource)
    class(actuatordisk_yaw), intent(inout) :: this
    real(rkind), intent(out) :: scalarSource
    real(rkind), intent(in) :: X,Y,Z
    real(rkind) :: delta_r = 0.8d0, smear_x = 1.5d0, delta, R, sumVal
    real(rkind) :: tmp, tmp2, diamFactor = 1.8d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Do everything in the i,j,k loop in the rhs call

    ! X,Y,Z are shifted to xc, yc, zx zero center as per AD location
    R = this%diam / 2.d0 * diamFactor ! Added to try and smear the forcing out more
    tmp = sqrt((Y/R)**2 + (Z/R)**2)
    tmp = (tmp-1.d0)/delta_r + 1.d0
    call Sfunc_point(tmp, tmp2)
    tmp = 1.d0 - tmp2
    delta = (this%dx)*smear_x
    scalarSource = tmp * (1.d0/(delta*sqrt(2.d0*pi))) * exp(-0.5d0*(X**2)/(delta**2))

end subroutine


subroutine AD_force(this, X, Y, Z, scalarSource)
    class(actuatordisk_yaw), intent(inout) :: this
    real(rkind), intent(inout), dimension(:,:,:) :: scalarSource
    real(rkind), intent(in), dimension(:,:,:) :: X, Y, Z
    real(rkind) :: delta_r = 0.22d0, smear_x = 1.5d0, delta, R, sumVal

    ! X,Y,Z are shifted to xc, yc, zx zero center as per AD location
    R = this%diam / 2.d0
    this%rbuff = sqrt((Y/R)**2 + (Z/R)**2)
    this%rbuff = (this%rbuff-1.d0)/delta_r + 1.d0
    !call Sfunc(this%rbuff, this%rbuff) << Can't operate in-place here because Sfunc uses in/out declarations.
    call Sfunc(this%rbuff, scalarSource)
    this%rbuff = 1.d0 - scalarSource
    delta = (this%dx)*smear_x
    scalarSource = this%rbuff * (1.d0/(delta*sqrt(2.d0*pi))) * exp(-0.5d0*(X**2)/(delta**2))
    sumVal = p_sum(scalarSource) * this%dx*this%dy*this%dz
    scalarSource = scalarSource / sumVal

end subroutine

pure subroutine Sfunc_point(x, val)
    !class(actuatordisk_yaw), intent(inout) :: this
    real(rkind), intent(in)  :: x
    real(rkind), intent(out) :: val
    
    val = 1.d0 / (1.d0 + exp(min(1.d0/(x-1.d0+1D-18) + 1.d0/(x + 1D-18), 50.d0)))
    if (x>1.d0) then
        val = 1.d0
    end if
    if (x<0.d0) then
        val = 0.d0
    end if

end subroutine

pure subroutine Sfunc(x, val)
    !class(actuatordisk_yaw), intent(inout) :: this
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

subroutine dumpPower(this, outputfile, tempname)
    class(actuatordisk_yaw), intent(inout) :: this
    character(len=*),    intent(in)            :: outputfile, tempname
    integer :: fid = 1234
    character(len=clen) :: fname

    ! Get power
    call this%get_power()
    ! Write power
    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname)
    !open(fid,file=trim(fname), form='unformatted',action='write',position='append')
    open(fid,file=trim(fname), form='formatted', action='write',position='append')
    write(fid, *) this%power
    close(fid)

end subroutine    

subroutine dumpPowerUpdate(this, outputfile, tempname, & 
                           powerUpdate, Phat, yaw, yawOld, meanP, kw, sigma, i)
    class(actuatordisk_yaw), intent(inout) :: this
    character(len=*),    intent(in)            :: outputfile, tempname
    integer :: fid = 1234
    integer, intent(in) :: i
    character(len=clen) :: fname, tempname2
    real(rkind), dimension(:), intent(in) :: powerUpdate, Phat, yaw, yawOld, meanP
    real(rkind), dimension(:), intent(in) :: kw, sigma

    ! Write power
    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname)
    !open(fid,file=trim(fname), form='unformatted',action='write',position='append')
    open(fid,file=trim(fname), form='formatted')
    write(fid, *) powerUpdate
    close(fid)
    ! Write Phat
    write(tempname2,"(A5,I3.3,A4)") "Phat_",i,".txt"
    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
    open(fid,file=trim(fname), form='formatted')
    write(fid, *) Phat
    close(fid)
    ! Write yaw opti
    write(tempname2,"(A8,I3.3,A4)") "YawOpti_",i,".txt"
    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
    open(fid,file=trim(fname), form='formatted')
    write(fid, *) yaw
    close(fid)
    ! Write yaw original in this time interval
    write(tempname2,"(A12,I3.3,A4)") "YawInterval_",i,".txt"
    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
    open(fid,file=trim(fname), form='formatted')
    write(fid, *) yawOld
    close(fid)
    ! Write mean power in this time interval
    write(tempname2,"(A6,I3.3,A4)") "MeanP_",i,".txt"
    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
    open(fid,file=trim(fname), form='formatted')
    write(fid, *) meanP
    close(fid)
    ! Write kw in this time interval
    write(tempname2,"(A3,I3.3,A4)") "kw_",i,".txt"
    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
    open(fid,file=trim(fname), form='formatted')
    write(fid, *) kw
    close(fid)
    ! Write sigma in this time interval
    write(tempname2,"(A6,I3.3,A4)") "sigma_",i,".txt"
    fname = outputfile(:len_trim(outputfile))//"/"//trim(tempname2)
    open(fid,file=trim(fname), form='formatted')
    write(fid, *) sigma
    close(fid)

end subroutine    

end module 
