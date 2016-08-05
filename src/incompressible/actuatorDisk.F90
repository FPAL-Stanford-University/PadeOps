module actuatorDiskmod
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
    public :: actuatorDisk
    
    real(rkind), parameter :: alpha_Smooth = 0.9d0 ! Exonential smoothing constant
    integer, parameter :: xReg = 4, yReg = 7, zReg = 7
    integer, parameter :: ntry = 10

    type :: actuatorDisk
        ! Actuator Disk Info
        integer :: xLoc_idx, ActutorDiskID
        integer(rkind), dimension(:,:), allocatable :: tag_face 
        real(rkind) :: yaw, tilt
        real(rkind) :: xLoc, yLoc, zLoc
        real(rkind) :: diam, cT, pfactor, normfactor, OneBydelSq
        real(rkind) :: uface = 0.d0, vface = 0.d0, wface = 0.d0
        integer :: totPointsOnFace
        real(rkind), dimension(:,:,:), allocatable :: eta_delta, dsq, xSmall, ySmall, zSmall, source
        real(rkind), dimension(:), allocatable :: xs, ys, zs
        integer :: xst, xen, yst, yen, zst, zen, xlen, ylen, zlen
        ! Grid Info
        integer :: nxLoc, nyLoc, nzLoc 
        real(rkind) :: delta ! Smearing size
        real(rkind) :: alpha_tau = 1.d0! Smoothing parameter (set to 1 for initialization) 
        real(rkind), dimension(:,:), allocatable :: rbuff
        real(rkind), dimension(:), allocatable ::  xline, yline, zline


        ! MPI communicator stuff
        logical :: Am_I_Active, Am_I_Split
        integer :: color, myComm, myComm_nproc, myComm_nrank

    contains
        procedure :: init
        procedure :: destroy
        procedure, private :: getMeanU
        procedure :: get_RHS
        procedure, private :: smear_this_source 
    end type


contains

subroutine init(this, inputDir, ActuatorDiskID, xG, yG, zG)
    class(actuatorDisk), intent(inout) :: this
    real(rkind), intent(in), dimension(:,:,:), target :: xG, yG, zG
    integer, intent(in) :: ActuatorDiskID
    character(len=*), intent(in) :: inputDir
    character(len=clen) :: tempname, fname
    integer :: ioUnit, tmpSum, totSum
    real(rkind) :: xLoc=1.d0, yLoc=1.d0, zLoc=0.1d0, diam=0.08d0, cT=0.65d0
    real(rkind) :: yaw=0.d0, tilt=0.d0, epsFact = 1.5d0, dx, dy, dz
    real(rkind), dimension(:,:), allocatable :: tmp
    integer, dimension(:,:), allocatable :: tmp_tag
    integer :: locator(1), ierr
    integer :: xLc(1), yLc(1), zLc(1)
    namelist /ACTUATOR_DISK/ xLoc, yLoc, zLoc, diam, cT, yaw, tilt
    
    ! Read input file for this turbine    
    write(tempname,"(A13,I3.3,A10)") "ActuatorDisk_", ActuatorDiskID, "_input.inp"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)

    ioUnit = 55
    open(unit=ioUnit, file=trim(fname), form='FORMATTED')
    read(unit=ioUnit, NML=ACTUATOR_DISK)
    close(ioUnit)
    
    this%xLoc = xLoc; this%yLoc = yLoc; this%zLoc = zLoc
    this%cT = cT; this%diam = diam; this%yaw = this%yaw
    dx=xG(2,1,1)-xG(1,1,1); dy=yG(1,2,1)-yG(1,1,1); dz=zG(1,1,2)-zG(1,1,1)
    this%nxLoc = size(xG,1); this%nyLoc = size(xG,2); this%nzLoc = size(xG,3)

    this%delta = epsFact * (dx*dy*dz)**(1.d0/3.d0)
    this%OneByDelSq = 1.d0/(this%delta**2)

    allocate(tmp(size(xG,2),size(xG,3)))
    allocate(tmp_tag(size(xG,2),size(xG,3)))
    allocate(this%tag_face(size(xG,2),size(xG,3)))
    allocate(this%xLine(size(xG,1)))
    allocate(this%yLine(size(xG,2)))
    allocate(this%zLine(size(xG,3)))
    
    this%xLine = xG(:,1,1); this%yLine = yG(1,:,1); this%zLine = zG(1,1,:)
    locator = minloc(abs(this%xLine - xLoc)); this%xLoc_idx = locator(1)

    tmp = sqrt((yG(1,:,:) - yLoc)**2 + (zG(1,:,:) - zLoc)**2)
    this%tag_face = 0
    tmp_tag = 0
    where(tmp < (diam/2.d0 + diam/2.d0))
        !this%tag_face = 1
        tmp_tag = 1
    end where
    where(tmp< (diam/2.d0))
            this%tag_face = 1
    end where

    if (sum(tmp_tag) > 0) then
        this%Am_I_Active = .true.
        this%color = ActuatorDiskID
        tmpSum = 1
    else
        this%Am_I_Active = .false.
        this%color = ActuatorDiskID*1000 
        tmpSum = 0
    end if 
    totSum = p_sum(tmpSum)
    if (totSum > 1) then
        this%Am_I_Split = .true.
    else
        this%Am_I_Split = .false. 
    end if 
    deallocate(tmp_tag)


    tmpSum = sum(this%tag_face)
    this%totPointsOnface = p_sum(tmpSum)
        
    this%pfactor = one/((this%delta**3)*(pi**(3.d0/2.d0)))
    if (this%Am_I_Split) then
        call MPI_COMM_SPLIT(mpi_comm_world, this%color, nrank, this%mycomm, ierr)
        call MPI_COMM_RANK( this%mycomm, this%myComm_nrank, ierr ) 
        call MPI_COMM_SIZE( this%mycomm, this%myComm_nproc, ierr )
    end if 

    
    if (this%Am_I_Active) then
        allocate(this%rbuff(size(xG,2),size(xG,3)))
        call sample_on_circle(diam/2.d0,yLoc,zLoc, this%ys,this%zs,ntry)
        allocate(this%xs(size(this%ys)))
        this%xs = xLoc
        this%normfactor = (1.d0/(real(size(this%xs),rkind)))*this%pfactor
        xLc = minloc(abs(xLoc - 1.5d0*diam/2.d0 - this%xline)); this%xst = max(1,xLc(1))
        yLc = minloc(abs(yLoc - diam - this%yline)); this%yst = max(1,yLc(1))
        zLc = minloc(abs(zLoc - diam - this%zline)); this%zst = max(1,zLc(1))

        xLc = minloc(abs(xLoc + 1.5d0*diam/2.d0 - this%xline)); this%xen = min(this%nxLoc,xLc(1))
        yLc = minloc(abs(yLoc + diam - this%yline)); this%yen = min(this%nyLoc,yLc(1))
        zLc = minloc(abs(zLoc + diam - this%zline)); this%zen = min(this%nzLoc,zLc(1))
        
        this%xlen=this%xen-this%xst+1
        this%ylen=this%yen-this%yst+1
        this%zlen=this%zen-this%zst+1
        
        allocate(this%dsq(this%xlen,this%ylen,this%zlen))    
        allocate(this%source(this%xlen,this%ylen,this%zlen))    
        allocate(this%eta_delta(this%xlen,this%ylen,this%zlen))    
        allocate(this%xSmall(this%xlen,this%ylen,this%zlen))    
        allocate(this%ySmall(this%xlen,this%ylen,this%zlen))    
        allocate(this%zSmall(this%xlen,this%ylen,this%zlen))    
        this%xSmall = xG(this%xst:this%xen,this%yst:this%yen,this%zst:this%zen)
        this%ySmall = yG(this%xst:this%xen,this%yst:this%yen,this%zst:this%zen)
        this%zSmall = zG(this%xst:this%xen,this%yst:this%yen,this%zst:this%zen)
    else
        deallocate(this%tag_face)
    end if 

    deallocate(tmp)
end subroutine 

subroutine destroy(this)
    class(actuatordisk), intent(inout) :: this

    if (Allocated(this%rbuff))  deallocate(this%rbuff)
    if (Allocated(this%tag_face))  deallocate(this%tag_face)
    
end subroutine 

subroutine getMeanU(this, u, v, w) 
    class(actuatordisk), intent(inout) :: this
    real(rkind), dimension(this%nxLoc,this%nyLoc,this%nzLoc), intent(in) :: u, v, w
    real(rkind) :: tmpSum, umn, vmn, wmn
    
    ! Get u face
    this%rbuff = u(this%xLoc_idx,:,:)
    this%rbuff = this%rbuff*this%tag_face
    if (this%AM_I_Split) then
        tmpSum = p_sum(sum(this%rbuff),this%myComm) 
        umn = tmpSum/real(this%totPointsOnFace,rkind)
    else
        umn = sum(this%rbuff)/real(this%totPointsOnFace,rkind)
    end if

    ! Get v face
    this%rbuff = v(this%xLoc_idx,:,:)
    this%rbuff = this%rbuff*this%tag_face
    if (this%AM_I_Split) then
        tmpSum = p_sum(sum(this%rbuff),this%myComm) 
        vmn = (tmpSum)/real(this%totPointsOnFace,rkind)
    else
        vmn = sum(this%rbuff)/real(this%totPointsOnFace,rkind)
    end if
    
    ! Get w face
    this%rbuff = w(this%xLoc_idx,:,:)
    this%rbuff = this%rbuff*this%tag_face
    if (this%AM_I_Split) then
        tmpSum = p_sum(sum(this%rbuff),this%myComm) 
        wmn = (tmpSum)/real(this%totPointsOnFace,rkind)
    else
        wmn = sum(this%rbuff)/real(this%totPointsOnFace,rkind)
    end if

    this%uface = this%alpha_tau*umn + (1.d0 - this%alpha_tau)*this%uface
    this%vface = this%alpha_tau*vmn + (1.d0 - this%alpha_tau)*this%vface
    this%wface = this%alpha_tau*wmn + (1.d0 - this%alpha_tau)*this%wface
    this%alpha_tau = alpha_smooth

end subroutine

subroutine get_RHS(this, u, v, w, rhsxvals, rhsyvals, rhszvals, inst_val)
    class(actuatordisk), intent(inout) :: this
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsxvals, rhsyvals, rhszvals
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(in)    :: u, v, w
    real(rkind), dimension(2),                                  intent(out)   :: inst_val
    integer :: j
    real(rkind) :: usp_sq, force

    if (this%Am_I_Active) then
        call this%getMeanU(u,v,w)
        usp_sq = this%uface**2 + this%vface**2 + this%wface**2
        force = this%normfactor * 0.5d0*this%cT*(pi*(this%diam**2)/4.d0)*usp_sq
        this%source = 0.d0
        do j =1,size(this%xs)
            this%dsq = (this%xSmall - this%xs(j))**2 + (this%ySmall - this%ys(j))**2 + (this%zSmall - this%zs(j))**2 
            this%dsq = exp(-this%dsq*this%oneByDelSq)
            this%eta_delta = (force) * this%dsq
            this%source = this%source + this%eta_delta
        end do
        rhsxvals(this%xst:this%xen,this%yst:this%yen,this%zst:this%zen) = this%source
        rhsyvals(this%xst:this%xen,this%yst:this%yen,this%zst:this%zen) = 0.d0
        rhszvals(this%xst:this%xen,this%yst:this%yen,this%zst:this%zen) = 0.d0
        inst_val(1) = force
        inst_val(2) = force*sqrt(usp_sq)
    end if 
end subroutine

subroutine smear_this_source(this, xC, yC, zC, valSource)
    class(actuatordisk), intent(inout) :: this
    real(rkind), intent(in) :: xC, yC, zC, valSource

    this%dsq = (this%xSmall - xC)**2 + (this%ySmall - yC)**2 + (this%zSmall - zC)**2 
    this%dsq = -this%dsq*this%oneByDelSq
    this%eta_delta = (valSource*this%normfactor) * exp(this%dsq)
    this%source = this%source + this%eta_delta

end subroutine 

subroutine sample_on_circle(R,xcen, ycen, xloc,yloc,np)
    use gridtools, only: linspace
    real(rkind), intent(in) :: R
    integer, intent(in) :: np
    real(rkind), intent(in) :: xcen, ycen
    integer, dimension(:), allocatable :: tag
    real(rkind), dimension(:), allocatable :: xline, yline
    real(rkind), dimension(:), allocatable, intent(out) :: xloc, yloc
    real(rkind), dimension(:), allocatable :: xtmp, ytmp, rtmp
    integer :: idx, i, j, nsz, iidx

    allocate(xline(np),yline(np))
    allocate(xtmp(np**2),ytmp(np**2), rtmp(np**2), tag(np**2))
    
    xline = linspace(-R,R,np+1)
    yline = linspace(-R,R,np+1)
    idx = 1
    do j = 1,np
        do i = 1,np
            xtmp(idx) = xline(i); ytmp(idx) = yline(j)
            idx = idx + 1
        end do 
    end do
    rtmp = sqrt(xtmp**2 + ytmp**2) 
    tag = 0
    where (rtmp < R) 
        tag = 1
    end where
    nsz = sum(tag)
    allocate(xloc(nsz), yloc(nsz))
    iidx = 1
    do idx = 1,size(tag)
        if (tag(idx) == 1) then
            xloc(iidx) = xtmp(idx)
            yloc(iidx) = ytmp(idx)
            iidx = iidx + 1
        end if
    end do

    xloc = xloc + xcen; yloc = yloc + ycen 
end subroutine

end module 
