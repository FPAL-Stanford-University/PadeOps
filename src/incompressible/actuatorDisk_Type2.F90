module actuatorDisk_T2mod
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
    public :: actuatorDisk_T2
    
    real(rkind), parameter :: alpha_Smooth = 0.9d0 ! Exonential smoothing constant
    integer, parameter :: xReg = 8, yReg = 8, zReg = 8

    type :: actuatorDisk_T2
        ! Actuator Disk_T2 Info
        integer :: xLoc_idx, ActutorDisk_T2ID
        integer, dimension(:,:), allocatable :: tag_face 
        real(rkind) :: yaw, tilt
        real(rkind) :: xLoc, yLoc, zLoc
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
        real(rkind), dimension(:,:), allocatable :: rbuff
        real(rkind), dimension(:), allocatable :: dline, xline, yline, zline
        real(rkind), dimension(:,:,:), pointer :: xG, yG, zG
        real(rkind), dimension(:,:,:), allocatable :: smearing_base

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

subroutine init(this, inputDir, ActuatorDisk_T2ID, xG, yG, zG)
    class(actuatorDisk_T2), intent(inout) :: this
    real(rkind), intent(in), dimension(:,:,:), target :: xG, yG, zG
    integer, intent(in) :: ActuatorDisk_T2ID
    character(len=*), intent(in) :: inputDir
    character(len=clen) :: tempname, fname
    integer :: ioUnit, tmpSum, totSum
    real(rkind) :: xLoc=1.d0, yLoc=1.d0, zLoc=0.1d0, diam=0.08d0, cT=0.65d0
    real(rkind) :: yaw=0.d0, tilt=0.d0, epsFact = 1.5d0, dx, dy, dz
    real(rkind), dimension(:,:), allocatable :: tmp
    integer, dimension(:,:), allocatable :: tmp_tag
    integer :: j, locator(1)
    !integer :: i, ylen, zlen
    integer :: xLc(1), yLc(1), zLc(1), xst, xen, yst, yen, zst, zen, ierr, xlen
    integer  :: ntry = 100
    real(rkind) :: time2initialize = 0, correction_factor, normfact_p

    namelist /ACTUATOR_DISK/ xLoc, yLoc, zLoc, diam, cT, yaw, tilt
    
    ! Read input file for this turbine    
    write(tempname,"(A13,I3.3,A10)") "ActuatorDisk_T2_", ActuatorDisk_T2ID, "_input.inp"
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
    allocate(this%dsq(2*xReg+1,2*yReg+1,2*zReg+1))
    allocate(this%dline(2*xReg+1))
    this%dsq = 0.d0
    this%xG => xG; this%yG => yG; this%zG => zG
    
    this%xLine = xG(:,1,1); this%yLine = yG(1,:,1); this%zLine = zG(1,1,:)
    locator = minloc(abs(this%xLine - xLoc)); this%xLoc_idx = locator(1)

    tmp = sqrt((yG(1,:,:) - yLoc)**2 + (zG(1,:,:) - zLoc)**2)
    this%tag_face = 0
    tmp_tag = 0
    where(tmp < (diam/2.d0 + max(xReg*dx,yReg*dy,zReg*dz)))
        !this%tag_face = 1
        tmp_tag = 1
    end where
    where(tmp< (diam/2.d0))
            this%tag_face = 1
    end where

    if (sum(tmp_tag) > 0) then
        this%Am_I_Active = .true.
        this%color = ActuatorDisk_T2ID
        tmpSum = 1
    else
        this%Am_I_Active = .false.
        this%color = ActuatorDisk_T2ID*1000 
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

     ntry = 2*ceiling(diam/min(dx, dy, dz))
    
    if (this%Am_I_Active) then
        allocate(this%rbuff(size(xG,2),size(xG,3)))
        call sample_on_circle(diam/2.d0,yLoc,zLoc, this%ys,this%zs,ntry)
        allocate(this%xs(size(this%ys)))
        this%xs = xLoc
        allocate(this%startEnds(7,size(this%xs)))
        this%startEnds = 0
        do j = 1,size(this%xs)
                xLc = minloc(abs(this%xs(j) - this%xline)); 
                yLc = minloc(abs(this%ys(j) - this%yline)); 
                zLc = minloc(abs(this%zs(j) - this%zline))

                xst = max(1, xLc(1) - xReg)
                yst = max(1, yLc(1) - yReg) 
                zst = max(1, zLc(1) - zReg)
                xen = min(this%nxLoc,xLc(1) + xReg) 
                yen = min(this%nyLoc,yLc(1) + yreg) 
                zen = min(this%nzLoc,zLc(1) + zreg)
                xlen = xen - xst + 1; !ylen = yen - yst + 1; zlen = zen - zst + 1
                this%startEnds(1,j) = xst; this%startEnds(2,j) = xen; this%startEnds(3,j) = yst
                this%startEnds(4,j) = yen; this%startEnds(5,j) = zst; this%startEnds(6,j) = zen
                this%startEnds(7,j) = xlen
        end do  
        this%normfactor = 1.d0/(real(size(this%xs),rkind))
    else
        deallocate(this%dsq, this%tag_face)
    end if 

    call message(1, "Initializing Actuator Disk (ADM Type=2) number", ActuatorDisk_T2ID)
    call tic()
    allocate(this%smearing_base(size(xG,1), size(xG,2), size(xG,3))) 
    this%smearing_base = 0.d0
    do j = 1,size(this%xs)
            call this%smear_this_source(this%smearing_base,this%xs(j),this%ys(j),this%zs(j), one, this%startEnds(1,j), &
                                this%startEnds(2,j),this%startEnds(3,j),this%startEnds(4,j), &
                                this%startEnds(5,j),this%startEnds(6,j),this%startEnds(7,j))
    end do 

    ! correction factor ::  required if this%smearing_base does not sum up to one because part of the cloud is outside the domain
    ! this can happen if turbine is close to the bottom wall, or, in the future, close to an immersed boundary
    normfact_p = p_maxval(this%normfactor)
    this%normfactor = normfact_p
    correction_factor = p_sum(this%smearing_base)*dx*dy*dz*this%pfactor*this%normfactor
    call message(2, "correction factor = ", correction_factor)
    !write(*,'(a,i4,e19.12,1x,e19.12)') '--', nrank, correction_factor, this%normfactor
    this%smearing_base = this%smearing_base/correction_factor

    call message(2, "Smearing grid parameter, ntry", ntry)
    call toc(mpi_comm_world, time2initialize)
    call message(2, "Time (seconds) to initialize", time2initialize)

    deallocate(tmp)
end subroutine 

subroutine destroy(this)
    class(actuatordisk_T2), intent(inout) :: this

    if (Allocated(this%rbuff))  deallocate(this%rbuff)
    if (Allocated(this%tag_face))  deallocate(this%tag_face)
    
    nullify(this%xG, this%yG, this%zG)
end subroutine 

subroutine getMeanU(this, u, v, w) 
    class(actuatordisk_T2), intent(inout) :: this
    real(rkind), dimension(this%nxLoc,this%nyLoc,this%nzLoc), intent(in) :: u, v, w
    real(rkind) :: tmpSum, umn, vmn, wmn
    
    if (this%AM_I_ACTIVE) then
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
    end if 
    this%alpha_tau = alpha_smooth

end subroutine

subroutine get_RHS(this, u, v, w, rhsxvals, rhsyvals, rhszvals)
!subroutine get_RHS(this, u, v, w, rhsvals)
    class(actuatordisk_T2), intent(inout) :: this
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsxvals, rhsyvals, rhszvals
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(in) :: u, v, w
    integer :: j
    real(rkind) :: usp_sq, force

    call this%getMeanU(u,v,w)
    usp_sq = this%uface**2 + this%vface**2 + this%wface**2
    force = -this%pfactor*this%normfactor*0.5d0*this%cT*(pi*(this%diam**2)/4.d0)*usp_sq
    !do j = 1,size(this%xs)
    !        call this%smear_this_source(rhsxvals,this%xs(j),this%ys(j),this%zs(j), force, this%startEnds(1,j), &
    !                            this%startEnds(2,j),this%startEnds(3,j),this%startEnds(4,j), &
    !                            this%startEnds(5,j),this%startEnds(6,j),this%startEnds(7,j))
    !end do 
    rhsxvals = rhsxvals + force*this%smearing_base 
    rhsyvals = zero
    rhszvals = zero

end subroutine

subroutine smear_this_source(this,rhsfield, xC, yC, zC, valSource, xst, xen, yst, yen, zst, zen, xlen)
    class(actuatordisk_T2), intent(inout) :: this
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsfield
    real(rkind), intent(in) :: xC, yC, zC, valSource
    integer, intent(in) :: xst, xen, yst, yen, zst, zen, xlen!, ylen, zlen
    integer :: jj, kk

    if (this%Am_I_Active) then
        do kk = zst,zen
            do jj = yst,yen
                this%dline(1:xlen) = (this%xG(xst:xen,jj,kk) - xC)**2 &
                                   + (this%yG(xst:xen,jj,kk) - yC)**2 & 
                                   + (this%zG(xst:xen,jj,kk) - zC)**2
                this%dline(1:xlen) = -this%dline(1:xlen)*this%oneBydelSq
                rhsfield(xst:xen,jj,kk) = rhsfield(xst:xen,jj,kk) + &
                                          valSource*exp(this%dline)

            end do 
        end do  
    end if 

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
