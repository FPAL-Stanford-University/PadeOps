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
    integer, parameter :: xReg = 8, yReg = 8, zReg = 8

    type :: actuatorDisk
        ! Actuator Disk Info
        integer :: xLoc_idx, ActutorDiskID
        integer(rkind), dimension(:,:), allocatable :: tag_face 
        real(rkind) :: yaw, tilt
        real(rkind) :: xLoc, yLoc, zLoc
        real(rkind) :: diam, cT, pfactor, normfactor
        real(rkind) :: uface = 0.d0, vface = 0.d0, wface = 0.d0
        integer :: totPointsOnFace
        real(rkind), dimension(:,:,:), allocatable :: eta_delta, dsq
        real(rkind), dimension(:,:), allocatable :: xp, yp, zp

        ! Grid Info
        integer :: nxLoc, nyLoc, nzLoc 
        real(rkind) :: delta ! Smearing size
        real(rkind) :: alpha_tau = 1.d0! Smoothing parameter (set to 1 for initialization) 
        real(rkind), dimension(:,:), allocatable :: rbuff
        real(rkind), dimension(:), allocatable :: xline, yline, zline
        real(rkind), dimension(:,:,:), pointer :: xG, yG, zG

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
    integer :: ioUnit, tmpSum, totSum, ierr
    integer :: n_rad = 10, n_theta = 20
    real(rkind) :: xLoc=1.d0, yLoc=1.d0, zLoc=0.1d0, diam=0.08d0, cT=0.65d0
    real(rkind) :: yaw=0.d0, tilt=0.d0, epsFact = 1.5d0, dx, dy, dz
    real(rkind), dimension(:,:), allocatable :: tmp
    integer, dimension(:,:), allocatable :: tmp_tag
    real(rkind), dimension(:), allocatable :: rEdge, rCell
    integer :: i, j, locator(1)
    namelist /ACTUATOR_DISK/ xLoc, yLoc, zLoc, diam, cT, yaw, tilt, epsFact, &
                             & n_rad, n_theta
    
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

    allocate(tmp(size(xG,2),size(xG,3)))
    allocate(tmp_tag(size(xG,2),size(xG,3)))
    allocate(this%tag_face(size(xG,2),size(xG,3)))
    allocate(this%xLine(size(xG,1)))
    allocate(this%yLine(size(xG,2)))
    allocate(this%zLine(size(xG,3)))
    allocate(this%dsq(2*xReg+1,2*yReg+1,2*zReg+1))
    this%dsq = 0.d0
    this%xG => xG; this%yG => yG; this%zG => zG
    
    this%xLine = xG(:,1,1); this%yLine = yG(1,:,1); this%zLine = zG(1,1,:)
    locator = minloc(abs(this%xLine - xLoc)); this%xLoc_idx = locator(1)

    tmp = sqrt((yG(1,:,:) - yLoc)**2 + (zG(1,:,:) - zLoc)**2)
    this%tag_face = 0
    tmp_tag = 0
    where(tmp < (diam/2.d0 + 8.d0*max(dz,dy,dz)))
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
        allocate(rEdge(n_rad+1), rCell(n_rad))
        rEdge = linspace(0.d0,diam/2.d0,n_rad+1)
        rCell = 0.5d0*(rEdge(1:n_rad) + rEdge(2:n_rad+1))
        allocate(this%xp(n_rad, n_theta), this%yp(n_rad, n_theta), this%zp(n_rad, n_theta))

        do j = 1,n_theta
            do i = 1,n_rad
                this%xp(i,j) = xLoc
                this%yp(i,j) = yLoc + rCell(i)*cos((j-1)*2.d0*pi/real(n_theta))
                this%zp(i,j) = zLoc + rCell(i)*sin((j-1)*2.d0*pi/real(n_theta))
            end do 
        end do
        this%normfactor = 1.d0/(real(n_theta)*real(n_rad)) 
    else
        deallocate(this%dsq, this%tag_face)
    end if 

    deallocate(tmp)
end subroutine 

subroutine destroy(this)
    class(actuatordisk), intent(inout) :: this

    if (Allocated(this%rbuff))  deallocate(this%rbuff)
    if (Allocated(this%tag_face))  deallocate(this%tag_face)
    
    nullify(this%xG, this%yG, this%zG)
end subroutine 

subroutine getMeanU(this, u, v, w) 
    class(actuatordisk), intent(inout) :: this
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

subroutine get_RHS(this, u, v, w, rhsvals)
    class(actuatordisk), intent(inout) :: this
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsvals
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(in) :: u, v, w
    integer :: i, j, ierr
    real(rkind) :: usp_sq, force

    call this%getMeanU(u,v,w)
    usp_sq = this%uface**2 + this%vface**2 + this%wface**2
    force = this%normfactor*0.5d0*this%cT*(pi*(this%diam**2)/4.d0)*usp_sq
    do j = 1,size(this%xP,2)
        do i = 1,size(this%xP,1)
            call this%smear_this_source(rhsvals,this%xp(i,j),this%yp(i,j),this%zp(i,j), force)
        end do  
    end do 
end subroutine

subroutine smear_this_source(this,rhsfield, xC, yC, zC, valSource)
    class(actuatordisk), intent(inout) :: this
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsfield
    real(rkind), intent(in) :: xC, yC, zC, valSource
    integer :: xLoc(1), yLoc(1), zLoc(1), xst, xen, yst, yen, zst, zen, ierr, xlen, ylen, zlen

    if (this%Am_I_Active) then
        xLoc = minloc(abs(xC - this%xline)); 
        yLoc = minloc(abs(yC - this%yline)); 
        zLoc = minloc(abs(zC - this%zline))

        xst = max(1, xLoc(1) - xReg)
        yst = max(1, yLoc(1) - yreg) 
        zst = max(1, zLoc(1) - zreg)
        xen = min(this%nxLoc,xLoc(1) + xReg) 
        yen = min(this%nyLoc,yLoc(1) + yreg) 
        zen = min(this%nzLoc,zLoc(1) + zreg)
        xlen = xen - xst + 1; ylen = yen - yst + 1; zlen = zen - zst + 1

        this%dsq(1:xlen,1:ylen,1:zlen) = (this%xG(xst:xen,yst:yen,zst:zen) - xC)**2  &
                                       + (this%yG(xst:xen,yst:yen,zst:zen) - yC)**2  &
                                       + (this%zG(xst:xen,yst:yen,zst:zen) - zC)**2

        this%dsq = -this%dsq/(this%delta**2)
        rhsfield(xst:xen,yst:yen,zst:zen) = rhsfield(xst:xen,yst:yen,zst:zen) + &
                                            (this%pfactor*valSource)*exp(this%dsq(1:xlen,1:ylen,1:zlen))
      

    end if 

end subroutine 

end module 
