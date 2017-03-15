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
    integer, parameter :: ntrymin = 10

    type :: CloudMod
       real(rkind) :: xTurbLoc, yTurbLoc
       real(rkind), allocatable, dimension(:,:,:) :: dsq, eta_delta, source, xSmall, ySmall, zSmall
    end type

    type :: actuatorDisk
        ! Actuator Disk Info
        integer :: xLoc_idx, ActutorDiskID
        integer, dimension(:,:), allocatable :: tag_face 
        real(rkind) :: yaw, tilt
        real(rkind) :: xLoc, yLoc, zLoc
        real(rkind) :: diam, cT, pfactor, normfactor, OneBydelSq
        real(rkind) :: uface = 0.d0, vface = 0.d0, wface = 0.d0
        integer :: totPointsOnFace
        !real(rkind), dimension(:,:,:), allocatable :: eta_delta, dsq, xSmall, ySmall, zSmall, source
        real(rkind), dimension(:), allocatable :: xs, ys, zs
        integer, allocatable, dimension(:) :: xst, xen, yst, yen, zst, zen, xlen, ylen, zlen, nlen
        ! Grid Info
        integer :: nxLoc, nyLoc, nzLoc 
        real(rkind) :: delta ! Smearing size
        real(rkind) :: alpha_tau = 1.d0! Smoothing parameter (set to 1 for initialization) 
        real(rkind), dimension(:,:), allocatable :: rbuff
        real(rkind), dimension(:), allocatable ::  xline, yline, zline

        type(CloudMod), allocatable, dimension(:) :: Cloud
        integer                                   :: max_num_clouds = 9

        ! MPI communicator stuff
        logical :: Am_I_Active, Am_I_Split
        integer :: color, myComm, myComm_nproc, myComm_nrank

    contains
        procedure :: init
        procedure :: destroy
        procedure, private :: getMeanU
        procedure :: get_RHS
        !procedure, private :: smear_this_source 
    end type


contains

subroutine init(this, inputDir, ActuatorDiskID, xG, yG, zG, gpC)
    class(actuatorDisk), intent(inout) :: this
    real(rkind), intent(in), dimension(:,:,:), target :: xG, yG, zG
    integer, intent(in) :: ActuatorDiskID
    character(len=*), intent(in) :: inputDir
    type(decomp_info), target :: gpC

    character(len=clen) :: tempname, fname, turboutfname
    integer :: ioUnit, tmpSum, totSum
    real(rkind) :: xLoc=1.d0, yLoc=1.d0, zLoc=0.1d0, diam=0.08d0, cT=0.65d0
    real(rkind) :: yaw=0.d0, tilt=0.d0, epsFact = 1.5d0, dx, dy, dz
    real(rkind) :: yLocGhLo=1.d0, yLocGhUp=1.d0, totProjRadius, Lx, Ly, xLocGhLo, xLocGhUp
    real(rkind) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(rkind), dimension(:,:), allocatable :: tmp,tmpGhlo,tmpGhUp
    integer, dimension(:,:), allocatable :: tmp_tag
    integer :: locator(1), ierr, stind, endind
    !integer :: xLc(1), yLc(1), zLc(1)
    integer :: icl
    logical :: periodicY = .true., periodicX = .true. ! hard coded for now
    integer :: ntry

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
    allocate(tmpGhLo(size(xG,2),size(xG,3)))
    allocate(tmpGhUp(size(xG,2),size(xG,3)))
    allocate(tmp_tag(size(xG,2),size(xG,3)))
    allocate(this%tag_face(size(xG,2),size(xG,3)))
    allocate(this%xLine(size(xG,1)))
    allocate(this%yLine(size(xG,2)))
    allocate(this%zLine(size(xG,3)))
    
    this%xLine = xG(:,1,1); this%yLine = yG(1,:,1); this%zLine = zG(1,1,:)
    locator = minloc(abs(this%xLine - xLoc)); this%xLoc_idx = locator(1)

    totProjRadius = diam/2.0_rkind + diam/2.0_rkind ! change the second diam/2.0 later
    tmp = sqrt((yG(1,:,:) - yLoc)**2 + (zG(1,:,:) - zLoc)**2)
    if(periodicX) then
      Lx = dx*gpC%xsz(1)
      xLocGhLo = xLoc - Lx; xLocGhUp = xLoc + Lx
    else
      xLocGhLo = xLoc; xLocGhUp = xLoc
    endif
    if(periodicY) then
      Ly = dy*gpC%ysz(2)
      yLocGhLo = yLoc - Ly
      tmpGhLo = sqrt((yG(1,:,:) - yLocGhLo)**2 + (zG(1,:,:) - zLoc)**2)

      yLocGhUp = yLoc + Ly
      tmpGhUp = sqrt((yG(1,:,:) - yLocGhLo)**2 + (zG(1,:,:) - zLoc)**2)
    else
      tmpGhLo = 10.0_rkind*totProjRadius
      tmpGhUp = 10.0_rkind*totProjRadius
    endif
    this%tag_face = 0
    tmp_tag = 0
    where((tmp < totProjRadius) .or. (tmpGhLo < totProjRadius) .or. (tmpGhUp < totProjRadius))
        !this%tag_face = 1
        tmp_tag = 1
    end where
    where((tmp < diam/2.d0) .or. (tmpGhLo < diam/2.0d0) .or. (tmpGhUp < diam/2.0d0))
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
    else
        this%myComm_nrank = -1
    end if 

    
    if (this%Am_I_Active) then
        ! Get the number of turbine points: 
        ntry = ntrymin!max(floor(diam/(sqrt(dy*dz))),ntrymin)
        allocate(this%rbuff(size(xG,2),size(xG,3)))
        call sample_on_circle(diam/2.d0,yLoc,zLoc, this%ys,this%zs,ntry)
        allocate(this%xs(size(this%ys)))
        this%xs = xLoc
        this%normfactor = (1.d0/(real(size(this%xs),rkind)))*this%pfactor
        write(turboutfname,'(a,i3.3,a)') 'TurbineInitLog_', ActuatorDiskID, '.out'
        ioUnit = 100+nrank
        open(unit=ioUnit, file=trim(turboutfname), form='FORMATTED')
        write(ioUnit,'(a,i3.3,1x,2(e19.12,1x),i5.5,1x,e19.12)') '----Turbine No.---', ActuatorDiskID, this%normfactor, this%pfactor, size(this%xs), this%delta

        allocate(this%Cloud(this%max_num_clouds))
        allocate(this%xst(this%max_num_clouds),  this%xen(this%max_num_clouds))
        allocate(this%yst(this%max_num_clouds),  this%yen(this%max_num_clouds))
        allocate(this%zst(this%max_num_clouds),  this%zen(this%max_num_clouds))
        allocate(this%xlen(this%max_num_clouds), this%ylen(this%max_num_clouds))
        allocate(this%zlen(this%max_num_clouds), this%nlen(this%max_num_clouds))
        do icl = 1, this%max_num_clouds
          if(icl==1) then
            this%Cloud(icl)%xTurbLoc = xLoc; this%Cloud(icl)%yTurbLoc = yLoc
          elseif(icl==2) then
            this%Cloud(icl)%xTurbLoc = xLoc; this%Cloud(icl)%yTurbLoc = yLocGhLo
          elseif(icl==3) then
            this%Cloud(icl)%xTurbLoc = xLoc; this%Cloud(icl)%yTurbLoc = yLocGhUp
          elseif(icl==4) then
            this%Cloud(icl)%xTurbLoc = xLocGhLo; this%Cloud(icl)%yTurbLoc = yLoc
          elseif(icl==5) then
            this%Cloud(icl)%xTurbLoc = xLocGhLo; this%Cloud(icl)%yTurbLoc = yLocGhLo
          elseif(icl==6) then
            this%Cloud(icl)%xTurbLoc = xLocGhLo; this%Cloud(icl)%yTurbLoc = yLocGhUp
          elseif(icl==7) then
            this%Cloud(icl)%xTurbLoc = xLocGhUp; this%Cloud(icl)%yTurbLoc = yLoc
          elseif(icl==8) then
            this%Cloud(icl)%xTurbLoc = xLocGhUp; this%Cloud(icl)%yTurbLoc = yLocGhLo
          elseif(icl==9) then
            this%Cloud(icl)%xTurbLoc = xLocGhUp; this%Cloud(icl)%yTurbLoc = yLocGhUp
          endif
          !xLc = minloc(abs(this%Cloud(icl)%xTurbLoc - 1.5d0*totProjRadius/2.d0 - this%xline)); this%xst(icl) = max(1,xLc(1))
          !yLc = minloc(abs(this%Cloud(icl)%yTurbLoc - totProjRadius - this%yline));            this%yst(icl) = max(1,yLc(1))
          !zLc = minloc(abs(zLoc     - totProjRadius - this%zline));            this%zst(icl) = max(1,zLc(1))

          !xLc = minloc(abs(this%Cloud(icl)%xTurbLoc + 1.5d0*totProjRadius/2.d0 - this%xline)); this%xen(icl) = min(this%nxLoc,xLc(1))
          !yLc = minloc(abs(this%Cloud(icl)%yTurbLoc + totProjRadius - this%yline));            this%yen(icl) = min(this%nyLoc,yLc(1))
          !zLc = minloc(abs(zLoc     + totProjRadius - this%zline));            this%zen(icl) = min(this%nzLoc,zLc(1))

          xmin = this%Cloud(icl)%xTurbLoc - 1.5d0*totProjRadius/2.0d0; xmax = this%Cloud(icl)%xTurbLoc + 1.5d0*totProjRadius/2.0d0
          ymin = this%Cloud(icl)%yTurbLoc - totProjRadius;             ymax = this%Cloud(icl)%yTurbLoc + totProjRadius
          zmin = zLoc - totProjRadius;             zmax = zLoc + totProjRadius

          call get_extents(xmin, xmax, xG(:,1,1), stind, endind); this%xst(icl) = stind; this%xen(icl) = endind
          call get_extents(ymin, ymax, yG(1,:,1), stind, endind); this%yst(icl) = stind; this%yen(icl) = endind
          call get_extents(zmin, zmax, zG(1,1,:), stind, endind); this%zst(icl) = stind; this%zen(icl) = endind
          
          this%xlen(icl)=this%xen(icl)-this%xst(icl)+1
          this%ylen(icl)=this%yen(icl)-this%yst(icl)+1
          this%zlen(icl)=this%zen(icl)-this%zst(icl)+1

          !if(this%Cloud(icl)%xTurbLoc - 1.5d0*totProjRadius/2.0d0 < this%xline(1))

          this%nlen(icl) = this%xlen(icl)*this%ylen(icl)*this%zlen(icl)
          if(this%nlen(icl) > 0) then
             allocate( this%Cloud(icl)%dsq   (this%xlen(icl),this%ylen(icl),this%zlen(icl)), this%Cloud(icl)%source(this%xlen(icl),this%ylen(icl),this%zlen(icl)), this%Cloud(icl)%eta_delta(this%xlen(icl),this%ylen(icl),this%zlen(icl)) )
             allocate( this%Cloud(icl)%xSmall(this%xlen(icl),this%ylen(icl),this%zlen(icl)), this%Cloud(icl)%ySmall(this%xlen(icl),this%ylen(icl),this%zlen(icl)), this%Cloud(icl)%zSmall   (this%xlen(icl),this%ylen(icl),this%zlen(icl)) )
             this%Cloud(icl)%xSmall = xG(this%xst(icl):this%xen(icl),this%yst(icl):this%yen(icl),this%zst(icl):this%zen(icl))
             this%Cloud(icl)%ySmall = yG(this%xst(icl):this%xen(icl),this%yst(icl):this%yen(icl),this%zst(icl):this%zen(icl))
             this%Cloud(icl)%zSmall = zG(this%xst(icl):this%xen(icl),this%yst(icl):this%yen(icl),this%zst(icl):this%zen(icl))
          endif
        enddo

        write(ioUnit,*) '-----New Turbine------' 
        write(ioUnit,*) 'xLoc, yLoc: ', this%Cloud(1)%xTurbLoc, this%Cloud(1)%yTurbLoc 
        do icl = 1, this%max_num_clouds
          write(ioUnit,*) '-------Cloud No.------' , icl
          write(ioUnit,*) 'nlen: ', this%nlen(icl)
          write(ioUnit,*) 'xTLoc, yTLoc: ', this%Cloud(icl)%xTurbLoc, this%Cloud(icl)%yTurbLoc
          write(ioUnit,'(a,4(e19.12,1x))') 'xLims: ', this%Cloud(icl)%xTurbLoc - 1.5d0*totProjRadius/2.0d0, this%Cloud(icl)%xTurbLoc + 1.5d0*totProjRadius/2.0d0,minval(xG(:,1,1)), maxval(xG(:,1,1))
          write(ioUnit,'(a,4(e19.12,1x))') 'yLims: ', this%Cloud(icl)%yTurbLoc - totProjRadius, this%Cloud(icl)%yTurbLoc + totProjRadius, minval(yG(1,:,1)), maxval(yG(1,:,1))
          write(ioUnit,'(a,4(e19.12,1x))') 'zLims: ', zLoc - totProjRadius, zLoc + totProjRadius, minval(zG(1,1,:)), maxval(zG(1,1,:))
          write(ioUnit,'(a,2(i5,1x))') 'xst: ', this%xst(icl), this%xen(icl)
          write(ioUnit,'(a,2(i5,1x))') 'yst: ', this%yst(icl), this%yen(icl)
          write(ioUnit,'(a,2(i5,1x))') 'zst: ', this%zst(icl), this%zen(icl)
          if(this%nlen(icl) > 0) then 
            write(ioUnit,'(a,2(e19.12,1x))') 'x: ', this%Cloud(icl)%xSmall(1,1,1), this%Cloud(icl)%xSmall(this%xlen(icl),1,1)
            write(ioUnit,'(a,2(e19.12,1x))') 'y: ', this%Cloud(icl)%ySmall(1,1,1), this%Cloud(icl)%ySmall(1,this%ylen(icl),1)
            write(ioUnit,'(a,2(e19.12,1x))') 'z: ', this%Cloud(icl)%zSmall(1,1,1), this%Cloud(icl)%zSmall(1,1,this%zlen(icl))
          endif
        enddo
        write(ioUnit,*) '-----Done Turbine------' 
        close(ioUnit)

        !allocate(this%dsq(this%xlen,this%ylen,this%zlen))    
        !allocate(this%source(this%xlen,this%ylen,this%zlen))    
        !allocate(this%eta_delta(this%xlen,this%ylen,this%zlen))    
        !allocate(this%xSmall(this%xlen,this%ylen,this%zlen))    
        !allocate(this%ySmall(this%xlen,this%ylen,this%zlen))    
        !allocate(this%zSmall(this%xlen,this%ylen,this%zlen))    
        !this%xSmall = xG(this%xst:this%xen,this%yst:this%yen,this%zst:this%zen)
        !this%ySmall = yG(this%xst:this%xen,this%yst:this%yen,this%zst:this%zen)
        !this%zSmall = zG(this%xst:this%xen,this%yst:this%yen,this%zst:this%zen)

    else
        deallocate(this%tag_face)
    end if 

    deallocate(tmp)
end subroutine 

subroutine destroy(this)
    class(actuatordisk), intent(inout) :: this
    integer :: icl

    if (Allocated(this%rbuff))  deallocate(this%rbuff)
    if (Allocated(this%tag_face))  deallocate(this%tag_face)
   
    do icl = 1, this%max_num_clouds
      if(this%nlen(icl) > 0) then
          deallocate(this%Cloud(icl)%dsq, this%Cloud(icl)%source, this%Cloud(icl)%eta_delta)
          deallocate(this%Cloud(icl)%xSmall, this%Cloud(icl)%ySmall, this%Cloud(icl)%zSmall)
      endif
    enddo
    deallocate(this%xst,this%yst,this%zst,this%xen,this%yen,this%zen,this%xlen,this%ylen,this%zlen,this%nlen,this%Cloud)
 
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
    class(actuatordisk), intent(inout), target :: this
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsxvals, rhsyvals, rhszvals
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(in)    :: u, v, w
    real(rkind), dimension(8),                                  intent(out), optional  :: inst_val
    integer :: j, icl
    real(rkind), dimension(:,:,:), pointer :: xCloudPtr, yCloudPtr, zCloudPtr, dsqPtr, eta_deltaPtr, sourcePtr
    real(rkind) :: usp_sq, force

    if (this%Am_I_Active) then
        call this%getMeanU(u,v,w)
        usp_sq = this%uface**2 !+ this%vface**2 + this%wface**2
        force = -this%normfactor * 0.5d0*this%cT*(pi*(this%diam**2)/4.d0)*usp_sq
        do icl = 1, this%max_num_clouds
          if(this%nlen(icl) < 1) cycle
          xCloudPtr => this%Cloud(icl)%xSmall; yCloudPtr    => this%Cloud(icl)%ySmall;    zCloudPtr => this%Cloud(icl)%zSmall
          dsqPtr    => this%Cloud(icl)%dsq;    eta_deltaPtr => this%Cloud(icl)%eta_delta; sourcePtr => this%Cloud(icl)%source
          sourcePtr = 0.d0
          do j =1,size(this%xs)
              dsqPtr = (xCloudPtr - this%xs(j))**2 + (yCloudPtr - this%ys(j))**2 + (zCloudPtr - this%zs(j))**2 
              dsqPtr = exp(-dsqPtr*this%oneByDelSq)
              eta_deltaPtr = (force) * dsqPtr
              sourcePtr = sourcePtr + eta_deltaPtr
          end do
          rhsxvals(this%xst(icl):this%xen(icl),this%yst(icl):this%yen(icl),this%zst(icl):this%zen(icl)) = sourcePtr
          rhsyvals(this%xst(icl):this%xen(icl),this%yst(icl):this%yen(icl),this%zst(icl):this%zen(icl)) = 0.d0
          rhszvals(this%xst(icl):this%xen(icl),this%yst(icl):this%yen(icl),this%zst(icl):this%zen(icl)) = 0.d0
          nullify(sourcePtr, eta_deltaPtr, dsqPtr, zCloudPtr, yCloudPtr, xCloudPtr)
        end do
        if (present(inst_val)) then
          if((this%Am_I_Split .and. this%myComm_nrank==0) .or. (.not. this%Am_I_Split)) then
            inst_val(1) = force
            inst_val(2) = force*sqrt(usp_sq)
            inst_val(3) = sqrt(usp_sq)
            inst_val(4) = usp_sq
            inst_val(5) = usp_sq*inst_val(3)
            inst_val(6) = this%uface
            inst_val(7) = this%vface
            inst_val(8) = this%wface
          end if
        end if 
    end if 
end subroutine

!subroutine smear_this_source(this, xC, yC, zC, valSource)
!    class(actuatordisk), intent(inout) :: this
!    real(rkind), intent(in) :: xC, yC, zC, valSource
!
!    this%dsq = (this%xSmall - xC)**2 + (this%ySmall - yC)**2 + (this%zSmall - zC)**2 
!    this%dsq = -this%dsq*this%oneByDelSq
!    this%eta_delta = (valSource*this%normfactor) * exp(this%dsq)
!    this%source = this%source + this%eta_delta
!
!end subroutine 

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

subroutine get_extents(xmin, xmax, procmesh, stind, endind)
    real(rkind), intent(in) :: xmin, xmax
    real(rkind), dimension(:), intent(in) :: procmesh
    integer, intent(out) :: stind, endind

    integer :: nloc, ind!, ist, iend

    nloc = size(procmesh)

    ! determine x extents
    if(xmax < procmesh(1)) then
       ! turbine i influence is exclusively to the left of this processor
       stind = 0; endind = -1;
    else if(xmin > procmesh(nloc)) then
       ! turbine i influence is exclusively to the right of this processor
       stind = 0; endind = -1;
    else
       ! turbine i influences this processor

       ! determine stind
       if(xmin < procmesh(1)) then
          stind = 1
       else
          stind = -1
          do ind = 1, nloc-1
            if(xmin >= procmesh(ind) .and. xmin < procmesh(ind+1)) exit
          enddo
          stind = ind+1
          if(ind == nloc) then
             write(*,*) 'Something wrong. Check details.', xmin, procmesh(1), procmesh(nloc)
             call GracefulExit("Exiting from min extents setup", 423)
          endif
       endif

       ! determine endind
       if(xmax > procmesh(nloc)) then
          endind = nloc
       else
          endind = -1
          do ind = 1, nloc-1
            if(xmax >= procmesh(ind) .and. xmax < procmesh(ind+1)) exit
          enddo
          endind = ind
          if(ind == nloc) then
             write(*,*) 'Something wrong. Check details.', xmax, procmesh(1), procmesh(nloc)
             call GracefulExit("Exiting from max extents setup", 423)
          endif
       endif
    endif
end subroutine


end module 
