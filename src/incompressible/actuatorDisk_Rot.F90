module actuatorDisk_Rotmod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa, twopi, rpm_to_radpersec, four, third, deg_to_radians
    use decomp_2d
    use exits, only: GracefulExit, message
    use mpi 
    use reductions, only: p_maxval, p_sum, p_minval
    use timer, only: tic, toc
    use Gridtools, only: linspace

    implicit none

    private
    public :: actuatorDisk_Rot
    
    real(rkind), parameter :: alpha_Smooth = 1.d0 ! 0.9d0 ! Exonential smoothing constant
    integer, parameter :: xReg = 8, yReg = 8, zReg = 8

    type :: actuatorDisk_Rot
        ! Actuator Disk_Rot Info
        integer :: xface_idx, xdisk_idx, ActuatorDisk_RotID, diagnostic_counter, NacelleRadInd, diagnostics_write_freq=100
        integer, dimension(:,:), allocatable :: tag_face 
        real(rkind) :: yaw, tilt, TSR, Omega, Lref
        real(rkind) :: xLoc, yLoc, zLoc, diam, blade_pitch, NacelleRad, DiskArea, CTp_eff
        real(rkind) :: pfactor, OneBydelSq, normfactor, sampleUpsDist
        real(rkind) :: uface = 0.d0, vface = 0.d0, wface = 0.d0 ! goes with xface_idx - can be upstream of disk
        real(rkind) :: udisk = 0.d0, vdisk = 0.d0, wdisk = 0.d0 ! goes with xdisk_idx - always closest to disk
        integer :: totPointsOnFace, RPMController=0, num_radial_elems=40, num_azimut_elems=40, num_blades=3
        real(rkind), dimension(:,:,:), allocatable :: eta_delta, dsq, clTable, cdTable
        real(rkind), dimension(:,:),   allocatable :: xp, yp, zp, diagnosticarr
        real(rkind), dimension(:),     allocatable :: xs, ys, zs, radius_elem, cos_azimut, sin_azimut
        real(rkind), dimension(:),     allocatable :: twist_elem, solidity_elem, chord_elem
        integer,     dimension(:,:),   allocatable :: startEnds
        integer,     dimension(:),     allocatable :: startEnds_global, airfoilID_elem, tableID_elem

        ! Grid Info
        integer :: nxLoc, nyLoc, nzLoc 
        real(rkind) :: delta ! Smearing size
        real(rkind) :: alpha_tau = 1.d0! Smoothing parameter (set to 1 for initialization) 
        real(rkind), dimension(:,:),     allocatable :: rbuff
        real(rkind), dimension(:),       allocatable :: dline, xline, yline, zline, force_x, force_y, force_z
        real(rkind), dimension(:,:,:),   pointer     :: xG, yG, zG
        real(rkind), dimension(:,:,:,:), allocatable :: smearing_base
        real(rkind), dimension(:,:,:),   allocatable :: force_local
        logical :: isFan = .false.

        ! MPI communicator stuff
        logical :: Am_I_Active, Am_I_Split
        integer :: color, myComm, myComm_nproc, myComm_nrank

    contains
        procedure :: init
        procedure :: destroy
        procedure :: get_RHS
        procedure, private :: getMeanU
        procedure, private :: smear_this_source 
        procedure, private :: sample_on_circle_rtheta
        procedure, private :: get_induction_factors
        procedure, private :: get_Omega
        procedure, private :: get_cl_cd
        procedure, private :: read_blade_properties
        procedure, private :: read_airfoil_properties
    end type


contains

subroutine init(this, inputDir, ActuatorDisk_RotID, xG, yG, zG, ntryparam)
    class(actuatorDisk_Rot), intent(inout) :: this
    real(rkind), intent(in), dimension(:,:,:), target :: xG, yG, zG
    integer, intent(in) :: ActuatorDisk_RotID
    integer, intent(in), optional :: ntryparam
    character(len=*), intent(in) :: inputDir
    character(len=clen) :: tempname, fname
    integer :: ioUnit, tmpSum, totSum, RPMController=0, num_blades=3, bladeTypeID=1, diagnostics_write_freq=100
    real(rkind) :: xLoc=1.d0, yLoc=1.d0, zLoc=0.1d0, diam=0.08d0
    real(rkind) :: TSR=8.d0, Omega = 1120.d0, blade_pitch, reference_length, NacelleRad
    real(rkind) :: yaw=0.d0, tilt=0.d0, epsFact = 1.5d0, dx, dy, dz
    real(rkind), dimension(:,:), allocatable :: tmp
    integer,     dimension(:,:), allocatable :: tmp_tag
    real(rkind), dimension(:),   allocatable :: corrfacarr
    integer :: i, j, locator(1)
    integer :: xLc(1), yLc(1), zLc(1), xst, xen, yst, yen, zst, zen, ierr, xlen
    integer  :: ntry = 100, xsize, ysize, zsize
    real(rkind) :: time2initialize = 0, correction_factor = one, sampleUpsDist = zero
    logical :: isFan = .false.

    namelist /ACTUATOR_DISK/ xLoc, yLoc, zLoc, diam, RPMController, TSR, Omega, yaw, tilt, blade_pitch, num_blades, reference_length, NacelleRad, bladeTypeID, sampleUpsDist, diagnostics_write_freq, isFan
    
    ! Read input file for this turbine    
    write(tempname,"(A17,I4.4,A10)") "ActuatorDisk_Rot_", ActuatorDisk_RotID, "_input.inp"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)

    ioUnit = 55
    open(unit=ioUnit, file=trim(fname), form='FORMATTED', action="read")
    read(unit=ioUnit, NML=ACTUATOR_DISK)
    close(ioUnit)
    
    this%xLoc = xLoc;  this%yLoc = yLoc;  this%zLoc       = zLoc
    this%TSR  = TSR;   this%diam = diam;  this%yaw        = this%yaw
    this%RPMController = RPMController;   this%tilt       = tilt
    this%blade_pitch   = blade_pitch;     this%num_blades = num_blades
    this%NacelleRad    = NacelleRad;      this%Omega      = Omega*rpm_to_radpersec
    this%sampleUpsDist = sampleUpsDist;   this%diagnostics_write_freq = diagnostics_write_freq
    this%Lref = reference_length;         this%isFan      = isFan

    dx=xG(2,1,1)-xG(1,1,1); dy=yG(1,2,1)-yG(1,1,1); dz=zG(1,1,2)-zG(1,1,1)
    this%nxLoc = size(xG,1); this%nyLoc = size(xG,2); this%nzLoc = size(xG,3)

    this%delta = epsFact * (dx*dy*dz)**(1.d0/3.d0)
    this%OneByDelSq = 1.d0/(this%delta**2)
    this%ActuatorDisk_RotID = ActuatorDisk_RotID

    this%DiskArea = pi * fourth * this%diam**2

    !do j = 0, nproc-1
    !  call mpi_barrier(mpi_comm_world, ierr)
    !  if(nrank==j) then
    !     write(*,*) '---------------------------------------------------'
    !     write(*,*) 'processor rank =', nrank
    !     write(*,*) 'nx, ny, nz : ', this%nxLoc, this%nyLoc, this%nzLoc
    !     write(*,*) 'xG         : ', minval(xG), maxval(xG)
    !     write(*,*) 'yG         : ', minval(yG), maxval(yG)
    !     write(*,*) 'zG         : ', minval(zG), maxval(zG)
    !     write(*,*) '---------------------------------------------------'
    !  endif
    !  call mpi_barrier(mpi_comm_world, ierr)
    !enddo

    allocate(tmp(size(xG,2),size(xG,3)))
    allocate(tmp_tag(size(xG,2),size(xG,3)))
    allocate(this%tag_face(size(xG,2),size(xG,3)))
    allocate(this%xLine(size(xG,1)))
    allocate(this%yLine(size(xG,2)))
    allocate(this%zLine(size(xG,3)))
    !allocate(this%dsq(2*xReg+1,2*yReg+1,2*zReg+1))
    allocate(this%dline(2*xReg+1))
    !this%dsq = 0.d0
    this%xG => xG; this%yG => yG; this%zG => zG
    
    this%xLine = xG(:,1,1); this%yLine = yG(1,:,1); this%zLine = zG(1,1,:)
    locator = minloc(abs(this%xLine - (xLoc-this%sampleUpsDist*this%diam))); this%xface_idx = locator(1)
    locator = minloc(abs(this%xLine - (xLoc)));                              this%xdisk_idx = locator(1)

    tmp = sqrt((yG(1,:,:) - yLoc)**2 + (zG(1,:,:) - zLoc)**2)
    this%tag_face = 0
    tmp_tag = 0
    where(tmp < (diam/2.d0 + max(0*xReg*dx,yReg*dy,zReg*dz)))
        !this%tag_face = 1
        tmp_tag = 1
    end where
    where(tmp< (diam/2.d0))
            this%tag_face = 1
    end where

    if (sum(tmp_tag) > 0) then
        this%Am_I_Active = .true.
        this%color = ActuatorDisk_RotID
        tmpSum = 1
    else
        this%Am_I_Active = .false.
        this%color = ActuatorDisk_RotID*1000 
        tmpSum = 0
    end if 
    totSum = p_sum(tmpSum)
    if (totSum > 1) then
        this%Am_I_Split = .true.
    else
        this%Am_I_Split = .false. 
    end if 
    !write(*,'(a,i4,1x,L1,1x,i7,1x,f12.6,1x,3(i5,1x,f12.6))') '---TTTT----', nrank, this%Am_I_Active, sum(tmp_tag), diam, xReg, dx, yReg, dy, zReg, dz
    deallocate(tmp_tag)


    tmpSum = sum(this%tag_face)
    this%totPointsOnface = p_sum(tmpSum)

    this%pfactor = one/((this%delta**3)*(pi**(3.d0/2.d0)))

    if (this%Am_I_Split) then
        call MPI_COMM_SPLIT(mpi_comm_world, this%color, nrank, this%mycomm, ierr)
        call MPI_COMM_RANK( this%mycomm, this%myComm_nrank, ierr ) 
        call MPI_COMM_SIZE( this%mycomm, this%myComm_nproc, ierr )
    end if 

    !print '(a,i4,1x,L2,1x,L2,a)', '((((', nrank, this%Am_I_Split, this%Am_I_Active, ')))))'

    if (present(ntryparam)) then
        ntry = ntryparam*ceiling(diam/min(dx, dy, dz))
    else
        ntry = 2*ceiling(diam/min(dx, dy, dz))
    endif
    ntry = p_maxval(ntry) ! prevents mismatch across processors due to roundoff

    if (this%Am_I_Active) then
        allocate(this%rbuff(size(xG,2),size(xG,3)))
        call this%sample_on_circle_rtheta(diam/2.d0,xLoc,yLoc,zLoc,ntry)
        allocate(this%startEnds(7,size(this%xs)))
        allocate(this%startEnds_global(7))
        this%startEnds = 0
        do j = 1,size(this%xs)
            xLc = minloc(abs(this%xs(j) - this%xline)); ! location of turbine
            yLc = minloc(abs(this%ys(j) - this%yline)); 
            zLc = minloc(abs(this%zs(j) - this%zline))

            xst = max(1, xLc(1) - xReg) ! grid starting point for turbine
            yst = max(1, yLc(1) - yReg) 
            zst = max(1, zLc(1) - zReg)
            xen = min(this%nxLoc,xLc(1) + xReg) 
            yen = min(this%nyLoc,yLc(1) + yreg) 
            zen = min(this%nzLoc,zLc(1) + zreg)
            xlen = xen - xst + 1; !ylen = yen - yst + 1; zlen = zen - zst + 1
            this%startEnds(1,j) = xst; this%startEnds(2,j) = yst; this%startEnds(3,j) = zst
            this%startEnds(4,j) = xen; this%startEnds(5,j) = yen; this%startEnds(6,j) = zen
            this%startEnds(7,j) = xlen
        end do
        this%startEnds_global(1:3) = minval(this%startEnds(1:3,:),2)
        this%startEnds_global(4:6) = maxval(this%startEnds(4:6,:),2)
        this%startEnds_global(7)   = this%startEnds_global(4) - this%startEnds_global(1) + 1

        this%normfactor = 1.d0/(real(this%num_azimut_elems,rkind))
    else
        !if(allocated(this%dsq))      deallocate(this%dsq)
        deallocate(this%tag_face)
    end if 

    call message(1, "Initializing Rotating Actuator Disk (ADM Type=3) number", ActuatorDisk_RotID)
    call message(2, "Sampling velocity at ", this%xLine(this%xface_idx))
    call tic()
    if(this%Am_I_Active) then
        !print *, '--', nrank, '---yy'
        xsize = this%startEnds_global(4)-this%startEnds_global(1) + 1
        ysize = this%startEnds_global(5)-this%startEnds_global(2) + 1
        zsize = this%startEnds_global(6)-this%startEnds_global(3) + 1
        allocate(this%smearing_base(xsize,ysize,zsize,size(this%xs)))
        allocate(this%force_local(xsize,ysize,zsize))
        allocate(this%force_x(size(this%xs)))
        allocate(this%force_y(size(this%xs)))
        allocate(this%force_z(size(this%xs)))
        this%smearing_base = 0.d0
        !do j = 0, nproc-1
        !  call mpi_barrier(this%myComm, ierr)
        !  if(nrank==j) then
        !    write(*,*) '-----------------------------------------------------------'
        !    write(*,*) 'processor rank =', nrank
        !    write(*,'(a,2(e19.12,1x))') 'smearing_base : ', maxval(this%smearing_base), minval(this%smearing_base)
        !    write(*,*) 'maxloc        : ', maxloc(this%smearing_base)
        !    write(*,*) 'minloc        : ', minloc(this%smearing_base)
        !    write(*,'(a,4(e19.12,1x))') 'dx dy dz pfac : ', dx, dy, dz, this%pfactor
        !    write(*,*) 'ist ien       : ', this%startEnds_global(1), this%startEnds_global(4), xsize
        !    write(*,*) 'jst jen       : ', this%startEnds_global(2), this%startEnds_global(5), ysize
        !    write(*,*) 'kst ken       : ', this%startEnds_global(3), this%startEnds_global(6), zsize
        !    write(*,*) '-----------------------------------------------------------'
        !    !write(*,*) 'correctionfac : ', correction_factor
        !  endif
        !  call mpi_barrier(this%myComm, ierr)
        !enddo
        do j = 1,size(this%xs)
            call this%smear_this_source(this%smearing_base(:,:,:,j),this%xs(j),this%ys(j),this%zs(j),one)
        end do
        !print *, '-------', nrank, '---yy'
    else
        !print *, '++', nrank, '+++++'
        ! makes it easier for reductions if this is allocated
        allocate(this%xs(ntry*ntry))
        allocate(this%smearing_base(1,1,1,size(this%xs)))
        allocate(this%force_local(1,1,1))
        allocate(this%force_x(1))
        allocate(this%force_y(1))
        allocate(this%force_z(1))

        ! initialize all arrays
        this%smearing_base = zero; this%force_local = zero;  this%force_x = zero;
        this%force_y = zero;  this%force_z = zero;
        !print *, '+++++', nrank, this%smearing_base(1,1,1,size(this%xs)), '+++++'
    endif
    ! correction factor :: required if this%smearing_base doesn't sum up to one 
    ! this can happen if part of the cloud is outside the domain, which can happen if
    ! turbine is close to the bottom wall, or, in the future, close to an immersed boundary
    ! note :: cannot be clubbed in the upper if block
    allocate(corrfacarr(size(this%xs)))
    corrfacarr = 0.0_rkind
    do j = 1,size(this%xs)
        correction_factor = p_sum(this%smearing_base(:,:,:,j),this%myComm)*dx*dy*dz*this%pfactor
        !write(*,'(a,i4,e19.12,1x,e19.12)') '--', nrank, correction_factor, this%normfactor
        this%smearing_base(:,:,:,j) = this%smearing_base(:,:,:,j)/correction_factor
        corrfacarr(j) = correction_factor
    enddo
    correction_factor = p_maxval(maxval(corrfacarr))
    
    call message(2, "correction factor max ", correction_factor)
    correction_factor = p_minval(minval(corrfacarr))
    call message(2, "correction factor min ", correction_factor)
    deallocate(corrfacarr)

    !if(nrank==5) then
    !   write(*,*) 'processor rank =', nrank
    !   write(*,*) 'smearing_base : ', maxval(this%smearing_base), minval(this%smearing_base)
    !   write(*,*) 'maxloc        : ', maxloc(this%smearing_base)
    !   write(*,*) 'minloc        : ', minloc(this%smearing_base)
    !   write(*,*) 'dx dy dz pfac : ', dx, dy, dz, this%pfactor
    !   write(*,*) 'correctionfac : ', correction_factor
    !endif

    call message(2, "Smearing grid parameter, ntry", ntry)

    this%NacelleRadInd = 0

    if(this%Am_I_Active) then
        ! setup twist, solidity, cl and cd tables
        allocate(this%twist_elem(this%num_radial_elems))
        allocate(this%chord_elem(this%num_radial_elems))
        allocate(this%solidity_elem(this%num_radial_elems))
        allocate(this%airfoilID_elem(this%num_radial_elems))
        allocate(this%tableID_elem(this%num_radial_elems))

        ! Read information about blade type used in this turbine
        ! This file does not use namelists, so order in which data is written in this file appears is important
        write(tempname,"(A10,I3.3,A10)") "BladeType_", bladeTypeID, "_input.inp"
        fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)

        call this%read_blade_properties(InputDir, fname, reference_length)
    endif

    call message(2, "Nacelle Radius Index ", p_maxval(this%NacelleRadInd))

    call toc(mpi_comm_world, time2initialize)
    call message(2, "Time (seconds) to initialize", time2initialize)

    deallocate(tmp)

    ! diagnostics
    if(this%Am_I_Active .and. ((this%Am_I_Split .and. this%myComm_nrank==0) .or. (.not. this%Am_I_Split))) then
         ioUnit = 100; write(tempname,"(a6,i4.4,a21)") 'ADRot_', ActuatorDisk_RotID, '_diagnostic_setup.dat'
         fname = trim(tempname)
         open(ioUnit,file=fname,form='formatted',action='write',status='unknown')
         write(ioUnit,*) 'VARIABLES="V1","V2","V3","V4","V5"'
         write(ioUnit,'(a,i5)') 'ZONE T="radius, twist, chord, airfoil", F=POINT, I=', this%num_radial_elems
         do j = 1, this%num_radial_elems
           write(ioUnit,'(4(e19.12,1x),i5)') this%radius_elem(j), this%twist_elem(j), this%chord_elem(j), this%solidity_elem(j), this%airfoilID_elem(j)
         enddo
         do i = 1, size(this%clTable,3)
           write(ioUnit,'(a,i5,a,i5)') 'ZONE T="clTable_',i,'", F=POINT, I=', size(this%clTable,1)
           do j = 1, size(this%clTable,1)
             write(ioUnit,'(5(e19.12,1x))') this%clTable(j,1:2,i),zero,zero,zero
           enddo
         enddo
         do i = 1, size(this%cdTable,3)
           write(ioUnit,'(a,i5,a,i5)') 'ZONE T="cdTable_',i,'", F=POINT, I=', size(this%cdTable,1)
           do j = 1, size(this%cdTable,1)
             write(ioUnit,'(5(e19.12,1x))') this%cdTable(j,1:2,i),zero,zero,zero
           enddo
         enddo
         close(ioUnit)
    endif

    allocate(this%diagnosticarr(10,this%num_radial_elems))
    this%diagnostic_counter = 0

end subroutine 

subroutine read_blade_properties(this, InputDir, fname, reference_length)
    class(actuatordisk_Rot), intent(inout) :: this
    character(len=*),    intent(in) :: InputDir
    character(len=clen), intent(in) :: fname
    real(rkind), intent(in) :: reference_length

    integer :: ioUnit, num_table_entries, k
    real(rkind), allocatable, dimension(:,:) :: twistAngTable, chordTable, airfoilIDTable
    integer,     allocatable, dimension(:)   :: unique_airfoilIDs

    ioUnit = 55
    open(unit=ioUnit, file=trim(fname), form='FORMATTED')
      ! first table for twist angle
      read(ioUnit, *); read(ioUnit, *) 
      read(ioUnit, *) num_table_entries
      allocate(twistAngTable(num_table_entries,2))
      do k = 1, num_table_entries
        read(ioUnit, *) twistAngTable(k,:)
      enddo
      ! normalize radial distance by reference length
      twistAngTable(:,1) = twistAngTable(:,1)/reference_length
      ! convert twist angle from degrees to radians
      twistAngTable(:,2) = twistAngTable(:,2) * deg_to_radians

      ! second table for chord
      read(ioUnit, *); read(ioUnit, *)
      read(ioUnit, *) num_table_entries
      allocate(chordTable(num_table_entries,2))
      do k = 1, num_table_entries
        read(ioUnit, *) chordTable(k,:)
      enddo
      ! normalize radial distance and chord length by reference length
      chordTable = chordTable/reference_length

      ! third table for airfoilIDs
      read(ioUnit, *); read(ioUnit, *)
      read(ioUnit, *) num_table_entries
      allocate(airfoilIDTable(num_table_entries,2))
      do k = 1, num_table_entries
        read(ioUnit, *) airfoilIDTable(k,:)
      enddo
      ! normalize radial distance by reference length
      airfoilIDTable(:,1) = airfoilIDTable(:,1)/reference_length
    close(ioUnit)


    call interp(this%radius_elem, this%twist_elem, twistAngTable)
    call interp(this%radius_elem, this%chord_elem, chordTable   )
    call interpint(this%radius_elem, this%airfoilID_elem, airfoilIDTable)

    call find_unique_entries(this%airfoilID_elem, unique_airfoilIDs)

    do k = 1, this%num_radial_elems
        if(this%radius_elem(k) < this%NacelleRad) then
          this%NacelleRadInd = k
        endif
    enddo

    this%solidity_elem = this%num_blades*this%chord_elem/twopi/this%radius_elem
    !where(this%radius_elem<this%NacelleRad)
    !    this%solidity_elem = one
    !end where
    do k = 1, this%num_radial_elems
        if(this%radius_elem(k) < this%NacelleRad) then
          this%solidity_elem(k) = one
        endif
    enddo

    deallocate(twistAngTable, chordTable, airfoilIDTable)

    call this%read_airfoil_properties(InputDir, unique_airfoilIDs)

end subroutine 

subroutine read_airfoil_properties(this, InputDir, airfoilIDs)
    class(actuatordisk_Rot), intent(inout) :: this
    character(len=*),        intent(in)    :: InputDir
    integer, dimension(:),   intent(in)    :: airfoilIDs

    character(len=clen) :: tempname, fname
    integer, allocatable, dimension(:) :: num_entries_cl, num_entries_cd
    integer :: i, j, size_table, ioUnit

    allocate(num_entries_cl(size(airfoilIDs)))
    allocate(num_entries_cd(size(airfoilIDs)))

    ! read number of entries in each airfoil data
    do i = 1, size(airfoilIDs)
        ! all cl files
        write(tempname,"(A12,I3.3,A7)") "AirfoilType_", airfoilIDs(i), "_CL.inp"
        fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
        ioUnit = 55
        open(unit=ioUnit,file=trim(fname), form='FORMATTED')
        read(ioUnit, *) num_entries_cl(i)
        close(ioUnit)

        ! all cd files
        write(tempname,"(A12,I3.3,A7)") "AirfoilType_", airfoilIDs(i), "_CD.inp"
        fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
        open(unit=ioUnit,file=trim(fname), form='FORMATTED')
        read(ioUnit, *) num_entries_cd(i)
        close(ioUnit)
    enddo
    ! allocate memory sufficient for the max
    size_table = maxval(num_entries_cl)
    allocate(this%clTable(size_table, 2, size(airfoilIDs)))
    size_table = maxval(num_entries_cd)
    allocate(this%cdTable(size_table, 2, size(airfoilIDs)))
    deallocate(num_entries_cl, num_entries_cd)

    ! fill all entries and pad where number of entries is smaller than max
    do i = 1, size(airfoilIDs)
        ! all cl files
        write(tempname,"(A12,I3.3,A7)") "AirfoilType_", airfoilIDs(i), "_CL.inp"
        fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
        open(unit=ioUnit,file=trim(fname), form='FORMATTED')
        read(ioUnit, *) size_table
        do j = 1, size_table
            read(ioUnit,*) this%clTable(j,:,i)
        enddo
        close(ioUnit)
        do j = size_table+1, size(this%clTable,1)
            this%clTable(j,:,i) = this%clTable(size_table,:,i)
        enddo
        this%clTable(:,1,i) = this%clTable(:,1,i)*deg_to_radians

        ! all cd files
        write(tempname,"(A12,I3.3,A7)") "AirfoilType_", airfoilIDs(i), "_CD.inp"
        fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
        open(unit=ioUnit,file=trim(fname), form='FORMATTED')
        read(ioUnit, *) size_table
        do j = 1, size_table
            read(ioUnit,*) this%cdTable(j,:,i)
        enddo
        close(ioUnit)
        do j = size_table+1, size(this%cdTable,1)
            this%cdTable(j,:,i) = this%cdTable(size_table,:,i)
        enddo
        this%cdTable(:,1,i) = this%cdTable(:,1,i)*deg_to_radians
    enddo

    ! store tableID for each element
    do i = 1, this%num_radial_elems
      do j = 1, size(airfoilIDs)
        if(airfoilIDs(j) == this%airfoilID_elem(i)) then
          this%tableID_elem(i) = j
        endif
      enddo
    enddo

end subroutine 


subroutine destroy(this)
    class(actuatordisk_Rot), intent(inout) :: this

    if (Allocated(this%rbuff))  deallocate(this%rbuff)
    if (Allocated(this%tag_face))  deallocate(this%tag_face)
    
    nullify(this%xG, this%yG, this%zG)
end subroutine 

subroutine getMeanU(this, u, v, w) 
    class(actuatordisk_Rot), intent(inout) :: this
    real(rkind), dimension(this%nxLoc,this%nyLoc,this%nzLoc), intent(in) :: u, v, w
    real(rkind) :: tmpSum, umn, vmn, wmn
    
    if (this%AM_I_ACTIVE) then
        ! Get u face
        this%rbuff = u(this%xface_idx,:,:)
        this%rbuff = this%rbuff*this%tag_face
        if (this%AM_I_Split) then
            tmpSum = p_sum(sum(this%rbuff),this%myComm) 
            umn = tmpSum/real(this%totPointsOnFace,rkind)
        else
            umn = sum(this%rbuff)/real(this%totPointsOnFace,rkind)
        end if

        ! Get v face
        this%rbuff = v(this%xface_idx,:,:)
        this%rbuff = this%rbuff*this%tag_face
        if (this%AM_I_Split) then
            tmpSum = p_sum(sum(this%rbuff),this%myComm) 
            vmn = (tmpSum)/real(this%totPointsOnFace,rkind)
        else
            vmn = sum(this%rbuff)/real(this%totPointsOnFace,rkind)
        end if
        
        ! Get w face
        this%rbuff = w(this%xface_idx,:,:)
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

        ! Get u disk
        this%rbuff = u(this%xdisk_idx,:,:); 
        this%rbuff = this%rbuff*this%tag_face
        if (this%AM_I_Split) then
            tmpSum = p_sum(sum(this%rbuff),this%myComm) 
            umn = tmpSum/real(this%totPointsOnFace,rkind)
        else
            umn = sum(this%rbuff)/real(this%totPointsOnFace,rkind)
        end if

        ! Get v disk
        this%rbuff = v(this%xdisk_idx,:,:)
        this%rbuff = this%rbuff*this%tag_face
        if (this%AM_I_Split) then
            tmpSum = p_sum(sum(this%rbuff),this%myComm) 
            vmn = (tmpSum)/real(this%totPointsOnFace,rkind)
        else
            vmn = sum(this%rbuff)/real(this%totPointsOnFace,rkind)
        end if
        
        ! Get w disk
        this%rbuff = w(this%xdisk_idx,:,:)
        this%rbuff = this%rbuff*this%tag_face
        if (this%AM_I_Split) then
            tmpSum = p_sum(sum(this%rbuff),this%myComm) 
            wmn = (tmpSum)/real(this%totPointsOnFace,rkind)
        else
            wmn = sum(this%rbuff)/real(this%totPointsOnFace,rkind)
        end if

        this%udisk = this%alpha_tau*umn + (1.d0 - this%alpha_tau)*this%udisk
        this%vdisk = this%alpha_tau*vmn + (1.d0 - this%alpha_tau)*this%vdisk
        this%wdisk = this%alpha_tau*wmn + (1.d0 - this%alpha_tau)*this%wdisk
    end if 
    this%alpha_tau = alpha_smooth

end subroutine

subroutine get_Omega(this)
    class(actuatordisk_Rot), intent(inout) :: this

    select case(this%RPMController)
    case(0)
      ! fixed TSR
      this%Omega = two*this%TSR*this%uface/this%diam  !; print '(2(i5,1x),a,e19.12)', nrank, this%ActuatorDisk_RotID, " Omega = ",this%Omega 
    case(1)
      ! fixed Omega read in at init, so do nothing
    case(2)
      ! adjust Omega to obtain a required CTp

      !! Step 1 :: Find current CTp
      !! !-----Already found and stored in this%CTp_eff

      !! Step 2 :: Drive Omega to the required CTp
      !! ! ---- Need to understand relation between CTp and TSR or Omega; If
      !! !they are proportional, then a PI controller can be implemented

      call GracefulExit("Only fixed TSR implemented for now. Set RPMController to 0", 111)
    case(3)
      ! variable based on some control algorithm to be implemented
      call GracefulExit("Only fixed TSR implemented for now. Set RPMController to 0", 111)
    end select

end subroutine

subroutine get_cl_cd(this,radind,aoa,cl,cd)
    class(actuatordisk_Rot), intent(in), target :: this
    integer,     intent(in)  :: radind
    real(rkind), intent(in)  :: aoa
    real(rkind), intent(out) :: cl, cd
    real(rkind) :: aoaloc(1), clloc(1), cdloc(1)
    real(rkind), dimension(:,:), pointer :: tablePtr

    aoaloc(1) = aoa

    tablePtr => this%clTable(1:size(this%clTable,1),1:2,this%tableID_elem(radind))
    call interp(aoaloc, clloc, tablePtr)

    tablePtr => this%cdTable(1:size(this%cdTable,1),1:2,this%tableID_elem(radind))
    call interp(aoaloc, cdloc, tablePtr)

    cl = clloc(1); cd = cdloc(1)

end subroutine


subroutine get_induction_factors(this,radind,dthrust,dtangen)
    class(actuatordisk_Rot), intent(inout) :: this
    integer,     intent(in)  :: radind
    real(rkind), intent(out) :: dthrust, dtangen
    real(rkind) :: radius, twist, fac, indfac, aprime, error, stopping_tol
    real(rkind) :: indfac_old, aprime_old, unrm, utan, urel, cosphi, sinphi
    real(rkind) :: phi, aoa, dfn, dft, drad, cl, cd
    integer     :: num_iter, max_iters

    dthrust = zero
    dtangen = zero

    radius = this%radius_elem(radind)
    twist = this%twist_elem(radind)
    fac = this%solidity_elem(radind)!*radius
    if(radind < this%num_radial_elems) then
      drad = this%radius_elem(radind+1)-this%radius_elem(radind)
    else
      drad = this%radius_elem(radind)-this%radius_elem(radind-1)
    endif

    if(radind <= this%NacelleRadInd) then
        urel   = this%uface
        phi    = half*pi
        sinphi = one
        cosphi = zero
        aprime = -one   ! dummy
        indfac = 0.1d0  ! dummy
        aoa    = 0.01d0 ! dummy
        call this%get_cl_cd(radind, aoa, cl, cd)
    else

        indfac = third; aprime = third
        dfn = 0.0d0; dft = 0.0d0; aoa = 0.0d0; cl =0.0d0; cd =0.0d0; phi = 0.0d0
        unrm = 0.0d0; utan = 0.0d0; urel = 0.0d0; cosphi = 0.0d0; sinphi = 0.0d0; 
        !if(nrank==4) then
        error = one; stopping_tol = 1.0d-4; num_iter = 0; max_iters = 100
        do while((error > stopping_tol) .and. (num_iter < max_iters))
            !if(nrank==0) then! .and. radind==3) then
            !  write(*,'(a,2(i5,1x),14(e19.12,1x))') '---deb---', radind, num_iter, unrm, utan, urel, &
            !               cosphi, sinphi, phi, aoa, cl, cd, dfn, dft, indfac, aprime, error
            !endif
            indfac_old = indfac; aprime_old = aprime

            unrm = this%uface*(one-indfac)
            utan = this%Omega*radius*(one+aprime)
            urel = sqrt(utan*utan + unrm*unrm)
            cosphi = utan/urel
            sinphi = unrm/urel
            phi = atan2(unrm, utan)
            aoa = phi - this%blade_pitch - twist
            call this%get_cl_cd(radind, aoa, cl, cd)

            !!-------first trial which did not work----------
            !dfn = fac*urel*urel*(cl*cosphi+cd*sinphi)
            !dft = fac*urel*urel*(cl*sinphi-cd*cosphi)
            !indfac = half*(one-sqrt(1-dfn/this%uface/this%uface))
            !aprime = dft/(four*this%uface*this%Omega*radius*(one-indfac))
            !!-------first trial which did not work----------

            !!-------second trial-----------------------------
            indfac = fac*(cl*cosphi+cd*sinphi)/(four*sinphi*sinphi + fac*(cl*cosphi+cd*sinphi))
            aprime = fac*(cl*sinphi-cd*cosphi)/(four*sinphi*cosphi - fac*(cl*sinphi-cd*cosphi)) 
            !!-------second trial-----------------------------

            ! realizability constraint
            if(indfac>half) then
              ! stop iterating and suppress the "error>stopping_tol" error message
              print '(a,i5,1x,5(e19.12,1x))', 'induction factor is beyond 0.5. Check details: ', radind, indfac, aprime, cl, cd, phi
              indfac_old = indfac; aprime_old = aprime 
              num_iter = max_iters + 2
              !stop
            endif

            ! stopping condition
            error = max(abs(indfac-indfac_old), abs(aprime-aprime_old))
            num_iter = num_iter + 1
        enddo

        if(error>stopping_tol) then
          write(*,'(a,i5,1x,3(e19.12,1x),2(i5,1x),2(e19.12,1x))') 'stopping: ', num_iter, error, indfac, &
                                                           aprime, radind, nrank, this%uface, this%Omega
        endif

    endif

    !!-------first trial which did not work----------
    !dthrust = dfn*pi*radius*drad
    !dtangen = dft*pi*radius*drad
    !!-------first trial which did not work----------

    !!-------second trial-----------------------------
    !fac = four*(one-indfac)*this%uface*pi*radius*drad
    !dthrust = fac*indfac*this%uface
    !dtangen = fac*aprime*this%Omega*radius

    fac = this%solidity_elem(radind)*pi*radius*drad*urel*urel
    dthrust = fac*(cl*cosphi+cd*sinphi)
    dtangen = fac*(cl*sinphi-cd*cosphi)
    !!-------second trial-----------------------------

    this%diagnosticarr(1, radind) = this%uface
    this%diagnosticarr(2, radind) = this%Omega
    this%diagnosticarr(3, radind) = indfac
    this%diagnosticarr(4, radind) = aprime
    this%diagnosticarr(5, radind) = phi
    this%diagnosticarr(6, radind) = aoa
    this%diagnosticarr(7, radind) = cl
    this%diagnosticarr(8, radind) = cd
    this%diagnosticarr(9, radind) = dthrust
    this%diagnosticarr(10,radind) = dtangen

end subroutine

subroutine get_RHS(this, u, v, w, rhsxvals, rhsyvals, rhszvals, inst_val)
    class(actuatordisk_Rot), intent(inout) :: this
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsxvals, rhsyvals, rhszvals
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(in)    :: u, v, w
    real(rkind), dimension(8),                                  intent(out)   :: inst_val
    real(rkind) :: dthrust, dtangen, total_thrust, total_torque
    integer     :: i, ist, ien, jst, jen, kst, ken, ioUnit, ierr
    character(len=clen) :: tempname, fname

    call this%getMeanU(u,v,w)
    call this%get_Omega()

    if(this%Am_I_Active) then
        total_thrust = zero
        total_torque = zero

        do i = 1, this%num_radial_elems
          call this%get_induction_factors(i,dthrust,dtangen)

          ist = (i-1)*this%num_azimut_elems + 1
          ien = i*this%num_azimut_elems
          this%force_x(ist:ien) = -dthrust*this%pfactor*this%normfactor
          this%force_y(ist:ien) = -dtangen*this%pfactor*this%normfactor*this%cos_azimut(1:this%num_azimut_elems)
          this%force_z(ist:ien) = -dtangen*this%pfactor*this%normfactor*this%sin_azimut(1:this%num_azimut_elems)

          total_thrust = total_thrust - dthrust
          total_torque = total_torque - dtangen*this%radius_elem(i)
        enddo

        if(this%isFan) then
            !fan_mult_fac = -one
            this%force_x = -this%force_x
            !this%force_y = -this%force_y
            !this%force_z = -this%force_z
            total_thrust = -total_thrust
            !total_torque = zero
        endif


        ist = this%startEnds_global(1);   ien = this%startEnds_global(4)
        jst = this%startEnds_global(2);   jen = this%startEnds_global(5)
        kst = this%startEnds_global(3);   ken = this%startEnds_global(6)

        this%force_local = zero
        do i = 1, size(this%xs)
          this%force_local = this%force_local + this%force_x(i)*this%smearing_base(:,:,:,i)
        enddo
        rhsxvals(ist:ien,jst:jen,kst:ken) = rhsxvals(ist:ien,jst:jen,kst:ken) + this%force_local

        this%force_local = zero
        do i = 1, size(this%xs)
          this%force_local = this%force_local + this%force_y(i)*this%smearing_base(:,:,:,i)
        enddo
        rhsyvals(ist:ien,jst:jen,kst:ken) = rhsyvals(ist:ien,jst:jen,kst:ken) + this%force_local

        this%force_local = zero
        do i = 1, size(this%xs)
          this%force_local = this%force_local + this%force_z(i)*this%smearing_base(:,:,:,i)
        enddo
        rhszvals(ist:ien,jst:jen,kst:ken) = rhszvals(ist:ien,jst:jen,kst:ken) + this%force_local

        !if (present(inst_val)) then
          if((this%Am_I_Split .and. this%myComm_nrank==0) .or. (.not. this%Am_I_Split)) then
            inst_val(1) = total_thrust
            inst_val(2) = total_thrust*sqrt(this%udisk**2+this%vdisk**2+this%wdisk**2)
            inst_val(3) = total_torque
            inst_val(4) = this%Omega
            inst_val(5) = total_torque*this%Omega
            inst_val(6) = this%udisk
            inst_val(7) = this%vdisk
            inst_val(8) = this%wdisk

            this%diagnostic_counter = this%diagnostic_counter + 1
            if(mod(this%diagnostic_counter, this%diagnostics_write_freq)==0) then
                ioUnit = 100; write(tempname,"(a6,i4.4,a12,i5.5,a4)") 'ADRot_', this%ActuatorDisk_RotID, '_diagnostic_', this%diagnostic_counter, '.dat'
                fname = trim(tempname)
                open(ioUnit,file=fname,form='formatted',action='write',status='unknown')
                do i = 1, this%num_radial_elems
                  write(ioUnit,'(10(e19.12,1x))') this%diagnosticarr(:,i)
                enddo
                close(ioUnit)
            end if
          end if
        !end if 

        ! calculate the effective CTp
        this%CTp_eff = -total_thrust/(half * (this%udisk**2 + this%vdisk**2 + this%wdisk**2) * this%DiskArea)

    endif

end subroutine

subroutine smear_this_source(this,rhsfield, xC, yC, zC, valSource)
    class(actuatordisk_Rot), intent(inout) :: this
    real(rkind), dimension(:,:,:), intent(inout) :: rhsfield
    real(rkind), intent(in) :: xC, yC, zC, valSource
    integer :: ist, ien, jst, jen, kst, ken, xlen, jj, kk, jloc, kloc

    ist = this%startEnds_global(1);   ien = this%startEnds_global(4)
    jst = this%startEnds_global(2);   jen = this%startEnds_global(5)
    kst = this%startEnds_global(3);   ken = this%startEnds_global(6)
    xlen = this%startEnds_global(7)

    do kk = kst, ken
        do jj = jst, jen
            kloc = kk-kst+1; jloc = jj-jst+1
            this%dline(1:xlen) = (this%xG(ist:ien,jj,kk) - xC)**2 + &
                                 (this%yG(ist:ien,jj,kk) - yC)**2 + &
                                 (this%zG(ist:ien,jj,kk) - zC)**2
            this%dline(1:xlen) = -this%dline(1:xlen)*this%oneBydelSq
            !write(*,'(12(i5,1x))') kk, jj, kst, ken, jst, jen, ist, ien, xlen, jloc, kloc, nrank
            rhsfield(1:xlen,jloc,kloc) = rhsfield(1:xlen,jloc,kloc) + valSource*exp(this%dline(1:xlen))
        end do
    end do  

end subroutine 

subroutine sample_on_circle_rtheta(this,R,xcen,ycen,zcen,np)
    use gridtools, only: linspace
    class(actuatordisk_Rot), intent(inout) :: this
    real(rkind), intent(in) :: R, xcen, ycen, zcen
    integer, intent(in) :: np
    !real(rkind), dimension(:), allocatable, intent(out) :: xloc, yloc
    real(rkind) :: drad, dthet, azimut
    integer :: i, j, iidx

    allocate(this%radius_elem(np),this%cos_azimut(np),this%sin_azimut(np))
    allocate(this%xs(np**2),this%ys(np**2),this%zs(np**2))

    this%num_radial_elems = np
    drad  = R/real(this%num_radial_elems,rkind)
    do i = 1, this%num_radial_elems
       this%radius_elem(i) = (real(i,rkind)-0.5_rkind)*drad
    enddo

    this%num_azimut_elems = np
    dthet = twopi/real(this%num_azimut_elems,rkind)
    do i = 1, this%num_azimut_elems
       azimut = (real(i,rkind)-0.5_rkind)*dthet
       this%cos_azimut(i) = cos(azimut)
       this%sin_azimut(i) = sin(azimut)
    enddo

    iidx = 0
    do j=1,this%num_radial_elems
     do i=1,this%num_azimut_elems
       iidx = iidx+1
       this%ys(iidx) = this%radius_elem(j)*this%cos_azimut(i)
       this%zs(iidx) = this%radius_elem(j)*this%sin_azimut(i)
     enddo
    enddo
    this%ys = this%ys + ycen; this%zs = this%zs + zcen 
    this%xs = xcen

end subroutine

subroutine interp(x, y, xyTable)
    !class(actuatorLine), intent(inout) :: this
    real(rkind), intent(in),  dimension(:)   :: x
    real(rkind), intent(out), dimension(:)   :: y
    real(rkind), intent(in),  dimension(:,:) :: xyTable

    integer :: i, locator(1), imin, i1, i2

       !write(*,'(a,3(i5,1x),4(f9.3,1x))') '++', nrank, size(xyTable,1), size(xyTable,2), maxval(xyTable(:,1)), minval(xyTable(:,1)), maxval(xyTable(:,2)), minval(xyTable(:,2))
       !  write(100+nrank,*) '---------'
       !do i = 1, size(xyTable,1)
       !  write(100+nrank,*) '++** ', xyTable(i,:)
       !enddo
       !  write(100+nrank,*) '---------'
    do i = 1, size(x)
       locator = minloc(abs(xyTable(:,1) - x(i))); imin = locator(1)
       if(xyTable(imin,1) < x(i)) then
         i1 = imin; i2 = imin+1
       else
         i1 = imin-1; i2 = imin
       endif
       if(i2 > size(xyTable,1)) then
         i1 = i1-1; i2 = i2-1
       endif
       if(i1 < 1) then
         i1 = i1+1; i2 = i2+1
       endif
       !write(*,'(a,4(i5,1x),5(f6.2,1x))') '--', nrank, i, i1, i2, xyTable(i2,2), xyTable(i1,2), xyTable(i2,1), xyTable(i1,1), x(i)
       y(i) = xyTable(i1,2) + (xyTable(i2,2)-xyTable(i1,2))/(xyTable(i2,1)-xyTable(i1,1))*(x(i) - xyTable(i1,1))
    enddo

end subroutine

subroutine interpint(x, y, xyTable)
    !class(actuatorLine), intent(inout) :: this
    real(rkind), intent(in),  dimension(:)   :: x
    integer,     intent(out), dimension(:)   :: y
    real(rkind), intent(in),  dimension(:,:) :: xyTable

    integer :: i, j

    do i = 1, size(x)
      do j = 1, size(xyTable,1)-1
        if(x(i)>=xyTable(j,1) .and. x(i)<=xyTable(j+1,1)) exit
      enddo
      y(i) = int(xyTable(j,2))
    enddo

end subroutine

subroutine find_unique_entries(array, unique_entries)
    integer, intent(in),  dimension(:) :: array
    integer, intent(out), allocatable, dimension(:) :: unique_entries

    integer, allocatable, dimension(:) :: temparr
    integer :: num_unique, j, k
    logical :: isEntryUnique

    allocate(temparr(size(array)))
    temparr = -1

    num_unique = 1
    temparr(num_unique) = array(1)

    do k = 2, size(array)
      isEntryUnique = .true.
      do j = 1, num_unique
        if(array(k) == temparr(j)) then
          isEntryUnique = .false.
          exit
        endif
      enddo
      if(isEntryUnique) then
        num_unique = num_unique+1
        temparr(num_unique) = array(k)
      endif
    enddo

    allocate(unique_entries(num_unique))
    unique_entries(1:num_unique) = temparr(1:num_unique)

    deallocate(temparr)

end subroutine


end module 
