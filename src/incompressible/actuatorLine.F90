module actuatorLineMod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa
    use decomp_2d
    use StaggOpsMod, only: staggOps  
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral  
    use mpi 
    use reductions, only: p_maxval, p_sum
    use timer, only: tic, toc

    implicit none

    private
    public :: actuatorLine

    ! default initializations
    integer :: num_turbines = 1, num_blades = 3, num_blade_points
    real(rkind) :: initial_yaw = zero, initial_azimuth = zero, initial_rpm = 8.0_rkind
    real(rkind) :: turb_xloc = pi, turb_yloc = pi, turb_zloc = 0.4_rkind
    real(rkind) :: tip_radius = 0.1_rkind, hub_radius = 0.05_rkind, epsfactor = 1.5d0
    real(rkind) :: nacelle_width = 0.1_rkind

    integer :: ioUnit

    !real(rkind) :: epsilon_sq, eps_pi_fac
    real(rkind), parameter :: degrees_to_radians = pi/180.0_rkind
    real(rkind), parameter :: rpm_to_radianspersec = two*pi/60.0_rkind

    type :: actuatorLine
        real(rkind) :: xLoc, yLoc, zLoc, tip_radius, hub_radius, hub_height
        integer     :: num_blades
        real(rkind) :: yaw_angle, blade_azimuth, rotspeed, nacelle_width
        logical     :: clockwise_rotation
        integer     :: num_blade_points  ! number of actuator points on each blade

        integer     :: nxLoc, nyLoc, nzLoc
        integer     :: ist, iend, jst, jend, kst, kend, xlen, ylen, zlen
        real(rkind) :: delta, OneByDelSq, turb_thrust, turb_torque, normfactor, invdxdy, invdz, distr_thrust
        real(rkind) :: uturbavg
        real(rkind) :: xRightPad, yRightPad, zLeftPad, zRightPad
        real(rkind), dimension(3) :: turbLoc, rotor_shaft, rotor_center
        real(rkind), dimension(:,:,:),   allocatable :: blade_points  ! number of actuator points on each blade
        real(rkind), dimension(:,:,:),   allocatable :: blade_forces  ! forces at actuator points
        real(rkind), dimension(:,:,:),   allocatable :: blade_forcesloc  ! needed if actuator points are split across processors
        real(rkind), dimension(:,:),     allocatable :: radial_dist, chord, twistAng  ! data at actuator points
        integer,     dimension(:,:),     allocatable :: airfoilID, airfoilIDIndex     ! data at actuator points
        real(rkind), dimension(:,:,:),   allocatable :: clTable, cdTable          ! CL and CD tables for each airfoil used
        integer,     dimension(:),       allocatable :: clTableSize, cdTableSize  ! Sizes of the CL and CD tables. All airfoils need not have same-sized tables
        real(rkind), dimension(:,:,:),   allocatable :: dsq, xSmall, ySmall, zSmall
        real(rkind), dimension(:,:,:,:), allocatable :: source
        real(rkind), dimension(:,:,:),   allocatable :: blade_statsloc, blade_stats  ! needed if actuator points are split across processors

        ! MPI communicator stuff
        logical :: Am_I_Active, Am_I_Split
        integer :: tag_proc, myComm, myComm_nproc, myComm_nrank

    contains

        procedure :: init
        procedure :: destroy
        procedure :: get_RHS 
        procedure, private :: get_extents
        procedure, private :: distribute_forces
        procedure, private :: get_blade_forces 
        procedure, private :: rotate_one_blade 
        procedure, private :: yaw_turbine
        procedure, private :: interp_velocity
        procedure, private :: interp_clcd
        procedure, private :: get_rotation_speed
        procedure, private :: update_turbine

    end type

contains

subroutine init(this, inputDir, ActuatorLineID, xG, yG, zG, xyzPads)
    use mpi
    class(actuatorLine), intent(inout) :: this
    character(len=*), intent(in)                           :: inputDir
    integer,          intent(in)                           :: ActuatorLineID
    real(rkind),      intent(in), dimension(:,:,:), target :: xG, yG, zG
    real(rkind),      intent(in), dimension(6)             :: xyzPads

    character(len=clen) :: tempname, fname
    integer :: j, k, ierr,  max_size_cltable, max_size_cdtable
    integer :: num_table_entries, turbineTypeID, num_airfoil_types
    integer :: num_blades, num_blade_points, bladeTypeID, size_table
    real(rkind) :: element_length, projection_radius, cloud_radius
    real(rkind) :: xmin, xmax, ymin, ymax, zmin, zmax 
    real(rkind) :: dx, dy, dz, xLoc, yLoc, zLoc, rotor_diam, hub_diam, hub_height, initial_yaw, reference_length, reference_velocity
    real(rkind) :: initial_azimuth, initial_rpm, nacelle_width, epsFact
    logical     :: isEntryUnique, clockwise_rotation

    real(rkind), allocatable, dimension(:,:) :: twistAngTable, chordTable, airfoilIDTable, airfoilIDIndexTable
    integer,     allocatable, dimension(:)   :: airfoilIDs

    namelist /ACTUATOR_LINE/ xLoc, yLoc, zLoc, turbineTypeID
    namelist /TURBINE_TYPE/ rotor_diam, hub_diam, hub_height, num_blades, initial_yaw, initial_azimuth, &
                        initial_rpm, nacelle_width, clockwise_rotation, num_blade_points, epsFact, bladeTypeID, &
                        reference_length, reference_velocity

    ! Read input file for this turbine    
    write(tempname,"(A13,I3.3,A10)") "ActuatorLine_", ActuatorLineID, "_input.inp"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)

    ioUnit = 55
    open(unit=ioUnit, file=trim(fname), form='FORMATTED')
    read(unit=ioUnit, NML=ACTUATOR_LINE)
    close(ioUnit)
    !call mpi_barrier(mpi_comm_world, ierr); call message(2,"Read ActuatorLine_input.inp"); call mpi_barrier(mpi_comm_world, ierr)

    ! Read information about this turbine type
    write(tempname,"(A12,I3.3,A10)") "TurbineType_", turbineTypeID, "_input.inp"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)

    ioUnit = 55
    open(unit=ioUnit, file=trim(fname), form='FORMATTED')
    read(unit=ioUnit, NML=TURBINE_TYPE)
    close(ioUnit)
    !call mpi_barrier(mpi_comm_world, ierr); call message(2,"Read TurbineType_input.inp"); call mpi_barrier(mpi_comm_world, ierr)

    ! Read information about blade type used in this turbine - identical blades per turbine for now
    ! This file does not use namelists, so order in which data is written in this file appears is important
    write(tempname,"(A10,I3.3,A10)") "BladeType_", bladeTypeID, "_input.inp"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)

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
      twistAngTable(:,2) = twistAngTable(:,2) * degrees_to_radians

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
    !call mpi_barrier(mpi_comm_world, ierr); call message(2,"Read BladeType_input.inp"); call mpi_barrier(mpi_comm_world, ierr)

    ! read in airfoil data
    ! first identify the number of airfoils used in this blade and the size of
    ! corresponding CL and CD table
    num_table_entries = size(airfoilIDTable,1)
    max_size_cltable = 0; max_size_cdtable = 0

    allocate(airfoilIDs(num_table_entries))    ! at most, all entries can be unique

    ! at least one airfoiltype is always present
    num_airfoil_types = 1
    airfoilIDs(num_airfoil_types) = airfoilIDTable(1,2)

    ! determine size of CL table
    write(tempname,"(A12,I3.3,A7)") "AirfoilType_", airfoilIDs(num_airfoil_types), "_CL.inp"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    open(unit=ioUnit,file=trim(fname), form='FORMATTED')
    read(ioUnit, *) size_table
    close(ioUnit)
    max_size_cltable = max(size_table, max_size_cltable)

    ! determine size of CD table
    write(tempname,"(A12,I3.3,A7)") "AirfoilType_", airfoilIDs(num_airfoil_types), "_CD.inp"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    open(unit=ioUnit,file=trim(fname), form='FORMATTED')
    read(ioUnit, *) size_table
    close(ioUnit)
    max_size_cdtable = max(size_table, max_size_cdtable)

    do k = 2, num_table_entries
      isEntryUnique = .true.
      do j = 1, k-1
        if(airfoilIDTable(k,2) == airfoilIDTable(j,2)) then
          isEntryUnique = .false.
          exit
        endif
      enddo
      if(isEntryUnique) then
        num_airfoil_types = num_airfoil_types + 1
        airfoilIDs(num_airfoil_types) = airfoilIDTable(k,2)

        ! determine size of CL table
        write(tempname,"(A12,I3.3,A7)") "AirfoilType_", airfoilIDs(num_airfoil_types), "_CL.inp"
        fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
        open(unit=ioUnit,file=trim(fname), form='FORMATTED')
        read(ioUnit, *) size_table
        close(ioUnit)
        max_size_cltable = max(size_table, max_size_cltable)
    
        ! determine size of CD table
        write(tempname,"(A12,I3.3,A7)") "AirfoilType_", airfoilIDs(num_airfoil_types), "_CD.inp"
        fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
        open(unit=ioUnit,file=trim(fname), form='FORMATTED')
        read(ioUnit, *) size_table
        close(ioUnit)
        max_size_cdtable = max(size_table, max_size_cdtable)
      endif
    enddo

    ! set airfoilIDIndexTable
    allocate(airfoilIDIndexTable(size(airfoilIDTable,1), 2))
    airfoilIDIndexTable = -1.0
    do k = 1, size(airfoilIDTable,1)
      do j = 1, num_airfoil_types
        if(nint(airfoilIDTable(k,2)) == airfoilIDs(j)) then
            airfoilIDIndexTable(k,1) = airfoilIDTable(k,1)
            airfoilIDIndexTable(k,2) = j
            exit
        endif
      enddo
      if(nint(airfoilIDIndexTable(k,1))==-1) then 
        write(*,*) 'airfoilIDTable(k,2): ', airfoilIDTable(k,2), ' not found in airfoilIDs:', airfoilIDs(:)
      endif
    enddo

    ! allocate memory and read in table entries for each airfoil type
    allocate(this%clTable(max_size_cltable, max_size_cltable, num_airfoil_types), this%clTableSize(num_airfoil_types))
    allocate(this%cdTable(max_size_cdtable, max_size_cdtable, num_airfoil_types), this%cdTableSize(num_airfoil_types))

    do k = 1, num_airfoil_types
        ! read CL table
        write(tempname,"(A12,I3.3,A7)") "AirfoilType_", airfoilIDs(num_airfoil_types), "_CL.inp"
        fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
        open(unit=ioUnit,file=trim(fname), form='FORMATTED')
        read(ioUnit, *) this%clTableSize(k)
        do j = 1, this%clTableSize(k)
          read(ioUnit, *) this%clTable(j,1,k), this%clTable(j,2,k)
          this%clTable(j,1,k) = this%clTable(j,1,k) * degrees_to_radians
        enddo
        close(ioUnit)
        !call mpi_barrier(mpi_comm_world, ierr); call message(2,"Read AirfoilType_CL_input.inp"); call mpi_barrier(mpi_comm_world, ierr)
       
        !if(nrank==0) then
        !  write(*,*) '-----ClTable:-----', k
        !  do j = 1, this%clTableSize(k)
        !    write(*, *) this%clTable(j,1,k), this%clTable(j,2,k)
        !  enddo
        !endif
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Done write CL table"); call mpi_barrier(mpi_comm_world, ierr)
        !stop 

        ! read CD table
        write(tempname,"(A12,I3.3,A7)") "AirfoilType_", airfoilIDs(num_airfoil_types), "_CD.inp"
        fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
        open(unit=ioUnit,file=trim(fname), form='FORMATTED')
        read(ioUnit, *) this%cdTableSize(k)
        do j = 1, this%cdTableSize(k)
          read(ioUnit, *) this%cdTable(j,1,k), this%cdTable(j,2,k)
          this%cdTable(j,1,k) = this%cdTable(j,1,k) * degrees_to_radians
        enddo
        close(ioUnit)
        !call mpi_barrier(mpi_comm_world, ierr); call message(2,"Read AirfoilType_CD_input.inp"); call mpi_barrier(mpi_comm_world, ierr)
    enddo
    deallocate(airfoilIDs)

    this%xLoc = xLoc; this%yLoc = yLoc; this%zLoc = zLoc
    this%tip_radius = half*rotor_diam/reference_length
    this%hub_radius = half*hub_diam  /reference_length
    this%hub_height = hub_height     /reference_length
    this%num_blades = num_blades;  
    this%yaw_angle  = initial_yaw * degrees_to_radians        ! initialize turbines with initial yaw zero (facing positive x direction)
    this%blade_azimuth = initial_azimuth * degrees_to_radians ! initialize first blade to zero azimuth; others are equally spaced in the 360 degree space
    this%rotspeed = initial_rpm * rpm_to_radianspersec * reference_length/reference_velocity
    this%nacelle_width = nacelle_width / reference_length
    this%clockwise_rotation = clockwise_rotation
    this%num_blade_points = num_blade_points

    dx=xG(2,1,1)-xG(1,1,1); dy=yG(1,2,1)-yG(1,1,1); dz=zG(1,1,2)-zG(1,1,1)
    this%nxLoc = size(xG,1); this%nyLoc = size(xG,2); this%nzLoc = size(xG,3)
    this%invdxdy = one/(dx*dy)
    this%invdz = one/dz

    ! factor for distributing blade forces to grid points
    this%delta = epsFact * (dx*dy*dz)**(1.d0/3.d0)
    this%OneByDelSq = 1.d0/(this%delta**2)
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Step 1"); call mpi_barrier(mpi_comm_world, ierr)

    allocate(this%blade_points(3, this%num_blade_points, this%num_blades))
    allocate(this%blade_forces(3, this%num_blade_points, this%num_blades))
    allocate(this%blade_forcesloc(3, this%num_blade_points, this%num_blades))
    allocate(this%blade_statsloc(28, this%num_blade_points, this%num_blades))
    allocate(this%blade_stats(28, this%num_blade_points, this%num_blades))

    allocate(this%radial_dist(this%num_blade_points, this%num_blades))
    allocate(this%twistAng(   this%num_blade_points, this%num_blades))
    allocate(this%chord(      this%num_blade_points, this%num_blades))
    allocate(this%airfoilID(  this%num_blade_points, this%num_blades))
    allocate(this%airfoilIDIndex(this%num_blade_points, this%num_blades))

    !do i = 1, this%nTurbines

    ! turbLoc is offset from (xLoc, yLoc, zLoc) by hub_height in z
    this%turbLoc(1) = this%xLoc
    this%turbLoc(2) = this%yLoc
    this%turbLoc(3) = this%zLoc + this%hub_height

    ! rotor center is offset from turbLoc by Nacelle width in x-y plane (only in x direction initially)
    this%rotor_center(1) = this%turbLoc(1) - this%nacelle_width
    this%rotor_center(2) = this%turbLoc(2)
    this%rotor_center(3) = this%turbLoc(3)

    ! rotor shaft points in the direction of the wind
    if(this%nacelle_width > 1.0d-10) then
        this%rotor_shaft = this%turbLoc - this%rotor_center
    else
        this%rotor_shaft(1) = one;  this%rotor_shaft(2:3) = zero
    endif
    this%rotor_shaft = this%rotor_shaft/sqrt(sum(this%rotor_shaft**2))

    ! length of each actuator line element
    element_length = (this%tip_radius - this%hub_radius) / this%num_blade_points
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Step 2"); call mpi_barrier(mpi_comm_world, ierr)

    do j = 1, this%num_blades
      do k = 1, this%num_blade_points
        ! set blade points beginning from the hub and going radially outward
        this%radial_dist(k,j) = this%hub_radius + element_length * (real(k-1, rkind) + half)

        ! initialize each blade with initial azimuth zero (vertically upwards)
        this%blade_points(:,k,j) = this%rotor_center
        this%blade_points(3,k,j) = this%rotor_center(3) + this%radial_dist(k,j)
      enddo
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Step 2.1"); call mpi_barrier(mpi_comm_world, ierr)

      ! rotate each blade to its correct azimuth
      call this%rotate_one_blade(j, this%blade_azimuth + real(j-1,rkind)*two*pi/real(this%num_blades, rkind))
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Step 2.2"); call mpi_barrier(mpi_comm_world, ierr)

      ! interpolate airfoil properties - twist angle, chord length, airfoilID - based on radial distance
      call interp(this%radial_dist(:,j), this%twistAng(:,j),  twistAngTable )
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Step 2.3"); call mpi_barrier(mpi_comm_world, ierr)
      call interp(this%radial_dist(:,j), this%chord(:,j),     chordTable    )
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Step 2.4"); call mpi_barrier(mpi_comm_world, ierr)
      call interpint(this%radial_dist(:,j), this%airfoilID(:,j), airfoilIDTable)
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Step 2.5"); call mpi_barrier(mpi_comm_world, ierr)
      call interpint(this%radial_dist(:,j), this%airfoilIDIndex(:,j), airfoilIDIndexTable)
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Step 2.6"); call mpi_barrier(mpi_comm_world, ierr)
    enddo
    deallocate(twistAngTable, chordTable, airfoilIDTable, airfoilIDIndexTable)

    ! now rotate all blades so that yaw angle is correct
    call this%yaw_turbine(this%yaw_angle)
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Step 2.7"); call mpi_barrier(mpi_comm_world, ierr)

    ! compute ist, iend, jst, jend, kst, kend based on turbine location, blade radius, projection radius and processor extents
    projection_radius = 1.1d0 * this%delta * sqrt(log(1.0D3))    ! distance where influence of actuator points reduces to 0.001 of the max value
    cloud_radius = this%tip_radius + projection_radius
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Step 2.8"); call mpi_barrier(mpi_comm_world, ierr)
 
    xmin = this%rotor_center(1) - cloud_radius; xmax = this%rotor_center(1) + cloud_radius
    ymin = this%rotor_center(2) - cloud_radius; ymax = this%rotor_center(2) + cloud_radius
    zmin = this%rotor_center(3) - cloud_radius; zmax = this%rotor_center(3) + cloud_radius
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Step 2.9"); call mpi_barrier(mpi_comm_world, ierr)

    call this%get_extents(1, xmin, xmax, xG(:,1,1))
    call this%get_extents(2, ymin, ymax, yG(1,:,1))
    call this%get_extents(3, zmin, zmax, zG(1,1,:))
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Step 2.10"); call mpi_barrier(mpi_comm_world, ierr)
 
    this%xlen = this%iend - this%ist + 1
    this%ylen = this%jend - this%jst + 1
    this%zlen = this%kend - this%kst + 1
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Step 3"); call mpi_barrier(mpi_comm_world, ierr)

    if(this%xlen * this%ylen * this%zlen > 0) then
      this%tag_proc = 1
      this%Am_I_Active = .true.
      !this%color = ActuatorLineID
    else
      this%tag_proc = 0
      this%Am_I_Active = .false.
      !this%color = ActuatorLineID*1000
    endif
    if(p_sum(this%tag_proc) > 1) then
      this%Am_I_Split = .true.
    else
      this%Am_I_Split = .false.
    endif

    if (this%Am_I_Split) then
        call MPI_COMM_SPLIT(mpi_comm_world, this%tag_proc, nrank, this%myComm, ierr)
        call MPI_COMM_RANK( this%myComm, this%myComm_nrank, ierr )
        call MPI_COMM_SIZE( this%myComm, this%myComm_nproc, ierr )
    end if

    ! allocate buffers ranging from min index to max index in each direction
    if(this%Am_I_Active) then
      allocate( this%dsq(   1:this%xlen, 1:this%ylen, 1:this%zlen   ))
      allocate( this%xSmall(1:this%xlen, 1:this%ylen, 1:this%zlen   ))
      allocate( this%ySmall(1:this%xlen, 1:this%ylen, 1:this%zlen   ))
      allocate( this%zSmall(1:this%xlen, 1:this%ylen, 1:this%zlen   ))
      allocate( this%source(1:this%xlen, 1:this%ylen, 1:this%zlen, 3))

      this%xSmall = xG(this%ist:this%iend, this%jst:this%jend, this%kst:this%kend)
      this%ySmall = yG(this%ist:this%iend, this%jst:this%jend, this%kst:this%kend)
      this%zSmall = zG(this%ist:this%iend, this%jst:this%jend, this%kst:this%kend)

    endif
        !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Step 4"); call mpi_barrier(mpi_comm_world, ierr)

    if(this%Am_I_Active) then
        this%normfactor = one/(this%delta**3 * pi**1.5d0)
        !this%normfactor = one/(real(this%xlen*this%ylen*this%zlen, rkind) * this%delta**3 * pi**1.5d0)
    endif

    ! x, y and z locations on the left and right. Used in conjunction with halo
    ! values (if required) in interpolation of velocities
    this%xRightPad = xyzPads(2);   this%yRightPad = xyzPads(4)
    this%zLeftPad  = xyzPads(5);   this%zRightPad = xyzPads(6)

    write(*,*) 'Done init subroutine'
    write(*,*) nrank, this%Am_I_Active, this%Am_I_Split
    if(this%Am_I_Active) then
      if((this%Am_I_Split .and. this%myComm_nrank==0) .or. (.not. this%Am_I_Split)) then
          write(*,*) 'Writing blade_init.dat from proc. no.', this%myComm_nrank, nrank
          open(13,file='blade_init.dat',status='replace',action='write')
          write(13,'(a300)') 'VARIABLES="i","x","y",z","radialdist","twistAng","chord","airfoilIDIndex"'
          do j = 1, this%num_blades
            write(13,'(a,i4.4,a,i6)') 'ZONE T="Blade ', k, '", F=POINT, I=', this%num_blade_points
            do k = 1, this%num_blade_points
              write(13,'(i6,1x,6(e19.12,1x),i6)') j, this%blade_points(:,k,j), this%radial_dist(k,j), this%twistAng(k,j), this%chord(k,j), this%airfoilIDIndex(k,j)
            enddo
          enddo
          close(13)
      endif
    endif

end subroutine


subroutine destroy(this)
    class(actuatorLine), intent(inout) :: this
    if(this%Am_I_Active) deallocate(this%dsq, this%xSmall, this%ySmall, this%zSmall, this%source)
    deallocate(this%airfoilIDIndex, this%airfoilID, this%chord, this%twistAng, this%radial_dist, this%blade_stats, this%blade_statsloc, this%blade_forcesloc, this%blade_forces, this%blade_points)
    deallocate(this%cdTableSize, this%cdTable, this%clTableSize, this%clTable)
end subroutine

subroutine get_extents(this, idir, xmin, xmax, procmesh)
    class(actuatorLine), intent(inout) :: this
    integer, intent(in) :: idir
    real(rkind), intent(in) :: xmin, xmax
    real(rkind), dimension(:), intent(in) :: procmesh

    integer :: nloc, ind, ist, iend

    nloc = size(procmesh)

    ! determine x extents
    if(xmax < procmesh(1)) then
       ! turbine i influence is exclusively to the left of this processor
       ist = 0; iend = -1;
    else if(xmin > procmesh(nloc)) then
       ! turbine i influence is exclusively to the right of this processor
       ist = 0; iend = -1;
    else
       ! turbine i influences this processor

       ! determine ist
       if(xmin < procmesh(1)) then
          ist = 1
       else
          ist = -1
          do ind = 1, nloc-1
            if(xmin >= procmesh(ind) .and. xmin < procmesh(ind+1)) exit
          enddo
          ist = ind+1
          if(ind == nloc) then
             write(*,*) 'Something wrong. Check details.', xmin, procmesh(1), procmesh(nloc)
             call GracefulExit("Exiting from min extents setup", 423)
          endif
       endif

       ! determine iend
       if(xmax > procmesh(nloc)) then
          iend = nloc
       else
          iend = -1
          do ind = 1, nloc-1
            if(xmax >= procmesh(ind) .and. xmax < procmesh(ind+1)) exit
          enddo
          iend = ind
          if(ind == nloc) then
             write(*,*) 'Something wrong. Check details.', xmax, procmesh(1), procmesh(nloc)
             call GracefulExit("Exiting from max extents setup", 423)
          endif
       endif
    endif

    if(idir==1) then
       this%ist = ist; this%iend = iend
    elseif(idir==2) then
       this%jst = ist; this%jend = iend
    elseif(idir==3) then
       this%kst = ist; this%kend = iend
    endif

end subroutine

subroutine yaw_turbine(this, angle)
    class(actuatorLine), intent(inout) :: this
    real(rkind), intent(in) :: angle

    real(rkind) :: axis(3), cosa, sina, onemcosa, onemsina, rot_matrix(3,3)
    integer :: ptID, blID

    ! set unit vector in the direction of axis of rotation
    axis(:) = zero; axis(3) = one

    ! cos and sin of angle 
    cosa = cos(angle); sina = sin(angle)
    onemcosa = one - cosa; onemsina = one - sina

    ! set rotation matrix (wikipedia, also same in SOWFA)
    rot_matrix(1,1) = axis(1)**2*onemcosa+cosa;  rot_matrix(1,2) = axis(1)*axis(2)*onemcosa-axis(3)*sina;  rot_matrix(1,3) = axis(1)*axis(3)*onemcosa+axis(2)*sina
    rot_matrix(2,2) = axis(2)**2*onemcosa+cosa;  rot_matrix(2,3) = axis(2)*axis(3)*onemcosa-axis(1)*sina;  rot_matrix(2,1) = axis(2)*axis(1)*onemcosa+axis(3)*sina
    rot_matrix(3,3) = axis(3)**2*onemcosa+cosa;  rot_matrix(3,1) = axis(3)*axis(1)*onemcosa-axis(2)*sina;  rot_matrix(3,2) = axis(3)*axis(2)*onemcosa+axis(1)*sina

    ! first rotate rotor_center
    this%rotor_center = this%rotor_center - this%turbLoc
    this%rotor_center = matmul(rot_matrix, this%rotor_center)
    this%rotor_center = this%rotor_center + this%turbLoc

    ! update shaft direction
    if(this%nacelle_width > 1.0d-10) then
        this%rotor_shaft = this%turbLoc - this%rotor_center
    else
        this%rotor_shaft(1) = one;  this%rotor_shaft(2:3) = zero
    endif
    this%rotor_shaft = this%rotor_shaft/sqrt(sum(this%rotor_shaft**2))

    ! next rotate each blade
    do blID = 1, this%num_blades
      do ptID = 1, this%num_blade_points
         ! shift origin to turbLoc
         this%blade_points(:, ptID, blID) = this%blade_points(:, ptID, blID) - this%turbLoc

         ! rotate point in this frame of reference
         this%blade_points(:, ptID, blID) = matmul(rot_matrix, this%blade_points(:, ptID, blID))

         ! revert to origin
         this%blade_points(:, ptID, blID) = this%blade_points(:, ptID, blID) + this%turbLoc
      enddo
    enddo
end subroutine

subroutine rotate_one_blade(this, blID, angle)
    class(actuatorLine), intent(inout) :: this
    integer, intent(in) :: blID
    real(rkind), intent(in) :: angle

    real(rkind) :: axis(3), cosa, sina, onemcosa, onemsina, rot_matrix(3,3)
    integer :: ptID

    ! set unit vector in the direction of axis of rotation
    axis(:) = this%rotor_shaft(:)
    axis = axis/sqrt(sum(axis**2))

    ! cos and sin of angle 
    cosa = cos(angle); sina = sin(angle)
    onemcosa = one - cosa; onemsina = one - sina

    ! set rotation matrix (wikipedia, also same in SOWFA)
    rot_matrix(1,1) = axis(1)**2*onemcosa+cosa;  rot_matrix(1,2) = axis(1)*axis(2)*onemcosa-axis(3)*sina;  rot_matrix(1,3) = axis(1)*axis(3)*onemcosa+axis(2)*sina
    rot_matrix(2,2) = axis(2)**2*onemcosa+cosa;  rot_matrix(2,3) = axis(2)*axis(3)*onemcosa-axis(1)*sina;  rot_matrix(2,1) = axis(2)*axis(1)*onemcosa+axis(3)*sina
    rot_matrix(3,3) = axis(3)**2*onemcosa+cosa;  rot_matrix(3,1) = axis(3)*axis(1)*onemcosa-axis(2)*sina;  rot_matrix(3,2) = axis(3)*axis(2)*onemcosa+axis(1)*sina

    do ptID = 1, this%num_blade_points
       ! shift origin to rotor_center
       this%blade_points(:, ptID, blID) = this%blade_points(:, ptID, blID) - this%rotor_center(:)

       ! rotate point in this frame of reference
       this%blade_points(:, ptID, blID) = matmul(rot_matrix, this%blade_points(:, ptID, blID))

       ! revert to origin
       this%blade_points(:, ptID, blID) = this%blade_points(:, ptID, blID) + this%rotor_center(:)
    enddo
end subroutine

subroutine distribute_forces(this, rhsx, rhsy, rhsz)
    class(actuatorLine), intent(inout) :: this
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsx, rhsy, rhsz

    integer :: j, k
    real(rkind) :: tmp_dsqsum, dsqsum, mindsqsum, maxdsqsum, tmp_distr_thrust

    mindsqsum = 1.0D30; maxdsqsum = -1.0D30
    if(this%Am_I_Active) then
      this%source = zero
      ! for each blade
      do j = 1, this%num_blades
        do k = 1, this%num_blade_points
          ! fx and fy first (C points)
          this%dsq = (this%xSmall-this%blade_points(1,k,j))**2 + (this%ySmall-this%blade_points(2,k,j))**2 + (this%zSmall-this%blade_points(3,k,j))**2
          this%dsq = this%normfactor * exp(-this%dsq*this%OneByDelSq)
          this%source(:,:,:,1) = this%source(:,:,:,1) + this%blade_forces(1, k, j)*this%dsq
          this%source(:,:,:,2) = this%source(:,:,:,2) + this%blade_forces(2, k, j)*this%dsq
          this%source(:,:,:,3) = this%source(:,:,:,3) + this%blade_forces(3, k, j)*this%dsq

          !! check if dsq sums to one
          !tmp_dsqsum = sum(this%dsq)
          !if(this%Am_I_Split) then
          !  dsqsum = p_sum(tmp_dsqsum, this%myComm)
          !else
          !  dsqsum = tmp_dsqsum
          !endif
          !mindsqsum = min(mindsqsum, dsqsum)
          !maxdsqsum = max(maxdsqsum, dsqsum)
        enddo
      enddo
      ! compute thrust from distributed forces -- this should equal turb_thrust
      tmp_distr_thrust = sum(this%source(:,:,:,1))*this%rotor_shaft(1) + &
                         sum(this%source(:,:,:,2))*this%rotor_shaft(2) + &
                         sum(this%source(:,:,:,3))*this%rotor_shaft(3)
      if(this%Am_I_Split) then
        this%distr_thrust = p_sum(tmp_distr_thrust, this%myComm)
      else
        this%distr_thrust = tmp_distr_thrust
      endif
      this%distr_thrust = this%distr_thrust / (this%invdxdy*this%invdz)
      !write(*,'(2(i3,1x),2(e19.12,1x))') nrank, this%myComm_nrank, mindsqsum, maxdsqsum
      !tmp_dsqsum = sum(this%source(:,:,:,1))
      !dsqsum = p_sum(tmp_dsqsum, this%myComm)

      rhsx(this%ist:this%iend, this%jst:this%jend, this%kst:this%kend) = rhsx(this%ist:this%iend, this%jst:this%jend, this%kst:this%kend) + this%source(:,:,:,1)
      rhsy(this%ist:this%iend, this%jst:this%jend, this%kst:this%kend) = rhsy(this%ist:this%iend, this%jst:this%jend, this%kst:this%kend) + this%source(:,:,:,2)
      rhsz(this%ist:this%iend, this%jst:this%jend, this%kst:this%kend) = rhsz(this%ist:this%iend, this%jst:this%jend, this%kst:this%kend) + this%source(:,:,:,3)
    endif
end subroutine

subroutine interp_velocity(this, ptID, blID, u, v, w, yRightHalo, zLeftHalo, zRightHalo, uloc, point_is_on_proc)
    class(actuatorLine), intent(inout) :: this
    integer,                                                    intent(in)  :: ptID, blID
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(in)  :: u, v, w
    real(rkind), dimension(this%nxLoc, this%nzLoc,   3), intent(in)         :: yRightHalo
    real(rkind), dimension(this%nxLoc, this%nyLoc+1, 3), intent(in)         :: zLeftHalo, zRightHalo
    real(rkind), dimension(3),                                  intent(out) :: uloc
    logical,                                                    intent(out) :: point_is_on_proc

    real(rkind) :: dx1, dy1, dx2, dy2, x1, y1, x2, y2, z1, invdz, blin1(3), blin2(3), acPt(3)
    real(rkind) :: f111(3), f121(3), f211(3), f221(3), f112(3), f122(3), f212(3), f222(3)
    integer     :: ind, jnd, knd, i, j, k
    logical     :: x_halo_reqd, y_halo_reqd, z_halo_reqd


    acPt = this%blade_points(:,ptID, blID)

    point_is_on_proc = .true.

    ! find x location of actuator point
    ind = -1; x_halo_reqd = .false.
    ! Check if point is between xSmall(xlen,1,1) and xRightPad
    if(acPt(1) >= this%xSmall(this%xlen,1,1) .and. acPt(1) <= this%xRightPad) then
        ind = this%xlen
        x_halo_reqd = .true.
    endif
    if(ind==-1) then
        do i=1, this%xlen-1
          if(acPt(1) >= this%xSmall(i,1,1) .and. acPt(1) <= this%xSmall(i+1,1,1)) then
              ind = i
              exit
          endif
        enddo
    endif

    ! find y location of actuator point
    jnd = -1; y_halo_reqd = .false.
    ! Check if point is between ySmall(1,ylen,1) and yRightPad
    if(acPt(2) >= this%ySmall(1,this%ylen,1) .and. acPt(2) <= this%yRightPad) then
        jnd = this%ylen
        y_halo_reqd = .true.
    endif
    if(jnd==-1) then
        do j=1, this%ylen-1
          if(acPt(2) >= this%ySmall(1,j,1) .and. acPt(2) <= this%ySmall(1,j+1,1)) then
              jnd = j
              exit
          endif
        enddo
    endif

    ! find z location of actuator point
    knd = -1; z_halo_reqd = .false.
    ! Check if point is between zLeftPad and zSmall(1,1,1)
    if(acPt(3) <= this%zSmall(1,1,1) .and. acPt(3) >= this%zLeftPad) then
        knd = 0
        z_halo_reqd = .true.
    endif
    ! Check if point is between zSmall(1,1,zlen) and zRightPad
    if(acPt(3) >= this%zSmall(1,1,this%zlen) .and. acPt(3) <= this%zRightPad) then
        knd = this%zlen
        z_halo_reqd = .true.
    endif
    if(knd==-1) then
        do k=1, this%zlen-1
          if(acPt(3) >= this%zSmall(1,1,k) .and. acPt(3) <= this%zSmall(1,1,k+1)) then
              knd = k
              exit
          endif
        enddo
    endif

    if(ind==-1 .or. jnd==-1 .or. knd==-1) then
       point_is_on_proc = .false.
       uloc = zero
       return
    endif

    if(x_halo_reqd) then
      if(y_halo_reqd) then
        if(z_halo_reqd) then
          ! T T T
          ! -- we are in x-decomp, so x-halo does not require communication
          x1 = this%xSmall(ind,1,1); y1 = this%ySmall(1,jnd,1);
          x2 = this%xRightPad;       y2 = this%yRightPad

          if(knd==0) then
              z1 = this%zLeftPad
              invdz = -this%invdz

              f111(1) = u(ind, jnd, knd+1);   f111(2) = v(ind, jnd, knd+1);   f111(3) = w(ind, jnd, knd+1)
              f211(1) = u(1,   jnd, knd+1);   f211(2) = v(1,   jnd, knd+1);   f211(3) = w(1,   jnd, knd+1)
              f221(1:3) = yRightHalo(1,   knd+1, 1:3)
              f121(1:3) = yRightHalo(ind, knd+1, 1:3)

              f112(1:3) = zLeftHalo(ind, jnd,   1:3)
              f212(1:3) = zLeftHalo(1,   jnd,   1:3)
              f222(1:3) = zLeftHalo(1,   jnd+1, 1:3)
              f122(1:3) = zLeftHalo(ind, jnd+1, 1:3)
          else ! knd == this%zlen
              z1 = this%zSmall(1,1,knd)
              invdz = this%invdz

              f111(1) = u(ind, jnd, knd);   f111(2) = v(ind, jnd, knd);   f111(3) = w(ind, jnd, knd)
              f211(1) = u(1,   jnd, knd);   f211(2) = v(1,   jnd, knd);   f211(3) = w(1,   jnd, knd)
              f221(1:3) = yRightHalo(1,   knd, 1:3)
              f121(1:3) = yRightHalo(ind, knd, 1:3)

              f112(1:3) = zRightHalo(ind, jnd,   1:3)
              f212(1:3) = zRightHalo(1,   jnd,   1:3)
              f222(1:3) = zRightHalo(1,   jnd+1, 1:3)
              f122(1:3) = zRightHalo(ind, jnd+1, 1:3)
          endif
        else
          ! T T F
          ! -- we are in x-decomp, so x-halo does not require communication
          x1 = this%xSmall(ind,1,1); y1 = this%ySmall(1,jnd,1); z1 = this%zSmall(1,1,knd)
          x2 = this%xRightPad;       y2 = this%yRightPad
          invdz = this%invdz
 
          f111(1) = u(ind, jnd, knd  ); f111(2) = v(ind, jnd, knd  ); f111(3) = w(ind, jnd, knd  )
          f112(1) = u(ind, jnd, knd+1); f112(2) = v(ind, jnd, knd+1); f112(3) = w(ind, jnd, knd+1)

          f211(1) = u(1, jnd, knd  );   f211(2) = v(1, jnd, knd  );   f211(3) = w(1, jnd, knd  )
          f212(1) = u(1, jnd, knd+1);   f212(2) = v(1, jnd, knd+1);   f212(3) = w(1, jnd, knd+1)

          f221(1:3) = yRighthalo(1, knd,   1:3)
          f222(1:3) = yRighthalo(1, knd+1, 1:3)

          f121(1:3) = yRighthalo(ind, knd,   1:3)
          f122(1:3) = yRighthalo(ind, knd+1, 1:3)
        endif
      else
        if(z_halo_reqd) then
          ! T F T
          ! -- we are in x-decomp, so x-halo does not require communication
          x1 = this%xSmall(ind,1,1); y1 = this%ySmall(1,jnd,1);
          x2 = this%xRightPad;       y2 = this%ySmall(1,jnd+1,1)

          if(knd==0) then
              z1 = this%zLeftPad
              invdz = -this%invdz

              f111(1) = u(ind, jnd,   knd+1);   f111(2) = v(ind, jnd,   knd+1);   f111(3) = w(ind, jnd,   knd+1)
              f211(1) = u(1,   jnd,   knd+1);   f211(2) = v(1,   jnd,   knd+1);   f211(3) = w(1,   jnd,   knd+1)
              f221(1) = u(1,   jnd+1, knd+1);   f221(2) = v(1,   jnd+1, knd+1);   f221(3) = w(1,   jnd+1, knd+1)
              f121(1) = u(ind, jnd+1, knd+1);   f121(2) = v(ind, jnd+1, knd+1);   f121(3) = w(ind, jnd+1, knd+1)

              f112(1:3) = zLeftHalo(ind, jnd,   1:3)
              f212(1:3) = zLeftHalo(1,   jnd,   1:3)
              f222(1:3) = zLeftHalo(1,   jnd+1, 1:3)
              f122(1:3) = zLeftHalo(ind, jnd+1, 1:3)
          else ! knd == this%zlen
              z1 = this%zSmall(1,1,knd)
              invdz = this%invdz

              f111(1) = u(ind, jnd,   knd);   f111(2) = v(ind, jnd,   knd);   f111(3) = w(ind, jnd,   knd)
              f211(1) = u(1,   jnd,   knd);   f211(2) = v(1,   jnd,   knd);   f211(3) = w(1,   jnd,   knd)
              f221(1) = u(1,   jnd+1, knd);   f221(2) = v(1,   jnd+1, knd);   f221(3) = w(1,   jnd+1, knd)
              f121(1) = u(ind, jnd+1, knd);   f121(2) = v(ind, jnd+1, knd);   f121(3) = w(ind, jnd+1, knd)

              f112(1:3) = zRightHalo(ind, jnd,   1:3)
              f212(1:3) = zRightHalo(1,   jnd,   1:3)
              f222(1:3) = zRightHalo(1,   jnd+1, 1:3)
              f122(1:3) = zRightHalo(ind, jnd+1, 1:3)
          endif
        else
          ! T F F
          ! -- we are in x-decomp, so x-halo does not require communication
          x1 = this%xSmall(ind,1,1); y1 = this%ySmall(1,jnd,  1); z1 = this%zSmall(1,1,knd)
          x2 = this%xRightPad;       y2 = this%ySmall(1,jnd+1,1)
          invdz = this%invdz
 
          f111(1) = u(ind, jnd, knd  ); f111(2) = v(ind, jnd, knd  ); f111(3) = w(ind, jnd, knd  )
          f112(1) = u(ind, jnd, knd+1); f112(2) = v(ind, jnd, knd+1); f112(3) = w(ind, jnd, knd+1)

          f211(1) = u(1, jnd, knd  );   f211(2) = v(1, jnd, knd  );   f211(3) = w(1, jnd, knd  )     !--periodic in x--
          f212(1) = u(1, jnd, knd+1);   f212(2) = v(1, jnd, knd+1);   f212(3) = w(1, jnd, knd+1)

          f221(1) = u(1, jnd+1, knd  );   f221(2) = v(1, jnd+1, knd  );   f221(3) = w(1, jnd+1, knd  )
          f222(1) = u(1, jnd+1, knd+1);   f222(2) = v(1, jnd+1, knd+1);   f222(3) = w(1, jnd+1, knd+1)

          f121(1) = u(ind, jnd+1, knd  );   f121(2) = v(ind, jnd+1, knd  );   f121(3) = w(ind, jnd+1, knd  )
          f122(1) = u(ind, jnd+1, knd+1);   f122(2) = v(ind, jnd+1, knd+1);   f122(3) = w(ind, jnd+1, knd+1)
        endif
      endif
    else
      if(y_halo_reqd) then
        if(z_halo_reqd) then
          ! F T T
          x1 = this%xSmall(ind,  1,1); y1 = this%ySmall(1,jnd,1)
          x2 = this%xSmall(ind+1,1,1); y2 = this%yRightPad
 
          if(knd==0) then
              z1 = this%zLeftPad
              invdz = -this%invdz

              f111(1) = u(ind,   jnd, knd+1);   f111(2) = v(ind,   jnd, knd+1);   f111(3) = w(ind,   jnd, knd+1)
              f211(1) = u(ind+1, jnd, knd+1);   f211(2) = v(ind+1, jnd, knd+1);   f211(3) = w(ind+1, jnd, knd+1)
              f221(1:3) = yRightHalo(ind+1, knd+1, 1:3)
              f121(1:3) = yRightHalo(ind,   knd+1, 1:3)

              f112(1:3) = zLeftHalo(ind,   jnd,   1:3)
              f212(1:3) = zLeftHalo(ind+1, jnd,   1:3)
              f222(1:3) = zLeftHalo(ind+1, jnd+1, 1:3)
              f122(1:3) = zLeftHalo(ind  , jnd+1, 1:3)
          else ! knd == this%zlen
              z1 = this%zSmall(1,1,knd)
              invdz = this%invdz

              f111(1) = u(ind,   jnd, knd);   f111(2) = v(ind,   jnd, knd);   f111(3) = w(ind,   jnd, knd)
              f211(1) = u(ind+1, jnd, knd);   f211(2) = v(ind+1, jnd, knd);   f211(3) = w(ind+1, jnd, knd)
              f221(1:3) = yRightHalo(ind+1, knd, 1:3)
              f121(1:3) = yRightHalo(ind,   knd, 1:3)

              f112(1:3) = zLeftHalo(ind,   jnd,   1:3)
              f112(1:3) = zRightHalo(ind,   jnd,   1:3)
              f212(1:3) = zRightHalo(ind+1, jnd,   1:3)
              f222(1:3) = zRightHalo(ind+1, jnd+1, 1:3)
              f122(1:3) = zRightHalo(ind  , jnd+1, 1:3)
          endif
        else
          ! F T F
          x1 = this%xSmall(ind,1,1);   y1 = this%ySmall(1,jnd,1); z1 = this%zSmall(1,1,knd)
          x2 = this%xSmall(ind+1,1,1); y2 = this%yRightPad
          invdz = this%invdz
 
          f111(1) = u(ind, jnd, knd  ); f111(2) = v(ind, jnd, knd  ); f111(3) = w(ind, jnd, knd  )
          f112(1) = u(ind, jnd, knd+1); f112(2) = v(ind, jnd, knd+1); f112(3) = w(ind, jnd, knd+1)

          f211(1) = u(ind+1, jnd, knd  );   f111(2) = v(ind+1, jnd, knd  );   f111(3) = w(ind+1, jnd, knd  )
          f212(1) = u(ind+1, jnd, knd+1);   f212(2) = v(ind+1, jnd, knd+1);   f212(3) = w(ind+1, jnd, knd+1)

          f221(1:3) = yRightHalo(ind+1, knd,  1:3)
          f222(1:3) = yRightHalo(ind+1, knd+1,1:3)

          f121(1:3) = yRightHalo(ind, knd,  1:3)
          f122(1:3) = yRightHalo(ind, knd+1,1:3)
        endif
      else
        if(z_halo_reqd) then
          !! F F T
          x1 = this%xSmall(ind,  1,1); y1 = this%ySmall(1,jnd,  1)
          x2 = this%xSmall(ind+1,1,1); y2 = this%ySmall(1,jnd+1,1)

          if(knd==0) then
              z1 = this%zLeftPad
              invdz = -this%invdz

              f111(1) = u(ind,   jnd,   knd+1);   f111(2) = v(ind,   jnd,   knd+1);   f111(3) = w(ind,   jnd,   knd+1)
              f211(1) = u(ind+1, jnd,   knd+1);   f211(2) = v(ind+1, jnd,   knd+1);   f211(3) = w(ind+1, jnd,   knd+1)
              f221(1) = u(ind+1, jnd+1, knd+1);   f221(2) = v(ind+1, jnd+1, knd+1);   f221(3) = w(ind+1, jnd+1, knd+1)
              f121(1) = u(ind,   jnd+1, knd+1);   f121(2) = v(ind,   jnd+1, knd+1);   f121(3) = w(ind,   jnd+1, knd+1)

              f112(1:3) = zLeftHalo(ind,   jnd,   1:3)
              f212(1:3) = zLeftHalo(ind+1, jnd,   1:3)
              f222(1:3) = zLeftHalo(ind+1, jnd+1, 1:3)
              f122(1:3) = zLeftHalo(ind,   jnd+1, 1:3)
          else ! knd == this%zlen
              z1 = this%zSmall(1,1,knd)
              invdz = this%invdz

              f111(1) = u(ind,   jnd,   knd);   f111(2) = v(ind,   jnd,   knd);   f111(3) = w(ind,   jnd,   knd)
              f211(1) = u(ind+1, jnd,   knd);   f211(2) = v(ind+1, jnd,   knd);   f211(3) = w(ind+1, jnd,   knd)
              f221(1) = u(ind+1, jnd+1, knd);   f221(2) = v(ind+1, jnd+1, knd);   f221(3) = w(ind+1, jnd+1, knd)
              f121(1) = u(ind,   jnd+1, knd);   f121(2) = v(ind,   jnd+1, knd);   f121(3) = w(ind,   jnd+1, knd)

              f112(1:3) = zRightHalo(ind,   jnd,   1:3)
              f212(1:3) = zRightHalo(ind+1, jnd,   1:3)
              f222(1:3) = zRightHalo(ind+1, jnd+1, 1:3)
              f122(1:3) = zRightHalo(ind,   jnd+1, 1:3)
          endif
        else
          ! F F F
          x1 = this%xSmall(ind,1,1);   y1 = this%ySmall(1,jnd,1);   z1 = this%zSmall(1,1,knd)
          x2 = this%xSmall(ind+1,1,1); y2 = this%ySmall(1,jnd+1,1);
          invdz = this%invdz
 
          f111(1) = u(ind, jnd, knd  ); f111(2) = v(ind, jnd, knd  ); f111(3) = w(ind, jnd, knd  )
          f112(1) = u(ind, jnd, knd+1); f112(2) = v(ind, jnd, knd+1); f112(3) = w(ind, jnd, knd+1)

          f211(1) = u(ind+1, jnd, knd  );   f211(2) = v(ind+1, jnd, knd  );   f211(3) = w(ind+1, jnd, knd  )
          f212(1) = u(ind+1, jnd, knd+1);   f212(2) = v(ind+1, jnd, knd+1);   f212(3) = w(ind+1, jnd, knd+1)

          f221(1) = u(ind+1, jnd+1, knd  );   f221(2) = v(ind+1, jnd+1, knd  );   f221(3) = w(ind+1, jnd+1, knd  )
          f222(1) = u(ind+1, jnd+1, knd+1);   f222(2) = v(ind+1, jnd+1, knd+1);   f222(3) = w(ind+1, jnd+1, knd+1)

          f121(1) = u(ind, jnd+1, knd  );   f121(2) = v(ind, jnd+1, knd  );   f121(3) = w(ind, jnd+1, knd  )
          f122(1) = u(ind, jnd+1, knd+1);   f122(2) = v(ind, jnd+1, knd+1);   f122(3) = w(ind, jnd+1, knd+1)
        endif
      endif
    endif

    dx1 = acPt(1) - x1; dx2 = x2 - acPt(1)
    dy1 = acPt(2) - y1; dy2 = y2 - acPt(2)
 
    ! trilinear interpolation is linear interpolation (in z) of two bilinear interpolations (in x and y)
    blin1 = this%invdxdy * (f111*dx2*dy2 + f211*dx1*dy2 + f121*dx2*dy1 + f221*dx1*dy1)
    blin2 = this%invdxdy * (f112*dx2*dy2 + f212*dx1*dy2 + f122*dx2*dy1 + f222*dx1*dy1)
    uloc = blin1 + invdz * (blin2 - blin1) * (acPt(3) - z1)

    point_is_on_proc = .true.

end subroutine

subroutine interp_clcd(this, ptID, blID, AOA, cl, cd)
    class(actuatorLine), intent(inout), target :: this
    integer, intent(in) :: ptID, blID
    real(rkind), intent(in) :: AOA
    real(rkind), intent(out) :: cl, cd

    integer :: airfIDInd, nsize, j
    real(rkind), dimension(:,:), pointer :: tablePtr
    real(rkind) :: aoaloc(1), clloc(1), cdloc(1)

    aoaloc(1) = AOA
    airfIDInd = this%airfoilIDIndex(ptID, blID)

    nsize = this%clTableSize(airfIDInd)
    tablePtr => this%clTable(1:nsize,1:nsize,airfIDInd)
    call interp(aoaloc, clloc, tablePtr)
    nullify(tablePtr)

    nsize = this%cdTableSize(airfIDInd)
    tablePtr => this%cdTable(1:nsize,1:nsize,airfIDInd)
    call interp(aoaloc, cdloc, tablePtr)
    !if(nrank==0) then
    !  write(*,'(a,3(e19.12,1x))') 'interpclcd: ', aoaloc(1), clloc, cdloc
    !endif
    nullify(tablePtr)

    cl = clloc(1); cd = cdloc(1)

end subroutine

subroutine get_blade_forces(this, u, v, w, yRightHalo, zLeftHalo, zRightHalo)
    use mpi
    use kind_parameters, only: mpirkind
    class(actuatorLine), intent(inout) :: this
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(in)  :: u, v, w
    real(rkind), dimension(this%nxLoc, this%nzLoc,   3),        intent(in)  :: yRightHalo
    real(rkind), dimension(this%nxLoc, this%nyLoc+1, 3),        intent(in)  :: zLeftHalo, zRightHalo

    integer :: j, k, ierr
    real(rkind) :: element_length, uloc(3), lift, drag
    real(rkind) :: e_rad(3), e_tan(3), e_nrm(3), urad, utan, unrm, umag, AOA, cl, cd
    real(rkind) :: dragVec(3), liftVec(3), tmp_thrust, tmp_torque, tmp_uturbavg
    logical     :: point_is_on_proc

    tmp_thrust = zero; tmp_torque = zero; tmp_uturbavg = zero

    if(this%Am_I_Active) then
      ! for each blade
      do j = 1, this%num_blades

        element_length = (this%tip_radius - this%hub_radius) / this%num_blade_points

        ! for each actuator point
        do k = 1, this%num_blade_points

          call this%interp_velocity(k, j, u, v, w, yRightHalo, zLeftHalo, zRightHalo, uloc, point_is_on_proc)

          if(.not. point_is_on_proc) then
            ! blade force at this point will be computed on a different processor
            this%blade_forcesloc(:,k,j) = zero
            this%blade_statsloc(:,k,j) = zero
          else
            ! this processor computes blade force at this point

            ! compute velocity in blade frame of reference
              ! first compute unit vectors
              ! radial direction,  radially outward if ccw; inward if cw
              e_rad = this%blade_points(:,k,j) - this%rotor_center
              if(this%clockwise_rotation) e_rad = -e_rad
              e_rad = e_rad/sqrt(sum(e_rad**2))

              ! tangential direction, facing the wind, e_rad x rotor_shaft
              e_tan(1) = e_rad(2)*this%rotor_shaft(3) - e_rad(3)*this%rotor_shaft(2)
              e_tan(2) = e_rad(3)*this%rotor_shaft(1) - e_rad(1)*this%rotor_shaft(3)
              e_tan(3) = e_rad(1)*this%rotor_shaft(2) - e_rad(2)*this%rotor_shaft(1)
              e_tan = e_tan/sqrt(sum(e_tan**2))
              
              ! normal, in the direction of wind, e_tan x e_rad
              e_nrm(1) = e_tan(2)*e_rad(3) - e_tan(3)*e_rad(2)
              e_nrm(2) = e_tan(3)*e_rad(1) - e_tan(1)*e_rad(3)
              e_nrm(3) = e_tan(1)*e_rad(2) - e_tan(2)*e_rad(1)
              e_nrm = e_nrm/sqrt(sum(e_nrm**2))

              ! now decompose local velocity in blade frame of reference
              urad = sum(uloc * e_rad);  utan = sum(uloc * e_tan);  unrm = sum(uloc * e_nrm)

            ! add rotation speed to tangential component: note e_tan is in the
            ! direction of radial motion, hence, negative sign for computing
            ! relative speed
            utan = utan - this%rotspeed*this%radial_dist(k,j)

            umag = sqrt(utan**2 + unrm**2)
            AOA = atan2(unrm, -utan)
            !AOA = atan2(utan, unrm)
            AOA = AOA - this%twistAng(k,j)
            call this%interp_clcd(k, j, AOA, cl, cd)
            lift = half * cl * umag**2 * this%chord(k,j) * element_length
            drag = half * cd * umag**2 * this%chord(k,j) * element_length

            ! drag vector in Cartesian frame
            dragVec = unrm * e_nrm + utan * e_tan
            dragVec = dragVec/sqrt(sum(dragVec**2))

            ! liftDir = radialDir x dragDir
            liftVec(1) = e_rad(2)*dragVec(3) - e_rad(3)*dragVec(2)
            liftVec(2) = e_rad(3)*dragVec(1) - e_rad(1)*dragVec(3)
            liftVec(3) = e_rad(1)*dragVec(2) - e_rad(2)*dragVec(1)

            this%blade_forcesloc(:,k,j) = - lift * liftVec(:) - drag * dragVec(:)

            ! compute sum over turbine of thrust and torque
            tmp_thrust = tmp_thrust - sum(this%blade_forcesloc(:,k,j) * this%rotor_shaft)
            tmp_torque = tmp_torque - sum(this%blade_forcesloc(:,k,j) * e_tan) * this%radial_dist(k,j)
            tmp_uturbavg = tmp_uturbavg + umag

            !if(j==1 .and. k==30) then
            !!write(100+nrank,'(7(e19.12,1x))') AOA, atan2(unrm, utan), this%radial_dist(k,j), this%twistAng(k,j), cl, cd
            !write(100+nrank,'(a,i4,1x,7(e19.12,1x))') 'e_rad   : ', nrank, e_rad
            !write(100+nrank,'(a,i4,1x,7(e19.12,1x))') 'rotshaft: ', nrank, this%rotor_shaft
            !write(100+nrank,'(a,i4,1x,7(e19.12,1x))') 'e_tan   : ', nrank, e_tan
            !write(100+nrank,'(a,i4,1x,7(e19.12,1x))') 'e_nrm   : ', nrank, e_nrm

            !write(100+nrank,'(a,i4,1x,7(e19.12,1x))') 'lift: ', nrank, lift, liftvec, lift*liftvec
            !write(100+nrank,'(a,i4,1x,7(e19.12,1x))') 'drag: ', nrank, drag, dragvec, drag*dragvec
            !write(100+nrank,'(3(i4,1x),4(e19.12,1x))') nrank, k, j, tmp_thrust, this%blade_forcesloc(:,k,j)
            !write(100+nrank,*) '---------'
            !endif

            !! if blade statistics are needed
            !this%blade_statsloc(1,k,j) = unrm;                    this%blade_statsloc(2,k,j) = utan;                  this%blade_statsloc(3,k,j) = urad
            !this%blade_statsloc(4,k,j) = umag;                    this%blade_statsloc(5,k,j) = aoa;                   this%blade_statsloc(6,k,j) = lift
            !this%blade_statsloc(7,k,j) = drag;                    this%blade_statsloc(8:10,k,j) = liftVec;            this%blade_statsloc(11:13,k,j) = dragVec
            !this%blade_statsloc(14,k,j) = this%radial_dist(k,j);  this%blade_statsloc(15:17,k,j) = this%rotor_shaft;  this%blade_statsloc(18:20,k,j) = e_nrm
            !this%blade_statsloc(21:23,k,j) = e_tan;               this%blade_statsloc(24:26,k,j) = e_rad;             this%blade_statsloc(27,k,j) = cl
            !this%blade_statsloc(28,k,j) = cd
          endif
        enddo
      enddo

      ! get values from actuator points that lie on other processors
      if(this%Am_I_Split) then
        call MPI_Allreduce(this%blade_forcesloc, this%blade_forces, 3*this%num_blade_points*this%num_blades, mpirkind, MPI_SUM, this%myComm, ierr)
        this%turb_thrust = p_sum(tmp_thrust, this%myComm)
        this%turb_torque = p_sum(tmp_torque, this%myComm)
        this%uturbavg = p_sum(tmp_uturbavg, this%myComm) / real(this%num_blades*this%num_blade_points, rkind)
        !call MPI_Allreduce(this%blade_statsloc, this%blade_stats, 28*this%num_blade_points*this%num_blades, mpirkind, MPI_SUM, this%myComm, ierr)
      else
        this%blade_forces = this%blade_forcesloc
        this%turb_thrust = tmp_thrust
        this%turb_torque = tmp_torque
        this%uturbavg =  tmp_uturbavg / real(this%num_blades*this%num_blade_points, rkind)
        !this%blade_stats = this%blade_statsloc
      endif
      !if(this%myComm_nrank==0) then
      !!  write(*,*) '---------'
      !!  do j = 1, this%num_blades
      !!   do k = 1, this%num_blade_points
      !!     write(*,'(2(i4,1x),3(e19.12,1x))') k, j, this%blade_forces(:,k,j)
      !!   enddo
      !!  enddo
      !!  write(*,*) '---------'
      !  open(13,file='blade_stats.dat',status='replace',action='write')
      !  write(13,'(a300)') 'VARIABLES="i","x","y",z","unrm","utan","urad","umag","aoa","lift","drag","liftvecx","liftvecy","liftvecz","dragvecx","gradvecy","dragvecz",&
      !                         "raddist","rotshaftx","rotshafty","rotshaftz","enrmx","enrmy","enrmz","etanx","etany","etanz","eradx","erady","eradz","cl","cd","bfx","bfy","bfz"' 
      !  do j = 1, this%num_blades
      !    write(13,'(a,i4.4,a,i6)') 'ZONE T="Blade ', j, '", F=POINT, I=', this%num_blade_points
      !    do k = 1, this%num_blade_points
      !      write(13,'(i6,1x,50(e19.12,1x))') k, this%blade_points(:,k,j), this%blade_stats(:,k,j), this%blade_forces(:,k,j)
      !    enddo
      !  enddo
      !  close(13)
      !  !k = 4; j = 1; write(*,'(i6,1x,50(e19.12,1x))') j, this%blade_points(:,k,j), this%blade_stats(:,k,j), this%blade_forces(:,k,j)
      !  !write(*,*) '---+++---'
      !endif
    endif

end subroutine

subroutine get_rotation_speed(this)
    class(actuatorLine), intent(inout), target :: this

    ! implement a controller here
    ! constant imposed rotation speed for now, so do nothing

    !if(this%Am_I_Active) then
    !endif

end subroutine

subroutine update_turbine(this, dt)
    class(actuatorLine), intent(inout), target :: this
    real(rkind), intent(in) :: dt

    integer :: j
    real(rkind) :: dAzimuth, dYaw

    if(this%Am_I_Active) then
      ! rotate blades
        ! set increment in azimuth
        dAzimuth = this%rotspeed * dt
        if( .NOT. this%clockwise_rotation) dAzimuth = -dAzimuth
  
        ! rotate blades
        do j = 1, this%num_blades
            ! rotate each blade to its correct azimuth
            call this%rotate_one_blade(j, dAzimuth)
        enddo
  
        ! update azimuth of first blade
        this%blade_azimuth = this%blade_azimuth + dAzimuth
        if(this%blade_azimuth > two*pi) then
            this%blade_azimuth = this%blade_azimuth - two*pi
        elseif(this%blade_azimuth < zero) then
            this%blade_azimuth = this%blade_azimuth + two*pi
        endif

      ! yaw rotor
        ! set increment in yaw
        dYaw = zero

        ! update turbine blades, rotor_center and shaft direction
        if(dYaw > 1.0d-10) then
          call this%yaw_turbine(dYaw)
          this%yaw_angle = this%yaw_angle + dYaw
          if(this%yaw_angle > two*pi) then
              this%yaw_angle = this%yaw_angle - two*pi
          elseif(this%yaw_angle < zero) then
              this%yaw_angle = this%yaw_angle + two*pi
          endif
        endif
    endif

end subroutine

subroutine get_RHS(this, dt, u, v, w, yRightHalo, zLeftHalo, zRightHalo, rhsxvals, rhsyvals, rhszvals, inst_val)
    class(actuatorLine),                                        intent(inout) :: this
    real(rkind),                                                intent(in)    :: dt
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(in)    :: u, v, w
    real(rkind), dimension(this%nxLoc, this%nzLoc,   3),        intent(in)    :: yRightHalo
    real(rkind), dimension(this%nxLoc, this%nyLoc+1, 3),        intent(in)    :: zLeftHalo, zRightHalo
    real(rkind), dimension(this%nxLoc, this%nyLoc, this%nzLoc), intent(inout) :: rhsxvals, rhsyvals, rhszvals
    real(rkind), dimension(8),                                  intent(out), optional  :: inst_val


    integer :: ierr

    ! update turbine point locations to account for nacelle yaw and blade rotation
    call this%update_turbine(dt)
    !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Done update_turbine"); call mpi_barrier(mpi_comm_world, ierr)

    ! compute turbine rotation speed based on a torque controller
    call this%get_rotation_speed()
    !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Done get_rotation_speed"); call mpi_barrier(mpi_comm_world, ierr)

    ! compute forces at ALM actuator points
    call this%get_blade_forces(u, v, w, yRightHalo, zLeftHalo, zRightHalo)
    !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Done get_blade_forces"); call mpi_barrier(mpi_comm_world, ierr)

    ! distribute forces from ALM points to Cartesian grid
    call this%distribute_forces(rhsxvals, rhsyvals, rhszvals)
    !call mpi_barrier(mpi_comm_world, ierr); call message(1,"Done distribute_forces"); call mpi_barrier(mpi_comm_world, ierr)

    ! some turbine-level quantities
    inst_val(1) = abs(this%turb_thrust)
    inst_val(2) = abs(this%distr_thrust)
    inst_val(3) = abs(this%turb_torque)
    inst_val(4) = abs(this%rotspeed)
    inst_val(5) = abs(this%turb_torque*this%rotspeed)
    inst_val(6) = abs(this%uturbavg)
    inst_val(7) = abs(this%uturbavg**2)
    inst_val(8) = abs(this%uturbavg**3)

end subroutine 

!subroutine write_turbineData(this)
!    class(actuatorLine),                                        intent(inout) :: this
!
!end subroutine

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
      y(i) = xyTable(j,2)
    enddo

end subroutine

subroutine interp(x, y, xyTable)
    !class(actuatorLine), intent(inout) :: this
    real(rkind), intent(in),  dimension(:)   :: x
    real(rkind), intent(out), dimension(:)   :: y
    real(rkind), intent(in),  dimension(:,:) :: xyTable

    integer :: i, locator(1), imin, i1, i2

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
       y(i) = xyTable(i1,2) + (xyTable(i2,2)-xyTable(i1,2))/(xyTable(i2,1)-xyTable(i1,1))*(x(i) - xyTable(i1,1))
    enddo

end subroutine

end module
