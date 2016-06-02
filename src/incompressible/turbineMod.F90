module turbineMod
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
    public :: TurbineArray

    real(rkind) :: epsilon_sq, eps_pi_fac
    real(rkind), parameter :: degrees_to_radians = pi/180.0_rkind

    integer, dimension(:), allocatable :: ist, iend, jst, jend, kst, kend, kstE, kendE  ! start and end indices for cube of cells that can be potentially affected by each turbine
    real(rkind), dimension(:,:,:,:), allocatable :: dist_sq, distE_sq, x_cloud, y_cloud, zC_cloud, zE_cloud

    type :: TurbineArray
        integer :: myProc
        integer, dimension(:), allocatable  :: xst, xen, yst, yen
        integer :: nTurbines
        type(decomp_info), pointer :: gpC, sp_gpC, gpE, sp_gpE
        type(spectral), pointer :: spectC, spectE
        integer :: myLeftNeigh, myRightNeigh, myTopNeigh, myBotNeigh
 
        integer, dimension(:),   allocatable :: num_cells_cloud                                     ! total number of cells in the cubic cloud around a turbine on this processor
        integer, dimension(:),   allocatable :: num_blades  ! number of blades
        integer, dimension(:,:), allocatable :: num_blade_points  ! number of actuator points on each blade
        real(rkind), dimension(:), allocatable :: yaw_angle, blade_azimuth, nacelle_width, hub_radius, tip_radius, turb_thrust, turb_torque, rotspeed
        real(rkind), dimension(:,:), allocatable :: rotor_center, turbLoc, rotor_shaft
        real(rkind), dimension(:,:,:,:), allocatable :: blade_points  ! number of actuator points on each blade
        real(rkind), dimension(:,:,:,:), allocatable :: blade_forces  ! forces at actuator points
        logical, dimension(:), allocatable :: clockwise_rotation

        real(rkind), dimension(:,:,:), pointer :: fx, fy, fz
        real(rkind), dimension(:,:,:,:), allocatable :: rbuffC, rbuffE
        complex(rkind), dimension(:,:,:,:), allocatable :: cbuffC, cbuffE

    contains

        procedure :: init
        procedure :: destroy
        procedure :: getForceRHS 
        procedure, private :: distribute_forces
        procedure, private :: get_blade_forces 
        procedure, private :: rotate_one_blade 
        procedure, private :: yaw_turbine
        procedure, private :: get_extents
        procedure, private :: interp_airfoil_props
        procedure, private :: interp_velocity
        procedure, private :: interp_clcd
        procedure, private :: get_rotation_speed
        procedure, private :: update_turbines

    end type

contains

subroutine init(this, inputFile, gpC, gpE, spectC, spectE, mesh, dx, dy, dz)
    class(TurbineArray), intent(inout), target :: this
    character(len=*), intent(in) :: inputFile
    type(spectral), target :: spectC, spectE
    type(decomp_info), target :: gpC, gpE
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), intent(in) :: dx, dy, dz

    integer :: i, j, k
    real(rkind) :: element_length, radial_dist, projection_radius
    real(rkind) :: xmin, xmax, ymin, ymax, zmin, zmax 

    print*, inputFile
    this%gpC => gpC
    this%spectC => this%spectC
    this%sp_gpC => this%spectC%spectdecomp

    this%gpE => gpE
    this%spectE => this%spectE
    this%sp_gpE => this%spectE%spectdecomp

    call GracefulExit("Wind Turbine stuff is incomplete", 423)

    allocate(this%num_cells_cloud(this%nTurbines), this%num_blades(this%nTurbines))
    allocate(ist(this%nTurbines), iend(this%nTurbines),  jst(this%nTurbines),  jend(this%nTurbines))
    allocate(kst(this%nTurbines), kend(this%nTurbines), kstE(this%nTurbines), kendE(this%nTurbines))

    allocate(this%cbuffC(this%sp_gpC%ysz(1), this%sp_gpC%ysz(2), this%sp_gpC%ysz(3), 1))
    allocate(this%cbuffE(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), this%sp_gpE%ysz(3), 1))

    allocate(this%rbuffC(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3), 2))
    allocate(this%rbuffE(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3), 1))

    this%fx => this%rbuffC(:,:,:,1); this%fy => this%rbuffC(:,:,:,2); this%fz => this%rbuffE(:,:,:,1)

    ! set number of turbines
    this%nTurbines = 1;

    allocate(this%num_blades(this%nTurbines),    this%yaw_angle(this%nTurbines))
    allocate(this%blade_azimuth(this%nTurbines), this%rotor_center(3, this%nTurbines))
    allocate(this%turbLoc(3, this%nTurbines),    this%nacelle_width(this%nTurbines))
    allocate(this%tip_radius(this%nTurbines),    this%hub_radius(this%nTurbines))
    allocate(this%turb_thrust(this%nTurbines),   this%turb_torque(this%nTurbines))
    allocate(this%rotor_shaft(3, this%nTurbines), this%rotspeed(this%nTurbines))
    allocate(this%clockwise_rotation(this%nTurbines))

    ! factor for distributing blade forces to grid points
    epsilon_sq = (min(1.5d0*dx, 1.5d0*dy, 1.5d0*dz))**2 ! adjust this factor later
    eps_pi_fac = (pi * epsilon_sq)**1.5d0

    do i = 1, this%nTurbines

      ! set number of blades for each turbine (identical for now)
      this%num_blades(i) = 3

      ! initialize turbines with initial yaw zero (facing positive x direction)
      this%yaw_angle(i) = zero * degrees_to_radians

      ! initialize first blade to zero azimuth; others are equally spaced in the 360 degree space
      this%blade_azimuth(i) = zero * degrees_to_radians

      ! rotor center is offset from (xLoc, yLoc, zLoc) by Nacelle width
      this%rotor_center(1,i) = this%turbLoc(1,i) - this%nacelle_width(i)
      this%rotor_center(2,i) = this%turbLoc(2,i); this%rotor_center(3,i) = this%turbLoc(3,i)

      this%rotor_shaft(:,i) = this%turbLoc(:,i) - this%rotor_center(:,i)

      ! set initial rotation speed of the turbines
      this%rotspeed(i) = 8.0d0*two*pi/60.0_rkind     ! 8 RPM

      allocate(this%num_blade_points(this%num_blades(i), this%nTurbines))

      do j = 1, this%num_blades(i)
        ! set number of actuator points along each blade (identical for now)
        this%num_blade_points(j,i) = 40

        ! length of each actuator line element
        element_length = (this%tip_radius(i) - this%hub_radius(i)) / this%num_blade_points(j,i)

        allocate(this%blade_points(3, this%num_blade_points(j,i), this%num_blades(i), this%nTurbines))
        allocate(this%blade_forces(3, this%num_blade_points(j,i), this%num_blades(i), this%nTurbines))

        do k = 1, this%num_blade_points(j,i)
          ! set blade points beginning from the hub and going radially outward
          radial_dist = this%hub_radius(i) + element_length * (real(k-1, rkind) + half)

          ! initialize each blade with initial azimuth zero (vertically upwards)
          this%blade_points(:,k,j,i) = this%rotor_center(:,i)
          this%blade_points(3,k,j,i) = this%rotor_center(3,i) + radial_dist
        enddo

        ! rotate each blade to its correct azimuth
        call this%rotate_one_blade(i, j, this%blade_azimuth(i) + real(j-1,rkind)*two*pi/real(this%num_blades(i), rkind))
      enddo

      ! now rotate all blades so that yaw angle is correct
      call this%yaw_turbine(i, this%yaw_angle(i))

      ! compute ist, iend, jst, jend, kst, kend based on turbine location, blade radius, projection radius and processor extents
      projection_radius = sqrt(epsilon_sq * log(1.0D3))    ! distance where influence of actuator points reduces to 0.001 of the max value
      radial_dist = this%tip_radius(i) + projection_radius

      xmin = this%rotor_center(1,i) - radial_dist; xmax = this%rotor_center(1,i) + radial_dist
      ymin = this%rotor_center(2,i) - radial_dist; ymax = this%rotor_center(2,i) + radial_dist
      zmin = this%rotor_center(3,i) - radial_dist; zmax = this%rotor_center(3,i) + radial_dist

      call this%get_extents(ist(i), iend(i), xmin, xmax, mesh(:,1,1,1))
      call this%get_extents(jst(i), jend(i), ymin, ymax, mesh(1,:,1,2))
      call this%get_extents(kst(i), kend(i), zmin, zmax, mesh(1,1,:,3))
      call this%get_extents(kstE(i), kendE(i), zmin, zmax, mesh(1,1,:,3))     ! this needs a zE or meshE array - not correct right now

      this%num_cells_cloud(i) = (iend(i)-ist(i)+1) * (jend(i)-jst(i)+1) * (kend(i)-kst(i)+1)

      if(this%num_cells_cloud(i) > 0) then
          allocate(dist_sq (ist(i):iend(i), jst(i):jend(i), kst(i):kend(i),   i))
          allocate(distE_sq(ist(i):iend(i), jst(i):jend(i), kstE(i):kendE(i), i))
          allocate( x_cloud(ist(i):iend(i), jst(i):jend(i), kst(i):kend(i),   i))
          allocate( y_cloud(ist(i):iend(i), jst(i):jend(i), kst(i):kend(i),   i))
          allocate(zC_cloud(ist(i):iend(i), jst(i):jend(i), kst(i):kend(i),   i))
          allocate(zE_cloud(ist(i):iend(i), jst(i):jend(i), kstE(i):kendE(i), i))
    
           x_cloud(ist(i):iend(i), jst(i):jend(i), kst(i):kend(i),   i) = mesh(ist(i):iend(i), jst(i):jend(i), kst(i):kend(i), 1)
           y_cloud(ist(i):iend(i), jst(i):jend(i), kst(i):kend(i),   i) = mesh(ist(i):iend(i), jst(i):jend(i), kst(i):kend(i), 2)
          zC_cloud(ist(i):iend(i), jst(i):jend(i), kst(i):kend(i),   i) = mesh(ist(i):iend(i), jst(i):jend(i), kst(i):kend(i), 3)
          zE_cloud(ist(i):iend(i), jst(i):jend(i), kstE(i):kendE(i), i) = zero ! is there a zE array?
      endif
    enddo

end subroutine


subroutine destroy(this)
    class(TurbineArray), intent(inout) :: this
    nullify(this%gpC, this%gpE, this%spectC, this%sp_gpC, this%fx, this%fy, this%fz)
    deallocate(kst, kend, jst, jend, ist, iend, kstE, kendE, this%num_cells_cloud)
    deallocate(this%cbuffC, this%cbuffE, this%num_cells_cloud, x_cloud, y_cloud, zC_cloud, zE_cloud)
    deallocate(this%blade_points, this%blade_forces)
end subroutine

subroutine get_extents(this, ist, iend, xmin, xmax, procmesh)
    class(TurbineArray), intent(in) :: this
    integer, intent(out) :: ist, iend
    real(rkind), intent(in) :: xmin, xmax
    real(rkind), dimension(:), intent(in) :: procmesh

    integer :: nloc, ind

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
end subroutine

subroutine yaw_turbine(this, turbID, angle)
    class(TurbineArray), intent(inout) :: this
    integer, intent(in) :: turbID
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
    this%rotor_center(:,turbID) = this%rotor_center(:,turbID) - this%turbLoc(:,turbID)
    this%rotor_center(:,turbID) = matmul(rot_matrix, this%rotor_center(:,turbID))
    this%rotor_center(:,turbID) = this%rotor_center(:,turbID) + this%turbLoc(:,turbID)

    ! update shaft direction
    this%rotor_shaft(:,turbID) = this%turbLoc(:,turbID) - this%rotor_center(:,turbID)

    ! next rotate each blade
    do blID = 1, this%num_blades(turbID)
      do ptID = 1, this%num_blade_points(blID, turbID)
         ! shift origin to turbLoc
         this%blade_points(:, ptID, blID, turbID) = this%blade_points(:, ptID, blID, turbID) - this%turbLoc(:,turbID)

         ! rotate point in this frame of reference
         this%blade_points(:, ptID, blID, turbID) = matmul(rot_matrix, this%blade_points(:, ptID, blID, turbID))

         ! revert to origin
         this%blade_points(:, ptID, blID, turbID) = this%blade_points(:, ptID, blID, turbID) + this%turbLoc(:,turbID)
      enddo
    enddo
end subroutine

subroutine rotate_one_blade(this, turbID, blID, angle)
    class(TurbineArray), intent(inout) :: this
    integer, intent(in) :: turbID, blID
    real(rkind), intent(in) :: angle

    real(rkind) :: axis(3), cosa, sina, onemcosa, onemsina, rot_matrix(3,3)
    integer :: ptID

    ! set unit vector in the direction of axis of rotation
    axis(:) = this%rotor_shaft(:,turbID)
    axis = axis/sqrt(sum(axis**2))

    ! cos and sin of angle 
    cosa = cos(angle); sina = sin(angle)
    onemcosa = one - cosa; onemsina = one - sina

    ! set rotation matrix (wikipedia, also same in SOWFA)
    rot_matrix(1,1) = axis(1)**2*onemcosa+cosa;  rot_matrix(1,2) = axis(1)*axis(2)*onemcosa-axis(3)*sina;  rot_matrix(1,3) = axis(1)*axis(3)*onemcosa+axis(2)*sina
    rot_matrix(2,2) = axis(2)**2*onemcosa+cosa;  rot_matrix(2,3) = axis(2)*axis(3)*onemcosa-axis(1)*sina;  rot_matrix(2,1) = axis(2)*axis(1)*onemcosa+axis(3)*sina
    rot_matrix(3,3) = axis(3)**2*onemcosa+cosa;  rot_matrix(3,1) = axis(3)*axis(1)*onemcosa-axis(2)*sina;  rot_matrix(3,2) = axis(3)*axis(2)*onemcosa+axis(1)*sina

    do ptID = 1, this%num_blade_points(blID, turbID)
       ! shift origin to rotor_center
       this%blade_points(:, ptID, blID, turbID) = this%blade_points(:, ptID, blID, turbID) - this%rotor_center(:,turbID)

       ! rotate point in this frame of reference
       this%blade_points(:, ptID, blID, turbID) = matmul(rot_matrix, this%blade_points(:, ptID, blID, turbID))

       ! revert to origin
       this%blade_points(:, ptID, blID, turbID) = this%blade_points(:, ptID, blID, turbID) + this%rotor_center(:,turbID)
    enddo
end subroutine

subroutine distribute_forces(this)
    class(TurbineArray), intent(inout) :: this

    integer :: i, j, k

    this%fx = zero; this%fy = zero; this%fz = zero

    ! for each turbine
    do i = 1, this%nTurbines
      if(this%num_cells_cloud(i) > 0) then
        ! for each blade
        do j = 1, this%num_blades(i)
          do k = 1, this%num_blade_points(j,i)
            ! fx and fy first (C points)
            dist_sq(:,:,:,i) = (x_cloud(:,:,:,i)  - this%blade_points(1,k,j,i))**2 &
                             + (y_cloud(:,:,:,i)  - this%blade_points(2,k,j,i))**2 &
                             + (zC_cloud(:,:,:,i) - this%blade_points(3,k,j,i))**2
            this%fx(ist(i):iend(i), jst(i):jend(i), kst(i):kend(i)) = &
            this%fx(ist(i):iend(i), jst(i):jend(i), kst(i):kend(i)) + this%blade_forces(1,k,j,i) * exp(-dist_sq(:,:,:,i)/epsilon_sq) / eps_pi_fac

            this%fy(ist(i):iend(i), jst(i):jend(i), kst(i):kend(i)) = &
            this%fy(ist(i):iend(i), jst(i):jend(i), kst(i):kend(i)) + this%blade_forces(2,k,j,i) * exp(-dist_sq(:,:,:,i)/epsilon_sq) / eps_pi_fac

            ! then fz (E points)
            distE_sq(:,:,:,i) = (x_cloud(:,:,:,i)  - this%blade_points(1,k,j,i))**2 &
                              + (y_cloud(:,:,:,i)  - this%blade_points(2,k,j,i))**2 &
                              + (zE_cloud(:,:,:,i) - this%blade_points(3,k,j,i))**2
            this%fz(ist(i):iend(i), jst(i):jend(i), kstE(i):kendE(i)) = &
            this%fz(ist(i):iend(i), jst(i):jend(i), kstE(i):kendE(i)) + this%blade_forces(3,k,j,i) * exp(-distE_sq(:,:,:,i)/epsilon_sq) / eps_pi_fac
            
          enddo
        enddo
      endif
    enddo

end subroutine

subroutine interp_airfoil_props(this, i, radial_dist, twistAng, chord, airfoilID)
    class(TurbineArray), intent(inout) :: this
    integer, intent(in) :: i
    real(rkind), intent(in) :: radial_dist
    real(rkind), intent(out) :: twistAng, chord
    integer, intent(out) :: airfoilID

end subroutine

subroutine interp_velocity(this, uloc, blPoint)
    class(TurbineArray), intent(inout) :: this
    real(rkind), dimension(3), intent(out) :: uloc
    real(rkind), dimension(3), intent(in) :: blPoint

end subroutine

subroutine interp_clcd(this, i, j, AOA, cl, cd)
    class(TurbineArray), intent(inout) :: this
    integer, intent(in) :: i, j
    real(rkind), intent(in) :: AOA
    real(rkind), intent(out) :: cl, cd

end subroutine
subroutine get_blade_forces(this)
    class(TurbineArray), intent(inout) :: this

    integer :: i, j, k, airfoilID
    real(rkind) :: element_length, radial_dist, twistAng, chord, uloc(3), lift, drag
    real(rkind) :: e_rad(3), e_tan(3), e_nrm(3), urad, utan, unrm, umag, AOA, cl, cd
    real(rkind) :: dragVec(3), liftVec(3)

    ! for each turbine
    do i = 1, this%nTurbines

      if(this%num_cells_cloud(i) > 0) then

        ! for each blade
        do j = 1, this%num_blades(i)

          element_length = (this%tip_radius(i) - this%hub_radius(i)) / this%num_blade_points(j,i)

          ! for each actuator point
          do k = 1, this%num_blade_points(j,i)

            radial_dist = this%hub_radius(i) + element_length * (real(k-1, rkind) + half)
            call this%interp_airfoil_props(i, radial_dist, twistAng, chord, airfoilID)

            call this%interp_velocity(uloc, this%blade_points(:,k,j,i))

            ! compute velocity in blade frame of reference
              ! first compute unit vectors
              ! radial direction,  radially outward if cw; inward if ccw
              e_rad = this%blade_points(:,k,j,i) - this%rotor_center(:,i)
              if(this%clockwise_rotation(i)) e_rad = -e_rad
              e_rad = e_rad/sqrt(sum(e_rad**2))

              ! tangential direction, facing the wind, e_rad x rotor_shaft
              e_tan(1) = e_rad(2)*this%rotor_shaft(3,i) - e_rad(3)*this%rotor_shaft(2,i)
              e_tan(2) = e_rad(3)*this%rotor_shaft(1,i) - e_rad(1)*this%rotor_shaft(3,i)
              e_tan(3) = e_rad(1)*this%rotor_shaft(2,i) - e_rad(2)*this%rotor_shaft(1,i)
              e_tan = e_tan/sqrt(sum(e_tan**2))
              
              ! normal, in the direction of wind, 
              e_nrm(1) = e_tan(2)*e_rad(3) - e_tan(3)*e_rad(2)
              e_nrm(2) = e_tan(3)*e_rad(1) - e_tan(1)*e_rad(3)
              e_nrm(3) = e_tan(1)*e_rad(2) - e_tan(2)*e_rad(1)
              e_nrm = e_nrm/sqrt(sum(e_nrm**2))

              ! now decompose local velocity in blade frame of reference
              urad = sum(uloc * e_rad);  utan = sum(uloc * e_tan);  unrm = sum(uloc * e_nrm)

            ! add rotation speed to tangential component
            utan = utan + this%rotspeed(i)*radial_dist

            umag = sqrt(utan**2 + unrm**2)
            AOA = atan2(unrm, utan)
            AOA = AOA - twistAng
            call this%interp_clcd(i, j, AOA, cl, cd)
            lift = half * cl * umag**2 * chord * element_length
            drag = half * cd * umag**2 * chord * element_length

            ! drag vector in Cartesian frame
            dragVec = urad * e_nrm + utan * e_tan

            ! liftDir = dragDir x radialDir
            liftVec(1) = dragVec(2)*e_rad(3) - dragVec(3)*e_rad(2)
            liftVec(2) = dragVec(3)*e_rad(1) - dragVec(1)*e_rad(3)
            liftVec(3) = dragVec(1)*e_rad(2) - dragVec(2)*e_rad(1)

            this%blade_forces(:,k,j,i) = - lift * liftVec(:) - drag * dragVec(:)

            ! compute sum over turbine of thrust and torque
            this%turb_thrust(i) = this%turb_thrust(i) - sum(this%blade_forces(:,k,j,i) * this%rotor_shaft(:,i))
            this%turb_torque(i) = this%turb_torque(i) + sum(this%blade_forces(:,k,j,i) * e_tan) * radial_dist
          enddo
        enddo
      endif
    enddo
end subroutine

subroutine get_rotation_speed(this)
    class(TurbineArray), intent(inout), target :: this

    integer :: i

    do i = 1, this%nTurbines
      if(this%num_cells_cloud(i) > 0) then
        ! implement a controller here
        ! constant imposed rotation speed for now, so do nothing
      endif
    enddo

end subroutine

subroutine update_turbines(this, dt)
    class(TurbineArray), intent(inout), target :: this
    real(rkind), intent(in) :: dt

    integer :: i, j
    real(rkind) :: dAzimuth, dYaw

    do i = 1, this%nTurbines
      if(this%num_cells_cloud(i) > 0) then
        ! rotate blades
          ! set increment in azimuth
          dAzimuth = this%rotspeed(i) * dt
          if( .NOT. this%clockwise_rotation(i)) dAzimuth = -dAzimuth
  
          ! rotate blades
          do j = 1, this%num_blades(i)
              ! rotate each blade to its correct azimuth
              call this%rotate_one_blade(i, j, dAzimuth)
          enddo
  
          ! update azimuth of first blade
          this%blade_azimuth(i) = this%blade_azimuth(i) + dAzimuth
          if(this%blade_azimuth(i) > two*pi) then
              this%blade_azimuth(i) = this%blade_azimuth(i) - two*pi
          elseif(this%blade_azimuth(i) < zero) then
              this%blade_azimuth(i) = this%blade_azimuth(i) + two*pi
          endif

        ! yaw rotor
          ! set increment in yaw
          dYaw = zero

          ! update turbine blades, rotor_center and shaft direction
          if(dYaw > 1.0d-10) then
            call this%yaw_turbine(i, dYaw)
            this%yaw_angle(i) = this%yaw_angle(i) + dYaw
            if(this%yaw_angle(i) > two*pi) then
                this%yaw_angle(i) = this%yaw_angle(i) - two*pi
            elseif(this%yaw_angle(i) < zero) then
                this%yaw_angle(i) = this%yaw_angle(i) + two*pi
            endif
          endif
      endif
    enddo

end subroutine

subroutine getForceRHS(this, dt, u, v, wC, urhs, vrhs, wrhs)
    class(TurbineArray), intent(inout), target :: this
    real(rkind),                                                                         intent(in) :: dt
    real(rkind),    dimension(this%gpC%xsz(1),   this%gpC%xsz(2),   this%gpC%xsz(3)),    intent(in) :: u, v, wC
    complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout) :: urhs, vrhs
    complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: wrhs 

    complex(rkind), dimension(:,:,:), pointer :: fChat, fEhat

    fChat => this%cbuffC(:,:,:,1); fEhat => this%cbuffE(:,:,:,1)

    call this%get_rotation_speed()          ! compute turbine rotation speed based on a torque controller
    call this%update_turbines(dt)           ! update turbine point locations to account for nacelle yaw and blade rotation
    call this%get_blade_forces()            ! compute forces at ALM actuator points
    call this%distribute_forces()           ! distribute forces from ALM points to Cartesian grid

    call this%spectC%fft(this%fx,fChat)
    urhs = urhs + fChat

    call this%spectC%fft(this%fy,fChat)
    vrhs = vrhs + fChat

    call this%spectC%fft(this%fz,fEhat)
    wrhs = wrhs + fEhat

    nullify(fChat, fEhat)

end subroutine 

end module
