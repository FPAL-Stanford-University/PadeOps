module particlesMod
    use kind_parameters, only: rkind, clen
    use decomp_2d,       only: nrank
    use constants,       only: zero, one, two, third, half
    use exits,           only: GracefulExit

    implicit none

    private
    public :: particles
    
    type :: particles
        integer :: np = 1, initflag_pos = 1, initflag_vel = 1
        real(rkind), allocatable, dimension(:,:) :: pos, vel, acc, fluidvel
        real(rkind), allocatable, dimension(:)   :: xline, yline, zline
        real(rkind) :: dx, dy, dz, Lx, Ly, Lz
        integer     :: nx, ny, nz
        logical :: tracer_particles = .true., inertial_particles = .false.
        logical :: periodicx = .true., periodicy = .true., periodicz = .false.

    contains
        procedure          :: init
        procedure          :: destroy
        procedure          :: update 
        procedure          :: write_viz
        procedure, private :: update_acc
        procedure, private :: update_vel
        procedure, private :: interp_fluidvel_to_particlepos
        procedure, private :: get_interp_factors
    end type

contains

subroutine init(this, inputfile, Lx, Ly, Lz, nx, ny, nz, u, v, w, xline, yline, zline, periodicx, periodicy, periodicz)
    class(particles), intent(inout) :: this
    character(len=clen),           intent(in) :: inputfile
    real(rkind),                   intent(in) :: Lx, Ly, Lz
    integer,                       intent(in) :: nx, ny, nz
    real(rkind), dimension(:,:,:), intent(in) :: u, v, w
    real(rkind), dimension(:),     intent(in) :: xline, yline, zline
    logical,                       intent(in) :: periodicx, periodicy, periodicz
    
    integer :: np = 1, initflag_pos = 1, initflag_vel = 1
    integer :: iounit, i, j, k, ip, np_x, np_y, np_z
    real(rkind) :: dx, dy, dz
    logical :: tracer_particles = .true., inertial_particles = .false.
    real(rkind), allocatable :: particle_locs(:, :)

    namelist /PARTICLES/ np, initflag_pos, initflag_vel, &
                         tracer_particles, inertial_particles

    ioUnit = 10
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PARTICLES)
    close(ioUnit)
    
    this%np = np
    this%inertial_particles = inertial_particles
    this%tracer_particles = tracer_particles
    this%periodicx = periodicx
    this%periodicy = periodicy
    this%periodicz = periodicz
    this%Lx = Lx; this%Ly = Ly; this%Lz = Lz

    if(nrank > 0) then
      call GracefulExit("More than one processor not permitted as of now", 11)
    endif

    if(this%inertial_particles) then
      call GracefulExit("Inertial particles not implemented as of now", 11)
    endif

    allocate(this%pos(this%np, 3))
    allocate(this%vel(this%np, 3))
    allocate(this%acc(this%np, 3))
    allocate(this%fluidvel(this%np, 3))
    allocate(this%xline(nx))
    allocate(this%yline(ny))
    allocate(this%zline(nz))

    this%xline = xline
    this%yline = yline
    this%zline = zline
    this%dx = xline(2) - xline(1)
    this%dy = yline(2) - yline(1)
    this%dz = zline(2) - zline(1)
    this%nx = nx
    this%ny = ny
    this%nz = nz

    if(initflag_pos==1) then
      ! random initialization
      if(this%np==1) then
        this%pos(1, 1) = 1
        this%pos(1, 2) = 6
        this%pos(1, 3) = 3
      else
         call GracefulExit("initflag_pos=1 (random) not implemented yet", 11)
      endif
    elseif(initflag_pos==2) then
      ! uniformly spaced initial locations
      np_x = nint(this%np**(third))
      np_y = np_x; np_z = np_x
      dx = Lx/nx; dy = Ly/ny; dz = Lz/nz
      ip = 0
      do k = 1, np_z
       do j = 1, np_y
        do i = 1, np_x
          ip = ip+1
          this%pos(ip, 1) = (i-half) * dx
          this%pos(ip, 2) = (j-half) * dy
          this%pos(ip, 3) = (k-half) * dz
        enddo
       enddo
      enddo
    elseif(initflag_pos==3) then
      ! read from a file
      call GracefulExit("initflag_pos=3 (read from a file) not implemented yet", 11)
    elseif(initflag_pos==4) then
      ! put particles at custom locations
        call initialize_particles(particle_locs, np)
        this%np = size(particle_locs, 2)

        do ip = 1, this%np
          this%pos(ip,1) = particle_locs(1,ip)
          this%pos(ip,2) = particle_locs(2,ip)
          this%pos(ip,3) = particle_locs(3,ip)
        end do
   endif

    if(initflag_vel==1) then
      ! random initialization
      call GracefulExit("initflag_vel=1 (random) not implemented yet", 11)
    elseif(initflag_vel==2) then
      ! interpolation from fluid velocities
      call this%interp_fluidvel_to_particlepos(u, v, w)
      this%vel = this%fluidvel
    endif

end subroutine

subroutine destroy(this)
    class(particles), intent(inout) :: this

    deallocate(this%zline)
    deallocate(this%yline)
    deallocate(this%xline)
    deallocate(this%fluidvel)
    deallocate(this%acc)
    deallocate(this%vel)
    deallocate(this%pos)

end subroutine

subroutine update(this, dt, u, v, w)
    class(particles),              intent(inout) :: this
    real(rkind),                   intent(in)    :: dt
    real(rkind), dimension(:,:,:), intent(in)    :: u, v, w

    call this%interp_fluidvel_to_particlepos(u, v, w)
    call this%update_acc(u, v, w)
    call this%update_vel(u, v, w)

    ! update positions of particles
    this%pos = this%pos + dt * this%vel

    ! apply boundary conditions
    ! along x
    if(this%periodicx) then
        where(this%pos(:,1) > this%Lx) this%pos(:,1) = this%pos(:,1) - this%Lx
        where(this%pos(:,1) < zero   ) this%pos(:,1) = this%pos(:,1) + this%Lx
    else
        where(this%pos(:,1) > this%Lx) this%pos(:,1) = two * this%Lx - this%pos(:,1)
        where(this%pos(:,1) < zero   ) this%pos(:,1) = -this%pos(:,1)
        !! velocity has to be changed: normal velocity component will be
        !reversed, tangential component will remain unchanged
    endif
    if( (maxval(this%pos(:,1)) > this%Lx) .or. (minval(this%pos(:,1)) < zero) ) then
        call GracefulExit("Time step is most likely too large. Check details.",11)
    endif

    ! along y
    if(this%periodicy) then
        where(this%pos(:,2) > this%Ly) this%pos(:,2) = this%pos(:,2) - this%Ly
        where(this%pos(:,2) < zero   ) this%pos(:,2) = this%pos(:,2) + this%Ly
    else
        where(this%pos(:,2) > this%Ly) this%pos(:,2) = two * this%Ly - this%pos(:,2)
        where(this%pos(:,2) < zero   ) this%pos(:,2) = -this%pos(:,2)
    endif
    if( (maxval(this%pos(:,2)) > this%Ly) .or. (minval(this%pos(:,2)) < zero) ) then
        call GracefulExit("Time step is most likely too large. Check details.",11)
    endif

    ! along z
    if(this%periodicz) then
        where(this%pos(:,3) > this%Lz) this%pos(:,3) = this%pos(:,3) - this%Lz
        where(this%pos(:,3) < zero   ) this%pos(:,3) = this%pos(:,3) + this%Lz
    else
        where(this%pos(:,3) > this%Lz) this%pos(:,3) = two * this%Lz - this%pos(:,3)
        where(this%pos(:,3) < zero   ) this%pos(:,3) = -this%pos(:,3)
    endif
    if( (maxval(this%pos(:,3)) > this%Lz) .or. (minval(this%pos(:,3)) < zero) ) then
        call GracefulExit("Time step is most likely too large. Check details.",11)
    endif

end subroutine

subroutine update_vel(this, u, v, w)
    class(particles),              intent(inout) :: this
    real(rkind), dimension(:,:,:), intent(in)    :: u, v, w

    if(this%tracer_particles) then
        this%vel = this%fluidvel
    elseif(this%inertial_particles) then
        call GracefulExit("Inertial particles not implemented yet",11)
    endif

end subroutine

subroutine interp_fluidvel_to_particlepos(this, u, v, w)
    class(particles),              intent(inout) :: this
    real(rkind), dimension(:,:,:), intent(in)    :: u, v, w

    integer :: ip, il, jl, kl, ilp1, jlp1, klp1
    real(rkind) :: facx, facy, facz, onemfacx, onemfacy, onemfacz

    do ip = 1, this%np
      ! along x
      call this%get_interp_factors(this%pos(ip,1), this%dx, this%xline, this%nx, il, ilp1, facx, onemfacx)

      ! along y
      call this%get_interp_factors(this%pos(ip,2), this%dy, this%yline, this%ny, jl, jlp1, facy, onemfacy)

      ! along z
      call this%get_interp_factors(this%pos(ip,3), this%dz, this%zline, this%nz, kl, klp1, facz, onemfacz)

      ! now interpolate all three velocity components
      this%fluidvel(ip,1) = u(il,   jl,   kl  ) *     facx *     facy *     facz + &
                            u(ilp1, jl,   kl  ) * onemfacx *     facy *     facz + &
                            u(ilp1, jlp1, kl  ) * onemfacx * onemfacy *     facz + &
                            u(il  , jlp1, kl  ) *     facx * onemfacy *     facz + &
                            u(il,   jl,   klp1) *     facx *     facy * onemfacz + &
                            u(ilp1, jl,   klp1) * onemfacx *     facy * onemfacz + &
                            u(ilp1, jlp1, klp1) * onemfacx * onemfacy * onemfacz + &
                            u(il  , jlp1, klp1) *     facx * onemfacy * onemfacz

      this%fluidvel(ip,2) = v(il,   jl,   kl  ) *     facx *     facy *     facz + &
                            v(ilp1, jl,   kl  ) * onemfacx *     facy *     facz + &
                            v(ilp1, jlp1, kl  ) * onemfacx * onemfacy *     facz + &
                            v(il  , jlp1, kl  ) *     facx * onemfacy *     facz + &
                            v(il,   jl,   klp1) *     facx *     facy * onemfacz + &
                            v(ilp1, jl,   klp1) * onemfacx *     facy * onemfacz + &
                            v(ilp1, jlp1, klp1) * onemfacx * onemfacy * onemfacz + &
                            v(il  , jlp1, klp1) *     facx * onemfacy * onemfacz

      this%fluidvel(ip,3) = w(il,   jl,   kl  ) *     facx *     facy *     facz + &
                            w(ilp1, jl,   kl  ) * onemfacx *     facy *     facz + &
                            w(ilp1, jlp1, kl  ) * onemfacx * onemfacy *     facz + &
                            w(il  , jlp1, kl  ) *     facx * onemfacy *     facz + &
                            w(il,   jl,   klp1) *     facx *     facy * onemfacz + &
                            w(ilp1, jl,   klp1) * onemfacx *     facy * onemfacz + &
                            w(ilp1, jlp1, klp1) * onemfacx * onemfacy * onemfacz + &
                            w(il  , jlp1, klp1) *     facx * onemfacy * onemfacz

    enddo

end subroutine

subroutine get_interp_factors(this, xloc, dx, xline, nx, il, ilp1, facx, onemfacx)
    class(particles),           intent(inout) :: this
    real(rkind),                intent(in)    :: xloc, dx
    integer,                    intent(in)    :: nx
    real(rkind), dimension(nx), intent(in)    :: xline
    integer,                    intent(out)   :: il, ilp1
    real(rkind),                intent(out)   :: facx, onemfacx

    il = minloc(abs(xloc-xline), 1)
    if(xline(il) > xloc) il = il-1
    il = min(max(il, 1), nx-1)
    ilp1 = il + 1

    facx = (xline(ilp1) - xloc) / dx
    onemfacx = one - facx

end subroutine

subroutine update_acc(this, u, v, w)
    class(particles),              intent(inout) :: this
    real(rkind), dimension(:,:,:), intent(in)    :: u, v, w

    if(this%tracer_particles) then
        this%acc = zero
    elseif(this%inertial_particles) then
        call GracefulExit("Inertial particles have not been implemented yet",11)
    endif

end subroutine

!subroutine write_viz(this, tviz, outputfile)
!    class(particles),    intent(in) :: this
!    real(rkind),         intent(in) :: tviz
!    character(len=clen), intent(in) :: outputfile
!
!    integer :: iounit, ip
!
!    ioUnit = 10
!    open(unit=ioUnit, file=trim(outputfile), status='unknown', action='write')
!    write(iounit, '(e21.15)') tviz
!    do ip = 1, this%np
!        write(iounit, '(i5,1x, 12(e21.15, 1x))') ip, this%pos(ip,:), this%vel(ip,:), this%acc(ip,:), this%fluidvel(ip,:)
!    enddo
!    close(ioUnit)
!
!end subroutine

subroutine write_viz(this, istep, tviz, runID)
    class(particles),    intent(in) :: this
    integer,             intent(in) :: istep        ! timestep index
    real(rkind),         intent(in) :: tviz         ! time at this step
    character(len=*),    intent(in) :: runID        ! e.g. "Run34"

    integer :: iounit, ip
    character(len=clen) :: filename
        do ip = 1, this%np
            ! File name: Run34_PARTICLE_00001.out
            write(filename,'(A,"_PARTICLE_",I5.5,".out")') trim(runID), ip
        
            ioUnit = 10 + ip   ! careful: avoid too large ioUnit numbers!
            open(unit=ioUnit, file=trim(filename), status='unknown', &
                 action='write', position='append')
        
            write(ioUnit,'(I0,1x,E21.15,1x,12(E21.15,1x))') &
                 istep, tviz, this%pos(ip,:), this%vel(ip,:), this%acc(ip,:), &
        this%fluidvel(ip,:)
        
            close(ioUnit)
        end do

end subroutine

!subroutine initialize_particles(particle_locs, np)
!    use kind_parameters, only: rkind
!    real(rkind), dimension(:,:), allocatable, intent(out) :: particle_locs
!    integer, intent(in) :: np
!
!    integer :: Ny, Nz, i, j, p
!    real(8) :: y_start, y_end, z_start, z_end, dy, dz
!
!    ! z parameters
!    z_start = 0.5d0
!    z_end   = 6.0d0
!    dz      = 0.5d0
!    Nz      = int((z_end - z_start)/dz) + 1   ! number of z locations
!
!    ! particles along y per z-level
!    Ny = np / Nz
!    if (Ny < 1) then
!        call GracefulExit("np too small for chosen z-locations",11)
!    end if
!
!    ! y range
!    y_start = 3.0d0
!    y_end   = 9.0d0
!    dy      = (y_end - y_start) / dble(Ny-1)
!
!    ! Allocate: 3 coords × np
!    allocate(particle_locs(3, Ny*Nz))
!
!    ! Assign positions
!    p = 0
!    do j = 1, Nz
!        do i = 1, Ny
!            p = p + 1
!            particle_locs(1,p) = 6.0d0                          ! x fixed
!            particle_locs(2,p) = y_start + (i-1)*dy             ! y varying
!            particle_locs(3,p) = z_start + (j-1)*dz             ! z varying
!        end do
!    end do
!
!end subroutine

subroutine initialize_particles(particle_locs, np)
    use kind_parameters, only: rkind
    real(rkind), dimension(:,:), allocatable, intent(out) :: particle_locs
    integer, intent(in) :: np

    integer :: Ny, Nz, i, j, p
    real(8) :: x_center, y_center, z_center
    real(8) :: side, y_start, y_end, z_start, z_end, dy, dz

    ! turbine center
    x_center = 6.0d0
    y_center = 6.0d0
    z_center = 1.365d0

    ! downstream release plane (0.5 downstream of turbine)
    x_center = x_center + 0.5d0

    ! square side
    side = 1.5d0

    ! choose Ny and Nz ~ sqrt(np)
    Ny = int(sqrt(dble(np)))
    Nz = np / Ny
    if (Ny * Nz /= np) then
        call GracefulExit("np must be a perfect rectangle (Ny×Nz)", 11)
    end if

    ! y range
    y_start = y_center - side/2.0d0
    y_end   = y_center + side/2.0d0
    dy      = (y_end - y_start) / dble(Ny-1)

    ! z range
    z_start = z_center - side/2.0d0
    z_end   = z_center + side/2.0d0
    dz      = (z_end - z_start) / dble(Nz-1)

    ! Allocate: 3 coords × np
    allocate(particle_locs(3, np))

    ! Assign positions
    p = 0
    do j = 1, Nz
        do i = 1, Ny
            p = p + 1
            particle_locs(1,p) = x_center                ! fixed downstream plane
            particle_locs(2,p) = y_start + (i-1)*dy      ! y grid
            particle_locs(3,p) = z_start + (j-1)*dz      ! z grid
        end do
    end do

end subroutine


end module
