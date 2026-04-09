module concurrentSimulation_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa 
    implicit none
    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.002_rkind ! 0.2% of the mean value
    integer :: nxg, nyg, nzg
    logical :: isPrecursor = .true.
    
    real(rkind), parameter :: xdim = 1000._rkind, udim = 0.45_rkind
    real(rkind), parameter :: timeDim = xdim/udim

end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use concurrentSimulation_parameters    
    use kind_parameters,  only: rkind
    use constants,        only: one,two
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    character(len=*),                intent(in)    :: inputfile
    integer :: i,j,k, ioUnit
    integer :: ix1, ixn, iy1, iyn, iz1, izn
    real(rkind)  :: Lx = one, Ly = one, Lz = one, z0init = 2.0d-4, ustarinit = 1.0d0
    namelist /concurrentSimulationINPUT/ Lx, Ly, Lz, z0init, ustarinit

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrentSimulationINPUT)
    close(ioUnit)    

    !Lx = two*pi; Ly = two*pi; Lz = one

    nxg = decomp%xsz(1); nyg = decomp%ysz(2); nzg = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%xst(1); iy1 = decomp%xst(2); iz1 = decomp%xst(3)
    ixn = decomp%xen(1); iyn = decomp%xen(2); izn = decomp%xen(3)
    
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = Lx/real(nxg,rkind)
        dy = Ly/real(nyg,rkind)
        dz = Lz/real(nzg,rkind)

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                    y(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 + k - 1, rkind ) * dz + dz/two
                end do
            end do
        end do

        ! Shift everything to the origin 
        x = x - dx
        y = y - dy
        z = z - dz 

    end associate

end subroutine

subroutine initfields_wallM(decompC, decompE, inputfile, mesh, fieldsC, fieldsE)
    use concurrentSimulation_parameters
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two, pi, half
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random
    use decomp_2d          
    use reductions,         only: p_maxval, p_minval
    use exits,              only: message_min_max
    implicit none
    type(decomp_info),               intent(in)    :: decompC
    type(decomp_info),               intent(in)    :: decompE
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    integer :: ioUnit
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, x, y, z
    real(rkind), dimension(:,:,:), allocatable :: randArr, ybuffC, ybuffE, zbuffC, zbuffE
    integer :: nz, nzE, k
    real(rkind)  :: Lx = one, Ly = one, Lz = one
    real(rkind) ::  z0init = 2.0d-4, epsnd, yperiods = 3.0d0, zpeak = 0.2d0, xperiods = 3.0d0, ustarinit = 1.0d0
    namelist /concurrentSimulationINPUT/ Lx, Ly, Lz, z0init, zpeak, ustarinit

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrentSimulationINPUT)
    close(ioUnit)    


    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)
 
    if(isPrecursor) then
       epsnd = 5.0d0 
       !epsnd = 0.0_rkind 
       u = (ustarinit/kappa)*log(z/z0init) + epsnd*cos(yperiods*two*pi*y/Ly)*exp(-half*(z/zpeak/Lz)**2)
       v = epsnd*(z/Lz)*cos(xperiods*two*pi*x/Lx)*exp(-half*(z/zpeak/Lz)**2)
       wC= zero  
    else
       epsnd = zero
       u = (ustarinit/kappa)*log(z/z0init) + epsnd*cos(yperiods*two*pi*y/Ly)*exp(-half*(z/zpeak/Lz)**2)
       v = epsnd*(z/Lz)*cos(xperiods*two*pi*x/Lx)*exp(-half*(z/zpeak/Lz)**2)
       wC= zero  
    endif
    isPrecursor = .false.

    !allocate(randArr(size(u,1),size(u,2),size(u,3)))
    !call gaussian_random(randArr,-one,one,seedu + 10*nrank)
    !do k = 1,size(randArr,3)
    !     u(:,:,k) = u(:,:,k) + 0.01*randArr(:,:,k)
    !end do
    !deallocate(randArr)

    !seedu = seedu + 100000

    call message_min_max(1,"Bounds for u:", p_minval(minval(u)), p_maxval(maxval(u)))
    call message_min_max(1,"Bounds for v:", p_minval(minval(v)), p_maxval(maxval(v)))
    call message_min_max(1,"Bounds for w:", p_minval(minval(w)), p_maxval(maxval(w)))
    
    !u = one!1.6d0*z*(2.d0 - z) 
    !v = zero;
    !w = zero;

    ! Interpolate wC to w
    allocate(ybuffC(decompC%ysz(1),decompC%ysz(2), decompC%ysz(3)))
    allocate(ybuffE(decompE%ysz(1),decompE%ysz(2), decompE%ysz(3)))

    allocate(zbuffC(decompC%zsz(1),decompC%zsz(2), decompC%zsz(3)))
    allocate(zbuffE(decompE%zsz(1),decompE%zsz(2), decompE%zsz(3)))
   
    nz = decompC%zsz(3)
    nzE = nz + 1

    call transpose_x_to_y(wC,ybuffC,decompC)
    call transpose_y_to_z(ybuffC,zbuffC,decompC)
    zbuffE = zero
    zbuffE(:,:,2:nzE-1) = half*(zbuffC(:,:,1:nz-1) + zbuffC(:,:,2:nz))
    call transpose_z_to_y(zbuffE,ybuffE,decompE)
    call transpose_y_to_x(ybuffE,w,decompE) 
    
    

    deallocate(ybuffC,ybuffE,zbuffC, zbuffE) 
  
      
    nullify(u,v,w,x,y,z)
   

    call message(0,"Velocity Field Initialized")

end subroutine


subroutine set_planes_io(xplanes, yplanes, zplanes)
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 3, nyplanes = 3, nzplanes = 3

    allocate(xplanes(nxplanes))
    allocate(yplanes(nyplanes))
    allocate(zplanes(nzplanes))

    xplanes = [32, 64, 128]
    yplanes = [15, 30, 45]
    zplanes = [15, 30, 45]

end subroutine

subroutine set_KS_planes_io(planesCoarseGrid, planesFineGrid)
    integer, dimension(:), allocatable,  intent(inout) :: planesFineGrid
    integer, dimension(:), allocatable,  intent(inout) :: planesCoarseGrid
    
    allocate(planesCoarseGrid(1), planesFineGrid(1))
    planesCoarseGrid = [8]
    planesFineGrid = [16]

end subroutine

subroutine setInhomogeneousNeumannBC_Temp(inputfile, wTh_surf)
    use kind_parameters,    only: rkind
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: wTh_surf
    real(rkind) :: ThetaRef, Lx, Ly, Lz, z0init
    integer :: iounit
    namelist /concurrentSimulationINPUT/ Lx, Ly, Lz, z0init
    
    wTh_surf = zero
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrentSimulationINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    real(rkind) :: ThetaRef, Lx, Ly, Lz, z0init
    integer :: iounit
    namelist /concurrentSimulationINPUT/ Lx, Ly, Lz, z0init
    
    Tsurf = zero; dTsurf_dt = zero; ThetaRef = one
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrentSimulationINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine


subroutine set_Reference_Temperature(inputfile, Tref)
    use kind_parameters,    only: rkind
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    real(rkind) :: Lx, Ly, Lz, z0init
    integer :: iounit
    
    namelist /concurrentSimulationINPUT/ Lx, Ly, Lz, z0init

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=concurrentSimulationINPUT)
    close(ioUnit)    
     
    Tref = 0.d0
    
    ! Do nothing really since this is an unstratified simulation

end subroutine

subroutine hook_probes(inputfile, probe_locs)
    use kind_parameters,    only: rkind
    real(rkind), dimension(:,:), allocatable, intent(inout) :: probe_locs
    character(len=*),                intent(in)    :: inputfile
    integer, parameter :: nprobes = 1 
    
    ! IMPORTANT : Convention is to allocate probe_locs(3,nprobes)
    ! Example: If you have at least 3 probes:
    ! probe_locs(1,3) : x -location of the third probe
    ! probe_locs(2,3) : y -location of the third probe
    ! probe_locs(3,3) : z -location of the third probe


    ! Add probes here if needed
    ! Example code: The following allocates 2 probes at (0.1,0.1,0.1) and
    ! (0.2,0.2,0.2)  
    print*, inputfile
    allocate(probe_locs(3,nprobes))
    probe_locs(1,1) = 5.0d0; probe_locs(2,1) = 6.0d0; probe_locs(3,1) = 1.365d0;
end subroutine

!subroutine hook_probes(inputfile, probe_locs)
!    use kind_parameters, only: rkind 
!    real(rkind), dimension(:,:), allocatable, intent(inout) :: probe_locs
!    character(len=*), intent(in) :: inputfile
!
!    integer :: nprobes, i
!    real(8) :: y_start, y_end, dy
!
!    ! Probe configuration
!    y_start = 0.5d0
!    y_end   = 7.0d0
!    dy      = 0.5d0
!
!    nprobes = int((y_end - y_start)/dy) + 1
!
!    ! Allocate and assign
!    allocate(probe_locs(3,nprobes))
!    do i = 1, nprobes
!        probe_locs(1,i) = 6.0d0            ! x
!        probe_locs(2,i) = 6.0d0  ! y
!        probe_locs(3,i) = y_start + (i-1)*dy          ! z
!    end do
!
!end subroutine


!subroutine hook_probes(inputfile, probe_locs)
!    use kind_parameters, only: rkind
!    real(rkind), dimension(:,:), allocatable, intent(out) :: probe_locs
!    character(len=*), intent(in) :: inputfile
!
!    integer :: nprobes, Ny, Nz, i, j, p
!    real(8) :: x_center, y_center, z_center
!    real(8) :: side, y_start, y_end, z_start, z_end, dy, dz
!
!    ! no of probes 
!    nprobes = 6400  
!
!    ! turbine center
!    x_center = 6.0d0
!    y_center = 6.0d0
!    z_center = 1.365d0
!
!    ! downstream release plane (0.5 downstream of turbine)
!    x_center = x_center + 2.5d0
!
!    ! square side
!    side = 1.5d0
!
!    ! grid dimensions (Ny x Nz)
!    Ny = int(sqrt(dble(nprobes)))
!    Nz = nprobes / Ny
!    if (Ny * Nz /= nprobes) then
!        call GracefulExit("nprobes must be a rectangular grid (Ny×Nz)", 12)
!    end if
!
!    ! y range
!    y_start = y_center - side/2.0d0
!    y_end   = y_center + side/2.0d0
!    dy      = (y_end - y_start) / dble(Ny-1)
!
!    ! z range
!    z_start = z_center - side/2.0d0
!    z_end   = z_center + side/2.0d0
!    dz      = (z_end - z_start) / dble(Nz-1)
!
!    ! Allocate: 3 coords × nprobes
!    allocate(probe_locs(3, nprobes))
!
!    ! Assign probe locations (grid in y–z plane)
!    p = 0
!    do j = 1, Nz
!        do i = 1, Ny
!            p = p + 1
!            probe_locs(1,p) = x_center                ! x fixed
!            probe_locs(2,p) = y_start + (i-1)*dy      ! y grid
!            probe_locs(3,p) = z_start + (j-1)*dz      ! z grid
!        end do
!    end do
!
!end subroutine



!subroutine hook_probes(inputfile, probe_locs)
!    use kind_parameters, only: rkind
!    real(rkind), dimension(:,:), allocatable, intent(inout) :: probe_locs
!    character(len=*), intent(in) :: inputfile
!
!    integer :: nbase, ny, nprobes, iprobe, ixz, iy
!    real(rkind), dimension(:,:), allocatable :: base_locs
!    real(rkind), dimension(:), allocatable   :: ylocs
!
!    ! base (x,z) locations (with dummy y=0.0, will be replaced by ylocs)
!    nbase = 28
!    allocate(base_locs(3,nbase))
!    base_locs(:,1)  = [1.0d0, 0.0d0, 0.1d0]
!    base_locs(:,2)  = [1.0d0, 0.0d0, 0.2d0]
!    base_locs(:,3)  = [1.0d0, 0.0d0, 0.5d0]
!    base_locs(:,4)  = [1.0d0, 0.0d0, 0.8d0]
!    base_locs(:,5)  = [2.5d0, 0.0d0, 0.1d0]
!    base_locs(:,6)  = [2.5d0, 0.0d0, 0.2d0]
!    base_locs(:,7)  = [2.5d0, 0.0d0, 0.5d0]
!    base_locs(:,8)  = [2.5d0, 0.0d0, 0.8d0]
!    base_locs(:,9)  = [3.0d0, 0.0d0, 0.1d0]
!    base_locs(:,10) = [3.0d0, 0.0d0, 0.2d0]
!    base_locs(:,11) = [3.0d0, 0.0d0, 0.5d0]
!    base_locs(:,12) = [3.0d0, 0.0d0, 0.8d0]
!    base_locs(:,13) = [3.5d0, 0.0d0, 0.1d0]
!    base_locs(:,14) = [3.5d0, 0.0d0, 0.2d0]
!    base_locs(:,15) = [3.5d0, 0.0d0, 0.5d0]
!    base_locs(:,16) = [3.5d0, 0.0d0, 0.8d0]
!    base_locs(:,17) = [4.0d0, 0.0d0, 0.1d0]
!    base_locs(:,18) = [4.0d0, 0.0d0, 0.2d0]
!    base_locs(:,19) = [4.0d0, 0.0d0, 0.5d0]
!    base_locs(:,20) = [4.0d0, 0.0d0, 0.8d0]
!    base_locs(:,21) = [7.0d0, 0.0d0, 0.1d0]
!    base_locs(:,22) = [7.0d0, 0.0d0, 0.2d0]
!    base_locs(:,23) = [7.0d0, 0.0d0, 0.5d0]
!    base_locs(:,24) = [7.0d0, 0.0d0, 0.8d0]
!    base_locs(:,25) = [12.0d0, 0.0d0, 0.1d0]
!    base_locs(:,26) = [12.0d0, 0.0d0, 0.2d0]
!    base_locs(:,27) = [12.0d0, 0.0d0, 0.5d0]
!    base_locs(:,28) = [12.0d0, 0.0d0, 0.8d0]
!
!    ! spanwise locations
!    allocate(ylocs(17))
!    ylocs = (/ 0.4d0, 0.45d0, 0.5d0, 0.55d0, 0.6d0, 0.65d0, 0.7d0, 0.75d0, 0.8d0, 0.85d0 ,0.9d0, 0.95d0, 1.0d0, 1.05d0, 1.1d0, 1.15d0, 1.2d0 /)
!
!    ny = size(ylocs)
!    nprobes = nbase * ny
!    allocate(probe_locs(3,nprobes))
!
!    ! fill probes
!    iprobe = 0
!    do iy = 1, ny
!        do ixz = 1, nbase
!            iprobe = iprobe + 1
!            probe_locs(1,iprobe) = base_locs(1,ixz)
!            probe_locs(2,iprobe) = ylocs(iy)
!            probe_locs(3,iprobe) = base_locs(3,ixz)
!        end do
!    end do
!
!!    print *, "Total probes placed = ", nprobes
!end subroutine
!

subroutine initScalar(decompC, inpDirectory, mesh, scalar_id, scalarField)
    use kind_parameters, only: rkind
    use decomp_2d,        only: decomp_info
    type(decomp_info),                                          intent(in)    :: decompC
    character(len=*),                intent(in)    :: inpDirectory
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    integer, intent(in)                            :: scalar_id
    real(rkind), dimension(:,:,:), intent(out)     :: scalarField

    scalarField = 0.d0
end subroutine 

subroutine setScalar_source(decompC, inpDirectory, mesh, scalar_id, scalarSource)
    use kind_parameters, only: rkind
    use decomp_2d,        only: decomp_info
    type(decomp_info),                                          intent(in)    :: decompC
    character(len=*),                intent(in)    :: inpDirectory
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    integer, intent(in)                            :: scalar_id
    real(rkind), dimension(:,:,:), intent(out)     :: scalarSource

    scalarSource = 0.d0
end subroutine 
