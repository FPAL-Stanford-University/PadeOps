module TaylorGreenPeriodic_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa 
    implicit none
    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.002_rkind ! 0.2% of the mean value
    integer :: nxg, nyg, nzg
   
    real(rkind), dimension(:,:,:), allocatable :: uexact, vexact, wexact, pexact

    real(rkind), parameter :: xdim = 1000._rkind, udim = 0.45_rkind
    real(rkind), parameter :: timeDim = xdim/udim
    integer :: direction
end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use TaylorGreenPeriodic_parameters    
    use kind_parameters,  only: rkind
    use constants,        only: one,two, pi
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k, ioUnit
    character(len=*),                intent(in)    :: inputfile
    integer :: ix1, ixn, iy1, iyn, iz1, izn, directionID
    real(rkind)  :: Lx = one, Ly = one, Lz = one
    namelist /TaylorGreenPeriodicINPUT/ directionID


    !Lx = two*pi; Ly = two*pi; Lz = one

    Lx = two*pi; Ly = two*pi; Lz = two*pi
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
    use TaylorGreenPeriodic_parameters
    use PadeDerOps, only: Pade6Stagg
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two, pi, half
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random
    use decomp_2d          
    use reductions,         only: p_maxval, p_minval
    use cd06staggstuff,     only: cd06stagg
    use exits,              only: gracefulExit,message_min_max
    implicit none
    type(decomp_info),               intent(in)    :: decompC
    type(decomp_info),               intent(in)    :: decompE
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    integer :: ioUnit
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, x, y, z
    real(rkind) :: dz
    real(rkind), dimension(:,:,:), allocatable :: xE, yE, zE, ybuffC, ybuffE, zbuffC, zbuffE
    type(cd06stagg), allocatable :: der
    integer :: directionID
    namelist /TaylorGreenPeriodicINPUT/ directionID

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=TaylorGreenPeriodicINPUT)
    close(ioUnit)    

    direction = directionID

    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)
 
    dz = z(1,1,2) - z(1,1,1)
    allocate(ybuffC(decompC%ysz(1),decompC%ysz(2), decompC%ysz(3)))
    allocate(ybuffE(decompE%ysz(1),decompE%ysz(2), decompE%ysz(3)))

    allocate(zbuffC(decompC%zsz(1),decompC%zsz(2), decompC%zsz(3)))
    allocate(zbuffE(decompE%zsz(1),decompE%zsz(2), decompE%zsz(3)))
    allocate(zE(decompE%xsz(1),decompE%xsz(2), decompE%xsz(3)))
    allocate(xE(decompE%xsz(1),decompE%xsz(2), decompE%xsz(3)))
    allocate(yE(decompE%xsz(1),decompE%xsz(2), decompE%xsz(3)))
    allocate(der)
    call der%init(decompC%zsz(3), dz, isTopEven = .false., isBotEven = .false., &
                             isTopSided = .true., isBotSided = .true.)
  

    call transpose_x_to_y(z,ybuffC,decompC)
    call transpose_y_to_z(ybuffC,zbuffC,decompC)
    call der%interpZ_C2E(zbuffC,zbuffE,size(zbuffC,1),size(zbuffC,2))                         
    call transpose_z_to_y(zbuffE,ybuffE,decompE)
    call transpose_y_to_x(ybuffE,zE,decompE) 
    
    call transpose_x_to_y(x,ybuffC,decompC)
    call transpose_y_to_z(ybuffC,zbuffC,decompC)
    call der%interpZ_C2E(zbuffC,zbuffE,size(zbuffC,1),size(zbuffC,2))                         
    call transpose_z_to_y(zbuffE,ybuffE,decompE)
    call transpose_y_to_x(ybuffE,xE,decompE) 

    call transpose_x_to_y(y,ybuffC,decompC)
    call transpose_y_to_z(ybuffC,zbuffC,decompC)
    call der%interpZ_C2E(zbuffC,zbuffE,size(zbuffC,1),size(zbuffC,2))                         
    call transpose_z_to_y(zbuffE,ybuffE,decompE)
    call transpose_y_to_x(ybuffE,yE,decompE) 
    
    select case(directionID)
    case(1) 
      u  =  sin(x)*cos(y)
      v  = -cos(x)*sin(y)
      w  = zero
      wC = zero 
    case(2)
      u  =  sin(x)*cos(z)
      v  =  zero
      w  = -cos(xE)*sin(zE)
      wC = -cos(x) *sin(z )
    case(3)
      u  =  zero
      v  =  sin(y)*cos(z)
      w  = -cos(yE)*sin(zE)
      wC = -cos(y )*sin(z )
    case default
      call gracefulExit("Invalid choice of directionID.",32)
    end select
  
    allocate(uexact(size(u,1),size(u,2), size(u,3)))
    allocate(vexact(size(v,1),size(v,2), size(v,3)))
    allocate(pexact(size(u,1),size(u,2), size(u,3)))
    allocate(wexact(size(wC,1),size(wC,2), size(wC,3)))


    call message_min_max(1,"Bounds for u:", p_minval(minval(u)), p_maxval(maxval(u)))
    call message_min_max(1,"Bounds for v:", p_minval(minval(v)), p_maxval(maxval(v)))
    call message_min_max(1,"Bounds for w:", p_minval(minval(wC)), p_maxval(maxval(wC)))
    
    deallocate(zE, xE, yE, ybuffC,ybuffE,zbuffC, zbuffE) 
    nullify(u,v,w,x,y,z)
   
    call message(0,"Velocity Field Initialized")

end subroutine


subroutine set_planes_io(xplanes, yplanes, zplanes)
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 1, nyplanes = 1, nzplanes = 1

    allocate(xplanes(nxplanes), yplanes(nyplanes), zplanes(nzplanes))

    xplanes = [64]
    yplanes = [64]
    zplanes = [20]

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
    use constants, only: one, zero 
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: wTh_surf
    real(rkind) :: ThetaRef
    integer :: iounit, directionID
    namelist /TaylorGreenPeriodicINPUT/ directionID
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=TaylorGreenPeriodicINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    real(rkind) :: ThetaRef
    integer :: iounit, directionID
    namelist /TaylorGreenPeriodicINPUT/ directionID
    
    Tsurf = zero; dTsurf_dt = zero; ThetaRef = one
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=TaylorGreenPeriodicINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine


subroutine set_Reference_Temperature(inputfile, Tref)
    use kind_parameters,    only: rkind
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    integer :: iounit, directionID
    
    namelist /TaylorGreenPeriodicINPUT/ directionID

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=TaylorGreenPeriodicINPUT)
    close(ioUnit)    
     
    Tref = 0.d0
    
    ! Do nothing really since this is an unstratified simulation

end subroutine

subroutine hook_probes(inputfile, probe_locs)
    use kind_parameters,    only: rkind
    real(rkind), dimension(:,:), allocatable, intent(inout) :: probe_locs
    character(len=*),                intent(in)    :: inputfile
    integer, parameter :: nprobes = 2
    
    ! IMPORTANT : Convention is to allocate probe_locs(3,nprobes)
    ! Example: If you have at least 3 probes:
    ! probe_locs(1,3) : x -location of the third probe
    ! probe_locs(2,3) : y -location of the third probe
    ! probe_locs(3,3) : z -location of the third probe


    ! Add probes here if needed
    ! Example code: The following allocates 2 probes at (0.1,0.1,0.1) and
    ! (0.2,0.2,0.2)  
    allocate(probe_locs(3,nprobes))
    probe_locs(1,1) = 0.1d0; probe_locs(2,1) = 0.1d0; probe_locs(3,1) = 0.1d0;
    probe_locs(1,2) = 0.2d0; probe_locs(2,2) = 0.2d0; probe_locs(3,2) = 0.2d0;


end subroutine

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

subroutine hook_source(gp,tsim,mesh,Re,urhs,vrhs,wrhs)
  use kind_parameters, only: rkind
  use decomp_2d,       only: decomp_info
  use constants,       only: zero
  type(decomp_info),               intent(in)         :: gp
  real(rkind), intent(in)                             :: tsim, Re
  real(rkind), dimension(:,:,:,:), intent(in), target :: mesh
  real(rkind), dimension(:,:,:),   intent(inout)      :: urhs, vrhs, wrhs
  real(rkind), dimension(:,:,:), pointer              :: x, y, z

end subroutine
