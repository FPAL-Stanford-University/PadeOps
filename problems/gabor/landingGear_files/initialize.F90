module HIT_Periodic_parameters
    use exits,            only: message
    use kind_parameters,  only: rkind
    use constants,        only: two, pi
    implicit none

    real(rkind) :: xmin = 0.d0, xmax = 1.d0
    real(rkind) :: ymin = 0.d0, ymax = 1.d0
    real(rkind) :: zmin = 0.d0, zmax = 1.d0

end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use HIT_Periodic_parameters    
    use kind_parameters,  only: rkind, clen
    use constants,        only: one
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k, ioUnit
    character(len=*),                intent(in)    :: inputfile
    integer :: ix1, ixn, iy1, iyn, iz1, izn, nxg, nyg, nzg
    namelist /landingGearINPUT/ xmin, xmax, ymin, ymax, zmin, zmax

    !Lx = two*pi; Ly = two*pi; Lz = one
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=landingGearINPUT)
    close(ioUnit)    

    !Lx = two*pi; Ly = two*pi; Lz = two*pi
    nxg = decomp%xsz(1); nyg = decomp%ysz(2); nzg = decomp%zsz(3)
    
    ! If base decomposition is in Y
    ix1 = decomp%xst(1); iy1 = decomp%xst(2); iz1 = decomp%xst(3)
    ixn = decomp%xen(1); iyn = decomp%xen(2); izn = decomp%xen(3)
    
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = (xmax - xmin)/real(nxg,rkind)
        dy = (ymax - ymin)/real(nyg,rkind)
        dz = (zmax - zmin)/real(nzg,rkind)

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = xmin + real( ix1 + i - 1, rkind ) * dx - dx/two
                    y(i,j,k) = ymin + real( iy1 + j - 1, rkind ) * dy - dy/two
                    z(i,j,k) = zmin + real( iz1 + k - 1, rkind ) * dz - dz/two
                end do
            end do
        end do
    end associate

end subroutine

subroutine initfields_wallM(decompC, decompE, inputfile, mesh, fieldsC, fieldsE)
    use HIT_Periodic_parameters
    use PadeDerOps, only: Pade6Stagg
    use kind_parameters,    only: rkind, clen 
    use constants,          only: zero, one, two, pi, half
    use gridtools,          only: alloc_buffs
    use decomp_2d          
    use decomp_2d_io
    use cd06staggstuff,     only: cd06stagg
    use exits,              only: gracefulExit,message_min_max
    use random

    implicit none
    type(decomp_info),               intent(in)    :: decompC
    type(decomp_info),               intent(in)    :: decompE
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    integer :: ioUnit
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, x, y, z

    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)

    u  = zero 
    v  = zero
    wC = zero 
    w  = zero 
      
    nullify(u,v,w,x,y,z)
   

end subroutine


subroutine set_planes_io(xplanes, yplanes, zplanes)
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 1, nyplanes = 1, nzplanes = 1

    allocate( xplanes(nxplanes))
    allocate( yplanes(nyplanes))
    allocate( zplanes(nzplanes))

    xplanes = [1]
    yplanes = [1]
    zplanes = [1]

end subroutine

subroutine set_KS_planes_io(planesCoarseGrid, planesFineGrid)
    integer, dimension(:), allocatable,  intent(inout) :: planesFineGrid
    integer, dimension(:), allocatable,  intent(inout) :: planesCoarseGrid
    
    allocate(planesCoarseGrid(1), planesFineGrid(1))
    planesCoarseGrid = [8]
    planesFineGrid = [16]

end subroutine

subroutine setInhomogeneousNeumannBC_Temp(inputfile, wTh_surf)
    use HIT_Periodic_parameters    
    use kind_parameters,    only: rkind, clen 
    use constants, only: one, zero
    implicit none
    real(rkind), intent(out) :: wTh_surf
    character(len=clen),                intent(in)    :: inputfile
    integer :: ioUnit 
    character(len=clen)  :: ufname, vfname, wfname 

    wTh_surf = zero 
    ! Do nothing really since this is an unstratified simulation
end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use HIT_Periodic_parameters    
    use kind_parameters,    only: rkind, clen 
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    real(rkind) :: ThetaRef
    integer :: iounit 
    character(len=clen)  :: ufname, vfname, wfname 
    
    Tsurf = zero; dTsurf_dt = zero; ThetaRef = one
    
    ! Do nothing really since this is an unstratified simulation
end subroutine


subroutine set_Reference_Temperature(inputfile, Tref)
    use HIT_Periodic_parameters    
    use kind_parameters,    only: rkind, clen 
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    integer :: iounit
     
    Tref = 0.d0
    
    ! Do nothing really since this is an unstratified simulation

end subroutine

subroutine hook_probes(inputfile, probe_locs)
    use kind_parameters,    only: rkind, clen 
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

subroutine hook_source(tsim,mesh,Re,urhs,vrhs,wrhs)
    use kind_parameters, only: rkind
    real(rkind),                     intent(in)    :: tsim, Re
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:),   intent(inout) :: urhs, vrhs, wrhs
    real(rkind) :: dx, t, R

    R = Re
    t = tsim
    dx = mesh(2,1,1,1) - mesh(1,1,1,1)

    urhs = urhs + 0.d0
    vrhs = vrhs + 0.d0
    wrhs = wrhs + 0.d0

end subroutine

subroutine getLargeScaleParams(KE,L,LES)
  use kind_parameters,    only: rkind
  use incompressibleGrid, only: igrid, readField3D
  use fortran_assert,     only: assert
  real(rkind), dimension(:,:,:), intent(out) :: KE, L
  type(igrid), intent(in) :: LES
  real(rkind), parameter :: kconst = 1.7d0
  integer, parameter :: model = 2

  select case (model)
  case (1)
    call readField3D(LES%runID, LES%step, LES%inputDir, "Diss", KE, LES%gpC) 
  
    ! Confirm dissipation is positive
    where (KE < 1.d-15) KE = 1.d-10
  
    ! Convert dissipation to KE
    KE = (KE**(2.d0/3.d0))*(L**(5.d0/3.d0))*1.03252d0*kconst
  case (2)
    call readField3D(LES%runID, LES%step, LES%inputDir, "tkeC", KE, LES%gpC) 
    
    ! Confirm KE is positive
    where (KE < 1.d-15) KE = 1.d-10
  end select

  call readField3D(LES%runID, LES%step, LES%inputDir, "Liso", L, LES%gpC) 
 
  
  where (L < 1.d-15) L = 1.d-10
  where(KE > 2.5d3) KE = 2.5d3

end subroutine

subroutine getDomainBoundaries(xDom,yDom,zDom,mesh)
  use kind_parameters, only: rkind
  use reductions,      only: p_minval, p_maxval
  real(rkind), dimension(2), intent(out) :: xDom, yDom, zDom
  real(rkind), dimension(:,:,:,:), intent(in) :: mesh
  integer :: nx, ny, nz
  real(rkind) :: dx, dy, dz
  
  nx = size(mesh,1)
  ny = size(mesh,2)
  nz = size(mesh,3)

  dx = mesh(2,1,1,1) - mesh(1,1,1,1)
  dy = mesh(1,2,1,2) - mesh(1,1,1,2)
  dz = mesh(1,1,2,3) - mesh(1,1,1,3)

  xDom(1) = p_minval(mesh(1 ,1 ,1 ,1)) - dx/2
  xDom(2) = p_maxval(mesh(nx,1 ,1 ,1)) + dx/2
  
  yDom(1) = p_minval(mesh(1 ,1 ,1 ,2)) - dy/2 
  yDom(2) = p_maxval(mesh(1 ,ny,1 ,2)) + dy/2
  
  zDom(1) = p_minval(mesh(1 ,1 ,1 ,3)) - dz/2 
  zDom(2) = p_maxval(mesh(1 ,1 ,nz,3)) + dz/2
end subroutine
