module AD_Coriolis_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa, pi
    implicit none
    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.002_rkind ! 0.2% of the mean value
    integer :: nxg, nyg, nzg
    
    real(rkind), parameter :: xdim = 1000._rkind, udim = 0.45_rkind
    real(rkind), parameter :: timeDim = xdim/udim
    real(rkind), dimension(:,:,:), allocatable :: utarget, vtarget, wtarget

    contains
subroutine init_fringe_targets(inputfile, mesh)
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two, pi, half
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random
    use decomp_2d          
    use reductions,         only: p_maxval, p_minval
    use exits,              only: message_min_max
    implicit none
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:), pointer :: z
    real(rkind) :: Lx, Ly, Lz, uInflow, vInflow, inflowAngle, alphaShear  
    real(rkind) :: InflowProfileAmplit, InflowProfileThick, zMid
    integer :: i,j,k, ioUnit
    integer :: InflowProfileType

    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, & 
                                InflowProfileAmplit, InflowProfileThick, InflowProfileType, inflowAngle, alphaShear

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=AD_CoriolisINPUT)
    close(ioUnit)    

    ! Initialize the velocity targets
    ! Put the same velocity profile in init fields and target
    ! Do something similar for v
    ! Compute the targets
    wtarget = zero 
    zMid = Lz / two
    z => mesh(:,:,:,3)
    !call get_u(uInflow, vInflow, InflowProfileAmplit, InflowProfileThick, mesh(:,:,:,3), zMid, InflowProfileType, utarget, vtarget)    
    select case(InflowProfileType)
      case(0)
          utarget = uInflow 
          vtarget = zero
      case(1)
          utarget = uInflow*(one  + InflowProfileAmplit*tanh((z-zMid)/InflowProfileThick))
          vtarget = zero
      case(2)
          utarget = uInflow*(one  + InflowProfileAmplit*tanh((z-zMid)/InflowProfileThick))
          vtarget = vInflow * tanh((z-zMid)/InflowProfileThick);
      case(3)
          utarget = uInflow*cos(inflowAngle*pi/180.d0) 
          vtarget = uInflow*sin(inflowAngle*pi/180.d0)
      case(4)
          utarget = uInflow*cos(inflowAngle*pi/180.d0) * (z/zMid) ** alphaShear
          vtarget = uInflow*sin(inflowAngle*pi/180.d0) * (z/zMid) ** alphaShear
    end select


    ! The velocity profile in z needs to go to slip wall at the top
    ! Both u and v need slip conditions

end subroutine

end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use AD_Coriolis_parameters    
    use kind_parameters,  only: rkind
    use constants,        only: one,two
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k, ioUnit
    character(len=*),                intent(in)    :: inputfile
    integer :: ix1, ixn, iy1, iyn, iz1, izn
    real(rkind)  :: Lx = one, Ly = one, Lz = one, G_alpha
    real(rkind) :: uInflow, vInflow, inflowAngle, alphaShear  
    real(rkind) :: InflowProfileAmplit, InflowProfileThick, zMid
    integer :: InflowProfileType
    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, & 
                                InflowProfileAmplit, InflowProfileThick, InflowProfileType, inflowAngle, alphaShear

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=AD_CoriolisINPUT)
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
    use AD_Coriolis_parameters
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
    integer :: nz, nzE
    real(rkind)  :: Lx = one, Ly = one, Lz = one, G_alpha
    real(rkind) :: uInflow, vInflow, inflowAngle, alphaShear
    real(rkind) :: InflowProfileAmplit, InflowProfileThick, zMid
    integer :: InflowProfileType
    
    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, & 
                                InflowProfileAmplit, InflowProfileThick, InflowProfileType, inflowAngle, alphaShear

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=AD_CoriolisINPUT)
    close(ioUnit)    


    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)

   
    !u = one * cos(G_alpha * pi / 180.d0)!z*(2 - z) + epsnd*(z/Lz)*cos(periods*2*pi*x/Lx)*sin(periods*2*pi*y/Lx)*exp(-0.5*(z/zpeak/Lz)**2) &
      !+ epsnd*((2-z)/Lz)*cos(periods*2*pi*x/Lx)*sin(periods*2*pi*y/Lx)*exp(-0.5*((2-z)/zpeak/Lz)**2)
    
    !v = one * sin(G_alpha * pi / 180.d0)!- epsnd*(z/Lz)*sin(periods*2*pi*x/Lx)*cos(periods*2*pi*y/Lx)*exp(-0.5*(z/zpeak/Lz)**2) &
        !- epsnd*((2-z)/Lz)*sin(periods*2*pi*x/Lx)*cos(periods*2*pi*y/Lx)*exp(-0.5*((2-z)/zpeak/Lz)**2)
    
    wC = zero
    zMid = Lz / 2.d0
    select case(InflowProfileType)
      case(0)
          u = uInflow 
          v = zero
      case(1)
          u = uInflow*(one  + InflowProfileAmplit*tanh((z-zMid)/InflowProfileThick))
          v = zero
      case(2)
          u = uInflow*(one  + InflowProfileAmplit*tanh((z-zMid)/InflowProfileThick))
          v = vInflow * tanh((z-zMid)/InflowProfileThick);
      case(3)
          u = uInflow*cos(inflowAngle*pi/180.d0)
          v = uInflow*sin(inflowAngle*pi/180.d0)
      case(4)
          u = uInflow*cos(inflowAngle*pi/180.d0) * (z/zMid) ** alphaShear
          v = uInflow*sin(inflowAngle*pi/180.d0) * (z/zMid) ** alphaShear
    end select

    !call get_u(uInflow, vInflow, InflowProfileAmplit, InflowProfileThick, z, zMid, InflowProfileType, u, v)    
    
    !allocate(randArr(size(u,1),size(u,2),size(u,3)))
    !call gaussian_random(randArr,-one,one,seedu + 10*nrank)
    !!do k = 1,size(randArr,3)
    !!     u(:,:,k) = u(:,:,k) + randscale*randArr(:,:,k)
    !!end do
    !deallocate(randArr)

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
    integer, parameter :: nxplanes = 3, nyplanes = 1, nzplanes = 1

    allocate(xplanes(nxplanes), yplanes(nyplanes), zplanes(nzplanes))

    xplanes = [128, 256, 384]
    yplanes = [256]
    zplanes = [96]

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
    integer :: ioUnit 
    real(rkind) :: ThetaRef, Lx, Ly, Lz
    logical :: initPurturbations = .false. 
    real(rkind) :: uInflow, vInflow, inflowAngle, alphaShear  
    real(rkind) :: InflowProfileAmplit, InflowProfileThick, zMid
    integer :: InflowProfileType
    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, & 
                                InflowProfileAmplit, InflowProfileThick, InflowProfileType, inflowAngle, alphaShear
     
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=AD_CoriolisINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    real(rkind) :: ThetaRef, Lx, Ly, Lz, G_alpha
    integer :: iounit
    real(rkind) :: uInflow, vInflow, inflowAngle, alphaShear  
    real(rkind) :: InflowProfileAmplit, InflowProfileThick, zMid
    integer :: InflowProfileType
    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, & 
                                InflowProfileAmplit, InflowProfileThick, InflowProfileType, inflowAngle, alphaShear
    
    Tsurf = zero; dTsurf_dt = zero; ThetaRef = one
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=AD_CoriolisINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine


subroutine set_Reference_Temperature(inputfile, Tref)
    use kind_parameters,    only: rkind
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    real(rkind) :: Lx, Ly, Lz, G_alpha
    integer :: iounit
    real(rkind) :: uInflow, vInflow, inflowAngle, alphaShear  
    real(rkind) :: InflowProfileAmplit, InflowProfileThick, zMid
    integer :: InflowProfileType
    
    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow, vInflow, & 
                                InflowProfileAmplit, InflowProfileThick, InflowProfileType, inflowAngle, alphaShear

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=AD_CoriolisINPUT)
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
    print*, inputfile
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

subroutine get_u(uInflow, vInflow, InflowProfileAmplit, InflowProfileThick, z, zMid, InflowProfileType, inflowAngle, alphaShear, u, v)
    use kind_parameters, only: rkind
    use constants,          only: zero, one, two, pi, half
    implicit none
    real(rkind), dimension(:,:,:), intent(inout) :: u, v
    real(rkind), dimension(:,:,:), intent(in) :: z
    real(rkind), intent(in) :: InflowProfileAmplit, InflowProfileThick, zMid, uInflow, vInflow, inflowAngle, alphaShear
    integer, intent(in) :: InflowProfileType

    select case(InflowProfileType)
      case(0)
          u = uInflow 
          v = zero
      case(1)
          u = uInflow*(one  + InflowProfileAmplit*tanh((z-zMid)/InflowProfileThick))
          v = zero
      case(2)
          u = uInflow*(one  + InflowProfileAmplit*tanh((z-zMid)/InflowProfileThick))
          v = vInflow * tanh((z-zMid)/InflowProfileThick);
      case(3)
          u = uInflow*cos(inflowAngle*pi/180.d0) 
          v = uInflow*sin(inflowAngle*pi/180.d0)
      case(4)
          u = uInflow*cos(inflowAngle*pi/180.d0) * (z/zMid) ** alphaShear
          v = uInflow*sin(inflowAngle*pi/180.d0) * (z/zMid) ** alphaShear
    end select

end subroutine


