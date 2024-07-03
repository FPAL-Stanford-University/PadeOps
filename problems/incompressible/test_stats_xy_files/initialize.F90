module test_stats_xy_parameters

    use exits, only: message
    use kind_parameters,  only: rkind, clen
    use constants, only: kappa, zero
    use fortran_assert,     only: assert
    implicit none
    integer :: nxSize = 256, nySize = 256, nzSize = 256
    real(rkind), dimension(3) :: wavenums
    real(rkind), dimension(8,3) :: coefs
    real(rkind), dimension(9) :: tauCoefs

contains

  pure subroutine Sfunc(x, val)
     real(rkind), dimension(:,:,:), intent(in) :: x
     real(rkind), dimension(:,:,:), intent(out) :: val
  
     val = 0.d0
     where (x>0.d0) 
        val = 1.d0/(1.d0 + exp(min(1.d0/(x - 1.d0 + 1.d-18) + 1.d0/(x + 1.d-18),50.d0)))
     end where
  
     where (x>1.d0) 
        val = 1.d0 
     end where
  
  end subroutine
  
end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use test_stats_xy_parameters    
    use kind_parameters,  only: rkind
    use constants,        only: one,two, pi
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k, ioUnit
    character(len=*),                intent(in)    :: inputfile
    integer :: nxg, nyg, nzg
    integer :: ix1, ixn, iy1, iyn, iz1, izn
    real(rkind)  :: Lx = one, Ly = one, Lz = one
    logical :: symmetricDomain = .true.
    real(rkind) :: zmin = -one, Tref
    character(len=clen) :: stats_info_dir
    integer :: num_stats_instances = 1

    namelist /SMinput/ Lx, Ly, Lz, symmetricDomain, zmin, Tref, stats_info_dir, num_stats_instances

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=SMinput)
    close(ioUnit)    

    nxg = decomp%xsz(1); nyg = decomp%ysz(2); nzg = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%xst(1); iy1 = decomp%xst(2); iz1 = decomp%xst(3)
    ixn = decomp%xen(1); iyn = decomp%xen(2); izn = decomp%xen(3)
    
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = Lx/real(nxg,rkind)
        dy = Ly/real(nyg,rkind)
        if (symmetricDomain) then
          ! z goes from -Lz to Lz
          dz = 2.d0*Lz/real(nzg,rkind)
          
          do k=1,size(mesh,3)
              do j=1,size(mesh,2)
                  do i=1,size(mesh,1)
                      x(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                      y(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                      z(i,j,k) = real( -nzg/2, rkind ) * dz + &
                        & real( iz1 + k - 1, rkind ) * dz + dz/two
                  end do
              end do
          end do
        else
          ! z goes from zmin to Lz where 
          dz = (Lz - zmin)/real(nzg,rkind)
          
          do k=1,size(mesh,3)
              do j=1,size(mesh,2)
                  do i=1,size(mesh,1)
                      x(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                      y(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                      z(i,j,k) = zmin + real( iz1 + k - 1, rkind ) * dz + dz/two
                  end do
              end do
          end do
        end if

        ! Shift everything to the origin 
        x = x - dx
        y = y - dy
        z = z - dz 

    end associate
   
    nxSize = nxg; nySize = nyg; nzSize = nzg

end subroutine

subroutine initfields_wallM(decompC, decompE, inputfile, meshC, meshE, fieldsC, fieldsE)
    use test_stats_xy_parameters
    use kind_parameters,    only: rkind, clen
    use constants,          only: zero, one, two, pi, half
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random
    use decomp_2d          
    use reductions,         only: p_maxval, p_minval
    use exits,              only: message_min_max, gracefulExit, message
    use cd06staggstuff,     only: cd06stagg
    use basic_io,           only: read_1d_ascii
    implicit none
    type(decomp_info),               intent(in)    :: decompC
    type(decomp_info),               intent(in)    :: decompE
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: meshC, meshE
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    real(rkind), dimension(:,:,:), pointer :: u => null(), v => null(), w => null(), &
      wC => null(), x => null(), y => null(), z => null(), T => null()
    real(rkind), dimension(:,:,:), pointer :: xE, yE, zE
    real(rkind), dimension(:,:,:), allocatable :: ybuffC, ybuffE, zbuffC, zbuffE
    real(rkind)  :: Lx = one, Ly = one, Lz = one
!    real(rkind)  :: zTop_cell, zBot_cell, zMid
    type(cd06stagg), allocatable :: der
    real(rkind) :: dz, zmin, Tref
    integer :: ioUnit, ierr
    integer :: i, j
    logical :: symmetricDomain
    character(len=clen) :: fname
    real(rkind), dimension(:), allocatable :: TinZ

    ! INPUT namelist variables (only inputdir is needed)
    integer :: nx, ny, nz, nsteps, nstepConstDt, nxS, nyS, nzS
    real(rkind) :: tstop, dt, CFL, CviscDT
    character(len=clen) :: inputdir, outputdir
    integer :: prow, pcol, restartFile_TID, restartFile_RID
    logical :: useRestartFile, restartFromDifferentGrid

    ! PHYSICS namelist variables (only isStratified is needed)
    logical :: isInviscid, useCoriolis, useExtraForcing, useMoisture, useSGS, &
      useforcedStratification, useGeostrophicForcing, assume_fplane, &
      useHITForcing,useScalars,useHITRealSpaceLinearForcing,addExtraSourceTerm, &
      useImmersedBodies
    real(rkind) :: Re, Ro, Pr, Fr, Ra, PrandtlFluid, BulkRichardson, &
      G_geostrophic,G_alpha,dpFdx,dpFdy,dpFdz,latitude,frameAngle, &
      HITForceTimeScale, immersed_taufact
    integer :: BuoyancyTermType,buoyancyDirection,numberOfImmersedBodies, localizedForceLayer
    logical :: isStratified = .false., removeMean = .false.
    character(len=clen) :: stats_info_dir
    integer :: num_stats_instances = 1

    namelist /SMinput/ Lx, Ly, Lz, symmetricDomain, zmin, Tref, stats_info_dir, num_stats_instances
    namelist /INPUT/ nx, ny, nz, tstop, dt, CFL, nsteps, inputdir, outputdir, prow, pcol, &
                    useRestartFile, restartFile_TID, restartFile_RID, CviscDT, &
                    nstepConstDt, restartFromDifferentGrid, nxS, nyS, nzS
    namelist /PHYSICS/isInviscid,useCoriolis,useExtraForcing,isStratified,&
      useMoisture,Re,Ro,Pr,Fr, Ra, useSGS, PrandtlFluid, BulkRichardson, &
      BuoyancyTermType,useforcedStratification, useGeostrophicForcing, &
      G_geostrophic, G_alpha, dpFdx,dpFdy,dpFdz,assume_fplane,latitude,&
      useHITForcing, useScalars, frameAngle, buoyancyDirection, &
      useHITRealSpaceLinearForcing, HITForceTimeScale, addExtraSourceTerm, &
      useImmersedBodies, numberOfImmersedBodies, immersed_taufact, &
      localizedForceLayer, removeMean

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=INPUT)
    read(unit=ioUnit, NML=PHYSICS)
    read(unit=ioUnit, NML=SMinput)
    close(ioUnit)    
   
    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1)
    T => fieldsC(:,:,:,7) 
      
    z => meshC(:,:,:,3)
    y => meshC(:,:,:,2)
    x => meshC(:,:,:,1)

    zE => meshE(:,:,:,3)
    yE => meshE(:,:,:,2)
    xE => meshE(:,:,:,1)

    ! Parameters defining the fields
    wavenums = [2.d0*pi/Lx, 6.d0*pi/Lx, 6.d0*pi/Lx]
    coefs(1,:) = [0.428751518293931d0,   0.04d0,   0.10d0]
    coefs(2,:) = [0.006509245148313d0,   0.02d0,   0.05d0]
    coefs(3,:) = [0.011821728885865d0,   0.02d0,   0.05d0]
    coefs(4,:) = [0.013463537687407d0,   0.017667033524145d0,   0.019068848073671d0]
    coefs(5,:) = [0.079220732955955d0,   0.095949242639290d0,   0.065574069915659d0]
    coefs(6,:) = [0.d0, 0.039222701953417d0,   0.065547789017756d0]
    coefs(7,:) = [0.d0, 0.017118668781156d0,   0.070604608801961d0]
    coefs(8,:) = [0.d0, 0.003183284637742d0,   0.027692298496089d0]
    tauCoefs   = [0.067913554086575d0, 0.039551521566859d0, 0.036743664854448d0, &
      0.098798200316163d0, 0.003773886623955d0, 0.088516800820248d0, &
      0.067873515485777d0,   0.075774013057833d0,   0.074313246812492d0]

    ! Set the value for the fields
    u = zero
    v = zero
    w = zero
    wC = zero
    T = zero
    !u  = 3.0d0 + coefs(1,1)*z  + &
    !  coefs(1,2)*(cos(wavenums(1)*x )*cos(wavenums(1)*y )*cos(wavenums(1)*(z-zmin-Lz/2.d0))) + &
    !  coefs(1,3)*(cos(wavenums(2)*y )*cos(wavenums(2)*x )*cos(wavenums(1)*(z-zmin-Lz/2.d0)))
    !v  = 2.5d0 + coefs(2,1)*z  + &
    !  coefs(2,2)*(cos(wavenums(1)*x )*cos(wavenums(1)*y )*cos(wavenums(1)*(z-zmin-Lz/2.d0))) + &
    !  coefs(2,3)*(cos(wavenums(2)*y )*cos(wavenums(2)*x )*cos(wavenums(1)*(z-zmin-Lz/2.d0)))
    !u  = 3.0d0 + coefs(1,1)*z  + coefs(1,2)*(cos(wavenums(1)*x ) + sin(wavenums(1)*y )) + coefs(1,3)*(sin(wavenums(2)*y ) + cos(wavenums(2)*x ))
    !v  = 2.5d0 + coefs(2,1)*z  + coefs(2,2)*(cos(wavenums(1)*x ) + sin(wavenums(1)*y )) + coefs(2,3)*(sin(wavenums(2)*y ) + cos(wavenums(2)*x ))
    !wC = 2.6d0 + coefs(3,1)*z  + coefs(3,2)*(cos(wavenums(1)*x ) + sin(wavenums(1)*y )) + coefs(3,3)*(sin(wavenums(2)*y ) + cos(wavenums(2)*x ))
    !w  = 2.6d0 + coefs(3,1)*zE + coefs(3,2)*(cos(wavenums(1)*xE) + sin(wavenums(1)*yE)) + coefs(3,3)*(sin(wavenums(2)*yE) + cos(wavenums(2)*xE))

    ! Allocated TinZ and read in Tinit data
    call message_min_max(1,"Bounds for u:", p_minval(minval(u)), p_maxval(maxval(u)))
    call message_min_max(1,"Bounds for v:", p_minval(minval(v)), p_maxval(maxval(v)))
    call message_min_max(1,"Bounds for w:", p_minval(minval(w)), p_maxval(maxval(w)))
    call message_min_max(1,"Bounds for T:", p_minval(minval(T)), p_maxval(maxval(T)))

    call message(0,"Velocity Field Initialized")
    nullify(T)
    nullify(u,v,w,wC,x,y,z)
end subroutine


subroutine set_planes_io(xplanes, yplanes, zplanes)
    use test_stats_xy_parameters
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 1, nyplanes = 1, nzplanes = 1

    allocate(yplanes(nyplanes))
    allocate(xplanes(nxplanes))
    allocate(zplanes(nzplanes))
    yplanes = [nySize/2]
    xPlanes = [nxSize/2] !800
    zplanes = [nzSize/2]
end subroutine

subroutine set_KS_planes_io(planesCoarseGrid, planesFineGrid)
    integer, dimension(:), allocatable,  intent(inout) :: planesFineGrid
    integer, dimension(:), allocatable,  intent(inout) :: planesCoarseGrid
    
    allocate(planesCoarseGrid(1), planesFineGrid(1))
    planesCoarseGrid = [8]
    planesFineGrid = [16]

end subroutine

subroutine setInhomogeneousNeumannBC_Temp(inputfile, wTh_surf)
    use kind_parameters,    only: rkind, clen
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: wTh_surf
    real(rkind) :: Lx, Ly, Lz, zmin, Tref
    logical :: symmetricDomain
    integer :: iounit
    character(len=clen) :: stats_info_dir
    integer :: num_stats_instances = 1

    namelist /SMinput/ Lx, Ly, Lz, symmetricDomain, zmin, Tref, stats_info_dir, num_stats_instances
    
    wTh_surf = zero
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=SMinput)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine

subroutine setDirichletBC_Temp(inputfile, Tfield, Tsurf, dTsurf_dt, whichSide)
    use kind_parameters,    only: rkind, clen
    use constants,          only: zero, one
    use reductions,         only: p_minval, p_maxval
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:), intent(in) :: Tfield
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    character(len=3), intent(in) :: whichSide
    integer :: iounit
    
    dTsurf_dt = zero
    if (whichSide == 'top') then
        Tsurf = p_maxval(maxval(Tfield))
    elseif (whichSide == 'bot') then
        Tsurf = p_minval(minval(Tfield))
    end if

end subroutine


subroutine set_Reference_Temperature(inputfile, Trefout)
    use kind_parameters,    only: rkind, clen
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Trefout
    real(rkind) :: Tref = 1.d0, Lx, Ly, Lz, zmin
    logical :: symmetricDomain
    integer :: iounit
    character(len=clen) :: stats_info_dir
    integer :: num_stats_instances = 1

    namelist /SMinput/ Lx, Ly, Lz, symmetricDomain, zmin, Tref, stats_info_dir, num_stats_instances

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=SMinput)
    close(ioUnit)    
    
    Trefout = Tref

end subroutine

subroutine hook_probes(inputfile, probe_locs)
    use kind_parameters,    only: rkind, clen
    use basic_io,           only: read_2d_ascii
    real(rkind), dimension(:,:), allocatable, intent(inout) :: probe_locs
    character(len=*),                intent(in)    :: inputfile
    character(len=clen) :: probefile
    logical :: useProbes = .false.
    integer :: ioUnit

    ! Unused namelist variables, but still need to declare them:
    integer :: vizDump_Schedule, t_restartDump, t_dataDump, ioType, runID
    integer :: t_planeDump, t_stop_planeDump, t_start_planeDump, t_start_pointProbe, t_stop_pointProbe, t_pointProbe
    logical :: dumpPlanes, dump_NU_SGS, dump_KAPPA_SGS
    real(rkind) :: deltaT_dump

    namelist /IO/ vizDump_Schedule, deltaT_dump, t_restartDump, t_dataDump, ioType, dumpPlanes, runID, useProbes, &
                & dump_NU_SGS, dump_KAPPA_SGS, t_planeDump, t_stop_planeDump, t_start_planeDump, t_start_pointProbe,&
                & t_stop_pointProbe, t_pointProbe
    namelist /Probes/ probefile

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=IO    )
    read(unit=ioUnit, NML=Probes)
    close(ioUnit)    
    
    ! IMPORTANT : Convention is to allocate probe_locs(3,nprobes)
    ! Example: If you have at least 3 probes:
    ! probe_locs(1,3) : x -location of the third probe
    ! probe_locs(2,3) : y -location of the third probe
    ! probe_locs(3,3) : z -location of the third probe


    ! Add probes here if needed
    ! Example code: The following allocates 2 probes at (0.1,0.1,0.1) and
    ! (0.2,0.2,0.2)  
    if (useProbes) then
        call read_2d_ascii(probe_locs,trim(probefile))
    end if

end subroutine

subroutine initScalar(decompC, inputfile, mesh, scalar_id, scalarField)
    use kind_parameters, only: rkind
    use decomp_2d,        only: decomp_info
    type(decomp_info),               intent(in)    :: decompC
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    integer, intent(in)                            :: scalar_id
    real(rkind), dimension(:,:,:), intent(out)     :: scalarField
    integer :: ioUnit = 11, ierr
    real(rkind) :: zlocScaInterface = 1.5d0
    real(rkind) :: interfaceThickness = 0.125d0
    real(rkind), dimension(:,:,:), pointer :: z => null()
    real(rkind) :: zmid, lf, kmin, kmax, tgtKE, tgtDissipation, gain, fringe_delta
    integer :: maskType
    logical :: dumpForce, projectDivergenceFree
    
    namelist /spectForceLayer/ zmid, lf, kmin, kmax, tgtKE, tgtDissipation, &
      gain, dumpForce, maskType, fringe_delta, projectDivergenceFree

    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
    read(unit=ioUnit, NML=spectForceLayer)
    close(ioUnit)

    associate( z => mesh(:,:,:,3))
        if (scalar_id < 4) then
            scalarField = 0.5d0*(1.d0 + tanh((z - zlocScaInterface)/interfaceThickness))
        elseif (scalar_id == 4) then
            scalarField = exp(-((z - zmid)/(lf/3.d0))**2.d0)
        end if
    end associate
end subroutine 

subroutine setScalar_source(decompC, inputfile, mesh, scalar_id, scalarSource)
    use kind_parameters, only: rkind
    use decomp_2d,        only: decomp_info
    use test_stats_xy_parameters, only: Sfunc
    use constants, only: pi
    use exits, only: message
    use reductions, only: p_sum

    type(decomp_info),               intent(in)    :: decompC
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    integer, intent(in)                            :: scalar_id
    real(rkind), dimension(:,:,:), intent(out)     :: scalarSource
    integer :: ioUnit = 11, ierr
    real(rkind), dimension(:,:,:), pointer :: z => null()
    real(rkind) :: zmid, lf, kmin, kmax, tgtKE, tgtDissipation, gain, fringe_delta
    integer :: maskType
    logical :: dumpForce, projectDivergenceFree
    
    namelist /spectForceLayer/ zmid, lf, kmin, kmax, tgtKE, tgtDissipation, &
      gain, dumpForce, maskType, fringe_delta, projectDivergenceFree

    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
    read(unit=ioUnit, NML=spectForceLayer)
    close(ioUnit)

    associate( z => mesh(:,:,:,3))
        scalarSource = exp(-((z - zmid)/(lf/3.d0))**2.d0)
    end associate

    call message(2,"Scalar source initialized") 
end subroutine 

subroutine hook_source(tsim,mesh,Re,urhs,vrhs,wrhs)
  use kind_parameters, only: rkind
  use decomp_2d,       only: decomp_info
  use constants,       only: zero
  real(rkind), intent(in)                             :: tsim, Re
  real(rkind), dimension(:,:,:,:), intent(in), target :: mesh
  real(rkind), dimension(:,:,:),   intent(inout)      :: urhs, vrhs, wrhs
  real(rkind), dimension(:,:,:), pointer              :: x => null(), y => null(), z => null()

  x => mesh(:,:,:,1)
  y => mesh(:,:,:,2)
  z => mesh(:,:,:,3)

  urhs = urhs + 0.d0 
  vrhs = vrhs + 0.d0 
  wrhs = wrhs + 0.d0

  nullify(x,y,z) 
end subroutine

