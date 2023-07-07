module shearLayerEntrainment_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa, zero
    use fortran_assert,     only: assert
    implicit none
    integer :: nxSize = 256, nySize = 256, nzSize = 256
    real(rkind) :: Umax = 1.d0
    !real(rkind) :: forceZmin = 0.d0, forceZmax = 0.1d0
    logical :: forceLayerInitialized = .false.
    !integer :: kstC, kenC, kstE, kenE
    !real(rkind), dimension(:), allocatable :: zC, zE
    !real(rkind), dimension(:,:,:), allocatable :: maskC, maskE
    real(rkind), dimension(:,:,:), allocatable :: UbC, UbE, sechSq
    real(rkind), dimension(:,:,:), allocatable, target :: buff
    real(rkind) :: h_u, h_T, onebyRe

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

  subroutine initializeForceLayer(sim)
      use incompressibleGrid, only: igrid
      use forcingLayerMod, only: onThisRank, getStEndIndices
      type(igrid), intent(in), target :: sim
      real(rkind), dimension(:,:,:), pointer :: zE, z

      zE => sim%zE
      z => sim%mesh(:,:,:,3)
     
      allocate(UbC(sim%gpC%xsz(1),sim%gpC%xsz(2),sim%gpC%xsz(3))) 
      allocate(sechSq(sim%gpC%xsz(1),sim%gpC%xsz(2),sim%gpC%xsz(3))) 
      allocate(UbE(sim%gpE%xsz(1),sim%gpE%xsz(2),sim%gpE%xsz(3))) 
      UbC = 1.d0 - 0.25*(1.d0 + tanh(z/h_u))
      UbE = 1.d0 - 0.25*(1.d0 + tanh(zE/h_u))
      sechsq = (1.d0/cosh(z/h_u))**2.d0

      !allocate(maskC(sim%gpC%xsz(1),sim%gpC%xsz(2),sim%gpC%xsz(3)))
      !allocate(maskE(sim%gpE%xsz(1),sim%gpE%xsz(2),sim%gpE%xsz(3)))
     
      !if (onThisRank(zC,forceZmin,forceZmax)) then
      !    call getStEndIndices(zC,forceZmin,forceZmax,kstC,kenC)
      !    call getStEndIndices(zE,forceZmin,forceZmax,kstE,kenE)
      !    maskC(:,:,kstC:kenC) = Umax
      !    maskE(:,:,kstE:kenE) = Umax
      !else
      !    maskC = 0.d0
      !    maskE = 0.d0
      !end if

      forceLayerInitialized = .true.

      !deallocate(zC,zE)
      nullify(zE, z)

  end subroutine

  subroutine finalizeProblem()
    deallocate(UbC, UbE, sechSq)
  end subroutine
  
end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use shearLayerEntrainment_parameters    
    use kind_parameters,  only: rkind
    use constants,        only: one,two, pi
    use decomp_2d,        only: decomp_info
    use forcingLayerMod, only: onThisRank, getStEndIndices
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
    real(rkind) :: zmin = -one

    namelist /problemInput/ Lx, Ly, Lz, symmetricDomain, zmin, Umax, h_u, h_T

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=problemInput)
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
    
    !allocate(zC(size(mesh,3)))
    !zC = mesh(1,1,:,3)

    !if (onThisRank(zC,forceZmin,forceZmax)) then
    !    call getStEndIndices(zC,forceZmin,forceZmax,kstC,kenC)
    !else
    !    kstC = 2
    !    kenC = 1    
    !end if

end subroutine

subroutine initfields_wallM(decompC, decompE, inputfile, mesh, fieldsC, fieldsE)
    use shearLayerEntrainment_parameters
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
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, x, y, z, T
    real(rkind), dimension(:,:,:), allocatable :: ybuffC, ybuffE, zbuffC, zbuffE
    real(rkind)  :: Lx = one, Ly = one, Lz = one
!    real(rkind)  :: zTop_cell, zBot_cell, zMid
    type(cd06stagg), allocatable :: der
    real(rkind) :: dz, zmin
    integer :: ioUnit, ierr
    integer :: i, j, n, Nmodes = 4
    logical :: symmetricDomain
    character(len=clen) :: fname
    real(rkind), dimension(:), allocatable :: TinZ
    real(rkind), dimension(4) :: phi, A, mode

    ! INPUT namelist variables (only inputdir is needed)
    integer :: nx, ny, nz, nsteps
    real(rkind) :: tstop, dt, CFL, CviscDT
    character(len=clen) :: inputdir, outputdir
    integer :: prow, pcol, restartFile_TID, restartFile_RID
    logical :: useRestartFile

    ! PHYSICS namelist variables (only isStratified is needed)
    logical :: isInviscid, useCoriolis, useExtraForcing, useMoisture, useSGS, &
      useforcedStratification, useGeostrophicForcing, assume_fplane, &
      useHITForcing,useScalars,useHITRealSpaceLinearForcing,addExtraSourceTerm, &
      useImmersedBodies, useLocalizedForceLayer
    real(rkind) :: Re, Ro, Pr, Fr, Ra, PrandtlFluid, BulkRichardson, &
      G_geostrophic,G_alpha,dpFdx,dpFdy,dpFdz,latitude,frameAngle, &
      HITForceTimeScale, immersed_taufact
    integer :: BuoyancyTermType,buoyancyDirection,numberOfImmersedBodies
    logical :: isStratified = .false.

    ! Read input file
    namelist /problemInput/ Lx, Ly, Lz, symmetricDomain, zmin, Umax, h_u, h_T
    namelist /INPUT/ nx, ny, nz, tstop, dt, CFL, nsteps, inputdir, outputdir, prow, pcol, &
                    useRestartFile, restartFile_TID, restartFile_RID, CviscDT
    namelist /PHYSICS/isInviscid,useCoriolis,useExtraForcing,isStratified,&
      useMoisture,Re,Ro,Pr,Fr, Ra, useSGS, PrandtlFluid, BulkRichardson, &
      BuoyancyTermType,useforcedStratification, useGeostrophicForcing, &
      G_geostrophic, G_alpha, dpFdx,dpFdy,dpFdz,assume_fplane,latitude,&
      useHITForcing, useScalars, frameAngle, buoyancyDirection, &
      useHITRealSpaceLinearForcing, HITForceTimeScale, addExtraSourceTerm, &
      useImmersedBodies, numberOfImmersedBodies, immersed_taufact, &
      useLocalizedForceLayer

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=INPUT)
    read(unit=ioUnit, NML=PHYSICS)
    read(unit=ioUnit, NML=problemInput)
    close(ioUnit)    
   
    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1)
    T => fieldsC(:,:,:,7) 
      
    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)
        
    u = zero
    v = zero 
    w = zero 
    wC = zero

    phi = [1.977467242167018d0, 2.549665936805474d0, -2.343710955248922d0, 2.597317105635469d0]
    A = 0.001443375672974d0
    mode = [5.d0, 6.d0, 7.d0, 8.d0]
    do n = 1,Nmodes
        u  = u  + A(n)*cos(mode(n)*2.d0*pi*x/Lx + mode(n)*2.d0*pi*y/Ly + phi(n))
        v  = v  + A(n)*sin(mode(n)*2.d0*pi*x/Lx + mode(n)*2.d0*pi*y/Ly + phi(n))
        wC = wC + A(n)*cos(mode(n)*2.d0*pi*x/Lx + mode(n)*2.d0*pi*y/Ly + phi(n))
    end do

    u  = u*exp(-25.d0*z*z)
    v  = v*exp(-25.d0*z*z)
    wC = wC*exp(-25.d0*z*z)

    ! Allocate buffers and interpolate wC to w
    allocate(ybuffC(decompC%ysz(1),decompC%ysz(2),decompC%ysz(3)))
    allocate(ybuffE(decompE%ysz(1),decompE%ysz(2),decompE%ysz(3)))
    allocate(zbuffC(decompC%zsz(1),decompC%zsz(2),decompC%zsz(3)))
    allocate(zbuffE(decompE%zsz(1),decompE%zsz(2),decompE%zsz(3)))
    allocate(der)
    call der%init(decompC%xsz(1),x(2,1,1)-x(1,1,1),.false.,.false.,.true.,.true.)
    call transpose_x_to_y(wC,ybuffC,decompC)
    call transpose_y_to_z(ybuffC,zbuffC,decompC)
    call der%interpZ_C2E(zbuffC,zbuffE,decompC%zsz(1),decompC%zsz(2))
    call transpose_z_to_y(zbuffE,ybuffE,decompE)
    call transpose_y_to_x(ybuffE,w,decompE)


    ! Allocated TinZ and read in Tinit data
    if (isStratified) then
        allocate(TinZ(decompC%zsz(3)))
        write(fname,'(A)')trim(inputdir)//'/Tinit.txt'
        call read_1d_ascii(TinZ,trim(fname))
        do j = 1,decompC%xsz(2)
            do i = 1,decompC%xsz(1)
                T(i,j,:) = TinZ(decompC%xst(3):decompC%xen(3))
            end do
        end do
        deallocate(TinZ)
    else
        T = one ! Make this a function of Richardson number. 
    end if
    
    call message_min_max(1,"Bounds for u:", p_minval(minval(u)), p_maxval(maxval(u)))
    call message_min_max(1,"Bounds for v:", p_minval(minval(v)), p_maxval(maxval(v)))
    call message_min_max(1,"Bounds for w:", p_minval(minval(w)), p_maxval(maxval(w)))
    call message_min_max(1,"Bounds for T:", p_minval(minval(T)), p_maxval(maxval(T)))

    call message(0,"Velocity Field Initialized")
    
    nullify(T)
    nullify(u,v,w,wC,x,y,z)
    deallocate(ybuffC,ybuffE,zbuffC,zbuffE)
    call der%destroy()
    deallocate(der)
end subroutine


subroutine set_planes_io(xplanes, yplanes, zplanes)
    use shearLayerEntrainment_parameters
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
    real(rkind) :: ThetaRef, Lx, Ly, Lz, zmin
    logical :: symmetricDomain
    integer :: iounit
    real(rkind) :: Umax, h_u, h_T

    namelist /problemInput/ Lx, Ly, Lz, symmetricDomain, zmin, Umax, h_u, h_T
    
    wTh_surf = zero
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=problemInput)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind, clen
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    real(rkind) :: ThetaRef, Lx, Ly, Lz, zmin
    logical :: symmetricDomain
    integer :: iounit
    real(rkind) :: Umax, h_u, h_T
    
    namelist /problemInput/ Lx, Ly, Lz, symmetricDomain, zmin, Umax, h_u, h_T
    
    Tsurf = zero; dTsurf_dt = zero; ThetaRef = one
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=problemInput)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine


subroutine set_Reference_Temperature(inputfile, Tref)
    use kind_parameters,    only: rkind, clen
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    real(rkind) :: ThetaRef, Lx, Ly, Lz, zmin
    logical :: symmetricDomain
    integer :: iounit
    real(rkind) :: Umax, h_u, h_T 
    
    namelist /problemInput/ Lx, Ly, Lz, symmetricDomain, zmin, Umax, h_u, h_T

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=problemInput)
    close(ioUnit)    
     
    Tref = 1.d0
    
    ! Do nothing really since this is an unstratified simulation

end subroutine

subroutine hook_probes(inputfile, probe_locs)
    use kind_parameters,    only: rkind
    real(rkind), dimension(:,:), allocatable, intent(inout) :: probe_locs
    character(len=*),                intent(in)    :: inputfile
    integer, parameter :: nprobes = 9
    
    ! IMPORTANT : Convention is to allocate probe_locs(3,nprobes)
    ! Example: If you have at least 3 probes:
    ! probe_locs(1,3) : x -location of the third probe
    ! probe_locs(2,3) : y -location of the third probe
    ! probe_locs(3,3) : z -location of the third probe


    ! Add probes here if needed
    ! Example code: The following allocates 2 probes at (0.1,0.1,0.1) and
    ! (0.2,0.2,0.2)  
    allocate(probe_locs(3,nprobes))
    
    ! Probe 1
    probe_locs(1,1) = 0.1d0; 
    probe_locs(2,1) = 3.141592653589d0; 
    probe_locs(3,1) = 3.141592653589d0; 
    
    ! Probe 2 
    probe_locs(1,2) = 5.783185307179d0; 
    probe_locs(2,2) = 3.141592653589d0; 
    probe_locs(3,2) = 3.141592653589d0; 

    ! Probe 3
    probe_locs(1,3) = 6.783185307179d0; 
    probe_locs(2,3) = 3.141592653589d0; 
    probe_locs(3,3) = 3.141592653589d0; 
    
    ! Probe 4
    probe_locs(1,4) = 10.28318530718d0; 
    probe_locs(2,4) = 3.141592653589d0; 
    probe_locs(3,4) = 3.141592653589d0; 

    ! Probe 5
    probe_locs(1,5) = 14.28318530718d0; 
    probe_locs(2,5) = 3.141592653589d0; 
    probe_locs(3,5) = 3.141592653589d0;

    ! Probe 6
    probe_locs(1,6) = 18.28318530718d0; 
    probe_locs(2,6) = 3.141592653589d0; 
    probe_locs(3,6) = 3.141592653589d0;

    ! Probe 7
    probe_locs(1,7) = 10.28318530718d0; 
    probe_locs(2,7) = 3.141592653589d0; 
    probe_locs(3,7) = 4.25d0; 

    ! Probe 8
    probe_locs(1,8) = 14.28318530718d0; 
    probe_locs(2,8) = 3.141592653589d0; 
    probe_locs(3,8) = 4.250;

    ! Probe 9
    probe_locs(1,9) = 18.28318530718d0; 
    probe_locs(2,9) = 3.141592653589d0; 
    probe_locs(3,9) = 4.25d0;
end subroutine

subroutine initScalar(decompC, inpDirectory, mesh, scalar_id, scalarField)
    use kind_parameters, only: rkind
    use decomp_2d,        only: decomp_info
    type(decomp_info),               intent(in)    :: decompC
    character(len=*),                intent(in)    :: inpDirectory
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    integer, intent(in)                            :: scalar_id
    real(rkind), dimension(:,:,:), intent(out)     :: scalarField

    scalarField = 0.d0
end subroutine 

subroutine setScalar_source(decompC, inputfile, mesh, scalar_id, scalarSource)
    use kind_parameters, only: rkind
    use decomp_2d,        only: decomp_info
    use shearLayerEntrainment_parameters, only: Sfunc
    use constants, only: pi
    use exits, only: message
    use reductions, only: p_sum

    type(decomp_info),               intent(in)    :: decompC
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    integer, intent(in)                            :: scalar_id
    real(rkind), dimension(:,:,:), intent(out)     :: scalarSource
    real(rkind), dimension(:,:,:), allocatable :: r, lambda, tmp
    real(rkind), dimension(:,:,:), pointer :: x, y, z
    real(rkind) :: xc = pi, yc = pi, zc = pi, rin = 0.75d0, rout = 1.25d0, delta_r = 0.22d0
    real(rkind) :: smear_x = 2.5d0, delta
    real(rkind) :: sumVal 

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)

    
    allocate(r(size(x,1),size(x,2),size(x,3)))
    allocate(lambda(size(x,1),size(x,2),size(x,3)))
    allocate(tmp(size(x,1),size(x,2),size(x,3)))

    r = sqrt((y - yc)**2 + (z - zc)**2)

    select case (scalar_id)
    case (1)
      tmp = (r - rout)/delta_r + 1 
      call Sfunc(tmp, lambda)
      lambda = -lambda
    case (2)
      tmp = (r - rin)/delta_r
      call Sfunc(tmp, lambda)
      lambda = 1.d0 - lambda
    end select 

    r = x - xc
    delta = (x(2,1,1) - x(1,1,1))*smear_x
    tmp = (1.d0/(delta*sqrt(2.d0*pi)))*exp(-0.5d0*(r**2)/(delta**2))
    scalarSource = tmp*lambda
    sumVal = p_sum(sum(scalarSource))*((x(2,1,1) - x(1,1,1))**3)

    call message(2,"Scalar source initialized with domain integrated value", sumVal)
    deallocate(r, lambda, tmp)

end subroutine 

subroutine hook_source(tsim,mesh,Re,w,urhs,vrhs,wrhs,duidxjC,duidxjE)
  use kind_parameters, only: rkind
  use decomp_2d,       only: decomp_info
  use constants,       only: zero
  use shearLayerEntrainment_parameters, only: forceLayerInitialized, UbC, UbE, &
    sechSq
  use fortran_assert, only: assert
  real(rkind), intent(in)                             :: tsim, Re
  real(rkind), dimension(:,:,:,:), intent(in), target :: mesh
  complex(rkind), dimension(:,:,:), intent(in)        :: w
  complex(rkind), dimension(:,:,:),   intent(inout)   :: urhs, vrhs, wrhs
  complex(rkind), dimension(:,:,:,:), intent(in)      :: duidxjC, duidxjE
  real(rkind), dimension(:,:,:), pointer              :: x, y, z

  x => mesh(:,:,:,1)
  y => mesh(:,:,:,2)
  z => mesh(:,:,:,3)
  
  ! The forcing terms are given here where uf and its derivatives are known
  ! analytically.
  ! fx = -w*ddz(uf) - uf*ddx(u)
  ! fy = -uf*ddx(v)
  ! fz = -uf*ddx(w)

  call assert(forceLayerInitialized,'forceLayerInitialized')
  urhs = urhs - (UbC*duidxjC(:,:,:,1) - w*0.25d0*sechSq/h_u)
  vrhs = vrhs - UbC*duidxjC(:,:,:,4)
  wrhs = wrhs - UbE*duidxjE(:,:,:,7)
  
  nullify(x,y,z)
end subroutine

