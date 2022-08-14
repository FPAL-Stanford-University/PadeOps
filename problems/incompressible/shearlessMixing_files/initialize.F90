module shearlessMixing_interact_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa, zero
    implicit none
    integer :: simulationID = 0
    integer :: nxSize = 256, nySize = 256, nzSize = 256
    character(len=1) :: streamWiseCoord = 'x'
    integer :: nxADSIM = 0, nyADSIM = 0, nzADSIM = 0, &
               nxHIT = 0,   nyHIT = 0,   nzHIT = 0

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

subroutine copyHITfieldsToSM(uhit,vhit,whit,uAD,vAD,wAD,hit,adsim,coord)
  use incompressibleGrid, only: igrid
  use decomp_2d,          only: transpose_x_to_y, transpose_y_to_x,&
                                transpose_y_to_z, transpose_z_to_y
  real(rkind), dimension(:,:,:), intent(in) :: uhit, vhit, whit
  real(rkind), dimension(:,:,:), intent(inout) :: uAD, vAD, wAD
  character(len=1), intent(in) :: coord
  type(igrid), intent(inout) :: hit, adsim
  integer :: ist, ien, jst, jen, kst, ken
  
  select case (coord)
    case ('x')
      ist = nxADSIM-nxHIT+1
      ien = nxADSIM
     
      ! Since everything is in x-decomposition anyway, we can just copy
      uAD(ist:ien,:,:) = uhit  
      vAD(ist:ien,:,:) = vhit  
      wAD(ist:ien,:,:) = whit  

    case ('y')
      jst = nyADSIM-nyHIT+1
      jen = nyADSIM

      ! We need to transpose from x->y then copy 
      call transpose_x_to_y(uhit, hit%rbuffyC(:,:,:,1), hit%gpC)
      adsim%rbuffyC(:,:,:,1) = zero
      adsim%rbuffyC(:,jst:jen,:,1) =  hit%rbuffyC(:,:,:,1)
      call transpose_y_to_x(adsim%rbuffyC(:,:,:,1), uAD, adSim%gpC)

      call transpose_x_to_y(vhit, hit%rbuffyC(:,:,:,1), hit%gpC)
      adsim%rbuffyC(:,:,:,1) = zero
      adsim%rbuffyC(:,jst:jen,:,1) =  hit%rbuffyC(:,:,:,1)
      call transpose_y_to_x(adsim%rbuffyC(:,:,:,1), vAD, adSim%gpC)
      
      call transpose_x_to_y(whit, hit%rbuffyE(:,:,:,1), hit%gpE)
      adsim%rbuffyE(:,:,:,1) = zero
      adsim%rbuffyE(:,jst:jen,:,1) =  hit%rbuffyE(:,:,:,1)
      call transpose_y_to_x(adsim%rbuffyE(:,:,:,1), wAD, adSim%gpE)

    case ('z')
      kst = nzADSIM-nzHIT+1
      ken = nzADSIM

      ! We need to transpose from x->y->z then copy 
      call transpose_x_to_y(uhit, hit%rbuffyC(:,:,:,1), hit%gpC)
      call transpose_y_to_z(hit%rbuffyC(:,:,:,1), hit%rbuffzC(:,:,:,1), hit%gpC)
      adsim%rbuffzC(:,:,:,1) = zero 
      adsim%rbuffzC(:,:,kst:ken,1) = hit%rbuffzC(:,:,:,1) 
      call transpose_z_to_y(adsim%rbuffzC(:,:,:,1),adsim%rbuffyC(:,:,:,1),adSim%gpC)
      call transpose_y_to_x(adsim%rbuffyC(:,:,:,1), uAD, adSim%gpC)

      call transpose_x_to_y(vhit, hit%rbuffyC(:,:,:,1), hit%gpC)
      call transpose_y_to_z(hit%rbuffyC(:,:,:,1), hit%rbuffzC(:,:,:,1), hit%gpC)
      adsim%rbuffzC(:,:,:,1) = zero 
      adsim%rbuffzC(:,:,kst:ken,1) = hit%rbuffzC(:,:,:,1) 
      call transpose_z_to_y(adsim%rbuffzC(:,:,:,1),adsim%rbuffyC(:,:,:,1),adSim%gpC)
      call transpose_y_to_x(adsim%rbuffyC(:,:,:,1), vAD, adSim%gpC)

      call transpose_x_to_y(whit, hit%rbuffyE(:,:,:,1), hit%gpE)
      call transpose_y_to_z(hit%rbuffyE(:,:,:,1), hit%rbuffzE(:,:,:,1), hit%gpE)
      adsim%rbuffzE(:,:,:,1) = zero 
      adsim%rbuffzE(:,:,kst:ken,1) = hit%rbuffzE(:,:,:,1) 
      call transpose_z_to_y(adsim%rbuffzE(:,:,:,1),adsim%rbuffyE(:,:,:,1),adSim%gpE)
      call transpose_y_to_x(adsim%rbuffyE(:,:,:,1), wAD, adSim%gpE)


  end select
end subroutine

end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use shearlessMixing_interact_parameters    
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
    real(rkind)  :: Lx = one, Ly = one, Lz = one, zero = one
    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, zero, InflowProfileType, &
            InflowProfileAmplit, InflowProfileThick, streamWiseCoord
    namelist /HIT_PeriodicINPUT/ Lx, Ly, Lz 

    select case (simulationID) 
    case (1) 
      ioUnit = 11
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
      read(unit=ioUnit, NML=AD_CoriolisINPUT)
      close(ioUnit)    
    case (2)
      ioUnit = 11
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
      read(unit=ioUnit, NML=HIT_PeriodicINPUT)
      close(ioUnit)    
    end select

    nxg = decomp%xsz(1); nyg = decomp%ysz(2); nzg = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%xst(1); iy1 = decomp%xst(2); iz1 = decomp%xst(3)
    ixn = decomp%xen(1); iyn = decomp%xen(2); izn = decomp%xen(3)
    
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = Lx/real(nxg,rkind)
        dy = Ly/real(nyg,rkind)
        ! z goes from -Lz to Lz
        dz = 2.d0*Lz/real(nzg,rkind)
        
        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                    y(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( -nzg/2, rkind ) * dz + real( ik1 + k - 1, rkind ) * dz + dz/two
                end do
            end do
        end do

        ! Shift everything to the origin 
        x = x - dx
        y = y - dy
        z = z - dz 

    end associate
   
    if (simulationID == 1) then
      nxSize = nxg; nySize = nyg; nzSize = nzg
    end if 

end subroutine

subroutine initfields_wallM(decompC, decompE, inputfile, mesh, fieldsC, fieldsE)
    use shearlessMixing_interact_parameters
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two, pi, half
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random
    use decomp_2d          
    use reductions,         only: p_maxval, p_minval
    use exits,              only: message_min_max, gracefulExit, message
    use fortran_assert,     only: assert
    implicit none
    type(decomp_info),               intent(in)    :: decompC
    type(decomp_info),               intent(in)    :: decompE
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, x, y, z
    real(rkind)  :: Lx = one, Ly = one, Lz = one
    real(rkind)  :: zTop_cell, zBot_cell, zMid
    integer :: ioUnit
    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, InflowProfileType, &
            InflowProfileAmplit, InflowProfileThick, streamWiseCoord
    
    if (simulationID == 1) then
      
      ioUnit = 11
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
      read(unit=ioUnit, NML=AD_CoriolisINPUT)
      close(ioUnit)    

      u  => fieldsC(:,:,:,1)
      v  => fieldsC(:,:,:,2)
      wC => fieldsC(:,:,:,3)
      w  => fieldsE(:,:,:,1)
      T  => fieldsC(:,:,:,7) 

      z => mesh(:,:,:,3)
      y => mesh(:,:,:,2)
      x => mesh(:,:,:,1)
 
      zTop_cell = p_maxval(mesh(:,:,:,3))
      zBot_cell = p_minval(mesh(:,:,:,3))
      zMid = half*(zTop_cell + zBot_cell)
 
      select case(InflowProfileType)
      T = one ! Make this a function of Richardson number. 
      u = zero
      v = zero 
      w = zero 
      wC = zero 

      call message_min_max(1,"Bounds for u:", p_minval(minval(u)), p_maxval(maxval(u)))
      call message_min_max(1,"Bounds for v:", p_minval(minval(v)), p_maxval(maxval(v)))
      call message_min_max(1,"Bounds for w:", p_minval(minval(w)), p_maxval(maxval(w)))
      

      nullify(u,v,w,x,y,z)
      call message(0,"Velocity Field for Simulation 1 Initialized")
    else 
      call gracefulExit("Only the Actuator disk simulation can be initialized like this. Check input file for HIT to ensure that it is restarted.",13)
    end if

end subroutine


subroutine set_planes_io(xplanes, yplanes, zplanes)
    use shearlessMixing_interact_parameters
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 1, nyplanes = 1, nzplanes = 1

    if (simulationID == 1) then
         !allocate(xplanes(nxplanes), yplanes(nyplanes), zplanes(nzplanes))
         allocate(yplanes(nyplanes))
         allocate(xplanes(nxplanes))
         !allocate(zplanes(nzplanes))
         !xplanes = [300,400,500,600,700]
         yplanes = [nySize/2]
         xPlanes = [5*nxSize/8] !800
         !xPlanes = [1, ceiling(nxSize/6.0), ceiling(nxSize/4.85), ceiling(nxSize/4.0), ceiling(nxSize/3.0), ceiling(nxSize/2.0), ceiling(nxSize/1.6)]
         !zplanes = [128]
    end if 
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
    real(rkind) :: ThetaRef, Lx, Ly, Lz, zero = 1.d0
    integer :: iounit
    namelist /shearlessMixing_interactINPUT/ Lx, Ly, Lz, zero 
    
    wTh_surf = zero
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=shearlessMixing_interactINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    real(rkind) :: ThetaRef, Lx, Ly, Lz, zero = 1.d0
    integer :: iounit
    namelist /shearlessMixing_interactINPUT/ Lx, Ly, Lz, zero 
    
    Tsurf = zero; dTsurf_dt = zero; ThetaRef = one
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=shearlessMixing_interactINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine


subroutine set_Reference_Temperature(inputfile, Tref)
    use kind_parameters,    only: rkind
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    real(rkind) :: Lx, Ly, Lz, zero = 1.d0
    integer :: iounit
    
    namelist /shearlessMixing_interactINPUT/ Lx, Ly, Lz, zero

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=shearlessMixing_interactINPUT)
    close(ioUnit)    
     
    Tref = 0.d0
    
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
    type(decomp_info),                                          intent(in)    :: decompC
    character(len=*),                intent(in)    :: inpDirectory
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    integer, intent(in)                            :: scalar_id
    real(rkind), dimension(:,:,:), intent(out)     :: scalarField

    scalarField = 0.d0
end subroutine 

subroutine setScalar_source(decompC, inputfile, mesh, scalar_id, scalarSource)
    use kind_parameters, only: rkind
    use decomp_2d,        only: decomp_info
    use shearlessMixing_interact_parameters, only: Sfunc
    use constants, only: pi
    use exits, only: message
    use reductions, only: p_sum

    type(decomp_info),                                          intent(in)    :: decompC
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
