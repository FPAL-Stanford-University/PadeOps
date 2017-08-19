module HIT_AD_interact_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa 
    implicit none
    integer :: simulationID = 0

end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use HIT_AD_interact_parameters    
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
    real(rkind)  :: Lx = one, Ly = one, Lz = one, uInflow = one
    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow
    
    select case (simulationID) 
    case (1) 
      ioUnit = 11
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
      read(unit=ioUnit, NML=AD_CoriolisINPUT)
      close(ioUnit)    
    case (2)
      Lx = two*pi
      Ly = two*pi
      Lz = two*pi
    end select

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
    use HIT_AD_interact_parameters
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two, pi, half
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random
    use decomp_2d          
    use reductions,         only: p_maxval, p_minval
    use exits,              only: message_min_max, gracefulExit, message
    implicit none
    type(decomp_info),               intent(in)    :: decompC
    type(decomp_info),               intent(in)    :: decompE
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, x, y, z
    real(rkind)  :: Lx = one, Ly = one, Lz = one, uInflow = one
    integer :: ioUnit
    namelist /AD_CoriolisINPUT/ Lx, Ly, Lz, uInflow
    
    if (simulationID == 1) then
      
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
   
      u = uInflow 
      v = zero
      wC= zero  
      w = zero

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
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 0, nyplanes = 1, nzplanes = 1

    !allocate(xplanes(nxplanes), yplanes(nyplanes), zplanes(nzplanes))
    allocate(yplanes(nyplanes))
    !allocate(zplanes(nzplanes))
    !xplanes = [300,400,500,600,700]
    yplanes = [128]
    !zplanes = [128]

end subroutine

subroutine set_KS_planes_io(planesCoarseGrid, planesFineGrid)
    integer, dimension(:), allocatable,  intent(inout) :: planesFineGrid
    integer, dimension(:), allocatable,  intent(inout) :: planesCoarseGrid
    
    allocate(planesCoarseGrid(1), planesFineGrid(1))
    planesCoarseGrid = [8]
    planesFineGrid = [16]

end subroutine


subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    real(rkind) :: ThetaRef, Lx, Ly, Lz, uInflow = 1.d0
    integer :: iounit
    namelist /HIT_AD_interactINPUT/ Lx, Ly, Lz, uInflow 
    
    Tsurf = zero; dTsurf_dt = zero; ThetaRef = one
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=HIT_AD_interactINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine


subroutine set_Reference_Temperature(inputfile, Tref)
    use kind_parameters,    only: rkind
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    real(rkind) :: Lx, Ly, Lz, uInflow = 1.d0
    integer :: iounit
    
    namelist /HIT_AD_interactINPUT/ Lx, Ly, Lz, uInflow

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=HIT_AD_interactINPUT)
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
