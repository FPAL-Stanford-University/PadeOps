module HIT_Periodic_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa 
    implicit none
    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.002_rkind ! 0.2% of the mean value
    integer :: nxg, nyg, nzg
 
    logical :: useBandpassFilter = .false. 
    real(rkind) :: k_bp_left, k_bp_right, uadvect = 10.0, x_shift 
    real(rkind), dimension(:,:,:), allocatable :: uTarget, vTarget, wTarget

end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use HIT_Periodic_parameters    
    use kind_parameters,  only: rkind, clen
    use constants,        only: one,two, pi
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k, ioUnit
    character(len=*),                intent(in)    :: inputfile
    integer :: ix1, ixn, iy1, iyn, iz1, izn
    real(rkind)  :: Lx = two*pi, Ly = two*pi, Lz = two*pi
    character(len=clen)  :: dir_init_files
    real(rkind) :: uadv = 1.d0, kleft = 10.d0, kright = 64.d0
    character(len=clen)  :: ufname, vfname, wfname
    logical :: BandpassFilterFields = .false. 
    integer :: initType = 0
    namelist /HIT_PeriodicINPUT/ ufname, vfname, wfname, uadv, kleft, kright, BandpassFilterFields, Lx, Ly, Lz, initType

    !Lx = two*pi; Ly = two*pi; Lz = one
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=HIT_PeriodicINPUT)
    close(ioUnit)    

    !Lx = two*pi; Ly = two*pi; Lz = two*pi
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

    k_bp_left = kleft
    k_bp_right = kright
    uadvect = uadv
    useBandPassFilter = BandpassFilterFields
    if (useBandPassFilter) then
      allocate(uTarget(size(mesh,1), size(mesh,2), size(mesh,3)))
      allocate(vTarget(size(mesh,1), size(mesh,2), size(mesh,3)))
      allocate(wTarget(size(mesh,1), size(mesh,2), size(mesh,3)))
      uTarget = 0.d0
      vTarget = 0.d0
      wTarget = 0.d0
    end if 
end subroutine

subroutine initfields_wallM(decompC, decompE, inputfile, mesh, fieldsC, fieldsE)
    use HIT_Periodic_parameters
    use PadeDerOps, only: Pade6Stagg
    use kind_parameters,    only: rkind, clen 
    use constants,          only: zero, one, two, pi, half
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random
    use decomp_2d          
    use decomp_2d_io
    use reductions,         only: p_maxval, p_minval, p_sum
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
    real(rkind) :: dz
    real(rkind), dimension(:,:,:), allocatable :: randArr, ybuffC, ybuffE, zbuffC, zbuffE
    type(cd06stagg), allocatable :: der
    integer :: nz, nzE, k
    character(len=clen)  :: ufname, vfname, wfname 
    real(rkind) :: uadv = 0.d0, kleft = 10.d0, kright = 64.d0, Lx, Ly, Lz
    logical :: BandpassFilterFields = .false. 
    integer :: initType = 0, seed = 23455
    namelist /HIT_PeriodicINPUT/ ufname, vfname, wfname, uadv, kleft, kright, BandpassFilterFields, Lx, Ly, Lz, initType

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=HIT_PeriodicINPUT)
    close(ioUnit)    


    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)
    
    if (initType == 0) then
        
       dz = z(1,1,2) - z(1,1,1)

       
       call decomp_2d_read_one(1,u ,ufname, decompC)
       call decomp_2d_read_one(1,v ,vfname, decompC)
       call decomp_2d_read_one(1,wC,wfname, decompC)
  
       call message_min_max(1,"Bounds for u:", p_minval(minval(u)), p_maxval(maxval(u)))
       call message_min_max(1,"Bounds for v:", p_minval(minval(v)), p_maxval(maxval(v)))
       call message_min_max(1,"Bounds for w:", p_minval(minval(wC)), p_maxval(maxval(wC)))
       
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
       allocate(der)
       call der%init(decompC%zsz(3), dz, isTopEven = .false., isBotEven = .false., &
                                isTopSided = .false., isBotSided = .false.)
       call der%interpZ_C2E(zbuffC,zbuffE,size(zbuffC,1),size(zbuffC,2))                         
       deallocate(der)
       call transpose_z_to_y(zbuffE,ybuffE,decompE)
       call transpose_y_to_x(ybuffE,w,decompE) 
      

       deallocate(ybuffC,ybuffE,zbuffC, zbuffE) 
    else

        call uniform_random(u,-5.d0,5.d0,seed+1234*nrank+54321)
        call uniform_random(v,-5.d0,5.d0,seed+25634*nrank+54321)
        call uniform_random(w,-5.d0,5.d0,seed+32454*nrank+54321)
        
        u = u - p_sum(u)/(decompC%xsz(1)*decompC%ysz(2)*decompC%zsz(3))
        v = v - p_sum(v)/(decompC%xsz(1)*decompC%ysz(2)*decompC%zsz(3))
        w = w - p_sum(w)/(decompE%xsz(1)*decompE%ysz(2)*decompE%zsz(3))
    end if 
      
    nullify(u,v,w,x,y,z)
   
    call message(0,"Velocity Field Initialized")

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

    xplanes = [64]
    yplanes = [64]
    zplanes = [64]

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
    use constants, only: one, zero
    implicit none
    real(rkind), intent(out) :: wTh_surf
    character(len=clen),                intent(in)    :: inputfile
    integer :: ioUnit 
    character(len=clen)  :: ufname, vfname, wfname 
    real(rkind) :: TI = 0.1, uadv = 1.d0, kleft = 10.d0, kright = 64.d0, Lx, Ly, Lz
    logical :: BandpassFilterFields = .false. 
    namelist /HIT_PeriodicINPUT/ ufname, vfname, wfname, TI, uadv, kleft, kright, BandpassFilterFields, Lx, Ly, Lz
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=HIT_PeriodicINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind, clen 
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    real(rkind) :: ThetaRef
    integer :: iounit 
    character(len=clen)  :: ufname, vfname, wfname 
    real(rkind) :: uadv = 1.d0, kleft = 10.d0, kright = 64.d0, Lx, Ly, Lz
    logical :: BandpassFilterFields = .false. 
    namelist /HIT_PeriodicINPUT/ ufname, vfname, wfname, uadv, kleft, kright, BandpassFilterFields, Lx, Ly, Lz
    
    Tsurf = zero; dTsurf_dt = zero; ThetaRef = one
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=HIT_PeriodicINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine


subroutine set_Reference_Temperature(inputfile, Tref)
    use kind_parameters,    only: rkind, clen 
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    integer :: iounit
    character(len=clen)  :: ufname, vfname, wfname 
    real(rkind) :: uadv = 1.d0, kleft = 10.d0, kright = 64.d0, Lx, Ly, Lz
    logical :: BandpassFilterFields = .false. 
    namelist /HIT_PeriodicINPUT/ ufname, vfname, wfname, uadv, kleft, kright, BandpassFilterFields, Lx, Ly, Lz
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=HIT_PeriodicINPUT)
    close(ioUnit)    
     
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
