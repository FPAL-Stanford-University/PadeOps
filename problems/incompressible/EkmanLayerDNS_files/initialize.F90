module EkmanDNS_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa 
    implicit none
    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.002_rkind ! 0.2% of the mean value
    integer :: nxg, nyg, nzg
 
    !character(len=clen) :: inputfname   
    real(rkind), parameter :: xdim = 1000._rkind, udim = 0.45_rkind
    real(rkind), parameter :: timeDim = xdim/udim

end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use EkmanDNS_parameters  
    use EkmanLayerDNS_IO, only: read_domain_info  
    use kind_parameters,  only: rkind, clen
    use constants,        only: one,two, pi
    use decomp_2d,        only: decomp_info
    use exits,            only: GracefulExit
    implicit none

    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    character(len=*),                intent(in)    :: inputfile
    integer :: i,j,k, ioUnit
    integer :: ix1, ixn, iy1, iyn, iz1, izn
    real(rkind)  :: Lx=26.d0, Ly=26.d0, Lz=24.d0, alphaRot = 0.d0, noiseAmp 
    character(len=clen) :: InitFileTag, InitFileDirectory!,inputfname
    logical :: do_laminar = .false.
    
    namelist /EkmanLayerDNS/Lx,Ly,Lz,alphaRot,noiseAmp,InitFileTag, InitFileDirectory, do_laminar
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=EkmanLayerDNS)
    close(ioUnit)    

    !inputfname = trim(inputfile)
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
    use EkmanDNS_parameters
    use kind_parameters,    only: rkind, clen
    use constants,          only: zero, one, two, pi, half
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random
    use decomp_2d          
    use reductions,         only: p_maxval, p_minval
    use exits,              only: message_min_max
    use EkmanLayerDNS_IO,   only: get_perturbations
    use cd06staggstuff,      only: cd06stagg
    use decomp_2d_io,        only: decomp_2d_read_one

    implicit none
    type(decomp_info),               intent(in)    :: decompC
    type(decomp_info),               intent(in)    :: decompE
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    integer :: ioUnit
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, x, y, z
    real(rkind) :: alphaRot=0.d0, dz,noiseAmp 
    character(len=clen) :: InitFileTag, InitFileDirectory
    real(rkind), dimension(:,:,:), allocatable :: randArr, ybuffC, ybuffE, zbuffC, zbuffE
    integer :: nz, nzE
    real(rkind), dimension(:,:,:), allocatable :: upurt, vpurt, wpurt
    type(cd06stagg), allocatable :: derW
    real(rkind)  :: Lx = 26.d0, Ly = 26.d0, Lz = 24.d0 
    logical :: do_laminar = .false.
    namelist /EkmanLayerDNS/Lx,Ly,Lz,alphaRot,noiseAmp,InitFileTag, InitFileDirectory, do_laminar

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=EkmanLayerDNS)
    close(ioUnit)    


    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)
 
    ! Laminar Ekman layer profile
    u  = cos(alphaRot*pi/180.d0) - exp(-z)*cos(z-(alphaRot*pi/180.d0)) 
    v  = sin(alphaRot*pi/180.d0) + exp(-z)*sin(z-(alphaRot*pi/180.d0)) 
    wC = zero

    ! Perturbations and noise
    allocate(upurt(size(u ,1),size(u ,2),size(u ,3)))
    allocate(vpurt(size(v ,1),size(v ,2),size(v ,3)))
    allocate(wpurt(size(wC,1),size(wC,2),size(wC,3)))
    if (do_laminar) then
      upurt = 0.d0
      vpurt = 0.d0
      wpurt = 0.d0
    else
      call get_perturbations(decompC, x, y, InitFileTag, InitFileDirectory, upurt, vpurt, wpurt)
    end if
    u  = u + upurt
    v  = v + vpurt
    wC = wC+ wpurt
    
    deallocate(upurt, vpurt, wpurt)
    allocate(randArr(size(wC,1),size(wC,2),size(wC,3)))
    
    call gaussian_random(randArr,zero,one,seedu + 100*nrank)
    u  = u + noiseAmp*randArr
    
    call gaussian_random(randArr,zero,one,seedv + 100*nrank)
    v  = v + noiseAmp*randArr

    deallocate(randArr)
    

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
    allocate(derW)
    dz = z(1,1,2) - z(1,1,1)
    call derW%init(decompC%zsz(3), dz, isTopEven = .false., isBotEven = .true., &
                                isTopSided = .false., isBotSided = .false.)

    call transpose_x_to_y(wC,ybuffC,decompC)
    call transpose_y_to_z(ybuffC,zbuffC,decompC)
    !zbuffE = zero
    !zbuffE(:,:,2:nzE-1) = half*(zbuffC(:,:,1:nz-1) + zbuffC(:,:,2:nz))
    call derW%interpz_C2E(zbuffC,zbuffE,size(zbuffC,1),size(zbuffC,2))
    zbuffE(:,:,1) = zero
    call transpose_z_to_y(zbuffE,ybuffE,decompE)
    call transpose_y_to_x(ybuffE,w,decompE) 
    call derW%destroy()
    deallocate(derW)

    deallocate(ybuffC,ybuffE,zbuffC, zbuffE) 
  

    nullify(u,v,w,x,y,z)
   

    call message(0,"Velocity Field Initialized")

end subroutine


subroutine set_planes_io(xplanes, yplanes, zplanes)
    use EkmanDNS_parameters  
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 0, nyplanes = 0, nzplanes = 4
    
    allocate(zplanes(nzplanes))
    zplanes =[1, 53, 106, 500]

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
    integer :: iounit, directionID
    namelist /ScalarFieldTestingINPUT/ directionID
     
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=ScalarFieldTestingINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    real(rkind) :: ThetaRef, Lx, Ly, Lz, alphaRot
    
    Tsurf = zero
    dTsurf_dt = zero
    ThetaRef = one
    
    ! Do nothing really since this is an unstratified simulation
end subroutine


subroutine set_Reference_Temperature(inputfile, Tref)
    use kind_parameters,    only: rkind
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
     
    Tref = 0.d0
    
    ! Do nothing really since this is an unstratified simulation

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
