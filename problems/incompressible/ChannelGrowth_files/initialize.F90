module EkmanDNS_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa 
    use decomp_2d
    use reductions, only: p_sum
    implicit none
    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.002_rkind ! 0.2% of the mean value
    integer :: nxg, nyg, nzg
    
    real(rkind), parameter :: xdim = 1000._rkind, udim = 0.45_rkind
    real(rkind), parameter :: timeDim = xdim/udim

end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use EkmanDNS_parameters    
    use kind_parameters,  only: rkind, clen
    use constants,        only: one,two, pi
    use decomp_2d,        only: decomp_info
    use exits,            only: GracefulExit
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k, ioUnit
    character(len=*),                intent(in)    :: inputfile
    integer :: ix1, ixn, iy1, iyn, iz1, izn
    character(len=clen) :: ModeFile
    real(rkind) ::  Re,Lx = two*pi,Ly = two*pi, Pert_Amp, Noise_Amp = 1.d-6
    integer :: ProblemMode = 1
    namelist /ChannelGrowth/ ModeFile, Re,Lx,Ly, Pert_Amp, Noise_Amp, ProblemMode

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=ChannelGrowth)
    close(ioUnit)    

    !Lx = two*pi; Ly = two*pi; Lz = one

    nxg = decomp%xsz(1); nyg = decomp%ysz(2); nzg = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%xst(1); iy1 = decomp%xst(2); iz1 = decomp%xst(3)
    ixn = decomp%xen(1); iyn = decomp%xen(2); izn = decomp%xen(3)


    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = Lx/real(nxg,rkind)
        dy = Ly/real(nyg,rkind)
        dz = two/real(nzg,rkind)

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
        z = z - dz - one 

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
    use ChannelGrowth_IO, only: get_perturbations 
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
    real(rkind) :: dz 
    character(len=clen) :: ModeFile
    real(rkind), dimension(:,:,:), allocatable :: randArr, ybuffC, ybuffE, zbuffC, zbuffE
    integer :: nz, nzE
    real(rkind), dimension(:,:,:), allocatable :: upurt, vpurt, wpurt
    real(rkind) :: Re,Lx,Ly, Pert_Amp, Noise_Amp = 1.d-6
    type(cd06stagg), allocatable :: derW
    integer :: ProblemMode = 1
    real(rkind) :: kx = 1.d0, ky = 0.d0 

    namelist /ChannelGrowth/ ModeFile, Re,Lx,Ly, Pert_Amp, Noise_Amp, ProblemMode

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=ChannelGrowth)
    close(ioUnit)    


    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)
 
  
    ! Laminar Ekman layer profile

    u  = one - z*z
    v  = zero 
    wC = zero


    allocate(upurt(size(u ,1),size(u ,2),size(u ,3)))
    allocate(vpurt(size(v ,1),size(v ,2),size(v ,3)))
    allocate(wpurt(size(wC,1),size(wC,2),size(wC,3)))

    call get_perturbations(decompC, x, y, kx, ky, ModeFile, upurt, vpurt, wpurt)

    u  = u  + pert_Amp*upurt
    v  = v  + pert_Amp*vpurt
    wC = wC + pert_Amp*wpurt
    
    deallocate(upurt, vpurt, wpurt)
    allocate(randArr(size(wC,1),size(wC,2),size(wC,3)))
    
    call gaussian_random(randArr,zero,one,seedu + 100*nrank)
    u  = u + Noise_Amp*randArr
    
    call gaussian_random(randArr,zero,one,seedv + 100*nrank)
    v  = v + Noise_Amp*randArr

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
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 1, nyplanes = 1, nzplanes = 6

    allocate(xplanes(nxplanes), yplanes(nyplanes), zplanes(nzplanes))

    xplanes = [64]
    yplanes = [64]
    zplanes = [5,15,30,50,80,150]

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
    real(rkind), intent(out) :: wTh_surf
    character(len=*),                intent(in)    :: inputfile
    integer :: ioUnit 
    real(rkind) :: Re,Lx,Ly, Pert_Amp, Noise_Amp = 1.d-6
    integer :: ProblemMode = 1
    namelist /ChannelGrowth/ Lx, Ly, Re,Lx,Ly, Pert_Amp, Noise_Amp, ProblemMode
     
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=ChannelGrowth)
    close(ioUnit)    

    wTh_surf = 0.d0 
    ! Do nothing really since this is an unstratified simulation
end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    integer :: iounit
    real(rkind) :: Re,Lx,Ly, Pert_Amp, Noise_Amp = 1.d-6, ThetaRef
    integer :: ProblemMode = 1
    namelist /ChannelGrowth/ Lx, Ly, Re,Lx,Ly, Pert_Amp, Noise_Amp, ProblemMode
    
    Tsurf = zero; dTsurf_dt = zero; ThetaRef = one
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=ChannelGrowth)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine


subroutine set_Reference_Temperature(inputfile, Tref)
    use kind_parameters,    only: rkind
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    integer :: iounit
    real(rkind) :: Re,Lx,Ly, Pert_Amp, Noise_Amp = 1.d-6
    integer :: ProblemMode = 1
    
    namelist /ChannelGrowth/ Lx, Ly, Re,Lx,Ly, Pert_Amp, Noise_Amp, ProblemMode

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=ChannelGrowth)
    close(ioUnit)    
     
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
    probe_locs(4,1) = 0.1d0; probe_locs(2,1) = 0.1d0; probe_locs(3,1) = 0.1d0;
    probe_locs(1,2) = 0.2d0; probe_locs(2,2) = 0.2d0; probe_locs(3,2) = 0.2d0;


end subroutine
