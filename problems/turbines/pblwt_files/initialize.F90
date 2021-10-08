module pblwt_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa 
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
    use pblwt_parameters    
    use kind_parameters,  only: rkind
    use constants,        only: one, two
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    real(rkind) :: z0init = 1.0d-4, ustarinit  = 1.0d0
    integer :: i,j,k, ioUnit
    character(len=*),                intent(in)    :: inputfile
    integer :: ix1, ixn, iy1, iyn, iz1, izn
    real(rkind)  :: Lx = one, Ly = one, Lz = one, zpeak
    namelist /PBLINPUT/ Lx, Ly, Lz, z0init, ustarinit, zpeak

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PBLINPUT)
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
    use pblwt_parameters
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two, pi, half
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random
    use decomp_2d          
    use reductions,         only: p_maxval, p_minval
    implicit none
    type(decomp_info),               intent(in)    :: decompC
    type(decomp_info),               intent(in)    :: decompE
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, x, y, z
    real(rkind), dimension(:,:,:), allocatable :: randArr
    real(rkind) :: z0init = 1.0d-4, epsnd = 0.1, sig, ustarinit = 1.0d0
    real(rkind), dimension(:,:,:), allocatable :: ybuffC, ybuffE, zbuffC, zbuffE
    integer :: nz, nzE, ioUnit, k
    real(rkind) :: Xperiods = 3.d0, Yperiods = 3.d0
    real(rkind) :: zpeak = 0.2d0
    real(rkind)  :: Lx = one, Ly = one, Lz = one
    namelist /PBLINPUT/ Lx, Ly, Lz, z0init, ustarinit, zpeak 

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PBLINPUT)
    close(ioUnit)    


    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)
 
    epsnd = 1.0d-5

    !u = (ustarinit/kappa)*log(z/z0init) + epsnd*cos(Yperiods*two*pi*y/Ly)*exp(-half*(z/zpeak/Lz)**2)
    u = 22.0d0*(one-(z-one)*(z-one)) + epsnd*cos(Yperiods*two*pi*y/Ly)*exp(-half*(z/zpeak/Lz)**2)
    v = epsnd*(z/Lz)*cos(Xperiods*two*pi*x/Lx)*exp(-half*(z/zpeak/Lz)**2)
    wC= zero  
   
    !Add random numbers
    randomScaleFact = 0.00d0
    allocate(randArr(size(u,1),size(u,2),size(u,3)))
    call gaussian_random(randArr,zero,one,seedu + 10*nrank)
    do k = 1,size(u,3)
        sig = randomScaleFact*(ustarinit/kappa)*log(z(1,1,k)/z0init)
        u(:,:,k) = u(:,:,k) + sig*randArr(:,:,k)
    end do  
    deallocate(randArr)
    
    allocate(randArr(size(v,1),size(v,2),size(v,3)))
    call gaussian_random(randArr,zero,one,seedv+ 10*nrank)
    do k = 1,size(v,3)
        sig = randomScaleFact*z(1,1,k)*exp(-half*(z(1,1,k)/zpeak/Lz)**2)
        v(:,:,k) = v(:,:,k) + sig*randArr(:,:,k)
    end do  
    deallocate(randArr)

    write(*,*) 'umax, umin: ', p_maxval(maxval(u)), p_minval(minval(u))
    write(*,*) 'vmax, vmin: ', p_maxval(maxval(v)), p_minval(minval(v))


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
    integer, parameter :: nxplanes = 1, nyplanes = 1, nzplanes = 3

    allocate(xplanes(nxplanes))
    allocate(yplanes(nxplanes))
    allocate(zplanes(nzplanes))

    xplanes = [144]
    yplanes = [48]
    zplanes = [5, 10, 15, 20]   ! approximately [0.1,0.2,0.3,0.4] for Lz=1, nz=48

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
    real(rkind) :: ThetaRef, Lx, Ly, Lz, z0init, zpeak
    integer :: iounit
    namelist /PBLINPUT/ Lx, Ly, Lz, z0init, zpeak
    
    wTh_surf = zero;
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PBLINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    real(rkind) :: ThetaRef, Lx, Ly, Lz, z0init, zpeak
    integer :: iounit
    namelist /PBLINPUT/ Lx, Ly, Lz, z0init, zpeak 
    
    Tsurf = zero; dTsurf_dt = zero; ThetaRef = one
    

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PBLINPUT)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine

subroutine set_Reference_Temperature(inputfile, Tref)
    use kind_parameters,    only: rkind
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    real(rkind) :: Lx, Ly, Lz, z0init, zpeak
    integer :: iounit
    
    namelist /PBLINPUT/ Lx, Ly, Lz, z0init , zpeak

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PBLINPUT)
    close(ioUnit)    
     
    Tref = 0.d0
    
    ! Do nothing really since this is an unstratified simulation

end subroutine

subroutine set_Reference_Temperatur(inputfile, Tref)
    use kind_parameters,    only: rkind
    use constants,          only: one, zero
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    real(rkind) :: Lx, Ly, Lz, z0init, zpeak
    integer :: iounit
    
    namelist /PBLINPUT/ Lx, Ly, Lz, z0init , zpeak

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PBLINPUT)
    close(ioUnit)    
     
    Tref = 0.d0
    
    ! Do nothing really since this is an unstratified simulation

end subroutine

subroutine hook_probes(inputfile, probe_locs)
    use kind_parameters,    only: rkind, clen
    use exits, only: message
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:), allocatable, intent(inout) :: probe_locs
    
    ! This code is specific to turbines where i am placing probes in front of
    ! each turbine. For a more generic implementation of probes check the
    ! initialization.F90 files for other igridWallM problems. 

    integer :: num_Turbines, nprobes = 2, ActuatorDiskID, ioUnit, ii, ADM_Type 
    logical :: useWindTurbines, ADM
    real(rkind) :: xloc, yloc, zloc, diam, ct, yaw, tilt
    real(rkind) :: upstreamdisplacement = 0.05d0 ! Place the probes this far upstream of the turbine centers
    character(len=clen) :: turbInfoDir, fname, tempname

    namelist /ACTUATOR_DISK/ xLoc, yLoc, zLoc, diam, cT, yaw, tilt
    namelist /WINDTURBINES/ useWindTurbines, num_turbines, ADM, turbInfoDir, ADM_Type
    
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=WINDTURBINES)
    close(ioUnit)    
    
   
    nprobes = 4*num_turbines
    allocate(probe_locs(3,nprobes))
    probe_locs = 0.d0   
   

    ! Set 1: Located at Hub height 
    do ActuatorDiskID = 1,num_turbines
        write(tempname,"(A13,I4.4,A10)") "ActuatorDisk_", ActuatorDiskID, "_input.inp"
        fname = turbInfoDir(:len_trim(turbInfoDir))//"/"//trim(tempname)

        ioUnit = 55
        open(unit=ioUnit, file=trim(fname), form='FORMATTED')
        read(unit=ioUnit, NML=ACTUATOR_DISK)
        close(ioUnit)
        
        probe_locs(1,ActuatorDiskID) = xLoc - upstreamdisplacement; 
        probe_locs(2,ActuatorDiskID) = yLoc; 
        probe_locs(3,ActuatorDiskID) = zLoc;
    end do 

    ! Set 2: Located 0.045 down from the hub 
    ii = 1 
    do ActuatorDiskID = num_turbines+1,2*num_turbines
        probe_locs(1,ActuatorDiskID) = probe_locs(1,ii)
        probe_locs(2,ActuatorDiskID) = probe_locs(2,ii)
        probe_locs(3,ActuatorDiskID) = probe_locs(3,ii) - 0.045d0
        ii = ii + 1
    end do  
        
    ! Set 3: Located 0.045 right from the hub 
    ii = 1 
    do ActuatorDiskID = 2*num_turbines+1,3*num_turbines
        probe_locs(1,ActuatorDiskID) = probe_locs(1,ii)
        probe_locs(2,ActuatorDiskID) = probe_locs(2,ii) + 0.045d0
        probe_locs(3,ActuatorDiskID) = probe_locs(3,ii) 
        ii = ii + 1
    end do  

    ! Set 4: Located 0.045 up from the hub 
    ii = 1 
    do ActuatorDiskID = 3*num_turbines+1,4*num_turbines
        probe_locs(1,ActuatorDiskID) = probe_locs(1,ii)
        probe_locs(2,ActuatorDiskID) = probe_locs(2,ii) 
        probe_locs(3,ActuatorDiskID) = probe_locs(3,ii) + 0.045d0 
        ii = ii + 1
    end do  

 
    call message(0,"Total number of probes desired:", nprobes)

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
