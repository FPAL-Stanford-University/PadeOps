module RBP_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa 
    use decomp_2d
    use reductions, only: p_sum
    use IncompressibleGrid, only: igrid
    implicit none
    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.002_rkind ! 0.2% of the mean value
    integer :: nxg, nyg, nzg
    
contains
    subroutine SetTemperatureBC_RBP(igp, inputfile)
        class(igrid), intent(inout) :: igp 
        character(len=*),                intent(in)    :: inputfile
        real(rkind) :: T_bottom = -0.5d0, T_top = 0.5d0
        integer :: ioUnit
        namelist /TEMPERATURE_BC/ T_bottom, T_top
        
        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=TEMPERATURE_BC)
        close(ioUnit)    
       
        call igp%ResetTemperatureBCs(T_top,T_bottom)
    end subroutine 
    
end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use RBP_parameters    
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
    character(len=clen) :: InitFileTag, InitFileDirectory
    character(len=clen) :: fname
    real(rkind)  :: Lx = one, Ly = one, Lz = one, Noise_Amp 
    integer :: ProblemMode = 1
    real(rkind) :: alpha = 1.0d0, phase = 0.0d0, x0 = 2.0d0, delta = 3.0d0, Pert_Amp = 1.d-4
    namelist /RBPinstability/ Lx, Ly, Noise_Amp, Pert_Amp, ProblemMode, alpha, phase, x0, delta

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=RBPinstability)
    close(ioUnit)    

    !Lx = two*pi; Ly = two*pi; Lz = one

    nxg = decomp%xsz(1); nyg = decomp%ysz(2); nzg = decomp%zsz(3)

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
        
        ! shift z to be centered about zero
        z = z - 0.5d0
    end associate

end subroutine

subroutine initfields_wallM(decompC, decompE, inputfile, mesh, fieldsC, fieldsE)
    use RBP_parameters
    use kind_parameters,    only: rkind, clen
    use constants,          only: zero, one, two, pi, half, imi
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random
    use decomp_2d          
    use reductions,         only: p_maxval, p_minval
    use exits,              only: message_min_max
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
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, x, y, z, T
    real(rkind) :: dz , speed, um, vm
    character(len=clen) :: InitFileTag, InitFileDirectory
    real(rkind), dimension(:,:,:), allocatable :: randArr, ybuffC, ybuffE, zbuffC, zbuffE
    integer :: nz, nzE
    real(rkind), dimension(:,:,:), allocatable :: upurt, vpurt, wpurt, Tpurt, xs, mask
    real(rkind) :: Noise_Amp = 0.d0, Lx, Ly
    type(cd06stagg), allocatable :: derW
    integer :: ProblemMode = 1
    character(len=clen) :: fname
    real(rkind) :: T_top, T_bottom
    real(rkind) :: alpha = 1.0d0, phase = 0.0d0, x0 = 2.0d0, delta = 3.0d0, Pert_Amp = 1.d-4
    namelist /RBPinstability/ Lx, Ly, Noise_Amp, Pert_Amp, ProblemMode, alpha, phase, x0, delta
    namelist /TEMPERATURE_BC/ T_bottom, T_top

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=RBPinstability)
    read(unit=ioUnit, NML=TEMPERATURE_BC)
    close(ioUnit)    
    


    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,3)
    T  => fieldsC(:,:,:,7) 
    w  => fieldsE(:,:,:,1)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)
 
  
    select case (ProblemMode)
    case (1)
         ! Laminar profile

         ! Base State
         u  = 1.d0 - 4.d0*z*z 
         v  = zero 
         wC = zero
         T  = T_bottom * (0.5 - z) + T_top * (z + 0.5) 

         allocate(upurt(size(u ,1),size(u ,2),size(u ,3)))
         allocate(vpurt(size(v ,1),size(v ,2),size(v ,3)))
         allocate(wpurt(size(wC,1),size(wC,2),size(wC,3)))
         allocate(Tpurt(size(T ,1),size(T ,2),size(T ,3)))
         allocate(xs   (size(x ,1),size(x ,2),size(x ,3)))
         allocate(mask (size(x ,1),size(x ,2),size(x ,3)))

         ! Set perturbation (wavepacket)
         xs = (x - x0) / delta
         mask = 0.5d0 - sign(0.5d0, abs(xs)-0.5d0)
         vpurt = 0.d0 
         Tpurt = 0.d0
         !upurt = (-0.25*imi*sqrt(pi)*delta) * exp(-0.25*alpha*(alpha*delta*delta + 4*imi*x0)) * (64.0*z*(z*z-0.25)) * exp(-imi*phase) * erf(-0.5*imi*alpha*delta - xs) - exp(2.0*imi*alpha*x0+imi*phase) * erf(0.5*imi*alpha*delta - xs)
         !wpurt = (16.0_rkind * (z*z - 0.25_rkind) ** 2_rkind) * exp(-xs*xs) * sin(alpha*x+phase)
         wpurt = (16.0_rkind * (z*z - 0.25_rkind) ** 2_rkind) * (0.5 + 0.5 * cos(2.0*pi*xs)) * sin(alpha*x+phase)
         upurt = (64.0_rkind * z * (z*z - 0.25_rkind)) * (0.25*delta/(2.0*pi-alpha*delta) * cos(alpha*x+phase-2.0*pi*xs) - 0.25*delta/(2.0*pi+alpha*delta) * cos(alpha*x+phase+2.0*pi*xs)) - 0.5/alpha * cos(alpha*x+phase)

         u  = u  + Pert_Amp*upurt*mask
         v  = v  + Pert_Amp*vpurt*mask
         wC = wC + Pert_Amp*wpurt*mask
         T  = T  + Pert_Amp*Tpurt*mask 

         deallocate(upurt, vpurt, wpurt, Tpurt, xs, mask)
         allocate(randArr(size(wC,1),size(wC,2),size(wC,3)))
         
         call gaussian_random(randArr,zero,one,seedu + 100*nrank)
         u  = u + Noise_Amp*randArr
         
         call gaussian_random(randArr,zero,one,seedv + 100*nrank)
         v  = v + Noise_Amp*randArr
        
         call gaussian_random(randArr,zero,one,seedv + seedu + 100*nrank)
         T  = T + Noise_Amp*randArr
         deallocate(randArr)
         

         call message_min_max(1,"Bounds for u:", p_minval(minval(u)), p_maxval(maxval(u)))
         call message_min_max(1,"Bounds for v:", p_minval(minval(v)), p_maxval(maxval(v)))
         call message_min_max(1,"Bounds for w:", p_minval(minval(wC)), p_maxval(maxval(wC)))
         call message_min_max(1,"Bounds for T:", p_minval(minval(T)), p_maxval(maxval(T)))
         
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
         call derW%interpz_C2E(zbuffC,zbuffE,size(zbuffC,1),size(zbuffC,2))
         zbuffE(:,:,1) = zero
         zbuffE(:,:,nzE) = zero
         call transpose_z_to_y(zbuffE,ybuffE,decompE)
         call transpose_y_to_x(ybuffE,w,decompE) 
         call derW%destroy()
         deallocate(derW)

         deallocate(ybuffC,ybuffE,zbuffC, zbuffE) 
  
    case (2)
        fname = InitFileDirectory(:len_trim(InitFileDirectory))//"/mystery_u.bin"
        call decomp_2d_read_one(1,u,fname, decompC)
        
        fname = InitFileDirectory(:len_trim(InitFileDirectory))//"/mystery_v.bin"
        call decomp_2d_read_one(1,v,fname, decompC)
        
        fname = InitFileDirectory(:len_trim(InitFileDirectory))//"/mystery_w.bin"
        call decomp_2d_read_one(1,w,fname, decompE)
        
        fname = InitFileDirectory(:len_trim(InitFileDirectory))//"/mystery_T.bin"
        call decomp_2d_read_one(1,T,fname, decompE)

    end select  
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
    real(rkind)  :: Lx = one, Ly = one, Lz = one
    real(rkind) :: Noise_Amp = 1.d-6 
    integer :: ProblemMode = 1
    real(rkind) :: alpha = 1.0d0, phase = 0.0d0, x0 = 2.0d0, delta = 3.0d0, Pert_Amp = 1.d-4
    namelist /RBPinstability/ Lx, Ly, Noise_Amp, Pert_Amp, ProblemMode, alpha, phase, x0, delta
     
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=RBPinstability)
    close(ioUnit)    

end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use constants,          only: zero, one
    implicit none

    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    real(rkind) :: ThetaRef, Lx, Ly 
    integer :: iounit
    real(rkind) :: Noise_Amp = 1.d-6
    integer :: ProblemMode = 1
    real(rkind) :: alpha = 1.0d0, phase = 0.0d0, x0 = 2.0d0, delta = 3.0d0, Pert_Amp = 1.d-4
    namelist /RBPinstability/ Lx, Ly, Noise_Amp, Pert_Amp, ProblemMode, alpha, phase, x0, delta
    
    Tsurf = zero; dTsurf_dt = zero; ThetaRef = one

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=RBPinstability)
    close(ioUnit)    

    ! Do nothing really since this is an unstratified simulation
end subroutine


subroutine set_Reference_Temperature(inputfile, Tref)
    use kind_parameters,    only: rkind
    implicit none 
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    real(rkind) :: Lx, Ly  
    integer :: iounit
    real(rkind) :: Noise_Amp = 1.d-6
    integer :: ProblemMode = 1
    real(rkind) :: alpha = 1.0d0, phase = 0.0d0, x0 = 2.0d0, delta = 3.0d0, Pert_Amp = 1.d-4
    
    namelist /RBPinstability/ Lx, Ly, Noise_Amp, Pert_Amp, ProblemMode, alpha, phase, x0, delta

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=RBPinstability)
    close(ioUnit)    
     
    Tref = 1.d0
    
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
