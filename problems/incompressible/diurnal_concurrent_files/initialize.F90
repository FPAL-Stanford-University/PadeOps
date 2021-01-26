module convective_igrid_parameters

      ! TAKE CARE OF TIME NON-DIMENSIONALIZATION IN THIS MODULE

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: zero, kappa 
    implicit none
    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.002_rkind ! 0.2% of the mean value
    integer :: nxg, nyg, nzg
    
    real(rkind), parameter :: xDim = 1000._rkind, uDim = 8._rkind !sqrt(3.0_rkind**2+9._rkind**2)
    real(rkind), parameter :: timeDim = xDim/uDim

end module     

module diurnalBCsmod
    use kind_parameters, only: rkind, clen
    use basic_io, only: read_2d_ascii 
    use interpolation, only: spline, ispline, binarysearch 
    real(rkind), dimension(:), allocatable :: t_geo, G_geo, t_flux, wT_flux, a_geo, b_geo, c_geo, a_flux, b_flux, c_flux
    real(rkind), dimension(:), allocatable :: t_galpha, G_alpha, a_galpha, b_galpha, c_galpha
    real(rkind) :: G_tolerance = 0.1d0
contains

    subroutine setup_diurnalBCs(inputfile)
        character(len=*),                intent(in)    :: inputfile
        real(rkind), dimension(:,:), allocatable :: data2read
        character(len=clen) :: fname_G, fname_wtheta, fname_galpha
        integer :: ioUnit
        namelist /DIURNAL_BCS/ fname_G, fname_wtheta, fname_galpha, G_tolerance

        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=DIURNAL_BCS)
        close(ioUnit)    

        ! Geostrophic wind speed
        call read_2d_ascii(data2read, fname_G)
        allocate(t_geo(size(data2read,1)))
        allocate(G_geo(size(data2read,1)))
        allocate(a_geo(size(data2read,1)))
        allocate(b_geo(size(data2read,1)))
        allocate(c_geo(size(data2read,1)))
        t_geo = data2read(:,1) * 60.d0*60.d0/1000.d0
        G_geo = data2read(:,2)
        call spline(t_geo, G_geo, a_geo, b_geo, c_geo, size(t_geo))
        deallocate(data2read)

        ! Surface sensible heat flux <w'theta'>
        call read_2d_ascii(data2read, fname_wtheta)
        allocate(t_flux(size(data2read,1)))
        allocate(wT_flux(size(data2read,1)))
        allocate(a_flux(size(data2read,1)))
        allocate(b_flux(size(data2read,1)))
        allocate(c_flux(size(data2read,1)))
        t_flux = data2read(:,1) * 60.d0*60.d0/1000.d0
        wT_flux = data2read(:,2)
        call spline(t_flux, wT_flux, a_flux, b_flux, c_flux, size(t_flux))
        deallocate(data2read)

        ! Geostrophic wind direction in computational domain g_alpha 
        ! NOTE: This is different from the frameAngle variable which describes
        ! the geostrophic wind direction in cardinal directions
        call read_2d_ascii(data2read, fname_galpha)
        allocate(t_galpha(size(data2read,1)))
        allocate(G_alpha(size(data2read,1)))
        allocate(a_galpha(size(data2read,1)))
        allocate(b_galpha(size(data2read,1)))
        allocate(c_galpha(size(data2read,1)))
        t_galpha = data2read(:,1) * 60.d0*60.d0/1000.d0
        G_alpha = data2read(:,2)
        call spline(t_galpha, G_alpha, a_galpha, b_galpha, c_galpha, size(t_galpha))
        deallocate(data2read)


    end subroutine

    subroutine diurnalBCs_CorrectnessCheck(time, G)
        use exits, only: GracefulExit
        real(rkind), intent(in) :: time, G
        real(rkind) :: Gtrue

        ! Test function
        !call linear_interp(size(t_geo,1),t_geo,G_geo,10.8d0,Gtrue)
        !Gtrue = ispline(10.8d0, t_geo, G_geo, a_geo, b_geo, c_geo, size(t_geo))
        !Gtrue = ispline(10.8d0, t_flux, wT_flux, a_flux, b_flux, c_flux, size(t_flux))
        !call linear_interp(size(t_flux,1),t_flux,wT_flux,10.8d0,Gtrue)
 
        ! Interpolate G: 
        !Gtrue = ispline(time, t_geo, G_geo, a_geo, b_geo, c_geo, size(t_geo))
        call linear_interp(size(t_geo,1),t_geo,G_geo,time,Gtrue)
        
        if (abs(Gtrue - G) .ge. G_tolerance) then 
            print*, "time:", time
            print*, "G entered:", G
            print*, "G from BC:", Gtrue
            call gracefulExit("Incorrect Geostrophic inputs",45)
        end if 
    end subroutine     
    
    subroutine get_diurnalBCs(time, G, wTheta, Gangle)
        real(rkind), intent(in) :: time
        real(rkind), intent(out) :: G, wTheta, Gangle

        ! Convert units for time? 
        ! This assumes that the input files' time vector is appropriately
        ! nondimensionalized       
         
        ! Interpolate G: 
        ! Only works for the diurnal case based on 1 hour of unstable
        ! initialization
        ! I directly modified the RESTART file to be zero time!!!
        !G = ispline(time, t_geo, G_geo, a_geo, b_geo, c_geo, size(t_geo))
        ! Interpolate wTheta: 
        !wTheta = ispline(time, t_flux, wT_flux, a_flux, b_flux, c_flux, size(t_flux))
 
        ! Linear inteprolation
        call linear_interp(size(t_geo,1),t_geo,G_geo,time,G)
        call linear_interp(size(t_flux,1),t_flux,wT_flux,time,wTheta)
        call linear_interp(size(t_galpha,1),t_galpha,G_alpha,time,Gangle)

    end subroutine 

    
    subroutine linear_interp(xlen,x,y,xv,yv)
        use kind_parameters, only: rkind
        use decomp_2d,        only: decomp_info
        implicit none
        integer, intent(in)                   :: xlen
        real(rkind), dimension(xlen), intent(in) :: x, y
        real(rkind), intent(in)               :: xv
        real(rkind), intent(out)              :: yv
        integer                               :: ind
     
        ! binary search
        call binarysearch(xlen, x, xv, ind, 1.d-6)
    
        ! Linear interpolation
        yv = y(ind) + (xv-x(ind)) * (y(ind+1)-y(ind)) / (x(ind+1)-x(ind))
    
    end subroutine

end module 

subroutine initfields_wallM(decompC, decompE, inputfile, mesh, fieldsC, fieldsE)
    use convective_igrid_parameters
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two, pi, half, three
    use gridtools,          only: alloc_buffs
    use random,             only: gaussian_random, uniform_random
    use decomp_2d          
    use reductions,         only: p_maxval
    use constants,          only: pi
    implicit none
    type(decomp_info),               intent(in)    :: decompC
    type(decomp_info),               intent(in)    :: decompE
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    integer :: ioUnit
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, T, x, y, z
    real(rkind), dimension(:,:,:), allocatable :: ybuffC, ybuffE, zbuffC, zbuffE, ztmp
    integer :: nz, nzE, k
    real(rkind) :: sig, ebar
    real(rkind)  :: Lx = one, Ly = one, Lz = one, Tref = zero, wTh_surf0 = one, dTsurf_dt = -0.05d0, z0init = 1.d-4, Tsurf0 = 290.d0
    real(rkind), dimension(:,:,:), allocatable :: randArr1, randArr2, randArr3
    
    namelist /PROBLEM_INPUT/ Lx, Ly, Lz, Tref, wTh_surf0, z0init, Tsurf0 

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBLEM_INPUT)
    close(ioUnit)    


    !!!!!!!!!!!!!!!!!!!!! DON'T CHANGE THE POINTERS / ALLOCATIONS !!!!!!!!!!!!!!!!!!!!!!
    u  => fieldsC(:,:,:,1); v  => fieldsC(:,:,:,2); wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1); T  => fieldsC(:,:,:,7) 
    z => mesh(:,:,:,3); y => mesh(:,:,:,2); x => mesh(:,:,:,1)
    !allocate(randArr(size(T,1),size(T,2),size(T,3)))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    u = uDim
    v = zero/uDim
    wC = zero/uDim

    ! From Kumar et al. (2006)
    allocate(ztmp(decompC%xsz(1),decompC%xsz(2),decompC%xsz(3)))
    ztmp = z*xDim
    do k = 1, decompC%xsz(3)
       if (ztmp(1,1,k)<800.d0) then
           T(:,:,k) = Tsurf0
       elseif (ztmp(1,1,k)>=800.d0 .and. ztmp(1,1,k)<1000.d0) then
           T(:,:,k) = Tsurf0 + (ztmp(1,1,k)-800.d0) * 6.d0 / 200.d0
       elseif (ztmp(1,1,k)>=1000.d0) then
           T(:,:,k) = Tsurf0 + (1000.d0 - 800.d0) * 6.d0 / 200.d0 + (ztmp(1,1,k)-1000.d0) * 1.d0 / 100.d0
       end if
    end do
    !do k = 1, decompC%xsz(3)
    !  if(ztmp(1,1,k) < 800.d0) then
    !    T(:,:,k) = 301.D0
    !  elseif(ztmp(1,1,k) < 850.d0) then
    !    T(:,:,k) = 286.0d0
    !  elseif(ztmp(1,1,k) < 900.d0) then
    !    T(:,:,k) = 286.0d0 + (ztmp(1,1,k)-850.d0)/50.0d0*2.0d0
    !  elseif(ztmp(1,1,k) < 1000.d0) then
    !    T(:,:,k) = 288.0d0 + (ztmp(1,1,k)-900.d0)/100.0d0*4.0d0
    !  else
    !    T(:,:,k) = 292.0d0 + (ztmp(1,1,k)-1000.d0)/1000.d0*8.0d0
    !  endif
    !enddo

    ! Add random numbers
    allocate(randArr1(size(T,1),size(T,2),size(T,3)))
    allocate(randArr2(size(T,1),size(T,2),size(T,3)))
    allocate(randArr3(size(T,1),size(T,2),size(T,3)))
    call uniform_random(randArr1,-half,half,seedu + 10*nrank)
    call uniform_random(randArr2,-half,half,seedv + 10*nrank)
    call uniform_random(randArr3,-half,half,seedw + 10*nrank)
    do k = 1,size(u,3)
        if(ztmp(1,1,k) .le. 800.0d0) then
           ebar = 0.5d0*(1.0d0-ztmp(1,1,k)/800.d0)
        else
           ebar = 0.0d0
        endif
        sig = sqrt(two*ebar/three) * sqrt(12.0d0) ! 1/sqrt(12) is stdev of uniform random distribution of width 1.
        sig = sig/uDim
        !write(*,'(5(e19.12,1x))') ztmp(1,1,k), ebar, sig, maxval(randArr1), minval(randArr1)
         u(:,:,k) =  u(:,:,k) + sig*randArr1(:,:,k)
         v(:,:,k) =  v(:,:,k) + sig*randArr2(:,:,k)
        wC(:,:,k) = wC(:,:,k) + sig*randArr3(:,:,k)
    end do
    deallocate(randArr3)
    deallocate(randArr2)
    deallocate(randArr1)
    deallocate(ztmp)

    !!!!!!!!!!!!!!!!!!!!! DON'T CHANGE ANYTHING UNDER THIS !!!!!!!!!!!!!!!!!!!!!!
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
    ! Deallocate local memory 
    deallocate(ybuffC,ybuffE,zbuffC, zbuffE) 
    nullify(u,v,w,x,y,z)
    call message(0,"Velocity Field Initialized")
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end subroutine

subroutine setInhomogeneousNeumannBC_Temp(inputfile, wTh_surf)
    use kind_parameters,    only: rkind
    use convective_igrid_parameters
    use constants, only: one, zero 
    use diurnalBCsmod
    implicit none
    real(rkind), intent(inout) :: wTh_surf
    character(len=*),                intent(in)    :: inputfile
    real(rkind) :: tmp, tmp2, time

    time = wTh_surf ! temporary workaround

    call get_diurnalBCs(time, tmp, wTh_surf, tmp2)  ! << will be an issue for restarts

end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use convective_igrid_parameters
    use constants, only: one, zero 
    implicit none
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    character(len=*),                intent(in)    :: inputfile
    integer :: ioUnit 
    real(rkind)  :: Lx = one, Ly = one, Lz = one, Tref = zero, wTh_surf0 = one, z0init = 1.d-4, Tsurf0 = 290.0d0
    namelist /PROBLEM_INPUT/ Lx, Ly, Lz, Tref, wTh_surf0, z0init, Tsurf0 
     
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBLEM_INPUT)
    close(ioUnit)    

    ! Do nothing here
    Tsurf = Tsurf0
    dTsurf_dt = zero
end subroutine

subroutine set_planes_io(xplanes, yplanes, zplanes)
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 1, nyplanes = 1, nzplanes = 1

    allocate(xplanes(nxplanes), yplanes(nyplanes), zplanes(nzplanes))

    xplanes = [64]
    yplanes = [64]
    zplanes = [256]

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! THE SUBROUTINES UNDER THIS DON'T TYPICALLY NEED TO BE CHANGED !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use convective_igrid_parameters    
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
    real(rkind)  :: Lx = one, Ly = one, Lz = one, Tref = zero, wTh_surf0 = one, dTsurf_dt = -0.05d0, z0init = 1.d-4, Tsurf0 = 290.d0
    namelist /PROBLEM_INPUT/ Lx, Ly, Lz, Tref, wTh_surf0, z0init, Tsurf0 

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBLEM_INPUT)
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

subroutine set_Reference_Temperature(inputfile, Thetaref)
    use kind_parameters,    only: rkind
    use constants, only: one, zero 
    implicit none
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Thetaref
    integer :: ioUnit 
    real(rkind)  :: Lx = one, Ly = one, Lz = one, Tref = zero, wTh_surf0 = one, dTsurf_dt = -0.05d0, z0init = 1.d-4, Tsurf0 = 290.d0
    namelist /PROBLEM_INPUT/ Lx, Ly, Lz, Tref, wTh_surf0, z0init, Tsurf0 
     
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBLEM_INPUT)
    close(ioUnit)    

    Thetaref = Tref
    ! This will set the value of Tref.     

end subroutine


subroutine set_KS_planes_io(planesCoarseGrid, planesFineGrid)
    integer, dimension(:), allocatable,  intent(inout) :: planesFineGrid
    integer, dimension(:), allocatable,  intent(inout) :: planesCoarseGrid
    
    allocate(planesCoarseGrid(1), planesFineGrid(1))
    planesCoarseGrid = [8]
    planesFineGrid = [16]

end subroutine

!subroutine initScalar(decompC, inpDirectory, mesh, scalar_id, scalarField)
!    use kind_parameters, only: rkind
!    use decomp_2d,        only: decomp_info
!    type(decomp_info),                                          intent(in)    :: decompC
!    character(len=*),                intent(in)    :: inpDirectory
!    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
!    integer, intent(in)                            :: scalar_id
!    real(rkind), dimension(:,:,:), intent(out)     :: scalarField
!
!    scalarField = 0.d0
!end subroutine 
!
!subroutine setScalar_source(decompC, inpDirectory, mesh, scalar_id, scalarSource)
!    use kind_parameters, only: rkind
!    use decomp_2d,        only: decomp_info
!    type(decomp_info),                                          intent(in)    :: decompC
!    character(len=*),                intent(in)    :: inpDirectory
!    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
!    integer, intent(in)                            :: scalar_id
!    real(rkind), dimension(:,:,:), intent(out)     :: scalarSource
!
!    scalarSource = 0.d0
!end subroutine 
!


module link_turbine_to_scalar
    use kind_parameters, only: rkind, clen 
    implicit none 
    complex(rkind), dimension(:,:,:), allocatable :: rhs_x, rhs_y, rhs_z
    real(rkind), dimension(:,:,:), allocatable :: utmp, vtmp, wtmp, scalar_source_0
    character(clen) :: OutputDir 
contains 

    subroutine setup_turb_scalar_source(primary)
        use decomp_2d_io
        use IncompressibleGrid, only: igrid
        use Exits, only : message
        class(igrid), intent(inout), target :: primary
        real(rkind), dimension(:,:,:), allocatable :: scalar_source_0, scalar_source_1
        real(rkind), dimension(10000) :: inst_horz_avg
        character(len=clen) :: fname

        
        allocate(utmp(primary%gpC%xsz(1),primary%gpC%xsz(2),primary%gpC%xsz(3)))
        allocate(vtmp(primary%gpC%xsz(1),primary%gpC%xsz(2),primary%gpC%xsz(3)))
        allocate(wtmp(primary%gpE%xsz(1),primary%gpE%xsz(2),primary%gpE%xsz(3)))
        
        allocate(scalar_source_0(primary%gpC%xsz(1),primary%gpC%xsz(2),primary%gpC%xsz(3)))
        allocate(scalar_source_1(primary%gpC%xsz(1),primary%gpC%xsz(2),primary%gpC%xsz(3)))
        
        allocate(rhs_y(primary%sp_gpC%ysz(1), primary%sp_gpC%ysz(2), primary%sp_gpC%ysz(3))) 
        allocate(rhs_z(primary%sp_gpE%ysz(1), primary%sp_gpE%ysz(2), primary%sp_gpE%ysz(3))) 
       
        utmp = 1.d0 
        vtmp = 0.d0 
        wtmp = 0.d0
        ! Scalar 3, primary%scalars(3)%source_hat contains spatial flags for the
        ! points of interest 
        call primary%spectC%ifft(primary%scalars(3)%source_hat,scalar_source_0)
        call primary%WindTurbineArr%getForceRHS( 1.d0, utmp, vtmp, wtmp, primary%scalars(1)%source_hat, rhs_y, rhs_z, .true., inst_horz_avg)
        call primary%spectC%ifft(primary%scalars(1)%source_hat,scalar_source_1)
        scalar_source_0 = -scalar_source_0 * scalar_source_1 / 1.0d4
        call primary%spectC%fft(scalar_source_0, primary%scalars(3)%source_hat)
        ! Scalar 4, primary%scalars(4)%source_hat contains spatial flags for the
        ! points of interest
        call primary%spectC%ifft(primary%scalars(4)%source_hat,scalar_source_0)
        call primary%WindTurbineArr%getForceRHS( 1.d0, utmp, vtmp, wtmp, primary%scalars(1)%source_hat, rhs_y, rhs_z, .true., inst_horz_avg)
        call primary%spectC%ifft(primary%scalars(1)%source_hat,scalar_source_1)
        scalar_source_0 = -scalar_source_0 * scalar_source_1 / 1.0d4
        call primary%spectC%fft(scalar_source_0, primary%scalars(4)%source_hat)
        ! Scalar 1 
        call primary%WindTurbineArr%getForceRHS( 1.d0, utmp, vtmp, wtmp, primary%scalars(1)%source_hat, rhs_y, rhs_z, .true., inst_horz_avg)
        primary%scalars(1)%source_hat = -primary%scalars(1)%source_hat / 1.0d4      
  
        call primary%spectC%ifft(primary%scalars(1)%source_hat,scalar_source_0)
        fname = primary%OutputDir(:len_trim(primary%OutputDir))//"/ScalarSource1.out"
        call decomp_2d_write_one(1,scalar_source_0,fname,primary%gpC)
        
        call primary%spectC%ifft(primary%scalars(2)%source_hat,scalar_source_0)
        fname = primary%OutputDir(:len_trim(primary%OutputDir))//"/ScalarSource2.out"
        call decomp_2d_write_one(1,scalar_source_0,fname,primary%gpC)
        
        call primary%spectC%ifft(primary%scalars(3)%source_hat,scalar_source_0)
        fname = primary%OutputDir(:len_trim(primary%OutputDir))//"/ScalarSource3.out"
        call decomp_2d_write_one(1,scalar_source_0,fname,primary%gpC)

        call primary%spectC%ifft(primary%scalars(4)%source_hat,scalar_source_0)
        fname = primary%OutputDir(:len_trim(primary%OutputDir))//"/ScalarSource4.out"
        call decomp_2d_write_one(1,scalar_source_0,fname,primary%gpC)
        
        call message(0,"SCALAR SOURCES WRITTEN TO DISK.")
    
        deallocate(utmp, vtmp, wtmp, scalar_source_0, rhs_y, rhs_z, scalar_source_1)

    end subroutine 

end module 


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
    use link_turbine_to_scalar, only: scalar_source_0

    type(decomp_info),                                          intent(in)    :: decompC
    character(len=*),                intent(in)    :: inpDirectory
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    integer, intent(in)                            :: scalar_id
    real(rkind), dimension(:,:,:), intent(out)     :: scalarSource
    real(rkind), dimension(:,:,:), pointer :: x, y, z
    real(rkind) :: dz, dy, dx, Lx, Ly

    z => mesh(:,:,:,3); y => mesh(:,:,:,2); x => mesh(:,:,:,1)
    dz = z(1,1,2) - z(1,1,1)
    dy = y(1,2,1) - y(1,1,1)
    dx = x(2,1,1) - x(1,1,1)
    Lx = 75.d0
    Ly = 35.d0
    select case (scalar_id)
    case (1) ! Turbine case
        scalarSource = 0.d0 ! Need to handle this case using init_turb2scalar_linker call 
    case (2) ! Above the turbine rows
        scalarSource = 0.d0
        scalarSource = exp(-(z-(0.25d0+1.5d0*126.d0/400.d0))**2 / 0.01) / 1.0d4 ! Mike implement this 
        do k = 1,size(z,3)
            do j = 1, size(z,2)
                do i = 1, size(z,1)
                    if ((x(i,j,k)>0.84*Lx) .or. (y(i,j,k)>0.84*Ly)) then
                        scalarSource(i,j,k) = 0.d0
                    end if
                end do
            end do
        end do
    case (3) ! Center turbinethis%blanks = 1.d0
        scalarSource = z*0.d0
        scalarSource = 1.d0 ! Need to handle this case using init_turb2scalar_linker call 
        do k = 1,size(z,3)
            do j = 1, size(y,2)
                do i = 1,size(x,1)
                    if ((x(i,j,k) < 18.d0) .or. (x(i,j,k)>19.d0)) then
                        scalarSource(i,j,k) = 0.d0
                    end if
                    if ((y(i,j,k) < 9.d0) .or. (y(i,j,k)>10.d0)) then
                        scalarSource(i,j,k) = 0.d0
                    end if
                end do
            end do
        end do
    case (4) ! Last turbine
        scalarSource = z*0.d0
        scalarSource = 1.d0 ! Need to handle this case using init_turb2scalar_linker call 
        do k = 1,size(z,3)
            do j = 1, size(y,2)
                do i = 1,size(x,1)
                    if (x(i,j,k) < 22.d0) then
                        scalarSource(i,j,k) = 0.d0
                    end if
                    if (y(i,j,k) < 11.d0) then
                        scalarSource(i,j,k) = 0.d0
                    end if
                end do
            end do
        end do
    !case (XXX) ! Turbine case
    !    scalarSource = 0.d0 ! Need to handle this case using init_turb2scalar_linker call 
    !    scalarSource = 0.d0 ! Need to handle this case using init_turb2scalar_linker call 
    !    scalarSource = z*0.d0
    !    scalarSource(1,:,:) = exp(-(z(1,:,:)-(dz))**2 / 0.01) ! Mike implement this 
    !    scalarSource(:,1,:) = exp(-(z(:,1,:)-(dz))**2 / 0.01) ! Mike implement this 
    end select

end subroutine




