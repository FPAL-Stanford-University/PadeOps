module randomJetHIT_parameters
    use exits, only: message
    use kind_parameters,  only: rkind, clen
    use constants, only: two, pi
    use fortran_assert, only: assert
    implicit none
    
    logical :: filterMask = .false.
    integer :: NjetsX = 0, NjetsY = 0
    real(rkind) :: jetXsz, jetYsz, jetZsz
    real(rkind) :: jetZst = 0.d0
    real(rkind) :: wJet = 1.d0
    
    ! GlobalDomainStuff
    real(rkind) :: Lx = two*pi, Ly = two*pi, Lz = two*pi
    integer :: nx, ny, nz

    ! IO stuff
    integer :: dumpMaskFreq = 100

    contains
      
end module     

subroutine initfields_wallM(decompC, decompE, inputfile, mesh, fieldsC, fieldsE)
    use randomJetHIT_parameters
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
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, T, x, y, z
    


    !!!!!!!!!!!!!!!!!!!!! DON'T CHANGE THE POINTERS / ALLOCATIONS !!!!!!!!!!!!!!!!!!!!!!
    u  => fieldsC(:,:,:,1); v  => fieldsC(:,:,:,2); wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1); T  => fieldsC(:,:,:,7) 
    z => mesh(:,:,:,3); y => mesh(:,:,:,2); x => mesh(:,:,:,1)
    !allocate(randArr(size(T,1),size(T,2),size(T,3)))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    u  = zero 
    v  = zero 
    wC = zero
    w  = zero 

    !!!!!!!!!!!!!!!!!!!!! DON'T CHANGE ANYTHING UNDER THIS !!!!!!!!!!!!!!!!!!!!!!
    nullify(u,v,w,x,y,z)
    call message(0,"Velocity Field Initialized")
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end subroutine

subroutine setInhomogeneousNeumannBC_Temp(inputfile, wTh_surf)
    use kind_parameters,    only: rkind
    use randomJetHIT_parameters
    use constants, only: one, zero 
    implicit none
    real(rkind), intent(inout) :: wTh_surf
    character(len=*),                intent(in)    :: inputfile
    real(rkind) :: tmp, tmp2, time

    wTh_surf = 0.d0 ! temporary workaround


end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use randomJetHIT_parameters
    use constants, only: one, zero 
    implicit none
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    character(len=*),                intent(in)    :: inputfile
    integer :: ioUnit 
    real(rkind)  :: Tref = zero, wTh_surf0 = one, z0init = 1.d-4, Tsurf0 = 290.0d0
    namelist /PROBLEM_INPUT/ Lx, Ly, Lz
     
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBLEM_INPUT)
    close(ioUnit)    

    ! Do nothing here
    Tsurf = Tsurf0
    dTsurf_dt = zero
end subroutine

subroutine set_planes_io(xplanes, yplanes, zplanes)
    use randomJetHIT_parameters, only: nx, ny, nz
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 1, nyplanes = 1, nzplanes = 3

    allocate(xplanes(nxplanes), yplanes(nyplanes), zplanes(nzplanes))

    xplanes = [nx/2]
    yplanes = [ny/2]
    zplanes = [12, 23, 64]

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
    use randomJetHIT_parameters    
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
    integer :: nxg, nyg, nzg
    
    namelist /PROBLEM_INPUT/ Lx, Ly, Lz, filterMask, dumpMaskFreq, &
      NjetsX, NjetsY, jetZst, jetXsz, jetYsz, jetZsz, wJet

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBLEM_INPUT)
    close(ioUnit)    

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

    nx = decomp%xsz(1)
    ny = decomp%ysz(2)
    nz = decomp%zsz(3)
end subroutine

subroutine set_Reference_Temperature(inputfile, Thetaref)
    use kind_parameters,    only: rkind
    use constants, only: one, zero 
    implicit none
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Thetaref
    integer :: ioUnit 
    real(rkind)  :: Lx = one, Ly = one, Lz = one, Tref = zero, wTh_surf0 = one, dTsurf_dt = -0.05d0, z0init = 1.d-4, Tsurf0 = 290.d0
    namelist /PROBLEM_INPUT/ Lx, Ly, Lz
     
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
subroutine hook_source(tsim,mesh,Re,urhs,vrhs,wrhs)
    use kind_parameters, only: rkind
    real(rkind),                     intent(in)    :: tsim, Re
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:),   intent(inout) :: urhs, vrhs, wrhs
    
    urhs = urhs + 0.d0
    vrhs = vrhs + 0.d0
    wrhs = wrhs + 0.d0
  end subroutine
