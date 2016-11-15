module ekmanNeutral_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa 
    implicit none
    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.002_rkind ! 0.2% of the mean value
    integer :: nxg, nyg, nzg

    real(rkind), parameter :: xdim = 400._rkind, udim = 8._rkind
    real(rkind), parameter :: timeDim = xdim/udim
end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use ekmanNeutral_parameters    
    use kind_parameters,  only: rkind
    use constants,        only: zero, one, two, three, four, pi
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    real(rkind) :: z0init = 1.d-4
    integer :: i,j,k, ioUnit
    character(len=*),                intent(in)    :: inputfile
    integer :: ix1, ixn, iy1, iyn, iz1, izn
    real(rkind)  :: Lx = one, Ly = one, Lz = one, Tref = zero
    namelist /EKMAN_NEUTRAL_INPUT/ Lx, Ly, Lz, z0init, Tref 

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=EKMAN_NEUTRAL_INPUT)
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
    use ekmanNeutral_parameters
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two, three, four, pi, half
    use gridtools,          only: alloc_buffs
    use IncompressibleGrid, only: u_index,v_index,w_index
    use random,             only: gaussian_random
    use decomp_2d          
    use reductions,         only: p_maxval
    implicit none
    type(decomp_info),               intent(in)    :: decompC
    type(decomp_info),               intent(in)    :: decompE
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    integer :: ioUnit, k
    real(rkind), dimension(:,:,:), pointer :: u, v, w, T, wC, x, y, z
    real(rkind), dimension(:,:,:), allocatable :: randArr, Tpurt
    real(rkind) :: z0init, sig!, epsnd
    real(rkind), dimension(:,:,:), allocatable :: ybuffC, ybuffE, zbuffC, zbuffE, eta
    integer :: nz, nzE
    !real(rkind) :: Xperiods = 3.d0, Yperiods = 3.d0, Zperiods = 1.d0
    real(rkind)  :: Lx = one, Ly = one, Lz = one, Tref = zero
    real(rkind), parameter :: a = 2.43_rkind, b = 0.027_rkind, thetam = 15._rkind + 273.15_rkind
    real(rkind), parameter :: c = 1.d0/3.d0, h0 = 0.41667d0, h1 = 0.45833d0, h2 = 0.5d0, D = h0/4.d0

    namelist /EKMAN_NEUTRAL_INPUT/ Lx, Ly, Lz, z0init, Tref 

    !z0init = 1.D-4
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=EKMAN_NEUTRAL_INPUT)
    close(ioUnit)    


    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1)
    T  => fieldsC(:,:,:,7)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)
 
    u = 1.d0!(1.d0 - exp(-z/D)*cos(z/D))
    v = zero!exp(-z/D)*sin(z/D)
    wC = zero 

    allocate(Tpurt(decompC%xsz(1),decompC%xsz(2),decompC%xsz(3)))
    allocate(eta(decompC%xsz(1),decompC%xsz(2),decompC%xsz(3)))
    eta = (z - h1)/(c*(h2 - h0))
    T = thetam + a*(tanh(eta) + 1.d0)/2.d0 + b*(log(2.d0*cosh(eta)) + eta)/2.d0 

    ! Add random numbers
    allocate(randArr(size(T,1),size(T,2),size(T,3)))
    call gaussian_random(randArr,zero,one,seedu + 10*nrank)
    do k = 1,size(u,3)
        sig = 0.08d0
        Tpurt(:,:,k) = sig*randArr(:,:,k)
    end do  
    deallocate(randArr)
    
    where (z > 0.25d0)
        Tpurt = zero
    end where
    T = T + Tpurt
   
    deallocate(eta, Tpurt)

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
    integer, parameter :: nxplanes = 1, nyplanes = 1, nzplanes = 1

    allocate(xplanes(nxplanes), yplanes(nyplanes), zplanes(nzplanes))

    xplanes = [128]
    yplanes = [128]
    zplanes = [20]

end subroutine

subroutine set_Reference_Temperatur(inputfile, Tref)
    use kind_parameters,    only: rkind
    use constants,          only: zero
    implicit none
    character(len=*),                intent(in)    :: inputfile
    real(rkind), intent(out) :: Tref
    real(rkind) :: Lx, Ly, Lz, z0init
    integer :: ioUnit 
    namelist /EKMAN_NEUTRAL_INPUT/ Lx, Ly, Lz, z0init, Tref 
     
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=EKMAN_NEUTRAL_INPUT)
    close(ioUnit)    

    ! This will have set the value of Tref.     

end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use constants,          only: zero
    use ekmanNeutral_parameters
    implicit none
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    real(rkind) :: Lx, Ly, Lz, z0init, Tref
    character(len=*),                intent(in)    :: inputfile
    integer :: ioUnit 
    namelist /EKMAN_NEUTRAL_INPUT/ Lx, Ly, Lz, z0init, Tref 
     
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=EKMAN_NEUTRAL_INPUT)
    close(ioUnit)    
    
    ! This subroutine is never called for this problem, since it used Neumann
    ! BC.

    Tsurf = 0.d0; dTsurf_dt = 0.d0
end subroutine

subroutine set_KS_planes_io(planesCoarseGrid, planesFineGrid)
    integer, dimension(:), allocatable,  intent(inout) :: planesFineGrid
    integer, dimension(:), allocatable,  intent(inout) :: planesCoarseGrid
    
    allocate(planesCoarseGrid(1), planesFineGrid(1))
    planesCoarseGrid = [8]
    planesFineGrid = [16]

end subroutine
