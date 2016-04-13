module pbl_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    
    implicit none
    real(rkind)  :: Lx, Ly, Lz, G = 15._rkind
    real(rkind)  :: ustar = 0.45_rkind, H = 1000._rkind, z0 = 0.1_rkind
    real(rkind)  :: f = 1.45d-4
    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind), parameter :: kappa = 0.41
    real(rkind) :: randomScaleFact = 0.01_rkind ! 5% of the mean value
    integer :: nxg, nyg, nzg

end module     

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use pbl_parameters    
    use kind_parameters,  only: rkind
    use constants,        only: one, two, three, four, pi
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k
    integer :: ix1, ixn, iy1, iyn, iz1, izn

    Lx = three; Ly = three; Lz = one

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

subroutine initfields_stagg(decompC, decompE, dx, dy, dz, inputfile, mesh, fieldsC, fieldsE, u_g, Ro)
    use pbl_parameters
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
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    real(rkind), intent(out), optional :: Ro, u_g
    integer :: ioUnit, k
    real(rkind), dimension(:,:,:), pointer :: u, v, w, wC, x, y, z
    real(rkind) :: mfactor, sig
    real(rkind), dimension(:,:,:), allocatable :: randArr
    real(rkind) :: z0nd, epsnd = 0.1
    real(rkind), dimension(:,:,:), allocatable :: ybuffC, ybuffE, zbuffC, zbuffE
    integer :: nz, nzE
    real(rkind) :: delta_Ek = 0.08, Uperiods = 3, Vperiods = 3 
    real(rkind) :: zpeak = 0.3
    namelist /PBLINPUT/ H, G, f, z0 


    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PBLINPUT)
    close(ioUnit)    

    ! No Coriolis right now
    u_g = one
    Ro = G/(H*f)
    
    z0nd = z0/H

    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,2)
    w  => fieldsE(:,:,:,1)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)
 
    ustar = kappa/log(Lz/z0nd);
    u = (ustar/kappa)*log(z/z0nd) + epsnd*cos(two*pi*x/Lx)*sin(two*pi*y/Ly)*cos(two*pi*z/Lz)
    v = epsnd*sin(two*pi*x/Lx)*cos(two*pi*y/Ly)*cos(two*pi*z/Lz)
    wC= epsnd*sin(two*pi*x/Lx)*sin(two*pi*y/Ly)*sin(two*pi*z/Lz)
     
    !u  = ustar/kappa * log(z/z0)
    !u = u + 1.*exp(0.5)*(z/Lz)*cos(Uperiods*2*pi*y/Ly)*exp(-0.5*(z/zpeak/Lz)**2);
    !v = exp(0.5)*(z/Lz)*cos(Vperiods*2*pi*x/Lx)*exp(-0.5*(z/zpeak/Lz)**2)

    !print*, maxval(v)
    !wC = zero
    ! Add random numbers
    !allocate(randArr(size(u,1),size(u,2),size(u,3)))
    !call gaussian_random(randArr,zero,one,seedu + 10*nrank)
    !do k = 1,size(u,3)
    !    sig = randomScaleFact*(1 - exp(-z(1,1,k)/delta_Ek)*cos(z(1,1,k)/delta_Ek))
    !    u(:,:,k) = u(:,:,k) + sig*randArr(:,:,k)
    !end do  
    !deallocate(randArr)
    !
    !allocate(randArr(size(v,1),size(v,2),size(v,3)))
    !call gaussian_random(randArr,zero,one,seedv+ 10*nrank)
    !do k = 1,size(v,3)
    !    sig = G*randomScaleFact*exp(-z(1,1,k)/delta_Ek)*sin(z(1,1,k)/delta_Ek)
    !    v(:,:,k) = v(:,:,k) + sig*randArr(:,:,k)
    !end do  
    !deallocate(randArr)


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


subroutine getForcing(inputfile, dpdx)
    use kind_parameters, only: rkind
    use constants, only: zero, one 
    use pbl_parameters    
    implicit none
    character(len=*), intent(in) :: inputfile 
    real(rkind), intent(out) :: dpdx
    integer :: ioUnit
    
    dpdx = zero
    

end subroutine

subroutine set_planes_io(xplanes, yplanes, zplanes)
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 2, nyplanes = 2, nzplanes = 2

    allocate(xplanes(nxplanes), yplanes(nyplanes), zplanes(nzplanes))

    xplanes = [2 , 10]
    yplanes = [2 , 10]
    zplanes = [2 , 10]

end subroutine


