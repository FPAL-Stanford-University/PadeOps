module gabls_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    use constants, only: kappa 
    implicit none
    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.002_rkind ! 0.2% of the mean value
    integer :: nxg, nyg, nzg

    real(rkind), parameter :: xdim = 400._rkind, udim = 8._rkind, Tref = 263.5_rkind
    real(rkind), parameter :: timeDim = xdim/udim
end module     

subroutine meshgen_wallM(decomp, dx, dy, dz, mesh, inputfile)
    use gabls_parameters    
    use kind_parameters,  only: rkind
    use constants,        only: zero, one, two, three, four, pi
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    real(rkind) :: z0init
    integer :: i,j,k, ioUnit
    character(len=*),                intent(in)    :: inputfile
    integer :: ix1, ixn, iy1, iyn, iz1, izn
    real(rkind)  :: Lx = one, Ly = one, Lz = one, Tsurf0 = zero, dTsurf_dt = zero
    namelist /GABLSINPUT/ Lx, Ly, Lz, z0init, Tsurf0, dTsurf_dt 

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=GABLSINPUT)
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
    use gabls_parameters
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
    real(rkind) :: mfactor, sig
    real(rkind), dimension(:,:,:), allocatable :: randArr, ztmp, Tpurt
    real(rkind) :: z0init, epsnd
    real(rkind), dimension(:,:,:), allocatable :: ybuffC, ybuffE, zbuffC, zbuffE
    integer :: nz, nzE
    real(rkind) :: delta_Ek = 0.08d0, Xperiods = 3.d0, Yperiods = 3.d0, Zperiods = 1.d0
    real(rkind) :: zpeak = 0.2d0, Tsurf0 = zero, dTsurf_dt = zero
    real(rkind)  :: Lx = one, Ly = one, Lz = one
    namelist /GABLSINPUT/ Lx, Ly, Lz, z0init, Tsurf0, dTsurf_dt 

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=GABLSINPUT)
    close(ioUnit)    


    u  => fieldsC(:,:,:,1)
    v  => fieldsC(:,:,:,2)
    wC => fieldsC(:,:,:,3)
    w  => fieldsE(:,:,:,1)
    T  => fieldsC(:,:,:,7)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)
 
    epsnd = 1.d-2

    u = one !+ epsnd*cos(Yperiods*two*pi*y/Ly)*exp(-half*(z/zpeak/Lz)**2.d0)
    v = zero!epsnd*(z/Lz)*cos(Xperiods*two*pi*x/Lx)*exp(-half*(z/zpeak/Lz)**2.d0)
    wC= zero 

    allocate(ztmp(decompC%xsz(1),decompC%xsz(2),decompC%xsz(3)))
    allocate(Tpurt(decompC%xsz(1),decompC%xsz(2),decompC%xsz(3)))
    ztmp = z*xDim
    T = 0.01d0*(ztmp - 100.d0) + 265.d0    
    where(ztmp < 100.d0)
        T = 265.d0
    end where
    T = T + 0.0001d0*ztmp
    Tpurt = 0.1d0*cos(12.d0*two*pi*x)*sin(12.d0*two*pi*y)*sin(12.d0*two*pi*(z - 50.d0/xDim))
    where (ztmp>50)
        Tpurt = zero
    end where
    T = T + Tpurt
    deallocate(ztmp, Tpurt)
    T = T/Tref

    ! Add random numbers
    !allocate(randArr(size(u,1),size(u,2),size(u,3)))
    !call gaussian_random(randArr,zero,one,seedu + 10*nrank)
    !do k = 1,size(u,3)
    !    sig = randomScaleFact*(one/kappa)*log(z(1,1,k)/z0nd)
    !    u(:,:,k) = u(:,:,k) + sig*randArr(:,:,k)
    !end do  
    !deallocate(randArr)
    !
    !allocate(randArr(size(v,1),size(v,2),size(v,3)))
    !call gaussian_random(randArr,zero,one,seedv+ 10*nrank)
    !do k = 1,size(v,3)
    !    sig = randomScaleFact*z(1,1,k)*exp(-half*(z(1,1,k)/zpeak/Lz)**2)
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


subroutine set_planes_io(xplanes, yplanes, zplanes)
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 1, nyplanes = 1, nzplanes = 1

    allocate(xplanes(nxplanes), yplanes(nyplanes), zplanes(nzplanes))

    xplanes = [64]
    yplanes = [64]
    zplanes = [30]

end subroutine

subroutine setDirichletBC_Temp(inputfile, Tsurf, dTsurf_dt)
    use kind_parameters,    only: rkind
    use constants,          only: zero
    use gabls_parameters
    implicit none
    real(rkind), intent(out) :: Tsurf, dTsurf_dt
    character(len=*),                intent(in)    :: inputfile
    real(rkind) :: Lx, Ly, Lz, z0init = zero, Tsurf0 = zero
    integer :: ioUnit 
    namelist /GABLSINPUT/ Lx, Ly, Lz, z0init, Tsurf0, dTsurf_dt 

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=GABLSINPUT)
    close(ioUnit)    

    ! First change dTsurf_dt to K/s from K/hr
    dTsurf_dt = dTsurf_dt / 3600._rkind
    
    ! Now Normalize: 
    dTsurf_dt = dTsurf_dt * timeDim / Tref 
    Tsurf = Tsurf0 / Tref

end subroutine

