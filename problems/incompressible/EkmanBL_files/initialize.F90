module EkmanBL_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    
    implicit none
    real(rkind)  :: Lx, Ly, Lz, D, G, ustar
    real(rkind)  :: nu = 1.14_rkind*(10**(-5._rkind))
    real(rkind)  :: f  = 1.454_rkind*(10**(-4._rkind))
    real(rkind)  :: Pr = 0.7_rkind
    real(rkind)  :: Ref = 400._rkind
    real(rkind)  :: ustar_by_G = 0.065_rkind
    real(rkind)  :: deltat_by_D = 13._rkind
    real(rkind)  :: dxP, dyP, dzP
    real(rkind)  :: Re_deltat, Re_tau
    real(rkind)  :: deltat

    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.02_rkind ! 5% of the mean value

    integer :: nperiods = 8
    real(rkind) :: oscScaleFact = 0.2_rkind    ! 20% of the mean value

end module     

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use EkmanBL_parameters    
    use kind_parameters,  only: rkind
    use constants,        only: two,pi
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn

    D = sqrt(two*nu/f)
    Lx = 26*D; Ly = 26*D; Lz = 24*D
    

    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    dx = Lx/real(nx,rkind)
    dy = Ly/real(ny,rkind)
    dz = Lz/real(nz,rkind)


    ! If base decomposition is in Y
    ix1 = decomp%xst(1); iy1 = decomp%xst(2); iz1 = decomp%xst(3)
    ixn = decomp%xen(1); iyn = decomp%xen(2); izn = decomp%xen(3)
    
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = Lx/real(nx,rkind)
        dy = Ly/real(ny,rkind)
        dz = Lz/real(nz,rkind)

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                    y(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 + k - 1, rkind ) * dz + dz/two
                end do
            end do
        end do

    end associate

end subroutine

subroutine initfields_stagg(decompC, decompE, dx, dy, dz, inputfile, mesh, fieldsC, fieldsE, u_g)
    use EkmanBL_parameters
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two,pi
    use gridtools,          only: alloc_buffs
    use IncompressibleGrid, only: u_index,v_index,w_index
    use random,             only: gaussian_random
    use decomp_2d,          only: decomp_info, nrank  
    implicit none
    type(decomp_info),               intent(in)    :: decompC
    type(decomp_info),               intent(in)    :: decompE
    character(len=*),                intent(in)    :: inputfile
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in), target    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsC
    real(rkind), dimension(:,:,:,:), intent(inout), target :: fieldsE
    integer :: ioUnit, runID, k
    real(rkind), dimension(:,:,:), pointer :: u, v, w, x, y, z
    real(rkind), intent(out), optional :: u_g
    logical :: useSGS 
    real(rkind) :: mfactor, sig
    real(rkind), dimension(:,:,:), allocatable :: randArr
    namelist /IINPUT/ runId, nu, Pr, useSGS, f, deltat_by_D, ustar_by_G, Ref 


    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=IINPUT)
    close(ioUnit)    

    u => fieldsC(:,:,:,1)
    v => fieldsC(:,:,:,2)
    w => fieldsE(:,:,:,1)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)
   
    ustar = f*D*deltat_by_D
    G = ustar/ustar_by_G
    deltat = D*deltat_by_D
    Re_deltat = G*deltat/nu
    Re_tau = ustar*deltat/nu 
    dxP =  dx*ustar/nu
    dyP =  dy*ustar/nu
    dzP =  dz*ustar/nu

    u_g = G
    u = G*(one - exp(-z/D)*cos(z/D))
    v = G*exp(-z/D)*sin(z/D)
    w = zero


    ! Add random numbers
    !allocate(randArr(size(u,1),size(u,2),size(u,3)))
    !call gaussian_random(randArr,zero,one,seedu + 10*nrank)
    !do k = 1,size(u,3)
    !    sig = randomScaleFact*abs(u(1,1,k))
    !    u(:,:,k) = u(:,:,k) + sig*randArr(:,:,1)
    !end do  
    !deallocate(randArr)
    !
    !allocate(randArr(size(v,1),size(v,2),size(v,3)))
    !call gaussian_random(randArr,zero,one,seedv+ 10*nrank)
    !do k = 1,size(v,3)
    !    sig = randomScaleFact*abs(v(1,1,k))
    !    v(:,:,k) = v(:,:,k) + sig*randArr(:,:,1)
    !end do  
    !deallocate(randArr)

    !allocate(randArr(size(w,1),size(w,2),size(w,3)))
    !call gaussian_random(randArr,zero,one,seedw+ 10*nrank)
    !do k = 2,size(w,3)-1
    !    sig = randomScaleFact*abs(v(1,1,k))
    !    w(:,:,k) = w(:,:,k) + sig*randArr(:,:,1)
    !end do  
    !deallocate(randArr)

    do k = 1,size(u,3)
        mfactor = u(2,2,k)
        u(:,:,k) = u(:,:,k) - mfactor*oscScaleFact*cos(nperiods*two*pi*x(:,:,k)/Lx)*sin(nperiods*two*pi*y(:,:,k)/Ly)
        mfactor = v(2,2,k)
        v(:,:,k) = v(:,:,k) + mfactor*oscScaleFact*sin(nperiods*two*pi*x(:,:,k)/Lx)*cos(nperiods*two*pi*y(:,:,k)/Ly)
    end do 

    do k = 1,size(u,3)
        mfactor = u(4,4,k)
        u(:,:,k) = u(:,:,k) + 0.5*mfactor*oscScaleFact*sin(4*nperiods*two*pi*x(:,:,k)/Lx)*cos(nperiods*two*pi*y(:,:,k)/Ly)
        mfactor = v(4,4,k)
        v(:,:,k) = v(:,:,k) - 0.5*mfactor*oscScaleFact*cos(4*nperiods*two*pi*x(:,:,k)/Lx)*sin(nperiods*two*pi*y(:,:,k)/Ly)
    end do 
    call message(0,"============================================================================")
    call message(0,"Initialized Velocity Fields (Lam. Ekman Solution + Perturbations)")
    call message(0,"Summary:")
    call message(1,"Geostrophic Velocity (x-direction):", G)
    call message(1,"Friction Velocity (ustar):", ustar)
    call message(1,"Friction Reynolds Number:", Re_tau)
    call message(1,"Grid Spacings:")
    call message(2,"dx_plus =", dxP) 
    call message(2,"dy_plus =", dyP) 
    call message(2,"dz_plus =", dzP) 
    call message(0,"============================================================================")

end subroutine



