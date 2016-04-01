module channel_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    
    implicit none
    real(rkind)  :: Lx, Ly, Lz, G, ustar

    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.000000005_rkind ! 5% of the mean value

    integer :: nperiods = 8
    real(rkind) :: oscScaleFact = 0.2_rkind    ! 20% of the mean value
    integer :: nxg, nyg, nzg

end module     

subroutine meshgen(decomp, dx, dy, dz, mesh)
    use channel_parameters    
    use kind_parameters,  only: rkind
    use constants,        only: two,pi,four,three
    use decomp_2d,        only: decomp_info
    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    integer :: i,j,k
    integer :: ix1, ixn, iy1, iyn, iz1, izn

    Lx = four*pi; Ly = four/three*pi; Lz = 2._rkind

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

subroutine getForcing(inputfile,dpdx)
    use constants, only: two
    use kind_parameters, only: rkind
    character(len=*), intent(in) :: inputfile
    real(rkind), intent(out) :: dpdx
    real(rkind) :: dpdxForcing = two/500._rkind
    integer :: ioUnit

    namelist /CHANNELINPUT/ dpdxForcing
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=CHANNELINPUT)
    close(ioUnit)

    dpdx = dpdxForcing 
    

end subroutine


subroutine initfields_stagg(decompC, decompE, dx, dy, dz, inputfile, mesh, fieldsC, fieldsE, u_g, Ro)
    use channel_parameters
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two,pi, half
    use gridtools,          only: alloc_buffs
    use IncompressibleGrid, only: u_index,v_index,w_index
    use random,             only: gaussian_random
    use decomp_2d,          only: decomp_info, nrank  
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
    integer :: ioUnit, i, j, k
    real(rkind), dimension(:,:,:), pointer :: u, v, w, x, y, z
    logical :: useSGS 
    real(rkind) :: mfactor, sig
    real(rkind), dimension(:,:,:), allocatable :: randArr
    real(rkind) :: epsfac, Uperiods, Vperiods , zpeak


    u_g = zero
    
    u => fieldsC(:,:,:,1)
    v => fieldsC(:,:,:,2)
    w => fieldsE(:,:,:,1)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)
 

    Uperiods = 4.0; Vperiods = 4.0; zpeak = 0.3;
    epsfac = 0.1d0;

    do k = 1,size(u,3)
        do j = 1,size(u,2)
            do i = 1,size(u,1)
                !u(i,j,k) = (1.0d0-exp(-z(i,j,k)*2.5d0)*cos(z(i,j,k)*2.50d0)+epsfac*exp(0.5d0)*(z(i,j,k)/Lz) &
                !        *cos(Uperiods*2.0d0*pi*y(i,j,k)/Ly)*exp(-0.5e0*(z(i,j,k)/zpeak/Lz)**2.0d0))
                !v(i,j,k) = (exp(-z(i,j,k)*2.5d0)*sin(z(i,j,k)*2.50d0)+epsfac*exp(0.5d0)*(z(i,j,k)/Lz)& 
                !            *cos(Vperiods*2.0d0*pi*x(i,j,k)/Lx)*exp(-0.5d0*(z(i,j,k)/zpeak/Lz)**2.0d0))
                !w(i,j,k) = zero  
                u(i,j,k) = (1.0d0-(z(i,j,k)-1.0d0)**8) + epsfac*0.5d0*Lx*dsin(pi*(z(i,j,k)-1.0d0))*dcos(2.0d0*two*pi/Lx*x(i,j,k))*dsin(two*pi/Ly*y(i,j,k))
                v(i,j,k) = -epsfac*0.5d0*Ly*dsin(pi*(z(i,j,k)-1.0d0))*dsin(2.0d0*two*pi/Lx*x(i,j,k))*dcos(two*pi/Ly*y(i,j,k))
                w(i,j,k) = -epsfac*(1.0d0+dcos(pi*(z(i,j,k)-1.0d0-dz/two)))*dsin(2.0d0*two*pi/Lx*x(i,j,k))*dsin(two*pi/Ly*y(i,j,k))

            end do 
        end do 
    end do 
   
    !u = half*z*(Lz - z) !one - exp(-z/delta_Ek)*cos(z/delta_Ek) &
    !    !    + half*exp(half)*(z/Lz)*cos(Uperiods*two*pi*y/Ly)*exp(-half*(z/zpeak/Lz)**2)
    !v = zero!exp(-z/delta_Ek)*sin(z/delta_Ek) + &
    !    !    + half*exp(half)*(z/Lz)*cos(Vperiods*two*pi*x/Lx)*exp(-half*(z/zpeak/Lz)**2)
    !!w = zero  
    !w = half*half*z*(Lz - z) 

    ! Add random numbers
    !allocate(randArr(size(u,1),size(u,2),size(u,3)))
    !call gaussian_random(randArr,zero,one,seedu + 10*nrank)
    !do k = 1,size(u,3)
    !    sig = G*randomScaleFact*(1 - exp(-z(1,1,k)/D)*cos(z(1,1,k)/D))
    !    u(:,:,k) = u(:,:,k) + sig*randArr(:,:,k)*0.5*(1 - cos(2*pi*y(:,:,k)/Ly))*0.5*(1 - cos(2*pi*x(:,:,k)/Lx))
    !end do  
    !deallocate(randArr)
    !
    !allocate(randArr(size(v,1),size(v,2),size(v,3)))
    !call gaussian_random(randArr,zero,one,seedv+ 10*nrank)
    !do k = 1,size(v,3)
    !    sig = G*randomScaleFact*exp(-z(1,1,k)/D)*sin(z(1,1,k)/D)
    !    v(:,:,k) = v(:,:,k) + sig*randArr(:,:,k)*0.5*(1 - cos(2*pi*y(:,:,k)/Ly))*0.5*(1 - cos(2*pi*x(:,:,k)/Lx))
    !end do  
    !deallocate(randArr)
    !write(*,*) '--x--', maxval(x), minval(x)
    !write(*,*) '--y--', maxval(y), minval(y)
    !write(*,*) '--z--', maxval(z), minval(z)
    !write(*,*) '--u--', maxval(u), minval(u)
    !write(*,*) '--v--', maxval(v), minval(v)
    !write(*,*) '--w--', maxval(w), minval(w)

    nullify(u,v,w,x,y,z)
    
    call message(0,"Velocity Field Initialized")

end subroutine




