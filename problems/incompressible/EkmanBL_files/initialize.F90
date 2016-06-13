module EkmanBL_parameters

    use exits, only: message
    use kind_parameters,  only: rkind
    
    implicit none
    real(rkind)  :: Lx, Ly, Lz, G, ustar
    real(rkind)  :: dxP, dyP, dzP
    real(rkind)  :: delta_Ek

    integer :: seedu = 321341
    integer :: seedv = 423424
    integer :: seedw = 131344
    real(rkind) :: randomScaleFact = 0.01_rkind ! 5% of the mean value

    integer :: nperiods = 8
    real(rkind) :: oscScaleFact = 0.2_rkind    ! 20% of the mean value
    integer :: nxg, nyg, nzg

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
    integer :: ix1, ixn, iy1, iyn, iz1, izn

    Lx = 1.5_rkind*pi; Ly = 1.5_rkind*pi; Lz = 1._rkind

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
    use EkmanBL_parameters
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
    integer :: ioUnit, k
    real(rkind), dimension(:,:,:), pointer :: u, v, w, x, y, z
    real(rkind) :: mfactor, sig
    real(rkind), dimension(:,:,:), allocatable :: randArr
    real(rkind) :: epsfac, Uperiods, Vperiods , zpeak

    namelist /EKMANINPUT/ Ro, delta_Ek  


    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=EKMANINPUT)
    close(ioUnit)    

    u_g = one
    
    u => fieldsC(:,:,:,1)
    v => fieldsC(:,:,:,2)
    w => fieldsE(:,:,:,1)

    z => mesh(:,:,:,3)
    y => mesh(:,:,:,2)
    x => mesh(:,:,:,1)
 

    !! CAREFUL - z in mesh is cell centers but W is at edges !!

    Uperiods = 3.0; Vperiods = 3.0; zpeak = 0.3;
    epsfac = 0.5d0;

    u = one - exp(-z/delta_Ek)*cos(z/delta_Ek) &
            + half*exp(half)*(z/Lz)*cos(Uperiods*two*pi*y/Ly)*exp(-half*(z/zpeak/Lz)**2)
    v = exp(-z/delta_Ek)*sin(z/delta_Ek) + &
            + half*exp(half)*(z/Lz)*cos(Vperiods*two*pi*x/Lx)*exp(-half*(z/zpeak/Lz)**2)
    w = zero  

    ! Add random numbers
    allocate(randArr(size(u,1),size(u,2),size(u,3)))
    call gaussian_random(randArr,zero,one,seedu + 10*nrank)
    do k = 1,size(u,3)
        sig = randomScaleFact*(1 - exp(-z(1,1,k)/delta_Ek)*cos(z(1,1,k)/delta_Ek))
        u(:,:,k) = u(:,:,k) + sig*randArr(:,:,k)
    end do  
    deallocate(randArr)
    
    allocate(randArr(size(v,1),size(v,2),size(v,3)))
    call gaussian_random(randArr,zero,one,seedv+ 10*nrank)
    do k = 1,size(v,3)
        sig = G*randomScaleFact*exp(-z(1,1,k)/delta_Ek)*sin(z(1,1,k)/delta_Ek)
        v(:,:,k) = v(:,:,k) + sig*randArr(:,:,k)
    end do  
    deallocate(randArr)

    nullify(u,v,w,x,y,z)
    
    call message(0,"Velocity Field Initialized")

end subroutine


subroutine getForcing(dpdx)
    use kind_parameters, only: rkind
    use constants, only: zero  
    real(rkind), intent(out) :: dpdx

    dpdx = zero
    

end subroutine


subroutine set_planes_io(xplanes, yplanes, zplanes)
    implicit none
    integer, dimension(:), allocatable,  intent(inout) :: xplanes
    integer, dimension(:), allocatable,  intent(inout) :: yplanes
    integer, dimension(:), allocatable,  intent(inout) :: zplanes
    integer, parameter :: nxplanes = 1, nyplanes = 1, nzplanes = 6

    allocate(xplanes(nxplanes), yplanes(nyplanes), zplanes(nzplanes))

    xplanes = 1
    yplanes = 1
    zplanes = [5 , 20, 50, 100, 200, 450]

end subroutine

