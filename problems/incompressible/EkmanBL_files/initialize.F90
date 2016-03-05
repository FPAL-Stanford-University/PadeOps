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
    real(rkind) :: randomScaleFact = 0.000000005_rkind ! 5% of the mean value

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

    D = sqrt(two*nu/f)
    Lx = 150*D; Ly = 75*D; Lz = 24*D
    

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

    end associate

end subroutine

subroutine initfields_stagg(decompC, decompE, dx, dy, dz, inputfile, mesh, fieldsC, fieldsE, u_g, f_corr)
    use EkmanBL_parameters
    use kind_parameters,    only: rkind
    use constants,          only: zero, one, two,pi
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
    real(rkind), intent(out), optional :: f_corr, u_g
    integer :: ioUnit, runID, k
    real(rkind), dimension(:,:,:), pointer :: u, v, w, x, y, z
    logical :: useSGS 
    real(rkind) :: mfactor, sig
    real(rkind), dimension(:,:,:), allocatable :: randArr
    real(rkind) :: epsfac, Uperiods, Vperiods , zpeak

    namelist /EKMANINPUT/ Ref, deltat_by_D, & 
                                ustar_by_G, f 


    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=EKMANINPUT)
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
    f_corr = f
    !u = G*(one - exp(-z/D)*cos(z/D))
    !v = G*exp(-z/D)*sin(z/D)
    !w = zero

    Uperiods = 4.0; Vperiods = 4.0; zpeak = D;
    epsfac = 0.5d0;

    u = 1.2*(1.0-exp(-z*2.5d0)*cos(z*2.50d0)+epsfac*exp(0.5d0)*(z/Lz)*cos(Uperiods*2.0d0**pi*y/Ly)*exp(-0.5e0*(z/zpeak/Lz)**2.0d0))
    v = 1.2*(exp(-z*2.5d0)*sin(z*2.50d0)+epsfac*exp(0.5d0)*(z/Lz)*cos(Vperiods*2.0d0**pi*x/Lx)*exp(-0.5d0*(z/zpeak/Lz)**2.0d0))
    w = zero

    u = u*G
    v = v*G

    ! Add random numbers
    allocate(randArr(size(u,1),size(u,2),size(u,3)))
    call gaussian_random(randArr,zero,one,seedu + 10*nrank)
    do k = 1,size(u,3)
        sig = G*randomScaleFact*(1 - exp(-z(1,1,k)/D)*cos(z(1,1,k)/D))
        u(:,:,k) = u(:,:,k) + sig*randArr(:,:,k)*0.5*(1 - cos(2*pi*y(:,:,k)/Ly))*0.5*(1 - cos(2*pi*x(:,:,k)/Lx))
    end do  
    deallocate(randArr)
    
    allocate(randArr(size(v,1),size(v,2),size(v,3)))
    call gaussian_random(randArr,zero,one,seedv+ 10*nrank)
    do k = 1,size(v,3)
        sig = G*randomScaleFact*exp(-z(1,1,k)/D)*sin(z(1,1,k)/D)
        v(:,:,k) = v(:,:,k) + sig*randArr(:,:,k)*0.5*(1 - cos(2*pi*y(:,:,k)/Ly))*0.5*(1 - cos(2*pi*x(:,:,k)/Lx))
    end do  
    deallocate(randArr)

    nullify(u,v,w,x,y,z)
    !allocate(randArr(size(w,1),size(w,2),size(w,3)))
    !call gaussian_random(randArr,zero,one,seedw+ 10*nrank)
    !do k = 2,size(w,3)-1
    !    sig = randomScaleFact*abs(v(1,1,k))
    !    w(:,:,k) = w(:,:,k) + sig*randArr(:,:,k)
    !end do  
    !deallocate(randArr)

    !do k = 1,size(u,3)
    !    mfactor = u(2,2,k)
    !    u(:,:,k) = u(:,:,k) - mfactor*oscScaleFact*cos(nperiods*two*pi*x(:,:,k)/Lx)*sin(nperiods*two*pi*y(:,:,k)/Ly)
    !    mfactor = v(2,2,k)
    !    v(:,:,k) = v(:,:,k) + mfactor*oscScaleFact*sin(nperiods*two*pi*x(:,:,k)/Lx)*cos(nperiods*two*pi*y(:,:,k)/Ly)
    !end do 

    !do k = 1,size(u,3)
    !    mfactor = u(4,4,k)
    !    u(:,:,k) = u(:,:,k) + 0.5*mfactor*oscScaleFact*sin(4*nperiods*two*pi*x(:,:,k)/Lx)*cos(nperiods*two*pi*y(:,:,k)/Ly)
    !    mfactor = v(4,4,k)
    !    v(:,:,k) = v(:,:,k) - 0.5*mfactor*oscScaleFact*cos(4*nperiods*two*pi*x(:,:,k)/Lx)*sin(nperiods*two*pi*y(:,:,k)/Ly)
    !end do 
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















module allStatistics
    use kind_parameters, only: rkind
    use decomp_2d
    use IncompressibleGridNP, only: igrid  
    implicit none

    real(rkind), dimension(:,:), allocatable, target :: zStats2dump, runningSum, TemporalMnNOW 
    real(rkind), dimension(:), pointer :: u_mean, v_mean, w_mean, uu_mean, vv_mean, ww_mean, uw_mean 
    real(rkind), dimension(:), pointer :: oxox_mean, oyoy_mean, ozoz_mean
    integer :: nxg, nyg, nzg
    integer :: tidSUM  
contains

    subroutine init_stats( ig)
        type(igrid), intent(in), target :: ig
        type(decomp_info), pointer  :: gpC

        gpC => ig%gpC
        nzg = gpC%zsz(3)
        nyg = gpC%ysz(2)
        nxg = gpC%xsz(1)
        tidSUM = 0

        allocate(zStats2dump(nzg,10))
        allocate(runningSum(nzg,10))
        allocate(TemporalMnNOW(nzg,10))
        u_mean => zStats2dump(:,1); v_mean => zStats2dump(:,2); w_mean => zStats2dump(:,3); 
        uu_mean => zStats2dump(:,4); vv_mean => zStats2dump(:,5); ww_mean => zStats2dump(:,6); 
        uw_mean => zStats2dump(:,7); oxox_mean => zStats2dump(:,8); oyoy_mean => zStats2dump(:,9); 
        ozoz_mean => zStats2dump(:,10)

        nullify(gpC)
    end subroutine



    subroutine dump_stats(ig)
        use basic_io, only: write_2d_ascii
        use kind_parameters, only: clen
        use mpi
        use exits, only: message
        type(igrid), intent(inout), target :: ig
        type(decomp_info), pointer :: gpC
        real(rkind), dimension(:,:,:), pointer :: rbuff1, rbuff2, rbuff3
        character(len=clen) :: fname
        character(len=clen) :: tempname
        character(len=clen) :: OutputDir
        integer :: tid

        rbuff1 => ig%rbuffxC(:,:,:,1); rbuff2 => ig%rbuffyC(:,:,:,1);
        rbuff3 => ig%rbuffzC(:,:,:,1); 
        gpC => ig%gpC
        tid = ig%step

        tidSUM = tidSUM + 1

        ! Compute u - mean 
        call transpose_x_to_y(ig%u,rbuff2,ig%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,ig%gpC)
        call compute_z_mean(rbuff3, u_mean)

        ! Compute v - mean 
        call transpose_x_to_y(ig%v,rbuff2,ig%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,ig%gpC)
        call compute_z_mean(rbuff3, v_mean)

        ! Compute w - mean 
        call transpose_x_to_y(ig%wC,rbuff2,ig%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,ig%gpC)
        call compute_z_mean(rbuff3, w_mean)

        ! Compute uu - mean 
        rbuff1 = ig%u*ig%u
        call transpose_x_to_y(rbuff1,rbuff2,ig%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,ig%gpC)
        call compute_z_mean(rbuff3, uu_mean)

        ! Compute vv - mean 
        rbuff1 = ig%v*ig%v
        call transpose_x_to_y(rbuff1,rbuff2,ig%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,ig%gpC)
        call compute_z_mean(rbuff3, vv_mean)

        ! Compute ww - mean 
        rbuff1 = ig%w*ig%w
        call transpose_x_to_y(rbuff1,rbuff2,ig%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,ig%gpC)
        call compute_z_mean(rbuff3, ww_mean)

        ! Compute uw - mean 
        rbuff1 = ig%u*ig%w
        call transpose_x_to_y(rbuff1,rbuff2,ig%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,ig%gpC)
        call compute_z_mean(rbuff3, uw_mean)

        call ig%compute_vorticity()
        
        ! Compute ox-ox mean
        rbuff1 = ig%ox*ig%ox
        call transpose_x_to_y(rbuff1,rbuff2,ig%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,ig%gpC)
        call compute_z_mean(rbuff3, oxox_mean)

        ! Compute oy-oy mean
        rbuff1 = ig%oy*ig%oy
        call transpose_x_to_y(rbuff1,rbuff2,ig%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,ig%gpC)
        call compute_z_mean(rbuff3, oyoy_mean)

        ! Compute oz-oz mean
        rbuff1 = ig%ox*ig%ox
        call transpose_x_to_y(rbuff1,rbuff2,ig%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,ig%gpC)
        call compute_z_mean(rbuff3, ozoz_mean)

        runningSum = runningSum + zStats2dump
        TemporalMnNOW = runningSum/tidSUM

        if (nrank == 0) then
            write(tempname,"(A3,I2.2,A2,I6.6,A4)") "Run", ig%RunID,"_t",tid,".stt"
            fname = ig%OutputDir(:len_trim(ig%OutputDir))//"/"//trim(tempname)
            call write_2d_ascii(TemporalMnNOW,fname)
        end if 
        call message(0, "Just dumped a .stt file")

    end subroutine 

    subroutine compute_z_mean(arr_in, vec_out)
        use reductions, only: P_SUM
        real(rkind), dimension(:,:,:), intent(in) :: arr_in
        real(rkind), dimension(:), intent(out) :: vec_out
        integer :: k

        do k = 1,nzg
            vec_out(k) = P_SUM(sum(arr_in(:,:,k)))/(real(nxg,rkind)*real(nyg,rkind))
        end do 

    end subroutine

    subroutine finalize_stats
        nullify( u_mean, v_mean, w_mean, uu_mean, vv_mean, ww_mean, uw_mean) 
        nullify( oxox_mean, oyoy_mean, ozoz_mean)
        deallocate(zStats2dump, runningSum, TemporalMnNOW)
    end subroutine 

end module 

