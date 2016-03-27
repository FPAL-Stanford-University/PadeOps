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
    real(rkind) :: randomScaleFact = 0.05_rkind ! 5% of the mean value

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

    Uperiods = 6.0; Vperiods = 6.0; zpeak = 0.3;
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

module allStatistics
    use kind_parameters, only: rkind
    use decomp_2d
    use IncompressibleGridNP, only: igrid 
    use constants, only: zero  
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

        allocate(zStats2dump(nzg,7))
        allocate(runningSum(nzg,7))
        allocate(TemporalMnNOW(nzg,7))
        u_mean => zStats2dump(:,1); v_mean => zStats2dump(:,2); w_mean => zStats2dump(:,3); 
        uu_mean => zStats2dump(:,4); vv_mean => zStats2dump(:,5); ww_mean => zStats2dump(:,6); 
        uw_mean => zStats2dump(:,7)

        runningSum = zero
        nullify(gpC)
    end subroutine

    subroutine compute_stats(ig)
        use mpi
        type(igrid), intent(inout), target :: ig
        type(decomp_info), pointer :: gpC
        real(rkind), dimension(:,:,:), pointer :: rbuff1, rbuff2, rbuff3, rbuff4

        rbuff1 => ig%rbuffxC(:,:,:,1); rbuff2 => ig%rbuffyC(:,:,:,1);
        rbuff3 => ig%rbuffzC(:,:,:,1); 
        rbuff4 => ig%rbuffzC(:,:,:,2); 
        gpC => ig%gpC

        tidSUM = tidSUM + 1

        ! Compute u - mean 
        call transpose_x_to_y(ig%u,rbuff2,ig%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,ig%gpC)
        call compute_z_mean(rbuff3, u_mean)
        ! uu mean
        call compute_z_fluct(rbuff3)
        rbuff4 = rbuff3*rbuff3
        call compute_z_mean(rbuff4, uu_mean)

        ! Compute w - mean 
        call transpose_x_to_y(ig%wC,rbuff2,ig%gpC)
        call transpose_y_to_z(rbuff2,rbuff4,ig%gpC)
        call compute_z_mean(rbuff4, w_mean)
        call compute_z_fluct(rbuff4)
        ! uw mean
        rbuff3 = rbuff3*rbuff4
        call compute_z_mean(rbuff3, uw_mean)
        ! ww mean 
        rbuff3 = rbuff4*rbuff4
        call compute_z_mean(rbuff3, ww_mean)
        
        ! Compute v - mean 
        call transpose_x_to_y(ig%v,rbuff2,ig%gpC)
        call transpose_y_to_z(rbuff2,rbuff3,ig%gpC)
        call compute_z_mean(rbuff3, v_mean)
        call compute_z_fluct(rbuff3)
        ! vv mean 
        rbuff4 = rbuff3*rbuff3
        call compute_z_mean(rbuff4, vv_mean)

        runningSum = runningSum + zStats2dump

    end subroutine 

    subroutine dump_stats(ig)
        use basic_io, only: write_2d_ascii, write_2D_binary
        use exits, only: message
        use kind_parameters, only: clen
        use mpi
        type(igrid), intent(in), target :: ig
        character(len=clen) :: fname
        character(len=clen) :: tempname
        character(len=clen) :: OutputDir
        integer :: tid

        TemporalMnNOW = runningSum/real(tidSUM,rkind)
        tid = ig%step

        if (nrank == 0) then
            write(tempname,"(A3,I2.2,A2,I6.6,A4)") "Run", ig%RunID,"_t",tid,".stt"
            fname = ig%OutputDir(:len_trim(ig%OutputDir))//"/"//trim(tempname)
            call write_2d_ascii(TemporalMnNOW,fname)
            !call write_2D_binary(TemporalMnNOW,fname)
        end if
        call message(1, "Just dumped a .stt file")
        call message(2, "Number ot tsteps averaged:",tidSUM)

    end subroutine

    subroutine compute_z_fluct(fin)
        use reductions, only: P_SUM
        real(rkind), dimension(:,:,:), intent(inout) :: fin
        integer :: k
        real(rkind) :: fmean

        do k = 1,size(fin,3)
            fmean = P_SUM(sum(fin(:,:,k)))/(real(nxg,rkind)*real(nyg,rkind))
            fin(:,:,k) = fin(:,:,k) - fmean
        end do 

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

