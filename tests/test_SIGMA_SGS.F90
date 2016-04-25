program test_SIGMA_SGS
    use kind_parameters, only: rkind, clen
    use constants, only: pi, half, two, one, imi, zero
    use reductions, only: p_maxval
    use timer, only: tic, toc
    use basic_io, only: read_2d_ascii 
    use decomp_2d
    use spectralMod, only: spectral
    use cd06staggstuff, only: cd06stagg 
    use PadePoissonMod, only: padepoisson
    use exits, only: message
    use cf90stuff, only: cf90
    use basic_io, only: write_binary
    use sgsmod, only: sgs
    use numerics, only: useCompactFD

    implicit none
   
    real(rkind), dimension(:,:,:), allocatable :: x, y, z, u, v, w, wE
    real(rkind), dimension(:,:,:), allocatable :: tmp
    type(decomp_info) :: gpC, gpE
    integer :: nx = 100, ny = 100, nz = 60
    integer :: ierr, prow = 0, pcol = 0
    real(rkind) :: dx, dy, dz
    real(rkind), dimension(:,:), allocatable :: temp 
    type(spectral), allocatable :: spect
    type(spectral), allocatable :: spectE
    type(padepoisson), allocatable :: poiss 
    real(rkind), dimension(:,:,:), allocatable :: divergence, rtmp1z, rtmp2z, rtmp1x, rtmp1y, rtmp1yE, rtmp1zE
    complex(rkind), dimension(:,:,:), allocatable :: urhs, vrhs, wrhs
    complex(rkind), dimension(:,:,:), allocatable :: uhat, vhat, what
    type(cd06stagg), allocatable :: derZu, derZv, derZw
    complex(rkind), dimension(:,:,:), allocatable :: ctmp1, ctmp2, ctmpz1,ctmpz2
    type(cf90) :: filzC, filzE
    real(rkind) :: Lx, Ly, Lz
    real(rkind) :: Uperiods, Vperiods, zpeak, epsfac, delta_Ek
    integer :: ix1, ixn, iy1, iyn, iz1, izn
    integer :: i,j,k
    real(rkind), dimension(:,:,:,:), allocatable :: duidxj
    complex(rkind), dimension(:,:,:,:), allocatable :: duidxjH
    integer :: dimTransform = 2
    character(len=clen) :: filename = "/home/aditya90/Codes/PadeOps/data/OpenFoam_AllData.txt" 
    type(sgs) :: sgsModel
    real(rkind) :: maxnuSGS
   
    integer :: ModelID = 0
    logical :: useEkmanInit = .false. 
    logical :: useDynamicProcedure = .true. 
    logical :: useClipping = .true. 
    logical :: useCompact = .true. 


    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, prow, pcol)

    call get_decomp_info(gpC)
    call decomp_info_init(nx, ny, nz+1, gpE)
   
    useCompactFD = useCompact

    allocate( x ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( y ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( z ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( u ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( v ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( w ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( wE ( gpE%xsz(1), gpE%xsz(2), gpE%xsz(3) ) )
    allocate( duidxj( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) , 9) )
    allocate( divergence ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
   
    allocate(spect, spectE)

    nx = gpC%xsz(1); ny = gpC%ysz(2); nz = gpC%zsz(3)

    if (useEkmanInit) then
        ! If base gpCosition is in Y
        ix1 = gpC%xst(1); iy1 = gpC%xst(2); iz1 = gpC%xst(3)
        ixn = gpC%xen(1); iyn = gpC%xen(2); izn = gpC%xen(3)
        
        Lx = 1.5_rkind*pi; Ly = 1.5_rkind*pi; Lz = 1._rkind

        dx = Lx/real(nx,rkind); dy = Ly/real(ny,rkind); dz = Lz/real(nz,rkind)

        do k=1,size(x,3)
            do j=1,size(x,2)
                do i=1,size(x,1)
                    x(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                    y(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                    z(i,j,k) = real( iz1 + k - 1, rkind ) * dz + dz/two
                end do
            end do
        end do
   
        x = x - dx
        y = y - dy
        z = z - dz 

        Uperiods = 2.0; Vperiods = 2.0; zpeak = 0.3;
        epsfac = 0.5d0; delta_Ek = 0.14071 

        u = one - exp(-z/delta_Ek)*cos(z/delta_Ek) &
                + half*exp(half)*(z/Lz)*cos(Uperiods*two*pi*y/Ly)*exp(-half*(z/zpeak/Lz)**2)
        v = exp(-z/delta_Ek)*sin(z/delta_Ek) + &
                + half*exp(half)*(z/Lz)*cos(Vperiods*two*pi*x/Lx)*exp(-half*(z/zpeak/Lz)**2)
        w = zero 

    else 
        call read_2d_ascii(temp,filename)
        allocate(tmp(nx,ny,nz))
        tmp = reshape(temp(:,4),[nx,ny,nz])
        u = tmp(gpC%xst(1):gpC%xen(1),gpC%xst(2):gpC%xen(2),gpC%xst(3):gpC%xen(3))

        tmp = reshape(temp(:,5),[nx,ny,nz])
        v = tmp(gpC%xst(1):gpC%xen(1),gpC%xst(2):gpC%xen(2),gpC%xst(3):gpC%xen(3))
        
        tmp = reshape(temp(:,6),[nx,ny,nz])
        w = tmp(gpC%xst(1):gpC%xen(1),gpC%xst(2):gpC%xen(2),gpC%xst(3):gpC%xen(3))
        
        tmp = reshape(temp(:,1),[nx,ny,nz])
        x = tmp(gpC%xst(1):gpC%xen(1),gpC%xst(2):gpC%xen(2),gpC%xst(3):gpC%xen(3))
        
        tmp = reshape(temp(:,2),[nx,ny,nz])
        y = tmp(gpC%xst(1):gpC%xen(1),gpC%xst(2):gpC%xen(2),gpC%xst(3):gpC%xen(3))
        
        tmp = reshape(temp(:,3),[nx,ny,nz])
        z = tmp(gpC%xst(1):gpC%xen(1),gpC%xst(2):gpC%xen(2),gpC%xst(3):gpC%xen(3))
        
        dx = x(2,1,1) - x(1,1,1)
        dy = y(1,2,1) - y(1,1,1)
        dz = z(1,1,2) - z(1,1,1)
   
        !deallocate(tmp, temp) 
    end if 

    call spect%init("x", nx, ny, nz, dx, dy, dz, "four", "2/3rd", dimTransform,.false.)
    call spectE%init("x", nx, ny, nz+1, dx, dy, dz, "four", "2/3rd", dimTransform,.false.)
    call spect%alloc_r2c_out(ctmp1)
    call spect%alloc_r2c_out(ctmp2)
    call spect%alloc_r2c_out(urhs)
    call spect%alloc_r2c_out(vrhs)
    call spectE%alloc_r2c_out(wrhs)
    call spect%alloc_r2c_out(uhat)
    call spect%alloc_r2c_out(vhat)
    call spect%alloc_r2c_out(what)
    call spect%alloc_r2c_out(duidxjH,9)
    
    urhs = zero
    vrhs = zero
    wrhs = zero

    uhat = zero
    vhat = zero
    what = zero
    
    allocate(derZu, derZv, derZw)
    call derZu%init(nz,dz,.true.,.false.,.false.,.false.)
    call derZv%init(nz,dz,.true.,.false.,.false.,.false.)
    call derZw%init(nz,dz,.false.,.false.,.true.,.true.)

    allocate(ctmpz1(spect%spectdecomp%zsz(1),spect%spectdecomp%zsz(2),spect%spectdecomp%zsz(3)))
    allocate(ctmpz2(spect%spectdecomp%zsz(1),spect%spectdecomp%zsz(2),spect%spectdecomp%zsz(3)))


    call spect%fft(u,ctmp1)
    call spect%mtimes_ik1_oop(ctmp1,ctmp2)
    call spect%ifft(ctmp2,duidxj(:,:,:,1))
    call spect%mtimes_ik2_oop(ctmp1,ctmp2)
    call spect%ifft(ctmp2,duidxj(:,:,:,2))
    call transpose_y_to_z(ctmp1,ctmpz1,spect%spectdecomp)
    call derZu%ddz_C2C(ctmpz1,ctmpz2, size(ctmpz1,1), size(ctmpz1,2))
    call transpose_z_to_y(ctmpz2,ctmp2,spect%spectdecomp)
    call spect%ifft(ctmp2,duidxj(:,:,:,3))

    call spect%fft(v,ctmp1)
    call spect%mtimes_ik1_oop(ctmp1,ctmp2)
    call spect%ifft(ctmp2,duidxj(:,:,:,4))
    call spect%mtimes_ik2_oop(ctmp1,ctmp2)
    call spect%ifft(ctmp2,duidxj(:,:,:,5))
    call transpose_y_to_z(ctmp1,ctmpz1,spect%spectdecomp)
    call derZv%ddz_C2C(ctmpz1,ctmpz2, size(ctmpz1,1), size(ctmpz1,2))
    call transpose_z_to_y(ctmpz2,ctmp2,spect%spectdecomp)
    call spect%ifft(ctmp2,duidxj(:,:,:,6))

    call spect%fft(w,ctmp1)
    call spect%mtimes_ik1_oop(ctmp1,ctmp2)
    call spect%ifft(ctmp2,duidxj(:,:,:,7))
    call spect%mtimes_ik2_oop(ctmp1,ctmp2)
    call spect%ifft(ctmp2,duidxj(:,:,:,8))
    call transpose_y_to_z(ctmp1,ctmpz1,spect%spectdecomp)
    call derZw%ddz_C2C(ctmpz1,ctmpz2, size(ctmpz1,1), size(ctmpz1,2))
    call transpose_z_to_y(ctmpz2,ctmp2,spect%spectdecomp)
    call spect%ifft(ctmp2,duidxj(:,:,:,9))

    do i = 1,9
        call spect%fft(duidxj(:,:,:,i),duidxjH(:,:,:,i))
    end do 

    call spect%fft(u,uhat)
    call spect%fft(v,vhat)
    call spect%fft(w,what)

    call sgsModel%init(modelID, spect, spectE, gpC, gpE, dx, dy, dz, useDynamicProcedure, useClipping, z, 1, 0.1d0, .false., .true. )
    !duidxj(2,3,4,1) = 1.d0; duidxj(2,3,4,2) = 2.d0; duidxj(2,3,4,3) = 3.d0
    !duidxj(2,3,4,4) = 4.d0; duidxj(2,3,4,5) = -1.d0; duidxj(2,3,4,6) = 6.d0

    do i = 1,9
        call spect%fft(duidxj(:,:,:,i),duidxjH(:,:,:,i))
    end do 

    call tic()
    call sgsModel%getRHS_SGS(duidxj, duidxjH, urhs, vrhs, wrhs, uhat, vhat, what, u, v, w, maxnuSGS)
    call toc()

    call spectE%ifft(wrhs,wE)
    
    !print*, wE(5,4,:)
    call message("Max val RHS:", p_maxval(u))
    call message("Max val nuSGS:", maxnuSGS)

    call derZu%destroy()
    call derZv%destroy()
    call derZw%destroy()
    call spect%destroy()
    deallocate(derZu, derZv, derZw)
    deallocate(spect)
end program 
