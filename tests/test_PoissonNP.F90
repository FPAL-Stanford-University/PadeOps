program testPoissonNP
    use mpi
    use kind_parameters, only : rkind
    use decomp_2d
    use constants, only: pi, two, one, imi, zero
    use spectralMod, only: spectral
    use poissonMod, only: poisson
    use reductions, only: p_maxval
    use timer, only: tic, toc

    implicit none
   
    real(rkind), dimension(:,:,:), allocatable :: f, fr, x, y, z, dfdx, rhs
    complex(rkind), dimension(:,:,:), allocatable :: fhat, phat, fhat_inZ, phat_inZ
    type(decomp_info) :: gp
    type(poisson) :: poiss 
    type(spectral), allocatable :: spect

    integer :: nx = 128, ny = 128, nz = 128
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k, ii, jj, kk
    real(rkind) :: dx, dy, dz, zbot, fmean
    integer :: dimTransform = 2

    call MPI_Init(ierr)
    
    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)
    

    dx = two*pi/nx; dy = two*pi/ny; dz = one/nz
    allocate(spect)
    call spect%init("x", nx, ny, nz, dx, dy, dz, "four", "2/3rd", dimTransform)
    call poiss%init(spect,.false.,dx, dy, dz)


    allocate( f ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( x ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( y ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( z ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( fr( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate(dfdx( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate(rhs( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    
    call spect%alloc_r2c_out(fhat)
    call spect%alloc_r2c_out(phat)
  
    call poiss%allocArrZ(fhat_inZ)
    call poiss%allocArrZ(phat_inZ)

    kk = 1
    do k = gp%xst(3),gp%xen(3)
        jj = 1
        do j = gp%xst(2),gp%xen(2)
            ii = 1
            do i = gp%xst(1),gp%xen(1)
                x(ii,jj,kk) = (i-1)*dx
                y(ii,jj,kk) = (j-1)*dy
                z(ii,jj,kk) = (k-1)*dz + dz/two
                ii = ii + 1
            end do 
            jj = jj + 1
        end do 
        kk = kk + 1
    end do 

    f = (- (z**2)*(2._rkind*z - 3._rkind)/6._rkind - one/12._rkind)!*cos(x)*sin(y)
    zbot = dz/two
    fmean = (- (zbot**2)*(2._rkind*zbot - 3._rkind)/6._rkind - one/12._rkind)
    print*, fmean
    !rhs = -cos(x)*sin(y)*(-(two*z**3)/3._rkind + z**2 + 2*z - 7._rkind/6._rkind)
    rhs = 1 - two*z

    call tic()
    call spect%fft(rhs,fhat)
    call transpose_y_to_z(fhat,fhat_inZ,poiss%sp_gp)
    call poiss%PoissonSolveZ(fhat_inZ,phat_inZ)
    call transpose_z_to_y(phat_inZ,phat,poiss%sp_gp)
    call spect%ifft(phat,fr)
    call toc()

    fr = fr + fmean
    print*, p_maxval(abs(fr - f))
    call spect%destroy()
    deallocate(spect)
    call MPI_Finalize(ierr)

end program 
