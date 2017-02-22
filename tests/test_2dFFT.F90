program test2dFFT
    use mpi
    use kind_parameters, only : rkind
    use decomp_2d
    use constants, only: pi, two, one, imi
    use spectralMod, only: spectral

    implicit none
   
    real(rkind), dimension(:,:,:), allocatable :: f, fr, x, y, z, dfdx, lapf
    complex(rkind), dimension(:,:,:), allocatable :: fhat
    type(decomp_info) :: gp
    
    type(spectral), allocatable :: spect

    integer :: nx = 32, ny = 32, nz = 32
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k, ii, jj, kk
    real(rkind) :: dx, dy, dz
    integer :: dimTransform = 2

    call MPI_Init(ierr)
    
    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)
    

    dx = two*pi/nx; dy = two*pi/ny; dz = two*pi/nz
    allocate(spect)
    !call spect%init("x", nx, ny, nz, dx, dy, dz, "four", "2/3rd", dimTransform)
    call spect%init("x", nx, ny, nz, dx, dy, dz, &
                "four", "2/3rd", 2 , .false.)

    allocate( f ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( x ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( y ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( z ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( fr( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate(dfdx( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate(lapf( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    
    call spect%alloc_r2c_out(fhat)
  

    kk = 1
    do k = gp%xst(3),gp%xen(3)
        jj = 1
        do j = gp%xst(2),gp%xen(2)
            ii = 1
            do i = gp%xst(1),gp%xen(1)
                x(ii,jj,kk) = (i-1)*dx
                y(ii,jj,kk) = (j-1)*dy
                z(ii,jj,kk) = (k-1)*dz
                ii = ii + 1
            end do 
            jj = jj + 1
        end do 
        kk = kk + 1
    end do 

    !f = cos(x)*sin(y)*(z*(two*pi - z))
    f = sin(x)*cos(y)
    dfdx = cos(x)*cos(y)
    lapf = -2.d0*cos(y)*sin(x)
    call sleep(nrank)

    call spect%fft(f,fhat)
    fhat = imi*spect%k1*fhat 
    call spect%ifft(fhat,fr)
    print*, "Done"
    print*, maxval(abs(fr - dfdx))
    
    call spect%fft(f,fhat)
    fhat = -spect%kabs_sq*fhat 
    call spect%ifft(fhat,fr)
    print*, maxval(abs(fr - lapf))
   

    call spect%destroy()
    deallocate(spect)


    call MPI_Finalize(ierr)

end program 
