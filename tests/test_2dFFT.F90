program test2dFFT
    use mpi
    use kind_parameters, only : rkind
    use decomp_2d
    use constants, only: pi, two, one
    use spectralMod, only: spectral

    implicit none
   
    real(rkind), dimension(:,:,:), allocatable :: f, fr
    complex(rkind), dimension(:,:,:), allocatable :: fhat
    type(decomp_info) :: gp
    
    type(spectral), allocatable :: spect

    integer :: nx = 128, ny = 128, nz = 128
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k
    real(rkind) :: dx = one, dy = one, dz = one
    integer :: dimTransform = 2

    call MPI_Init(ierr)
    
    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)

    allocate(spect)
    call spect%init("x", nx, ny, nz, dx, dy, dz, "four", "2/3rd", dimTransform)

    allocate( f ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( fr( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    call spect%alloc_r2c_out(fhat)
   
    do k = gp%xst(3),gp%xen(3)
        do j = gp%xst(2),gp%xen(2)
            do i = gp%xst(1),gp%xen(1)
                f(i,j,k) = real(i + j + k,rkind)
            end do 
        end do 
    end do 

    call spect%fft(f,fhat)
    !print*, fhat(1:5,1,1)
    call spect%ifft(fhat,fr)
    print*, maxval(abs(fr - f))
    call spect%destroy()
    deallocate(spect)


    call MPI_Finalize(ierr)

end program 
