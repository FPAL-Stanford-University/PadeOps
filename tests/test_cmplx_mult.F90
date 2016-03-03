module complx_mults
    use kind_parameters, only: rkind
    use constants, only: imi

    interface mult_by_ik
        module procedure mult_by_ik_oop, mult_by_ik_ip
    end interface 
contains
    pure subroutine mult_by_ik_oop(kmat,fhat,dfhat)
        real(rkind), dimension(:,:,:), intent(in) :: kmat
        complex(rkind), dimension(:,:,:), intent(in) :: fhat
        complex(rkind), dimension(:,:,:), intent(out) :: dfhat

        integer :: i, j, k
        real(rkind) :: rpart, ipart

        do k = 1,size(fhat,3)
            do j = 1,size(fhat,2)
                do i = 1,size(fhat,1)
                    rpart = -kmat(i,j,k)*dimag(fhat(i,j,k))
                    ipart = kmat(i,j,k)*dreal(fhat(i,j,k))
                    dfhat(i,j,k) = dcmplx(rpart,ipart)
                end do 
            end do
        end do  
    end subroutine

    subroutine mult_by_ik_ip(kmat,fhat)
        real(rkind), dimension(:,:,:), intent(in) :: kmat
        complex(rkind), dimension(:,:,:), intent(inout) :: fhat

        integer :: i, j, k
        real(rkind) :: rpart, ipart

        do k = 1,size(fhat,3)
            do j = 1,size(fhat,2)
                do i = 1,size(fhat,1)
                    rpart = -kmat(i,j,k)*dimag(fhat(i,j,k))
                    ipart = kmat(i,j,k)*dreal(fhat(i,j,k))
                    fhat(i,j,k) = dcmplx(rpart,ipart)
                end do 
            end do
        end do  
    end subroutine

    subroutine mult_by_i(fhat)
        complex(rkind), dimension(:,:,:), intent(inout) :: fhat

        integer :: i, j, k
        real(rkind) :: rpart, ipart

        do k = 1,size(fhat,3)
            do j = 1,size(fhat,2)
                do i = 1,size(fhat,1)
                    rpart = -dimag(fhat(i,j,k))
                    ipart = dreal(fhat(i,j,k))
                    fhat(i,j,k) = dcmplx(rpart,ipart)
                end do 
            end do
        end do  
    end subroutine
end module 



program test_cmplx_mult
    use kind_parameters, only: rkind
    use timer, only: tic, toc
    use spectralMod, only: spectral
    use decomp_2d
    use constants, only: imi, two, pi, one, zero
    use exits, only: message
    use complx_mults, only: mult_by_ik, mult_by_i

    implicit none
   
    real(rkind), dimension(:,:,:), allocatable :: dfdx, fr, f, x, y, z
    complex(rkind), dimension(:,:,:), allocatable :: fhat, fhat2

    type(decomp_info) :: gp
    
    type(spectral), allocatable :: spect

    integer :: nx = 256, ny = 256, nz = 256
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k, ii, jj, kk
    real(rkind) :: dx, dy, dz
    integer :: dimTransform = 2

    call MPI_Init(ierr)
    
    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)
    

    dx = two*pi/nx; dy = two*pi/ny; dz = two*pi/nz
    allocate(spect)
    call spect%init("x", nx, ny, nz, dx, dy, dz, "four", "2/3rd", dimTransform)

    allocate( f ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( x ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( y ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( z ( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate( fr( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    allocate(dfdx( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) )
    
    call spect%alloc_r2c_out(fhat)
    call spect%alloc_r2c_out(fhat2)
  

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

    f = cos(x)*sin(y)*(z*(two*pi - z))
    dfdx = cos(x)*cos(y)*(z*(two*pi - z))


    call spect%fft(f,fhat)

    call tic()
    fhat = imi*spect%k2*fhat 
    call toc()
    
    !print*, fhat(1:5,1,1)
    call spect%ifft(fhat,fr)


    call spect%fft(f,fhat)
    call tic()
    call mult_by_ik(spect%k2, fhat) 
    call toc()
    call spect%ifft(fhat,fr)
    
    call message("Maximum Error:", maxval(abs(fr - dfdx)))

    call tic()
    fhat = imi*fhat
    call toc()

    call tic()
    call mult_by_i(fhat)
    call toc()

    call spect%destroy()
    deallocate(spect)



    call MPI_Finalize(ierr)

end program 
