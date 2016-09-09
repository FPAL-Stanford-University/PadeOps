module wavenums
    use kind_parameters, only : rkind
    implicit none
contains
    pure function GetWaveNums(nx,dx) result(k)
        use constants, only: pi, two
        integer, intent(in) :: nx
        real(rkind), intent(in) :: dx
        real(rkind), dimension(nx) :: k

        integer :: i,dummy

        dummy = nx - MOD(nx,2)

        do i = 1,nx
            k(i) = ( -pi + (i-1)*two*pi/real(dummy,rkind) ) / dx
        end do

        k = ifftshift(k)

    end function
    
    pure function ifftshift(k) result(kshift)

        real(rkind), dimension(:), intent(in) :: k
        real(rkind), dimension(SIZE(k)) :: kshift
        integer :: n

        n = SIZE(k)

        select case ( MOD(n,2) )
        case (0)
            kshift(1:n/2) = k(n/2+1:n)
            kshift(n/2+1:n) = k(1:n/2)
        case (1)
            kshift(1:(n+1)/2) = k((n+1)/2:n)
            kshift((n+1)/2+1:n) = k(1:(n-1)/2)
        end select

    end function


end module


program test_fft3d_allreals
    use timer, only: tic, toc
    use kind_parameters, only : rkind
    use mpi
    use constants, only: pi
    use wavenums

    implicit none 

    include "fftw3.f"

    integer :: nx = 256, ny = 256, nz = 256
    real(rkind), dimension(:,:,:), allocatable :: x,y,z, f, fnew, fhat, fhatnew, g
    complex(rkind), dimension(:,:,:), allocatable :: tmp1, tmp2
    real(rkind), dimension(:), allocatable :: kx
    integer :: i, j, k
    real(rkind) :: dx, dy, dz, normfact
    integer(kind=8) :: plan_fftx, plan_ifftx, plan_ffty, plan_iffty, plan_fftz, plan_ifftz
    real(rkind) :: start, finish1, finish2, finish3

    allocate(f(nx,ny,nz), fhat(nx + 2, ny, nz), fhatnew(nx+2,ny,nz), g(nx, ny, nz), x(nx, ny,nz), y(nx, ny, nz), z(nx, ny, nz), fnew(nx, ny, nz))

    dx = 2.d0*pi/real(nx,rkind)
    dy = 2.d0*pi/real(ny,rkind)
    dz = 2.d0*pi/real(nz,rkind)
    normfact = 1.d0/real(nx,rkind)
    fhat = 0.d0
    do k = 1,nz
        do j = 1,ny
            do i = 1,nx
                x(i,j,k) = (i - 1)*dx
                y(i,j,k) = (j - 1)*dy
                z(i,j,k) = (k - 1)*dz
            end do 
        end do 
    end do 

    allocate(kx(nx))
    kx = getWavenums(nx,dx)

    f = cos(x)*cos(y)*cos(z)
    call dfftw_plan_many_dft_r2c(plan_fftx, 1, nx, &
         & ny*nz, f, nx, 1, nx, fhat, nx/2+1, 1, nx/2+1, FFTW_EXHAUSTIVE)

    call dfftw_plan_many_dft_c2r(plan_ifftx, 1, nx, &
         & ny*nz, fhat, nx/2+1, 1, nx/2+1, f, nx, 1, nx, FFTW_EXHAUSTIVE)

    start = MPI_WTIME()
    call dfftw_execute_dft_r2c(plan_fftx, f, fhat)  
    finish1 = MPI_WTIME()
    call mtimes_ik1_ip(fhat)
    !call mtimes_ik1_ip(fhat, fhatnew)
    finish2 = MPI_WTIME()
    call dfftw_execute_dft_c2r(plan_ifftx, fhat, fnew)
    !call dfftw_execute_dft_c2r(plan_ifftx, fhatnew, fnew)
    fnew = fnew*normfact 
    finish3 = MPI_WTIME()
    print*, "---------------"
    print*, "max error:", maxval(fnew + sin(x)*cos(y)*cos(z))
    print*, "Time for fft:", finish1 - start
    print*, "Time for mtimes_ik:", finish2 - finish1
    print*, "Time for ifft:", finish3 - finish2



    deallocate( f,fhat, g, x, y, z)


contains
    subroutine mtimes_ik1_ip(fhat)
        real(rkind), dimension(:, :, :), intent(inout) :: fhat 
        integer :: i, j, k
        real(rkind) :: tmp

        do k = 1,size(fhat,3)
            do j = 1,size(fhat,2)
                do i = 1, size(fhat,1)/2
                    tmp = fhat(2*i-1,j,k) 
                    fhat(2*i-1,j,k)  = - fhat(2*i,j,k)*kx(i)
                    fhat(2*i  ,j,k)  =   tmp          *kx(i)
                end do 
            end do 
        end do 


    end subroutine

    subroutine mtimes_ik1_oop(fhat, dfhatdx)
        real(rkind), dimension(:, :, :), intent(in) :: fhat 
        real(rkind), dimension(:, :, :), intent(out) :: dfhatdx 
        integer :: i, j, k

        do k = 1,size(fhat,3)
            do j = 1,size(fhat,2)
                do i = 1, size(fhat,1)/2
                    dfhatdx(2*i-1,j,k)  = - kx(i)*fhat(2*i  ,j,k)
                    dfhatdx(2*i  ,j,k)  =   kx(i)*fhat(2*i-1,j,k)
                end do 
            end do 
        end do 


    end subroutine
end program 
