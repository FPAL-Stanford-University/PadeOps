program test_ffts

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two,pi
    use timer,           only: tic,toc
    use fftstuff,        only: ffts
    implicit none

    integer, parameter :: nx=512, ny=128, nz=128
    logical, parameter :: periodic = .TRUE.

    type( ffts ) :: myffts
    real(rkind), dimension(nx,ny,nz) :: x,f,df,df_exact
    real(rkind) :: dx

    integer :: i,j,k,ierr

    dx = two*pi/real(nx,rkind)

    do i=1,nx
        x(i,:,:) = real(i-1,rkind)*two*pi/real(nx,rkind)
        f(i,:,:) = sin(4._rkind * x(i,:,:))
        df_exact(i,:,:) = 4._rkind * cos( 4._rkind * x(i,:,:))
    end do

    ierr = myffts%init( nx, dx )

    call tic()
    do k=1,nz
        do j=1,ny
            df(:,j,k) = myffts % fourder1(f(:,j,k))
        end do
    end do
    call toc("Time to get x derivatives")

    print*, "Maximum error = ", MAXVAL( ABS(df - df_exact))

end program
