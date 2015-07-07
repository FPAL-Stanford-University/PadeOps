program test_ffts

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two,pi
    use fftstuff,        only: ffts
    implicit none

    integer, parameter :: nx = 16
    logical, parameter :: periodic = .TRUE.

    type( ffts ) :: myffts
    real(rkind), dimension(nx) :: x,f,df,df_exact
    real(rkind) :: dx

    integer :: i,ierr

    dx = two*pi/real(nx,rkind)

    do i=1,nx
        x(i) = real(i-1,rkind)*two*pi/real(nx,rkind)
        f(i) = sin(4._rkind * x(i))
        df_exact(i) = 4._rkind * cos( 4._rkind * x(i))
    end do

    ierr = myffts%init( nx, dx )

    df = myffts % fourder1(f)

    print*, "Maximum error = ", MAXVAL( ABS(df - df_exact))

end program
