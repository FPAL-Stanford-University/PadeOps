program test_cd06

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two,pi
    use cd06stuff,       only: cd06
    implicit none

    integer, parameter :: nx = 128
    logical, parameter :: periodic = .TRUE.

    type( cd06 ) :: mycd06
    real(rkind), dimension(nx) :: x,f,df,df_exact
    real(rkind) :: dx

    integer :: i,ierr

    dx = two*pi/real(nx,rkind)

    do i=1,nx
        x(i) = real(i-1,rkind)*two*pi/real(nx,rkind)
        f(i) = sin(4._rkind * x(i))
        df_exact(i) = 4._rkind * cos( 4._rkind * x(i))
    end do

    ierr = mycd06%init( nx, dx, periodic, 0, 0)

    df = mycd06 % cd06der1(f)

    print*, "Maximum error = ", MAXVAL( ABS(df - df_exact))

end program
