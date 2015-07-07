program test_dcts

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two,pi
    use dctstuff,        only: dcts
    implicit none

    integer, parameter :: nx = 2049
    logical, parameter :: periodic = .TRUE.

    type( dcts ) :: mydcts
    real(rkind), dimension(nx) :: x,f,df,df_exact

    integer :: i,ierr

    do i=1,nx
        x(i) = - cos( real((i-1),rkind) * pi / real((nx-1),rkind) )
        f(i) = cos( x(i) )
        df_exact(i) = -sin( x(i) )
    end do

    ierr = mydcts%init( nx )

    df = mydcts % chebder1(f)

    print*, "Maximum error = ", MAXVAL( ABS(df - df_exact))

end program
