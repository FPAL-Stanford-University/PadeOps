program upsample1d
    use kind_parameters, only: rkind
    use gridtools, only: upsample_Periodic_1d 
    use constants, only: two, pi

    integer, parameter :: nx = 16
    real(rkind), dimension(nx) :: yS, xS
    real(rkind), dimension(2*nx) :: yL
    integer :: i
    real(rkind) :: dx

    dx = two*pi/real(nx,rkind)
    do i = 2,nx
        xS(i) = xS(i-1) + dx
    end do 

    yS = cos(xS)

    call upsample_Periodic_1d(yS,yL, nx)

    print*, yL



end program 
