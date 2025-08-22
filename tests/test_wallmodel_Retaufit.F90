program test_wallmodel_Retau
    implicit none
    integer, parameter :: rkind = selected_real_kind(15, 307)
    real(rkind) :: Re_del, Re_tauf
    integer :: i, npts
    real(rkind) :: Re_min, Re_max, logRe, dlogRe
    character(len=100) :: filename
    open(unit=10, file='Re_del_vs_Re_tauf.dat', status='replace', action='write')

    ! Logarithmic spacing between 1e-2 and 1e7
    npts = 10000
    Re_min = 1.0e-2_rkind
    Re_max = 1.0e+7_rkind
    dlogRe = (log10(Re_max) - log10(Re_min)) / real(npts - 1, rkind)

    do i = 1, npts
        logRe = log10(Re_min) + real(i - 1, rkind) * dlogRe
        Re_del = 10.0_rkind ** logRe
        Re_tauf = get_Retaufit(Re_del)
        write(10, '(E15.7E3,1X,E15.7E3)') Re_del, Re_tauf
    end do

    close(10)
    print *, 'Output written to Re_del_vs_Re_tauf.dat'

contains

    function get_Retaufit(Re_del) result(Re_tauf)
        real(rkind), intent(in) :: Re_del
        real(rkind) :: b1, b2, k3, k4, p1, p2, p3
        real(rkind) :: Re_tauf

        b1 = (1.0d0 + 1.55d-1 * Re_del**(-3.0d-2))**(-1.0d0)
        b2 = 1.7d0 - (1.0d0 + 3.6d1 * Re_del**(-7.5d-1))**(-1.0d0)

        k3 = 5.0d-3
        k4 = k3 ** (b1 - 5.0d-1)
        p1 = k4 * Re_del ** b1
        p2 = 1.0d0 + (k3 * Re_del) ** (-b2)
        p3 = (b1 - 5.0d-1) / b2

        Re_tauf = p1 * p2 ** p3
    end function

end program test_wallmodel_Retau

