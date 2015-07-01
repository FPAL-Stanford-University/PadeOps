! Routines for Fourier derivatives

function InitFFT() result(ierr)
    integer :: ierr
    ierr = xfft % init(nx_,dx_)
    ierr = yfft % init(ny_,dy_)
    ierr = zfft % init(nz_,dz_)
    k1 = xfft % k
    k2 = yfft % k
    k3 = zfft % k
end function

function ddxFOUR(f), result(df)
#ifdef TEST_FOUR
    use kind_parameters, only: rkind
    use constants, only: imi
#endif
    real(rkind), dimension(nx), intent(in) :: f
    real(rkind), dimension(nx) :: df

    associate ( k => xfft%
    df = xfft%ifft(imi*xfft




function InitFOUR(direction) result (ierr)
    character(len=1), intent(in) :: direction
    integer :: ierr

    ierr = 0
end function
