! Routines for Fourier derivatives

function InitFFT() result(ierr)
    integer :: ierr
    ierr = xfft % init(nx,dx)
    ierr = yfft % init(ny,dy)
    ierr = zfft % init(nz,dz)
    k1 = xfft % k
    k2 = yfft % k
    k3 = zfft % k

    ierr = 0
end function

function d1xFOUR(f), result(df)
#ifdef TEST_FOUR
    use kind_parameters, only: rkind
    use constants, only: imi
#endif
    real(rkind), dimension(nx), intent(in) :: f
    real(rkind), dimension(nx) :: df

    df = xfft%ifft(imi*k1*xfft%fft(f))

end function 

function d1yFOUR(f), result(df)
#ifdef TEST_FOUR
    use kind_parameters, only: rkind
    use constants, only: imi
#endif
    real(rkind), dimension(ny), intent(in) :: f
    real(rkind), dimension(ny) :: df

    df = yfft%ifft(imi*k2*yfft%fft(f))

end function 

function d1zFOUR(f), result(df)
#ifdef TEST_FOUR
    use kind_parameters, only: rkind
    use constants, only: imi
#endif
    real(rkind), dimension(nz), intent(in) :: f
    real(rkind), dimension(nz) :: df

    df = zfft%ifft(imi*k3*zfft%fft(f))

end function 

function d2xFOUR(f), result(df)
#ifdef TEST_FOUR
    use kind_parameters, only: rkind
    use constants, only: imi
#endif
    real(rkind), dimension(nx), intent(in) :: f
    real(rkind), dimension(nx) :: df

    df = xfft%ifft(-k1*k1*xfft%fft(f))

end function 

function d2yFOUR(f), result(df)
#ifdef TEST_FOUR
    use kind_parameters, only: rkind
    use constants, only: imi
#endif
    real(rkind), dimension(ny), intent(in) :: f
    real(rkind), dimension(ny) :: df

    df = yfft%ifft(-k2*k2*yfft%fft(f))

end function 

function d2zFOUR(f), result(df)
#ifdef TEST_FOUR
    use kind_parameters, only: rkind
    use constants, only: imi
#endif
    real(rkind), dimension(nz), intent(in) :: f
    real(rkind), dimension(nz) :: df

    df = zfft%ifft(-k3*k3*zfft%fft(f))

end function 

function InitFOUR(direction) result (ierr)
    character(len=1), intent(in) :: direction
    integer :: ierr

    ierr = 0
end function
