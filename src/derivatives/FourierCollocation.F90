! Routines specific to Fourier Collocation scheme
! This file is included in the Derivatives module

subroutine ddx_four(this, f, df) 
    
    class( derivatives ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%nx) :: wrkArr
    integer :: i, j, k

    select case (this%xoprank)
    case (1)
        if (.not. fourUseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % xfft % fourder1( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % xfft % fourder1( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. fourUseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % xfft % fourder1( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % xfft % fourder1( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. fourUseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % xfft % fourder1( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % xfft % fourder1( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 

subroutine ddy_four(this, f, df) 
    
    class( derivatives ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%ny) :: wrkArr
    integer :: i, j, k

    select case (this%yoprank)
    case (1)
        if (.not. fourUseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % yfft % fourder1( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % yfft % fourder1( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. fourUseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % yfft % fourder1( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % yfft % fourder1( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. fourUseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % yfft % fourder1( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % yfft % fourder1( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 


subroutine ddz_four(this, f, df) 
    
    class( derivatives ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%nz) :: wrkArr
    integer :: i, j, k

    select case (this%zoprank)
    case (1)
        if (.not. fourUseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % zfft % fourder1( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % zfft % fourder1( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. fourUseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % zfft % fourder1( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % zfft % fourder1( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. fourUseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % zfft % fourder1( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % zfft % fourder1( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 


subroutine d2dx2_four(this, f, df) 
    
    class( derivatives ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%nx) :: wrkArr
    integer :: i, j, k

    select case (this%xoprank)
    case (1)
        if (.not. fourUseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % xfft % fourder2( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % xfft % fourder2( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. fourUseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % xfft % fourder2( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % xfft % fourder2( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. fourUseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % xfft % fourder2( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % xfft % fourder2( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 

subroutine d2dy2_four(this, f, df) 
    
    class( derivatives ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%ny) :: wrkArr
    integer :: i, j, k

    select case (this%yoprank)
    case (1)
        if (.not. fourUseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % yfft % fourder2( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % yfft % fourder2( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. fourUseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % yfft % fourder2( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % yfft % fourder2( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. fourUseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % yfft % fourder2( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % yfft % fourder2( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 


subroutine d2dz2_four(this, f, df) 
    
    class( derivatives ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%nz) :: wrkArr
    integer :: i, j, k

    select case (this%zoprank)
    case (1)
        if (.not. fourUseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % zfft % fourder2( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % zfft % fourder2( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. fourUseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % zfft % fourder2( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % zfft % fourder2( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. fourUseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % zfft % fourder2( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % zfft % fourder2( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 


