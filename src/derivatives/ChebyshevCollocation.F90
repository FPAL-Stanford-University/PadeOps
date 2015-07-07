! Routines specific to Chebyshev Collocation scheme
! This file is included in the Derivatives module

subroutine ddx_cheb(this, f, df) 
    
    class( derivative ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%nx) :: wrkArr
    integer :: i, j, k

    select case (this%xoprank)
    case (1)
        if (.not. chebUseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % xdct % chebder1( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % xdct % chebder1( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. chebUseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % xdct % chebder1( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % xdct % chebder1( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. chebUseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % xdct % chebder1( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % xdct % chebder1( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 

subroutine ddy_cheb(this, f, df) 
    
    class( derivative ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%ny) :: wrkArr
    integer :: i, j, k

    select case (this%yoprank)
    case (1)
        if (.not. chebUseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % ydct % chebder1( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % ydct % chebder1( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. chebUseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % ydct % chebder1( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % ydct % chebder1( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. chebUseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % ydct % chebder1( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % ydct % chebder1( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 


subroutine ddz_cheb(this, f, df) 
    
    class( derivative ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%nz) :: wrkArr
    integer :: i, j, k

    select case (this%zoprank)
    case (1)
        if (.not. chebUseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % zdct % chebder1( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % zdct % chebder1( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. chebUseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % zdct % chebder1( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % zdct % chebder1( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. chebUseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % zdct % chebder1( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % zdct % chebder1( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 


subroutine d2dx2_cheb(this, f, df) 
    
    class( derivative ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%nx) :: wrkArr
    integer :: i, j, k

    select case (this%xoprank)
    case (1)
        if (.not. chebUseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % xdct % chebder2( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % xdct % chebder2( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. chebUseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % xdct % chebder2( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % xdct % chebder2( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. chebUseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % xdct % chebder2( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % xdct % chebder2( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 

subroutine d2dy2_cheb(this, f, df) 
    
    class( derivative ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%ny) :: wrkArr
    integer :: i, j, k

    select case (this%yoprank)
    case (1)
        if (.not. chebUseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % ydct % chebder2( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % ydct % chebder2( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. chebUseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % ydct % chebder2( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % ydct % chebder2( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. chebUseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % ydct % chebder2( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % ydct % chebder2( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 


subroutine d2dz2_cheb(this, f, df) 
    
    class( derivative ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%nz) :: wrkArr
    integer :: i, j, k

    select case (this%zoprank)
    case (1)
        if (.not. chebUseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % zdct % chebder2( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % zdct % chebder2( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. chebUseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % zdct % chebder2( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % zdct % chebder2( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. chebUseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % zdct % chebder2( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % zdct % chebder2( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 


