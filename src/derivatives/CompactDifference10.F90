! Routines specific to 10th Order Compact Difference scheme
! This file is included in the Derivatives module

subroutine ddx_cd10(this, f, df) 
    
    class( derivative ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%nx) :: wrkArr
    integer :: i, j, k

    select case (this%xoprank)
    case (1)
        if (.not. cd10UseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % xcd10 % cd10der1( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % xcd10 % cd10der1( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. cd10UseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % xcd10 % cd10der1( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % xcd10 % cd10der1( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. cd10UseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % xcd10 % cd10der1( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % xcd10 % cd10der1( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 

subroutine ddy_cd10(this, f, df) 
    
    class( derivative ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%ny) :: wrkArr
    integer :: i, j, k

    select case (this%yoprank)
    case (1)
        if (.not. cd10UseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % ycd10 % cd10der1( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % ycd10 % cd10der1( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. cd10UseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % ycd10 % cd10der1( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % ycd10 % cd10der1( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. cd10UseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % ycd10 % cd10der1( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % ycd10 % cd10der1( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 


subroutine ddz_cd10(this, f, df) 
    
    class( derivative ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%nz) :: wrkArr
    integer :: i, j, k

    select case (this%zoprank)
    case (1)
        if (.not. cd10UseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % zcd10 % cd10der1( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % zcd10 % cd10der1( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. cd10UseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % zcd10 % cd10der1( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % zcd10 % cd10der1( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. cd10UseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % zcd10 % cd10der1( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % zcd10 % cd10der1( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 


subroutine d2dx2_cd10(this, f, df) 
    
    class( derivative ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%nx) :: wrkArr
    integer :: i, j, k

    select case (this%xoprank)
    case (1)
        if (.not. cd10UseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % xcd10 % cd10der2( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % xcd10 % cd10der2( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. cd10UseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % xcd10 % cd10der2( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % xcd10 % cd10der2( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. cd10UseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % xcd10 % cd10der2( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % xcd10 % cd10der2( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 

subroutine d2dy2_cd10(this, f, df) 
    
    class( derivative ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%ny) :: wrkArr
    integer :: i, j, k

    select case (this%yoprank)
    case (1)
        if (.not. cd10UseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % ycd10 % cd10der2( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % ycd10 % cd10der2( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. cd10UseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % ycd10 % cd10der2( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % ycd10 % cd10der2( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. cd10UseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % ycd10 % cd10der2( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % ycd10 % cd10der2( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 


subroutine d2dz2_cd10(this, f, df) 
    
    class( derivative ), intent(in) :: this
    real(rkind), dimension(:,:,:), intent(in) :: f
    real(rkind), dimension(:,:,:), intent(out) :: df
    real(rkind), dimension(this%nz) :: wrkArr
    integer :: i, j, k

    select case (this%zoprank)
    case (1)
        if (.not. cd10UseWrkArr1) then
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    df(:, j, k) = this % zcd10 % cd10der2( f(:, j, k) ) 
                end do 
            end do 
        else
            do k = 1,size(f,3)
                do j = 1,size(f,2)
                    wrkArr = f(:, j, k)
                    df(:, j, k) = this % zcd10 % cd10der2( wrkArr ) 
                end do 
            end do 
        end if 

    case (2)
        if (.not. cd10UseWrkArr2) then
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    df(i, :, k) = this % zcd10 % cd10der2( f(i, :, k) ) 
                end do 
            end do 
        else 
            do k = 1,size(f,3)
                do i = 1,size(f,1)
                    wrkArr = f(i, :, k)
                    df(i, :, k) = this % zcd10 % cd10der2( wrkArr ) 
                end do 
            end do 
        end if 

    case (3)
        if (.not. cd10UseWrkArr3) then
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % zcd10 % cd10der2( f(i, j, :) ) 
                end do 
            end do 
        else
            do j = 1,size(f,2)
                do i = 1,size(f,1)
                    df(i, j, :) = this % zcd10 % cd10der2( f(i, j, :) ) 
                end do 
            end do 
        end if

    end select 

end subroutine 


