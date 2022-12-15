module interpolatorMod
    use kind_parameters, only: rkind
    use decomp_2d
    use exits, only: GracefulExit 

    implicit none 
    private 
    public :: interpolator

    type :: interpolator
        private  
        real(rkind), dimension(:,:,:), allocatable :: fx_X, fx_Y, fxy_Y, fxy_Z, fxyz_Z, fxyz_Y
        type(decomp_info), pointer :: gpSource, gpDest 
        type(decomp_info) :: gpSX, gpSXY
        integer, dimension(:), allocatable  :: xInd, yInd, zInd  
        real(rkind), dimension(:), allocatable :: wx, wy, wz
        contains
            procedure :: init 
            procedure :: destroy 
            procedure :: LinInterp3D 
    end type 

contains 

subroutine init(this, gpSource, gpDest, xSource, ySource, zSource, xDest, yDest, zDest)
    class(interpolator), intent(inout) :: this
    type(decomp_info), intent(in), target :: gpSource, gpDest 
    real(rkind), dimension(:), intent(in) :: xSource, ySource, zSource, xDest, yDest, zDest
    integer :: nxS, nyS, nzS, nxD, nyD, nzD, idx 
    real(rkind) :: delta, start

    this%gpSource => gpSource
    this%gpDest => gpDest

    ! Safeguards
    !if (xSource(1) > xDest(1)) then
    !    call GracefulExit("Low bound of x-axis in out of bounds (interpolator)",34)
    !end if
    !if (xSource(size(xSource)) < xDest(size(xDest))) then
    !    print*, xSource(size(xSource)), xDest(size(xDest))
    !    call GracefulExit("High bound of x-axis in out of bounds (interpolator)",34)
    !end if  
    !if (ySource(1) > yDest(1)) then
    !    call GracefulExit("Low bound of y-axis in out of bounds (interpolator)",34)
    !end if
    !if (ySource(size(ySource)) < yDest(size(yDest))) then
    !    print*, ySource(size(ySource)),yDest(size(yDest)) 
    !    call GracefulExit("High bound of y-axis in out of bounds (interpolator)",34)
    !end if  
    !if (zSource(1) > zDest(1)) then
    !    call GracefulExit("Low bound of z-axis in out of bounds (interpolator)",34)
    !end if
    !if (zSource(size(zSource)) < zDest(size(zDest))) then
    !    call GracefulExit("High bound of z-axis in out of bounds (interpolator)",34)
    !end if  

    allocate(this%xInd(size(xDest)))
    allocate(this%wx(size(xDest)))
    
    allocate(this%yInd(size(yDest)))
    allocate(this%wy(size(yDest)))
    
    allocate(this%zInd(size(zDest)))
    allocate(this%wz(size(zDest)))

    nxS = gpSource%xsz(1)
    nyS = gpSource%ysz(2)
    nzS = gpSource%zsz(3)

    nxD = gpDest%xsz(1)
    nyD = gpDest%ysz(2)
    nzD = gpDest%zsz(3)
    
   
    ! Get interpolation indices and weights
    delta = xSource(2) - xSource(1)
    start = xSource(1)
    do idx = 1,size(this%wx)
        this%xInd(idx) = ceiling((xDest(idx) - start)/delta)
        if (this%xInd(idx) < 1) then 
            this%wx(idx) = 1.d0
            this%xInd(idx) = 1
        elseif (this%xInd(idx) > size(xSource) - 1) then
            this%wx(idx) = 0.d0
            this%xInd(idx) = size(xSource) - 1 
        else
            this%wx(idx) = (xSource(this%xInd(idx) + 1) - xDest(idx))/delta 
        endif  
    end do 

    delta = ySource(2) - ySource(1)
    start = ySource(1)
    do idx = 1,size(this%wy)
        this%yInd(idx) = ceiling((yDest(idx) - start)/delta)
        if (this%yInd(idx) < 1) then 
            this%wy(idx) = 1.d0
            this%yInd(idx) = 1
        elseif (this%yInd(idx) > size(ySource) - 1) then
            this%wy(idx) = 0.d0
            this%yInd(idx) = size(ySource) - 1 
        else
            this%wy(idx) = (ySource(this%yInd(idx) + 1) - yDest(idx))/delta 
        end if 
    end do 

    delta = zSource(2) - zSource(1)
    start = zSource(1)
    do idx = 1,size(this%wz)
        this%zInd(idx) = ceiling((zDest(idx) - start)/delta)
        if (this%zInd(idx) < 1) then 
            this%wz(idx) = 1.d0
            this%zInd(idx) = 1
        elseif (this%zInd(idx) > size(zSource) - 1) then
            this%wz(idx) = 0.d0
            this%zInd(idx) = size(zSource) - 1 
        else
            this%wz(idx) = (zSource(this%zInd(idx) + 1) - zDest(idx))/delta 
        end if 
    end do 

    ! Create 2 intermediate transposers and buffer arrays    
    call decomp_info_init(nxD,nyS,nzS,this%gpSX)
    call decomp_info_init(nxD,nyD,nzS,this%gpSXY)
    allocate(this%fx_X  (this%gpSX %xsz(1),this%gpSX %xsz(2),this%gpSX %xsz(3)))
    allocate(this%fx_Y  (this%gpSX %ysz(1),this%gpSX %ysz(2),this%gpSX %ysz(3)))
    allocate(this%fxy_Y (this%gpSXY%ysz(1),this%gpSXY%ysz(2),this%gpSXY%ysz(3)))
    allocate(this%fxy_Z (this%gpSXY%zsz(1),this%gpSXY%zsz(2),this%gpSXY%zsz(3)))
    allocate(this%fxyz_Z(this%gpDest%zsz(1),this%gpDest%zsz(2),this%gpDest%zsz(3)))
    allocate(this%fxyz_Y(this%gpDest%ysz(1),this%gpDest%ysz(2),this%gpDest%ysz(3)))

end subroutine 

subroutine destroy(this)
    class(interpolator), intent(inout) :: this
    deallocate(this%wx, this%wy, this%wz, this%xInd, this%yInd, this%zInd)
    deallocate(this%fx_X, this%fx_Y, this%fxy_Y, this%fxy_Z)
end subroutine

subroutine LinInterp3D(this, fS, fD, buffer)
    use constants, only: one 
    class(interpolator), intent(inout) :: this 
    real(rkind), dimension(this%gpSource%xsz(1),this%gpSource%xsz(2), this%gpSource%xsz(3)), intent(in) :: fS
    real(rkind), dimension(this%gpDest%xsz(1),this%gpDest%xsz(2), this%gpDest%xsz(3)), intent(inout) :: fD
    real(rkind), dimension(this%gpDest%xsz(1),this%gpDest%xsz(2), this%gpDest%xsz(3)), intent(out), optional :: buffer

    integer :: i, j, k 
    
    ! interpolate in x
    do k = 1,this%gpSX%xsz(3)
        do j = 1,this%gpSX%xsz(2)
            do i = 1,this%gpSX%xsz(1)
                this%fx_X(i,j,k) = this%wx(i)*fS(this%xInd(i),j,k) + (one - this%wx(i))*fS(this%xInd(i)+1,j,k)
            end do 
        end do 
    end do
    call transpose_x_to_y(this%fx_X, this%fx_Y, this%gpSX)

    ! interpolate in y 
    do k = 1,this%gpSXY%ysz(3)
        do j = 1,this%gpSXY%ysz(2)
            do i = 1,this%gpSXY%ysz(1)
                this%fxy_Y(i,j,k) = this%wy(j)*this%fx_Y(i,this%yInd(j),k) + (one - this%wy(j))*this%fx_Y(i,this%yInd(j)+1,k) 
            end do 
        end do 
    end do
    call transpose_y_to_z(this%fxy_Y,this%fxy_Z, this%gpSXY) 

    ! intepolate in z 
    do k = 1,this%gpDest%zsz(3)
        do j = 1,this%gpDest%zsz(2)
            do i = 1,this%gpDest%zsz(1)
                this%fxyz_Z(i,j,k) = this%wz(k)*this%fxy_Z(i,j,this%Zind(k)) + (one - this%wz(k))*this%fxy_Z(i,j,this%Zind(k)+1)
            end do 
        end do 
    end do 

    ! Finally, transpose back to x decomposition 
    call transpose_z_to_y(this%fxyz_Z,this%fxyz_Y,this%gpDest)
    !call transpose_y_to_x(this%fxyz_Y,fD, this%gpDest) ! DONE!

    if (present(buffer)) then 
        call transpose_y_to_x(this%fxyz_Y, buffer, this%gpDest) 
        fD = fD + buffer ! DONE!
    else
        call transpose_y_to_x(this%fxyz_Y,fD, this%gpDest) ! DONE!
    end if 

end subroutine 

end module 
