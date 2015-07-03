! Routines for Chebyshev derivatives
! Assumes Gauss-Lobatto points

function InitDCT() result(ierr)
    integer :: ierr
    ierr = xdct % init(nx)
    ierr = ydct % init(ny)
    ierr = zdct % init(nz)

    ierr = 0
end function


pure function chebder(fhat) result(dfhat)
    real(rkind), dimension(:), intent(in) :: fhat
    real(rkind), dimension(size(fhat)) :: dfhat
    integer :: n, i  

    n = size(fhat) - 1
    dfhat(n+1) = zero
    dfhat((n-1) + 1) = two*n*fhat(n+1)
    do i = (n - 1),2,-1
        dfhat((i - 1) + 1) = two*i*fhat(i + 1) + dfhat((i+1)+1)
    end do 

    dfhat(1) = fhat(2) + half*fhat(3) 

end function 


function d1xCHEB(f), result(df)
    real(rkind), dimension(nx,ny,nz), intent(in) :: f
    real(rkind), dimension(nx,ny,nz) :: df
    integer :: i, j

     do j = 1,nz
        do i = 1,ny 
            df(:,i,j) = xdct % idct( chebder(xdct % dct(f(:,i,j)))
        end do 
     end do 
end function

function d1yCHEB(f), result(df)
    real(rkind), dimension(nx,ny,nz), intent(in) :: f
    real(rkind), dimension(nx,ny,nz) :: df
    integer :: i, j

     do j = 1,nz
        do i = 1,nx 
            df(i,:,j) = ydct % idct( chebder(ydct % dct(f(i,:,j)))
        end do 
     end do 
end function

function d1zCHEB(f), result(df)
    real(rkind), dimension(nx,ny,nz), intent(in) :: f
    real(rkind), dimension(nx,ny,nz) :: df
    integer :: i, j

     do j = 1,ny
        do i = 1,nx 
            df(i,j,:) = zdct % idct( chebder(zdct % dct(f(i,j,:)))
        end do 
     end do 
end function

function d2xCHEB(f), result(df)
    real(rkind), dimension(nx,ny,nz), intent(in) :: f
    real(rkind), dimension(nx,ny,nz) :: df
    integer :: i, j

     do j = 1,nz
        do i = 1,ny 
            df(:,i,j) = xdct % idct( chebder ( chebder(xdct % dct(f(:,i,j))))
        end do 
     end do 
end function

function d2yCHEB(f), result(df)
    real(rkind), dimension(nx,ny,nz), intent(in) :: f
    real(rkind), dimension(nx,ny,nz) :: df
    integer :: i, j

     do j = 1,nz
        do i = 1,nx 
            df(i,:,j) = ydct % idct( chebder ( chebder(ydct % dct(f(i,:,j))))
        end do 
     end do 
end function

function d2zCHEB(f), result(df)
    real(rkind), dimension(nx,ny,nz), intent(in) :: f
    real(rkind), dimension(nx,ny,nz) :: df
    integer :: i, j

     do j = 1,ny
        do i = 1,nx 
            df(i,j,:) = zdct % idct( chebder ( chebder(zdct % dct(f(i,j,:))))
        end do 
     end do 
end function

function InitCHEB(direction) result (ierr)
    character(len=1), intent(in) :: direction
    integer :: ierr

    ierr = 0
end function
