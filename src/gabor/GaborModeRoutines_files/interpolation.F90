subroutine interpToLocation(datIn,datOut,gp,dx,dy,dz,x,y,z)
  ! Tri-linear interpolation to a given location (x,y,z)
  ! Inputs:
  !     datIn    --> 3D array of input data
  !     gp       --> Grid partition corresponding to datIn
  !     dx,dy,dz --> Mesh spacing corresponding to datIn
  !     x,y,z    --> Coordinates of interpolation location
  ! Output:
  !     datOut --> Scalar value of interpolated data at (x,y,z)
  use decomp_2d, only: decomp_info
  real(rkind), dimension(:,:,:), intent(in) :: datIn
  real(rkind), intent(out) :: datOut
  type(decomp_info), intent(in) :: gp
  real(rkind), intent(in) :: dx, dy, dz, x, y, z
  integer :: xlo, xhi, ylo, yhi, zlo, zhi
  real(rkind), dimension(8) :: weights
  integer, dimension(8) :: xid, yid, zid
  integer :: ist, jst, kst
  real(rkind) :: xst, yst, zst

  ist = gp%xst(1)
  jst = gp%xst(2)
  kst = gp%xst(3)

  xst = (real(ist,rkind)-1.d0)*dx
  yst = (real(jst,rkind)-1.d0)*dy
  zst = (real(kst,rkind)-1.d0)*dz

  call getXYZneighbors(x,y,z,ist,jst,kst,dx,dy,dz,xlo,xhi,ylo,yhi,zlo,zhi)
  call getWeightsForLinInterp(x,y,z,xlo,xhi,ylo,yhi,zlo,zhi,&
    xst,yst,zst,dx,dy,dz,weights,xid,yid,zid)
  !call getWeightsForLinInterp(idx,xlo,xhi,ylo,yhi,zlo,zhi,&
  !  x(xlo),y(ylo),z(zlo),dx,dy,dz,weights,xid,yid,zid)
  call interpDat(datOut,datIn,weights,xid,yid,zid)
end subroutine

subroutine interpDat(datOut, datIn, weights, xid, yid, zid)
    real(rkind), intent(out) :: datOut
    real(rkind), dimension(:,:,:), intent(in) :: datIn
    real(rkind), dimension(8), intent(in) :: weights
    integer :: idx
    integer, dimension(8), intent(in) :: xid, yid, zid
    
    ! Get interpolated value
    datOut  = weights(1)*datIn(xid(1),yid(1),zid(1))
    do idx = 2,8
      datOut  = datOut  + weights(idx)*datIn(xid(idx),yid(idx),zid(idx))
    end do
end subroutine

pure subroutine getWeightsForLinInterp(x, y, z, xlo, xhi, ylo, yhi, zlo, zhi,&
                                    xLESlo, yLESlo, zLESlo, dx, dy, dz, weights, &
                                    xid, yid, zid)
    ! Computes the weights for tri-linear interplation
    ! Inputs:
    !   x, y, z --> Physical location of interpolated (output) data
    !   xlo, xhi, etc. --> The LES edge mesh index for the neighboring LES nodes
    !                      (neighbor to Gabor mode location)
    !   xLESlo, etc. --> Physical location of LES neighbor to left
    !
    ! Outputs:
    !   weights --> Weights for interpolation
    !   xid, etc. --> arrays corresponding to node indices of neighbors

    real(rkind), intent(in) :: x, y, z, xLESlo, yLESlo, zLESlo, dx, dy, dz
    real(rkind), dimension(8), intent(out) :: weights
    real(rkind) :: xweight, yweight, zweight
    integer, intent(in) :: xlo, xhi, ylo, yhi, zlo, zhi
    integer, dimension(8), intent(out) :: xid, yid, zid

    ! Step2: assemble xid, yid, zid
    xid(1) = xlo; xid(2) = xlo; xid(3) = xlo; xid(4) = xlo
    xid(5) = xhi; xid(6) = xhi; xid(7) = xhi; xid(8) = xhi
    
    yid(1) = ylo; yid(2) = ylo; yid(3) = yhi; yid(4) = yhi
    yid(5) = ylo; yid(6) = ylo; yid(7) = yhi; yid(8) = yhi
    
    zid(1) = zlo; zid(2) = zhi; zid(3) = zlo; zid(4) = zhi
    zid(5) = zlo; zid(6) = zhi; zid(7) = zlo; zid(8) = zhi

    ! Step3: Compute weights for each of the 8 points
    xweight = (x - xLESlo)/dx
    yweight = (y - yLESlo)/dy
    zweight = (z - zLESlo)/dz

    weights(1) = (1-yweight) * (1-zweight) * (1-xweight)
    weights(2) = (1-yweight) * zweight     * (1-xweight)

    weights(3) = yweight     * (1-zweight) * (1-xweight)
    weights(4) = yweight     * zweight     * (1-xweight)

    weights(5) = (1-yweight) * (1-zweight) * xweight
    weights(6) = (1-yweight) * zweight     * xweight

    weights(7) = yweight     * (1-zweight) * xweight
    weights(8) = yweight     * zweight     * xweight

end subroutine
  
pure subroutine getXYZneighbors(x,y,z,ist,jst,kst,dx,dy,dz,xlo,xhi,ylo,yhi,zlo,zhi)
  ! Get indices of neighboring LES edge nodes 
  ! Inputs:
  !     x, y, z       --> Physical location of desired interoplated data
  !     dx, dy, dz    --> Mesh spacing of input data
  !     ist, jst, kst --> Global indices for the start of MPI partition arrays
  ! Outputs:
  !     xlo, xhi, etc. --> Local indices of arrays relative to MPI partition

    integer, intent(in) :: ist, jst, kst
    real(rkind), intent(in) :: dx, dy, dz, x, y, z
    integer, intent(out) :: xlo, xhi, ylo, yhi, zlo, zhi
    
    ! Step1: compute xlo, xhi, ylo, etc... for the 8 nearest points

    xlo = floor(x/dx) + 1 - ist + 1
    xhi = xlo + 1
    
    ylo = floor(y/dy) + 1 - jst + 1
    yhi = ylo + 1
    
    zlo = floor(z/dz) + 1 - kst + 1
    zhi = zlo + 1
    
end subroutine 

integer function findMeshIdx(xloc,x)
  ! Returns the index of cell x where xloc is found

  real(rkind), intent(in) :: xloc
  real(rkind), dimension(:), intent(in) :: x
  integer :: i, N

  N = size(x)-1
  findMeshIdx = 0
  do i = 1,N
    if (xloc > x(i) .and. xloc <= x(i+1)) then
      findMeshIdx = i
    end if
  end do
end function

subroutine getNearestNeighborValue(datIn, datOut, gp, dx, dy, dz, x, y, z)
  use decomp_2d, only: decomp_info
  real(rkind), dimension(:,:,:), intent(in) :: datIn 
  real(rkind), intent(out) :: datOut
  type(decomp_info), intent(in) :: gp
  real(rkind), intent(in) :: dx, dy, dz, x, y, z
  integer, dimension(8) :: xid, yid, zid
  real(rkind), dimension(8) :: weights
  integer :: ist, jst, kst
  integer :: xlo, xhi, ylo, yhi, zlo, zhi, idx, neighIdx
  real(rkind) :: minDist, dist, xst, yst, zst, xneigh, yneigh, zneigh 

  ist = gp%xst(1)
  jst = gp%xst(2)
  kst = gp%xst(3)
  
  xst = (real(ist,rkind)-1.d0)*dx
  yst = (real(jst,rkind)-1.d0)*dy
  zst = (real(kst,rkind)-1.d0)*dz

  call getXYZneighbors(x,y,z,ist,jst,kst,dx,dy,dz,xlo,xhi,ylo,yhi,zlo,zhi)
  call getWeightsForLinInterp(x,y,z,xlo,xhi,ylo,yhi,zlo,zhi,&
    xst,yst,zst,dx,dy,dz,weights,xid,yid,zid)

  minDist = 1.d12
  neighIdx = 0
  do idx = 1,8
    xneigh = dx*(xid(idx) - 1)
    yneigh = dy*(yid(idx) - 1)
    zneigh = dz*(zid(idx) - 1)

    dist = sqrt((xneigh - x)**2.d0 + (yneigh - y)**2.d0 + (zneigh - z)**2.d0)
    if (dist < minDist) then
      minDist = dist 
      neighIdx = idx 
    end if 
  end do 

  datOut = datIn(xid(neighIdx), yid(neighIdx), zid(neighIdx))
end subroutine