subroutine interpToMode(datIn,datOut,idx,gp,dx,dy,dz,x,y,z)
  real(rkind), dimension(:,:,:), intent(in) :: datIn
  real(rkind), intent(out) :: datOut
  integer, intent(in) :: idx
  type(decomp_info), intent(in) :: gp
  real(rkind), intent(in) :: dx, dy, dz 
  real(rkind), dimension(:), intent(in) :: x, y, z
  integer :: xlo, xhi, ylo, yhi, zlo, zhi
  real(rkind), dimension(8) :: weights
  integer, dimension(8) :: xid, yid, zid
  integer :: ist, ien, jst, jen, kst, ken, isz, jsz, ksz

  call getStartAndEndIndices(gp,ist,ien,jst,jen,kst,ken,isz,jsz,ksz)
  call getXYZneighbors(idx,ist,jst,kst,dx,dy,dz,xlo,xhi,ylo,yhi,zlo,zhi)
  call getWeightsForLinInterp(idx,xlo,xhi,ylo,yhi,zlo,zhi,&
    x(xlo),y(ylo),z(zlo),dx,dy,dz,weights,xid,yid,zid)
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

pure subroutine getWeightsForLinInterp(gmID, xlo, xhi, ylo, yhi, zlo, zhi,&
                                    xLESlo, yLESlo, zLESlo, dx, dy, dz, weights, &
                                    xid, yid, zid)
    ! Computes the weights for tri-linear interplation
    ! Inputs:
    !   gmID --> Gabor mode index
    !   xlo, xhi, etc. --> The LES edge mesh index for the neighboring LES nodes
    !                      (neighbor to Gabor mode location)
    !   xLESlo, etc. --> Physical location of LES neighbor to left
    !
    ! Outputs:
    !   weights --> Weights for interpolation
    !   xid, etc. --> arrays corresponding to node indices of neighbors

    real(rkind), intent(in) :: xLESlo, yLESlo, zLESlo, dx, dy, dz
    real(rkind), dimension(8), intent(out) :: weights
    real(rkind) :: xweight, yweight, zweight, onebydx, onebydy, onebydz
    integer, intent(in) :: gmID, xlo, xhi, ylo, yhi, zlo, zhi
    integer, dimension(8), intent(out) :: xid, yid, zid

    ! Step2: assemble xid, yid, zid
    xid(1) = xlo; xid(2) = xlo; xid(3) = xlo; xid(4) = xlo;
    xid(5) = xhi; xid(6) = xhi; xid(7) = xhi; xid(8) = xhi;
    
    yid(1) = ylo; yid(2) = ylo; yid(3) = yhi; yid(4) = yhi;
    yid(5) = ylo; yid(6) = ylo; yid(7) = yhi; yid(8) = yhi;
    
    zid(1) = zlo; zid(2) = zhi; zid(3) = zlo; zid(4) = zhi;
    zid(5) = zlo; zid(6) = zhi; zid(7) = zlo; zid(8) = zhi;

    ! Step3: Compute weights for each of the 8 points
    onebydx = 1/dx; onebydy = 1/dy; onebydz = 1/dz
    xweight = (gmxloc(gmID) - xLESlo)*onebydx
    yweight = (gmyloc(gmID) - yLESlo)*onebydy
    zweight = (gmzloc(gmID) - zLESlo)*onebydz

    weights(1) = (1-yweight) * (1-zweight) * (1-xweight);
    weights(2) = (1-yweight) * zweight     * (1-xweight);

    weights(3) = yweight     * (1-zweight) * (1-xweight);
    weights(4) = yweight     * zweight     * (1-xweight);

    weights(5) = (1-yweight) * (1-zweight) * xweight;
    weights(6) = (1-yweight) * zweight     * xweight;

    weights(7) = yweight     * (1-zweight) * xweight;
    weights(8) = yweight     * zweight     * xweight;

end subroutine
  
pure subroutine getXYZneighbors(idx,ist,jst,kst,dx,dy,dz,xlo,xhi,ylo,yhi,zlo,zhi)
  ! Get indices of neighboring LES edge nodes 
  ! Inputs:
  !     idx --> Gabor mode index
  !     ist, jst, kst --> Global indices for the start of MPI partition arrays
  ! Outputs:
  !     xlo, xhi, etc. --> Local indices of arrays relative to MPI partition

    integer, intent(in) :: idx, ist, jst, kst
    real(rkind), intent(in) :: dx, dy, dz
    integer, intent(out) :: xlo, xhi, ylo, yhi, zlo, zhi
    real(rkind) :: onebydx, onebydy, onebydz
    
    ! Step1: compute xlo, xhi, ylo, etc... for the 8 nearest points
    onebydx = 1.d0/dx
    onebydy = 1.d0/dy
    onebydz = 1.d0/dz

    xlo = floor(  gmxloc(idx)*onebydx) + 1 - ist + 1
    xhi = xlo + 1
    !xhi = ceiling(gmxloc(idx)*onebydx) + 1 - ist + 1
    
    ylo = floor(  gmyloc(idx)*onebydy) + 1 - jst + 1
    yhi = ylo + 1
    !yhi = ceiling(gmyloc(idx)*onebydy) + 1 - jst + 1
    
    zlo = floor(  gmzloc(idx)*onebydz) + 1 - kst + 1
    zhi = zlo + 1
    !zhi = ceiling(gmzloc(idx)*onebydz) + 1 - kst + 1
    
end subroutine 
