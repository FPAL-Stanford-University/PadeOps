subroutine renderVelocityXYperiodic()
  use domainSetup, only: dxF, dyF, dzF, xFh, yFh, zFh, nzF, nxsupp, nysupp, &
    nzsupp, gpF
  use exits, only: message
  integer :: n
  real(rkind), dimension(:,:,:), allocatable :: u, v, w
  real(rkind) :: small = 1.d-14
  character(len=clen) :: mssg
  real(rkind) :: kdotx, kdotx2, kdotx3
  integer :: i, j, k, iist, iien, jjst, jjen, kkst, kken
  integer :: istF, ienF, jstF, jenF, kstF, kenF, iszF, jszF, kszF
  real(rkind) :: cs, ss, fx, fy, fz, f
  real(rkind) :: wxSupport, wySupport, wzSupport

  ! Deine global indices for the halo-padded arrays (i.e. xFh, u, etc.)
  call getStartAndEndIndices(gpF,istF,ienF,jstF,jenF,kstF,kenF,iszF,jszF,kszF)
  istF = istF - nxsupp/2; ienF = ienF + nxsupp/2
  jstF = jstF - nysupp/2; jenF = jenF + nysupp/2
  kstF = max(1,kstF-nzsupp/2)
  kenF = min(nzF+1,kenF+nzsupp/2)

  allocate(u(istF:ienF,jstF:jenF,kstF:kenF))
  allocate(v(istF:ienF,jstF:jenF,kstF:kenF))
  allocate(w(istF:ienF,jstF:jenF,kstF:kenF))

  wxSupport = real(nxsupp+1,rkind)*dxF
  wySupport = real(nysupp+1,rkind)*dyF
  wzSupport = real(nzsupp+1,rkind)*dzF

  do n = 1,nmodes
    ! NOTE: These are global indices
    iist = ceiling((gmxloc(n)+small)/dxF) - nxsupp/2
    iien = floor((gmxloc(n)+small)/dxF) + nxsupp/2

    jjst = ceiling((gmyloc(n)+small)/dyF) - nxsupp/2
    jjen = floor((gmyloc(n)+small)/dyF) + nysupp/2

    kkst = max(1,ceiling((gmzloc(n)+small)/dzF) - nzsupp/2)
    kken = min(nzF,floor((gmzloc(n)+small)/dzF) + nzsupp/2)

    do k = kkst,kken
      kdotx3 = kz(n)*(zFh(k)-gmzloc(n))
      fz = cos(pi*(zFh(k)-gmzloc(n))/wzSupport)
      do j = jjst,jjen
        kdotx2 = ky(n)*(yFh(j)-gmyloc(n))
        fy = cos(pi*(yFh(j)-gmyloc(n))/wySupport)
        do i = iist,iien
          kdotx = kdotx2 + kdotx3 + kz(n)*(xFh(i)-gmxloc(n))

          cs = cos(kdotx)
          ss = sin(kdotx)

          fx = cos(pi*(xFh(i)-gmxloc(n))/wxSupport)
          f = fx*fy*fz

          u(i,j,k) = u(i,j,k) + f*(2.d0*uhatR(n)*cs - 2.d0*uhatI(n)*ss)
          v(i,j,k) = v(i,j,k) + f*(2.d0*vhatR(n)*cs - 2.d0*vhatI(n)*ss)
          w(i,j,k) = w(i,j,k) + f*(2.d0*whatR(n)*cs - 2.d0*whatI(n)*ss)
        end do
      end do
    end do
    if (mod(n,1000) == 0) then
      write(mssg,'(F7.4,A10)')real(n,rkind)/real(nmodes,rkind)*100.d0,'% Complete'
      call message(trim(mssg))
    end if
  end do

  !TODO: Add periodic contribution
  !call update_halo(

  deallocate(u,v,w)
end subroutine

subroutine meshgrid(x,y,z,x1D,y1D,z1D)
  real(rkind), dimension(:,:,:), intent(inout) :: x, y, z
  real(rkind), dimension(:), intent(in) :: x1D, y1D, z1D
  integer :: i, j, k
  call assert(size(x1D) == size(x,1),'Mismatched array sizes - meshgrid')
  call assert(size(y1D) == size(y,1),'Mismatched array sizes - meshgrid')
  call assert(size(z1D) == size(z,1),'Mismatched array sizes - meshgrid')

  do i = 1,size(x1D)
    x(i,:,:) = x1D(i)
  end do
  do j = 1,size(y1D)
    y(:,j,:) = y1D(j)
  end do
  do k = 1,size(z1D)
    z(:,:,k) = z1D(k)
  end do
end subroutine
