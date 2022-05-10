subroutine renderVelocityXYperiodic()
  use exits, only: message
  integer :: n
  !real(rkind), dimension(:,:,:), allocatable :: buff
  real(rkind) :: small = 1.d-14
  character(len=clen) :: mssg
  real(rkind) :: kdotx, kdotx2, kdotx3
  integer :: i, j, k, iist, iien, jjst, jjen, kkst, kken
  !integer :: iisz, jjsz, kksz
  !integer :: ii, jj, kk
  real(rkind) :: cs, ss, fx, fy, fz, f
  real(rkind) :: wxSupport, wySupport, wzSupport

  call assert(isotropicModesInitialized,"Gabor modes aren't initialized. "//&
    "Cannot render the velocity field")
  call message("Rendering the Gabor-inducded velocity field")
  
  wxSupport = real(nxsupp+1,rkind)*dxF
  wySupport = real(nysupp+1,rkind)*dyF
  wzSupport = real(nzsupp+1,rkind)*dzF
  
  do n = 1,nmodes  
    ! NOTE: These are global indices of the physical domain
    iist =           ceiling((gmxloc(n)+small)/dxF) - nxsupp/2
    iien =           floor(  (gmxloc(n)+small)/dxF) + nxsupp/2

    jjst =           ceiling((gmyloc(n)+small)/dyF) - nysupp/2
    jjen =           floor(  (gmyloc(n)+small)/dyF) + nysupp/2

    kkst = max(1    ,ceiling((gmzloc(n)+small)/dzF) - nzsupp/2)
    kken = min(nzF+1,floor(  (gmzloc(n)+small)/dzF) + nzsupp/2)

    do k = kkst,kken
      kdotx3 = kz(n)*(zFh(k)-gmzloc(n))
      fz = cos(pi*(zFh(k)-gmzloc(n))/wzSupport)
      do j = jjst,jjen
        kdotx2 = ky(n)*(yFh(j)-gmyloc(n))
        fy = cos(pi*(yFh(j)-gmyloc(n))/wySupport)
        do i = iist,iien
          kdotx = kdotx2 + kdotx3 + kx(n)*(xFh(i)-gmxloc(n))

          cs = cos(kdotx)
          ss = sin(kdotx)

          fx = cos(pi*(xFh(i)-gmxloc(n))/wxSupport)
          f = fx*fy*fz

           !call globalToGlobalIdx(i,j,k,ii,jj,kk)
           !uFh(ii,jj,kk) = uFh(ii,jj,kk) + f*(2.d0*uhatR(n)*cs - 2.d0*uhatI(n)*ss)
           !vFh(ii,jj,kk) = vFh(ii,jj,kk) + f*(2.d0*vhatR(n)*cs - 2.d0*vhatI(n)*ss)
           !wFh(ii,jj,kk) = wFh(ii,jj,kk) + f*(2.d0*whatR(n)*cs - 2.d0*whatI(n)*ss)
           uFh(i,j,k) = uFh(i,j,k) + f*(2.d0*uhatR(n)*cs - 2.d0*uhatI(n)*ss)
           vFh(i,j,k) = vFh(i,j,k) + f*(2.d0*vhatR(n)*cs - 2.d0*vhatI(n)*ss)
           wFh(i,j,k) = wFh(i,j,k) + f*(2.d0*whatR(n)*cs - 2.d0*whatI(n)*ss)
        end do
      end do
    end do
    if (mod(n,1000) == 0) then
      write(mssg,'(F7.4,A10)')real(n,rkind)/real(nmodes,rkind)*100.d0,'% Complete'
      call message(trim(mssg))
    end if
  end do
  
  ! Add periodic contribution
!  call assert(nxsupp == nysupp,'Window supports must be isotropic for "//&
!    MPI halo exchange in renderVelocityXYperiodic -- GaborRoutines.F90')
!  call assert(nxsupp == nzsupp,'Window supports must be isotropic for "//&
!    MPI halo exchange in renderVelocityXYperiodic -- GaborRoutines.F90')
!  call assert(nysupp == nzsupp,'Window supports must be isotropic for "//&
!    MPI halo exchange in renderVelocityXYperiodic -- GaborRoutines.F90')
!  call getStartAndEndIndices(gpFh,iist,iien,jjst,jjen,kkst,kken,iisz,jjsz,kksz)
!  call update_halo(uFh,buff,level=nxsupp/2,opt_decomp=gpFh,opt_global=.true.)
!  call haloExchange(uFh,buff)
!  
!  call update_halo(vFh,buff,level=nxsupp/2,opt_decomp=gpFh,opt_global=.true.)
!  call haloExchange(vFh,buff)
!  
!  call update_halo(wFh,buff,level=nxsupp/2,opt_decomp=gpFh,opt_global=.true.)
!  call haloExchange(wFh,buff)
  
!  deallocate(buff)
  
  ! TODO: Uses custom implementation of halo_exchange. This does not require isotropic
  ! window supports  
  call haloExchangeMPI(uFh)
  call haloExchangeMPI(vFh)
  call haloExchangeMPI(wFh)

  ! Copy data to final velocity arrays
!  call getStartAndEndIndices(gpFh,iist,iien,jjst,jjen,kkst,kken,iisz,jjsz,kksz)
!  if (kkst == 1 .and. kken == nzF+1) then
!    uG = uFh(iist+nxsupp/2:iien-nxsupp/2,jjst+nysupp/2:jjen-nysupp/2,:)
!    vG = vFh(iist+nxsupp/2:iien-nxsupp/2,jjst+nysupp/2:jjen-nysupp/2,:)
!    wG = wFh(iist+nxsupp/2:iien-nxsupp/2,jjst+nysupp/2:jjen-nysupp/2,:)
!  elseif (kkst == 1) then
!    uG = uFh(iist+nxsupp/2:iien-nxsupp/2,jjst+nysupp/2:jjen-nysupp/2,1:kken-nzsupp/2)
!    vG = vFh(iist+nxsupp/2:iien-nxsupp/2,jjst+nysupp/2:jjen-nysupp/2,1:kken-nzsupp/2)
!    wG = wFh(iist+nxsupp/2:iien-nxsupp/2,jjst+nysupp/2:jjen-nysupp/2,1:kken-nzsupp/2)
!  elseif (kken == nzF+1) then
!    uG = uFh(iist+nxsupp/2:iien-nxsupp/2,jjst+nysupp/2:jjen-nysupp/2,kkst+nzsupp/2:kken)
!    vG = vFh(iist+nxsupp/2:iien-nxsupp/2,jjst+nysupp/2:jjen-nysupp/2,kkst+nzsupp/2:kken)
!    wG = wFh(iist+nxsupp/2:iien-nxsupp/2,jjst+nysupp/2:jjen-nysupp/2,kkst+nzsupp/2:kken)
!  else
!    uG = uFh(iist+nxsupp/2:iien-nxsupp/2,jjst+nysupp/2:jjen-nysupp/2,kkst+nzsupp/2:kken-nzsupp/2)
!    vG = vFh(iist+nxsupp/2:iien-nxsupp/2,jjst+nysupp/2:jjen-nysupp/2,kkst+nzsupp/2:kken-nzsupp/2)
!    wG = wFh(iist+nxsupp/2:iien-nxsupp/2,jjst+nysupp/2:jjen-nysupp/2,kkst+nzsupp/2:kken-nzsupp/2)
!  end if
  
  !TODO: Interpolate velocity to dzF/2:dzF:Lz-dzF/2 mesh
  call message("WARNING: Gabor-induced velocity has not been interpolated "//&
    "to cell centers")
  
end subroutine

subroutine haloExchange(u,buff)
  real(rkind), dimension(:,:,:), intent(inout) :: u
  real(rkind), dimension(:,:,:), intent(in) :: buff
  integer :: iszU, jszU, kszU, iszB, jszB, kszB
  
  iszU = size(u,1)
  jszU = size(u,2)
  kszU = size(u,3)

  ! For x-pencil: shape(buff) = iszU        X jszU+nysupp X kszU+nzsupp
  ! For y-pencil: shape(buff) = iszU+nxsupp X jszU        X kszU+nzsupp
  ! For z-pencil: shape(buff) = iszU+nxsupp X jszU+nzsupp X kszU
  iszB = size(buff,1)
  jszB = size(buff,2)
  kszB = size(buff,3)
call assert(maxval(abs(u - buff(:,nysupp/2+1:jszB-nysupp/2,nzsupp/2+1:kszB-nzsupp/2)))<1.d-12,'A')!DEBUG

print*, "u(1,1,1):", u(1,1,1) !DEBUG
print*, "buff(1,1+nysupp/2,1+nzsupp/2):", buff(1,1+nysupp/2,1+nzsupp/2) !DEBUG
print*, "u(1,2,1):", u(1,2,1) !DEBUG
print*, "buff(1,jszB-nysupp/2+1,1+nzsupp/2):", buff(1,jszB-nysupp/2+1,1+nzsupp/2) !DEBUG
print*, "buff(1,jszB,1+nzsupp/2):", buff(1,jszB,1+nzsupp/2) !DEBUG
print*, "u(1,jszU,1):", u(1,jszU,1) !DEBUG
print*, "buff(1,nysupp/2,1+nzsupp/2):", buff(1,nysupp/2,1+nzsupp/2) !DEBUG
call assert(.false.)
print*, maxval(abs(u(:,1,:) - buff(:,jszB-nysupp/2+1,nzsupp/2+1:kszB-nzsupp/2))) !DEBUG
call assert(maxval(abs(u(:,1,:) - buff(:,jszB-nysupp/2+1,nzsupp/2+1:kszB-nzsupp/2))) < 1.d-12,'C') !DEBUG

print*, maxval(abs(buff(:,nysupp/2,1+nzsupp/2:kszB-nzsupp/2) - u(:,jszU-nysupp/2,:))) !DEBUG
call assert(maxval(abs(buff(:,nysupp/2,1+nzsupp/2:kszB-nzsupp/2) - u(:,jszU-nysupp/2,:)))<1.d-14,'B') !DEBUG
  select case (decomp2Dpencil)
  case ('x') 
     if (periodic(1)) then
       ! Left halo exchange 
       u(1+nxsupp/2:nxsupp,:,:) = u(1+nxsupp/2:nxsupp,:,:) + &
         u(iszU-nxsupp/2+1:iszU,:,:)
       u(:,1+nysupp/2:nysupp,:) = u(:,1+nysupp/2:nysupp,:) + &
         buff(:,1:nysupp/2,1+nzsupp/2:kszB-nzsupp/2)
       u(:,:,1+nzsupp/2:nzsupp) = u(:,:,1+nzsupp/2:nzsupp) + &
         buff(:,1+nysupp/2:jszB-nysupp/2,1:nzsupp/2)

       ! Right halo exchange
       u(iszU-nxsupp+1:iszU-nxsupp/2,:,:) = u(iszU-nxsupp+1:iszU-nxsupp/2,:,:) + &
         u(1:nxsupp/2,:,:)
       u(:,jszU-nysupp+1:jszU-nysupp/2,:) = u(:,jszU-nysupp+1:jszU-nysupp/2,:) + &
         buff(:,jszB-nysupp/2+1:jszB,1+nzsupp/2:kszB-nzsupp/2)
       u(:,:,kszU-nzsupp+1:kszU-nzsupp/2) = u(:,:,kszU-nzsupp+1:kszU-nzsupp/2) + &
         buff(:,1+nysupp/2:jszB-nysupp/2,kszB-nzsupp/2+1:kszB)
     else
       call assert(.false.,'x direction assumed periodic. Other BCs not'//&
         ' supported -- GaborModeRoutines.F90')
     end if
  case ('y')
    call assert(.false.,'y pencil not implemented yet')
  case ('z')
    call assert(.false.,'z pencil not implemented yet')
  end select
end subroutine

pure subroutine globalToGlobalIdx(iin,jin,kin,iout,jout,kout)
  ! Converts a global index from the physical domain to a global index in the
  ! halo padded arrays
  ! Inputs:
  !     iin --> global index in physical domain. Ranges from [1-nxsupp/2:nx+nxsupp/2]
  !             Similar definitions for jin and kin
  ! Ouptputs:
  !     iout --> global index in halo padded array. Ranges from
  !              [1:nx + nxsupp + nprocX*nxsupp] for periodic directions
  !              and [1:nx + nprocX*nxsupp] for non-periodic directions.
  !              Similar definitions for jout and kout.

  integer, intent(in) :: iin, jin, kin
  integer, intent(out) :: iout, jout, kout
  iout = iin + nrankX*nxsupp + nxsupp/2
  jout = jin + nrankY*nysupp + nysupp/2
!  kout = kin + nrankZ*nzsupp + nzsupp/2
  
!  if (.not. periodic(1)) iout = iin + nrankX*nxsupp
!  if (.not. periodic(2)) jout = jin + nrankY*nysupp
!  if (.not. periodic(3)) kout = kin + nrankZ*nzsupp
  kout = kin + nrankZ*nzsupp
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
