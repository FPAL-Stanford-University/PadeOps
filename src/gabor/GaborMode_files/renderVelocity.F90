subroutine renderVelocityXYperiodic()
  use exits, only: message
  use omp_lib, only: omp_get_thread_num, omp_get_num_threads
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
  integer :: tid
  real(rkind), dimension(:,:,:,:), allocatable :: utmp,vtmp,wtmp

  call assert(isotropicModesInitialized,"Gabor modes aren't initialized. "//&
    "Cannot render the velocity field")
  call message("Rendering the Gabor-induced velocity field")
  
  wxSupport = real(nxsupp+1,rkind)*dxF
  wySupport = real(nysupp+1,rkind)*dyF
  wzSupport = real(nzsupp+1,rkind)*dzF
   
  ! Allocate velocity arrays for each thread  
  allocate(utmp(istFh:ienFh,jstFh:jenFh,kstFh:kenFh,nthreads))
  allocate(vtmp(istFh:ienFh,jstFh:jenFh,kstFh:kenFh,nthreads))
  allocate(wtmp(istFh:ienFh,jstFh:jenFh,kstFh:kenFh,nthreads))
 
  !$OMP PARALLEL SHARED(utmp,vtmp,wtmp,uFh,vFh,wFh) &
  !$OMP PRIVATE(tid,n,i,j,k,kdotx,kdotx3,kdotx2,fx,fy,fz,f) &
  !$OMP PRIVATE(cs,ss,iist,iien,jjst,jjen,kkst,kken)
  tid = omp_get_thread_num()
  if (tid == 0 .and. nrank == 0) print*, "# of threads spawned:", omp_get_num_threads()
  !$OMP DO
  do n = 1,nmodes  
    ! NOTE: These are global indices of the physical domain
    iist =           ceiling((gmxloc(n)+small)/dxF) - nxsupp/2
    iien =           floor(  (gmxloc(n)+small)/dxF) + nxsupp/2

    jjst =           ceiling((gmyloc(n)+small)/dyF) - nysupp/2
    jjen =           floor(  (gmyloc(n)+small)/dyF) + nysupp/2

    kkst = max(1    ,ceiling((gmzloc(n)+small)/dzF) - nzsupp/2)
    kken = min(nzF+1,floor(  (gmzloc(n)+small)/dzF) + nzsupp/2)

    !!$OMP PARALLEL SHARED(uFh, vFh, wFh) PRIVATE(tid,i,j,k)
    !tid = omp_get_thread_num()
    !if (n == 1 .and. tid == 0) print*, "# of threads spawned:", omp_get_num_threads()
    !!$OMP DO
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
          
          utmp(i,j,k,tid+1) = utmp(i,j,k,tid+1) + f*(2.d0*uhatR(n)*cs - 2.d0*uhatI(n)*ss)
          vtmp(i,j,k,tid+1) = vtmp(i,j,k,tid+1) + f*(2.d0*vhatR(n)*cs - 2.d0*vhatI(n)*ss)
          wtmp(i,j,k,tid+1) = wtmp(i,j,k,tid+1) + f*(2.d0*whatR(n)*cs - 2.d0*whatI(n)*ss)
          !uFh(i,j,k) = uFh(i,j,k) + f*(2.d0*uhatR(n)*cs - 2.d0*uhatI(n)*ss)
          !vFh(i,j,k) = vFh(i,j,k) + f*(2.d0*vhatR(n)*cs - 2.d0*vhatI(n)*ss)
          !wFh(i,j,k) = wFh(i,j,k) + f*(2.d0*whatR(n)*cs - 2.d0*whatI(n)*ss)
        end do
      end do
    end do

    if (mod(n,100000) == 0 .and. tid == 0) then
      write(mssg,'(F7.4,A10)')real(n,rkind)/real(nmodes,rkind)*100.d0*32.d0,'% Complete'
      call message(trim(mssg))
    end if
  end do
  !$OMP END DO
  !$OMP CRITICAL
  uFh = uFh + utmp(:,:,:,tid+1)
  vFh = vFh + vtmp(:,:,:,tid+1)
  wFh = wFh + wtmp(:,:,:,tid+1)
  !$OMP END CRITICAL
  !$OMP END PARALLEL
  deallocate(utmp,vtmp,wtmp)
 
  ! Send and recieve halo velocity information 
  call haloExchangeMPI(uFh)
  call haloExchangeMPI(vFh)
  call haloExchangeMPI(wFh)

  ! TODO: Enforce no-penetration boundary condition
  !call enforceNoPenBCxyPeriodic()

  ! Copy data to final velocity arrays
  call getOutputFields()
  
  ! Interpolate velocity to dzF/2:dzF:Lz-dzF/2 mesh
  call interp%interpz_E2C(uG,uGout,0,1)
  call interp%interpz_E2C(vG,vGout,0,1)
  call interp%interpz_E2C(wG,wGout,0,-1)
  call message("Finished interpolating Gabor-induced velocity "//&
    "to cell centers")
  
end subroutine

subroutine enforceNoPenBCxyPeriodic()
  use domainSetup, only: zFE, Lx, Ly, Lz
  integer :: nx, ny, nz
  nx = size(uG,1)
  ny = size(uG,2)
  nz = size(uG,3)-1

  call assert(size(zFE) == nz+1,'size(zFE) == nz+1 -- GaborModeRoutines.F90')


end subroutine

subroutine getOutputFields()
  if (kstFh == 1 .and. kenFh == nzF+1) then
    uG = uFh(istFh+nxsupp/2:ienFh-nxsupp/2,jstFh+nysupp/2:jenFh-nysupp/2,:)
    vG = vFh(istFh+nxsupp/2:ienFh-nxsupp/2,jstFh+nysupp/2:jenFh-nysupp/2,:)
    wG = wFh(istFh+nxsupp/2:ienFh-nxsupp/2,jstFh+nysupp/2:jenFh-nysupp/2,:)
  elseif (kstFh == 1) then
    uG = uFh(istFh+nxsupp/2:ienFh-nxsupp/2,jstFh+nysupp/2:jenFh-nysupp/2,&
      1:kenFh-nzsupp/2)
    vG = vFh(istFh+nxsupp/2:ienFh-nxsupp/2,jstFh+nysupp/2:jenFh-nysupp/2,&
      1:kenFh-nzsupp/2)
    wG = wFh(istFh+nxsupp/2:ienFh-nxsupp/2,jstFh+nysupp/2:jenFh-nysupp/2,&
      1:kenFh-nzsupp/2)
  elseif (kenFh == nzF+1) then
    uG = uFh(istFh+nxsupp/2:ienFh-nxsupp/2,jstFh+nysupp/2:jenFh-nysupp/2,&
      kstFh+nzsupp/2:kenFh)
    vG = vFh(istFh+nxsupp/2:ienFh-nxsupp/2,jstFh+nysupp/2:jenFh-nysupp/2,&
      kstFh+nzsupp/2:kenFh)
    wG = wFh(istFh+nxsupp/2:ienFh-nxsupp/2,jstFh+nysupp/2:jenFh-nysupp/2,&
      kstFh+nzsupp/2:kenFh)
  else
    uG = uFh(istFh+nxsupp/2:ienFh-nxsupp/2,jstFh+nysupp/2:jenFh-nysupp/2,&
      kstFh+nzsupp/2:kenFh-nzsupp/2)
    vG = vFh(istFh+nxsupp/2:ienFh-nxsupp/2,jstFh+nysupp/2:jenFh-nysupp/2,&
      kstFh+nzsupp/2:kenFh-nzsupp/2)
    wG = wFh(istFh+nxsupp/2:ienFh-nxsupp/2,jstFh+nysupp/2:jenFh-nysupp/2,&
      kstFh+nzsupp/2:kenFh-nzsupp/2)
  end if
end subroutine
