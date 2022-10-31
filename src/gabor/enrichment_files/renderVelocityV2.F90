subroutine renderLocalVelocity(this)
    use exits,     only: message
    use omp_lib,   only: omp_get_thread_num, omp_get_num_threads
    use decomp_2d, only: nrank
    
    class(enrichmentOperator), intent(inout), target :: this
    integer :: n
    real(rkind) :: small = 1.d-14
    character(len=clen) :: mssg
    integer :: i, j, k, iist, iien, jjst, jjen, kkst, kken
    integer :: ist, ien, jst, jen, kst, ken
    real(single_kind), dimension(:,:,:), pointer :: cs, ss, f, xF, yF, zF
    real(single_kind), dimension(:,:,:), pointer :: kdotx
    real(single_kind) :: dx, dy, dz, x, y, z
    real(single_kind) :: uR, uI, vR, vI, wR, wI
    real(single_kind) :: kx, ky, kz
    real(single_kind) :: wxSupport, wySupport, wzSupport
    integer :: tid, idx

    call message("Rendering the Gabor-induced velocity field")
  
    ! Make single precision versions of various parameters
    dx = castSingle(this%smallScales%dx)
    dy = castSingle(this%smallScales%dy)
    dz = castSingle(this%smallScales%dz)

    wxSupport = castSingle(this%nxsupp+1)*dx
    wySupport = castSingle(this%nysupp+1)*dy
    wzSupport = castSingle(this%nzsupp+1)*dz
     
    ! Allocate velocity arrays for each thread 
    ist = this%smallScales%gpC%xst(1) 
    ien = this%smallScales%gpC%xen(1) 
    jst = this%smallScales%gpC%xst(2) 
    jen = this%smallScales%gpC%xen(2) 
    kst = this%smallScales%gpC%xst(3) 
    ken = this%smallScales%gpC%xen(3)

    ! Zero the arrays
    this%utmp = 0.e0
    this%vtmp = 0.e0
    this%wtmp = 0.e0

    this%smallScales%u  = 0.d0
    this%smallScales%v  = 0.d0
    this%smallScales%wC = 0.d0

    !$OMP PARALLEL &
    !$OMP PRIVATE(tid,n,i,j,k,kdotx,f,xF,yF,zF) &
    !$OMP PRIVATE(cs,ss,iist,iien,jjst,jjen,kkst,kken)
    tid = omp_get_thread_num()
    if (tid == 0 .and. nrank == 0) print*, "# of threads spawned:", omp_get_num_threads()
    !$OMP DO
    do n = 1,this%nmodes  
      ! NOTE: These are global indices of the physical domain
      ! NOTE: The contribution of Gabor modes on neighboring processes is not
      ! accounted for here, nor is the periodic contribution for periodic
      ! directions whose data resides exlusively on the process (e.g. in x)
      iist = max(ceiling((this%x(n)+small)/this%smallScales%dx) - this%nxsupp/2, ist)
      iien = min(floor(  (this%x(n)+small)/this%smallScales%dx) + this%nxsupp/2, ien)

      jjst = max(ceiling((this%y(n)+small)/this%smallScales%dy) - this%nysupp/2, jst)
      jjen = min(floor(  (this%y(n)+small)/this%smallScales%dy) + this%nysupp/2, jen)

      kkst = max(ceiling((this%z(n)+small)/this%smallScales%dz) - this%nzsupp/2, kst)
      kken = min(floor(  (this%z(n)+small)/this%smallScales%dz) + this%nzsupp/2, ken)

      xF    => this%buff(iist:iien,jjst:jjen,kkst:kken,tid+1,1)
      yF    => this%buff(iist:iien,jjst:jjen,kkst:kken,tid+1,2)
      zF    => this%buff(iist:iien,jjst:jjen,kkst:kken,tid+1,3)

      kdotx => this%buff(iist:iien,jjst:jjen,kkst:kken,tid+1,4)
      cs    => this%buff(iist:iien,jjst:jjen,kkst:kken,tid+1,5)
      ss    => this%buff(iist:iien,jjst:jjen,kkst:kken,tid+1,6)
      f     => this%buff(iist:iien,jjst:jjen,kkst:kken,tid+1,7)

      ! Creat single-precision version of variables
      uR = castSingle(this%uhatR(n))
      uI = castSingle(this%uhatI(n))
      vR = castSingle(this%vhatR(n))
      vI = castSingle(this%vhatI(n))
      wR = castSingle(this%whatR(n))
      wI = castSingle(this%whatI(n))
      
      kx = castSingle(this%kx(n))
      ky = castSingle(this%ky(n))
      kz = castSingle(this%kz(n))

      x  = castSingle(this%x(n))
      y  = castSingle(this%y(n))
      z  = castSingle(this%z(n))

      idx = 1
      do k = kkst,kken
        zF(:,:,idx) = dz*castSingle(k - 1)
        idx = idx + 1
      end do
      idx = 1
      do j = jjst,jjen
        yF(:,idx,:) = dy*castSingle(j - 1)
        idx = idx + 1
      end do
      idx = 1
      do i = iist,iien
        xF(idx,:,:) = dx*castSingle(i - 1)
        idx = idx + 1
      end do
      
      kdotx = kx*(xF - x) + ky*(yF - y) + kz*(zF - z)
      cs = cos(kdotx)
      ss = sin(kdotx)

      f =     cos(pi*(xF - x)/wxSupport)
      f = f * cos(pi*(yF - y)/wySupport)
      f = f * cos(pi*(zF - z)/wzSupport)

      this%utmp(iist:iien,jjst:jjen,kkst:kken,tid+1) = &
        this%utmp(iist:iien,jjst:jjen,kkst:kken,tid+1) + &
        f*(2.d0*uR*cs - 2.d0*uI*ss) 
      
      this%vtmp(iist:iien,jjst:jjen,kkst:kken,tid+1) = &
        this%vtmp(iist:iien,jjst:jjen,kkst:kken,tid+1) + &
        f*(2.d0*vR*cs - 2.d0*vI*ss)
      
      this%wtmp(iist:iien,jjst:jjen,kkst:kken,tid+1) = &
        this%wtmp(iist:iien,jjst:jjen,kkst:kken,tid+1) + &
        f*(2.d0*wR*cs - 2.d0*wI*ss)

      nullify(xF, yF, zF, kdotx, cs, ss, f)
      if (mod(n,100000) == 0 .and. tid == 0) then
        write(mssg,'(F7.4,A10)')real(n,rkind)/real(this%nmodes,rkind)*100.d0,'% Complete'
        call message(trim(mssg))
      end if
    end do
    !$OMP END DO
    !$OMP CRITICAL
    this%smallScales%u  = this%smallScales%u  + real(this%utmp(:,:,:,tid+1),rkind)
    this%smallScales%v  = this%smallScales%v  + real(this%vtmp(:,:,:,tid+1),rkind)
    this%smallScales%wC = this%smallScales%wC + real(this%wtmp(:,:,:,tid+1),rkind)
    !$OMP END CRITICAL
    !$OMP END PARALLEL
end subroutine
