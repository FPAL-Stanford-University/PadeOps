subroutine renderLocalVelocity(this,x,y,z,kx,ky,kz,uR,uI,vR,vI,wR,wI)
    use exits,     only: message
    use omp_lib,   only: omp_get_thread_num, omp_get_num_threads
    use decomp_2d, only: nrank
    ! Render the velocity local to the MPI rank
    ! INPUTS:
    !   x, y, z     --> Gabor mode location
    !   kx, ky,kz   --> Gabor mode wave-vector components 
    !   uR, uI, ... --> Gabor mode velocity amplitudes
    
    class(enrichmentOperator), intent(inout) :: this
    real(rkind), dimension(:), intent(in) :: uR, uI, vR, vI, wR, wI, x, y, z, &
      kx, ky, kz
    integer :: n
    real(rkind) :: small = 1.d-14
    character(len=clen) :: mssg
    real(single_kind) :: kdotx, kdotx2, kdotx3
    integer :: i, j, k, iist, iien, jjst, jjen, kkst, kken
    integer :: ist, ien, jst, jen, kst, ken
    real(single_kind) :: cs, ss, fx, fy, fz, f, xF, yF, zF, kxs, kys, kzs
    real(single_kind) :: uRs, uIs, vRs, vIs, wRs, wIs, xs, ys, zs, dx, dy, dz
    real(single_kind) :: wxSupport, wySupport, wzSupport, du, dv, dw
    integer :: tid
    real(single_kind), parameter :: pi_single = 4.0*atan(1.0)

    call message(1,"Rendering the Gabor-induced velocity field")
    
    this%utmp = 0.e0
    this%vtmp = 0.e0
    this%wtmp = 0.e0

    ! Cast variables to single precision
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

    !$OMP PARALLEL &
    !$OMP PRIVATE(tid,n,i,j,k,kdotx,kdotx3,kdotx2,fx,fy,fz,f) &
    !$OMP PRIVATE(cs,ss,iist,iien,jjst,jjen,kkst,kken)
    tid = omp_get_thread_num()
    if (tid == 0 .and. nrank == 0) print*, "# of threads spawned:", omp_get_num_threads()
    !$OMP DO
    do n = 1,size(x) 
      ! NOTE: These are global indices of the physical domain
      ! NOTE: The contribution of Gabor modes on neighboring processes is not
      ! accounted for here, nor is the periodic contribution for periodic
      ! directions whose data resides exlusively on the process (e.g. in x)
      iist = max(ceiling((x(n)+small)/this%smallScales%dx) - this%nxsupp/2, ist)
      !iien = min(ceiling((x(n)+small)/this%smallScales%dx) + this%nxsupp/2, ien)
      iien = min(floor(  (x(n)+small)/this%smallScales%dx) + this%nxsupp/2, ien)

      jjst = max(ceiling((y(n)+small)/this%smallScales%dy) - this%nysupp/2, jst)
      !jjen = min(ceiling((y(n)+small)/this%smallScales%dy) + this%nysupp/2, jen)
      jjen = min(floor(  (y(n)+small)/this%smallScales%dy) + this%nysupp/2, jen)

      kkst = max(ceiling((z(n)+small)/this%smallScales%dz) - this%nzsupp/2, kst)
      !kken = min(ceiling((z(n)+small)/this%smallScales%dz) + this%nzsupp/2, ken)
      kken = min(floor(  (z(n)+small)/this%smallScales%dz) + this%nzsupp/2, ken)

      if (iien < iist .or. jjen < jjst .or. kken < kkst) then
        continue
      else
        ! Cast variables to single precision
        kxs = castSingle(kx(n))
        kys = castSingle(ky(n))
        kzs = castSingle(kz(n))
        
        xs  = castSingle(x(n))
        ys  = castSingle(y(n))
        zs  = castSingle(z(n))

        uRs = castSingle(uR(n))
        uIs = castSingle(uI(n))
        vRs = castSingle(vR(n))
        vIs = castSingle(vI(n))
        wRs = castSingle(wR(n))
        wIs = castSingle(wI(n))

        do k = kkst,kken
          zF = dz*castSingle(k - 1)
          kdotx3 = kzs*(zF - zs)
          fz = cos(pi_single*(zF - zs)/wzSupport)
          do j = jjst,jjen
            yF = dy*castSingle(j - 1)
            kdotx2 = kys*(yF - ys)
            fy = cos(pi_single*(yF - ys)/wySupport)
            do i = iist,iien
              xF = dx*castSingle(i - 1)
              kdotx = kdotx2 + kdotx3 + kxs*(xF - xs)

              cs = cos(kdotx)
              ss = sin(kdotx)

              fx = cos(pi_single*(xF - xs)/wxSupport)
              f = fx*fy*fz

              du = 2*f*(uRs*cs - uIs*ss)
              dv = 2*f*(vRs*cs - vIs*ss)
              dw = 2*f*(wRs*cs - wIs*ss)

              this%utmp(i,j,k,tid+1) = this%utmp(i,j,k,tid+1) + du 
              this%vtmp(i,j,k,tid+1) = this%vtmp(i,j,k,tid+1) + dv 
              this%wtmp(i,j,k,tid+1) = this%wtmp(i,j,k,tid+1) + dw 
            end do
          end do
        end do
      end if
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
