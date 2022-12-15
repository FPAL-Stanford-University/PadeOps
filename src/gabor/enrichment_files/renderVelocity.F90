subroutine renderLocalVelocity(this,x,y,z,kx,ky,kz,uR,uI,vR,vI,wR,wI)
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
    character(len=clen) :: mssg
    real(single_kind) :: kdotx, kdotx2, kdotx3
    integer :: i, j, k, iist, iien, jjst, jjen, kkst, kken
    integer :: ist, ien, jst, jen, kst, ken
    real(single_kind) :: cs, ss, fx, fy, fz, f, xF, yF, zF, kxs, kys, kzs
    real(single_kind) :: uRs, uIs, vRs, vIs, wRs, wIs, xs, ys, zs, dx, dy, dz
    real(single_kind) :: wxSupport, wySupport, wzSupport, du, dv, dw
    integer :: tid
    real(single_kind), parameter :: pi_single = 4.0*atan(1.0)
    
    this%utmp = 0.e0
    this%vtmp = 0.e0
    this%wtmp = 0.e0

    ! Cast variables to single precision
    dx = castSingle(this%smallScales%dx) 
    dy = castSingle(this%smallScales%dy) 
    dz = castSingle(this%smallScales%dz) 

    wxSupport = castSingle(this%nxsupp)*dx
    wySupport = castSingle(this%nysupp)*dy
    wzSupport = castSingle(this%nzsupp)*dz
     
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
    !$OMP DO
    do n = 1,size(x) 
      ! NOTE: These are global indices of the physical domain
      ! NOTE: The contribution of Gabor modes on neighboring processes is not
      ! accounted for here, nor is the periodic contribution for periodic
      ! directions whose data resides exlusively on the process (e.g. in x)
      iist = max(ceiling((x(n) - xDom(1))/this%smallScales%dx) - this%nxsupp/2 - 1, ist)
      iien = min(floor(  (x(n) - xDom(1))/this%smallScales%dx) + this%nxsupp/2 + 1, ien)

      jjst = max(ceiling((y(n) - yDom(1))/this%smallScales%dy) - this%nysupp/2 - 1, jst)
      jjen = min(floor(  (y(n) - yDom(1))/this%smallScales%dy) + this%nysupp/2 + 1, jen)

      kkst = max(ceiling((z(n) - zDom(1))/this%smallScales%dz) - this%nzsupp/2 - 1, kst)
      kken = min(floor(  (z(n) - zDom(1))/this%smallScales%dz) + this%nzsupp/2 + 1, ken)
      
      !iist = max(floor(  x(n)/this%smallScales%dx) - this%nxsupp/2, ist)
      !iien = min(ceiling(x(n)/this%smallScales%dx) + this%nxsupp/2, ien)

      !jjst = max(floor(  y(n)/this%smallScales%dy) - this%nysupp/2, jst)
      !jjen = min(ceiling(y(n)/this%smallScales%dy) + this%nysupp/2, jen)
      !
      !kkst = max(floor(  z(n)/this%smallScales%dz) - this%nzsupp/2, kst)
      !kken = min(ceiling(z(n)/this%smallScales%dz) + this%nzsupp/2, ken)
!if (n == 1) print*, jjst, jjen
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
          zF = zDom(1) + dz*castSingle(k - 1)
          kdotx3 = kzs*(zF - zs)
          fz = max(cos(pi_single*(zF - zs)/wzSupport), 0.e0)
          do j = jjst,jjen
            yF = yDom(1) + dy*castSingle(j - 1)
            kdotx2 = kys*(yF - ys)
            fy = max(cos(pi_single*(yF - ys)/wySupport), 0.e0)
            do i = iist,iien
              xF = xDom(1) + dx*castSingle(i - 1)
              kdotx = kdotx2 + kdotx3 + kxs*(xF - xs)

              cs = cos(kdotx)
              ss = sin(kdotx)

              fx = max(cos(pi_single*(xF - xs)/wxSupport), 0.e0)
              f = fx*fy*fz

              du = 2*f*(uRs*cs - uIs*ss)
              dv = 2*f*(vRs*cs - vIs*ss)
              dw = 2*f*(wRs*cs - wIs*ss)

              this%utmp(i,j,k,tid+1) = this%utmp(i,j,k,tid+1) + du 
              this%vtmp(i,j,k,tid+1) = this%vtmp(i,j,k,tid+1) + dv 
              this%wtmp(i,j,k,tid+1) = this%wtmp(i,j,k,tid+1) + dw
              
              !if (i == 11 .and. j == 13 .and. k == 16) then
              !  print*, this%utmp(i,j,k,tid+1), xs, ys, zs
              !end if
            end do
          end do
        end do
      end if
      if (mod(n,100000) == 0 .and. tid == 0) then
        write(mssg,'(F7.4,A10)')real(n,rkind)/real(size(x),rkind)*100.d0,'% Complete'
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
