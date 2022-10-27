subroutine renderLocalVelocity(this)
    use exits,     only: message
    use omp_lib,   only: omp_get_thread_num, omp_get_num_threads
    use decomp_2d, only: nrank
    
    class(enrichmentOperator), intent(inout) :: this
    integer :: n
    real(rkind) :: small = 1.d-14
    character(len=clen) :: mssg
    real(rkind) :: kdotx, kdotx2, kdotx3
    integer :: i, j, k, iist, iien, jjst, jjen, kkst, kken
    integer :: ist, ien, jst, jen, kst, ken
    real(rkind) :: cs, ss, fx, fy, fz, f, xF, yF, zF
    real(rkind) :: wxSupport, wySupport, wzSupport, du, dv, dw
    integer :: tid

    call message("Rendering the Gabor-induced velocity field")
  
    wxSupport = real(this%nxsupp+1,rkind)*this%smallScales%dx
    wySupport = real(this%nysupp+1,rkind)*this%smallScales%dy
    wzSupport = real(this%nzsupp+1,rkind)*this%smallScales%dz
     
    ! Allocate velocity arrays for each thread 
    ist = this%smallScales%gpC%xst(1) 
    ien = this%smallScales%gpC%xen(1) 
    jst = this%smallScales%gpC%xst(2) 
    jen = this%smallScales%gpC%xen(2) 
    kst = this%smallScales%gpC%xst(3) 
    ken = this%smallScales%gpC%xen(3)

    ! Zero the arrays
    this%utmp = 0.d0
    this%vtmp = 0.d0
    this%wtmp = 0.d0

    this%smallScales%u  = 0.d0
    this%smallScales%v  = 0.d0
    this%smallScales%wC = 0.d0

    !$OMP PARALLEL &
    !$OMP PRIVATE(tid,n,i,j,k,kdotx,kdotx3,kdotx2,fx,fy,fz,f) &
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

      do k = kkst,kken
        zF = this%smallScales%dz*real(kkst - 1, rkind)
        kdotx3 = this%kz(n)*(zF - this%z(n))
        fz = cos(pi*(zF - this%z(n))/wzSupport)
        do j = jjst,jjen
          yF = this%smallScales%dy*real(jjst - 1, rkind)
          kdotx2 = this%ky(n)*(yF - this%y(n))
          fy = cos(pi*(yF - this%y(n))/wySupport)
          do i = iist,iien
            xF = this%smallScales%dx*real(iist - 1, rkind)
            kdotx = kdotx2 + kdotx3 + this%kx(n)*(xF - this%x(n))

            cs = cos(kdotx)
            ss = sin(kdotx)

            fx = cos(pi*(xF - this%x(n))/wxSupport)
            f = fx*fy*fz

            du = f*(2.d0*this%uhatR(n)*cs - 2.d0*this%uhatI(n)*ss)
            dv = f*(2.d0*this%vhatR(n)*cs - 2.d0*this%vhatI(n)*ss)
            dw = f*(2.d0*this%whatR(n)*cs - 2.d0*this%whatI(n)*ss)

            this%utmp(i,j,k,tid+1) = this%utmp(i,j,k,tid+1) + du 
            this%vtmp(i,j,k,tid+1) = this%vtmp(i,j,k,tid+1) + dv 
            this%wtmp(i,j,k,tid+1) = this%wtmp(i,j,k,tid+1) + dw 
          end do
        end do
      end do

      if (mod(n,100000) == 0 .and. tid == 0) then
        write(mssg,'(F7.4,A10)')real(n,rkind)/real(this%nmodes,rkind)*100.d0,'% Complete'
        call message(trim(mssg))
      end if
    end do
    !$OMP END DO
    !$OMP CRITICAL
    this%smallScales%u  = this%smallScales%u  + this%utmp(:,:,:,tid+1)
    this%smallScales%v  = this%smallScales%v  + this%vtmp(:,:,:,tid+1)
    this%smallScales%wC = this%smallScales%wC + this%wtmp(:,:,:,tid+1)
    !$OMP END CRITICAL
    !$OMP END PARALLEL
end subroutine
