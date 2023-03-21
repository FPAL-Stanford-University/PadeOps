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
    integer, dimension(2) :: shift
    real(rkind) :: xtgt, ytgt, ztgt, distance
    integer :: nmodesPrintStatus = 500000
    
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

    shift = [1,1]

    call tic()
    if (this%genModesOnUniformGrid) then
      do n = 1,size(x) 
        ! NOTE: These are global indices of the physical domain
        ! NOTE: The contribution of Gabor modes on neighboring processes is not
        ! accounted for here, nor is the periodic contribution for periodic
        ! directions whose data resides exlusively on the process (e.g. in x)
        call getStartEndGlobalIndices(x(n), xDom(1), wxSupport, this%smallScales%dx, &
          shift, ist, ien, iist, iien)
        call getStartEndGlobalIndices(y(n), yDom(1), wySupport, this%smallScales%dy, &
          shift, jst, jen, jjst, jjen)
        call getStartEndGlobalIndices(z(n), zDom(1), wzSupport, this%smallScales%dz, &
          shift, kst, ken, kkst, kken)
        
        if (iien < iist .or. jjen < jjst .or. kken < kkst) then
          continue
        else

          do k = kkst,kken
            ztgt = (k-1)*dz
            do j = jjst,jjen
              ytgt = (j-1)*dy
              do i = iist,iien
                xtgt = (i-1)*dx
                distance = sqrt( (x(n) - xtgt)**2 + (y(n) - ytgt)**2 + (z(n) - ztgt)**2)
                du = ceiling(wxSupport/2 - distance) ! 1 for modes inside and 0 for modes outside
                this%utmp(i,j,k) = this%utmp(i,j,k) + du 
              end do
            end do
          end do
        end if
        if (mod(n,nmodesPrintStatus) == 0) then
          write(mssg,'(F7.4,A10)')real(n,rkind)/real(size(x),rkind)*100.d0,'% Complete'
          call message(trim(mssg))
        end if
      end do
    else if (this%useFastTrigFunctions) then
      do n = 1,size(x) 
        ! NOTE: These are global indices of the physical domain
        ! NOTE: The contribution of Gabor modes on neighboring processes is not
        ! accounted for here, nor is the periodic contribution for periodic
        ! directions whose data resides exlusively on the process (e.g. in x)
        
        call getStartEndGlobalIndices(x(n), xDom(1), wxSupport, this%smallScales%dx, &
          shift, ist, ien, iist, iien)
        call getStartEndGlobalIndices(y(n), yDom(1), wySupport, this%smallScales%dy, &
          shift, jst, jen, jjst, jjen)
        call getStartEndGlobalIndices(z(n), zDom(1), wzSupport, this%smallScales%dz, &
          shift, kst, ken, kkst, kken)
        
        if (iien < iist .or. jjen < jjst .or. kken < kkst) then
          continue
        else
          ! Cast variables to single precision
          kxs = castSingle(kx(n)); kys = castSingle(ky(n)); kzs = castSingle(kz(n))
          xs  = castSingle(x(n));  ys  = castSingle(y(n)); zs  = castSingle(z(n))

          uRs = castSingle(uR(n)); vRs = castSingle(vR(n)); wRs = castSingle(wR(n))
          uIs = castSingle(uI(n)); vIs = castSingle(vI(n)); wIs = castSingle(wI(n))

          do k = kkst,kken
            zF = real(zDom(1),kind=4) + dz*castSingle(k - 1)
            kdotx3 = kzs*(zF - zs)
            fz = max(fastcos(pi_single*(zF - zs)/wzSupport), 0.e0)
            do j = jjst,jjen
              yF = real(yDom(1),kind=4) + dy*castSingle(j - 1)
              kdotx2 = kys*(yF - ys)
              fy = max(fastcos(pi_single*(yF - ys)/wySupport), 0.e0)
              do i = iist,iien
                xF = real(xDom(1),kind=4) + dx*castSingle(i - 1)
                kdotx = kdotx2 + kdotx3 + kxs*(xF - xs)

                cs = fastcos(kdotx)
                ss = fastsin(kdotx)

                fx = max(fastcos(pi_single*(xF - xs)/wxSupport), 0.e0)
                f = fx*fy*fz

                du = 2*f*(uRs*cs - uIs*ss)
                dv = 2*f*(vRs*cs - vIs*ss)
                dw = 2*f*(wRs*cs - wIs*ss)

                this%utmp(i,j,k) = this%utmp(i,j,k) + du 
                this%vtmp(i,j,k) = this%vtmp(i,j,k) + dv 
                this%wtmp(i,j,k) = this%wtmp(i,j,k) + dw
                
              end do
            end do
          end do
        end if
        if (mod(n,nmodesPrintStatus) == 0) then
          write(mssg,'(F7.4,A10)')real(n,rkind)/real(size(x),rkind)*100.d0,'% Complete'
          call message(trim(mssg))
        end if
      end do
    else
      do n = 1,size(x) 
        ! NOTE: These are global indices of the physical domain
        ! NOTE: The contribution of Gabor modes on neighboring processes is not
        ! accounted for here, nor is the periodic contribution for periodic
        ! directions whose data resides exlusively on the process (e.g. in x)
        
        call getStartEndGlobalIndices(x(n), xDom(1), wxSupport, this%smallScales%dx, &
          shift, ist, ien, iist, iien)
        call getStartEndGlobalIndices(y(n), yDom(1), wySupport, this%smallScales%dy, &
          shift, jst, jen, jjst, jjen)
        call getStartEndGlobalIndices(z(n), zDom(1), wzSupport, this%smallScales%dz, &
          shift, kst, ken, kkst, kken)
        
        if (iien < iist .or. jjen < jjst .or. kken < kkst) then
          continue
        else
          ! Cast variables to single precision
          kxs = castSingle(kx(n)); kys = castSingle(ky(n)); kzs = castSingle(kz(n))
          xs  = castSingle(x(n));  ys  = castSingle(y(n)); zs  = castSingle(z(n))

          uRs = castSingle(uR(n)); vRs = castSingle(vR(n)); wRs = castSingle(wR(n))
          uIs = castSingle(uI(n)); vIs = castSingle(vI(n)); wIs = castSingle(wI(n))

          do k = kkst,kken
            zF = real(zDom(1),kind=4) + dz*castSingle(k - 1)
            kdotx3 = kzs*(zF - zs)
            fz = max(cos(pi_single*(zF - zs)/wzSupport), 0.e0)
            do j = jjst,jjen
              yF = real(yDom(1),kind=4) + dy*castSingle(j - 1)
              kdotx2 = kys*(yF - ys)
              fy = max(cos(pi_single*(yF - ys)/wySupport), 0.e0)
              do i = iist,iien
                xF = real(xDom(1),kind=4) + dx*castSingle(i - 1)
                kdotx = kdotx2 + kdotx3 + kxs*(xF - xs)

                cs = cos(kdotx)
                ss = sin(kdotx)

                fx = max(cos(pi_single*(xF - xs)/wxSupport), 0.e0)
                f = fx*fy*fz

                du = 2*f*(uRs*cs - uIs*ss)
                dv = 2*f*(vRs*cs - vIs*ss)
                dw = 2*f*(wRs*cs - wIs*ss)

                this%utmp(i,j,k) = this%utmp(i,j,k) + du 
                this%vtmp(i,j,k) = this%vtmp(i,j,k) + dv 
                this%wtmp(i,j,k) = this%wtmp(i,j,k) + dw
                
              end do
            end do
          end do
        end if
        if (mod(n,nmodesPrintStatus) == 0) then
          write(mssg,'(F7.4,A10)')real(n,rkind)/real(size(x),rkind)*100.d0,'% Complete'
          call message(trim(mssg))
        end if
      end do
    end if
    this%smallScales%u  = this%smallScales%u  + real(this%utmp,rkind)
    this%smallScales%v  = this%smallScales%v  + real(this%vtmp,rkind)
    this%smallScales%wC = this%smallScales%wC + real(this%wtmp,rkind)
    call toc('Finished rendering velocity')
end subroutine
  
subroutine renderVelocity(this)
  class(enrichmentOperator), intent(inout), target :: this
  real(rkind) :: Lx, Ly, Lz
  real(rkind), dimension(:,:), pointer :: haloBuffY, haloBuffZ 
  
  call message(1,"Rendering the Gabor-induced velocity field")
 
  haloBuffY => null()
  haloBuffZ => null()

  Lx = xDom(2) - xDom(1)
  Ly = yDom(2) - yDom(1)
  Lz = zDom(2) - zDom(1)
  
  ! Zero the velocity arrays
  this%smallScales%u  = 0.d0
  this%smallScales%v  = 0.d0
  this%smallScales%wC = 0.d0
   
  ! STEP 1: Exchange Gabor modes from neighbors that have influence on your domain 
  if (allocated(this%renderModeData)) deallocate(this%renderModedata)
  allocate(this%renderModeData(this%nmodes, size(this%ModeData,2)))
  this%renderModeData = this%modeData
  call this%sendRecvHaloModes(this%renderModeData)
  
  ! Step 2: Render velocity 
    call this%renderLocalVelocity(this%renderModeData(:,1), &
      this%renderModeData(:,2),  this%renderModeData(:,3), &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))

  ! Step 3: Add x-periodic contribution
  include "enrichment_files/periodicContribution.F90"

  ! TODO:
  ! Step 2.b: interpolate wC to w. Do we need to get uhat, vhat, what from u,
  ! v, w? 

  ! Aditya
  ! STEP 3: Impose boundary condition (no-penetration BC)
  if (this%imposeNoPenetrationBC) then
    call this%smallScales%projectToFixBC()
  end if

  ! Aditya 
  ! STEP 4: Compute pressure (in needed)
  if (this%renderPressure) then 
    call this%smallScales%computePressure()
  end if 

  nullify(haloBuffY,haloBuffZ)
end subroutine  

subroutine getstartEndGlobalIndices(loc, domst, wSupport, delta, shift, st, en, stg, eng)
  real(rkind), intent(in) :: loc, delta, domst
  real(single_kind), intent(in) :: wSupport
  integer, dimension(2), intent(in) :: shift
  integer, intent(in) :: st, en
  integer, intent(out) :: stg, eng

  stg = max(ceiling((loc - domst - wSupport/2 - delta/2)/delta) + shift(1), st)
  eng = min(ceiling((loc - domst + wSupport/2 - delta/2)/delta) + shift(2), en)
end subroutine

recursive function fastsin(xin) result (y)
  real(single_kind), intent(in) :: xin
  real(single_kind) :: y, x

  x = xin

  ! Shift all negative x-values to the positive x-axis
  x = abs(x + 0.5e0*pi_single) - 0.5e0*pi_single;

  ! Wrap x-values that are greater than pi
  if (x > 2.e0*pi_single) then
      y = fastsin(x - 2.e0*pi_single);
  else if (x > pi_single) then
      y = -fastsin(x - pi_single);
  else
      y = 1.2732e0*x - 0.4053e0*x*abs(x);
      y = 0.225e0*(y*abs(y) - y) + y;
  end if
end function

function fastcos(xin) result(y)
  real(single_kind), intent(in) :: xin
  real(single_kind) :: y, x
  x = xin
  y = fastsin(x + 0.5e0*pi_single)
end function
