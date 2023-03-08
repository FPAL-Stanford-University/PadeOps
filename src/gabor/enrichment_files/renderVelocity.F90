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
    real(single_kind), parameter :: pi_single = 4.0*atan(1.0)
    integer, dimension(2) :: shift
    
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

    if (this%genModesOnUniformGrid) then
      shift = [2,0]
    else
      shift = [1,1]
    end if

    do n = 1,size(x) 
      ! NOTE: These are global indices of the physical domain
      ! NOTE: The contribution of Gabor modes on neighboring processes is not
      ! accounted for here, nor is the periodic contribution for periodic
      ! directions whose data resides exlusively on the process (e.g. in x)
      !iist = max(ceiling((x(n) - xDom(1))/this%smallScales%dx) - this%nxsupp/2 - 1, ist)
      !iien = min(floor(  (x(n) - xDom(1))/this%smallScales%dx) + this%nxsupp/2 + 1, ien)

      !jjst = max(ceiling((y(n) - yDom(1))/this%smallScales%dy) - this%nysupp/2 - 1, jst)
      !jjen = min(floor(  (y(n) - yDom(1))/this%smallScales%dy) + this%nysupp/2 + 1, jen)

      !kkst = max(ceiling((z(n) - zDom(1))/this%smallScales%dz) - this%nzsupp/2 - 1, kst)
      !kken = min(floor(  (z(n) - zDom(1))/this%smallScales%dz) + this%nzsupp/2 + 1, ken)
      
      iist = max(nint((x(n) - xDom(1) - wxSupport/2 - this%smallScales%dx/2)/this%smallScales%dx) + shift(1), ist)
      iien = min(nint((x(n) - xDom(1) + wxSupport/2 - this%smallScales%dx/2)/this%smallScales%dx) + shift(2), ien)

      jjst = max(nint((y(n) - yDom(1) - wySupport/2 - this%smallScales%dy/2)/this%smallScales%dy) + shift(1), jst)
      jjen = min(nint((y(n) - yDom(1) + wySupport/2 - this%smallScales%dy/2)/this%smallScales%dy) + shift(2), jen)

      kkst = max(nint((z(n) - zDom(1) - wzSupport/2 - this%smallScales%dz/2)/this%smallScales%dz) + shift(1), kst)
      kken = min(nint((z(n) - zDom(1) + wzSupport/2 - this%smallScales%dz/2)/this%smallScales%dz) + shift(2), ken)
      
      !iist = max(floor(  x(n)/this%smallScales%dx) - this%nxsupp/2, ist)
      !iien = min(ceiling(x(n)/this%smallScales%dx) + this%nxsupp/2, ien)

      !jjst = max(floor(  y(n)/this%smallScales%dy) - this%nysupp/2, jst)
      !jjen = min(ceiling(y(n)/this%smallScales%dy) + this%nysupp/2, jen)
      !
      !kkst = max(floor(  z(n)/this%smallScales%dz) - this%nzsupp/2, kst)
      !kken = min(ceiling(z(n)/this%smallScales%dz) + this%nzsupp/2, ken)
      
      if (iien < iist .or. jjen < jjst .or. kken < kkst) then
        continue
      else
        ! Cast variables to single precision
        kxs = 0.e0!castSingle(kx(n))
        kys = 0.e0!castSingle(ky(n))
        kzs = 0.e0!castSingle(kz(n))
        
        xs  = castSingle(x(n))
        ys  = castSingle(y(n))
        zs  = castSingle(z(n))

        uRs = 0.5e0!castSingle(uR(n))
        uIs = 0.e0!castSingle(uI(n))
        vRs = 0.e0!castSingle(vR(n))
        vIs = 0.e0!castSingle(vI(n))
        wRs = 0.e0!castSingle(wR(n))
        wIs = 0.e0!castSingle(wI(n))


        do k = kkst,kken
          zF = real(zDom(1),kind=4) + dz*castSingle(k - 1)
          kdotx3 = kzs*(zF - zs)
          !fz = max(cos(pi_single*(zF - zs)/wzSupport), 0.e0)
          fz = 1.e0
          do j = jjst,jjen
            yF = real(yDom(1),kind=4) + dy*castSingle(j - 1)
            kdotx2 = kys*(yF - ys)
            !fy = max(cos(pi_single*(yF - ys)/wySupport), 0.e0)
            fy = 1.e0
            do i = iist,iien
              xF = real(xDom(1),kind=4) + dx*castSingle(i - 1)
              kdotx = kdotx2 + kdotx3 + kxs*(xF - xs)

              cs = cos(kdotx)
              ss = sin(kdotx)

              !fx = max(cos(pi_single*(xF - xs)/wxSupport), 0.e0)
              fx = 1.e0
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
      if (mod(n,100000) == 0) then
        write(mssg,'(F7.4,A10)')real(n,rkind)/real(size(x),rkind)*100.d0,'% Complete'
        call message(trim(mssg))
      end if
    end do

    this%smallScales%u  = this%smallScales%u  + real(this%utmp,rkind)
    this%smallScales%v  = this%smallScales%v  + real(this%vtmp,rkind)
    this%smallScales%wC = this%smallScales%wC + real(this%wtmp,rkind)
end subroutine
  
subroutine renderVelocity(this)
  class(enrichmentOperator), intent(inout), target :: this
  real(rkind) :: Lx
  real(rkind), dimension(:,:), pointer :: haloBuffY, haloBuffZ 
  
  call message(1,"Rendering the Gabor-induced velocity field")
 
  haloBuffY => null()
  haloBuffZ => null()

  Lx = xDom(2) - xDom(1)
  
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
  if (periodicBCs(1)) then
    call this%renderLocalVelocity(this%renderModeData(:,1) + Lx, &
      this%renderModeData(:,2),  this%renderModeData(:,3), &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
    call this%renderLocalVelocity(this%renderModeData(:,1) - Lx, &
      this%renderModeData(:,2),  this%renderModeData(:,3), &
      this%renderModeData(:,4),  this%renderModeData(:,5), &
      this%renderModeData(:,6),  this%renderModeData(:,7), &
      this%renderModeData(:,8),  this%renderModeData(:,9), &
      this%renderModeData(:,10), this%renderModeData(:,11), &
      this%renderModeData(:,12))
  end if

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

