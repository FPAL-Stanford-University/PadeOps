subroutine destroy(this)
  class(sgs_igrid), intent(inout) :: this
  nullify(this%gpC, this%gpE, this%spectC, this%spectE, this%sp_gpC, this%sp_gpE, this%fxC)
  nullify(this%cbuffyC, this%cbuffzC, this%rbuffxC, this%Tsurf, this%fyC, this%fzE, this%wTh_surf)
  if (this%isEddyViscosityModel) call this%destroyMemory_EddyViscosity()
  if (this%DynamicProcedureType .ne. 0) call this%destroyMemory_DynamicProcedure()
  select case (this%mid)
  case (0)
     call this%destroy_smagorinsky()
  case (1)
     call this%destroy_sigma()
  end select
  nullify(this%tau_11, this%tau_12, this%tau_22, this%tau_33)
  deallocate(this%tau_13, this%tau_23)
  if(this%is_z0_varying) then
      deallocate(this%ustarsqvar, this%z0var, this%Uxvar, this%Uyvar, this%WallMFactorvar, this%lamfact, this%deli, this%mask_upstream, this%alpfac, this%Uxspan, this%Uyspan, this%ustarsqspan) !!!!
  endif

end subroutine


subroutine link_pointers(this, nuSGS, tauSGS_ij, tau13, tau23, q1, q2, q3, kappaSGS, kappa_bounding)
   class(sgs_igrid), intent(in), target :: this
   real(rkind), dimension(:,:,:)  , pointer, intent(inout) :: nuSGS
   real(rkind), dimension(:,:,:)  , pointer, intent(inout) :: tau13, tau23
   real(rkind), dimension(:,:,:,:), pointer, intent(inout) :: tauSGS_ij
   real(rkind), dimension(:,:,:)  , pointer, intent(inout) :: q1, q2, q3, kappaSGS
   real(rkind), dimension(:,:,:)  , pointer, optional, intent(inout)  :: kappa_bounding 

   nuSGS => this%nu_sgs_C
   tau13 => this%tau_13
   tau23 => this%tau_23

   tauSGS_ij => this%tau_ij
  
   if (this%isStratified) then
      q1 => this%q1C
      q2 => this%q2C
      q3 => this%q3E
      kappaSGS => this%kappa_sgs_C
   end if

   if (this%useScalarBounding) then
      if(present(kappa_bounding)) then
         kappa_bounding => this%kappa_boundingC
      end if 
   end if
end subroutine 

subroutine init(this, gpC, gpE, spectC, spectE, dx, dy, dz, inputfile, Lx, Ly, xMesh, zMeshE, zMeshC, fBody_x, fBody_y, fBody_z, computeFbody, PadeDer, cbuffyC, cbuffzC, cbuffyE, cbuffzE, rbuffxC, rbuffyC, rbuffzC, rbuffyE, rbuffzE, Tsurf, ThetaRef, wTh_surf, Fr, Re, isInviscid, isStratified, botBC_temp, useShiftedPeriodicBC, Lxeff, initSpinUp)
  class(sgs_igrid), intent(inout), target :: this
  class(decomp_info), intent(in), target :: gpC, gpE
  class(spectral), intent(in), target :: spectC, spectE
  real(rkind), intent(in) :: dx, dy, dz, ThetaRef, Fr, Re, Lx, Ly, Lxeff
  real(rkind), intent(in), target :: Tsurf, wTh_surf
  character(len=*), intent(in) :: inputfile
  real(rkind), dimension(:), intent(in) :: zMeshE, zMeshC, xMesh
  real(rkind), dimension(:,:,:), intent(in), target :: fBody_x, fBody_y, fBody_z
  logical, intent(out) :: computeFbody
  real(rkind), dimension(:,:,:,:), intent(in), target :: rbuffxC, rbuffyE, rbuffzE, rbuffyC, rbuffzC
  complex(rkind), dimension(:,:,:,:), intent(in), target :: cbuffyC, cbuffzC, cbuffyE, cbuffzE
  type(Pade6stagg), target, intent(in) :: PadeDer
  logical, intent(in) :: isInviscid, isStratified, useShiftedPeriodicBC
  integer, intent(in) :: botBC_temp
  logical, intent(in), optional :: initSpinUp

  ! Input file variables
  logical :: DomainAveraged_DynProc = .false., useWallDamping = .false., useSGSDynamicRestart = .false., useVerticalTfilter = .false.
  integer :: DynamicProcedureType = 0, SGSmodelID = 0, WallModelType = 0, DynProcFreq = 1
  real(rkind) :: ncWall = 1.d0, Csgs = 0.17d0, z0 = 0.01d0, deltaRatio = 2.d0, turbPrandtl = 0.4d0, Cy = 100.d0 
  real(rkind) :: z0t = 0.001d0
  character(len=clen) :: SGSDynamicRestartFile, fname, tempname
  logical :: explicitCalcEdgeEddyViscosity = .false., UseDynamicProcedureScalar = .false., useScalarBounding = .false.
  logical :: usePrSGS = .false. 
  integer :: ierr, WM_matchingIndex = 1, i, j
  real(rkind) :: lowbound = 0.d0 , highbound = 1.d0, xpos = 0.0_rkind, epssmall = 1.0d-6
  real(rkind), dimension(gpC%xsz(1)) :: sp_map, x1, x2, S1, S2, deli
  logical :: is_z0_varying = .false., filter_for_heterog = .true., use_const_alpfac = .true.
  real(rkind) :: z0r, z0s, spx, spy, rpx, rpy, totpx, xl, xlmod, rpstart = -1.0d0, spx_delta = 1.0d0, spy_delta = 1.0d0
  real(rkind) :: Mfactor, dele, excludedist = 3.0_rkind, const_alpfac = 0.027_rkind, blht_for_alpfac = 3.0_rkind, betfac = 0.15_rkind
  integer :: spnumx, spnumy

  namelist /SGS_MODEL/ DynamicProcedureType, SGSmodelID, z0, z0t, &
                 useWallDamping, ncWall, Csgs, WallModelType, usePrSGS, &
                 DynProcFreq, useSGSDynamicRestart, useVerticalTfilter,&
                 DomainAveraged_DynProc, SGSDynamicRestartFile, &
                 explicitCalcEdgeEddyViscosity, &
                 UseDynamicProcedureScalar, deltaRatio, turbPrandtl, &
                 useScalarBounding, Cy, lowbound, highbound, WM_matchingIndex, &
                 is_z0_varying

  namelist /Z0VARYING/ spx, spy, rpx, rpy, spnumx, spnumy, z0s, z0r, rpstart, spx_delta, spy_delta, &
                       filter_for_heterog, excludedist, use_const_alpfac, const_alpfac, blht_for_alpfac, betfac

  this%gpC => gpC
  this%gpE => gpE
  this%spectC => spectC
  this%spectE => spectE
  this%sp_gpC => spectC%spectdecomp
  this%sp_gpE => spectE%spectdecomp
  this%Tsurf => Tsurf
  this%wTh_surf => wTh_surf
  this%Fr = Fr
  this%Re = Re
  !this%Pr = Pr
  this%Pr = turbPrandtl
  this%ThetaRef = ThetaRef
  this%PadeDer => PadeDer
  this%fxC => fBody_x
  this%fyC => fBody_y
  this%fzE => fBody_z
  this%isStratified = isStratified
  this%usePrSGS = usePrSGS
  !if (present(botBC_Temp)) 
  this%botBC_Temp = botBC_Temp

  ! compute nxeff for horz averaging if useShiftedPeriodicBC is true
  this%useShiftedPeriodicBC = useShiftedPeriodicBC
  if(this%useShiftedPeriodicBC) then
    this%nxeff = minloc(abs(xMesh-Lxeff*Lx), 1)
    this%nxeff = this%nxeff-4   ! 4 is arbitrary here
    print *, 'Lxeff = ', Lxeff*Lx
    print *, 'nxeff = ', this%nxeff
  endif

  if(this%useShiftedPeriodicBC) then
      this%meanFact = one/(real(this%nxeff,rkind) * real(gpC%ysz(2),rkind))
  else 
      this%meanFact = one/(real(gpC%xsz(1),rkind) * real(gpC%ysz(2),rkind))
  endif
  print *, "meanFact = ", this%meanFact

  this%dx = dx
  this%dy = dy
  this%dz = dz
  this%DomainAveraged_DynProc = DomainAveraged_DynProc

  allocate(this%tau_ij(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3),6))
  this%tau_11   => this%tau_ij(:,:,:,1)
  this%tau_12   => this%tau_ij(:,:,:,2)
  this%tau_13C  => this%tau_ij(:,:,:,3) ! This always going to be zero unless populated
  this%tau_22   => this%tau_ij(:,:,:,4)
  this%tau_23C  => this%tau_ij(:,:,:,5) ! This always going to be zero unless populated
  this%tau_33   => this%tau_ij(:,:,:,6)
  allocate(this%tau_13(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
  allocate(this%tau_23(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
  this%tau_ij = zero

  if (present(initSpinUp)) then
      this%initSpinUp = initSpinUp
  end if

  if (this%isStratified .or. this%initSpinUp) then
     allocate(this%q1C(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
     allocate(this%q2C(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
     allocate(this%q3E(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
     this%q1C = zero; this%q2C = zero; this%q3E = zero
  end if 

  this%isPeriodic = this%PadeDer%isPeriodic


  ! Link buffers
  this%cbuffyC => cbuffyC
  this%cbuffzC => cbuffzC
  this%cbuffyE => cbuffyE
  this%cbuffzE => cbuffzE
  this%rbuffxC => rbuffxC
  this%rbuffyE => rbuffyE
  this%rbuffzE => rbuffzE
  this%rbuffyC => rbuffyC
  this%rbuffzC => rbuffzC

  open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
  read(unit=123, NML=SGS_MODEL)
  close(123)

  this%useScalarBounding = useScalarBounding 
  this%Cy = Cy  
  this%lowbound = lowbound
  this%highbound = highbound 
  this%UseDynamicProcedureScalar = UseDynamicProcedureScalar
  this%explicitCalcEdgeEddyViscosity = explicitCalcEdgeEddyViscosity
  this%mid = SGSmodelID
  this%z0  = z0
  this%z0t = z0t
  this%DynamicProcedureType = DynamicProcedureType
  this%DynProcFreq = DynProcFreq
  this%useVerticalTfilter = useVerticalTfilter

  this%isInviscid = isInviscid

  this%is_z0_varying  = is_z0_varying
  if(this%is_z0_varying) then
    ! read more stuff from the input file
    open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
    read(unit=123, NML=Z0VARYING)
    close(123)

    if((abs(spy-Ly)>1.0d-9) .or. (abs(rpy-Ly)>1.0d-9)) then
        call GracefulExit("Only stripe heterogeneity allowed. Change spy.", 11)
    endif

    if((abs(spnumx*(spx+rpx)-Lx)>1.0d-9)) then
        call GracefulExit("Inconsistent spx, rpx, spnumx and Lx. Check details.", 11)
    endif

    allocate(this%z0var(this%gpC%xsz(1), this%gpC%xsz(2)))
    allocate(this%ustarsqvar(this%gpC%xsz(1), this%gpC%xsz(2)))
    allocate(this%Uxvar(this%gpC%xsz(1), this%gpC%xsz(2)))
    allocate(this%alpfac(this%gpC%xsz(1), this%gpC%xsz(2)))
    allocate(this%Uyvar(this%gpC%xsz(1), this%gpC%xsz(2)))
    allocate(this%WallMFactorvar(this%gpC%xsz(1), this%gpC%xsz(2)))
    allocate(this%lamfact(this%gpC%xsz(1), this%gpC%xsz(2)))
    allocate(this%mask_upstream(this%gpC%xsz(1), this%gpC%xsz(2)))
    allocate(this%deli(this%gpC%xsz(1), this%gpC%xsz(2)))
    allocate(this%Uxspan(this%gpC%xsz(1), this%gpC%xsz(2)))
    allocate(this%Uyspan(this%gpC%xsz(1), this%gpC%xsz(2)))
    allocate(this%ustarsqspan(this%gpC%xsz(1), this%gpC%xsz(2))) !!!!!

    ! set all arrays to zero :: this is important
    this%z0var         = 0.0_rkind;   this%ustarsqvar     = 0.0_rkind;   this%Uxvar   = 0.0_rkind
    this%Uyvar         = 0.0_rkind;   this%WallMFactorvar = 0.0_rkind;   this%lamfact = 0.0_rkind
    this%mask_upstream = 0.0_rkind;   this%deli           = 0.0_rkind;   sp_map       = 0.0_rkind
    this%Uxspan        = 0.0_rkind;   this%Uyspan         = 0.0_rkind;   this%ustarsqspan = 0.0_rkind

    this%z0s = z0s; this%z0r = z0r
    this%filter_for_heterog = filter_for_heterog
    this%betfac = betfac
    totpx = spx + rpx
    if(rpstart < zero) then
        ! overwrite rough patch start with length of smooth patch
        rpstart = spx
    elseif(rpstart > totpx) then
        !call message(1, "!!!!!WARNING!!!!! rpstart is larger than totpx. Ensure this is a developing case.")
        call GracefulExit("rpstart is larger than totpx. Check details.", 11)
    endif

    if(spx_delta < epssmall) then
        do i = 1, this%gpC%xsz(1)
          xl = xMesh(i)
          xlmod = mod(xl, totpx)
          if((xlmod >= (rpstart+rpx)) .or. (xlmod < rpstart)) then
              sp_map(i) = 0.0_rkind !this%z0var(i,j) = z0s
          else
              sp_map(i) = 1.0_rkind !this%z0var(i,j) = z0r
          endif
        enddo
    else
        ! handle xpos = (rpstart-spx) separately
        xpos = rpstart-spx; x1 = (xMesh - (xpos + 0.5_rkind*spx_delta))/spx_delta + 1.0_rkind
        call S_fringe(x1, S1)
        sp_map = 1.0_rkind-S1

        !xpos = rpstart - spx; sp_map = 0.0_rkind

        ! then calculate map for all stripes
        do i = 1, spnumx
          xpos = xpos + spx; x1 = (xMesh - (xpos - 0.5_rkind*spx_delta))/spx_delta
          xpos = xpos + rpx; x2 = (xMesh - (xpos + 0.5_rkind*spx_delta))/spx_delta + 1.0_rkind
          call S_fringe(x1, S1)
          call S_fringe(x2, S2)
          sp_map = sp_map + S1 - S2
        enddo
    endif

    do j = 1, this%gpC%xsz(2)
      this%z0var(:,j) = z0s + sp_map * (z0r - z0s)
    enddo

    ! NOTE :: rpstart must be figured out before this block for setting alpfac
    this%alpfac = const_alpfac
    if(.not. use_const_alpfac) then
        do j = 1, this%gpC%xsz(2)
          do i = 1, this%gpC%xsz(1)
              xl = xMesh(i) - rpstart
              if(xl > epssmall) then
                  this%alpfac(i,j) = 0.004551_rkind * (xl/blht_for_alpfac)**(-0.3593)
              endif
          enddo
        enddo
    endif

    !!! only for Abkar-PA model. Think where to put this later
    !!! determine lamfact
    !!! smooth patch before rpstart and after fringe. Rough between these.
    Mfactor = 0.75_rkind-0.03_rkind*log(this%z0r/this%z0s)
    deli = 0.0_rkind
    do j = 1, this%gpC%xsz(2)
      do i = 1, this%gpC%xsz(1)
        if((xMesh(i) < rpstart) .or. (xMesh(i) > rpstart+spx)) then
            ! upstream smooth patch
            this%lamfact(i,j) = 1.1d0 ! set some number above 1.0
        else
            ! downstream rough or smooth region
            xl = xMesh(i) - rpstart
            deli(i) = this%z0r*Mfactor*(xl/this%z0r)**0.8_rkind
            this%lamfact(i,j) = -log(half*this%dz/(deli(i)*this%alpfac(i,j)+1.0d-18)) / log(this%alpfac(i,j))
            this%lamfact(i,j) = min(this%lamfact(i,j), one)     ! note :: must be here, not outside the loop
            this%lamfact(i,j) = max(this%lamfact(i,j), zero)    ! note :: must be here, not outside the loop
            write(100+nrank,'(2(i5,1x),9(e19.12,1x))') i, j, xl, deli(i), this%alpfac(i,j), half*this%dz, half*this%dz/(deli(i)*this%alpfac(i,j)+1.0d-18), -log(this%alpfac(i,j)), this%lamfact(i,j)
        endif
      enddo
    enddo

    if(this%gpC%xst(3)==1) then
        do j = 1, this%gpC%xsz(2)
          do i = 1, this%gpC%xsz(1)
            if((xMesh(i) < (rpstart-excludedist)) .or. (xMesh(i) > (rpstart+spx+excludedist))) then
                this%mask_upstream(i,j) = one
            endif
          enddo
        enddo
    endif
    this%mask_normfac = p_sum(sum(this%mask_upstream))

    do j=1, this%gpC%xsz(2)
      this%deli(:,j) = deli(:)
    enddo

    call message(1, "Printing debug info about heterog ")
    if(nrank==0) then
      !print *, 'nx = ', this%gpC%xsz(1)
      open(10, file='debug_heterog.dat',action='write',status='unknown')
      do i=1,this%gpC%xsz(1)
        write(10,'(5(e19.12,1x))') xMesh(i), sp_map(i), this%z0var(i,1), this%lamfact(i,1), deli(i)
      enddo
      close(10)
    endif

    ! write all data as 2D arrays
    this%rbuffxC(:,:,1,1) = this%z0var(:,:)
    this%rbuffxC(:,:,2,1) = this%mask_upstream(:,:)
    this%rbuffxC(:,:,3,1) = this%lamfact(:,:)
    this%rbuffxC(:,:,4,1) = this%deli(:,:)
    write(tempname,"(A)") "z0var_setup.dat"
    fname = "./"//trim(tempname)
    call decomp_2d_write_one(1, this%rbuffxC(:,:,:,1), fname, gpC)
    call message(1, "Done printing debug info about heterog ")

  endif

  this%WallModel  = WallModelType
  this%WM_matchingIndex = WM_matchingIndex
  if (this%WallModel .ne. 0) then
      if (this%PadeDer%isPeriodic) then
         call GracefulExit("You cannot use a wall model if the problem is periodic in Z",12)
      else
         call this%initWallModel()
      end if 
  else
      this%useWallModel = .false. 
  end if


  allocate(this%cmodelC(size(zMeshC)))
  allocate(this%cmodelE(size(zMeshE)))

  select case (SGSmodelID)
  case (0)
     call this%init_smagorinsky(dx,dy,dz,Csgs,ncWall,z0,useWallDamping,zMeshC, zMeshE)
  case (1)
     call this%init_sigma(dx, dy, dz, Csgs)
  case (2)
     call this%init_AMD(dx, dy, dz, Csgs)
  case default
     call GracefulExit("Incorrect choice for SGS model ID.", 213)
  end select

  if (this%isEddyViscosityModel) call this%allocateMemory_EddyViscosity()
  
  if (this%useScalarBounding) then 
      allocate(this%kappa_boundingC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
      allocate(this%kappa_boundingE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
      ierr = this%gaussianX%init(gpC%xsz(1), .true.)
      ierr = this%gaussianY%init(gpC%ysz(2), .true.)
      ierr = this%gaussianZ%init(gpC%zsz(3), this%isPeriodic)
      this%lowbound_PotT = lowbound
      this%highbound_PotT = highbound 
      this%Cy_PotT        = Cy 
  end if 

  if (this%isStratified) then
      this%TurbPrandtlNum_PotT = turbPrandtl
  end if

  if (DynamicProcedureType .ne. 0) then
      call this%allocateMemory_DynamicProcedure(computeFbody, deltaRatio)
      if (DynamicProcedureType .ne. 1) then
         call gracefulExit("Only planar averaged dynamic procedure is allowed; &
                  & all other code has been temporarily redacted. Contact Aditya & 
                  & (aditya90@stanford.edu) if you want to use any other method.",4234)
      end if
      if(useSGSDynamicRestart) then
         call this%readSGSDynamicRestart(SGSDynamicRestartFile)
      endif

  endif

  if ((UseDynamicProcedureScalar) .and. (DynamicProcedureType == 0)) then
      call gracefulExit("You cannot use dynamic procedure for a scalar without & 
               & using dynamic procedure for momentum.",123)
  end if 
  if ((this%mid .ne. 0) .and. (UseDynamicProcedureScalar)) then
      call gracefulExit("Dynamic procedure for scalar  only supported for smagorinsky",312)
  end if

end subroutine

pure subroutine S_fringe(x, output)
   real(rkind), dimension(:), intent(in)    :: x
   real(rkind), dimension(:), intent(out)   :: output
   integer :: i
   real(rkind) :: exparg

   do i = 1,size(x)
     if (x(i) .le. 0.d0) then
        output(i) = 0.d0
     else if (x(i) .ge. 1.d0) then
        output(i) = 1.d0
     else
        exparg = 1.d0/(x(i) - 1.d0 + 1.0D-32) + 1.d0/(x(i) + 1.0D-32)
        exparg = min(exparg,708.0d0) ! overflows if exparg > 709. need a better fix for this
        output(i) = 1.d0/(1.d0 + exp(exparg))
     end if
   end do

end subroutine

subroutine readSGSDynamicRestart(this,SGSDynamicRestartFile)
  class(sgs_igrid), intent(inout) :: this
  character(len=clen), intent(in) :: SGSDynamicRestartFile 
  integer :: ierr
  integer :: oldDynProcType, oldSGSmodel

  call GracefulExit("Restarting dynamic procedure using previous history has &
      & been temporarily redacted. Contact ADITYA if you want to know why.", 214)
  ! Open the restart file 
  open(unit=123, file=trim(SGSDynamicRestartFile), form='FORMATTED', status='old', action='read', iostat=ierr)
   
  ! Read in the mstep
  read(123,*) this%mstep

  ! Read in the old DynamicProcedure type 
  read(123,*) oldSGSmodel
  
  ! Read in the old DynamicProcedure type 
  read(123,*) oldDynProcType

  if (oldSGSmodel .ne. this%mid) then
      call message(1,"WARNING: Mismatch in the SGS model type between the restart file and new execution. &
                      & Will regenerate model constants at first new time step.")
      this%mstep = 0
  end if

  if (oldDynProcType .ne. this%DynamicProcedureType) then
      call message(1, "WARNING: Mismatch in the dynamic procedure type between the restart file and new execution.&
                       & Will regenerate model constants at first new time step.")
      this%mstep = 0
  else ! Good everything is matching
      select case( this%DynamicProcedureType ) 
      case (1) ! Planar averaged standard dynamic procedure
         ! This one is tricky because of the parallel processors (will complete
         ! later)
      case (2)
         read(123,*) this%mstep
      end select
  end if
  close(123)
end subroutine

subroutine dumpSGSDynamicRestart(this, SGSDynamicRestartFile)
  class(sgs_igrid), intent(in) :: this
  character(len=clen), intent(in) :: SGSDynamicRestartFile 
  integer :: ierr

  open(unit=123, file=trim(SGSDynamicRestartFile), form='FORMATTED', status='replace', iostat=ierr)
  write(123,*) this%mstep
  write(123,*) this%mid
  write(123,*) this%DynamicProcedureType

  select case (this%DynamicProcedureType) 
  case (1)
    ! This one is tricky because of the parallel processors (will complete
    ! later)
  case (2)
      write(123,*) this%cmodel_global
  end select
  close(123)

end subroutine
