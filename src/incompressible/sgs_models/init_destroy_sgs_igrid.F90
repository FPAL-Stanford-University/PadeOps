subroutine destroy(this)
  class(sgs_igrid), intent(inout) :: this
  nullify(this%gpC, this%gpE, this%spectC, this%spectE, this%sp_gpC, this%sp_gpE, this%fxC)
  nullify(this%cbuffyC, this%cbuffzC, this%rbuffxC, this%Tsurf, this%fyC, this%fzE)
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
end subroutine


subroutine link_pointers(this, nuSGS, tauSGS_ij, tau13, tau23, q1, q2, q3)
   class(sgs_igrid), intent(in), target :: this
   real(rkind), dimension(:,:,:)  , pointer, intent(inout) :: nuSGS
   real(rkind), dimension(:,:,:)  , pointer, intent(inout) :: tau13, tau23
   real(rkind), dimension(:,:,:,:), pointer, intent(inout) :: tauSGS_ij
   real(rkind), dimension(:,:,:)  , pointer, intent(inout) :: q1, q2, q3

   nuSGS => this%nu_sgs_C
   tau13 => this%tau_13
   tau23 => this%tau_23

   tauSGS_ij => this%tau_ij
  
   if (this%isStratified) then
      q1 => this%q1C
      q2 => this%q2C
      q3 => this%q3E
   end if 
end subroutine 

subroutine init(this, gpC, gpE, spectC, spectE, dx, dy, dz, inputfile, zMeshE, zMeshC, fBody_x, fBody_y, fBody_z, computeFbody, PadeDer, cbuffyC, cbuffzC, cbuffyE, cbuffzE, rbuffxC, rbuffyC, rbuffzC, rbuffyE, rbuffzE, Tsurf, ThetaRef, Fr, Re, Pr, isInviscid, isStratified, botBC_temp, initSpinUp)
  class(sgs_igrid), intent(inout), target :: this
  class(decomp_info), intent(in), target :: gpC, gpE
  class(spectral), intent(in), target :: spectC, spectE
  real(rkind), intent(in) :: dx, dy, dz, ThetaRef, Fr, Re, Pr
  real(rkind), intent(in), target :: Tsurf
  character(len=*), intent(in) :: inputfile
  real(rkind), dimension(:), intent(in) :: zMeshE, zMeshC
  real(rkind), dimension(:,:,:), intent(in), target :: fBody_x, fBody_y, fBody_z
  logical, intent(out) :: computeFbody
  real(rkind), dimension(:,:,:,:), intent(in), target :: rbuffxC, rbuffyE, rbuffzE, rbuffyC, rbuffzC
  complex(rkind), dimension(:,:,:,:), intent(in), target :: cbuffyC, cbuffzC, cbuffyE, cbuffzE
  type(Pade6stagg), target, intent(in) :: PadeDer
  logical, intent(in) :: isInviscid, isStratified
  integer, intent(in) :: botBC_temp
  logical, intent(in), optional :: initSpinUp

  ! Input file variables
  logical :: useWallDamping = .false., useSGSDynamicRestart = .false., useVerticalTfilter = .false.
  integer :: DynamicProcedureType = 0, SGSmodelID = 0, WallModelType = 0, DynProcFreq = 1 
  real(rkind) :: ncWall = 1.d0, Csgs = 0.17d0, z0 = 0.01d0
  character(len=clen) :: SGSDynamicRestartFile
  logical :: explicitCalcEdgeEddyViscosity = .false.
  integer :: ierr
  
  namelist /SGS_MODEL/ DynamicProcedureType, SGSmodelID, z0,  &
                 useWallDamping, ncWall, Csgs, WallModelType, &
                 DynProcFreq, useSGSDynamicRestart, useVerticalTfilter,           &
                 SGSDynamicRestartFile,explicitCalcEdgeEddyViscosity


  this%gpC => gpC
  this%gpE => gpE
  this%spectC => spectC
  this%spectE => spectE
  this%sp_gpC => spectC%spectdecomp
  this%sp_gpE => spectE%spectdecomp
  this%dz = dz
  this%Tsurf => Tsurf
  this%Fr = Fr
  this%Re = Re
  this%Pr = Pr
  this%ThetaRef = ThetaRef
  this%PadeDer => PadeDer
  this%fxC => fBody_x
  this%fyC => fBody_y
  this%fzE => fBody_z
  this%meanfact = one/(real(gpC%xsz(1),rkind) * real(gpC%ysz(2),rkind))
  this%isStratified = isStratified
  !if (present(botBC_Temp)) 
  this%botBC_Temp = botBC_Temp

  allocate(this%tau_ij(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3),6))
  this%tau_11   => this%tau_ij(:,:,:,1)
  this%tau_12   => this%tau_ij(:,:,:,2)
  this%tau_13C  => this%tau_ij(:,:,:,3) ! This always going to be zero
  this%tau_22   => this%tau_ij(:,:,:,4)
  this%tau_23C  => this%tau_ij(:,:,:,5) ! This always going to be zero 
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

  this%explicitCalcEdgeEddyViscosity = explicitCalcEdgeEddyViscosity
  this%mid = SGSmodelID
  this%z0 = z0
  this%DynamicProcedureType = DynamicProcedureType
  this%DynProcFreq = DynProcFreq
  this%useVerticalTfilter = useVerticalTfilter
  
  this%isInviscid = isInviscid

  this%WallModel  = WallModelType
  
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
  
  if (DynamicProcedureType .ne. 0) then
      call this%allocateMemory_DynamicProcedure(computeFbody)
      if(useSGSDynamicRestart) then
         call this%readSGSDynamicRestart(SGSDynamicRestartFile)
      endif
  endif


end subroutine

subroutine readSGSDynamicRestart(this,SGSDynamicRestartFile)
  class(sgs_igrid), intent(inout) :: this
  character(len=clen), intent(in) :: SGSDynamicRestartFile 
  integer :: ierr
  integer :: oldDynProcType, oldSGSmodel

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
