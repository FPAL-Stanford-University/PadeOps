subroutine destroy(this)
  class(sgs_igrid), intent(inout) :: this
  nullify(this%gpC, this%gpE, this%spectC, this%spectE, this%sp_gpC, this%sp_gpE, this%fiC)
  nullify(this%cbuffyC, this%cbuffzC, this%rbuffxC, this%Tsurf)
  if (this%isEddyViscosityModel) call this%destroyMemory_EddyViscosity()
  if (this%DynamicProcedureType .ne. 0) call this%destroyMemory_DynamicProcedure()
  select case (this%mid)
  case (0)
     call this%destroy_smagorinsky()
  case (1)
     call this%destroy_sigma()
  end select
  deallocate(this%tau_11, this%tau_12, this%tau_13, this%tau_22, this%tau_23, this%tau_33)
end subroutine

subroutine init(this, gpC, gpE, spectC, spectE, dx, dy, dz, inputfile, zMeshE, zMeshC, fBody, computeFbody, PadeDer, cbuffyC, cbuffzC, cbuffyE, cbuffzE, rbuffxC, rbuffyC, rbuffzC, rbuffyE, rbuffzE, Tsurf, ThetaRef, Fr, Re, isInviscid, isStratified)
  class(sgs_igrid), intent(inout) :: this
  class(decomp_info), intent(in), target :: gpC, gpE
  class(spectral), intent(in), target :: spectC, spectE
  real(rkind), intent(in) :: dx, dy, dz, ThetaRef, Fr, Re
  real(rkind), intent(in), target :: Tsurf
  character(len=*), intent(in) :: inputfile
  real(rkind), dimension(:), intent(in) :: zMeshE, zMeshC
  real(rkind), dimension(:,:,:,:), intent(in), target :: fBody
  logical, intent(out) :: computeFbody
  real(rkind), dimension(:,:,:,:), intent(in), target :: rbuffxC, rbuffyE, rbuffzE, rbuffyC, rbuffzC
  complex(rkind), dimension(:,:,:,:), intent(in), target :: cbuffyC, cbuffzC, cbuffyE, cbuffzE
  type(Pade6stagg), target, intent(in) :: PadeDer
  logical, intent(in) :: isInviscid, isStratified

  ! Input file variables
  logical :: useWallDamping = .false., useSGSDynamicRestart
  integer :: DynamicProcedureType = 0, SGSmodelID = 0, WallModelType = 0, DynProcFreq = 1 
  real(rkind) :: ncWall = 1.d0, Csgs = 0.17d0, z0 = 0.01d0
  character(len=clen) :: SGSDynamicRestartFile
  namelist /SGS_MODEL/ DynamicProcedureType, SGSmodelID, z0,  &
                 useWallDamping, ncWall, Csgs, WallModelType, &
                 DynProcFreq, useSGSDynamicRestart, SGSDynamicRestartFile

  integer :: ierr

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
  this%ThetaRef = ThetaRef
  this%PadeDer => PadeDer
  this%fiC => fBody
  this%meanfact = one/(real(gpC%xsz(1),rkind) * real(gpC%ysz(2),rkind))
  allocate(this%tau_11(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
  allocate(this%tau_12(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
  allocate(this%tau_22(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
  allocate(this%tau_33(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
  allocate(this%tau_13(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
  allocate(this%tau_23(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))

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

  this%mid = SGSmodelID
  this%z0 = z0
  this%DynamicProcedureType = DynamicProcedureType
  this%DynProcFreq = DynProcFreq
  this%WallModel        = WallModelType
  if (this%WallModel .ne. 0) call this%initWallModel()

  select case (SGSmodelID)
  case (0)
     call this%init_smagorinsky(dx,dy,dz,Csgs,ncWall,z0,useWallDamping,zMeshC, zMeshE)
  case (1)
     call this%init_sigma(dx, dy, dz, Csgs)
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

  this%isInviscid = isInviscid
  this%isStratified = isStratified

end subroutine

subroutine readSGSDynamicRestart(this,SGSDynamicRestartFile)
  class(sgs_igrid), intent(inout) :: this
  character(len=clen), intent(in) :: SGSDynamicRestartFile 
  integer :: ierr

  if(this%DynamicProcedureType==1) then
     ! read mstep
     open(unit=123, file=trim(SGSDynamicRestartFile), form='FORMATTED', status='old', action='read', iostat=ierr)
     read(123,*) this%mstep
     ! read cmodelC, cmodelE
     !! figure this out later
     !do k=1,size(cmodelC)
     !read(123,*)
     close(123)
  else
     ! read cmodel_global, mstep
     open(unit=123, file=trim(SGSDynamicRestartFile), form='FORMATTED', status='old', action='read', iostat=ierr)
     read(123,*) this%mstep
     read(123,*) this%cmodel_global
     close(123)
  endif

end subroutine

subroutine dumpSGSDynamicRestart(this, SGSDynamicRestartFile)
  class(sgs_igrid), intent(in) :: this
  character(len=clen), intent(in) :: SGSDynamicRestartFile 
  integer :: ierr

  if(this%DynamicProcedureType==1) then
     ! write mstep
     open(unit=123, file=trim(SGSDynamicRestartFile), form='FORMATTED', status='replace', iostat=ierr)
     write(123,*) this%mstep
     ! write cmodelC, cmodelE
     !! figure this out later
     !do k=1,size(cmodelC)
     !write(123,*)
     close(123)
  else
     ! write cmodel_global, mstep
     open(unit=123, file=trim(SGSDynamicRestartFile), form='FORMATTED', status='replace', iostat=ierr)
     write(123,*) this%mstep
     write(123,*) this%cmodel_global
     close(123)
  endif

end subroutine
