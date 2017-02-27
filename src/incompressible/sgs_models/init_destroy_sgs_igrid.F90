subroutine destroy(this)
  class(sgs_igrid), intent(inout) :: this
  nullify(this%gpC, this%gpE, this%spectC, this%spectE, this%sp_gpC, this%sp_gpE, this%fi)
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

subroutine init(this, gpC, gpE, spectC, spectE, dx, dy, dz, inputfile, zMeshE, zMeshC, fBody, computeFbody, cbuffyC, cbuffzC, rbuffxC, rbuffyE, rbuffzE, Tsurf, ThetaRef, Fr)
  class(sgs_igrid), intent(inout) :: this
  class(decomp_info), intent(in), target :: gpC, gpE
  class(spectral), intent(in), target :: spectC, spectE
  real(rkind), intent(in) :: dx, dy, dz, ThetaRef, Fr
  real(rkind), intent(in), target :: Tsurf
  character(len=*), intent(in) :: inputfile
  real(rkind), dimension(:), intent(in) :: zMeshE, zMeshC
  real(rkind), dimension(:,:,:,:), intent(in), target :: fBody
  logical, intent(out) :: computeFbody
  integer :: ierr
  real(rkind), dimension(:,:,:,:), intent(in), target :: rbuffxC, rbuffyE, rbuffzE
  complex(rkind), dimension(:,:,:,:), intent(in), target :: cbuffyC, cbuffzC

  ! Input file variables
  logical :: useWallDamping = .false.
  integer :: DynamicProcedureType = 0, SGSmodelID = 0, WallModelType = 0
  real(rkind) :: ncWall = 1.d0, Csgs = 0.17d0, z0 = 0.01d0
  namelist /SGS_MODEL/ DynamicProcedureType, SGSmodelID, z0, &
                 useWallDamping, ncWall, Csgs, WallModelType 


  this%gpC => gpC
  this%gpE => gpE
  this%spectC => spectC
  this%spectE => spectE
  this%sp_gpC => spectC%spectdecomp
  this%sp_gpE => spectE%spectdecomp
  this%dz = dz
  this%Tsurf => Tsurf
  this%Fr = Fr
  this%ThetaRef = ThetaRef

  open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
  read(unit=123, NML=SGS_MODEL)
  close(123)

  this%mid = SGSmodelID
  this%z0 = z0
  this%DynamicProcedureType = DynamicProcedureType
  this%WallModelType        = WallModelType
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
  
  if (DynamicProcedureType .ne. 0) call this%allocateMemory_DynamicProcedure(computeFbody)

  allocate(this%tau_11(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
  allocate(this%tau_12(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
  allocate(this%tau_22(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
  allocate(this%tau_33(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
  allocate(this%tau_13(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
  allocate(this%tau_23(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))

  this%fi => fBody
  this%meanfact = one/(real(gpC%xsz(1),rkind) * real(gpC%ysz(2),rkind))

  ! Link buffers
  this%cbuffyC => cbuffyC
  this%cbuffzC => cbuffzC
  this%rbuffxC => rbuffxC
  this%rbuffyE => rbuffyE
  this%rbuffzE => rbuffzE
end subroutine
