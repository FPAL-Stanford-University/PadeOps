!!!!!NOTES FOR ADITYA: 
!!!!! 1. `type sgs_scalar', `subroutine init' need to be finished or trimmed. 
!!!!! 2. subroutine destropy needs to be added
!!!!! 3. arguments to all subroutine calls need fixing
!!!!! 4. linking to igrid


module sgsmod_scalar
    use kind_parameters, only: rkind, clen
    use constants, only: imi, pi, zero,one,two,three,half, four,eight, nine, six, kappa 
    use decomp_2d
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral  
    use mpi 
    use reductions, only: p_maxval, p_sum, p_minval
    use numerics, only: useCompactFD 
    !use StaggOpsMod, only: staggOps  
    use gaussianstuff, only: gaussian
    use lstsqstuff, only: lstsq
    use PadeDerOps, only: Pade6stagg
    implicit none

    private
    public :: sgs_scalar

    complex(rkind), parameter :: zeroC = zero + imi*zero
    logical :: useVerticalTfilter = .false. 
    real(rkind), parameter :: beta_h = 7.8_rkind, beta_m = 4.8_rkind

    type :: sgs_scalar !----this needs cleaning up--------
        private !-- how do I make only DynamicProcedureType and cmodel_gloal
        !public, so that it can be accessed from problem_files/temporalHooks.F90 ???
        type(decomp_info), pointer :: gpC, gpE
        type(spectral), pointer :: spectC, spectE
        type(decomp_info), pointer :: sp_gpC, sp_gpE
        integer :: mid, DynamicProcedureType, WallModel, DynProcFreq
        real(rkind), dimension(:), allocatable :: cmodelC, cmodelE
        real(rkind) :: cmodel_global, cmodel_global_x, cmodel_global_y, cmodel_global_z
        real(rkind), dimension(:,:,:), allocatable :: nu_sgs_C, nu_sgs_E
        real(rkind), dimension(:,:,:), allocatable :: kappa_sgs_C, kappa_sgs_E
        logical :: isEddyViscosityModel = .false.
        real(rkind), dimension(:,:,:,:), allocatable :: tau_ij
        real(rkind), dimension(:,:,:), pointer :: tau_11, tau_12, tau_22, tau_33, tau_13C, tau_23C
        real(rkind), dimension(:,:,:), allocatable :: tau_13, tau_23
        real(rkind), dimension(:,:,:,:), allocatable :: S_ij_C, S_ij_E
        real(rkind), dimension(:,:,:,:), pointer :: rbuffxE, rbuffxC, rbuffzC, rbuffyC, rbuffyE, rbuffzE
        complex(rkind), dimension(:,:,:,:), pointer :: cbuffyC, cbuffzC, cbuffyE, cbuffzE
        type(Pade6stagg), pointer :: PadeDer
        logical :: explicitCalcEdgeEddyViscosity = .false.
        real(rkind), dimension(:,:,:), allocatable :: q1C, q2C, q3E 
        logical :: initspinup = .false., isPeriodic = .false.  

        ! Wall model
        real(rkind), dimension(:,:,:,:), allocatable :: tauijWM
        complex(rkind),dimension(:,:,:,:), allocatable :: tauijWMhat_inZ, tauijWMhat_inY
        real(rkind), dimension(:,:,:), allocatable :: filteredSpeedSq
        complex(rkind), dimension(:,:,:), allocatable :: Tfilhat, Tfilhatz1, Tfilhatz2
        logical :: useWallModel = .false. 
        integer :: botBC_temp = 1
        real(rkind), public :: ustar = 1.d0, InvObLength = 0.d0
        real(rkind) :: umn = 1.d0, vmn = 1.d0, uspmn = 1.d0, Tmn = 1.d0, wTh_surf = 0.d0
        real(rkind) :: dz, z0, meanfact, ThetaRef, Fr, WallMfactor, Re, Pr
        real(rkind), pointer :: Tsurf
        complex(rkind), dimension(:,:), allocatable :: q3HAT_AtWall

        ! for dynamic procedures - all are at edges
        type(gaussian) :: gaussianTestFilterZ, gfiltx, gfilty, gfiltz
        real(rkind), dimension(:,:,:,:), allocatable :: Mij, Lij, Sij_Filt, alphaij_Filt, tauijWM_Filt
        real(rkind), dimension(:,:,:,:), allocatable :: fi_Filt, ui_Filt, fiE
        real(rkind), dimension(:,:,:),   allocatable :: Dsgs_Filt, buff1, buff2, cmodelC_local, cmodelE_local
        real(rkind), dimension(:,:,:),   pointer     :: fxC, fyC, fzE
        real(rkind), dimension(:,:,:),   pointer     :: Dsgs
        logical :: isInviscid, isStratified, useDynamicProcedure, useVerticalTfilter = .false. 
        real(rkind), dimension(:), allocatable :: cmodel_allZ
        real(rkind) :: invRe, deltaRat
        integer :: mstep, averageType = 0

        ! model constant values/properties
        real(rkind) :: camd_x, camd_y, camd_z
        !logical :: useCglobal = .false.
        integer :: modelConstType = 0 ! 0 :: scalar; 1 :: z-vector; 2 :: xyz array

        contains 
            procedure          :: init
            procedure          :: destroy
            procedure          :: getRHS_SGS_Scalar
            procedure, private :: getQjSGS
            procedure, private :: GetPrSGSStdDyn
            procedure, private :: GetPrSGSGlobDyn

    end type 

contains
!#include "sgs_models/init_destroy_sgs_igrid.F90"
!#include "sgs_models/smagorinsky.F90"
!#include "sgs_models/sigma.F90"
!#include "sgs_models/AMD.F90"
!#include "sgs_models/eddyViscosity.F90"
!#include "sgs_models/dynamicProcedure_sgs_igrid.F90"
!#include "sgs_models/standardDynamicProcedure.F90"
!#include "sgs_models/globalDynamicProcedure.F90"
!#include "sgs_models/wallmodel.F90"
!#include "sgs_models/accessors.F90"

subroutine init(this, gpC, gpE, spectC, spectE, dx, dy, dz, inputfile, zMeshE, zMeshC, fBody_x, fBody_y, fBody_z, computeFbody, PadeDer, cbuffyC, cbuffzC, cbuffyE, cbuffzE, rbuffxC, rbuffyC, rbuffzC, rbuffxE, rbuffyE, rbuffzE, Tsurf, ThetaRef, Fr, Re, Pr, isInviscid, isStratified, botBC_temp, initSpinUp, momSgsModel)
  class(sgs_igrid), intent(inout), target :: this
  class(decomp_info), intent(in), target :: gpC, gpE
  class(spectral), intent(in), target :: spectC, spectE
  real(rkind), intent(in) :: dx, dy, dz, ThetaRef, Fr, Re, Pr
  real(rkind), intent(in), target :: Tsurf
  character(len=*), intent(in) :: inputfile
  real(rkind), dimension(:), intent(in) :: zMeshE, zMeshC
  real(rkind), dimension(:,:,:), intent(in), target :: fBody_x, fBody_y, fBody_z
  logical, intent(out) :: computeFbody
  real(rkind), dimension(:,:,:,:), intent(in), target :: rbuffxC, rbuffyE, rbuffzE, rbuffyC, rbuffzC, rbuffxE
  complex(rkind), dimension(:,:,:,:), intent(in), target :: cbuffyC, cbuffzC, cbuffyE, cbuffzE
  type(Pade6stagg), target, intent(in) :: PadeDer
  logical, intent(in) :: isInviscid, isStratified
  integer, intent(in) :: botBC_temp
  logical, intent(in), optional :: initSpinUp
  type(sgsmod_igrid), intent(inout), target :: momSgsModel

  ! Input file variables
  logical :: useWallDamping = .false., useSGSDynamicRestart = .false., useVerticalTfilter = .false.
  integer :: DynamicProcedureType = 0, SGSmodelID = 0, WallModelType = 0, DynProcFreq = 1 
  real(rkind) :: ncWall = 1.d0, Csgs = 0.17d0, z0 = 0.01d0
  character(len=clen) :: SGSDynamicRestartFile
  logical :: explicitCalcEdgeEddyViscosity = .false.
  integer :: ierr, averageType = 0
  
  namelist /SGS_MODEL/ DynamicProcedureType, SGSmodelID, z0,  &
                 useWallDamping, ncWall, Csgs, WallModelType, &
                 DynProcFreq, useSGSDynamicRestart, useVerticalTfilter,  &
                 SGSDynamicRestartFile,explicitCalcEdgeEddyViscosity,    &
                 averageType


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
  this%averageType = averageType

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
  this%rbuffxE => rbuffxE
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

  if(this%useVerticalTFilter .and. this%isPeriodic) then
      call GracefulExit("You cannot use separate vertical test filter if the problem is periodic in Z",12)
  endif

 
  this%isInviscid = isInviscid

  this%WallModel  = WallModelType
 
  ! --- need to work on this-- 
  !if (this%WallModel .ne. 0) then
  !    if (this%PadeDer%isPeriodic) then
  !       call GracefulExit("You cannot use a wall model if the problem is periodic in Z",12)
  !    else
  !       call this%initWallModel()
  !    end if 
  !else
  !    this%useWallModel = .false. 
  !end if

  if (this%isEddyViscosityModel) then
      allocate(this%kappa_sgs_C(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
      allocate(this%kappa_sgs_E(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)))
  end if
  
  if (DynamicProcedureType .ne. 0) then
      call this%allocateMemory_DynamicProcedure_scalar(computeFbody, dx, dy, dz)
      !if(useSGSDynamicRestart) then
      !   call this%readSGSDynamicRestart(SGSDynamicRestartFile)
      !endif
  else
      !if(useWallDamping) then
      !   call this%initWallDamping(dx, dy, dz, Csgs,z0,ncWall,zMeshC, zMeshE) 
      !endif
  endif

end subroutine

subroutine allocateMemory_DynamicProcedure_scalar(this)
   class(sgs_scalar), intent(inout), target :: this

   this%useDynamicProcedure = .true.
   select case( this%dynamicProcedureType) 
   case (1)
      this%Pj => this%momSgsModel%alphaij_Filt(:,:,:,1:3)
      this%Rj => this%momSgsModel%alphaij_Filt(:,:,:,4:6)
      this%TE_Filt => this%momSgsModel%alphaij_Filt(:,:,:,7)
      this%buff1 => this%momSgsModel%alphaij_Filt(:,:,:,8)
      this%buff2 => this%momSgsModel%alphaij_Filt(:,:,:,9)
      this%ui_Filt => this%momSgsModel%ui_Filt
      this%nusgs_E => this%momSgsModel%nusgs_E
      this%nusgs_E_Filt => this%momSgsModel%Mij(:,:,:,1)

      select case(this%AverageType)
      case(0) ! Volume 
         this%PrSGSType = 0
         if(this%isPeriodic) then
           call this%spectC%init_TestFilt(this%deltaRat, dx, dy, dz)
         else
           call this%spectE%init_TestFilt(this%deltaRat, dx, dy, dz)
         endif
      case(1) ! Planar
         if(this%isPeriodic) then
           call GracefulExit("You cannot use planar averaging if the problem is periodic in Z",12)
         endif
         this%PrSGSType = 1
         allocate(this%InvPrSGSZ_allZ(this%gpE%zsz(3)))
         allocate(this%InvPrSGSZ_C(this%gpC%xsz(3)))
         allocate(this%InvPrSGSZ_E(this%gpE%xsz(3)))
         call this%spectE%init_TestFilt(this%deltaRat, dx, dy, dz)
      case(2)
         this%PrSGSType = 2
         allocate(this%InvPrSGSC_local(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
         allocate(this%InvPrSGSE_local(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)))
         ierr = this%gfiltx%init(this%gpE%xsz(1), .true.)
         ierr = this%gfilty%init(this%gpE%ysz(2), .true.)
         ierr = this%gfiltz%init(this%gpE%zsz(3), this%isPeriodic)
         if(this%isPeriodic) then
           call this%spectC%init_TestFilt(this%deltaRat, dx, dy, dz)
         else
           call this%spectE%init_TestFilt(this%deltaRat, dx, dy, dz)
         endif
      end select
   case (2)
         this%PrSGSType = 0
   end select

end subroutine


subroutine getRHS_SGS_Scalar(this, Trhs, dTdxC, dTdyC, dTdzE)
   class(sgs_scalar), intent(inout), target :: this
   complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout) :: Trhs
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: dTdxC, dTdyC
   real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)), intent(in) :: dTdzE
   complex(rkind), dimension(:,:,:), pointer :: cbuffy1, cbuffy2, cbuffz1, cbuffz2

   cbuffy1 => this%cbuffyC(:,:,:,1); cbuffy2 => this%cbuffyE(:,:,:,1); 
   cbuffz1 => this%cbuffzC(:,:,:,1); cbuffz2 => this%cbuffzE(:,:,:,1) 

   ! First get qj's 
   call this%getQjSGS(dTdxC, dTdyC, dTdzE)

   ! ddx(q1)
   call this%spectC%fft(this%q1C, cbuffy1)
   call this%spectC%mtimes_ik1_ip(cbuffy1)
   Trhs = Trhs - cbuffy1

   ! ddy(q2)
   call this%spectC%fft(this%q2C, cbuffy1)
   call this%spectC%mtimes_ik2_ip(cbuffy1)
   Trhs = Trhs - cbuffy1

   ! ddz(q3)
   call this%spectE%fft(this%q3E, cbuffy2)
   call transpose_y_to_z(cbuffy2, cbuffz2, this%sp_gpE)
   call this%PadeDer%ddz_E2C(cbuffz2, cbuffz1, 0, 0)
   call transpose_z_to_y(cbuffz1,cbuffy1,this%sp_gpC)
   Trhs = Trhs - cbuffy1

end subroutine

subroutine getQjSGS(this,dTdxC, dTdyC, dTdzE)
   class(sgs_scalar), intent(inout), target :: this
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: dTdxC, dTdyC
   real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)), intent(in) :: dTdzE

   if (this%useWallModel) call this%computeWall_PotTFlux()

   if (this%isEddyViscosityModel) then

      if(newTimeStep .and. this%useDynamicProcedure) then
          if (mod(this%mstep, this%DynProcFreq) == 0) then

              select case(this%dynamicProcedureType) 
              case (1) ! Standard (planar averaged) dynamic procedure
                 call this%StandardDynamicProcedure_scalar()     !uE, vE, wE, uhatE, vhatE, whatE, duidxjEhat)
              case (2) ! Global Dynamic Procedure
                 call this%GlobalDynamicProcedure_scalar()     !uhatE, vhatE, whatE, uE, vE, wE, duidxjEhat, duidxjE) ! Pass in the relevant stuff, finish the procedure implementation
              end select 
          endif
      endif

      select case (this%PrSGSType) 
      case(0) ! Constant or global dynamic procedure
        this%kappa_sgs_C = this%InvPrSGS*this%nu_sgs_C
        this%kappa_sgs_E = this%InvPrSGS*this%nu_sgs_E
      case(1) ! vector
        do k = 1,size(this%kappa_sgs_C,3)
           this%kappa_sgs_C(:,:,k) = this%InvPrSGSZ_C(k)*this%nu_sgs_C(:,:,k)
        end do 
        do k = 1,size(this%kappa_sgs_E,3)
           this%kappa_sgs_E(:,:,k) = this%InvPrSGSZ_E(k)*this%nu_sgs_E(:,:,k)
        end do 
      case(2) ! 3D field
        this%kappa_sgs_E = this%InvPrSGSE_local*this%nu_sgs_E
        this%kappa_sgs_C = this%InvPrSGSC_local*this%nu_sgs_E
      end select

      this%q1C = -this%kappa_sgs_C*dTdxC
      this%q2C = -this%kappa_sgs_C*dTdyC
      this%q3E = -this%kappa_sgs_E*dTdzE
   end if

   if (this%useWallModel) call this%embed_WM_PotTflux()
 
end subroutine


subroutine GlobalDynamicProcedure_scalar(this)     !, uE, vE, wE, uhatE, vhatE, whatE, duidxjEhat)
   class(sgs_igrid), intent(inout) :: this
 
end subroutine

subroutine StandardDynamicProcedure_scalar(this)     !, uE, vE, wE, uhatE, vhatE, whatE, duidxjEhat)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)), intent(in) :: uE, vE, wE
   complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3))  , intent(in) :: uhatE, vhatE, whatE
   complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3),9), intent(in) :: duidxjEhat
   integer :: idx
   real(rkind) :: tmp1, tmp2

   ! STEP 1: Test filter velocities (required for Lij)
   call this%TestFilter_Cmplx_to_Real(ThatE, this%TE_Filt)

   ! STEP 2: Compute Pj
   this%buff1 = uE*TE
   call this%TestFilter_real_to_real(this%buff1, this%buff2)
   this%Pj(:,:,:,1) = -this%buff2 + this%ui_Filt(:,:,:,1)*this%TE_Filt
   
   this%buff1 = vE*TE
   call this%TestFilter_real_to_real(this%buff1, this%buff2)
   this%Pj(:,:,:,2) = -this%buff2 + this%ui_Filt(:,:,:,2)*this%TE_Filt

   this%buff1 = wE*TE
   call this%TestFilter_real_to_real(this%buff1, this%buff2)
   this%Pj(:,:,:,3) = -this%buff2 + this%ui_Filt(:,:,:,3)*this%TE_Filt

   ! STEP 3: Compute R_j
   ! Part a: Compute \tilde{dTdxj}
   call this%TestFilter_cmplx_to_real(dTdxEhat,this%Rj(:,:,:,1))
   call this%TestFilter_cmplx_to_real(dTdyEhat,this%Rj(:,:,:,2))
   call this%TestFilter_cmplx_to_real(dTdzEhat,this%Rj(:,:,:,3))

   ! Part b: Compute \tilde{nusgs_E}
   call this%TestFilter_real_to_real(this%nusgs_E, this%nusgs_E_Filt)

   ! Part c: Compute the rest
   this%buff1 = this%nusgs_E*dTdxE
   call this%TestFilter_real_to_real(this%buff1, this%buff2)
   this%Rj(:,:,:,1) = this%Rj(:,:,:,1)*this%nusgs_E_filt - this%buff2

   this%buff1 = this%nusgs_E*dTdyE
   call this%TestFilter_real_to_real(this%buff1, this%buff2)
   this%Rj(:,:,:,2) = this%Rj(:,:,:,2)*this%nusgs_E_filt - this%buff2

   this%buff1 = this%nusgs_E*dTdzE
   call this%TestFilter_real_to_real(this%buff1, this%buff2)
   this%Rj(:,:,:,3) = this%Rj(:,:,:,3)*this%nusgs_E_filt - this%buff2

   ! STEP 4: Compute the numerator
   this%buff1 = this%Pj(:,:,:,1)*this%Rj(:,:,:,1)
   do idx = 2,3
      this%buff1 = this%buff1 + this%Pj(:,:,:,idx)*this%Rj(:,:,:,idx)
   end do 

   ! STEP 5: Compute the denominator
   this%buff1 = this%Rj(:,:,:,1)*this%Rj(:,:,:,1)
   do idx = 2,3
      this%buff1 = this%buff1 + this%Rj(:,:,:,idx)*this%Rj(:,:,:,idx)
   end do 

   ! STEP 6: Get the planar average and interpolate
   select case(this%PrSGSType)
   case(0)
      tmp1 = p_sum(sum(this%buff1))
      tmp2 = p_sum(sum(this%buff2))
      this%InvPrSGS = tmp1/(tmp2+1.0D-14)
   case(1) ! Planar Average
     call this%planarAverageAndInterpolateToCells_scalar(this%buff1, this%buff2, this%rbuffxC(:,:,:,1))
     this%InvPrSGSZ_E = this%buff1(1,1,:)
     this%InvPrSGSZ_C = this%rbuffxC(1,1,:,1)
   case(2) ! Gaussian Filter
     where(this%buff1 < zero)
       this%buff1 = zero
     end where
     
     call this%gaussFilter3D_scalar(this%buff1)
     call this%gaussFilter3D_scalar(this%buff2)
     this%InvPrSGSE_local = this%buff1/(this%buff2 + 1.0D-14)

     call transpose_x_to_y(this%InvPrSGSE_local, this%rbuffyE(:,:,:,1), this%gpE)
     call transpose_y_to_z(this%rbuffyE(:,:,:,1), this%rbuffzE(:,:,:,1), this%gpE)
     this%rbuffzC(:,:,1:this%gpC%zsz(3),1) = 0.5d0*(this%rbuffzE(:,:,1:this%gpC%zsz(3),1)+this%rbuffzE(:,:,2:this%gpC%zsz(3)+1,1))
     call transpose_z_to_y(this%rbuffzC(:,:,:,1), this%rbuffyC(:,:,:,1), this%gpC)
     call transpose_y_to_x(this%rbuffyC(:,:,:,1), this%InvPrSGSC_local, this%gpC)
   end select
end subroutine


subroutine gaussFilter3D_scalar(this, fin)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)), intent(inout) :: fin

   call this%gfiltx%filter1(fin,                   this%rbuffxE(:,:,:,1), this%gpE%xsz(2), this%gpE%xsz(3), 0, 0)
   call transpose_x_to_y   (this%rbuffxE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
   call this%gfilty%filter2(this%rbuffyE(:,:,:,1), this%rbuffyE(:,:,:,2), this%gpE%ysz(1), this%gpE%ysz(3), 0, 0)
   call transpose_y_to_z   (this%rbuffyE(:,:,:,2), this%rbuffzE(:,:,:,1), this%gpE)
   call this%gfiltz%filter3(this%rbuffzE(:,:,:,1), this%rbuffzE(:,:,:,2), this%gpE%zsz(1), this%gpE%zsz(2), 0, 0)

   call transpose_z_to_y   (this%rbuffzE(:,:,:,2), this%rbuffyE(:,:,:,1), this%gpE)
   call transpose_y_to_x   (this%rbuffyE(:,:,:,1), fin,                   this%gpE)
end subroutine 


subroutine planarAverageAndInterpolateToCells_scalar(this, numE, denE, ratC)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)), intent(inout) :: numE
   real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)), intent(in)    :: denE
   real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)), intent(out)   :: ratC
   integer :: idx
   real(rkind) :: tmp1, tmp2  
 
   call transpose_x_to_y(numE, this%rbuffyE(:,:,:,1), this%gpE)
   call transpose_y_to_z(this%rbuffyE(:,:,:,1), this%rbuffzE(:,:,:,1), this%gpE)
   do idx = 1,this%gpE%zsz(3)
      this%rbuffzE(:,:,idx,1) = max(p_sum(sum(this%rbuffzE(:,:,idx,1)))*this%meanfact, zero)
   end do 
   
   call transpose_x_to_y(denE, this%rbuffyE(:,:,:,1), this%gpE)
   call transpose_y_to_z(this%rbuffyE(:,:,:,1), this%rbuffzE(:,:,:,2), this%gpE)
   do idx = 1,this%gpE%zsz(3)
      this%rbuffzE(:,:,idx,2) = p_sum(sum(this%rbuffzE(:,:,idx,2)))*this%meanfact
   end do
   this%rbuffzE(:,:,:,1) = this%rbuffzE(:,:,:,1)/(this%rbuffzE(:,:,:,2) + 1.d-14)
   this%cmodel_allZ = this%rbuffzE(1,1,:,1)

   this%rbuffzC(:,:,1:this%gpC%zsz(3),1) = 0.5d0*(this%rbuffzE(:,:,1:this%gpC%zsz(3),1)+this%rbuffzE(:,:,2:this%gpC%zsz(3)+1,1))
   call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
   call transpose_y_to_x(this%rbuffyE(:,:,:,1), numE, this%gpE)

   call transpose_z_to_y(this%rbuffzC(:,:,:,1), this%rbuffyC(:,:,:,1), this%gpC)
   call transpose_y_to_x(this%rbuffyC(:,:,:,1), ratC, this%gpC)
end subroutine


end module 
