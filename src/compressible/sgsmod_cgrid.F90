module sgsmod_cgrid
    use kind_parameters, only: rkind, clen
    use constants,       only: imi, zero, two, eight, eps, third, one, two, three, half, four,eight, nine, six, kappa, pi
    use decomp_2d
    use gridtools,       only: alloc_buffs, destroy_buffs
    use exits, only: GracefulExit, message
    use mpi 
    use reductions, only: p_maxval, p_sum, p_minval
    implicit none

    private
    public :: sgs_cgrid

    complex(rkind), parameter :: zeroC = zero + imi*zero
    logical :: useVerticalTfilter = .false. 

    type :: sgs_cgrid
        private 
        class(decomp_info), pointer :: decomp
        integer :: SGSmodelID, mid, DynamicProcedureType, DynProcFreq
        real(rkind) :: cmodel_global
        logical :: isEddyViscosityModel, isEddyDiffmodel, isTurbPrandtlconst
        real(rkind), dimension(:,:,:,:), allocatable :: S_ij
        real(rkind), dimension(:,:,:),   allocatable :: nusgs, kapsgs
        logical ::  isPeriodic = .false.

        integer :: nxL, nyL, nzL
        real(rkind) :: dx, dy, dz, deltaLES

        real(rkind) :: Csgs, ncWall, PrSGS, Cp, Pr
        logical :: useWallDamping, useSGSDynamicRestartFile
        logical :: DomainAveraged_DynProc, useDynamicProcedureScalar
        character(len=clen) :: SGSDynamicRestartFile

        !! model constant values/properties
        real(rkind) :: camd_x, camd_y, camd_z, cmgm_x, cmgm_y, cmgm_z, c1_mgm , c2_mgm, PrCpfac
        logical :: useCglobal = .false. 
        contains 
            !! ALL INIT PROCEDURES
            procedure          :: init
            procedure, private :: init_smagorinsky
            procedure, private :: init_sigma
            procedure, private :: init_amd
            procedure, private :: init_mgm
            !! ALL SGS SOURCE/SINK PROCEDURES
            procedure          :: getQjSGS
            procedure          :: getTauSGS
            procedure, private :: get_SGS_kernel
            procedure, private :: get_smagorinsky_kernel
            procedure, private :: get_sigma_kernel
            procedure, private :: get_amd_kernel
            procedure, private :: get_amd_Dkernel
            procedure, private :: get_tausgs_mgm
            procedure, private :: get_Qjsgs_mgm
            procedure, private :: multiply_by_model_constant 
            procedure, private :: get_Sij_from_duidxj
            procedure, private :: get_modS_Sii
            !! ALL DESTROY PROCEDURES
            procedure          :: destroy
            procedure, private :: destroy_smagorinsky 
            procedure, private :: destroy_sigma
            procedure, private :: destroy_amd
            procedure, private :: destroy_mgm
            !! ACCESSORS 
            procedure          :: getMax_nusgs
            procedure          :: getMax_kapsgs
    end type 

contains

!!!#include "sgs_models/init_destroy_sgs_cgrid.F90"
#include "sgs_models/smagorinsky.F90"
#include "sgs_models/sigma.F90"
#include "sgs_models/AMD.F90"
#include "sgs_models/MGM.F90"
#include "sgs_models/eddyViscosity.F90"
#include "sgs_models/accessors.F90"


subroutine init(this, decomp, Cp, Pr, dx, dy, dz, inputfile)
  class(sgs_cgrid), intent(inout) :: this
  class(decomp_info), intent(in), target    :: decomp
  real(rkind) , intent(in) :: Cp, Pr, dx, dy, dz
  character(len=clen), intent(in) :: inputfile

  integer :: SGSmodelID = 0, DynProcFreq = 1, DynamicProcedureType = 0, ierr
  real(rkind) :: Csgs = 1.0_rkind, ncWall = 1.0_rkind, PrSGS = 1.0_rkind
  logical :: useWallDamping = .false., useSGSDynamicRestartFile = .false. 
  logical :: DomainAveraged_DynProc = .false., useDynamicProcedureScalar = .false.
  logical :: isEddyViscosityModel = .false., isEddyDiffmodel = .true., isTurbPrandtlconst = .true.
  character(len=clen) :: SGSDynamicRestartFile

   namelist /SGS_MODEL/ DynamicProcedureType, SGSmodelID,        &
                  useWallDamping, ncWall, Csgs, DynProcFreq,     &
                  useSGSDynamicRestartFile,                      &
                  DomainAveraged_DynProc, SGSDynamicRestartFile, &
                  useDynamicProcedureScalar, PrSGS,              &
                  isEddyViscosityModel, isEddyDiffmodel, isTurbPrandtlconst

   open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=123, NML=SGS_MODEL)
   close(123)

   this%decomp => decomp
   this%Cp = Cp
   this%Pr = Pr

   this%nxL = this%decomp%ysz(1)
   this%nyL = this%decomp%ysz(2)
   this%nzL = this%decomp%ysz(3)

   this%dx = dx;    this%dy = dy;    this%dz = dz

   if (.not. this%isPeriodic) then
      this%deltaLES = (1.5d0*this%dx*1.5d0*this%dy*this%dz)**(1.d0/3.d0)
   else
      this%deltaLES =  (1.5d0*this%dx*1.5d0*this%dy*1.5d0*this%dz)**(1.d0/3.d0)
   end if 
    
   allocate(  this%S_ij(this%nxL, this%nyL, this%nzL, 6))
   allocate( this%nusgs(this%nxL, this%nyL, this%nzL)   )
   allocate(this%kapsgs(this%nxL, this%nyL, this%nzL)   )

   this%SGSmodelID = SGSmodelID;             this%Csgs = Csgs;
   this%useWallDamping = useWallDamping;     this%ncWall = ncWall
   this%PrSGS = PrSGS

   this%DynProcFreq = DynProcFreq;                               this%useSGSDynamicRestartFile = useSGSDynamicRestartFile
   this%DomainAveraged_DynProc = DomainAveraged_DynProc;         this%SGSDynamicRestartFile = SGSDynamicRestartFile
   this%useDynamicProcedureScalar = useDynamicProcedureScalar;   this%DynamicProcedureType = DynamicProcedureType
   this%isEddyViscosityModel = isEddyViscosityModel;             this%isEddyDiffModel = isEddyDiffModel
   this%isTurbPrandtlconst = isTurbPrandtlconst

   select case (this%SGSmodelID)
   case (0)
      call this%init_smagorinsky()
   case (1)
      call this%init_sigma()
   case (2)
      call this%init_AMD()
   case (3)
      call this%init_mgm()
   case default
      call GracefulExit("Incorrect choice for SGS model ID.", 213)
   end select

    !!!! Safeguards against wrong inputs !!!!!
    if(this%isEddyDiffModel) then

      if( (.not. this%isEddyViscosityModel) ) then
          call GracefulExit("It can be an EddyDiff Model only if it is an EddyViscosityModel", 213)
      endif

      if( (.not. this%isTurbPrandtlconst) .and. (this%SGSmodelID .ne. 2) ) then
          call GracefulExit("To use EddyDiffModel with non-const TurbPrandtl, SGSModel must be AMD", 213)
      endif

    else

      if( this%SGSmodelID .ne. 3 ) then
          call GracefulExit("To use non-EddyDiffModel, SGSModel must be MGM", 213)
      endif

    endif
end subroutine 

subroutine destroy(this)
   class(sgs_cgrid), intent(inout) :: this

   select case (this%SGSmodelID)
   case (0)
      call this%destroy_smagorinsky()
   case (1)
      call this%destroy_sigma()
   case (2)
      call this%destroy_AMD()
   case (3)
      call this%destroy_mgm()
   case default
      call GracefulExit("Incorrect choice for SGS model ID.", 213)
   end select

   if(allocated(this%kapsgs)) deallocate(this%kapsgs)
   if(allocated(this%nusgs) ) deallocate(this%nusgs )
   if(allocated(this%S_ij)  ) deallocate(this%S_ij  )
   nullify(this%decomp)

end subroutine 

subroutine getTauSGS(this, duidxj, rho, tausgs)
   use constants, only : third, twothird 
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL, this%nyL, this%nzL, 9), intent(in)  :: duidxj
   real(rkind), dimension(this%nxL, this%nyL, this%nzL   ), intent(in)  :: rho
   real(rkind), dimension(this%nxL, this%nyL, this%nzL, 6), intent(out) :: tausgs

   real(rkind) :: S, CI = 0.003_rkind, deltaLES
   real(rkind), dimension(this%nxL, this%nyL, this%nzL)  :: modS_sq, q, Sii
   integer :: i,j,k

   if (this%isEddyViscosityModel) then

      ! Step 0: Compute Sij
      call this%get_Sij_from_duidxj(duidxj) 
      
      ! Step 1: Get mag-square and trace of S_ij
      call this%get_modS_Sii(modS_sq, Sii)

      ! Step 2: Get TKE (q)
      q = twothird * CI * (this%deltaLES**2) * rho * modS_sq

      ! Step 4: Get nusgs
      call this%get_SGS_kernel(duidxj, modS_sq)

      !!! Step 5: Dynamic Procedure ?

      ! Step 6: Multiply by model constant
      call this%multiply_by_model_constant()
     
      !! Subsumed in Step 1 
      !!do k = 1, this%nzL
      !!   do j = 1, this%nyL
      !!      do i = 1, this%nxL
      !!         S = this%S_ij(i,j,k,1)*this%S_ij(i,j,k,1) ! S11*S11
      !!         S = S + 2.d0*(this%S_ij(i,j,k,2)*this%S_ij(i,j,k,2)) ! S12*S12 + S21*S21
      !!         S = S + 2.d0*(this%S_ij(i,j,k,3)*this%S_ij(i,j,k,3)) ! S13*S13 + S31*S31
      !!         S = S + (this%S_ij(i,j,k,4)*this%S_ij(i,j,k,4)) ! S22*S22
      !!         S = S + 2.d0*(this%S_ij(i,j,k,5)*this%S_ij(i,j,k,5)) ! S23*S23 + S32*S32
      !!         S = S + (this%S_ij(i,j,k,6)*this%S_ij(i,j,k,6)) ! S33*S33
      !!      ! Now do modS = sqrt(2* S_ij*S_ij)
      !!         S = 2.d0*S
      !!         modS(i,j,k) = sqrt(S)
      !!      end do 
      !!   end do 
      !!end do

      ! Step 7: Get tau_sgs    
      tausgs(:,:,:,1) = -two * rho * this%nusgs * (this%S_ij(:,:,:,1)-Sii) + q
      tausgs(:,:,:,2) = -two * rho * this%nusgs * this%S_ij(:,:,:,2)
      tausgs(:,:,:,3) = -two * rho * this%nusgs * this%S_ij(:,:,:,3)
      tausgs(:,:,:,4) = -two * rho * this%nusgs * (this%S_ij(:,:,:,4)-Sii) + q
      tausgs(:,:,:,5) = -two * rho * this%nusgs * this%S_ij(:,:,:,5)
      tausgs(:,:,:,6) = -two * rho * this%nusgs * (this%S_ij(:,:,:,6)-Sii) + q
  
   else
     !! call MGM model from here
     call this%get_Sij_from_duidxj(duidxj)
     call this%get_tausgs_mgm(rho, duidxj, tausgs)
    
     
   end if

end subroutine

subroutine getQjSGS(this, duidxj, rho, gradT, Qjsgs)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL, this%nyL, this%nzL,9), intent(in)  :: duidxj
   real(rkind), dimension(this%nxL, this%nyL, this%nzL  ), intent(in)  :: rho
   real(rkind), dimension(this%nxL, this%nyL, this%nzL,3), intent(in)  :: gradT
   real(rkind), dimension(this%nxL, this%nyL, this%nzL,3), intent(out) :: Qjsgs

   integer :: k

   if (this%isEddyDiffModel) then
      if (this%isTurbPrandtlconst) then
          this%kapsgs = this%nusgs/this%PrSGS
      else
          call this%get_amd_Dkernel(duidxj, gradT)
      endif
      
      do k = 1, 3
         Qjsgs(:,:,:,k) = - rho * this%Cp * this%kapsgs * gradT(:,:,:,k)
      end do
      
   else
      !! MGM model
      call this%get_Sij_from_duidxj(duidxj) 
      call this%get_Qjsgs_mgm(rho, duidxj, gradT, Qjsgs)
   end if
end subroutine

end module 
