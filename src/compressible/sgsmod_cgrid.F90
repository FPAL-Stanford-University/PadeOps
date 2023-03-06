module sgsmod_cgrid
    use kind_parameters, only: rkind, clen
    use constants,       only: imi, zero, two, eight, eps, third, one, two, three, half, four,eight, nine, six, kappa, pi, twothird
    use decomp_2d
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
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
        type(derivatives), pointer :: der
        type(decomp_info), pointer :: decomp
        type(filters), allocatable :: testfil
        integer :: SGSmodelID, mid, DynamicProcedureType, DynProcFreq, mstep = 0
        logical :: isEddyViscosityModel, isEddyDiffmodel, isTurbPrandtlconst
        real(rkind), dimension(:,:,:,:), allocatable :: S_ij
        real(rkind), dimension(:,:,:,:), allocatable :: ybuflocal
        real(rkind), dimension(:,:,:),   allocatable :: nusgs, kapsgs
        real(rkind), dimension(:),       allocatable :: cmodel_local, cmodel_local_Qjsgs, cmodel_local_tke
        real(rkind), dimension(:,:,:,:), pointer     :: xbuf, ybuf, zbuf, duiFildxj, SFil_ij, tausgsFil, gradTFil, QjsgsFil
        real(rkind), dimension(:,:,:  ), pointer     :: uFil, vFil, wFil, rhoFil, TFil, numer, denom
        real(rkind), dimension(:,:,:  ), pointer     :: duFildx, duFildy, duFildz, Lij, Mij
        real(rkind), dimension(:,:,:  ), pointer     :: dvFildx, dvFildy, dvFildz
        real(rkind), dimension(:,:,:  ), pointer     :: dwFildx, dwFildy, dwFildz
        real(rkind), dimension(:,:,:  ), pointer     :: dTFildx, dTFildy, dTFildz, nusgsFil, SiiFil, modSFil_sq, qtkeFil! , kapsgsFil
        real(rkind) :: cmodel_global, cmodel_global_Qjsgs, deltaRatioSq
        logical ::  isPeriodic = .false., periodicx = .true., periodicy = .true., periodicz = .true.
        integer, dimension(2) :: x_bc, y_bc, z_bc
        logical :: filter_in_x = .true., filter_in_y = .false., filter_in_z = .true., preComputed_SFil_duFil

        integer :: nxL, nyL, nzL
        real(rkind) :: dx, dy, dz, deltaLES

        real(rkind) :: Csgs, Ctke, ncWall, PrSGS, Cp, Pr
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
            procedure, private :: init_dynamic_procedure
            !! ALL SGS SOURCE/SINK PROCEDURES
            procedure          :: getQjSGS
            procedure          :: getTauSGS
            procedure          :: get_qtke
            procedure, private :: get_SGS_kernel
            procedure, private :: get_smagorinsky_kernel
            procedure, private :: get_sigma_kernel
            procedure, private :: get_amd_kernel
            procedure, private :: get_amd_Dkernel
            !!procedure, private :: get_tausgs_mgm
            !!procedure, private :: get_Qjsgs_mgm
            procedure, private :: get_mgm_kernel
            procedure, private :: get_Qjsgs_mgm_kernel
            procedure, private :: get_Qjsgs_eddy_kernel
            procedure, private :: multiply_by_model_coefficient_mgm
            procedure, private :: multiply_by_model_coefficient_mgm_Qjsgs
            procedure, private :: multiply_by_model_coefficient_eddy_Qjsgs
            procedure, private :: multiply_by_model_constant 
            procedure, private :: get_Sij_from_duidxj
            procedure, private :: get_modS_Sii
            procedure, private :: doLocalDynamicProcedure_mgm
            procedure, private :: doLocalDynamicProcedure_mgm_Qjsgs
            procedure, private :: doLocalDynamicProcedure_eddy
            procedure, private :: doLocalDynamicProcedure_eddy_Qjsgs
            procedure, private :: filter_xyz
            procedure, private :: gradient
            !! ALL DESTROY PROCEDURES
            procedure          :: destroy
            procedure, private :: destroy_smagorinsky 
            procedure, private :: destroy_sigma
            procedure, private :: destroy_amd
            procedure, private :: destroy_mgm
            procedure, private :: destroy_dynamic_procedure
            !! ACCESSORS 
            procedure          :: getMax_nusgs
            procedure          :: getMax_kapsgs
            procedure          :: get_LocalDynamicProcedure_Coeff
            procedure          :: get_LocalDynamicProcedure_Coeff_Qjsgs
            procedure          :: get_LocalDynamicProcedure_Coeff_tke
            procedure          :: get_Max_LocalDynamicProcedure_Coeff
            procedure          :: get_Max_LocalDynamicProcedure_Coeff_Qjsgs
            procedure          :: get_Max_LocalDynamicProcedure_Coeff_tke
            procedure          :: get_Min_LocalDynamicProcedure_Coeff
            procedure          :: get_Min_LocalDynamicProcedure_Coeff_Qjsgs
            procedure          :: get_Min_LocalDynamicProcedure_Coeff_tke
    end type 

contains

!!!#include "sgs_models/init_destroy_sgs_cgrid.F90"
#include "sgs_models/smagorinsky.F90"
#include "sgs_models/sigma.F90"
#include "sgs_models/AMD.F90"
#include "sgs_models/MGM.F90"
#include "sgs_models/eddyViscosity.F90"
#include "sgs_models/accessors.F90"
#include "sgs_models/dynamicProcedures.F90"


subroutine init(this, der, decomp, Cp, Pr, dx, dy, dz, inputfile, xbuf, ybuf, zbuf, periodicx, periodicy, periodicz, x_bc1, x_bcn, y_bc1, y_bcn, z_bc1, z_bcn)
  class(sgs_cgrid), intent(inout) :: this
  class(derivatives), intent(in), target :: der
  class(decomp_info), intent(in), target :: decomp
  real(rkind) , intent(in) :: Cp, Pr, dx, dy, dz
  character(len=clen), intent(in) :: inputfile
  real(rkind), dimension(:,:,:,:), intent(in), target :: xbuf, ybuf, zbuf
  logical, intent(in) ::  periodicx, periodicy, periodicz
  integer, intent(in) :: x_bc1, x_bcn
  integer, intent(in) :: y_bc1, y_bcn
  integer, intent(in) :: z_bc1, z_bcn

  integer :: SGSmodelID = 0, DynProcFreq = 1, DynamicProcedureType = 0, ierr
  real(rkind) :: Csgs = 1.0_rkind, Ctke = 0.003_rkind, ncWall = 1.0_rkind, PrSGS = 1.0_rkind, deltaRatio = 2.0_rkind
  logical :: useWallDamping = .false., useSGSDynamicRestartFile = .false. 
  logical :: DomainAveraged_DynProc = .false., useDynamicProcedureScalar = .false.
  logical :: isEddyViscosityModel = .false., isEddyDiffmodel = .true., isTurbPrandtlconst = .true.
  character(len=clen) :: SGSDynamicRestartFile
  character(len=clen) :: testfilter_x = "box2", testfilter_y = "box2", testfilter_z = "box2"
  logical :: filter_in_x = .true., filter_in_y = .false., filter_in_z = .true.

   namelist /SGS_MODEL/ DynamicProcedureType, SGSmodelID,        &
                  useWallDamping, ncWall, Csgs, Ctke, DynProcFreq,     &
                  useSGSDynamicRestartFile,                      &
                  DomainAveraged_DynProc, SGSDynamicRestartFile, &
                  useDynamicProcedureScalar, PrSGS,              &
                  filter_in_x, filter_in_y, filter_in_z,         &
                  testfilter_x, testfilter_y, testfilter_z, deltaRatio, &
                  isEddyViscosityModel, isEddyDiffmodel, isTurbPrandtlconst

   open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=123, NML=SGS_MODEL)
   close(123)

   this%der    => der   
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
    
   this%xbuf => xbuf
   this%ybuf => ybuf
   this%zbuf => zbuf

   this%periodicx = periodicx;  this%periodicy = periodicy;  this%periodicz = periodicz
   this%x_bc = [x_bc1, x_bcn]
   this%y_bc = [y_bc1, y_bcn]
   this%z_bc = [z_bc1, z_bcn]
   this%filter_in_x = filter_in_x;   this%filter_in_y = filter_in_y;   this%filter_in_z = filter_in_z

   allocate(  this%S_ij(this%nxL, this%nyL, this%nzL, 6))
   allocate( this%nusgs(this%nxL, this%nyL, this%nzL)   )
   allocate(this%kapsgs(this%nxL, this%nyL, this%nzL)   )

   this%SGSmodelID = SGSmodelID;             this%Csgs = Csgs;   this%Ctke = Ctke
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

   call this%init_dynamic_procedure(testfilter_x, testfilter_y, testfilter_z, deltaRatio)

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

    if(this%DynamicProcedureType==2) then
        call GracefulExit("Global-Dynamic Procedure not enabled as of now", 213)
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

   call this%destroy_dynamic_procedure()

   if(allocated(this%kapsgs)) deallocate(this%kapsgs)
   if(allocated(this%nusgs) ) deallocate(this%nusgs )
   if(allocated(this%S_ij)  ) deallocate(this%S_ij  )

   nullify(this%zbuf)
   nullify(this%ybuf)
   nullify(this%xbuf)
   nullify(this%der)
   nullify(this%decomp)

end subroutine 

subroutine getTauSGS(this, newTimeStep, duidxj, rho, u, v, w, tausgs)
   use constants, only : third, twothird 
   class(sgs_cgrid), intent(inout) :: this
   logical, intent(in) :: newTimeStep
   real(rkind), dimension(this%nxL, this%nyL, this%nzL, 9), intent(in)  :: duidxj
   real(rkind), dimension(this%nxL, this%nyL, this%nzL   ), intent(in)  :: rho, u, v, w
   real(rkind), dimension(this%nxL, this%nyL, this%nzL, 6), intent(out) :: tausgs

   !real(rkind) :: S, CI = 0.003_rkind, deltaLES
   real(rkind), dimension(this%nxL, this%nyL, this%nzL)  :: modS_sq, qtke, Sii
   integer :: i,j,k

   if (this%isEddyViscosityModel) then

      ! Step 0: Compute Sij
      call this%get_Sij_from_duidxj(duidxj, this%S_ij) 
      
      ! Step 1: Get mag-square and trace of S_ij
      call this%get_modS_Sii(this%S_ij,modS_sq, Sii)

      ! Step 2: Get TKE (q)
      call this%get_qtke(rho, modS_sq, qtke)
      !qtke = twothird * (this%deltaLES**2) * rho * modS_sq

      ! Step 4: Get nusgs
      call this%get_SGS_kernel(duidxj, this%S_ij, modS_sq, this%nusgs)

      ! Step 5: Dynamic Procedure 
      if( (newTimeStep) .and. (mod(this%mstep, this%DynProcFreq) == 0) ) then
        if(this%DynamicProcedureType==1) then
            call this%DoLocalDynamicProcedure_eddy(rho, u, v, w, this%nusgs, this%S_ij, Sii, modS_sq, qtke)
        elseif(this%DynamicProcedureType==2) then
            !! not yet implemented
        endif
      endif
      
      ! Step 6: Multiply by model constant     
      call this%multiply_by_model_constant(qtke)
     
      ! Step 7: Get tau_sgs    
      tausgs(:,:,:,1) = -two * rho * this%nusgs * (this%S_ij(:,:,:,1)-Sii) + qtke
      tausgs(:,:,:,2) = -two * rho * this%nusgs * this%S_ij(:,:,:,2)
      tausgs(:,:,:,3) = -two * rho * this%nusgs * this%S_ij(:,:,:,3)
      tausgs(:,:,:,4) = -two * rho * this%nusgs * (this%S_ij(:,:,:,4)-Sii) + qtke
      tausgs(:,:,:,5) = -two * rho * this%nusgs * this%S_ij(:,:,:,5)
      tausgs(:,:,:,6) = -two * rho * this%nusgs * (this%S_ij(:,:,:,6)-Sii) + qtke
  
   else
     !!!! call MGM model from here
     !!call this%get_Sij_from_duidxj(duidxj)
     !!call this%get_tausgs_mgm(rho, duidxj, tausgs)
   
     !!! new way of writing MGM !!!
     call this%get_Sij_from_duidxj(duidxj, this%S_ij)
     call this%get_mgm_kernel(rho, duidxj, this%S_ij, tausgs)
 
     !! Dynam Procedure
     if( (newTimeStep) .and. (mod(this%mstep, this%DynProcFreq) == 0) ) then
       if(this%DynamicProcedureType==1) then
           call this%DoLocalDynamicProcedure_mgm(rho, u, v, w, tausgs)
       elseif(this%DynamicProcedureType==2) then
           !! not yet implemented
       endif
     endif

     call this%multiply_by_model_coefficient_MGM(tausgs)
 
   end if

   if(newTimeStep) this%mstep = this%mstep + 1

end subroutine

subroutine getQjSGS(this, newTimeStep, duidxj, rho, u, v, w, T, gradT, Qjsgs)
   class(sgs_cgrid), intent(inout) :: this
   logical, intent(in) :: newTimeStep
   real(rkind), dimension(this%nxL, this%nyL, this%nzL,9), intent(in)  :: duidxj
   real(rkind), dimension(this%nxL, this%nyL, this%nzL  ), intent(in)  :: rho, u, v, w, T
   real(rkind), dimension(this%nxL, this%nyL, this%nzL,3), intent(in)  :: gradT
   real(rkind), dimension(this%nxL, this%nyL, this%nzL,3), intent(out) :: Qjsgs

   integer :: k

   if (this%isEddyDiffModel) then
     ! if (this%isTurbPrandtlconst) then
     !     this%kapsgs = this%nusgs ! /this%PrSGS
     ! else
     !     call this%get_amd_Dkernel(duidxj, gradT, this%kapsgs)
     ! endif
     ! 
     ! do k = 1, 3
     !    Qjsgs(:,:,:,k) = - rho * this%Cp * this%kapsgs * gradT(:,:,:,k)
     ! end do
     !!! Above copied here
     call this%get_Qjsgs_eddy_kernel(rho, this%nusgs, duidxj, gradT, Qjsgs)

     !! Dynam Procedure
     if( (newTimeStep) .and. (mod(this%mstep, this%DynProcFreq) == 0) ) then
       if(this%DynamicProcedureType==1) then
           call this%DoLocalDynamicProcedure_eddy_Qjsgs(rho, u, v, w, T, Qjsgs)
       elseif(this%DynamicProcedureType==2) then
           !! not yet implemented
       endif
     endif

     call this%multiply_by_model_coefficient_eddy_Qjsgs(Qjsgs) 
      
   else
     !!!! MGM model
     !!call this%get_Sij_from_duidxj(duidxj) 
     !!call this%get_Qjsgs_mgm(rho, duidxj, gradT, Qjsgs)
   
     !!! new way of writing MGM !!!
     call this%get_Sij_from_duidxj(duidxj, this%S_ij)
     call this%get_Qjsgs_mgm_kernel(rho, duidxj, this%S_ij, gradT, Qjsgs)
 
     !! Dynam Procedure
     if( (newTimeStep) .and. (mod(this%mstep, this%DynProcFreq) == 0) ) then
       if(this%DynamicProcedureType==1) then
           call this%DoLocalDynamicProcedure_mgm_Qjsgs(rho, u, v, w, T, Qjsgs)
       elseif(this%DynamicProcedureType==2) then
           !! not yet implemented
       endif
     endif

     call this%multiply_by_model_coefficient_mgm_Qjsgs(Qjsgs)
 
   end if
end subroutine

end module 
