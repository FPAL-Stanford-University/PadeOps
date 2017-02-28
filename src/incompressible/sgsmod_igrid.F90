module sgsmod_igrid
    use kind_parameters, only: rkind, clen
    use constants, only: imi, pi, zero,one,two,three,half, four,eight, nine, six, kappa 
    use decomp_2d
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral  
    use mpi 
    use reductions, only: p_maxval, p_sum, p_minval
    use numerics, only: useCompactFD 
    use StaggOpsMod, only: staggOps  
    use gaussianstuff, only: gaussian
    use lstsqstuff, only: lstsq
    use PadeDerOps, only: Pade6stagg
    implicit none

    private
    public :: sgs_igrid

    complex(rkind), parameter :: zeroC = zero + imi*zero
    logical :: useVerticalTfilter = .false. 
    real(rkind), parameter :: beta_h = 7.8_rkind, beta_m = 4.8_rkind

    type :: sgs_igrid
        private
        type(decomp_info), pointer :: gpC, gpE
        type(spectral), pointer :: spectC, spectE
        type(decomp_info), pointer :: sp_gpC, sp_gpE
        integer :: mid, DynamicProcedureType, WallModel, DynProcFreq
        real(rkind), dimension(:), allocatable :: cmodelC, cmodelE
        real(rkind) :: cmodel_global, cmodel_global_x, cmodel_global_y, cmodel_global_z
        real(rkind), dimension(:,:,:), allocatable :: nu_sgs_C, nu_sgs_E
        logical :: isEddyViscosityModel = .false. 
        real(rkind), dimension(:,:,:), allocatable :: tau_11, tau_12, tau_13, tau_22, tau_23, tau_33
        real(rkind), dimension(:,:,:,:), allocatable :: S_ij_C, S_ij_E
        real(rkind), dimension(:,:,:,:), pointer :: rbuffxC, rbuffzC, rbuffyC, rbuffxE, rbuffyE, rbuffzE
        complex(rkind), dimension(:,:,:,:), pointer :: cbuffyC, cbuffzC, cbuffyE, cbuffzE
        type(Pade6stagg), pointer :: PadeDer

        ! Wall model
        real(rkind), dimension(:,:,:,:), allocatable :: tauijWM
        complex(rkind),dimension(:,:,:,:), allocatable :: tauijWMhat_inZ, tauijWMhat_inY
        real(rkind), dimension(:,:,:), allocatable :: filteredSpeedSq
        complex(rkind), dimension(:,:,:), allocatable :: Tfilhat, Tfilhatz1, Tfilhatz2
        logical :: useWallModel
        integer :: botBC_temp = 1
        real(rkind) :: ustar, InvObLength, umn, vmn, uspmn, Tmn, wTh_surf
        real(rkind) :: dz, z0, meanfact, ThetaRef, Fr, WallMfactor, Re
        real(rkind), pointer :: Tsurf

        ! for dynamic procedures - all are at edges
        type(gaussian) :: gaussianTestFilterZ
        real(rkind), dimension(:,:,:,:), allocatable :: Mij, Lij, Sij_Filt, alphaij_Filt, tauijWM_Filt
        real(rkind), dimension(:,:,:,:), allocatable :: fi_Filt, ui_Filt, fiE
        real(rkind), dimension(:,:,:),   allocatable :: Dsgs_Filt, buff1, buff2
        real(rkind), dimension(:,:,:,:), pointer     :: fiC
        real(rkind), dimension(:,:,:),   pointer     :: Dsgs
        logical :: isInviscid, isStratified, useDynamicProcedure, useVerticalTfilter = .false. 
        real(rkind) :: invRe, deltaRat
        integer :: mstep
        contains 
            !! ALL INIT PROCEDURES
            procedure          :: init
            procedure, private :: init_smagorinsky 
            procedure, private :: init_sigma
            procedure, private :: allocateMemory_EddyViscosity

            !! ALL WALL MODEL PROCEDURE
            procedure, private :: initWallModel
            procedure, private :: destroyWallModel
            procedure, private :: getfilteredSpeedSqAtWall
            procedure, private :: computeWallStress
            procedure, private :: compute_and_bcast_surface_Mn
            procedure, private :: getSurfaceQuantities
         
            !! ALL DYNAMIC PROCEDURE SUBROUTINES
            procedure, private :: allocateMemory_DynamicProcedure
            procedure, private :: destroyMemory_DynamicProcedure
            procedure, private :: applyDynamicProcedure 
            procedure, private :: TestFilter_Cmplx_to_Real
            procedure, private :: TestFilter_Real_to_Real
            procedure, private :: TestFilter_Real_to_Real_ip
            procedure, private :: planarAverageAndInterpolateToCells 
            procedure, private :: interp_bForce_CellToEdge
            procedure, private :: DoStandardDynamicProcedure
            procedure, private :: DoGlobalDynamicProcedure
            procedure, private :: readSGSDynamicRestart
            procedure, private :: dumpSGSDynamicRestart

            !! ALL GET_TAU PROCEDURES
            procedure          :: getTauSGS
            procedure          :: getRHS_SGS
            procedure, private :: get_SGS_kernel
            procedure, private :: multiply_by_model_constant 
            
            !! ALL DESTROY PROCEDURES
            procedure          :: destroy
            procedure, private :: destroy_smagorinsky 
            procedure, private :: destroy_sigma
            procedure, private :: destroyMemory_EddyViscosity
    end type 

contains
#include "sgs_models/init_destroy_sgs_igrid.F90"
#include "sgs_models/smagorinsky.F90"
#include "sgs_models/sigma.F90"
#include "sgs_models/eddyViscosity.F90"
#include "sgs_models/dynamicProcedure_sgs_igrid.F90"
#include "sgs_models/standardDynamicProcedure.F90"
#include "sgs_models/globalDynamicProcedure.F90"
#include "sgs_models/wallmodel.F90"



subroutine getRHS_SGS(this, urhs, vrhs, wrhs, duidxjC, duidxjE, duidxjEhat, uhatE, vhatE, whatE, uhatC, vhatC, ThatC, uC, vC, uE, vE, wE, newTimeStep)
   class(sgs_igrid), intent(inout), target :: this
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),9), intent(in) :: duidxjC
   real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),9), intent(in) :: duidxjE
   complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3),9), intent(in) :: duidxjEhat
   complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(in) :: uhatE, vhatE, whatE
   complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(in) :: uhatC, vhatC, ThatC
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: uC, vC
   real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)), intent(in) :: uE, vE, wE
   complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout) :: urhs, vrhs
   complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: wrhs
   logical, intent(in) :: newTimeStep

   complex(rkind), dimension(:,:,:), pointer :: cbuffy1, cbuffy2, cbuffy3, cbuffz1, cbuffz2
   
   call this%getTauSGS(duidxjC, duidxjE, duidxjEhat, uhatE, vhatE, whatE, uhatC, vhatC, ThatC, uC, vC, uE, vE, wE, newTimeStep)

   cbuffy1 => this%cbuffyC(:,:,:,1); cbuffy2 => this%cbuffyE(:,:,:,1); 
   cbuffz1 => this%cbuffzC(:,:,:,1); cbuffz2 => this%cbuffzE(:,:,:,1) 
   cbuffy3 => this%cbuffyC(:,:,:,2); 

   ! ddx(tau11)
   call this%spectC%fft(this%tau_11, cbuffy1)
   call this%spectC%mtimes_ik1_ip(cbuffy1)
   urhs = urhs - cbuffy1

   ! ddy(tau22)
   call this%spectC%fft(this%tau_22, cbuffy1)
   call this%spectC%mtimes_ik2_ip(cbuffy1)
   vrhs = vrhs - cbuffy1

   ! ddz(tau33)
   call this%spectC%fft(this%tau_33, cbuffy1)
   call transpose_y_to_z(cbuffy1, cbuffz1, this%sp_gpC)
   call this%PadeDer%ddz_C2E(cbuffz1, cbuffz2, 0, 0)
   call transpose_z_to_y(cbuffz2, cbuffy2, this%sp_gpE)
   wrhs = wrhs - cbuffy2

   ! ddy(tau12) for urhs, and ddx(tau12) for vrhs
   call this%spectC%fft(this%tau_12, cbuffy1)
   call this%spectC%mtimes_ik1_oop(cbuffy1, cbuffy3)
   vrhs = vrhs - cbuffy3
   call this%spectC%mtimes_ik2_ip(cbuffy1)
   urhs = urhs - cbuffy1

   ! ddz(tau13) for urhs, ddx(tau13) for wrhs
   call this%spectE%fft(this%tau_13, cbuffy2)
   if (this%useWallModel) cbuffy2 = cbuffy2 + this%tauijWMhat_inY(:,:,:,1)
   call transpose_y_to_z(cbuffy2, cbuffz2, this%sp_gpE)
   call this%PadeDer%ddz_E2C(cbuffz2, cbuffz1, 0, 0)
   call transpose_z_to_y(cbuffz1, cbuffy1, this%sp_gpC)
   urhs = urhs - cbuffy1
   call this%spectE%mtimes_ik1_ip(cbuffy2)
   wrhs = wrhs - cbuffy2

   ! ddz(tau23) for vrhs, ddy(tau23) for wrhs
   call this%spectE%fft(this%tau_23, cbuffy2)
   if (this%useWallModel) cbuffy2 = cbuffy2 + this%tauijWMhat_inY(:,:,:,2)
   call transpose_y_to_z(cbuffy2, cbuffz2, this%sp_gpE)
   call this%PadeDer%ddz_E2C(cbuffz2, cbuffz1, 0, 0)
   call transpose_z_to_y(cbuffz1, cbuffy1, this%sp_gpC)
   vrhs = vrhs - cbuffy1
   call this%spectE%mtimes_ik2_ip(cbuffy2)
   wrhs = wrhs - cbuffy2

end subroutine

subroutine getTauSGS(this, duidxjC, duidxjE, duidxjEhat, uhatE, vhatE, whatE, uhatC, vhatC, ThatC, uC, vC, uE, vE, wE, newTimeStep)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),9), intent(in) :: duidxjC
   real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),9), intent(in) :: duidxjE
   complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3),9), intent(in) :: duidxjEhat
   complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(in) :: uhatE, vhatE, whatE
   complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(in) :: uhatC, vhatC, ThatC
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: uC, vC
   real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)), intent(in) :: uE, vE, wE
   logical, intent(in) :: newTimeStep


   if(newTimeStep) this%mstep = this%mstep + 1

   if (this%useWallModel) call this%computeWallStress( uC, vC, uhatC, vhatC, ThatC) 

   if (this%isEddyViscosityModel) then

      ! Step 0: Compute Sij
      call get_Sij_from_duidxj(duidxjC, this%S_ij_C, this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)) 
      call get_Sij_from_duidxj(duidxjE, this%S_ij_E, this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)) 
      
      ! Step 1: Get nuSGS
      call this%get_SGS_kernel(duidxjC, duidxjE)

      ! Step 2: Dynamic Procedure ?
      if(newTimeStep .and. this%useDynamicProcedure) then
          if (mod(this%mstep, this%DynProcFreq ==0)) call this%applyDynamicProcedure(uE, vE, wE, uhatE, vhatE, whatE, duidxjE, duidxjEhat)
      endif

      ! Step 3: Multiply by model constant
      call this%multiply_by_model_constant()

      ! Step 2: Get tau_sgs
      this%tau_11 = -two*this%nu_sgs_C*this%S_ij_C(:,:,:,1)
      this%tau_12 = -two*this%nu_sgs_C*this%S_ij_C(:,:,:,2)
      this%tau_13 = -two*this%nu_sgs_E*this%S_ij_E(:,:,:,3)
      this%tau_22 = -two*this%nu_sgs_C*this%S_ij_C(:,:,:,4)
      this%tau_23 = -two*this%nu_sgs_E*this%S_ij_E(:,:,:,5)
      this%tau_33 = -two*this%nu_sgs_C*this%S_ij_C(:,:,:,6)

   end if
   
end subroutine
    

end module 
