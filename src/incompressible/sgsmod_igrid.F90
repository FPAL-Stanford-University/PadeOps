module sgsmod_igrid
    use kind_parameters, only: rkind, clen
    use constants, only: imi, pi, zero,one,two,three,half, four,eight, nine, six, kappa 
    use decomp_2d
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral  
    use mpi 
    use cd06staggstuff, only: cd06stagg
    use reductions, only: p_maxval, p_sum, p_minval
    use numerics, only: useCompactFD 
    use StaggOpsMod, only: staggOps  
    use gaussianstuff, only: gaussian
    use lstsqstuff, only: lstsq
    implicit none

    private
    public :: sgs_igrid

    complex(rkind), parameter :: zeroC = zero + imi*zero
    logical :: useVerticalTfilter = .false. 

    type :: sgs_igrid
        private
        type(decomp_info), pointer :: gpC, gpE
        type(spectral), pointer :: spectC, spectE
        type(decomp_info), pointer :: sp_gpC, sp_gpE
        integer :: mid, DynamicProcedureType, WallModelType
        real(rkind), dimension(:), allocatable :: cmodelC, cmodelE
        real(rkind) :: cmodel_global, cmodel_global_x, cmodel_global_y, cmodel_global_z
        real(rkind), dimension(:,:,:), allocatable :: nu_sgs_C, nu_sgs_E
        logical :: isEddyViscosityModel = .false. 
        real(rkind), dimension(:,:,:), allocatable :: tau_11, tau_12, tau_13, tau_22, tau_23, tau_33
        real(rkind), dimension(:,:,:,:), allocatable :: S_ij_C, S_ij_E
        ! for dyanamic procedures - all are at edges
        real(rkind), dimension(:,:,:,:), allocatable :: Sij_Filt, alphaij_Filt, tauijWM_Filt, fi_Filt, ui_Filt
        real(rkind), dimension(:,:,:),   allocatable :: Dsgs_Filt, buff1, buff2
        real(rkind), dimension(:,:,:,:), pointer     :: fi
        real(rkind), dimension(:,:,:),   pointer     :: Dsgs
        logical :: isInviscid
        real(rkind) :: invRe
        contains 
            !! ALL INIT PROCEDURES
            procedure          :: init
            procedure, private :: init_smagorinsky 
            procedure, private :: init_sigma
            procedure, private :: allocateMemory_EddyViscosity
            procedure, private :: allocateMemory_DynamicProcedure

            !! ALL GET_TAU PROCEDURES
            procedure, private :: getTauSGS
            procedure, private :: get_SGS_kernel
            procedure, private :: multiply_by_model_constant 
            
            !! ALL DESTROY PROCEDURES
            procedure          :: destroy
            procedure, private :: destroy_smagorinsky 
            procedure, private :: destroy_sigma
            procedure, private :: destroyMemory_EddyViscosity
            procedure, private :: destroyMemory_DynamicProcedure
    end type 

contains
#include "sgs_models/init_destroy_sgs_igrid.F90"
#include "sgs_models/smagorinsky.F90"
#include "sgs_models/sigma.F90"
#include "sgs_models/eddyViscosity.F90"
#include "sgs_models/dynamicProcedure_sgs_igrid.F90"


subroutine getTauSGS(this, duidxjC, duidxjE)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: duidxjC
   real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)), intent(in) :: duidxjE

   if (this%isEddyViscosityModel) then

      ! Step 0: Compute Sij
      call get_Sij_from_duidxj(duidxjC, this%S_ij_C, this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)) 
      call get_Sij_from_duidxj(duidxjE, this%S_ij_E, this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)) 
      
      ! Step 1: Get nuSGS
      call this%get_SGS_kernel(duidxjC, duidxjE)

      ! Step 2: Dynamic Procedure ?
      

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
