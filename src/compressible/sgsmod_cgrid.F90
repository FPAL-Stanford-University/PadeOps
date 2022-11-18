module sgsmod_cgrid
    use kind_parameters, only: rkind, clen
    use constants,       only: imi, zero !,one,two,three,half, four,eight, nine, six, kappa, piby2 
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
        integer :: SGSmodelID, mid, DynamicProcedureType
        real(rkind) :: cmodel_global,csgs,turbPrandtl = 0.7_rkind
        logical :: isEddyViscosityModel = .false., isEddyDiffmodel = .true., isTurbPrandtlconst = .true.
        real(rkind), dimension(:,:,:,:), allocatable :: S_ij
        real(rkind), dimension(:,:,:),   allocatable :: nusgs, kapsgs
        logical ::  isPeriodic = .false.

        !! model constant values/properties
        real(rkind) :: camd_x, camd_y, camd_z, cmgm_x, cmgm_y, cmgm_z, c1_mgm , c2_mgm, PrCpfac
        logical :: useCglobal = .false. 
        contains 
            !! ALL INIT PROCEDURES
            procedure          :: init
            !! ALL SGS SOURCE/SINK PROCEDURES
            procedure          :: getQjSGS
            procedure          :: getTauSGS
            procedure, private :: get_SGS_kernel
            procedure, private :: multiply_by_model_constant 
            !! ALL DESTROY PROCEDURES
            procedure          :: destroy
            procedure, private :: destroy_smagorinsky 
            procedure, private :: destroy_sigma
            procedure, private :: destroy_amd
            procedure, private :: destroy_mgm
            !! ACCESSORS 
            procedure          :: getMax_nuSGS
            procedure          :: getMax_kapSGS
    end type 

contains

#include "sgs_models/init_destroy_sgs_cgrid.F90"
#include "sgs_models/smagorinsky.F90"
#include "sgs_models/sigma.F90"
#include "sgs_models/AMD.F90"
#include "sgs_models/MGM.F90"
#include "sgs_models/eddyViscosity.F90"
#include "sgs_models/accessors.F90"


subroutine init(this, decomp, Cp, Pr)
   class(sgs_cgrid), intent(inout) :: this
   class(decomp_info), intent(in), target    :: decomp
   real(rkind) , intent(in) :: Cp, Pr
   this%decomp => decomp
   this%Cp = Cp
   this%Pr = Pr
   allocate(  this%S_ij(this%decomp%ysz(1), this%decomp%ysz(2), this%decomp%ysz(3), 6))
   allocate( this%nusgs(this%decomp%ysz(1), this%decomp%ysz(2), this%decomp%ysz(3))   )
   allocate(this%kapsgs(this%decomp%ysz(1), this%decomp%ysz(2), this%decomp%ysz(3))   )

end subroutine 

subroutine destroy(this)
   class(sgs_cgrid), intent(inout) :: this

   if(allocated(this%kapsgs)) deallocate(this%kapsgs)
   if(allocated(this%nusgs) ) deallocate(this%nusgs )
   if(allocated(this%S_ij)  ) deallocate(this%S_ij  )
   nullify(this%decomp)

end subroutine 

subroutine getTauSGS(this, duidxj, rho, tausgs, nxL, nyL, nzL, dx, dy, dz) 
   use constants, only : third, twothird 
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), intent(in) :: nxL, nyL, nzL, dx, dy, dz
   real(rkind), dimension(nxL,nyL,nzL,9), intent(in)  :: duidxj
   real(rkind), dimension(nxL, nyL, nzL), intent(in)  :: rho
   real(rkind), dimension(nxL,nyL,nzL,6), intent(out) :: tausgs
   real(rkind) :: S, CI = 0.003_rkind, deltaLES
   real(rkind), dimension(nxL,nyL,nzL)  :: modS, Sii, q
   integer :: i,j,k

   if (this%isEddyViscosityModel) then

      ! Step 0: Compute Sij
      call get_Sij_from_duidxj(duidxj, this%S_ij, nxL, nyL, nzL) 
      
      ! Step 1: Get nuSGS
      call this%get_SGS_kernel(duidxj)

      !!! Step 2: Dynamic Procedure ?

      ! Step 3: Multiply by model constant
      call this%multiply_by_model_constant()
      
      do k = 1,nzL
         do j = 1,nyL
            do i = 1,nxL
               S = this%S_ij(i,j,k,1)*this%S_ij(i,j,k,1) ! S11*S11
               S = S + 2.d0*(this%S_ij(i,j,k,2)*this%S_ij(i,j,k,2)) ! S12*S12 + S21*S21
               S = S + 2.d0*(this%S_ij(i,j,k,3)*this%S_ij(i,j,k,3)) ! S13*S13 + S31*S31
               S = S + (this%S_ij(i,j,k,4)*this%S_ij(i,j,k,4)) ! S22*S22
               S = S + 2.d0*(this%S_ij(i,j,k,5)*this%S_ij(i,j,k,5)) ! S23*S23 + S32*S32
               S = S + (this%S_ij(i,j,k,6)*this%S_ij(i,j,k,6)) ! S33*S33
            ! Now do modS = sqrt(2* S_ij*S_ij)
               S = 2.d0*S
               modS(i,j,k) = sqrt(S)
            end do 
         end do 
      end do
     
      if (.not. this%isPeriodic) then
         deltaLES = (1.5d0*dx*1.5d0*dy*dz)**(1.d0/3.d0)
      else
         deltaLES =  (1.5d0*dx*1.5d0*dy*1.5d0*dz)**(1.d0/3.d0)
      end if 
    
      q   = twothird * CI * rho * (deltaLES**2) * (modS**2)
      Sii = (this%S_ij(:,:,:,1) + this%S_ij(:,:,:,4) + this%S_ij(:,:,:,6)) * third
      ! Step 2: Get tau_sgs    
      tausgs(:,:,:,1) = -two * rho * this%nusgs * (this%S_ij(:,:,:,1)-Sii) + q
      tausgs(:,:,:,2) = -two * rho * this%nusgs * this%S_ij(:,:,:,2) 
      tausgs(:,:,:,3) = -two * rho * this%nusgs * this%S_ij(:,:,:,3) 
      tausgs(:,:,:,4) = -two * rho * this%nusgs * (this%S_ij(:,:,:,4)-Sii) + q
      tausgs(:,:,:,5) = -two * rho * this%nusgs * this%S_ij(:,:,:,5) 
      tausgs(:,:,:,6) = -two * rho * this%nusgs * (this%S_ij(:,:,:,6)-Sii) + q
   else
     !! call MGM model from here
      call get_Sij_from_duidxj(duidxj, this%S_ij, nxL, nyL, nzL) 
      get_tausgs_mgm(this%cmgm_x, this%cmgm_y, this%cmgm_z, this%c1_mgm, rho, duidxj, this%S_ij, nxL, nyL, nzL, tausgs)      
   end if

end subroutine


subroutine getQjSGS(this,duidxj, rho, gradT, Qjsgs, nxL, nyL, nzL, dx, dy, dz)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), intent(in) :: nxL, nyL, nzL, dx, dy, dz
   real(rkind), dimension(nxL,nyL,nzL,9), intent(in)  :: duidxj
   real(rkind), dimension(nxL, nyL, nzL), intent(in)  :: rho
   real(rkind), dimension(nxL,nyL,nzL,3), intent(in)  :: gradT
   real(rkind), dimension(nxL,nyL,nzL,3), intent(out) :: Qjsgs
   real(rkind) :: k
!! not yet completed
   if (this%isEddyDiffModel) then
      if (this%IsPrandtlconstant) then
         this%kapsgs = this%nusgs/this%TurbPrandtl
      elseif () then 
         get_amd_Dkernel(this%kapsgs,this%camd_x,this%camd_y,this%camd_z, duidxj, gradT, nxL, nyL, nzL)
      else
         call message ()
      end if
      
      do k=1,3
         Qjsgs(:,:,:,k) = - rho * this%Cp * this%kapsgs * gradT(:,:,:,k)
      end do
      
   else   
      !! MGM model
      call get_Sij_from_duidxj(duidxj, this%S_ij, nxL, nyL, nzL) 
      get_Qjsgs_mgm(this%cmgm_x, this%cmgm_y, this%cmgm_z, this%c2_mgm, this%PrCpfac, rho, duidxj,gradT, this%S_ij, nxL, nyL, nzL,Qjsgs)
   end if
end subroutine

end module 
