module sgsmod
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
    public :: sgs

    real(rkind) :: c_sigma = 1.3_rkind
    real(rkind) :: c_smag = 0.17_rkind
    real(rkind) :: c_mgm = 1.0_rkind
    integer, parameter     :: mgm_option = 2
    real(rkind), parameter :: deltaRatio = 2.0_rkind
    complex(rkind), parameter :: zeroC = zero + imi*zero
    logical :: useVerticalTfilter = .false. 
    integer :: applyDynEvery  = 10
    real(rkind), parameter :: cx_amd = 1.d0/sqrt(12.d0), cy_amd = 1.d0/sqrt(12.d0), cz_amd = 1.d0/sqrt(3.d0)

    real(rkind), parameter :: bm = 4.8_rkind, bh = 7.8_rkind
    type :: sgs
        private
        type(spectral), pointer :: spectC, spectE 
        real(rkind), allocatable, dimension(:,:,:,:) :: rbuff, rbuffE, SIGMAbuffs, Lij, Mij, MGMbuffs, MGMbuffsE
        real(rkind), pointer, dimension(:,:,:) :: cSMAG_WALL, nuSGS, nuSGSfil,nuSGSE, nuSCA, nuSCAE
        real(rkind) :: deltaFilter, mconst, deltaTFilter 
        type(decomp_info), pointer :: sp_gp, gpC
        type(decomp_info), pointer :: sp_gpE, gpE

        complex(rkind), allocatable, dimension(:,:,:,:) :: cbuff
        complex(rkind), allocatable, dimension(:,:,:) :: ctmpCz, ctmpEz, ctmpEy, ctmpCz2
        complex(rkind), pointer, dimension(:,:,:) :: nuSGShat
        
        real(rkind), dimension(:,:,:), allocatable :: uf, vf, wCf, rtmpY, rtmpZ, rtmpZ2, rtmpZE2, rtmpYE, rtmpZE, dTdzC_diff
        real(rkind), dimension(:,:,:), allocatable :: SSI_xbuff1, SSI_xbuff2, q1, q2, q3, AMD_NUM, AMD_DEN

        type(cd06stagg), allocatable :: derZ_SS, derTAU33
        type(staggOps), allocatable :: Ops2ndOrder, OpsNU

        logical :: useWallFunction = .false. 
        logical :: useDynamicProcedure = .false.
        logical :: useClipping = .false., CompStokesP = .false. 
        logical :: eddyViscModel = .true., isStratified = .false. 
        
        real(rkind) :: meanFact, Pr, Theta0, Fr

        logical :: useWallModel = .false.  
        real(rkind) :: z0, dz, dxsq, dysq, dzsq, cx, cy, cz 
        integer :: nz, ntimeAvgQs = 3

        type(gaussian) :: Gfilz 
        type(lstsq) :: Tfilz 
        integer :: SGSmodel 

        integer :: WallModel = 1
        real(rkind) :: WallMfactor
        integer :: mstep = 0
        contains 
            procedure :: init
            procedure :: destroy
            procedure :: setStratificationConstants
            procedure, private :: DynamicProcedure
            procedure, private :: planarAverage 
            procedure, private :: planarAverage_oop 
            procedure, private :: get_SMAG_Op
            procedure, private :: get_AMD_Op
            procedure, private :: get_AMD_Op_Strat
            procedure, private :: get_ShearImpSMAG_Op
            procedure, private :: get_SIGMA_Op
            procedure, private :: get_MGM_Op
            procedure, private :: testFilter_ip
            procedure, private :: testFilter_oop
            procedure, private :: testFilter_oop_C2R
            procedure, private :: filtCoeff
            procedure, private :: getLocalWallModel
            procedure :: getRHS_SGS
            procedure :: getRHS_SGS_WallM
            procedure :: getRHS_SGS_Scalar_WallM
            procedure :: link_pointers
    end type 

contains

#include "sgs_models/initialize.F90"
#include "sgs_models/sigma_model_get_nuSGS.F90"
#include "sgs_models/amd_model_get_nuSGS.F90"
#include "sgs_models/mgm_model.F90"

    subroutine setStratificationConstants(this, Fr, Theta0)
        class(sgs), intent(inout) :: this
        real(rkind), intent(in) :: Fr, Theta0
    
        this%Fr = Fr; this%theta0 = theta0
    end subroutine

    subroutine link_pointers(this,nuSGS,c_SGS, tauSGS_ij, tau13, tau23, q1, q2, q3)
        class(sgs), intent(in), target :: this
        real(rkind), dimension(:,:,:)  , pointer, intent(inout) :: c_SGS, nuSGS
        real(rkind), dimension(:,:,:)  , pointer, intent(inout), optional :: tau13, tau23
        real(rkind), dimension(:,:,:,:), pointer, intent(inout) :: tauSGS_ij
        real(rkind), dimension(:,:,:)  , pointer, intent(inout) :: q1, q2, q3

        if (allocated(this%rbuff)) then
            if (this%useDynamicProcedure) then
                c_SGS => this%Lij(:,:,:,1)
            else
                c_SGS => this%rbuff(:,:,:,8) !Locations that are certainly zero
            end if 
            nuSGS => this%rbuff(:,:,:,7)
            tauSGS_ij => this%rbuff(:,:,:,1:6)
            if ((present(tau13)) .and. (present(tau23))) then
                if (this%useWallModel) then
                    tau13 => this%rbuffE(:,:,:,1)
                    tau23 => this%rbuffE(:,:,:,2)
                else
                    tau13 => this%rbuff(:,:,:,3)
                    tau23 => this%rbuff(:,:,:,5)
                end if 
            end if 
            if (this%isStratified) then
                q1 => this%q1; q2 => this%q2; q3 => this%q3
            end if 
        else
            call gracefulExit("You have called SGS%LINK_POINTERS before &
                & initializing SGS",324)
        end if 


    end subroutine


    subroutine get_SMAG_Op(this, nuSGS, S11, S22, S33, S12, S13, S23)
        class(sgs), intent(inout) :: this
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: S11, S22, S33
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: S12, S13, S23
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(out) :: nuSGS

        nuSGS = S12*S12 + S13*S13 + S23*S23
        nuSGS = two*nuSGS
        nuSGS = nuSGS + S11*S11 + S22*S22 + S33*S33
        nuSGS = two*nuSGS
        nuSGS = sqrt(nuSGS)

    end subroutine

    subroutine get_ShearImpSMAG_Op(this, nuSGS, S11, S22, S33, S12, S13, S23)
        class(sgs), intent(inout) :: this
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: S11, S22, S33
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: S12, S13, S23
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(out) :: nuSGS

        nuSGS = S12*S12 + S13*S13 + S23*S23
        nuSGS = two*nuSGS
        nuSGS = nuSGS + S11*S11 + S22*S22 + S33*S33
        nuSGS = two*nuSGS
        nuSGS = sqrt(nuSGS)

        call this%planarAverage_oop(S12,this%SSI_xbuff1)
        this%SSI_xbuff2 = this%SSI_xbuff1*this%SSI_xbuff1
        call this%planarAverage_oop(S13,this%SSI_xbuff1)
        this%SSI_xbuff2 = this%SSI_xbuff2 + this%SSI_xbuff1*this%SSI_xbuff1
        call this%planarAverage_oop(S23,this%SSI_xbuff1)
        this%SSI_xbuff2 = this%SSI_xbuff2 + this%SSI_xbuff1*this%SSI_xbuff1
        this%SSI_xbuff2 = two*this%SSI_xbuff2

        call this%planarAverage_oop(S11,this%SSI_xbuff1)
        this%SSI_xbuff2 = this%SSI_xbuff2 + this%SSI_xbuff1*this%SSI_xbuff1
        call this%planarAverage_oop(S22,this%SSI_xbuff1)
        this%SSI_xbuff2 = this%SSI_xbuff2 + this%SSI_xbuff1*this%SSI_xbuff1
        call this%planarAverage_oop(S33,this%SSI_xbuff1)
        this%SSI_xbuff2 = this%SSI_xbuff2 + this%SSI_xbuff1*this%SSI_xbuff1
        this%SSI_xbuff2 = two*this%SSI_xbuff2
        this%SSI_xbuff2 = sqrt(this%SSI_xbuff2)

        nuSGS = nuSGS - this%SSI_xbuff2

    end subroutine

    subroutine getRHS_SGS(this, duidxj, duidxjhat, urhs, vrhs, wrhs, uhat, vhat, wChat, u, v, wC, max_nuSGS)
        class(sgs), intent(inout), target :: this
        real(rkind)   , dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),9), intent(inout), target :: duidxj
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3),9), intent(inout) :: duidxjhat
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(inout) :: urhs, vrhs
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(in) :: uhat, vhat, wChat
        complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: wrhs
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(inout) :: u, v, wC

        real(rkind), dimension(:,:,:), pointer :: tau11, tau12, tau13, tau22, tau23, tau33
        real(rkind), dimension(:,:,:), pointer :: S11, S12, S13, S22, S23, S33
        complex(rkind), dimension(:,:,:), pointer :: tauhat, tauhat2    
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
        real(rkind), intent(out), optional :: max_nuSGS
        integer :: k

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3)
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6)
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9)

        ! STEP 0: Associate the pointers to buffers
        tau11 => this%rbuff(:,:,:,1); tau12 => this%rbuff(:,:,:,2); tau13 => this%rbuff(:,:,:,3)
        tau22 => this%rbuff(:,:,:,4); tau23 => this%rbuff(:,:,:,5); tau33 => this%rbuff(:,:,:,6)
        S11 => this%rbuff(:,:,:,1); S12 => this%rbuff(:,:,:,2); S13 => this%rbuff(:,:,:,3)
        S22 => this%rbuff(:,:,:,4); S23 => this%rbuff(:,:,:,5); S33 => this%rbuff(:,:,:,6)


        ! STEP 1: Compute S_ij 
        S11 = dudx; S22 = dvdy; S33 = dwdz   
        S12 = half*(dvdx + dudy); S13 = half*(dwdx + dudz); S23 = half*(dvdz + dwdy)
        
        ! STEP 2: Call the SGS model operator
        select case (this%SGSmodel)
        case (0)     
            call this%get_SMAG_Op(this%nuSGS, S11,S22,S33,S12,S13,S23)
        case (1)
            call this%get_SIGMA_Op(this%nuSGS, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
        case (2)
            call this%get_SMAG_Op(this%nuSGS, S11,S22,S33,S12,S13,S23)
        end select
        
        if (this%useDynamicProcedure) then
            call this%DynamicProcedure(u,v,wC,uhat, vhat, wChat, duidxjhat) 
        else
            if (this%useWallFunction) then
                this%nuSGS = this%cSMAG_WALL*this%nuSGS
            else
                this%nuSGS = this%mconst*this%nuSGS
            end if 
        end if 

        
        do k = 1,size(this%nuSGS,3)
            !print*, sum(this%nuSGS(:,:,k))/(size(this%nuSGS,1)*size(this%nuSGS,2))
        end do 

        tau11 = -two*this%nuSGS*S11; tau12 = -two*this%nuSGS*S12; tau13 = -two*this%nuSGS*S13
        tau22 = -two*this%nuSGS*S22; tau23 = -two*this%nuSGS*S23; tau33 = -two*this%nuSGS*S33

        tauhat => this%cbuff(:,:,:,1); tauhat2 => this%cbuff(:,:,:,2)       


        if (useCompactFD) then
#include "sgs_models/getRHS_common_CD.F90"
        else
#include "sgs_models/getRHS_common_FD.F90"
        end if
        
        if (present(max_nuSGS)) then
            max_nuSGS = p_maxval(maxval(this%nuSGS))
        end if 

    end subroutine


    subroutine getRHS_SGS_WallM(this, duidxjC, duidxjE, duidxjChat, urhs, vrhs, wrhs, uhat, vhat, wChat, u, v, wC, ustar, Umn, Vmn, Uspmn, filteredSpeedSq, InvObLength, max_nuSGS, inst_horz_avg, dTdx, dTdy, dTdzHC)    
        class(sgs), intent(inout), target :: this
        real(rkind)   , dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),9), intent(in), target :: duidxjC
        real(rkind)   , dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),4), intent(in), target :: duidxjE
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3),9), intent(inout) :: duidxjChat
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(inout) :: urhs, vrhs
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(in) :: uhat, vhat, wChat
        complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: wrhs
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: u, v, wC
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: filteredSpeedSq
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(inout), optional :: dTdzHC
        real(rkind)   , dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in), optional :: dTdx, dTdy
        real(rkind), dimension(this%ntimeAvgQs), intent(out) :: inst_horz_avg
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz
        real(rkind), dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
        real(rkind), dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
        real(rkind), dimension(:,:,:), pointer :: dvdzC, dudzC
        real(rkind), dimension(:,:,:), pointer :: dwdxC, dwdyC
        real(rkind), dimension(:,:,:), pointer :: tau11, tau12, tau13, tau13C, tau22, tau23, tau23C, tau33
        real(rkind), dimension(:,:,:), pointer :: S11, S12, S13, S13C, S22, S23, S23C, S33
        real(rkind), intent(in) :: ustar, Umn, Vmn, Uspmn, InvObLength
        real(rkind), intent(out), optional :: max_nuSGS
        complex(rkind), dimension(:,:,:), pointer :: tauhat, tauhat2    


        dudx  => duidxjC(:,:,:,1); dudy  => duidxjC(:,:,:,2); dudzC => duidxjC(:,:,:,3); 
        dvdx  => duidxjC(:,:,:,4); dvdy  => duidxjC(:,:,:,5); dvdzC => duidxjC(:,:,:,6); 
        dwdxC => duidxjC(:,:,:,7); dwdyC => duidxjC(:,:,:,8); dwdz  => duidxjC(:,:,:,9); 
        dwdx => duidxjE(:,:,:,1); dwdy => duidxjE(:,:,:,2);
        dudz => duidxjE(:,:,:,3); dvdz => duidxjE(:,:,:,4);

        tau11 => this%rbuff(:,:,:,1); tau12 => this%rbuff(:,:,:,2) ; tau13C => this%rbuff(:,:,:,3)
        tau22 => this%rbuff(:,:,:,4); tau23C => this%rbuff(:,:,:,5); tau33 => this%rbuff(:,:,:,6)
        S11 => this%rbuff(:,:,:,1); S12 => this%rbuff(:,:,:,2); S13C => this%rbuff(:,:,:,3)
        S22 => this%rbuff(:,:,:,4); S23C => this%rbuff(:,:,:,5); S33 => this%rbuff(:,:,:,6)

        tau13 => this%rbuffE(:,:,:,1); tau23 => this%rbuffE(:,:,:,2)
        S13 => this%rbuffE(:,:,:,1); S23 => this%rbuffE(:,:,:,2)
        tauhat => this%cbuff(:,:,:,1); tauhat2 => this%cbuff(:,:,:,2)       

        ! STEP 1: Compute Sij
        S11 = dudx; S22 = dvdy; S33 = dwdz
        S13 = half*(dudz + dwdx); S23 = half*(dvdz + dwdy); S12 = half*(dudy + dvdx)

        S13C = half*(dudzC + dwdxC); S23C = half*(dvdzC + dwdyC)

        ! STEP 2: Call the SGS model operator
        select case (this%SGSmodel)
        case (0)     
            call this%get_SMAG_Op(this%nuSGS, S11,S22,S33,S12,S13C,S23C)
        case (1)
            call this%get_SIGMA_Op(this%nuSGS, dudx, dudy, dudzC, dvdx, dvdy, dvdzC, dwdxC, dwdyC, dwdz)
        case (2)
            call this%get_ShearImpSMAG_Op(this%nuSGS, S11,S22,S33,S12,S13C,S23C)
        case (3)
            call this%get_MGM_Op(duidxjC, duidxjE)!, maxDissp)
        case (4)
            if (this%isStratified) then
                if ((present(dTdx)).and.(present(dTdy)).and.(present(dTdzHC))) then
                    if (this%spectC%carryingZeroK) then
                        dTdzHC(1,1,:) = cmplx(0.d0,0.d0)
                    end if
                    call this%spectC%ifft(dTdzHC,this%dTdzC_diff)
                    call this%get_AMD_Op_strat(this%nuSGS,dudx, dudy, dudzC, dvdx, dvdy, dvdzC, dwdxC, dwdyC, &
                                               dwdz,S11,S22,S33,S12,S13C,S23C, dTdx, dTdy)
                else ! This should never happen if called from igridWallM
                    this%nuSGS = zero        
                end if
            else
                call this%get_AMD_Op(this%nuSGS,dudx, dudy, dudzC, dvdx, dvdy, dvdzC, dwdxC, dwdyC, dwdz, S11,S22,S33,S12,S13C,S23C)
            end if
        end select

        !print*, this%nuSGS(3,2,3)
        if (this%useDynamicProcedure) then
            if (mod(this%mstep,ApplyDynEvery) == 0) then
                call this%DynamicProcedure(u,v,wC,uhat, vhat, wChat, duidxjChat) 
            end if
            this%nuSGS = this%Lij(:,:,:,1) * (this%deltafilter * this%deltafilter) * this%nuSGS  
        elseif (this%eddyViscModel) then
            if (this%useWallFunction) then
                !print*, this%deltaFIlter, this%cSMAG_WALL(3,2,3)
                this%nuSGS = this%cSMAG_WALL*this%nuSGS
            else
                this%nuSGS = this%mconst*this%nuSGS
            end if 
        end if


        this%mstep = this%mstep + 1
        if (present(max_nuSGS)) then
            max_nuSGS = p_maxval(maxval(this%nuSGS))
        end if 

        ! STEP 3: Compute TAUij
        if(this%eddyViscModel) then
          this%nuSGS = -two*this%nuSGS
          tau11 = this%nuSGS*S11; tau22 = this%nuSGS*S22; tau33 = this%nuSGS*S33
          tau12 = this%nuSGS*S12

          call transpose_x_to_y(this%nuSGS,this%rtmpY,this%gpC)
          call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
          if (useCompactFD) then
            call this%derZ_SS%InterpZ_C2E(this%rtmpZ,this%rtmpZE,size(this%rtmpZ,1),size(this%rtmpZ,2))
          else
            call this%OpsNU%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
            this%rtmpzE(:,:,this%nz+1) = two*this%rtmpZ(:,:,this%nz) - this%rtmpzE(:,:,this%nz)
          end if
          call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
          call transpose_y_to_x(this%rtmpYE,this%nuSGSE,this%gpE)
          !this%nuSGSE = max(zero,this%nuSGSE) !-- bug here. nuSGSE here should be non-positive to avoid backscatter
          this%nuSGSE = min(zero,this%nuSGSE) !-- bug here. nuSGSE here should be non-positive to avoid backscatter
          tau13 = this%nuSGSE * S13
          tau23 = this%nuSGSE * S23
        endif
     

        !if(nrank==0) then
        !  print *, "delta", this%deltaFilter, this%mconst
        !  print*, "duidxj:", dudx(3,2,3), dudy(3,2,3), dudzC(3,2,3)
        !  print*, "duidxj:", dvdx(3,2,3), dvdy(3,2,3), dvdzC(3,2,3)
        !  print*, "duidxj:", dwdxC(3,2,3), dwdyC(3,2,3), dwdz(3,2,3)
        !  print*, "tau11:", tau11(3,2,3)
        !  print*, "tau22:", tau22(3,2,3)
        !  print*, "tau33:", tau33(3,2,3)
        !  print*, "tau13:", tau13(3,2,3)
        !  print*, "tau23:", tau23(3,2,3)
        !  print*, "tau12:", tau12(3,2,3)
        !endif

        select case(this%wallModel) 
        case(0)
            this%WallMFactor = -ustar*ustar/(Uspmn + 1.D-13)
            
            ! STEP 4: tau13 -> ddz() in urhs, ddx in wrhs
            call transpose_y_to_z(uhat,this%ctmpCz,this%sp_gp) ! <- send uhat to z decomp (wallmodel)
            call this%spectE%fft(tau13,this%ctmpEy)
            call transpose_y_to_z(this%ctmpEy,this%ctmpEz,this%sp_gpE)
            
            this%ctmpEz(:,:,1) = this%WallMFactor*this%ctmpCz(:,:,1)
            
            if(nrank==0) inst_horz_avg(2) = real(this%ctmpEz(1,1,1), rkind)*this%MeanFact
            if(useCompactFD) then
                call this%derZ_SS%ddz_E2C(this%ctmpEz,this%ctmpCz,size(this%ctmpEz,1),size(this%ctmpEz,2))
            else
                call this%Ops2ndOrder%ddz_E2C(this%ctmpEz,this%ctmpCz)
            end if 
            call transpose_z_to_y(this%ctmpCz,tauhat,this%sp_gp)
            urhs = urhs - tauhat
            if (.not. this%CompStokesP) then
                this%ctmpEz(:,:,1) = zeroC; this%ctmpEz(:,:,this%nz+1) = zeroC
            end if 
            call transpose_z_to_y(this%ctmpEz,this%ctmpEy,this%sp_gpE)
            call this%spectE%mtimes_ik1_ip(this%ctmpEy)
            wrhs = wrhs - this%ctmpEy
            
            ! STEP 5: tau23 -> ddz() in vrhs, ddy in wrhs
            call transpose_y_to_z(vhat,this%ctmpCz,this%sp_gp) ! <- send vhat to z decomp (wallmodel)

            call this%spectE%fft(tau23,this%ctmpEy)
            call transpose_y_to_z(this%ctmpEy,this%ctmpEz,this%sp_gpE)
            
            this%ctmpEz(:,:,1) = this%WallMFactor*this%ctmpCz(:,:,1)
            
            if(nrank==0) inst_horz_avg(3) = real(this%ctmpEz(1,1,1), rkind)*this%MeanFact
            if(useCompactFD) then
                call this%derZ_SS%ddz_E2C(this%ctmpEz,this%ctmpCz,size(this%ctmpEz,1),size(this%ctmpEz,2))
            else
                call this%Ops2ndOrder%ddz_E2C(this%ctmpEz,this%ctmpCz)
            end if 
            call transpose_z_to_y(this%ctmpCz,tauhat,this%sp_gp)
            vrhs = vrhs - tauhat
            if (.not. this%CompStokesP) then
                this%ctmpEz(:,:,1) = zeroC; this%ctmpEz(:,:,this%nz+1) = zeroC
            end if 
            call transpose_z_to_y(this%ctmpEz,this%ctmpEy,this%sp_gpE)
            call this%spectE%mtimes_ik2_ip(this%ctmpEy)
            wrhs = wrhs - this%ctmpEy
        case(1)
            call this%getLocalWallModel(filteredSpeedSq, this%ctmpCz)
            this%WallMFactor = -(kappa/(log(this%dz/(two*this%z0)) + bm*InvObLength*this%dz/two))**2 
            
            call this%spectE%fft(tau13,this%ctmpEy)
            call transpose_y_to_z(this%ctmpEy,this%ctmpEz,this%sp_gpE)
            this%ctmpEz(:,:,1) = (this%WallMFactor*Umn/Uspmn)*this%ctmpCz(:,:,1) 
            if(nrank==0) inst_horz_avg(2) = real(this%ctmpEz(1,1,1), rkind)*this%MeanFact
            if(useCompactFD) then
                call this%derZ_SS%ddz_E2C(this%ctmpEz,this%ctmpCz,size(this%ctmpEz,1),size(this%ctmpEz,2))
            else
                call this%Ops2ndOrder%ddz_E2C(this%ctmpEz,this%ctmpCz)
            end if 
            call transpose_z_to_y(this%ctmpCz,tauhat,this%sp_gp)
            urhs = urhs - tauhat
            if (.not. this%CompStokesP) then
                this%ctmpEz(:,:,1) = zeroC; this%ctmpEz(:,:,this%nz+1) = zeroC
            end if 
            call transpose_z_to_y(this%ctmpEz,this%ctmpEy,this%sp_gpE)
            call this%spectE%mtimes_ik1_ip(this%ctmpEy)
            wrhs = wrhs - this%ctmpEy
            
            ! STEP 5: tau23 -> ddz() in vrhs, ddy in wrhs
            call this%spectE%fft(tau23,this%ctmpEy)
            call transpose_y_to_z(this%ctmpEy,this%ctmpEz,this%sp_gpE)
            this%ctmpEz(:,:,1) = (this%WallMFactor*Vmn/Uspmn)*this%ctmpCz(:,:,1) 
            if(nrank==0) inst_horz_avg(3) = real(this%ctmpEz(1,1,1), rkind)*this%MeanFact
            if(useCompactFD) then
                call this%derZ_SS%ddz_E2C(this%ctmpEz,this%ctmpCz,size(this%ctmpEz,1),size(this%ctmpEz,2))
            else
                call this%Ops2ndOrder%ddz_E2C(this%ctmpEz,this%ctmpCz)
            end if 
            call transpose_z_to_y(this%ctmpCz,tauhat,this%sp_gp)
            vrhs = vrhs - tauhat
            if (.not. this%CompStokesP) then
                this%ctmpEz(:,:,1) = zeroC; this%ctmpEz(:,:,this%nz+1) = zeroC
            end if 
            call transpose_z_to_y(this%ctmpEz,this%ctmpEy,this%sp_gpE)
            call this%spectE%mtimes_ik2_ip(this%ctmpEy)
            wrhs = wrhs - this%ctmpEy
        end select

        ! STEP 6: tau12 -> ddy in urhs, ddx in vrhs
        call this%spectC%fft(tau12,tauhat)
        call this%spectC%mtimes_ik1_oop(tauhat,tauhat2)
        vrhs = vrhs - tauhat2
        call this%spectC%mtimes_ik2_ip(tauhat)
        urhs = urhs - tauhat

        ! STEP 7: tau11 -> ddx() in urhs 
        call this%spectC%fft(tau11,tauhat)
        call this%spectC%mtimes_ik1_ip(tauhat)
        urhs = urhs - tauhat

        ! STEP 8: tau22 -> ddy() in vrhs
        call this%spectC%fft(tau22,tauhat)
        call this%spectC%mtimes_ik2_ip(tauhat)
        vrhs = vrhs - tauhat

        ! STEP 9: tau33 -> ddz() in wrhs
        call this%spectC%fft(tau33,tauhat)
        call transpose_y_to_z(tauhat,this%ctmpCz,this%sp_gp)
        if (useCompactFD) then
            call this%derTAU33%ddz_C2E(this%ctmpCz,this%ctmpEz,size(this%ctmpCz,1),size(this%ctmpCz,2))
        else
            call this%Ops2ndOrder%ddz_C2E(this%ctmpCz,this%ctmpEz,.true.,.true.)
        end if 
        call transpose_z_to_y(this%ctmpEz,this%ctmpEy,this%sp_gpE)
        wrhs = wrhs - this%ctmpEy

    end subroutine

#include "sgs_models/dynamicprocedure.F90"

    subroutine getRHS_SGS_Scalar_WallM(this, duidxjC, dTdxC, dTdyC, dTdzE, dTdzC, T_rhs, wTh_surf)
        class(sgs), intent(inout), target :: this
        real(rkind)   , dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),9), intent(in), target :: duidxjC
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(inout) :: T_rhs
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: dTdxC,dTdyC, dTdzC 
        real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)), intent(in) :: dTdzE
        real(rkind), intent(in) :: wTh_surf
        complex(rkind), dimension(:,:,:), pointer :: tauhat, tauhat2    
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz
        real(rkind), dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
        real(rkind), dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
        !real(rkind) :: wallMfactor  
        
        tauhat => this%cbuff(:,:,:,1); tauhat2 => this%cbuff(:,:,:,2)       
        dudx  => duidxjC(:,:,:,1); dudy  => duidxjC(:,:,:,2); dudz  => duidxjC(:,:,:,3); 
        dvdx  => duidxjC(:,:,:,4); dvdy  => duidxjC(:,:,:,5); dvdz  => duidxjC(:,:,:,6); 
        dwdx  => duidxjC(:,:,:,7); dwdy  => duidxjC(:,:,:,8); dwdz  => duidxjC(:,:,:,9); 

        select case (this%SGSmodel)
        case(0)
            this%nuSCA = this%nuSGS/(two*this%Pr); this%nuSCAE = this%nuSGSE/(two*this%Pr)
        case(1)
            this%nuSCA = this%nuSGS/(two*this%Pr); this%nuSCAE = this%nuSGSE/(two*this%Pr)
        case(2)
            this%nuSCA = this%nuSGS/(two*this%Pr); this%nuSCAE = this%nuSGSE/(two*this%Pr)
        case(3)
            this%nuSCA = this%nuSGS/(two*this%Pr); this%nuSCAE = this%nuSGSE/(two*this%Pr)
        case(4)
            this%nuSCA = (this%cx**2)*(dudx*dTdxC*dTdxC + dvdx*dTdxC*dTdyC + dwdx*dTdxC*dTdzC)
            this%nuSCA = this%nuSCA + (this%cy**2)*(dudy*dTdyC*dTdxC + dvdy*dTdyC*dTdyC + dwdy*dTdyC*dTdzC)
            this%nuSCA = this%nuSCA + (this%cz**2)*(dudz*dTdzC*dTdxC + dvdz*dTdzC*dTdyC + dwdz*dTdzC*dTdzC)
            this%amd_DEN = dTdxC*dTdxC + dTdyC*dTdyC + dTdzC*dTdzC
            this%nuSCA = -this%nuSCA/this%amd_DEN
            this%nuSCA = max(this%nuSCA,zero)
            call transpose_x_to_y(this%nuSGS,this%rtmpY,this%gpC)
            call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
            if (useCompactFD) then
              call this%derZ_SS%InterpZ_C2E(this%rtmpZ,this%rtmpZE,size(this%rtmpZ,1),size(this%rtmpZ,2))
            else
              call this%OpsNU%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              this%rtmpzE(:,:,this%nz+1) = two*this%rtmpZ(:,:,this%nz) - this%rtmpzE(:,:,this%nz)
            end if
            call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
            call transpose_y_to_x(this%rtmpYE,this%nuSCAE,this%gpE)
        end select

        this%q1 = this%nuSCA  * dTdxC
        this%q2 = this%nuSCA  * dTdyC
        this%q3 = this%nuSCAE * dTdzE
        
        call this%spectC%fft(this%q1,tauhat)
        call this%spectC%mtimes_ik1_ip(tauhat)
        T_rhs = T_rhs - tauhat 

        call this%spectC%fft(this%q2,tauhat)
        call this%spectC%mtimes_ik2_ip(tauhat)
        T_rhs = T_rhs - tauhat 

        call this%spectE%fft(this%q3,this%ctmpEy)
        call transpose_y_to_z(this%ctmpEy,this%ctmpEz,this%sp_gpE)
        this%ctmpEz(:,:,1) = zeroC
        if (nrank == 0) then
            this%ctmpEz(1,1,1) = cmplx(wTh_surf/this%MeanFact,zero,rkind)
        end if
        if (useCompactFD) then
            call this%derZ_SS%ddz_E2C(this%ctmpEz,this%ctmpCz,size(this%ctmpCz,1),size(this%ctmpCz,2))
        else 
            call this%Ops2ndOrder%ddz_E2C(this%ctmpEz,this%ctmpCz)
        end if 
        call transpose_z_to_y(this%ctmpCz,tauhat,this%sp_gp)
        T_rhs = T_rhs - tauhat


    end subroutine  

    subroutine getLocalWallModel(this, filteredSpeedSq, tauWallH)
        class(sgs), intent(inout), target :: this
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) ::filteredSpeedSq
        complex(rkind), dimension(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)), intent(out) :: tauWallH

        complex(rkind), dimension(:,:,:), pointer :: cbuffy

        cbuffy => this%cbuff(:,:,:,1)
        call this%spectC%fft(filteredSpeedSq,cbuffy)
        call transpose_y_to_z(cbuffy,tauWallH,this%sp_gp)
    end subroutine  


end module 
