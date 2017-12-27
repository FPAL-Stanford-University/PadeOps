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
        call this%get_SMAG_Op(this%nuSGS, S11,S22,S33,S12,S13,S23)
        
        if (this%useDynamicProcedure) then
            call this%DynamicProcedure(u,v,wC,uhat, vhat, wChat, duidxj,duidxjhat) 
        else
            if (this%useWallFunction) then
                this%nuSGS = this%cSMAG_WALL*this%nuSGS
            else
                this%nuSGS = this%mconst*this%nuSGS
            end if 
        end if 

        tau11 = -two*this%nuSGS*S11; tau12 = -two*this%nuSGS*S12; tau13 = -two*this%nuSGS*S13
        tau22 = -two*this%nuSGS*S22; tau23 = -two*this%nuSGS*S23; tau33 = -two*this%nuSGS*S33

        tauhat => this%cbuff(:,:,:,1); tauhat2 => this%cbuff(:,:,:,2)       


        if (useCompactFD) then
#include "getRHS_common_CD.F90"
        else
#include "getRHS_common_FD.F90"
        end if
        
        if (present(max_nuSGS)) then
            max_nuSGS = p_maxval(maxval(this%nuSGS))
        end if 

    end subroutine
