module sgsmod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, pi, zero,one,two,three,half, four,eight, nine, six  
    use decomp_2d
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral  
    use mpi 
    use cd06staggstuff, only: cd06stagg
    use reductions, only: p_maxval, p_sum
    use numerics, only: useCompactFD 
    use StaggOpsMod, only: staggOps  
    use wallModelmod, only: wallmodel

    implicit none

    private
    public :: sgs

    real(rkind) :: c_sigma = 1.5_rkind
    real(rkind) :: c_smag = 0.17_rkind
    real(rkind) :: kappa = 0.41_rkind
    real(rkind), parameter :: deltaRatio = two
    complex(rkind), parameter :: zeroC = zero + imi*zero

    type :: sgs
        private
        type(spectral), pointer :: spect, spectE 
        real(rkind), allocatable, dimension(:,:,:,:) :: rbuff, Lij, Mij
        real(rkind), pointer, dimension(:,:,:) :: cSMAG_WALL, nuSGS, nuSGSfil
        real(rkind) :: deltaFilter, mconst, deltaTFilter 
        type(decomp_info), pointer :: sp_gp, gp
        type(decomp_info), pointer :: sp_gpE, gpE

        complex(rkind), allocatable, dimension(:,:,:,:) :: cbuff
        complex(rkind), allocatable, dimension(:,:,:) :: ctmpCz, ctmpEz, ctmpEy, ctmpCz2
        complex(rkind), pointer, dimension(:,:,:) :: nuSGShat
        
        real(rkind), dimension(:,:,:), allocatable :: rtmpY, rtmpZ

        type(cd06stagg), allocatable :: derZ_EE, derZ_OO
        type(staggOps), allocatable :: Ops2ndOrder

        logical :: useWallFunction = .false. 
        logical :: useDynamicProcedure = .false.
        logical :: useClipping = .false. 
        !type(moengWall), allocatable :: wallModel
        
        real(rkind) :: meanFact

        logical :: useWallModel = .false.  
        type(wallmodel), pointer :: moengWall
        real(rkind) :: UmeanAtWall

        integer :: SGSmodel ! 0: Standard Smag, 1: Sigma Model 

        contains 
            procedure :: init
            procedure :: destroy
            procedure, private :: DynamicProcedure
            procedure, private :: planarAverage 
            procedure, private :: BroadcastMeanAtWall 
            procedure, private :: get_SMAG_Op
            procedure, private :: testFilter_ip
            procedure, private :: testFilter_oop
            procedure, private :: testFilter_oop_C2R
            procedure :: getRHS_SGS
            procedure :: link_pointers
    end type 

contains

#include "sgs_models/initialize.F90"

    subroutine link_pointers(this,nuSGS,c_SGS, tauSGS_ij)
        class(sgs), intent(in), target :: this
        real(rkind), dimension(:,:,:)  , pointer, intent(inout) :: c_SGS, nuSGS
        real(rkind), dimension(:,:,:,:), pointer, intent(inout) :: tauSGS_ij

        if (allocated(this%rbuff)) then
            if (this%useDynamicProcedure) then
                c_SGS => this%rbuff(:,:,:,8)
            else
                c_SGS => this%rbuff(:,:,:,19) !Locations that are certainly zero
            end if 
            nuSGS => this%rbuff(:,:,:,7)
            tauSGS_ij => this%rbuff(:,:,:,13:18)
        else
            call gracefulExit("You have called SGS%LINK_POINTERS before &
                & initializing SGS",324)
        end if 

    end subroutine


    subroutine get_SMAG_Op(this, nuSGS, S11, S22, S33, S12, S13, S23)
        class(sgs), intent(inout) :: this
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in) :: S11, S22, S33
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in) :: S12, S13, S23
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: nuSGS

        nuSGS = S12*S12 + S13*S13 + S23*S23
        nuSGS = two*nuSGS
        nuSGS = nuSGS + S11*S11 + S22*S22 + S33*S33
        nuSGS = two*nuSGS
        nuSGS = sqrt(nuSGS)

    end subroutine

    subroutine getRHS_SGS(this, duidxj, duidxjhat, urhs, vrhs, wrhs, uhat, vhat, wChat, u, v, wC, max_nuSGS)
        class(sgs), intent(inout), target :: this
        real(rkind)   , dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3),9), intent(inout), target :: duidxj
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3),9), intent(inout) :: duidxjhat
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(inout) :: urhs, vrhs
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(in) :: uhat, vhat, wChat
        complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: wrhs
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(inout) :: u, v, wC

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


        if (this%useWallModel) then
            call this%BroadcastMeanAtWall(u,v)
        end if 

#include "sgs_models/getRHS_common.F90"

        if (present(max_nuSGS)) then
            max_nuSGS = p_maxval(maxval(this%nuSGS))
        end if 

    end subroutine

    subroutine BroadcastMeanAtWall(this, u,v)
        use kind_parameters, only: mpirkind
        class(sgs), intent(inout), target :: this
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in) :: u, v
        integer :: ierr
        real(rkind), dimension(:,:,:), pointer :: rbuff
        complex(rkind), dimension(:,:,:), pointer :: cbuff

        rbuff => this%rbuff(:,:,:,1)
        cbuff => this%cbuff(:,:,:,1)

        rbuff = sqrt(u*u + v*v)
        call this%spect%fft(rbuff,cbuff)

        if (nrank == 0) then
            this%UmeanAtWall = real(cbuff(1,1,1),rkind)*this%meanFact
        end if

        call mpi_bcast(this%UmeanAtWall,1,mpirkind,0,mpi_comm_world,ierr)
      
    end subroutine

    subroutine testFilter_ip(this, field)
        class(sgs), intent(inout), target :: this
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(inout) :: field
        complex(rkind), dimension(:,:,:), pointer :: ctmpY 
        
        ctmpY => this%cbuff(:,:,:,1)
        call this%spect%fft(field,ctmpY)
        call this%spect%testFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,field)
    end subroutine
   
    subroutine testFilter_oop(this, fin, fout)
        class(sgs), intent(inout), target :: this
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: fin
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: fout
        complex(rkind), dimension(:,:,:), pointer :: ctmpY 
        
        ctmpY => this%cbuff(:,:,:,1)
        call this%spect%fft(fin,ctmpY)
        call this%spect%testFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,fout)
    end subroutine
   

    subroutine testFilter_oop_C2R(this, fhat, fout)
        class(sgs), intent(in), target :: this
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(in) :: fhat 
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: fout
        complex(rkind), dimension(:,:,:), pointer :: ctmpY 

        ctmpY => this%cbuff(:,:,:,1)
        call this%spect%TestFilter_oop(fhat, ctmpY)
        call this%spect%ifft(ctmpY,fout)

    end subroutine

    subroutine DynamicProcedure(this,u,v,wC,uhat,vhat, wChat, duidxj,duidxjhat)
        class(sgs), intent(inout), target :: this
        real(rkind)   , dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3),9), intent(inout), target :: duidxj
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3),9), intent(inout), target :: duidxjhat
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(inout) :: u, v, wC
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(in) :: uhat, vhat, wChat

        real(rkind), dimension(:,:,:), pointer :: L11, L12, L13, L22, L23, L33
        real(rkind), dimension(:,:,:), pointer :: M11, M12, M13, M22, M23, M33
        complex(rkind), dimension(:,:,:), pointer :: dudxH, dudyH, dudzH 
        complex(rkind), dimension(:,:,:), pointer :: dvdxH, dvdyH, dvdzH
        complex(rkind), dimension(:,:,:), pointer :: dwdxH, dwdyH, dwdzH
        complex(rkind), dimension(:,:,:), pointer :: ctmpY 
        real(rkind), dimension(:,:,:), pointer :: rtmpX
        real(rkind), dimension(:,:,:), pointer :: numerator, denominator        
        integer :: idx
        real(rkind), dimension(:,:,:), pointer :: C_DYN 


        ! STEP 0: Allocate pointers
        C_DYN => this%rbuff(:,:,:,8); numerator => this%Lij(:,:,:,1); denominator => this%Mij(:,:,:,1)
        
        M11 => this%Mij(:,:,:,1); M12 => this%Mij(:,:,:,2); M13 => this%Mij(:,:,:,3);
        M22 => this%Mij(:,:,:,4); M23 => this%Mij(:,:,:,5); M33 => this%Mij(:,:,:,6);
        
        L11 => this%Lij(:,:,:,1); L12 => this%Lij(:,:,:,2); L13 => this%Lij(:,:,:,3);
        L22 => this%Lij(:,:,:,4); L23 => this%Lij(:,:,:,5); L33 => this%Lij(:,:,:,6);
        
        rtmpX => this%rbuff(:,:,:,8); ctmpY => this%cbuff(:,:,:,2)
        
        dudxH => duidxjhat(:,:,:,1); dudyH => duidxjhat(:,:,:,2); dudzH => duidxjhat(:,:,:,3); 
        dvdxH => duidxjhat(:,:,:,4); dvdyH => duidxjhat(:,:,:,5); dvdzH => duidxjhat(:,:,:,6); 
        dwdxH => duidxjhat(:,:,:,7); dwdyH => duidxjhat(:,:,:,8); dwdzH => duidxjhat(:,:,:,9); 


        ! STEP 1: Compute Mij = delta^2 * \testF{Ssmag*Sij)
        do idx = 1,6
            this%Mij(:,:,:,idx) = this%nuSGS * this%rbuff(:,:,:,idx)
            call this%testFilter_ip(this%Mij(:,:,:,idx))
        end do 
        this%Mij = ( this%deltaFilter * this%deltaFilter ) * this%Mij 
     
        ! STEP 2: Compute Mij = Mij - deltaF^2 * \testF{Ssmag}*\testF{Sij} 
        !         NOTE: Lij is used as a temporary storage
         
        ! Step 2a: Compute \testF{Sij} - Stored in Lij
        call this%testFilter_oop_C2R(dudxH,L11)
        call this%testFilter_oop_C2R(dvdyH,L22)
        call this%testFilter_oop_C2R(dwdzH,L33)
        ctmpY = half*(dudyH + dvdxH)
        call this%testFilter_oop_C2R(ctmpY,L12)
        ctmpY = half*(dudzH + dwdxH)
        call this%testFilter_oop_C2R(ctmpY,L13)
        ctmpY = half*(dvdzH + dwdyH)
        call this%testFilter_oop_C2R(ctmpY,L23)

        ! Step 2b: Compute \testF{Ssmag}
        call this%get_SMAG_Op(rtmpX,L11,L22,L33,L12,L13,L23)
       
        ! Step 2c: multiply and add!
        rtmpX = (this%deltaTfilter * this%deltaTfilter) * rtmpX
        do idx = 1,6
            this%Lij(:,:,:,idx) = rtmpx*this%Lij(:,:,:,idx)
        end do 
        this%Mij = this%Mij - this%Lij
        this%Mij = two*this%Mij 


        ! STEP 3: Compute Lij = \testF{ui*uj}
        L11 = u*u ; L12 = u*v ; L13 = u *wC
        L22 = v*v ; L23 = v*wC; L33 = wC*wC
        do idx = 1,6
            call this%testFilter_ip(this%Lij(:,:,:,idx))
        end do 


        ! STEP 4: Compute Lij = Lij - \testF{ui}*\testF{uj}
        !call this%testFilter_ip(u); call this%testFilter_ip(v); call this%testFilter_ip(wC)
        call this%testFilter_oop_C2R(uhat,u)
        call this%testFilter_oop_C2R(vhat,v)
        call this%testFilter_oop_C2R(wChat,wC)
        
        L11 = L11 - u*u ; L12 = L12 - u*v  ; L13 = L13 - u *wC
        L22 = L22 - v*v ; L23 = L23 - v*wC ; L33 = L33 - wC*wC

        ! STEP 5: Compute Numerator = Mij * Lij 
        this%Lij = this%Lij * this%Mij 
        do idx = 2,6
            numerator = numerator + this%Lij(:,:,:,idx)
        end do 
        numerator = numerator + L12 + L13 + L23
        print*, numerator(4,5,3)

        ! STEP 6: Compute Denominator = Mij * Mij 
        this%Mij = this%Mij * this%Mij 
        do idx = 2,6
            denominator = denominator + this%Mij(:,:,:,idx)
        end do 
        denominator = denominator + M12 + M13 + M23

        ! STEP 7: Filter and clip 
        call this%planarAverage(numerator,this%useClipping)
        call this%planarAverage(denominator,.false.)
        denominator = denominator + 1d-14
        
        ! STEP 8: Compute the SMAG constant and the nuSGS
        numerator = numerator/denominator
        this%nuSGS =  numerator*(this%deltaFilter * this%deltaFilter) * this%nuSGS
    end subroutine

    subroutine planarAverage(this,f, useClipping)
        class(sgs), intent(inout) :: this
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(inout) :: f
        logical, intent(in) :: useClipping
        integer :: k
        real(rkind) :: mnVal

        call transpose_x_to_y(f,this%rtmpY,this%gp)
        call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gp)

        do k = 1,this%gp%zsz(3)
            mnVal = P_SUM(sum(this%rtmpZ(:,:,k)))*this%meanFact
            if (useClipping) then
                this%rtmpZ(:,:,k) = max(mnVal, zero)
            else
                this%rtmpZ(:,:,k) = mnVal
            end if 
        end do 

        call transpose_z_to_y(this%rtmpZ,this%rtmpY,this%gp)
        call transpose_y_to_x(this%rtmpY,f,this%gp)

    end subroutine
end module 
