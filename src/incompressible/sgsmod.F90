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

    real(rkind) :: c_sigma = 1.35_rkind
    real(rkind) :: c_smag = 0.165_rkind
    real(rkind), parameter :: deltaRatio = four**(two/three)
    complex(rkind), parameter :: zeroC = zero + imi*zero

    type :: sgs
        private
        type(spectral), pointer :: spect, spectE 
        real(rkind), allocatable, dimension(:,:,:,:) :: rbuff, Lij, Mkl
        real(rkind), pointer, dimension(:,:,:) :: nuSGS, nuSGSfil
        real(rkind) :: deltaFilter, mconst
        type(decomp_info), pointer :: sp_gp, gp
        type(decomp_info), pointer :: sp_gpE, gpE

        complex(rkind), allocatable, dimension(:,:,:,:) :: cbuff
        complex(rkind), allocatable, dimension(:,:,:) :: ctmpCz, ctmpEz, ctmpEy, ctmpCz2
        complex(rkind), pointer, dimension(:,:,:) :: nuSGShat
        
        real(rkind), dimension(:,:,:), allocatable :: rtmpY, rtmpZ

        type(cd06stagg), allocatable :: derZ_EE, derZ_OO
        type(staggOps), allocatable :: Ops2ndOrder

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
            procedure, private :: get_nuSGS
            procedure, private :: DynamicProcedure
            procedure, private :: planarAverage 
            procedure, private :: BroadcastMeanAtWall 
            procedure :: getRHS_SGS
            procedure :: link_pointers
    end type 

contains


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

    subroutine init(this, ModelID, spectC, spectE, gpC, gpE, dx, dy, dz, useDynamicProcedure, useClipping, moengWall)
        class(sgs), intent(inout), target :: this
        type(spectral), intent(in), target :: spectC, spectE
        type(decomp_info), intent(in), target :: gpC, gpE
        real(rkind), intent(in) :: dx, dy, dz
        logical, intent(in) :: useDynamicProcedure, useClipping
        integer, intent(in) :: modelID
        type(wallModel), intent(in), target, optional :: moengWall

        if (present(moengWall)) then
            this%moengWall => moengWall
            this%useWallModel = .true. 
        end if

        this%SGSmodel = modelID 
        this%useDynamicProcedure = useDynamicProcedure
        this%useClipping = useClipping

        allocate(this%rbuff(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),21))
        this%rbuff = zero  
        this%spect => spectC
        this%spectE => spectE
        this%deltaFilter = ((1.5_rkind**2)*dx*dy*dz)**(one/three)

        select case (this%SGSmodel)
        case(0)
            this%mconst = (this%deltaFilter*c_smag)**2
        case(1)
            this%mconst = (this%deltaFilter*c_sigma)**2
        case default 
            call GracefulExit("Invalid choice for SGS model.",2013)
        end select

        this%nuSGS => this%rbuff(:,:,:,7)
        this%nuSGSfil => this%rbuff(:,:,:,19)
        this%gp => gpC
        this%gpE => gpE
        this%sp_gp => this%spect%spectdecomp
        this%sp_gpE => this%spectE%spectdecomp
        
        allocate(this%cbuff(this%sp_gp%ysz(1), this%sp_gp%ysz(2), this%sp_gp%ysz(3),3))
        this%nuSGShat => this%cbuff(:,:,:,1)

        allocate(this%ctmpCz(this%sp_gp%zsz(1), this%sp_gp%zsz(2), this%sp_gp%zsz(3)))
        allocate(this%ctmpEz(this%sp_gpE%zsz(1), this%sp_gpE%zsz(2), this%sp_gpE%zsz(3)))
        allocate(this%ctmpEy(this%sp_gpE%ysz(1), this%sp_gpE%ysz(2), this%sp_gpE%ysz(3)))
        
        allocate(this%ctmpCz2(this%sp_gp%zsz(1), this%sp_gp%zsz(2), this%sp_gp%zsz(3)))

        if (useCompactFD) then
            allocate(this%derZ_EE, this%derZ_OO)
            call this%derZ_EE%init( this%sp_gp%zsz(3), dz, isTopEven = .true., isBotEven = .true., & 
                             isTopSided = .true., isBotSided = .true.) 
            call this%derZ_OO%init( this%sp_gp%zsz(3), dz, isTopEven = .false., isBotEven = .false., & 
                             isTopSided = .true., isBotSided = .true.) 
        else
            allocate(this%Ops2ndOrder)
            call this%Ops2ndOrder%init(gpC,gpE,0,dx,dy,dz,spectC%spectdecomp,spectE%spectdecomp)
        end if 
        if (this%useDynamicProcedure) then
           allocate(this%Lij(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),6))
           allocate(this%Mkl(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),6))
        end if 

        this%meanFact = one/(real(gpC%xsz(1))*real(gpC%ysz(2)))
        allocate(this%rtmpY(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3)))
        allocate(this%rtmpZ(gpC%zsz(1),gpC%zsz(2),gpC%zsz(3)))

    end subroutine

    subroutine destroy(this)
        class(sgs), intent(inout) :: this

        if (useCompactFD) then
            call this%derZ_OO%destroy()
            call this%derZ_EE%destroy()
            deallocate(this%derZ_OO, this%derZ_EE)
        else
            call this%Ops2ndOrder%destroy()
            deallocate(this%Ops2ndOrder)
        end if     
        nullify( this%nuSGS)
        if(this%useDynamicProcedure) deallocate(this%Lij, this%Mkl)
        nullify(this%sp_gp, this%gp, this%spect)
        deallocate(this%rbuff, this%cbuff)
        deallocate(this%rtmpZ, this%rtmpY)
    end subroutine

    subroutine get_nuSGS(this,duidxj,nuSGSfil)
        class(sgs), intent(inout), target :: this
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3),9), intent(in), target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
        real(rkind), dimension(:,:,:), pointer :: I1, I2, I3, I1sq, I1cu
        real(rkind), pointer, dimension(:,:,:) :: G11, G12, G13, G22, G23, G33
        real(rkind), dimension(:,:,:), pointer :: alpha1, alpha2, alpha3, alpha1sqrt
        real(rkind), dimension(:,:,:), pointer :: alpha1tmp
        real(rkind), dimension(:,:,:), pointer :: sigma1, sigma2, sigma3, sigma1sq 
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out), target, optional :: nuSGSfil
        real(rkind), pointer, dimension(:,:,:) :: S11, S12, S13, S22, S23, S33

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3)
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6)
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9)


        select case (this%SGSmodel)
        case (0) ! Standard Smagorinsky
#include "sgs_models/smagorinsky_model_get_nuSGS.F90"
        case (1) ! Standard Sigma
#include "sgs_models/sigma_model_get_nuSGS.F90"
        case default
            this%nuSGSfil = zero
            if (present(nuSGSfil)) nuSGSfil = zero
        end select

    end subroutine

    subroutine getRHS_SGS(this, duidxj, urhs, vrhs, wrhs, uhat, vhat, wChat, u, v, wC, max_nuSGS)
        class(sgs), intent(inout), target :: this
        real(rkind)   , dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3),9), intent(inout), target :: duidxj
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(inout) :: urhs, vrhs
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(in) :: uhat, vhat, wChat
        complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: wrhs
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(inout) :: u, v, wC

        real(rkind), dimension(:,:,:), pointer :: tau11, tau12, tau13, tau22, tau23, tau33
        complex(rkind), dimension(:,:,:), pointer :: tauhat, tauhat2    
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
        real(rkind), intent(out), optional :: max_nuSGS

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3)
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6)
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9)
    

        tau11 => this%rbuff(:,:,:,13); tau12 => this%rbuff(:,:,:,14); tau13 => this%rbuff(:,:,:,15)
        tau22 => this%rbuff(:,:,:,16); tau23 => this%rbuff(:,:,:,17); tau33 => this%rbuff(:,:,:,18)


        ! COmpute S_ij (here denoted as tau_ij)
        tau11 = dudx
        tau12 = half*(dvdx + dudy)
        tau13 = half*(dwdx + dudz)
        tau22 = dvdy
        tau23 = half*(dvdz + dwdy)
        tau33 = dwdz   
        
        call this%get_nuSGS(duidxj)
        
        if (this%useDynamicProcedure) then
            call this%DynamicProcedure(uhat,vhat,wChat,u,v,wC,duidxj) 
        else
            this%nuSGS = this%mconst*this%nuSGS
        end if 

        tau11 = two*this%nuSGS*tau11
        tau12 = two*this%nuSGS*tau12
        tau13 = two*this%nuSGS*tau13
        tau22 = two*this%nuSGS*tau22
        tau23 = two*this%nuSGS*tau23
        tau33 = two*this%nuSGS*tau33

        tauhat => this%cbuff(:,:,:,1); tauhat2 => this%cbuff(:,:,:,2)       


        if (this%useWallModel) then
            call this%BroadcastMeanAtWall(uhat)
            call this%moengWall%updateWallStress(tau13, tau23, u, v, this%UmeanAtWall)
        end if 

#include "sgs_models/getRHS_common.F90"

        if (present(max_nuSGS)) then
            max_nuSGS = p_maxval(maxval(this%nuSGS))
        end if 

    end subroutine

    subroutine BroadcastMeanAtWall(this, uhat)
        use kind_parameters, only: mpirkind
        class(sgs), intent(inout) :: this
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(in) :: uhat
        integer :: ierr

        if (nrank == 0) then
            this%UmeanAtWall = real(uhat(1,1,1),rkind)*this%meanFact
        end if

        call mpi_bcast(this%UmeanAtWall,1,mpirkind,0,mpi_comm_world,ierr)
      
    end subroutine

    subroutine DynamicProcedure(this,uhat,vhat,wChat,u,v,wC,duidxj)
        class(sgs), intent(inout), target :: this
        real(rkind)   , dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3),9), intent(inout), target :: duidxj
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(in) :: uhat, vhat, wChat
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(inout) :: u, v, wC

        real(rkind), dimension(:,:,:), pointer :: tau11, tau12, tau13, tau22, tau23, tau33
        real(rkind), dimension(:,:,:), pointer :: Sf11, Sf12, Sf13, Sf22, Sf23, Sf33
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
        real(rkind), dimension(:,:,:), pointer :: L11, L12, L13, L22, L23, L33
        real(rkind), dimension(:,:,:), pointer :: M11, M12, M13, M22, M23, M33
        complex(rkind), dimension(:,:,:), pointer :: ctmpY 
        real(rkind), dimension(:,:,:), pointer :: numerator, denominator        


        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3)
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6)
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9)
    
        Sf11 => this%rbuff(:,:,:,1); Sf12 => this%rbuff(:,:,:,2); Sf13 => this%rbuff(:,:,:,3)
        Sf22 => this%rbuff(:,:,:,4); Sf23 => this%rbuff(:,:,:,5); Sf33 => this%rbuff(:,:,:,6)

        tau11 => this%rbuff(:,:,:,13); tau12 => this%rbuff(:,:,:,14); tau13 => this%rbuff(:,:,:,15)
        tau22 => this%rbuff(:,:,:,16); tau23 => this%rbuff(:,:,:,17); tau33 => this%rbuff(:,:,:,18)

        L11 => this%Lij(:,:,:,1); L12 => this%Lij(:,:,:,2); L13 => this%Lij(:,:,:,3);
        L22 => this%Lij(:,:,:,4); L23 => this%Lij(:,:,:,5); L33 => this%Lij(:,:,:,6);
        
        M11 => this%Mkl(:,:,:,1); M12 => this%Mkl(:,:,:,2); M13 => this%Mkl(:,:,:,3);
        M22 => this%Mkl(:,:,:,4); M23 => this%Mkl(:,:,:,5); M33 => this%Mkl(:,:,:,6);
        
        ctmpY => this%cbuff(:,:,:,1)
        numerator => this%rbuff(:,:,:,8); denominator => this%rbuff(:,:,:,9)

        ! NOTE: duidxj, u, v, wC are corrupted fields after this subroutine
        ! call. 

        ! Step 1: Compute Lij : \tide{ui uj}
        L11 = u*u
        call this%spect%fft(L11,ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,L11)

        L12 = u*v
        call this%spect%fft(L12,ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,L12)

        L13 = u*wC
        call this%spect%fft(L13,ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,L13)

        L22 = v*v
        call this%spect%fft(L22,ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,L22)

        L23 = v*wC
        call this%spect%fft(L23,ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,L23)

        L33 = wC*wC
        call this%spect%fft(L33,ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,L33)

        ! Step 2: Compute Lij : \tilde{ui}*\tilde{uj}
        call this%spect%TestFilter_oop(uhat,ctmpY)
        call this%spect%ifft(ctmpY,u)  ! <= u is not corrupted 

        call this%spect%TestFilter_oop(vhat,ctmpY)
        call this%spect%ifft(ctmpY,v)  ! <= v is not corrupted 
        
        call this%spect%TestFilter_oop(wChat,ctmpY)
        call this%spect%ifft(ctmpY,wC)  ! <= wC is not corrupted 
       
        L11 = L11 - u*u; L12 = L12 - u*v ; L13 = L13 -  u*wC
        L22 = L22 - v*v; L23 = L23 - v*wC; L33 = L33 - wC*wC
        
        ! Step 3: Compute Mkl: \tilde{Dm * S_kl}
        ! Note that we have already computed S_ij (stored in Tau_ij buffers)
        ! and also computed Dm which is stored in this%nuSGS
        M11 = -this%nuSGS * tau11
        call this%spect%fft(M11,ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,M11)

        M12 = -this%nuSGS * tau12
        call this%spect%fft(M12,ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,M12)

        M13 = -this%nuSGS * tau13
        call this%spect%fft(M13,ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,M13)

        M22 = -this%nuSGS * tau22
        call this%spect%fft(M22,ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,M22)

        M23 = -this%nuSGS * tau23
        call this%spect%fft(M23,ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,M23)

        M33 = -this%nuSGS * tau33
        call this%spect%fft(M33,ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,M33)

        ! Step 4: Compute Mkl: get \tide{duidxj}
        call this%spect%fft(dudx, ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,dudx)  ! <= dudx is corrupted

        call this%spect%fft(dudy, ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,dudy)  ! <= dudy is corrupted

        call this%spect%fft(dudz, ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,dudz)  ! <= dudz is corrupted

        call this%spect%fft(dvdx, ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,dvdx)  ! <= dvdx is corrvpted

        call this%spect%fft(dvdy, ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,dvdy)  ! <= dvdy is corrvpted

        call this%spect%fft(dvdz, ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,dvdz)  ! <= dvdz is corrvpted

        call this%spect%fft(dwdx, ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,dwdx)  ! <= dwdx is corrwpted

        call this%spect%fft(dwdy, ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,dwdy)  ! <= dwdy is corrwpted

        call this%spect%fft(dwdz, ctmpY)
        call this%spect%TestFilter_ip(ctmpY)
        call this%spect%ifft(ctmpY,dwdz)  ! <= dwdz is corrwpted
    

        ! Step 5: Compute Mkl: Compute \tilde{Skl}
        Sf11 = dudx
        Sf12 = half*(dvdx + dudy)
        Sf13 = half*(dwdx + dudz)
        Sf22 = dvdy
        Sf23 = half*(dvdz + dwdy)
        Sf33 = dwdz   

        ! Step 6: Compute Mkl: get \tilde{Dm}
        call this%get_nuSGS(duidxj,this%nuSGSfil)
        this%nuSGSfil = deltaRatio * this%nuSGSfil
        
        ! Step 7: Add deltaRat*\tilde{Dm}*\tilde{Skl} to Mkl
        M11 = M11 + this%nuSGSfil*Sf11
        M12 = M12 + this%nuSGSfil*Sf12
        M13 = M13 + this%nuSGSfil*Sf13
        M22 = M22 + this%nuSGSfil*Sf22
        M23 = M23 + this%nuSGSfil*Sf23
        M33 = M33 + this%nuSGSfil*Sf33


        ! Step 8: Compute the numerator: Lij * Mij
        numerator = M11*L11 + M22*L22 + M33*L33
        numerator = numerator + two*M12*L12 + two*M13*L13 + two*M23*L23

        ! Step 9: Compute the denominator: Mij * Mij
        denominator = M11*M11 + M22*M22 + M33*M33
        denominator = denominator + two*M12*M12 + two*M13*M13 + two*M23*M23

        ! Step 10: Compute planar averages
        call this%planarAverage(numerator, this%useClipping)
        call this%planarAverage(denominator,.false.)

        ! Step 11: Compute the true nuSGS
        numerator = -half*numerator/(denominator + 1d-14)
        this%nuSGS = numerator*this%nuSGS
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
                this%rtmpZ(:,:,k) = min(mnVal, zero)
            else
                this%rtmpZ(:,:,k) = mnVal
            end if 
        end do 

        call transpose_z_to_y(this%rtmpZ,this%rtmpY,this%gp)
        call transpose_y_to_x(this%rtmpY,f,this%gp)

    end subroutine
end module 
