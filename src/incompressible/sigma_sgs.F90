module sigmaSGSmod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, pi, zero,one,two,three,half, four, nine, six  
    use decomp_2d
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral  
    use mpi 
    use cd06staggstuff, only: cd06stagg
    use reductions, only: p_maxval, p_sum

    implicit none

    private
    public :: sigmaSGS

    real(rkind) :: c_sigma = 1.35_rkind
    real(rkind), parameter :: deltaRatio = four**(two/three)

    type :: sigmaSGS
        private
        type(spectral), pointer :: spect, spectE 
        real(rkind), allocatable, dimension(:,:,:,:) :: rbuff, Lij, Mkl
        real(rkind), pointer, dimension(:,:,:) :: G11, G12, G13, G22, G23, G33
        real(rkind), pointer, dimension(:,:,:) :: nuSGS, nuSGSfil
        real(rkind) :: deltaFilter, mconst
        type(decomp_info), pointer :: sp_gp, gp
        type(decomp_info), pointer :: sp_gpE, gpE

        complex(rkind), allocatable, dimension(:,:,:,:) :: cbuff
        complex(rkind), allocatable, dimension(:,:,:) :: ctmpCz, ctmpEz, ctmpEy, ctmpCz2
        complex(rkind), pointer, dimension(:,:,:) :: nuSGShat
        
        real(rkind), dimension(:,:,:), allocatable :: rtmpY, rtmpZ

        type(cd06stagg), allocatable :: derZ

        logical :: useWallModel = .false.
        logical :: useDynamicProcedure = .false.
        logical :: useClipping = .false. 
        !type(moengWall), allocatable :: wallModel
        
        real(rkind) :: meanFact

        contains 
            procedure :: init
            procedure :: destroy
            procedure, private :: get_nuSGS
            procedure, private :: DynamicProcedure
            procedure, private :: planarAverage 
            procedure :: getRHS_SGS
    end type 

contains

    subroutine init(this, spectC, spectE, gpC, gpE, dx, dy, dz, useDynamicProcedure, useClipping)!, WallModel)
        class(sigmaSGS), intent(inout), target :: this
        type(spectral), intent(in), target :: spectC, spectE
        type(decomp_info), intent(in), target :: gpC, gpE
        real(rkind), intent(in) :: dx, dy, dz
        logical, intent(in) :: useDynamicProcedure, useClipping
       ! type(moengWall), intent(in), target, optional :: WallModel

        !if (present(WallModel)) then
        !    this%useWallModel = .true.
        !end if 
        
        this%useDynamicProcedure = useDynamicProcedure
        this%useClipping = useClipping

        allocate(this%rbuff(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),19))
        
        this%spect => spectC
        this%spectE => spectE
        this%deltaFilter = ((1.5_rkind**2)*dx*dy*dz)**(one/three)
        this%mconst = (this%deltaFilter*c_sigma)**2

        this%G11 => this%rbuff(:,:,:,1)
        this%G12 => this%rbuff(:,:,:,2)
        this%G13 => this%rbuff(:,:,:,3)
        this%G22 => this%rbuff(:,:,:,4)
        this%G23 => this%rbuff(:,:,:,5)
        this%G33 => this%rbuff(:,:,:,6)
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

        allocate(this%derZ)
        call this%derZ%init( this%sp_gp%zsz(3), dz, isTopEven = .true., isBotEven = .true., & 
                             isTopSided = .true., isBotSided = .true.) 

         if (this%useDynamicProcedure) then
            allocate(this%Lij(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),6))
            allocate(this%Mkl(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),6))
         end if 

         this%meanFact = one/(gpC%xsz(1)*gpC%ysz(2))
         allocate(this%rtmpY(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3)))
         allocate(this%rtmpZ(gpC%zsz(1),gpC%zsz(2),gpC%zsz(3)))

    end subroutine

    subroutine destroy(this)
        class(sigmaSGS), intent(inout) :: this

        call this%derZ%destroy()
        deallocate(this%derZ)
        nullify(this%G11, this%G13, this%G13, this%G22, this%G23, this%G33, this%nuSGS)
        deallocate(this%Lij, this%Mkl)
        nullify(this%sp_gp, this%gp, this%spect)
        deallocate(this%rbuff, this%cbuff)
        deallocate(this%rtmpZ, this%rtmpY)
    end subroutine

    subroutine get_nuSGS(this,duidxj,nuSGSfil)
        class(sigmaSGS), intent(inout), target :: this
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3),9), intent(in), target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
        real(rkind), dimension(:,:,:), pointer :: I1, I2, I3, I1sq, I1cu
        real(rkind), dimension(:,:,:), pointer :: alpha1, alpha2, alpha3, alpha1sqrt
        real(rkind), dimension(:,:,:), pointer :: alpha1tmp
        real(rkind), dimension(:,:,:), pointer :: sigma1, sigma2, sigma3, sigma1sq 
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out), optional :: nuSGSfil

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3)
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6)
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9)

        
        I1 => this%rbuff(:,:,:,8); I2 => this%rbuff(:,:,:,9); I3 => this%rbuff(:,:,:,10)
        I1sq => this%rbuff(:,:,:,11); I1cu => this%rbuff(:,:,:,12)

        alpha1 => this%rbuff(:,:,:,1); alpha2 => this%rbuff(:,:,:,2); alpha3 => this%rbuff(:,:,:,3)

        this%G11 = dudx*dudx + dvdx*dvdx + dwdx*dwdx 
        this%G12 = dudx*dudy + dvdx*dvdy + dwdx*dwdy 
        this%G13 = dudx*dudz + dvdx*dvdz + dwdx*dwdz
        this%G22 = dudy*dudy + dvdy*dvdy + dwdy*dwdy 
        this%G23 = dudy*dudz + dvdy*dvdz + dwdy*dwdz 
        this%G33 = dudz*dudz + dvdz*dvdz + dwdz*dwdz

        I1 = this%G11 + this%G22 + this%G33
        I1sq = I1*I1
        I1cu = I1sq*I1
        
        I2 = -this%G11*this%G11 - this%G22*this%G22 - this%G33*this%G33
        I2 = I2 - two*this%G12*this%G12 - two*this%G13*this%G13
        I2 = I2 - two*this%G23*this%G23
        I2 = I2 + I1sq
        I2 = half*I2

        I3 = this%G11*(this%G22*this%G33 - this%G23*this%G23) 
        I3 = I3 + this%G12*(this%G13*this%G23 - this%G12*this%G33) 
        I3 = I3 + this%G13*(this%G12*this%G23 - this%G22*this%G13)  

        alpha1 = I1sq/nine - I2/three
        alpha1 = max(alpha1,zero)
        
        alpha2 = I1cu/27._rkind - I1*I2/six + I3/two
       
        alpha1sqrt => this%rbuff(:,:,:,4)
        alpha1sqrt = sqrt(alpha1)    
        alpha1tmp => this%rbuff(:,:,:,5) 
        alpha1tmp = alpha1*alpha1sqrt
        alpha1tmp = alpha2/(alpha1tmp)
        alpha1tmp = min(alpha1tmp,one)
        alpha1tmp = max(alpha1tmp,-one)
        alpha1tmp = acos(alpha1tmp)
        alpha3 = (one/three)*(alpha1tmp)
          
  
        sigma1 => this%rbuff(:,:,:,9); sigma2 => this%rbuff(:,:,:,10); sigma3 => this%rbuff(:,:,:,11)
        sigma1sq => this%rbuff(:,:,:,12)

        sigma1sq = I1/three + two*alpha1sqrt*cos(alpha3)
        sigma1sq = max(sigma1sq,zero)
        sigma1 = sqrt(sigma1sq)

        sigma2 = pi/three + alpha3
        sigma2 = (-two)*alpha1sqrt*cos(sigma2)
        sigma2 = sigma2 + I1/three
        sigma2 = max(sigma2,zero)
        sigma2 = sqrt(sigma2)
        
        sigma3 = pi/three - alpha3
        sigma3 = (-two)*alpha1sqrt*cos(sigma3)
        sigma3 = sigma3 + I1/three
        sigma3 = max(sigma3,zero)
        sigma3 = sqrt(sigma3)

        if (present(nuSGSfil)) then
            nuSGSfil = sigma3*(sigma1 - sigma2)*(sigma2 - sigma3)/(sigma1sq + 1.d-15)
        else
            this%nuSGS = sigma3*(sigma1 - sigma2)*(sigma2 - sigma3)/(sigma1sq + 1.d-15)
            call this%spect%fft(this%nuSGS,this%nuSGShat)
            call this%spect%dealias(this%nuSGShat)
            call this%spect%ifft(this%nuSGShat,this%nuSGS)  
        end if 
        

    end subroutine

    subroutine getRHS_SGS(this, duidxj, urhs, vrhs, wrhs, uhat, vhat, wChat, u, v, wC, max_nuSGS)
        class(sigmaSGS), intent(inout), target :: this
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

        call this%get_nuSGS(duidxj)

        ! COmpute S_ij (here denoted as tau_ij)
        tau11 = dudx
        tau12 = half*(dvdx + dudy)
        tau13 = half*(dwdx + dudz)
        tau22 = dvdy
        tau23 = half*(dvdz + dwdy)
        tau33 = dwdz   
        
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
#include "sgs_models/getRHS_common.F90"

        if (present(max_nuSGS)) then
            max_nuSGS = p_maxval(maxval(this%nuSGS))
        end if 

    end subroutine

    subroutine DynamicProcedure(this,uhat,vhat,wChat,u,v,wC,duidxj)
        class(sigmaSGS), intent(inout), target :: this
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
    
        ! Step 5: Compute Mkl: get \tilde{Dm}
        call this%get_nuSGS(duidxj,this%nuSGSfil)
        this%nuSGSfil = deltaRatio * this%nuSGSfil

        ! Step 6: Compute Mkl: Compute \tilde{Skl}
        Sf11 = dudx
        Sf12 = half*(dvdx + dudy)
        Sf13 = half*(dwdx + dudz)
        Sf22 = dvdy
        Sf23 = half*(dvdz + dwdy)
        Sf33 = dwdz   

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
        class(sigmasgs), intent(inout) :: this
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
