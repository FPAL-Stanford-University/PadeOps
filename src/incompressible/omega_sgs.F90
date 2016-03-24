module omegaSGSmod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, pi, zero,one,two,three,half, nine, six  
    use decomp_2d
    use exits, only: GracefulExit, message
    use spectralMod, only: spectral  
    use mpi 


    implicit none

    private
    public :: omegaSGS

    real(rkind) :: c_sigma = 1.35_rkind

    type :: omegaSGS
        private
        type(spectral), pointer :: spect, spectE 
        real(rkind), allocatable, dimension(:,:,:,:) :: rbuff
        real(rkind), pointer, dimension(:,:,:) :: G11, G12, G13, G22, G23, G33
        real(rkind), pointer, dimension(:,:,:) :: nuSGS
        real(rkind) :: deltaFilter, mconst
        type(decomp_info), pointer :: sp_gp, gp
        type(decomp_info), pointer :: sp_gpE, gpE

        complex(rkind), allocatable, dimension(:,:,:,:) :: cbuff
        complex(rkind), allocatable, dimension(:,:,:) :: ctmpC, ctmpE
        complex(rkind), pointer, dimension(:,:,:) :: nuSGShat

        contains 
            procedure :: init
            procedure :: destroy
            procedure, private :: get_nuSGS
            procedure :: getRHS_SGS
    end type 

contains

    subroutine init(this, spectC, spectE, gpC, gpE, dx, dy, dz)
        class(omegaSGS), intent(inout), target :: this
        type(spectral), intent(in), target :: spectC, spectE
        type(decomp_info), intent(in), target :: gpC, gpE
        real(rkind), intent(in) :: dx, dy, dz

        allocate(this%rbuff(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3),15))
        
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
        this%gp => gpC
        this%gpE => gpE
        this%sp_gp => this%spect%spectdecomp
        this%sp_gpE => this%spectE%spectdecomp
        
        allocate(this%cbuff(this%sp_gp%ysz(1), this%sp_gp%ysz(2), this%sp_gp%ysz(3),3))
        this%nuSGShat => this%cbuff(:,:,:,1)

        allocate(this%ctmpC(this%sp_gp%zsz(1), this%sp_gp%zsz(2), this%sp_gp%zsz(3)))
        allocate(this%ctmpE(this%sp_gpE%zsz(1), this%sp_gpE%zsz(2), this%sp_gpE%zsz(3)))

    end subroutine

    subroutine destroy(this)
        class(omegaSGS), intent(inout) :: this

        nullify(this%G11, this%G13, this%G13, this%G22, this%G23, this%G33, this%nuSGS)
        nullify(this%sp_gp, this%gp, this%spect)
        deallocate(this%rbuff, this%cbuff)
    end subroutine

    subroutine get_nuSGS(this,duidxj)
        class(omegaSGS), intent(inout), target :: this
        real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3),9), intent(in), target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
        real(rkind), dimension(:,:,:), pointer :: I1, I2, I3, I1sq, I1cu
        real(rkind), dimension(:,:,:), pointer :: alpha1, alpha2, alpha3, alpha1sqrt
        real(rkind), dimension(:,:,:), pointer :: tmp1, tmp2, tmp3, alpha1tmp
        real(rkind), dimension(:,:,:), pointer :: sigma1, sigma2, sigma3, sigma1sq 
        integer :: i, j, k, nx, ny, nz, info

        real(rkind), dimension(3,3) :: Gloc
        real(rkind), dimension(3) :: lambda

        integer, parameter :: lwork = 8
        real(rkind), dimension(lwork) :: work 

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
        alpha1tmp = alpha2/alpha1tmp
        alpha1tmp = min(alpha1tmp,one)
        alpha1tmp = max(alpha1tmp,-one)
        alpha3 = (one/three)*acos(alpha1tmp)
          
  
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

        this%nuSGS = sigma3*(sigma1 - sigma2)*(sigma2 - sigma3)/sigma1sq
        this%nuSGS = this%mconst*this%nuSGS

        call this%spect%fft(this%nuSGS,this%nuSGShat)
        call this%spect%dealias(this%nuSGShat)
        call this%spect%ifft(this%nuSGShat,this%nuSGS)  

    end subroutine

    subroutine getRHS_SGS(this,uhat,vhat,what, duidxj, urhs, vrhs, wrhs)
        class(omegaSGS), intent(inout) :: this
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(in) :: uhat, vhat
        complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(in) :: what
        real(rkind)   , dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3),9), intent(in) :: duidxj

        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(inout) :: urhs, vrhs
        complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: wrhs
        
        real(rkind), dimension(:,:,:), pointer :: tau11, tau12, tau13, tau22, tau23, tau33
        complex(rkind), dimension(:,:,:), pointer :: tauhat                

        !dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3)
        !dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6)
        !dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9)
    
        !tauhat => this%cbuff(:,:,:,1)        

        !tau11 => this%rbuff(:,:,:,1); tau12 => this%rbuff(:,:,:,2); tau13 => this%rbuff(:,:,:,3)
        !tau22 => this%rbuff(:,:,:,4); tau23 => this%rbuff(:,:,:,5); tau33 => this%rbuff(:,:,:,6)

        !call this%get_nuSGS(duidxj)
        !
        !tau11 = two*this%nuSGS*dudx
        !tau12 = this%nuSGS*(dvdx + dudy)
        !tau13 = this%nuSGS*(dwdx + dudz)
        !
        !tau22 = two*this%nuSGS*dvdy
        !tau23 = this%nuSGS*(dvdz + dwdy)

        !tau33 = two*this%nuSGS*dwdz    
       
        !call this%spect%fft(tau11,tauhat)
        !call this%spect%mtimes_ik1_ip(tauhat)
        !urhs = urhs + tauhat

        !call this%spect%fft(tau22,tauhat)
        !call this%spect%mtimes_ik2_ip(tauhat)
        !vrhs = vrhs + tauhat
        !
        !call this%spect%fft(tau33,tauhat)
        !call transpose_y_to_z(tauhat,this%ctmpC,this%sp_gp)


 
    end subroutine



end module 
