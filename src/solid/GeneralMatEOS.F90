module GeneralMatEOS

    use kind_parameters, only: rkind
    use constants,       only: half,one,zero,third,twothird,sixth,three,two
    use EOSMod,          only: eos
    use decomp_2d,       only: decomp_info

    implicit none

    !type, extends(GenSolFluidEos) :: generaleos
    type :: generaleos

        real(rkind) :: K0   = 138.0d9
        real(rkind) :: Kp   = 4.96_rkind
        real(rkind) :: G0   = 46.9d9
        real(rkind) :: Gp   = 0.57_rkind
        real(rkind) :: beta = 0.0_rkind
        real(rkind) :: T0   = 300.0_rkind
        real(rkind) :: Cv   = 3.9d-4
        real(rkind) :: gam0 = 1.96_rkind
        real(rkind) :: qpar = 1.0_rkind
        
        real(rkind) :: invCv

        real(rkind), allocatable, dimension(:,:,:,:) :: finger
        real(rkind), allocatable, dimension(:,:,:,:) :: fingersq
        real(rkind), allocatable, dimension(:,:,:)   :: Inv1,Inv2,Inv3,dedI1fac,dedI2fac,dedI3fac,GI3,GpI3

    contains

        procedure :: init
        procedure :: get_p_devstress
        procedure :: get_T
        procedure :: get_sos
        procedure :: destroy

    end type

contains

    subroutine init(this,decomp,eosparams)
        class(generaleos),         intent(inout) :: this
        type(decomp_info),         intent(in)    :: decomp
        real(rkind), dimension(:), intent(in)    :: eosparams

        this%K0   = eosparams(1);     this%Kp   = eosparams(2)
        this%G0   = eosparams(3);     this%Gp   = eosparams(4)
        this%beta = eosparams(5);     this%T0   = eosparams(6)
        this%Cv   = eosparams(7);     this%gam0 = eosparams(8)
        this%qpar = eosparams(9)

        this%invCv = one/this%Cv

        allocate( this%finger  (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),6) )
        allocate( this%fingersq(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3),6) )
        allocate( this%Inv1    (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
        allocate( this%Inv2    (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
        allocate( this%Inv3    (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
        allocate( this%dedI1fac(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
        allocate( this%dedI2fac(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
        allocate( this%dedI3fac(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
        allocate( this%GI3     (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )
        allocate( this%GpI3    (decomp%ysz(1),decomp%ysz(2),decomp%ysz(3))   )

    end subroutine

    subroutine destroy(this)
        class(generaleos),         intent(inout) :: this
 
        deallocate(this%finger,this%fingersq,this%Inv1,this%Inv2,this%Inv3,this%dedI1fac,this%dedI2fac,this%dedI3fac,this%GI3,this%GpI3)
    end subroutine

    pure subroutine get_p_devstress(this,rho0,g,rho,entr,p,devstress)
        class(generaleos), intent(inout) :: this
        real(rkind),                     intent(in)  :: rho0
        real(rkind), dimension(:,:,:,:), intent(in)  :: g
        real(rkind), dimension(:,:,:),   intent(in)  :: rho, entr
        real(rkind), dimension(:,:,:),   intent(out) :: p
        real(rkind), dimension(:,:,:,:), intent(out) :: devstress

        integer :: ii

        p = zero!(this%gam-one)*rho*e - this%gam*this%PInf
        devstress = zero!(this%gam-one)*rho*e - this%gam*this%PInf

        associate ( g11 => g(:,:,:,1), g12 => g(:,:,:,2), g13 => g(:,:,:,3), &
                    g21 => g(:,:,:,4), g22 => g(:,:,:,5), g23 => g(:,:,:,6), &
                    g31 => g(:,:,:,7), g32 => g(:,:,:,8), g33 => g(:,:,:,9)  )

            ! compute g*gT first
            this%finger(:,:,:,1) = g11*g11 + g12*g12 + g13*g13
            this%finger(:,:,:,2) = g11*g21 + g12*g22 + g13*g23
            this%finger(:,:,:,3) = g11*g31 + g12*g32 + g13*g33
            this%finger(:,:,:,4) = g21*g21 + g22*g22 + g23*g23
            this%finger(:,:,:,5) = g21*g31 + g22*g32 + g23*g33
            this%finger(:,:,:,6) = g31*g31 + g32*g32 + g33*g33

            ! store tr(finger) in Inv2
            this%Inv2 = this%finger(:,:,:,1) + this%finger(:,:,:,4) + this%finger(:,:,:,6)
            

            ! compute inverse of (g*gT) and store in finger. This is right Cauchy-Green tensor, C ( = g^(-T)*g^(-1) = FT*F) )

            ! compute invariants of C
            this%Inv1 = this%finger(:,:,:,1) + this%finger(:,:,:,4) + this%finger(:,:,:,6)

            this%Inv3 = this%finger(:,:,:,1) * (this%finger(:,:,:,4)*this%finger(:,:,:,6) - this%finger(:,:,:,5)*this%finger(:,:,:,5)) &
                      + this%finger(:,:,:,2) * (this%finger(:,:,:,3)*this%finger(:,:,:,5) - this%finger(:,:,:,6)*this%finger(:,:,:,2)) &
                      + this%finger(:,:,:,3) * (this%finger(:,:,:,5)*this%finger(:,:,:,2) - this%finger(:,:,:,3)*this%finger(:,:,:,4)) 

            this%Inv2 = this%Inv2*this%Inv3

            ! compute gT*g first
            this%finger(:,:,:,1) = g11*g11 + g21*g21 + g31*g31
            this%finger(:,:,:,2) = g11*g12 + g21*g22 + g31*g32
            this%finger(:,:,:,3) = g11*g13 + g21*g23 + g31*g33
            this%finger(:,:,:,4) = g12*g12 + g22*g22 + g32*g32
            this%finger(:,:,:,5) = g12*g13 + g22*g23 + g32*g33
            this%finger(:,:,:,6) = g13*g13 + g23*g23 + g33*g33

        end associate

        ! compute inverse of (gT*g) and store in finger. This is left Cauchy-Green tensor, b ( = g^(-1)*g^(-T) = F*FT)
            

        associate ( GG11 => this%finger(:,:,:,1), GG12 => this%finger(:,:,:,2), GG13 => this%finger(:,:,:,3), &
                    GG21 => this%finger(:,:,:,2), GG22 => this%finger(:,:,:,4), GG23 => this%finger(:,:,:,5), &
                    GG31 => this%finger(:,:,:,3), GG32 => this%finger(:,:,:,5), GG33 => this%finger(:,:,:,6)  )

                ! compute square of b and store in fingersq
                this%fingersq(:,:,:,1) = GG11*GG11 + GG12*GG21 + GG13*GG31
                this%fingersq(:,:,:,2) = GG11*GG12 + GG12*GG22 + GG13*GG32
                this%fingersq(:,:,:,3) = GG11*GG13 + GG12*GG23 + GG13*GG33
                this%fingersq(:,:,:,4) = GG21*GG12 + GG22*GG22 + GG23*GG32
                this%fingersq(:,:,:,5) = GG21*GG13 + GG22*GG23 + GG23*GG33
                this%fingersq(:,:,:,6) = GG31*GG13 + GG32*GG23 + GG33*GG33
        end associate

        this%dedI1fac = this%beta * this%GI3 * this%Inv3**(-third)
        this%dedI2fac = (one-this%beta) * this%GI3 * this%Inv3**(-twothird)

        ! 2*rho*dedI3*I3
        this%dedI3fac = three*this%K0*exp(-1.5D0*(this%Kp-one)*(this%Inv3**(sixth)-one))*(this%Inv3**(-sixth)-this%Inv3**(-third)) &
                 - (two*rho0*this%Cv*this%T0*this%gam0**2/this%qpar) * this%Inv3**(half*(one-this%qpar)) * (one - this%Inv3**this%qpar) * &
                   (exp(entr)-one) * exp(this%gam0/this%qpar*(one-this%Inv3**this%qpar))      &
                 + this%GpI3 * (this%beta*this%Inv1*this%Inv3**(-third) + (one-this%beta)*this%Inv2*this%Inv3**(-twothird) - three) &
                 + this%GI3/this%Inv3 * (this%beta*sixth*this%Inv1*this%Inv3**(-third) - 1.5_rkind*(one-this%beta)*this%Inv2*this%Inv3**(-twothird)-three)

        do ii = 1, 6
            devstress(:,:,:,ii) = -this%dedI2fac * this%fingersq(:,:,:,ii) + (this%dedI1fac + this%Inv1*this%dedI2fac)*this%finger(:,:,:,ii)
        enddo
        p = -third*(devstress(:,:,:,1) + devstress(:,:,:,4) + devstress(:,:,:,6))
        devstress(:,:,:,1) = devstress(:,:,:,1) + p 
        devstress(:,:,:,4) = devstress(:,:,:,4) + p 
        devstress(:,:,:,6) = devstress(:,:,:,6) + p 

        p = p + this%Inv3*this%dedI3fac

    end subroutine

    pure subroutine get_T(this,e,T)
        class(generaleos), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: e
        real(rkind), dimension(:,:,:), intent(out) :: T

        T = this%invCv*e

    end subroutine

    pure subroutine get_sos(this,rho0,rho,devstress,sos)
        class(generaleos), intent(in) :: this
        real(rkind),                     intent(in)  :: rho0
        real(rkind), dimension(:,:,:),   intent(in)  :: rho
        real(rkind), dimension(:,:,:,:), intent(in)  :: devstress
        real(rkind), dimension(:,:,:),   intent(out) :: sos

        sos = zero!sqrt(this%gam*(p+this%PInf)/rho)

    end subroutine

!    pure subroutine get_e_from_p(this,rho,p,e)
!        class(generaleos), intent(in) :: this
!        real(rkind), dimension(:,:,:), intent(in)  :: rho,p
!        real(rkind), dimension(:,:,:), intent(out) :: e
!
!        e = zero!(p + this%gam*this%PInf) * this%onebygam_m1 / rho
!
!    end subroutine

end module
