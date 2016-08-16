module GeneralMatEOS

    use kind_parameters, only: rkind
    use constants,       only: half,one,zero
    use EOSMod,          only: eos

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

    contains

        procedure :: init
        procedure :: get_p_devstress
        procedure :: get_T
        procedure :: get_sos

    end type

contains

    subroutine init(this,eosparams)
        class(generaleos), intent(inout) :: this
        real(rkind), dimension(:) :: eosparams

        this%K0   = eosparams(1);     this%Kp   = eosparams(2)
        this%G0   = eosparams(3);     this%Gp   = eosparams(4)
        this%beta = eosparams(5);     this%T0   = eosparams(6)
        this%Cv   = eosparams(7);     this%gam0 = eosparams(8)
        this%qpar = eosparams(9)

        this%invCv = one/this%Cv

    end subroutine

    pure subroutine get_p_devstress(this,rho0,g,rho,entr,p,devstress)
        class(generaleos), intent(in) :: this
        real(rkind),                     intent(in)  :: rho0
        real(rkind), dimension(:,:,:,:), intent(in)  :: g
        real(rkind), dimension(:,:,:),   intent(in)  :: rho, entr
        real(rkind), dimension(:,:,:),   intent(out) :: p
        real(rkind), dimension(:,:,:,:), intent(out) :: devstress

        p = zero!(this%gam-one)*rho*e - this%gam*this%PInf
        devstress = zero!(this%gam-one)*rho*e - this%gam*this%PInf

        real(rkind), dimension(:,:,:,:), intent(in)  :: g
        real(rkind), dimension(:,:,:,:), intent(out) :: finger
        real(rkind), dimension(:,:,:,:), intent(out), optional :: fingersq
        real(rkind), dimension(:,:,:),   intent(out), optional :: trG, trG2, detG

        associate ( g11 => g(:,:,:,1), g12 => g(:,:,:,2), g13 => g(:,:,:,3), &
                    g21 => g(:,:,:,4), g22 => g(:,:,:,5), g23 => g(:,:,:,6), &
                    g31 => g(:,:,:,7), g32 => g(:,:,:,8), g33 => g(:,:,:,9)  )
            finger(:,:,:,1) = g11*g11 + g21*g21 + g31*g31
            finger(:,:,:,2) = g11*g12 + g21*g22 + g31*g32
            finger(:,:,:,3) = g11*g13 + g21*g23 + g31*g33
            finger(:,:,:,4) = g12*g12 + g22*g22 + g32*g32
            finger(:,:,:,5) = g12*g13 + g22*g23 + g32*g33
            finger(:,:,:,6) = g13*g13 + g23*g23 + g33*g33
        end associate

        associate ( GG11 => finger(:,:,:,1), GG12 => finger(:,:,:,2), GG13 => finger(:,:,:,3), &
                    GG21 => finger(:,:,:,2), GG22 => finger(:,:,:,4), GG23 => finger(:,:,:,5), &
                    GG31 => finger(:,:,:,3), GG32 => finger(:,:,:,5), GG33 => finger(:,:,:,6)  )
            if(present(detG)) then
                detG = GG11*(GG22*GG33-GG23*GG32) - GG12*(GG21*GG33-GG31*GG23) + GG13*(GG21*GG32-GG31*GG22)
            endif
            if(present(fingersq)) then
                fingersq(:,:,:,1) = GG11*GG11 + GG12*GG21 + GG13*GG31
                fingersq(:,:,:,2) = GG11*GG12 + GG12*GG22 + GG13*GG32
                fingersq(:,:,:,3) = GG11*GG13 + GG12*GG23 + GG13*GG33
                fingersq(:,:,:,4) = GG21*GG12 + GG22*GG22 + GG23*GG32
                fingersq(:,:,:,5) = GG21*GG13 + GG22*GG23 + GG23*GG33
                fingersq(:,:,:,6) = GG31*GG13 + GG32*GG23 + GG33*GG33
            endif
        end associate

        ! compute inverse of finger
        ! compute inverse of fingersq

        Inv1 = finger(:,:,:,1) + finger(:,:,:,4) + finger(:,:,:,6)
        Inv3 = 
        dedI1fac = this%beta * GI3 * Inv3**(-third)
        dedI2fac = (one-this%beta) * GI3 * Inv3**(-twothird)

        ! 2*rho*dedI3*I3
        dedI3fac = three*this%K0*exp(-1.5D0*(this%Kp-one)*(Inv3**(sixth)-one))*(Inv3**(-sixth)-Inv3**(-third)) &
                 - (two*this%rho0*this%Cv*this%T0*this%gam0**2/this%qpar) * Inv3**(half*(one-this%qpar)) * (one - Inv3**this%qpar) * &
                   (exp(entr)-one) * exp(this%gam0/this%qpar*(one-Inv3**this%qpar))      &
                 + GpI3 * (this%beta*Inv1*Inv3**(-third) + (one-this%beta)*Inv2*Inv3**(-twothird) - three) &
                 + GI3/Inv3 * (this%beta*sixth*Inv1*Inv3**(-third) - 1.5_rkind*(one-this%beta)*Inv2*Inv3**(-twothird)-three)

        devstress = -dedI2fac * fingersq + (dedI1fac + Inv11*dedI2fac)*finger
        p = -third*(devstress(:,:,:,1) + devstress(:,:,:,4) + devstress(:,:,:,6))
        devstress(:,:,:,1) = devstress(:,:,:,1) + p 
        devstress(:,:,:,4) = devstress(:,:,:,4) + p 
        devstress(:,:,:,6) = devstress(:,:,:,6) + p 

        p = p + Inv3*dedI3fac

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
