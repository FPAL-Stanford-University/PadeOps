module FingerEOSMod

    use kind_parameters, only: rkind,clen
    use constants,       only: one
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit
    use AbstractEOSMod,  only: abstracteos, get_invariants

    implicit none

    ! This is an abstract that holds all the logic for EOS that depend on invariants of the Finger tensor
    type, abstract, extends(abstracteos) :: fingereos

    contains

        procedure :: get_e_from_rho_g_T
        procedure :: get_p_devstress_T_sos2
        procedure :: get_pT_derivatives_wrt_energyVF
        procedure :: get_pT_from_energyVF

        procedure :: get_finger

        procedure(get_energy_derivatives_interface),      deferred :: get_energy_derivatives
        procedure(get_e_from_rho_invariants_T_interface), deferred :: get_e_from_rho_invariants_T

    end type

    abstract interface

        ! Subroutine to compute the dervatives of the internal energy with respect to the invariants of the finger tensor and
        ! entropy (in the form of the temperature)
        subroutine get_energy_derivatives_interface(this,rho,e,I1,I2,I3,dedI1,dedI2,dedI3,T)
            import :: fingereos
            import :: rkind
            class(fingereos),                                            intent(in)  :: this
            real(rkind), dimension(:,:,:),                               intent(in)  :: rho
            real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)), intent(in)  :: e, I1, I2, I3
            real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)), intent(out) :: dedI1, dedI2, dedI3, T
        end subroutine

        ! Subroutine to get the internal energy from the temperature and material density and invariants of the finger tensor
        subroutine get_e_from_rho_invariants_T_interface(this,rho,I1,I2,I3,T,e)
            import :: fingereos
            import :: rkind
            class(fingereos),                                              intent(in)  :: this
            real(rkind), dimension(:,:,:),                                 intent(in)  :: rho
            real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)),   intent(in)  :: I1,I2,I3,T
            real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)),   intent(out) :: e
        end subroutine

    end interface

    interface fingereos
        module procedure init
    end interface

contains

    function init() result(this)
        class(fingereos), pointer :: this
        this%rho0 = one
        this%usegTg = .TRUE.
    end function

    subroutine get_e_from_rho_g_T(this,rho,g,T,e)
        class(fingereos),                                      intent(in)  :: this
        real(rkind), dimension(:,:,:,:),                       intent(in)  :: g
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)), intent(in)  :: rho,T
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)), intent(out) :: e
        
        real(rkind), dimension(size(g,1),size(g,2),size(g,3))   :: I1, I2, I3
        real(rkind), dimension(size(g,1),size(g,2),size(g,3),6) :: finger

        ! Get the finger tensor and square of the finger tensor
        call this%get_finger(g,finger)

        ! Get invariants of the Finger tensor
        call get_invariants(finger,I1,I2,I3)
        
        ! Get energy from density, invariants and T now (deffered to concrete class)
        call this%get_e_from_rho_invariants_T(rho,I1,I2,I3,T,e)
    end subroutine

    subroutine get_p_devstress_T_sos2(this,g,rho,e,p,T,devstress,sos_sq)
        use constants, only: two, three, four
        class(fingereos),                                        intent(in)  :: this
        real(rkind), dimension(:,:,:,:),                         intent(in)  :: g
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)),   intent(in)  :: rho, e
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)),   intent(out) :: p, T, sos_sq
        real(rkind), dimension(size(g,1),size(g,2),size(g,3),6), intent(out) :: devstress

        real(rkind), dimension(size(g,1),size(g,2),size(g,3))   :: I1, I2, I3
        real(rkind), dimension(size(g,1),size(g,2),size(g,3))   :: dedI1, dedI2, dedI3
        real(rkind), dimension(size(g,1),size(g,2),size(g,3),6) :: finger, fingersq

        ! Get the finger tensor and square of the finger tensor
        call this%get_finger(g,finger,fingersq)

        ! Get invariants of the Finger tensor
        call get_invariants(finger,I1,I2,I3)

        ! Get the derivatives of the internal energy wrt invariants of the finger tensor
        call this%get_energy_derivatives(rho,e,I1,I2,I3,dedI1,dedI2,dedI3,T)

        ! Compute the Cauchy stress tensor
        devstress(:,:,:,1) = -two*rho*( (dedI1 + I1*dedI2)*finger(:,:,:,1) - dedI2*fingersq(:,:,:,1) + I3*dedI3 )
        devstress(:,:,:,2) = -two*rho*( (dedI1 + I1*dedI2)*finger(:,:,:,2) - dedI2*fingersq(:,:,:,2)            )
        devstress(:,:,:,3) = -two*rho*( (dedI1 + I1*dedI2)*finger(:,:,:,3) - dedI2*fingersq(:,:,:,3)            )
        devstress(:,:,:,4) = -two*rho*( (dedI1 + I1*dedI2)*finger(:,:,:,4) - dedI2*fingersq(:,:,:,4) + I3*dedI3 )
        devstress(:,:,:,5) = -two*rho*( (dedI1 + I1*dedI2)*finger(:,:,:,5) - dedI2*fingersq(:,:,:,5)            )
        devstress(:,:,:,6) = -two*rho*( (dedI1 + I1*dedI2)*finger(:,:,:,6) - dedI2*fingersq(:,:,:,6) + I3*dedI3 )

        ! Get pressure = -tr(stress)/3
        p = -( devstress(:,:,:,1) + devstress(:,:,:,4) + devstress(:,:,:,6) ) / three

        ! Get deviatoric part of the stress
        devstress(:,:,:,1) = devstress(:,:,:,1) + p
        devstress(:,:,:,4) = devstress(:,:,:,4) + p
        devstress(:,:,:,6) = devstress(:,:,:,6) + p

        ! NEED TO CHANGE THIS. CURRENTLY WRONG
        sos_sq = real(1.528D7,rkind) + (four/three)*real(4.41D6,rkind)  ! Hardcoded to linear speed for material in IVPs in Barton's thesis

    end subroutine

    subroutine get_finger(this,g,finger,fingersq)
        class(fingereos),                                        intent(in)  :: this
        real(rkind), dimension(:,:,:,:), target,                 intent(in)  :: g
        real(rkind), dimension(size(g,1),size(g,2),size(g,3),6), intent(out) :: finger
        real(rkind), dimension(size(g,1),size(g,2),size(g,3),6), intent(out), optional :: fingersq

        if (this%usegTg) then
            associate ( g11 => g(:,:,:,1), g12 => g(:,:,:,2), g13 => g(:,:,:,3), &
                        g21 => g(:,:,:,4), g22 => g(:,:,:,5), g23 => g(:,:,:,6), &
                        g31 => g(:,:,:,7), g32 => g(:,:,:,8), g33 => g(:,:,:,9)  )
                finger(:,:,:,1) = g11
                finger(:,:,:,2) = g12
                finger(:,:,:,3) = g13
                finger(:,:,:,4) = g22
                finger(:,:,:,5) = g23
                finger(:,:,:,6) = g33
            end associate
        else
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
        end if

        if ( present(fingersq) ) then
            associate ( GG11 => finger(:,:,:,1), GG12 => finger(:,:,:,2), GG13 => finger(:,:,:,3), &
                        GG21 => finger(:,:,:,2), GG22 => finger(:,:,:,4), GG23 => finger(:,:,:,5), &
                        GG31 => finger(:,:,:,3), GG32 => finger(:,:,:,5), GG33 => finger(:,:,:,6)  )
                fingersq(:,:,:,1) = GG11*GG11 + GG12*GG21 + GG13*GG31
                fingersq(:,:,:,2) = GG11*GG12 + GG12*GG22 + GG13*GG32
                fingersq(:,:,:,3) = GG11*GG13 + GG12*GG23 + GG13*GG33
                fingersq(:,:,:,4) = GG21*GG12 + GG22*GG22 + GG23*GG32
                fingersq(:,:,:,5) = GG21*GG13 + GG22*GG23 + GG23*GG33
                fingersq(:,:,:,6) = GG31*GG13 + GG32*GG23 + GG33*GG33
            end associate
        end if

    end subroutine

    subroutine get_pT_from_energyVF(this, VF0, g0, energy, VF, p, T)
        use constants, only: zero
        class(fingereos),          intent(in)  :: this
        real(rkind), dimension(9), intent(in)  :: g0
        real(rkind),               intent(in)  :: VF0, VF, energy
        real(rkind),               intent(out) :: p, T
        
        p = zero

        T = zero
    end subroutine

    subroutine get_pT_derivatives_wrt_energyVF(this, VF0, g0, energy, VF, dpde, dpdVF, dTde, dTdVF)
        use constants, only: zero
        class(fingereos),          intent(in)  :: this
        real(rkind), dimension(9), intent(in)  :: g0
        real(rkind),               intent(in)  :: VF0, VF, energy
        real(rkind),               intent(out) :: dpde, dpdVF, dTde, dTdVF
        
        dpde = zero

        dpdVF = zero

        dTde = zero

        dTdVF = zero
    end subroutine

end module
