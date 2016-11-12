module CauchyGreenEOSMod

    use kind_parameters, only: rkind,clen
    use constants,       only: one
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit
    use AbstractEOSMod,  only: abstracteos, get_invariants

    implicit none

    ! This is an abstract that holds all the logic for EOS that depend on invariants of the Finger tensor
    type, abstract, extends(abstracteos) :: cauchygreeneos

    contains

        procedure :: get_e_from_rho_g_T
        procedure :: get_p_devstress_T_sos2

        procedure :: get_cauchygreen

        procedure(get_energy_derivatives_interface),      deferred :: get_energy_derivatives
        procedure(get_e_from_rho_invariants_T_interface), deferred :: get_e_from_rho_invariants_T

    end type

    abstract interface

        ! Subroutine to compute the dervatives of the internal energy with respect to the invariants of the Cauchy-Green tensor and
        ! entropy (in the form of the temperature)
        subroutine get_energy_derivatives_interface(this,rho,e,I1,I2,I3,dedI1,dedI2,dedI3,T)
            import :: cauchygreeneos
            import :: rkind
            class(cauchygreeneos),                                       intent(in)  :: this
            real(rkind), dimension(:,:,:),                               intent(in)  :: rho
            real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)), intent(in)  :: e, I1, I2, I3
            real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)), intent(out) :: dedI1, dedI2, dedI3, T
        end subroutine

        ! Subroutine to get the internal energy from the temperature and material density and invariants of the Cauchy-Green tensor
        subroutine get_e_from_rho_invariants_T_interface(this,rho,I1,I2,I3,T,e)
            import :: cauchygreeneos
            import :: rkind
            class(cauchygreeneos),                                         intent(in)  :: this
            real(rkind), dimension(:,:,:),                                 intent(in)  :: rho
            real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)),   intent(in)  :: I1,I2,I3,T
            real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)),   intent(out) :: e
        end subroutine

    end interface

    interface cauchygreeneos
        module procedure init
    end interface

contains

    function init() result(this)
        class(cauchygreeneos), pointer :: this
        this%rho0 = one
        this%usegTg = .TRUE.
    end function

    subroutine get_e_from_rho_g_T(this,rho,g,T,e)
        class(cauchygreeneos),                                 intent(in)  :: this
        real(rkind), dimension(:,:,:,:),                       intent(in)  :: g
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)), intent(in)  :: rho,T
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)), intent(out) :: e
        
        real(rkind), dimension(size(g,1),size(g,2),size(g,3))   :: I1, I2, I3
        real(rkind), dimension(size(g,1),size(g,2),size(g,3),6) :: cauchygreen

        ! Get the Cauchy-Green tensor and square of the Cauchy-Green tensor
        call this%get_cauchygreen(g,cauchygreen)

        ! Get invariants of the Finger tensor
        call get_invariants(cauchygreen,I1,I2,I3)
        
        ! Get energy from density, invariants and T now (deffered to concrete class)
        call this%get_e_from_rho_invariants_T(rho,I1,I2,I3,T,e)
    end subroutine

    subroutine get_p_devstress_T_sos2(this,g,rho,e,p,T,devstress,sos_sq)
        use constants, only: zero, two, three
        class(cauchygreeneos),                                   intent(in)  :: this
        real(rkind), dimension(:,:,:,:),                         intent(in)  :: g
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)),   intent(in)  :: rho, e
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)),   intent(out) :: p, T, sos_sq
        real(rkind), dimension(size(g,1),size(g,2),size(g,3),6), intent(out) :: devstress

        real(rkind), dimension(size(g,1),size(g,2),size(g,3))   :: I1, I2, I3
        real(rkind), dimension(size(g,1),size(g,2),size(g,3))   :: dedI1, dedI2, dedI3
        real(rkind), dimension(size(g,1),size(g,2),size(g,3),6) :: cauchygreen, cauchygreensq

        ! Get the Cauchy-Green tensor and square of the Cauchy-Green tensor
        call this%get_cauchygreen(g,cauchygreen,cauchygreensq)

        ! Get invariants of the Finger tensor
        call get_invariants(cauchygreen,I1,I2,I3)

        ! Get the derivatives of the internal energy wrt invariants of the Cauchy-Green tensor
        call this%get_energy_derivatives(rho,e,I1,I2,I3,dedI1,dedI2,dedI3,T)!,sos_sq)

        ! Compute the Cauchy stress tensor
        devstress(:,:,:,1) = two*rho*( (dedI1 + I1*dedI2)*cauchygreen(:,:,:,1) - dedI2*cauchygreensq(:,:,:,1) + I3*dedI3 )
        devstress(:,:,:,2) = two*rho*( (dedI1 + I1*dedI2)*cauchygreen(:,:,:,2) - dedI2*cauchygreensq(:,:,:,2)            )
        devstress(:,:,:,3) = two*rho*( (dedI1 + I1*dedI2)*cauchygreen(:,:,:,3) - dedI2*cauchygreensq(:,:,:,3)            )
        devstress(:,:,:,4) = two*rho*( (dedI1 + I1*dedI2)*cauchygreen(:,:,:,4) - dedI2*cauchygreensq(:,:,:,4) + I3*dedI3 )
        devstress(:,:,:,5) = two*rho*( (dedI1 + I1*dedI2)*cauchygreen(:,:,:,5) - dedI2*cauchygreensq(:,:,:,5)            )
        devstress(:,:,:,6) = two*rho*( (dedI1 + I1*dedI2)*cauchygreen(:,:,:,6) - dedI2*cauchygreensq(:,:,:,6) + I3*dedI3 )

        ! Get pressure = -tr(stress)/3
        p = -( devstress(:,:,:,1) + devstress(:,:,:,4) + devstress(:,:,:,6) ) / three

        ! Get deviatoric part of the stress
        devstress(:,:,:,1) = devstress(:,:,:,1) + p
        devstress(:,:,:,4) = devstress(:,:,:,4) + p
        devstress(:,:,:,6) = devstress(:,:,:,6) + p

        ! NEED TO CHANGE THIS. CURRENTLY WRONG
        sos_sq = zero

    end subroutine

    subroutine get_cauchygreen(this,g,cauchygreen,cauchygreensq)
        class(cauchygreeneos),                                   intent(in)  :: this
        real(rkind), dimension(:,:,:,:), target,                 intent(in)  :: g
        real(rkind), dimension(size(g,1),size(g,2),size(g,3),6), intent(out) :: cauchygreen
        real(rkind), dimension(size(g,1),size(g,2),size(g,3),6), intent(out), optional :: cauchygreensq

        real(rkind), dimension(size(g,1),size(g,2),size(g,3),6) :: finger
        real(rkind), dimension(size(g,1),size(g,2),size(g,3)  ) :: invdet

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

        ! compute inverse of Finger tensor to get Cauchy-Green tensor
        invdet = one / (finger(:,:,:,1)*(finger(:,:,:,4)*finger(:,:,:,6)-finger(:,:,:,5)*finger(:,:,:,5)) &
                     -  finger(:,:,:,2)*(finger(:,:,:,2)*finger(:,:,:,6)-finger(:,:,:,5)*finger(:,:,:,3)) &
                     +  finger(:,:,:,3)*(finger(:,:,:,2)*finger(:,:,:,5)-finger(:,:,:,4)*finger(:,:,:,3)) )

        cauchygreen(:,:,:,1) = invdet * (finger(:,:,:,4) * finger(:,:,:,6) - finger(:,:,:,5) * finger(:,:,:,5))
        cauchygreen(:,:,:,2) = invdet * (finger(:,:,:,3) * finger(:,:,:,5) - finger(:,:,:,6) * finger(:,:,:,2))
        cauchygreen(:,:,:,3) = invdet * (finger(:,:,:,2) * finger(:,:,:,5) - finger(:,:,:,3) * finger(:,:,:,4))
        cauchygreen(:,:,:,4) = invdet * (finger(:,:,:,1) * finger(:,:,:,6) - finger(:,:,:,3) * finger(:,:,:,3))
        cauchygreen(:,:,:,5) = invdet * (finger(:,:,:,3) * finger(:,:,:,2) - finger(:,:,:,1) * finger(:,:,:,5))
        cauchygreen(:,:,:,6) = invdet * (finger(:,:,:,1) * finger(:,:,:,4) - finger(:,:,:,2) * finger(:,:,:,2))

        if ( present(cauchygreensq) ) then
            associate ( GG11 => cauchygreen(:,:,:,1), GG12 => cauchygreen(:,:,:,2), GG13 => cauchygreen(:,:,:,3), &
                        GG21 => cauchygreen(:,:,:,2), GG22 => cauchygreen(:,:,:,4), GG23 => cauchygreen(:,:,:,5), &
                        GG31 => cauchygreen(:,:,:,3), GG32 => cauchygreen(:,:,:,5), GG33 => cauchygreen(:,:,:,6)  )
                cauchygreensq(:,:,:,1) = GG11*GG11 + GG12*GG21 + GG13*GG31
                cauchygreensq(:,:,:,2) = GG11*GG12 + GG12*GG22 + GG13*GG32
                cauchygreensq(:,:,:,3) = GG11*GG13 + GG12*GG23 + GG13*GG33
                cauchygreensq(:,:,:,4) = GG21*GG12 + GG22*GG22 + GG23*GG32
                cauchygreensq(:,:,:,5) = GG21*GG13 + GG22*GG23 + GG23*GG33
                cauchygreensq(:,:,:,6) = GG31*GG13 + GG32*GG23 + GG33*GG33
            end associate
        end if

    end subroutine

end module
