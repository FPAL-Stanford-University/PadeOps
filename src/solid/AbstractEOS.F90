module AbstractEOSMod

    use kind_parameters, only: rkind,clen
    use constants,       only: one
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit

    implicit none

    type, abstract :: abstracteos


        real(rkind) :: rho0 = one      ! Material density in undeformed configuration
        logical     :: usegTg = .TRUE. ! Use g^Tg formulation?

    contains

        procedure(get_e_from_rho_g_T_interface),              deferred :: get_e_from_rho_g_T
        procedure(get_p_devstress_T_sos2_interface),          deferred :: get_p_devstress_T_sos2
        procedure(get_pT_derivatives_wrt_energyVF_interface), deferred :: get_pT_derivatives_wrt_energyVF
        procedure(get_pT_from_energyVF_interface),            deferred :: get_pT_from_energyVF

    end type

    abstract interface

        ! Subroutine to get the internal energy from the material density, inverse deformation gradient tensor (of Finger tensor)
        ! and temperature
        subroutine get_e_from_rho_g_T_interface(this,rho,g,T,e)
            import :: abstracteos
            import :: rkind
            class(abstracteos),                                    intent(in)  :: this
            real(rkind), dimension(:,:,:,:),                       intent(in)  :: g
            real(rkind), dimension(size(g,1),size(g,2),size(g,3)), intent(in)  :: rho,T
            real(rkind), dimension(size(g,1),size(g,2),size(g,3)), intent(out) :: e
        end subroutine

        subroutine get_p_devstress_T_sos2_interface(this,g,rho,e,p,T,devstress,sos_sq)
            import :: abstracteos
            import :: rkind
            class(abstracteos),                                      intent(in)  :: this
            real(rkind), dimension(:,:,:,:),                         intent(in)  :: g
            real(rkind), dimension(size(g,1),size(g,2),size(g,3)),   intent(in)  :: rho, e
            real(rkind), dimension(size(g,1),size(g,2),size(g,3)),   intent(out) :: p, T, sos_sq
            real(rkind), dimension(size(g,1),size(g,2),size(g,3),6), intent(out) :: devstress
        end subroutine

        subroutine get_pT_derivatives_wrt_energyVF_interface(this, VF0, g0, energy, VF, dpde, dpdVF, dTde, dTdVF)
            import :: abstracteos
            import :: rkind
            class(abstracteos),        intent(in)  :: this
            real(rkind), dimension(9), intent(in)  :: g0
            real(rkind),               intent(in)  :: VF0, VF, energy
            real(rkind),               intent(out) :: dpde, dpdVF, dTde, dTdVF
        end subroutine

        subroutine get_pT_from_energyVF_interface(this, VF0, g0, energy, VF, p, T)
            import :: abstracteos
            import :: rkind
            class(abstracteos),        intent(in)  :: this
            real(rkind), dimension(9), intent(in)  :: g0
            real(rkind),               intent(in)  :: VF0, VF, energy
            real(rkind),               intent(out) :: p, T
        end subroutine

    end interface

    interface abstracteos
        module procedure init
    end interface

    interface get_invariants
        module procedure get_invariants_field, get_invariants_point
    end interface

contains

    function init() result(this)
        class(abstracteos), pointer :: this
        this%rho0 = one
        this%usegTg = .TRUE.
    end function

    subroutine get_invariants_field(tensor,I1,I2,I3)
        real(rkind), dimension(:,:,:,:), target,                              intent(in)  :: tensor
        real(rkind), dimension(size(tensor,1),size(tensor,2),size(tensor,3)), intent(out) :: I1, I2, I3

        if (size(tensor,4) == 6) then ! Symmetric tensor
            associate ( G11 => tensor(:,:,:,1), G12 => tensor(:,:,:,2), G13 => tensor(:,:,:,3), &
                        G21 => tensor(:,:,:,2), G22 => tensor(:,:,:,4), G23 => tensor(:,:,:,5), &
                        G31 => tensor(:,:,:,3), G32 => tensor(:,:,:,5), G33 => tensor(:,:,:,6)  )
                I1 = G11 + G22 + G33 
                I2 = G11*G22 + G22*G33 + G33*G11 - G12*G21 - G23*G32 - G13*G31
                I3 = G11*(G22*G33-G23*G32) - G12*(G21*G33-G31*G23) + G13*(G21*G32-G31*G22)
            end associate
        else if (size(tensor,4) == 9) then ! General tensor
            associate ( G11 => tensor(:,:,:,1), G12 => tensor(:,:,:,2), G13 => tensor(:,:,:,3), &
                        G21 => tensor(:,:,:,4), G22 => tensor(:,:,:,5), G23 => tensor(:,:,:,6), &
                        G31 => tensor(:,:,:,7), G32 => tensor(:,:,:,8), G33 => tensor(:,:,:,9)  )
                I1 = G11 + G22 + G33 
                I2 = G11*G22 + G22*G33 + G33*G11 - G12*G21 - G23*G32 - G13*G31
                I3 = G11*(G22*G33-G23*G32) - G12*(G21*G33-G31*G23) + G13*(G21*G32-G31*G22)
            end associate
        end if

    end subroutine

    subroutine get_invariants_point(tensor,I1,I2,I3)
        real(rkind), dimension(:), target, intent(in)  :: tensor
        real(rkind),                       intent(out) :: I1, I2, I3

        if (size(tensor,1) == 6) then ! Symmetric tensor
            associate ( G11 => tensor(1), G12 => tensor(2), G13 => tensor(3), &
                        G21 => tensor(2), G22 => tensor(4), G23 => tensor(5), &
                        G31 => tensor(3), G32 => tensor(5), G33 => tensor(6)  )
                I1 = G11 + G22 + G33 
                I2 = G11*G22 + G22*G33 + G33*G11 - G12*G21 - G23*G32 - G13*G31
                I3 = G11*(G22*G33-G23*G32) - G12*(G21*G33-G31*G23) + G13*(G21*G32-G31*G22)
            end associate
        else if (size(tensor,1) == 9) then ! General tensor
            associate ( G11 => tensor(1), G12 => tensor(2), G13 => tensor(3), &
                        G21 => tensor(4), G22 => tensor(5), G23 => tensor(6), &
                        G31 => tensor(7), G32 => tensor(8), G33 => tensor(9)  )
                I1 = G11 + G22 + G33 
                I2 = G11*G22 + G22*G33 + G33*G11 - G12*G21 - G23*G32 - G13*G31
                I3 = G11*(G22*G33-G23*G32) - G12*(G21*G33-G31*G23) + G13*(G21*G32-G31*G22)
            end associate
        end if

    end subroutine

end module
