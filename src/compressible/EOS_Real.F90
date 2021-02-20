module EOSMod_Real

    use kind_parameters, only: rkind,clen
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit

    type, abstract :: eos_real

    contains

        procedure(get_p_interface),         deferred :: get_p
        procedure(get_T_interface),         deferred :: get_T
        procedure(get_e_from_T_interface),  deferred :: get_e_from_T
        procedure(get_sos_interface),       deferred :: get_sos

    end type

    abstract interface

        pure subroutine get_p_interface(this,rho,T,p)
            import :: eos_real
            import :: rkind
            class(eos_real),               intent(inout)  :: this 
            real(rkind), dimension(:,:,:), intent(in)  :: rho,T
            real(rkind), dimension(:,:,:), intent(out) :: p
        end subroutine

        subroutine get_T_interface(this,rho,e,T)
            import :: eos_real
            import :: rkind
            class(eos_real),               intent(inout)  :: this
            real(rkind), dimension(:,:,:), intent(in)  :: rho,e
            real(rkind), dimension(:,:,:), intent(out) :: T
        end subroutine

        pure subroutine get_e_from_T_interface(this,rho,T,e)
            import :: eos_real
            import :: rkind
            class(eos_real),               intent(inout)  :: this
            real(rkind), dimension(:,:,:), intent(in)  :: rho,T
            real(rkind), dimension(:,:,:), intent(out) :: e
        end subroutine

        pure subroutine get_sos_interface(this,rho,T,sos)
            import :: eos_real
            import :: rkind
            class(eos_real),               intent(inout)  :: this
            real(rkind), dimension(:,:,:), intent(in)  :: rho,T
            real(rkind), dimension(:,:,:), intent(out) :: sos
        end subroutine

    end interface

contains

    ! function init() result(this)
    !     class(eos), pointer :: this
    ! end function

end module
