module EOSMod

    use kind_parameters, only: rkind,clen
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit

    type, abstract :: eos

    contains

        procedure(get_p_interface),         deferred :: get_p
        procedure(get_T_interface),         deferred :: get_T
        procedure(get_e_from_p_interface),  deferred :: get_e_from_p
        procedure(get_sos_interface),       deferred :: get_sos

    end type

    abstract interface

        pure subroutine get_p_interface(this,rho,e,p)
            import :: eos
            import :: rkind
            class(eos),                    intent(in)  :: this 
            real(rkind), dimension(:,:,:), intent(in)  :: rho,e
            real(rkind), dimension(:,:,:), intent(out) :: p
        end subroutine

        pure subroutine get_T_interface(this,e,T)
            import :: eos
            import :: rkind
            class(eos),                    intent(in)  :: this
            real(rkind), dimension(:,:,:), intent(in)  :: e
            real(rkind), dimension(:,:,:), intent(out) :: T
        end subroutine

        pure subroutine get_e_from_p_interface(this,rho,p,e)
            import :: eos
            import :: rkind
            class(eos),                    intent(in)  :: this
            real(rkind), dimension(:,:,:), intent(in)  :: rho,p
            real(rkind), dimension(:,:,:), intent(out) :: e
        end subroutine

        pure subroutine get_sos_interface(this,rho,p,sos)
            import :: eos
            import :: rkind
            class(eos),                    intent(in)  :: this
            real(rkind), dimension(:,:,:), intent(in)  :: rho,p
            real(rkind), dimension(:,:,:), intent(out) :: sos
        end subroutine

    end interface

    ! interface eos
    !     module procedure init
    ! end interface

contains

    ! function init() result(this)
    !     class(eos), pointer :: this
    ! end function

end module
