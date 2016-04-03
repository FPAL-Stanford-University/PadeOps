module MixtureEOSMod

    use kind_parameters, only: rkind,clen
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit
    use EOSMod,          only: eos
    use IdealGasEOS,     only: idealgas

    integer, parameter :: max_materials = 99

    type :: material_eos
        class(eos), pointer :: mat
    end type

    type, abstract :: mixture

        integer :: ns
        integer :: nxp, nyp, nzp
        type(material_eos), dimension(:), allocatable :: material
        real(rkind), dimension(:,:,:), allocatable :: gam, Cp, Cv

    contains

        procedure(init_interface),          deferred :: init
        procedure(update_interface),        deferred :: update
        procedure(get_p_interface),         deferred :: get_p
        procedure(get_T_interface),         deferred :: get_T
        procedure(get_e_from_p_interface),  deferred :: get_e_from_p
        procedure(get_sos_interface),       deferred :: get_sos
        procedure(destroy_interface),       deferred :: destroy

    end type

    abstract interface

        subroutine init_interface(this,decomp,ns,gam,Rgas,eostype)
            import :: mixture
            import :: rkind
            import :: decomp_info
            class(mixture),                                intent(inout) :: this
            type(decomp_info),                             intent(in)    :: decomp
            integer,                                       intent(in)    :: ns
            real(rkind), dimension(max_materials),         intent(in)    :: gam, Rgas
            character(len=clen), dimension(max_materials), intent(in)    :: eostype
        end subroutine

        pure subroutine update_interface(this,Ys)
            import :: mixture
            import :: rkind
            class(mixture),                                             intent(inout) :: this 
            real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns), intent(in)    :: Ys
        end subroutine

        pure subroutine get_p_interface(this,rho,e,Ys,p)
            import :: mixture
            import :: rkind
            class(mixture),                                             intent(in)  :: this 
            real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(in)  :: rho,e
            real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns), intent(in)  :: Ys
            real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(out) :: p
        end subroutine

        pure subroutine get_T_interface(this,e,Ys,T)
            import :: mixture
            import :: rkind
            class(mixture),                                             intent(in)  :: this
            real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(in)  :: e
            real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns), intent(in)  :: Ys
            real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(out) :: T
        end subroutine

        pure subroutine get_e_from_p_interface(this,rho,p,Ys,e)
            import :: mixture
            import :: rkind
            class(mixture),                                             intent(in)  :: this
            real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(in)  :: rho,p
            real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns), intent(in)  :: Ys
            real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(out) :: e
        end subroutine

        pure subroutine get_sos_interface(this,rho,p,Ys,sos)
            import :: mixture
            import :: rkind
            class(mixture),                                             intent(in)  :: this
            real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(in)  :: rho,p
            real(rkind), dimension(this%nxp,this%nyp,this%nzp,this%ns), intent(in)  :: Ys
            real(rkind), dimension(this%nxp,this%nyp,this%nzp),         intent(out) :: sos
        end subroutine

        subroutine destroy_interface(this)
            import :: mixture
            class(mixture), intent(inout)  :: this
        end subroutine

    end interface

end module
