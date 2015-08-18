module EOSMod

    use kind_parameters, only: rkind,clen
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use FiltersMod,      only: filters
    use exits,           only: GracefulExit

    type, abstract :: eos

    contains

        procedure(get_primitive_interface), deferred :: get_primitive
        procedure(get_conserved_interface), deferred :: get_conserved
        procedure(get_p_interface),         deferred :: get_p
        procedure(get_T_interface),         deferred :: get_T
        procedure(get_e_from_p_interface),  deferred :: get_e_from_p

    end type

    abstract interface

        pure subroutine get_primitive_interface(this,rho,rhou,rhov,rhow,TE,u,v,w,e,p,T)
            import :: eos
            import :: rkind
            class(eos),                    intent(in)  :: this 
            real(rkind), dimension(:,:,:), intent(in)  :: rho,rhou,rhov,rhow,TE
            real(rkind), dimension(:,:,:), intent(out) :: u,v,w,e,p,T
        end subroutine

        pure subroutine get_conserved_interface(this,rho,u,v,w,p,rhou,rhov,rhow,TE)
            import :: eos
            import :: rkind
            class(eos),                    intent(in)  :: this 
            real(rkind), dimension(:,:,:), intent(in)  :: rho,u,v,w,p
            real(rkind), dimension(:,:,:), intent(out) :: rhou,rhov,rhow,TE
        end subroutine

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

    end interface

end module
