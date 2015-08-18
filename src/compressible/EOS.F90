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

        pure elemental subroutine get_primitive(rho,rhou,rhov,rhow,TE,u,v,w,e)
            import :: eos
            real(rkind), intent(in)  :: rho,rhou,rhov,rhow,TE
            real(rkind), intent(out) :: u,v,w,e
        end subroutine

        pure elemental subroutine get_conserved(rho,u,v,w,e,rhou,rhov,rhow,TE)
            import :: eos
            real(rkind), intent(in)   :: rho,u,v,w,e
            real(rkind), intent(out)  :: rhou,rhov,rhow,TE
        end subroutine

        pure elemental subroutine get_p(rho,u,v,w,e,p)
            import :: eos
            real(rkind), intent(in)  :: rho,u,v,w,e
            real(rkind), intent(out) :: p
        end subroutine

        pure elemental subroutine get_T(rho,u,v,w,e,T)
            import :: eos
            real(rkind), intent(in)  :: rho,u,v,w,e
            real(rkind), intent(out) :: T
        end subroutine

        pure elemental subroutine get_e_from_p(rho,u,v,w,p,e)
            import :: eos
            real(rkind), intent(in)  :: rho,u,v,w,p
            real(rkind), intent(out) :: e
        end subroutine

    end interface

end module
