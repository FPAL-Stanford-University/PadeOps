module IdealGasEOS

    use kind_parameters, only: rkind
    use constants,       only: one
    use EOSMod,          only: eos

    type, extends(eos) :: idealgas

        real(rkind) :: gam    ! Ratio of specific heats
        real(rkind) :: Rgas   ! Gas constant

    contains

        procedure :: init
        procedure :: get_p
        procedure :: get_e
        procedure :: get_T
        procedure :: get_e_from_p

    end type

contains

    pure elemental subroutine get_primitive(rho,rhou,rhov,rhow,TE,u,v,w,e)
        import :: eos
        real(rkind), intent(in)  :: rho,rhou,rhov,rhow,TE
        real(rkind), intent(out) :: u,v,w,e
        real(rkind) :: onebyrho

        onebyrho = one/rho
        u = rhou * onebyrho
        v = rhov * onebyrho
        w = rhow * onebyrho
        e = (TE*onebyrho) - ( u*u + v*v + w*w )

    end subroutine

    pure elemental subroutine get_p(rho,u,v,w,e,p)
        class(idealgas), intent(in) :: this
        real(rkind), intent(in)  :: rho,u,v,w,e
        real(rkind), intent(out) :: p

        p = (this%gam-one)*rho*e

    end subroutine

end module
