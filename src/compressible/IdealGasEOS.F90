module IdealGasEOS

    use kind_parameters, only: rkind
    use constants,       only: half,one
    use EOSMod,          only: eos

    implicit none

    type, extends(eos) :: idealgas

        real(rkind) :: gam = 1.4_rkind                     ! Ratio of specific heats
        real(rkind) :: onebygam_m1 = one/(1.4_rkind-one)   ! 1/(gamma-1)
        
        real(rkind) :: Rgas = one           ! Gas constant
        real(rkind) :: onebyRgas = one      ! 1/Rgas

    contains

        ! procedure :: init
        procedure :: get_p
        procedure :: get_e_from_p
        procedure :: get_T
        procedure :: get_sos
        procedure :: get_enthalpy

    end type

    interface idealgas
        module procedure init
    end interface

contains

    function init(gam,Rgas) result(this)
        type(idealgas) :: this
        real(rkind), intent(in) :: gam, Rgas

        this%gam = gam
        this%onebygam_m1 = one/(this%gam-one)

        this%Rgas = Rgas
        this%onebyRgas = one/this%Rgas
    end function

    pure subroutine get_p(this,rho,e,p)
        class(idealgas), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,e
        real(rkind), dimension(:,:,:), intent(out) :: p

        p = (this%gam-one)*rho*e

    end subroutine

    pure subroutine get_e_from_p(this,rho,p,e)
        class(idealgas), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,p
        real(rkind), dimension(:,:,:), intent(out) :: e

        e = p * this%onebygam_m1 / rho

    end subroutine

    pure subroutine get_T(this,e,T)
        class(idealgas), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: e
        real(rkind), dimension(:,:,:), intent(out) :: T

        T = (this%gam-one)*e*this%onebyRgas

    end subroutine

    pure subroutine get_sos(this,rho,p,sos)
        class(idealgas), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,p
        real(rkind), dimension(:,:,:), intent(out) :: sos

        sos = sqrt(this%gam*p/rho)

    end subroutine

    pure subroutine get_enthalpy(this,T,h)
        class(idealgas), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: T
        real(rkind), dimension(:,:,:), intent(out) :: h

        h = this%gam * this%onebygam_m1 * this%Rgas * T

    end subroutine

end module
