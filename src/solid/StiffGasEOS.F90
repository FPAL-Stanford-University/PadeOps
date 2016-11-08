module StiffGasEOS

    use kind_parameters, only: rkind
    use constants,       only: half,one
    use EOSMod,          only: eos

    implicit none

    type, extends(eos) :: stiffgas

        real(rkind) :: gam = 1.4_rkind                     ! Ratio of specific heats
        real(rkind) :: onebygam_m1 = one/(1.4_rkind-one)   ! 1/(gamma-1)
        
        real(rkind) :: Rgas = one           ! Gas constant
        real(rkind) :: onebyRgas = one      ! 1/Rgas

        real(rkind) :: PInf = one           ! Gas constant
        real(rkind) :: Cv = one           ! Gas constant

    contains

        procedure :: init
        procedure :: get_p
        procedure :: get_e_from_p
        procedure :: get_e_from_T
        procedure :: get_T
        procedure :: get_enthalpy
        procedure :: get_sos
        procedure :: get_sos2

    end type

    ! interface stiffgas
    !     module procedure init
    ! end interface

contains

    ! function init(gam_,Rgas_,PInf_) result(this)
    subroutine init(this,gam_,Rgas_,PInf_)
        class(stiffgas), intent(inout) :: this
        real(rkind),     intent(in)    :: gam_
        real(rkind),     intent(in)    :: Rgas_
        real(rkind),     intent(in)    :: PInf_

        this%gam = gam_
        this%onebygam_m1 = one/(this%gam-one)

        this%Rgas = Rgas_
        this%onebyRgas = one/this%Rgas

        this%PInf = PInf_
        this%Cv = this%Rgas/(this%gam-one)

    end subroutine
    ! end function

    pure subroutine get_p(this,rho,e,p)
        class(stiffgas), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,e
        real(rkind), dimension(:,:,:), intent(out) :: p

        p = (this%gam-one)*rho*e - this%gam*this%PInf

    end subroutine

    pure subroutine get_e_from_p(this,rho,p,e)
        class(stiffgas), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,p
        real(rkind), dimension(:,:,:), intent(out) :: e

        e = (p + this%gam*this%PInf) * this%onebygam_m1 / rho

    end subroutine

    pure subroutine get_e_from_T(this,rho,T,e)
        class(stiffgas), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,T
        real(rkind), dimension(:,:,:), intent(out) :: e

        e = this%Cv * T + this%PInf/rho

    end subroutine

    !pure subroutine get_Cv(this,Cv_)
    !    class(stiffgas), intent(in) :: this
    !    real(rkind), intent(out)  :: Cv_

    !    Cv_ = this%Rgas / (this%gam-one)

    !end subroutine

    pure subroutine get_T(this,e,T,rho)
        class(stiffgas), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: e
        real(rkind), dimension(:,:,:), intent(out) :: T
        real(rkind), dimension(:,:,:), intent(in), optional  :: rho

        T =  (e - this%PInf/rho)/this%Cv

    end subroutine

    pure subroutine get_enthalpy(this,T,enthalpy)
        class(stiffgas), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: T
        real(rkind), dimension(:,:,:), intent(out) :: enthalpy

        enthalpy = this%gam*this%Cv*T
    end subroutine

    pure subroutine get_sos(this,rho,p,sos)
        class(stiffgas), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,p
        real(rkind), dimension(:,:,:), intent(out) :: sos

        sos = sqrt(this%gam*(p+this%PInf)/rho)

    end subroutine

    pure subroutine get_sos2(this,rho,p,sos)
        class(stiffgas), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,p
        real(rkind), dimension(:,:,:), intent(out) :: sos

        sos = this%gam*(p+this%PInf)/rho

    end subroutine

end module
