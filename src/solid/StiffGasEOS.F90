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

    contains

        procedure :: init
        procedure :: get_p
        procedure :: get_e_from_p
        procedure :: get_T
        procedure :: get_sos

    end type

contains

    subroutine init(this,gam_,Rgas_,PInf_)
        class(idealgas), intent(inout) :: this
        real(rkind) :: gam_
        real(rkind) :: Rgas_
        real(rkind) :: PInf_

        this%gam = gam_
        this%onebygam_m1 = one/(this%gam-one)

        this%Rgas = Rgas_
        this%onebyRgas = one/this%Rgas

        this%PInf = PInf_

    end subroutine

    pure subroutine get_p(this,rho,e,p)
        class(idealgas), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,e
        real(rkind), dimension(:,:,:), intent(out) :: p

        p = (this%gam-one)*rho*e - this%gam*this%PInf

    end subroutine

    pure subroutine get_e_from_p(this,rho,p,e)
        class(idealgas), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: rho,p
        real(rkind), dimension(:,:,:), intent(out) :: e

        e = (p + this%gam*this%PInf) * this%onebygam_m1 / rho

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

        sos = sqrt(this%gam*(p+this%PInf)/rho)

    end subroutine

end module
