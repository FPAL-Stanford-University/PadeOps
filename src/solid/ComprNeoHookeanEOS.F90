module ComprNeoHookeanEOSMod

    use kind_parameters,   only: rkind,clen
    use constants,         only: one
    use exits,             only: GracefulExit
    use CauchyGreenEOSMod, only: cauchygreeneos

    implicit none

    ! This class extends the cauchygreeneos class and implements the Godunov Romenskii EOS as defined in Barton (2009), Section 2.6.
    type, extends(cauchygreeneos) :: comprneohookeos

        real(rkind) :: Cv
        real(rkind) :: T0
        real(rkind) :: mu
        real(rkind) :: gam

    contains

        procedure :: get_energy_derivatives
        procedure :: get_e_from_rho_invariants_T

    end type

    interface comprneohookeos
        module procedure init
    end interface

contains

    function init(rho0,Cv,T0,mu,gam,usegTg) result(this)
        type(comprneohookeos)     :: this
        real(rkind),   intent(in) :: rho0, Cv, T0, mu, gam
        logical,       intent(in) :: usegTg

        this%rho0 = rho0
        this%Cv = Cv
        this%T0 = T0
        this%mu = mu
        this%gam = gam
        this%usegTg = usegTg

    end function

    ! Subroutine to compute the dervatives of the internal energy with respect to the invariants of the cauchygreen tensor and
    ! entropy (in the form of the temperature)
    subroutine get_energy_derivatives(this,rho,e,I1,I2,I3,dedI1,dedI2,dedI3,T)
        use constants, only: zero, half, one, two, three, four
        class(comprneohookeos),                                      intent(in)  :: this
        real(rkind), dimension(:,:,:),                               intent(in)  :: rho
        real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)), intent(in)  :: e, I1, I2, I3
        real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)), intent(out) :: dedI1, dedI2, dedI3, T
        
        real(rkind) :: invCv, muBy2Rho0, OneMGamBy2Cv, GamGm1Cv, muByRho0

        invCv = one/this%Cv
        muBy2Rho0 = half*this%mu/this%rho0
        OneMGamBy2Cv = half*(one - this%gam) * this%Cv
        GamGm1Cv = this%gam * (this%gam - one) * this%Cv
        muByRho0 = this%mu/this%rho0

        T = invCv * (e - muBy2Rho0 * I1 )

        dedI1 = muBy2Rho0
        dedI2 = zero
        dedI3 = OneMGamBy2Cv * T / I3

        !sos_sq = GamGm1Cv * T + muByRho0 * I3

    end subroutine

    ! Subroutine to get the internal energy from the temperature and material density and invariants of the cauchygreen tensor
    subroutine get_e_from_rho_invariants_T(this,rho,I1,I2,I3,T,e)
        use constants, only: half, two, three
        class(comprneohookeos),                                        intent(in)  :: this
        real(rkind), dimension(:,:,:),                                 intent(in)  :: rho
        real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)),   intent(in)  :: I1,I2,I3,T
        real(rkind), dimension(size(rho,1),size(rho,2),size(rho,3)),   intent(out) :: e

        real(rkind) :: muBy2Rho0

        muBy2Rho0 = half*this%mu/this%rho0
        e = muBy2Rho0 * I1 + this%Cv * T

    end subroutine

end module
