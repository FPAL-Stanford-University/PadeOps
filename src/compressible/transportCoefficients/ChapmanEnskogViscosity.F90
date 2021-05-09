module ChapmanEnskogViscosityMod

    use kind_parameters,   only: rkind
    use constants,         only: zero,one
    use ShearViscosityMod, only: shearViscosity

    implicit none

    type, extends(shearViscosity) :: chapmanEnskogViscosity

        real(rkind) :: const
        real(rkind) :: molwt
        real(rkind) :: eps_by_k
        real(rkind) :: sigma2
        real(rkind) :: A, B, C, D, E, F

    contains

        procedure :: get_mu

    end type

    interface chapmanEnskogViscosity
        module procedure init
    end interface

contains

    function init(const, molwt, eps_by_k, sigma, A, B, C, D, E, F) result(this)
        type(chapmanEnskogViscosity) :: this
        real(rkind),      intent(in) :: const, molwt, eps_by_k, sigma, A, B, C, D, E, F

        this%const    = const
        this%molwt    = molwt
        this%eps_by_k = eps_by_k
        this%sigma2   = sigma*sigma
        this%const    = const

        this%A = A
        this%B = B
        this%C = C
        this%D = D
        this%E = E
        this%F = F

    end function

    pure subroutine get_mu(this, T, simtime, mu)
        class(chapmanEnskogViscosity), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: T
        real(rkind),                   intent(in)  :: simtime
        real(rkind), dimension(:,:,:), intent(out) :: mu

        real(rkind), dimension(size(T,1), size(T,2), size(T,3)) :: T_star, Omega

        T_star = T / this%eps_by_k
        Omega = this%A * (T_star)**this%B + this%C * exp(this%D * T_star) + this%E * exp(this%F * T_star)

        mu = this%const * sqrt( T * this%molwt ) / (this%sigma2 * Omega)

    end subroutine

end module
