module ConstRatioBulkViscosityMod

    use kind_parameters,  only: rkind
    use constants,        only: zero
    use BulkViscosityMod, only: bulkViscosity

    implicit none

    type, extends(bulkViscosity) :: constRatioBulkViscosity

        real(rkind) :: ratio = zero

    contains

        procedure :: get_beta

    end type

    interface constRatioBulkViscosity
        module procedure init
    end interface

contains

    function init(ratio) result(this)
        type(constRatioBulkViscosity) :: this
        real(rkind),           intent(in) :: ratio

        this%ratio = ratio

    end function

    pure subroutine get_beta(this, T, mu, beta)
        class(constRatioBulkViscosity), intent(in) :: this
        real(rkind), dimension(:,:,:), intent(in)  :: T, mu
        real(rkind), dimension(:,:,:), intent(out) :: beta

        beta = this%ratio * mu

    end subroutine

end module
