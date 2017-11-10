module BulkViscosityMod

    use kind_parameters, only: rkind

    implicit none

    type, abstract :: bulkViscosity

    contains

        procedure(get_beta_interface), deferred :: get_beta

    end type

    abstract interface

        pure subroutine get_beta_interface(this, T, mu, beta)
            import :: bulkViscosity
            import :: rkind
            class(bulkViscosity),           intent(in) :: this
            real(rkind), dimension(:,:,:), intent(in)  :: T, mu
            real(rkind), dimension(:,:,:), intent(out) :: beta
        end subroutine

    end interface

contains

end module
