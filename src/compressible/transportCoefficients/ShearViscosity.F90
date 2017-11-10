module ShearViscosityMod

    use kind_parameters, only: rkind

    implicit none

    type, abstract :: shearViscosity

    contains

        procedure(get_mu_interface), deferred :: get_mu

    end type

    abstract interface

        pure subroutine get_mu_interface(this, T, mu)
            import :: shearViscosity
            import :: rkind
            class(shearViscosity),         intent(in)  :: this
            real(rkind), dimension(:,:,:), intent(in)  :: T
            real(rkind), dimension(:,:,:), intent(out) :: mu
        end subroutine

    end interface

contains

end module
