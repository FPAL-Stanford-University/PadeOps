module MassDiffusivityMod

    use kind_parameters, only: rkind

    implicit none

    type, abstract :: massDiffusivity

    contains

        procedure(get_diff_interface), deferred :: get_diff

    end type

    abstract interface

        ! pure subroutine get_diff_interface(this, p, T, Xs, diff)
        subroutine get_diff_interface(this, p, T, Xs, simtime, diff)
            import :: massDiffusivity
            import :: rkind
            class(massDiffusivity),  intent(in) :: this
            real(rkind), dimension(:,:,:),   intent(in)  :: p, T
            real(rkind), dimension(:,:,:,:), intent(in)  :: Xs
            real(rkind),                     intent(in)  :: simtime
            real(rkind), dimension(:,:,:,:), intent(out) :: diff
        end subroutine

    end interface

contains

end module
