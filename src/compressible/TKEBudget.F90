module TKEBudgetMod

    use mpi
    use decomp_2d,       only: decomp_info, nrank
    use kind_parameters, only: rkind, mpirkind
    use constants,       only: half,one
    use exits,           only: GracefulExit
    use AveragingMod,    only: averaging

    implicit none

    type :: tkeBudget
        private

        type(decomp_info), pointer :: gp
        type(averaging) :: avg
        logical, dimension(3) :: averaging_directions = [.true., .true., .true.]

    contains

        procedure          :: reynolds_avg
        procedure          :: reynolds_avg_and_fluct
        procedure          :: favre_avg
        procedure          :: favre_avg_and_fluct
        final              :: destroy

    end type

    interface tkeBudget
        module procedure init
    end interface

contains

    function init(gp, pencil, averaging_directions) result(this)
        type(tkeBudget)                        :: this
        class(decomp_info), target, intent(in) :: gp
        integer,                    intent(in) :: pencil
        logical, dimension(3),      intent(in) :: averaging_directions

        integer, dimension(3) :: st
        integer :: ierr

        this%gp => gp
        this%averaging_directions = averaging_directions

        this%avg = averaging(this%gp, pencil, this%averaging_directions)

    end function

    impure elemental subroutine destroy(this)
        type(tkeBudget), intent(inout) :: this

        nullify(this%gp)

    end subroutine

    subroutine reynolds_avg(this, f, f_bar)
        class(tkeBudget),                                                                       intent(in)  :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),       intent(in)  :: f
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)), intent(out) :: f_bar

        call this%avg%get_average(f, f_bar)

    end subroutine

    subroutine reynolds_avg_and_fluct(this, f, f_bar, f_prime)
        class(tkeBudget),                                                                       intent(in)  :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),       intent(in)  :: f
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)), intent(out) :: f_bar
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),       intent(out) :: f_prime

        call this%reynolds_avg(f, f_bar)
        call this%avg%get_fluctuations(f, f_bar, f_prime)

    end subroutine

    subroutine favre_avg(this, rho, f, f_tilde)
        class(tkeBudget),                                                                       intent(in)  :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),       intent(in)  :: rho, f
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)), intent(out) :: f_tilde

        this%avg%get_weighted_average(rho, f, f_tilde)

    end subroutine

    subroutine favre_avg_and_fluct(this, rho, f, f_tilde, f_pprime)
        class(tkeBudget),                                                                       intent(in)  :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),       intent(in)  :: rho, f
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)), intent(out) :: f_tilde
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),       intent(out) :: f_pprime

        call this%favre_avg(rho, f, f_tilde)
        call this%avg%get_fluctuations(f, f_tilde, f_pprime)

    end subroutine

end module
