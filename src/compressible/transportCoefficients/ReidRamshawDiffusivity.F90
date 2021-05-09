module ReidRamshawDiffusivityMod

    use kind_parameters,    only: rkind
    use constants,          only: zero,epssmall,one,two
    use exits,              only: GracefulExit
    use MassDiffusivityMod, only: massDiffusivity

    implicit none

    type, extends(massDiffusivity) :: reidRamshawDiffusivity

        ! Compute constant diffusivity for all species
        ! based on Ramshaw's effective binary diffusion model
        ! and Reid's model for the mass diffusion of a binary mixture

        integer     :: ns

        real(rkind) :: const                                   ! Constant in compatible units
        real(rkind) :: A, B, C, D, E, F, G, H                  ! Constants in collision integral
        real(rkind), dimension(:), allocatable :: molwt        ! Molecular weights
        real(rkind), dimension(:), allocatable :: sigma        ! Collision diameter
        real(rkind), dimension(:), allocatable :: eps_by_k     ! Lennard-Jones energy parameter

    contains

        procedure :: reidBinaryDiffusion
        procedure :: get_diff
        final     :: destroy

    end type

    interface reidRamshawDiffusivity
        module procedure init
    end interface

contains

    function init(ns, const, molwt, sigma, eps_by_k, A, B, C, D, E, F, G, H) result(this)
        type(reidRamshawDiffusivity)           :: this
        integer,                    intent(in) :: ns
        real(rkind),                intent(in) :: const, A, B, C, D, E, F, G, H
        real(rkind), dimension(ns), intent(in) :: molwt, sigma, eps_by_k

        if (ns < 1) call GracefulExit("Cannot have less than 1 species", 4529)
        this%ns = ns

        if (allocated(this%molwt))    deallocate(this%molwt);    allocate(this%molwt(this%ns))
        if (allocated(this%sigma))    deallocate(this%sigma);    allocate(this%sigma(this%ns))
        if (allocated(this%eps_by_k)) deallocate(this%eps_by_k); allocate(this%eps_by_k(this%ns))

        this%const    = const
        this%molwt    = molwt
        this%sigma    = sigma
        this%eps_by_k = eps_by_k

        this%A = A
        this%B = B
        this%C = C
        this%D = D
        this%E = E
        this%F = F
        this%G = G
        this%H = H

    end function

    subroutine destroy(this)
        type(reidRamshawDiffusivity), intent(inout) :: this

        if (allocated(this%molwt))    deallocate(this%molwt)
        if (allocated(this%sigma))    deallocate(this%sigma)
        if (allocated(this%eps_by_k)) deallocate(this%eps_by_k)

    end subroutine

    pure subroutine reidBinaryDiffusion(this, p, T, diff, i, j)
        class(reidRamshawDiffusivity),   intent(in)  :: this
        real(rkind), dimension(:,:,:),   intent(in)  :: p, T
        integer,                         intent(in)  :: i, j
        real(rkind), dimension(:,:,:),   intent(out) :: diff

        real(rkind), dimension(size(T,1),size(T,2),size(T,3)) :: T_star, Omega_ij
        real(rkind) :: sigma_ij, M_ij

        T_star = T / sqrt(this%eps_by_k(i) * this%eps_by_k(j))
        Omega_ij = this%A * T_star**(this%B) + this%C * exp(this%D * T_star) &
                 + this%E * exp(this%F * T_star) + this%G * exp(this%H * T_star)

        M_ij = two / (one/this%molwt(i) + one/this%molwt(j))
        sigma_ij = (this%sigma(i) + this%sigma(j)) / two

        diff = (this%const/(sqrt(M_ij)*sigma_ij*sigma_ij)) * (sqrt(T))**3 / (Omega_ij * p)
        
    end subroutine

    !pure subroutine get_diff(this, p, T, Xs, diff)
    subroutine get_diff(this, p, T, Xs, simtime, diff)
        class(reidRamshawDiffusivity),   intent(in)  :: this
        real(rkind), dimension(:,:,:),   intent(in)  :: p, T
        real(rkind), dimension(:,:,:,:), intent(in)  :: Xs
        real(rkind),                     intent(in)  :: simtime
        real(rkind), dimension(:,:,:,:), intent(out) :: diff

        real(rkind), dimension(size(T,1),size(T,2),size(T,3)) :: diff_ij
        integer :: i, j

        diff = zero
        do i = 1,this%ns
            ! do j = 1,this%ns
            do j = 1,i-1
                ! if (j == i) cycle
                call this%reidBinaryDiffusion(p, T, diff_ij, i, j)
                diff(:,:,:,i) = diff(:,:,:,i) + Xs(:,:,:,j) / diff_ij
                diff(:,:,:,j) = diff(:,:,:,j) + Xs(:,:,:,i) / diff_ij
            end do
        end do

        diff = (one - Xs) / (diff + epssmall)

    end subroutine

end module
