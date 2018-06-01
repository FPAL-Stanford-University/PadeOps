module TKEBudgetMod

    use mpi
    use decomp_2d,        only: decomp_info, nrank
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero,half,one
    use exits,            only: GracefulExit
    use AveragingMod,     only: averaging
    use io_hdf5_stuff,    only: io_hdf5
    use DerivativesMod,   only: derivatives
    use mytranspose2DMod, only: mytranspose2D
    use operators,        only: gradient

    implicit none

    type :: tkeBudget
        ! private

        type(decomp_info), pointer :: gp
        type(averaging)            :: avg
        type(derivatives), pointer :: der
        
        type(derivatives)   :: der_avg
        type(mytranspose2D) :: gp_avg !!! Only support 2D for now

        type(io_hdf5)       :: viz

        logical, dimension(3) :: averaging_directions = [.true., .true., .true.]
        integer, dimension(2) :: x_bc, y_bc, z_bc

    contains

        procedure          :: init
        procedure          :: reynolds_avg
        procedure          :: reynolds_avg_and_fluct
        procedure          :: favre_avg
        procedure          :: favre_avg_and_fluct
        procedure          :: get_tke
        procedure          :: get_reynolds_stress
        procedure          :: get_duidxj_avg
        procedure          :: get_gradp_avg
        procedure          :: get_div_tau_avg
        procedure          :: get_production
        procedure          :: get_p_dil
        procedure          :: get_dissipation
        procedure          :: tke_budget
        final              :: destructor

    end type

    ! interface tkeBudget
    !     module procedure init
    ! end interface

contains

    ! function init(gp, der, mesh, dx, dy, dz, averaging_directions, outputdir, x_bc, y_bc, z_bc, reduce_precision) result(this)
    subroutine init(this, gp, der, mesh, dx, dy, dz, averaging_directions, outputdir, x_bc, y_bc, z_bc, reduce_precision)
        ! type(tkeBudget)                             :: this
        class(tkeBudget)                            :: this
        class(decomp_info), target,      intent(in) :: gp
        class(derivatives), target,      intent(in) :: der
        real(rkind), dimension(:,:,:,:), intent(in) :: mesh
        real(rkind),                     intent(in) :: dx, dy, dz
        logical, dimension(3),           intent(in) :: averaging_directions
        character(len=clen),             intent(in) :: outputdir
        integer, dimension(2),           intent(in) :: x_bc, y_bc, z_bc
        logical, optional,               intent(in) :: reduce_precision

        logical :: reduce_precision_
        integer, dimension(3) :: st
        integer :: ierr

        this%gp => gp
        this%der => der
        this%averaging_directions = averaging_directions

        this%x_bc = x_bc
        this%y_bc = y_bc
        this%z_bc = z_bc

        reduce_precision_ = .true.
        if (present(reduce_precision)) reduce_precision_ = reduce_precision

        ! Initialize averaging object in the y  pencil
        this%avg = averaging(this%gp, 2, this%averaging_directions)

        if ( (this%avg%avg_dim /= 0) .and. (this%avg%avg_dim /= 2) ) then
            call GracefulExit("Only all direction or 1 direction averaging is supported for now.", 4589)
        end if

        if (this%avg%avg_dim == 2) then
            call this%der_avg%init( this%gp, &
                                    dx, dy, dz, &
                                    averaging_directions(1), averaging_directions(2), averaging_directions(3), &
                                    der%getMethodx(), der%getMethody(), der%getMethodz() )

            if (averaging_directions(1)) then
                ! call this%der_avg%set_ysz( [1, this%gp%ysz(2), this%gp%ysz(3)] )                 ! Set ysz to make arrays 2D
                ! call this%der_avg%set_zsz( [1, this%gp%zsz(2), this%gp%zsz(3)] )                 ! Set zsz to make arrays 2D

                call GracefulExit("Averaging in X is supported for budgets! :(", 4589)
            end if

            if (averaging_directions(2)) then
                call GracefulExit("Averaging in Y is supported for budgets! :(", 4589)
            end if

            if (averaging_directions(3)) then
                call this%der_avg%set_xsz( [this%gp%xsz(1), this%gp%xsz(2), 1] )                 ! Set xsz to make arrays 2D
                call this%der_avg%set_ysz( [this%gp%ysz(1), this%gp%ysz(2), 1] )                 ! Set ysz to make arrays 2D

                call this%gp_avg%init(this%gp%xsz(1),this%gp%ysz(2),this%avg%xy_comm)

                ! Initialize HDF5 output object (only z index of 1 since we're outputting a 1D field)
                call this%viz%init(mpi_comm_world, this%gp, 'y', outputdir, 'TKEBudget', reduce_precision=reduce_precision_, &
                                   read_only=.false., subdomain_lo=[1,1,1], subdomain_hi=[this%gp%xsz(1),this%gp%ysz(2),1], &
                                   jump_to_last=.true.)

                ! Write the coordinates of subdomain out
                call this%viz%write_coords(mesh)
            end if

        end if

    ! end function
    end subroutine

    ! impure elemental subroutine destructor(this)
    subroutine destructor(this)
        type(tkeBudget), intent(inout) :: this

        nullify(this%gp)
        nullify(this%der)

        call this%viz%destroy()

    end subroutine

    subroutine reynolds_avg(this, f, f_bar)
        class(tkeBudget),                                                                       intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),       intent(in)    :: f
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)), intent(out)   :: f_bar

        call this%avg%get_average(f, f_bar)

    end subroutine

    subroutine reynolds_avg_and_fluct(this, f, f_bar, f_prime)
        class(tkeBudget),                                                                       intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),       intent(in)    :: f
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)), intent(out)   :: f_bar
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),       intent(out)   :: f_prime

        call this%reynolds_avg(f, f_bar)
        call this%avg%get_fluctuations(f, f_bar, f_prime)

    end subroutine

    subroutine favre_avg(this, rho, f, f_tilde)
        class(tkeBudget),                                                                       intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),       intent(in)    :: rho, f
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)), intent(out)   :: f_tilde

        call this%avg%get_weighted_average(rho, f, f_tilde)

    end subroutine

    subroutine favre_avg_and_fluct(this, rho, f, f_tilde, f_pprime)
        class(tkeBudget),                                                                       intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),       intent(in)    :: rho, f
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)), intent(out)   :: f_tilde
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),       intent(out)   :: f_pprime

        call this%favre_avg(rho, f, f_tilde)
        call this%avg%get_fluctuations(f, f_tilde, f_pprime)

    end subroutine

    subroutine get_tke(this, rho, u, v, w, tke)
        class(tkeBudget),                                                                       intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),       intent(in)    :: rho, u, v, w
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)), intent(out)   :: tke

        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3))       :: tke3d, tmp
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: dum, rho_bar

        call this%reynolds_avg(rho, rho_bar)

        call this%favre_avg_and_fluct(rho, u, dum, tmp) ! tmp has u_pprime
        tke3d = tmp*tmp ! u_pprime^2

        call this%favre_avg_and_fluct(rho, v, dum, tmp) ! tmp has v_pprime
        tke3d = tke3d + tmp*tmp ! Add v_pprime^2

        call this%favre_avg_and_fluct(rho, w, dum, tmp) ! tmp has w_pprime
        tke3d = tke3d + tmp*tmp ! Add w_pprime^2

        ! Get tke as favre averaged tke3d
        call this%favre_avg(rho, tke3d, tke)

        tke = half*rho_bar*tke

    end subroutine

    subroutine get_reynolds_stress(this, rho, u_pprime, v_pprime, w_pprime, Rij)
        class(tkeBudget),                                                                          intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),          intent(in)    :: rho
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),          intent(in)    :: u_pprime, &
                                                                                                                    v_pprime, &
                                                                                                                    w_pprime
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3), 6), intent(out)   :: Rij

        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3)) :: tmp

        ! R11
        tmp = u_pprime*u_pprime
        call this%favre_avg(rho, tmp, Rij(:,:,:,1))

        ! R12
        tmp = u_pprime*v_pprime
        call this%favre_avg(rho, tmp, Rij(:,:,:,2))

        ! R13
        tmp = u_pprime*w_pprime
        call this%favre_avg(rho, tmp, Rij(:,:,:,3))

        ! R22
        tmp = v_pprime*v_pprime
        call this%favre_avg(rho, tmp, Rij(:,:,:,4))

        ! R23
        tmp = v_pprime*w_pprime
        call this%favre_avg(rho, tmp, Rij(:,:,:,5))

        ! R33
        tmp = w_pprime*w_pprime
        call this%favre_avg(rho, tmp, Rij(:,:,:,6))

    end subroutine

    subroutine get_duidxj_avg(this, u, v, w, duidxj_avg)
        class(tkeBudget),                                                                                  intent(in)  :: this
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),            intent(in)  :: u, v, w
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3), 9), target, intent(out) :: duidxj_avg

        real(rkind), dimension(this%gp_avg%nx,this%gp_avg%ay,1) :: xbuf1, xbuf2 !!! HACK !!! Only works for Z averaging
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        if (this%avg%avg_dim == 0) then
            duidxj_avg = zero

        else if (this%avg%avg_dim == 2) then
            dudx => duidxj_avg(:,:,:,1); dudy => duidxj_avg(:,:,:,2); dudz => duidxj_avg(:,:,:,3);
            dvdx => duidxj_avg(:,:,:,4); dvdy => duidxj_avg(:,:,:,5); dvdz => duidxj_avg(:,:,:,6);
            dwdx => duidxj_avg(:,:,:,7); dwdy => duidxj_avg(:,:,:,8); dwdz => duidxj_avg(:,:,:,9);

            ! Set Z derivatives to zero
            dudz = zero; dvdz = zero; dwdz = zero

            ! dudx
            call this%gp_avg%transpose_y_to_x(u(:,:,1),xbuf1(:,:,1))
            call this%der_avg%ddx(xbuf1,xbuf2,-this%x_bc(1),-this%x_bc(2))
            call this%gp_avg%transpose_x_to_y(xbuf2(:,:,1),dudx(:,:,1))
            
            ! dvdx
            call this%gp_avg%transpose_y_to_x(v(:,:,1),xbuf1(:,:,1))
            call this%der_avg%ddx(xbuf1,xbuf2, this%x_bc(1), this%x_bc(2))
            call this%gp_avg%transpose_x_to_y(xbuf2(:,:,1),dvdx(:,:,1))
            
            ! dwdx
            call this%gp_avg%transpose_y_to_x(w(:,:,1),xbuf1(:,:,1))
            call this%der_avg%ddx(xbuf1,xbuf2, this%x_bc(1), this%x_bc(2))
            call this%gp_avg%transpose_x_to_y(xbuf2(:,:,1),dwdx(:,:,1))

            call this%der_avg%ddy(u,dudy, this%y_bc(1), this%y_bc(2))
            call this%der_avg%ddy(v,dvdy,-this%y_bc(1),-this%y_bc(2))
            call this%der_avg%ddy(w,dwdy, this%y_bc(1), this%y_bc(2))
        end if
        
    end subroutine

    subroutine get_gradp_avg(this,p_avg,gradp_avg)
        class(tkeBudget),                                                                                  intent(in)  :: this
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),            intent(in)  :: p_avg
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3), 3), target, intent(out) :: gradp_avg

        real(rkind), dimension(this%gp_avg%nx,this%gp_avg%ay,1) :: xbuf1, xbuf2 !!! HACK !!! Only works for Z averaging
        real(rkind), dimension(:,:,:), pointer :: dpdx, dpdy, dpdz

        dpdx => gradp_avg(:,:,:,1); dpdy => gradp_avg(:,:,:,2); dpdz => gradp_avg(:,:,:,3);

        ! Set Z derivatives to zero
        dpdz = zero

        ! dpdx
        call this%gp_avg%transpose_y_to_x(p_avg(:,:,1),xbuf1(:,:,1))
        call this%der_avg%ddx(xbuf1,xbuf2, this%x_bc(1), this%x_bc(2))
        call this%gp_avg%transpose_x_to_y(xbuf2(:,:,1),dpdx(:,:,1))
        
        call this%der_avg%ddy(p_avg,dpdy, this%y_bc(1), this%y_bc(2))
        
    end subroutine

    subroutine get_div_tau_avg(this,tau_avg,tau_avg_div)
        class(tkeBudget),                                                                                 intent(in)  :: this
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),6), target, intent(in)  :: tau_avg
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),3), target, intent(out) :: tau_avg_div

        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: tmp
        real(rkind), dimension(this%gp_avg%nx,this%gp_avg%ay,1) :: xbuf1, xbuf2 !!! HACK !!! Only works for Z averaging
        real(rkind), dimension(:,:,:), pointer :: tauxx, tauxy, tauxz, tauyy, tauyz, tauzz, divx, divy, divz

        tauxx => tau_avg(:,:,:,1); tauxy => tau_avg(:,:,:,2); tauxz => tau_avg(:,:,:,3);
                                   tauyy => tau_avg(:,:,:,4); tauyz => tau_avg(:,:,:,5);
                                                              tauzz => tau_avg(:,:,:,6);

        divx => tau_avg_div(:,:,:,1); divy => tau_avg_div(:,:,:,2); divz => tau_avg_div(:,:,:,3);

        ! \frac{ \partial \tau_{1j} }{ \partial x_j }
        call this%gp_avg%transpose_y_to_x(tauxx(:,:,1),xbuf1(:,:,1))
        call this%der_avg%ddx(xbuf1,xbuf2, this%x_bc(1), this%x_bc(2))
        call this%gp_avg%transpose_x_to_y(xbuf2(:,:,1),divx(:,:,1))
        
        call this%der_avg%ddy(tauxy(:,:,1),tmp(:,:,1),-this%y_bc(1),-this%y_bc(2))

        divx = divx + tmp ! Z derivatives are zero
        
        ! \frac{ \partial \tau_{2j} }{ \partial x_j }
        call this%gp_avg%transpose_y_to_x(tauxy(:,:,1),xbuf1(:,:,1))
        call this%der_avg%ddx(xbuf1,xbuf2,-this%x_bc(1),-this%x_bc(2))
        call this%gp_avg%transpose_x_to_y(xbuf2(:,:,1),divy(:,:,1))
        
        call this%der_avg%ddy(tauyy(:,:,1),tmp(:,:,1), this%y_bc(1), this%y_bc(2))

        divy = divy + tmp ! Z derivatives are zero
        
        ! \frac{ \partial \tau_{3j} }{ \partial x_j }
        call this%gp_avg%transpose_y_to_x(tauxz(:,:,1),xbuf1(:,:,1))
        call this%der_avg%ddx(xbuf1,xbuf2,-this%x_bc(1),-this%x_bc(2))
        call this%gp_avg%transpose_x_to_y(xbuf2(:,:,1),divz(:,:,1))
        
        call this%der_avg%ddy(tauyz(:,:,1),tmp(:,:,1),-this%y_bc(1),-this%y_bc(2))

        divz = divz + tmp ! Z derivatives are zero
        
    end subroutine

    subroutine get_production(this, rho_bar, Rij, grad_u_tilde, production)
        class(tkeBudget),                                                                                 intent(in)  :: this
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),           intent(in)  :: rho_bar
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),6), target, intent(in)  :: Rij
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),9), target, intent(in)  :: grad_u_tilde
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),           intent(out) :: production
        
        real(rkind), dimension(:,:,:), pointer :: R11, R12, R13, R22, R23, R33
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz

        R11 => Rij(:,:,:,1); R12 => Rij(:,:,:,2); R13 => Rij(:,:,:,3);
                             R22 => Rij(:,:,:,4); R23 => Rij(:,:,:,5);
                                                  R33 => Rij(:,:,:,6);

        dudx => grad_u_tilde(:,:,:,1); dudy => grad_u_tilde(:,:,:,2); dudz => grad_u_tilde(:,:,:,3);
        dvdx => grad_u_tilde(:,:,:,4); dvdy => grad_u_tilde(:,:,:,5); dvdz => grad_u_tilde(:,:,:,6);
        dwdx => grad_u_tilde(:,:,:,7); dwdy => grad_u_tilde(:,:,:,8); dwdz => grad_u_tilde(:,:,:,9);

        production = R11*dudx + R12*dudy + R13*dudz &
                   + R12*dvdx + R22*dvdy + R23*dvdz &
                   + R13*dwdx + R23*dwdy + R33*dwdz

        production = -rho_bar * production
    end subroutine
    
    subroutine get_p_dil(this, p, p_prime, u_pprime_bar, grad_p_bar, grad_u_pprime, p_dil_fluct, baropycnal, fluct_p_dil)
        class(tkeBudget),                                                                                 intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),                 intent(in)    :: p, p_prime
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),3),         intent(in)    :: u_pprime_bar
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),3),         intent(in)    :: grad_p_bar
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3),      9), target, intent(in)    :: grad_u_pprime
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),           intent(out)   :: p_dil_fluct
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),           intent(out)   :: baropycnal
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),           intent(out)   :: fluct_p_dil

        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3)) :: tmp
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz

        dudx => grad_u_pprime(:,:,:,1); dudy => grad_u_pprime(:,:,:,2); dudz => grad_u_pprime(:,:,:,3);
        dvdx => grad_u_pprime(:,:,:,4); dvdy => grad_u_pprime(:,:,:,5); dvdz => grad_u_pprime(:,:,:,6);
        dwdx => grad_u_pprime(:,:,:,7); dwdy => grad_u_pprime(:,:,:,8); dwdz => grad_u_pprime(:,:,:,9);

        ! Pressure fluctuation-dilatation correlation \overline{ p \frac{\partial u_i^{\prime \prime}}{\partial x_i} }
        tmp = p * (dudx + dvdy + dwdz)
        call this%reynolds_avg(tmp, p_dil_fluct)

        ! Fluctuating pressure fluctuation-dilatation correlation \overline{ p^\prime \frac{\partial u_i^{\prime \prime}}{\partial x_i} }
        tmp = p_prime*(dudx + dvdy + dwdz)
        call this%reynolds_avg(tmp, fluct_p_dil)
        
        ! Baropycnal term -\overline(u_i^{\prime\prime}) \frac{\partial \overline{p}}{\partial x_i}
        baropycnal = - ( u_pprime_bar(:,:,:,1)*grad_p_bar(:,:,:,1) &
                       + u_pprime_bar(:,:,:,2)*grad_p_bar(:,:,:,2) &
                       + u_pprime_bar(:,:,:,3)*grad_p_bar(:,:,:,3) )
    end subroutine

    subroutine get_dissipation(this, tauij, tau_bar, tau_prime, u_pprime_bar, grad_u_pprime, dissipation, diss_mass_flux, diss_fluct)
        class(tkeBudget),                                                                                 intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3),      6), target, intent(in)    :: tauij, tau_prime
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),6), target, intent(in)    :: tau_bar
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),3),         intent(in)    :: u_pprime_bar
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3),      9), target, intent(in)    :: grad_u_pprime
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),           intent(out)   :: dissipation
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),           intent(out)   :: diss_mass_flux
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),           intent(out)   :: diss_fluct

        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3)) :: tmp
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),3) :: tau_bar_div
        real(rkind), dimension(:,:,:), pointer :: tauxx, tauxy, tauxz, tauyy, tauyz, tauzz
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz

        tauxx => tauij(:,:,:,1); tauxy => tauij(:,:,:,2); tauxz => tauij(:,:,:,3);
                                 tauyy => tauij(:,:,:,4); tauyz => tauij(:,:,:,5);
                                                          tauzz => tauij(:,:,:,6);

        dudx => grad_u_pprime(:,:,:,1); dudy => grad_u_pprime(:,:,:,2); dudz => grad_u_pprime(:,:,:,3);
        dvdx => grad_u_pprime(:,:,:,4); dvdy => grad_u_pprime(:,:,:,5); dvdz => grad_u_pprime(:,:,:,6);
        dwdx => grad_u_pprime(:,:,:,7); dwdy => grad_u_pprime(:,:,:,8); dwdz => grad_u_pprime(:,:,:,9);

        tmp = tauxx*dudx + tauxy*dudy + tauxz*dudz &
            + tauxy*dvdx + tauyy*dvdy + tauyz*dvdz &
            + tauxz*dwdx + tauyz*dwdy + tauzz*dwdz
        call this%reynolds_avg(tmp, dissipation)

        tauxx => tau_bar(:,:,:,1); tauxy => tau_bar(:,:,:,2); tauxz => tau_bar(:,:,:,3);
                                   tauyy => tau_bar(:,:,:,4); tauyz => tau_bar(:,:,:,5);
                                                              tauzz => tau_bar(:,:,:,6);

        ! Need to get tau_bar divergence
        call this%get_div_tau_avg(tau_bar,tau_bar_div)

        ! Mass flux term -\overline(u_i^{\prime\prime}) \frac{\partial \overline{\tau_{ij}}}{\partial x_j}
        diss_mass_flux = ( u_pprime_bar(:,:,:,1)*tau_bar_div(:,:,:,1) &
                         + u_pprime_bar(:,:,:,2)*tau_bar_div(:,:,:,2) &
                         + u_pprime_bar(:,:,:,3)*tau_bar_div(:,:,:,3) )

        tauxx => tau_prime(:,:,:,1); tauxy => tau_prime(:,:,:,2); tauxz => tau_prime(:,:,:,3);
                                     tauyy => tau_prime(:,:,:,4); tauyz => tau_prime(:,:,:,5);
                                                                  tauzz => tau_prime(:,:,:,6);

        ! Fluctuating stress dissipation
        tmp = tauxx*dudx + tauxy*dudy + tauxz*dudz &
            + tauxy*dvdx + tauyy*dvdy + tauyz*dvdz &
            + tauxz*dwdx + tauyz*dwdy + tauzz*dwdz
        call this%reynolds_avg(tmp, diss_fluct)

    end subroutine

    subroutine tke_budget(this, rho, u, v, w, p, tauij, tke_old, tke_prefilter, tsim, dt)
        class(tkeBudget),                                                                         intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),         intent(in)    :: rho, u, v, w, p
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3),      6), intent(in)    :: tauij
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),   intent(in)    :: tke_old, tke_prefilter
        real(rkind),                                                                              intent(in)    :: tsim, dt

        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: ddt_tke
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: production
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: p_dil_fluct
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: fluct_p_dil
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: baropycnal
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: dissipation
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: diss_mass_flux
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: diss_fluct
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: dissipation_num

        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3))   :: rho_bar, p_bar
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3))   :: u_tilde, v_tilde, w_tilde
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3))   :: tke
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),6) :: Rij
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),6) :: tau_bar
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),9) :: grad_u_tilde
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),3) :: grad_p_bar
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),3) :: u_pprime_bar

        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3))   :: u_pprime, v_pprime, w_pprime, p_prime
        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3),9) :: grad_u_pprime
        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3),6) :: tau_prime
        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3))   :: tmp

        integer :: i

        ! Get density average
        call this%reynolds_avg(rho, rho_bar)

        ! Get velocity favre average and fluctuations
        call this%favre_avg_and_fluct(rho, u, u_tilde, u_pprime)
        call this%favre_avg_and_fluct(rho, v, v_tilde, v_pprime)
        call this%favre_avg_and_fluct(rho, w, w_tilde, w_pprime)

        ! Get pressure average and fluctuations
        call this%reynolds_avg_and_fluct(p, p_bar, p_prime)

        ! Get tke
        tmp = half*rho*( u_pprime*u_pprime + v_pprime*v_pprime + w_pprime*w_pprime )
        call this%reynolds_avg(tmp, tke)

        ! Get tke rate of change
        ddt_tke = (tke - tke_old)/dt

        ! Get numerical dissipation (from filtering)
        dissipation_num = (tke_prefilter - tke)/dt

        ! Get Reynolds stresses
        call this%get_reynolds_stress(rho, u_pprime, v_pprime, w_pprime, Rij)

        ! Get mean velocity gradients
        call this%get_duidxj_avg(u_tilde, v_tilde, w_tilde, grad_u_tilde)

        ! Get production
        call this%get_production(rho_bar, Rij, grad_u_tilde, production)

        ! Get mean pressure gradients
        call this%get_gradp_avg(p_bar, grad_p_bar)

        ! Get turbulent mass fluxes
        call this%reynolds_avg(u_pprime, u_pprime_bar(:,:,:,1))
        call this%reynolds_avg(v_pprime, u_pprime_bar(:,:,:,2))
        call this%reynolds_avg(w_pprime, u_pprime_bar(:,:,:,3))

        ! Get fluctuating velocity gradients
        call gradient(this%gp, this%der, u_pprime, grad_u_pprime(:,:,:,1), grad_u_pprime(:,:,:,2), grad_u_pprime(:,:,:,3),&
                      -this%x_bc, this%y_bc, this%z_bc)
        call gradient(this%gp, this%der, v_pprime, grad_u_pprime(:,:,:,4), grad_u_pprime(:,:,:,5), grad_u_pprime(:,:,:,6),&
                       this%x_bc,-this%y_bc, this%z_bc)
        call gradient(this%gp, this%der, w_pprime, grad_u_pprime(:,:,:,7), grad_u_pprime(:,:,:,8), grad_u_pprime(:,:,:,9),&
                       this%x_bc, this%y_bc,-this%z_bc)


        ! Get pressure dilatation correlation terms (includes baropycnal work term)
        call this%get_p_dil(p, p_prime, u_pprime_bar, grad_p_bar, grad_u_pprime, p_dil_fluct, baropycnal, fluct_p_dil)

        ! Get mean and fluctuating shear stresses
        do i = 1, 6
            call this%reynolds_avg_and_fluct(tauij(:,:,:,i), tau_bar(:,:,:,i), tau_prime(:,:,:,i))
        end do

        ! Get dissipation terms
        call this%get_dissipation(tauij, tau_bar, tau_prime, u_pprime_bar, grad_u_pprime, dissipation, diss_mass_flux, diss_fluct)

        ! Write out data to output file
        call this%viz%start_viz(tsim)

        call this%viz%write_variable(rho_bar, 'rho_bar')
        call this%viz%write_variable(u_tilde, 'u_tilde')
        call this%viz%write_variable(v_tilde, 'v_tilde')
        call this%viz%write_variable(w_tilde, 'w_tilde')
        call this%viz%write_variable(p_bar,   'p_bar')

        call this%viz%write_variable(tke, 'TKE')
        call this%viz%write_variable(Rij(:,:,:,1), 'R_11')
        call this%viz%write_variable(Rij(:,:,:,2), 'R_12')
        call this%viz%write_variable(Rij(:,:,:,3), 'R_13')
        call this%viz%write_variable(Rij(:,:,:,4), 'R_22')
        call this%viz%write_variable(Rij(:,:,:,5), 'R_23')
        call this%viz%write_variable(Rij(:,:,:,6), 'R_33')

        call this%viz%write_variable(tau_bar(:,:,:,1), 'tau_bar_11')
        call this%viz%write_variable(tau_bar(:,:,:,2), 'tau_bar_12')
        call this%viz%write_variable(tau_bar(:,:,:,3), 'tau_bar_13')
        call this%viz%write_variable(tau_bar(:,:,:,4), 'tau_bar_22')
        call this%viz%write_variable(tau_bar(:,:,:,5), 'tau_bar_23')
        call this%viz%write_variable(tau_bar(:,:,:,6), 'tau_bar_33')

        call this%viz%write_variable(u_pprime_bar(:,:,:,1), 'u_pprime_bar')
        call this%viz%write_variable(u_pprime_bar(:,:,:,2), 'v_pprime_bar')
        call this%viz%write_variable(u_pprime_bar(:,:,:,3), 'w_pprime_bar')

        call this%viz%write_variable(ddt_tke,         'TKE_rate')
        call this%viz%write_variable(production,      'production')
        call this%viz%write_variable(p_dil_fluct,     'p_dil_fluct')
        call this%viz%write_variable(fluct_p_dil,     'fluct_p_dil')
        call this%viz%write_variable(baropycnal,      'baropycnal')
        call this%viz%write_variable(dissipation,     'dissipation')
        call this%viz%write_variable(diss_mass_flux,  'diss_mass_flux')
        call this%viz%write_variable(diss_fluct,      'diss_fluct')
        call this%viz%write_variable(dissipation_num, 'dissipation_num')

        call this%viz%end_viz()

    end subroutine

end module
