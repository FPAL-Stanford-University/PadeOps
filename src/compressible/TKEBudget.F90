module TKEBudgetMod

    use mpi
    use decomp_2d,        only: decomp_info, nrank
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero,half,one,two
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

        type(io_hdf5)       :: tke_viz
        type(io_hdf5)       :: mix_viz

        logical, dimension(3) :: averaging_directions = [.true., .true., .true.]
        integer, dimension(2) :: x_bc, y_bc, z_bc

        integer :: ns

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
        
        procedure          :: get_production
        procedure          :: get_dissipation
        procedure          :: get_transport
        procedure          :: get_p_dil
        procedure          :: get_mass_flux
        procedure          :: get_baropycnal
        procedure          :: tke_budget
        
        procedure          :: get_rhoPsi_bar
        procedure          :: mixing_budget
        final              :: destructor

    end type

    ! interface tkeBudget
    !     module procedure init
    ! end interface

contains

    subroutine init(this, gp, der, mesh, dx, dy, dz, ns, averaging_directions, outputdir, x_bc, y_bc, z_bc, reduce_precision)
        ! type(tkeBudget)                             :: this
        class(tkeBudget)                            :: this
        class(decomp_info), target,      intent(in) :: gp
        class(derivatives), target,      intent(in) :: der
        real(rkind), dimension(:,:,:,:), intent(in) :: mesh
        real(rkind),                     intent(in) :: dx, dy, dz
        integer,                         intent(in) :: ns
        logical, dimension(3),           intent(in) :: averaging_directions
        character(len=clen),             intent(in) :: outputdir
        integer, dimension(2),           intent(in) :: x_bc, y_bc, z_bc
        logical, optional,               intent(in) :: reduce_precision

        logical :: reduce_precision_
        ! integer, dimension(3) :: st
        ! integer :: ierr

        this%gp => gp
        this%der => der
        this%averaging_directions = averaging_directions

        this%x_bc = x_bc
        this%y_bc = y_bc
        this%z_bc = z_bc

        this%ns = ns

        reduce_precision_ = .true.
        if (present(reduce_precision)) reduce_precision_ = reduce_precision

        ! Initialize averaging object in the y  pencil
        this%avg = averaging(this%gp, 2, this%averaging_directions)

        ! if ( (this%avg%avg_dim /= 0) .and. (this%avg%avg_dim /= 2) ) then
        !     call GracefulExit("Only all direction or 1 direction averaging is supported for now.", 4589)
        ! end if

        ! Averaging in all directions
        if (this%avg%avg_dim == 0) then
            ! Initialize HDF5 output object (only x,y,z indices of 1 since we're outputting a 0D field)
            call this%tke_viz%init(mpi_comm_world, this%gp, 'y', outputdir, 'TKEBudget', reduce_precision=reduce_precision_, &
                               read_only=.false., subdomain_lo=[1,1,1], subdomain_hi=[1,1,1], &
                               jump_to_last=.true.)

            if (this%ns > 1) then
                ! Initialize HDF5 output object (only x,y,z indices of 1 since we're outputting a 0D field)
                call this%mix_viz%init(mpi_comm_world, this%gp, 'y', outputdir, 'MixBudget', reduce_precision=reduce_precision_, &
                                   read_only=.false., subdomain_lo=[1,1,1], subdomain_hi=[1,1,1], &
                                   jump_to_last=.true.)
            end if
        end if

        ! Averaging in X-Z directions
        if (this%avg%avg_dim == 1) then
            if (averaging_directions(2)) then
                call GracefulExit("Averaging in Y is not supported for 1D budgets! :(", 4589)
            end if

            call this%der_avg%init( this%gp, &
                                    dx, dy, dz, &
                                    averaging_directions(1), averaging_directions(2), averaging_directions(3), &
                                    der%getMethodx(), der%getMethody(), der%getMethodz() )

            call this%der_avg%set_ysz( [1, this%gp%ysz(2), 1] )                 ! Set ysz to make arrays 1D

            ! Initialize HDF5 output object (only x,z index of 1 since we're outputting a 1D field)
            call this%tke_viz%init(mpi_comm_world, this%gp, 'y', outputdir, 'TKEBudget', reduce_precision=reduce_precision_, &
                               read_only=.false., subdomain_lo=[1,1,1], subdomain_hi=[1,this%gp%ysz(2),1], &
                               jump_to_last=.true.)

            if (this%ns > 1) then
                ! Initialize HDF5 output object (only x,z index of 1 since we're outputting a 2D field)
                call this%mix_viz%init(mpi_comm_world, this%gp, 'y', outputdir, 'MixBudget', reduce_precision=reduce_precision_, &
                                   read_only=.false., subdomain_lo=[1,1,1], subdomain_hi=[1,this%gp%ysz(2),1], &
                                   jump_to_last=.true.)
            end if
        end if

        if (this%avg%avg_dim == 2) then
            call this%der_avg%init( this%gp, &
                                    dx, dy, dz, &
                                    averaging_directions(1), averaging_directions(2), averaging_directions(3), &
                                    der%getMethodx(), der%getMethody(), der%getMethodz() )

            if (averaging_directions(1)) then
                ! call this%der_avg%set_ysz( [1, this%gp%ysz(2), this%gp%ysz(3)] )                 ! Set ysz to make arrays 2D
                ! call this%der_avg%set_zsz( [1, this%gp%zsz(2), this%gp%zsz(3)] )                 ! Set zsz to make arrays 2D

                call GracefulExit("Averaging in X is not supported for 2D budgets! :(", 4589)
            end if

            if (averaging_directions(2)) then
                call GracefulExit("Averaging in Y is not supported for 2D budgets! :(", 4589)
            end if

            if (averaging_directions(3)) then
                call this%der_avg%set_xsz( [this%gp%xsz(1), this%gp%xsz(2), 1] )                 ! Set xsz to make arrays 2D
                call this%der_avg%set_ysz( [this%gp%ysz(1), this%gp%ysz(2), 1] )                 ! Set ysz to make arrays 2D

                call this%gp_avg%init(this%gp%xsz(1),this%gp%ysz(2),this%avg%xy_comm)

                ! Initialize HDF5 output object (only z index of 1 since we're outputting a 2D field)
                call this%tke_viz%init(mpi_comm_world, this%gp, 'y', outputdir, 'TKEBudget', reduce_precision=reduce_precision_, &
                                   read_only=.false., subdomain_lo=[1,1,1], subdomain_hi=[this%gp%xsz(1),this%gp%ysz(2),1], &
                                   jump_to_last=.true.)

                if (this%ns > 1) then
                    ! Initialize HDF5 output object (only z index of 1 since we're outputting a 2D field)
                    call this%mix_viz%init(mpi_comm_world, this%gp, 'y', outputdir, 'MixBudget', reduce_precision=reduce_precision_, &
                                       read_only=.false., subdomain_lo=[1,1,1], subdomain_hi=[this%gp%xsz(1),this%gp%ysz(2),1], &
                                       jump_to_last=.true.)
                end if
            end if

        end if

        ! Write the coordinates of subdomain out
        call this%tke_viz%write_coords(mesh)

        if (this%ns > 1) then
            ! Write the coordinates of subdomain out
            call this%mix_viz%write_coords(mesh)
        end if
    end subroutine

    subroutine destructor(this)
        type(tkeBudget), intent(inout) :: this

        nullify(this%gp)
        nullify(this%der)

        call this%tke_viz%destroy()
        call this%mix_viz%destroy()
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

        dudx => duidxj_avg(:,:,:,1); dudy => duidxj_avg(:,:,:,2); dudz => duidxj_avg(:,:,:,3);
        dvdx => duidxj_avg(:,:,:,4); dvdy => duidxj_avg(:,:,:,5); dvdz => duidxj_avg(:,:,:,6);
        dwdx => duidxj_avg(:,:,:,7); dwdy => duidxj_avg(:,:,:,8); dwdz => duidxj_avg(:,:,:,9);

        ! Set Z derivatives to zero
        dudz = zero; dvdz = zero; dwdz = zero

        if (this%averaging_directions(1)) then
            dudx = zero
            dvdx = zero
            dwdx = zero
        else
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
        end if

        if (this%averaging_directions(2)) then
            dudy = zero
            dvdy = zero
            dwdy = zero
        else
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
        if (this%averaging_directions(1)) then
            dpdx = zero
        else
            call this%gp_avg%transpose_y_to_x(p_avg(:,:,1),xbuf1(:,:,1))
            call this%der_avg%ddx(xbuf1,xbuf2, this%x_bc(1), this%x_bc(2))
            call this%gp_avg%transpose_x_to_y(xbuf2(:,:,1),dpdx(:,:,1))
        end if
        
        if (this%averaging_directions(2)) then
            dpdy = zero
        else
            call this%der_avg%ddy(p_avg,dpdy, this%y_bc(1), this%y_bc(2))
        end if
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
  
    subroutine get_dissipation(this, tauij_prime, grad_u_pprime, dissipation)
        class(tkeBudget),                                                                                 intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3),      6), target, intent(in)    :: tauij_prime
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3),      9), target, intent(in)    :: grad_u_pprime
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),           intent(out)   :: dissipation

        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3)) :: tmp
        real(rkind), dimension(:,:,:), pointer :: tauxx, tauxy, tauxz, tauyy, tauyz, tauzz
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz

        ! Store the tauij_prime fluctutations
        tauxx => tauij_prime(:,:,:,1); tauxy => tauij_prime(:,:,:,2); tauxz => tauij_prime(:,:,:,3);
                                        tauyy => tauij_prime(:,:,:,4); tauyz => tauij_prime(:,:,:,5);
                                                                        tauzz => tauij_prime(:,:,:,6);

        ! Store ddxi_upprime
        dudx => grad_u_pprime(:,:,:,1); dudy => grad_u_pprime(:,:,:,2); dudz => grad_u_pprime(:,:,:,3);
        dvdx => grad_u_pprime(:,:,:,4); dvdy => grad_u_pprime(:,:,:,5); dvdz => grad_u_pprime(:,:,:,6);
        dwdx => grad_u_pprime(:,:,:,7); dwdy => grad_u_pprime(:,:,:,8); dwdz => grad_u_pprime(:,:,:,9);

        ! Compute resolved dissipation
        tmp = tauxx*dudx + tauxy*dudy + tauxz*dudz &
            + tauxy*dvdx + tauyy*dvdy + tauyz*dvdz &
            + tauxz*dwdx + tauyz*dwdy + tauzz*dwdz
        call this%reynolds_avg(-tmp, dissipation)
    end subroutine

    subroutine get_transport(this, rho, pp, up, vp, wp, upp, vpp, wpp, &
            tauij_prime, tke, transport)
        class(tkeBudget), intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3)),& 
            intent(in)    :: rho, pp, up, vp, wp, upp, vpp, wpp, tke
        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3), 6),&
            target, intent(in)  :: tauij_prime
        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3)),& 
            intent(out) :: transport
        
        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3)) &
            :: tmp1, tmp2
        real(rkind), dimension(:,:,:), pointer :: tauxy, tauyy, tauzy !primes

        tauxy => tauij_prime(:,:,:,2); 
        tauyy => tauij_prime(:,:,:,4); 
        tauzy => tauij_prime(:,:,:,5); 
        
        ! x terms: set zero
        !tmp1 = rho*tke*upp + pp*up - tauxx
        !call this%gp_avg%transpose_y_to_x(tmp1,tmp2)
        !call this%der_avg%ddx(tmp2,tmp1, this%x_bc(1), this%x_bc(2))
        !call this%gp_avg%transpose_x_to_y(tmp1,tmp2)
        !tmpsum = tmp2

        ! y terms:
        tmp1 = rho*tke*vpp + pp*vp - (tauxy*upp + tauyy*vpp + tauzy*wpp)
        call this%der_avg%ddy(tmp1,tmp2, this%y_bc(1), this%y_bc(2))
        
        ! z terms: set to zero? 
        !tmp1 = rho*tke*wpp + pp*wp - tauzz
        !call this%gp_avg%transpose_y_to_z(tmp1,tmp2)
        !call this%der_avg%ddz(tmp2,tmp1, this%z_bc(1), this%z_bc(2))
        !call this%gp_avg%transpose_z_to_y(tmp1,tmp2)
        !tmpsum = tmpsum + tmp2

        call this%reynolds_avg(-tmp2, transport)
    end subroutine

    subroutine get_p_dil(this, p_prime, grad_u_pprime, p_dil )
        class(tkeBudget),                                                                                 intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),                 intent(in)    :: p_prime
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3),      9), target, intent(in)    :: grad_u_pprime
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),           intent(out)   :: p_dil

        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3)) :: tmp
        real(rkind), dimension(:,:,:), pointer :: dudx, dvdy, dwdz

        dudx => grad_u_pprime(:,:,:,1);
        dvdy => grad_u_pprime(:,:,:,5);
        dwdz => grad_u_pprime(:,:,:,9);

        ! Pressure dilatation (pressure strain contracted)
        tmp = p_prime*(dudx + dvdy + dwdz)
        call this%reynolds_avg(tmp, p_dil)
    end subroutine

    subroutine get_mass_flux(this, tau_bar, u_pprime_bar, diss_mass_flux )
        class(tkeBudget),                                                                                 intent(inout) :: this
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),6), target, intent(in)    :: tau_bar
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),3), target, intent(in)    :: u_pprime_bar
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),           intent(out)   :: diss_mass_flux 

        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3)) :: tmp1,tmp2,tmpsum
        real(rkind), dimension(:,:,:), pointer :: tauxx, tauxy, tauxz, tauyy, tauyz, tauzz, upp,vpp,wpp

        upp => u_pprime_bar(:,:,:,1);
        vpp => u_pprime_bar(:,:,:,2);
        wpp => u_pprime_bar(:,:,:,3);
        tauxx => tau_bar(:,:,:,1); tauxy => tau_bar(:,:,:,2); tauxz => tau_bar(:,:,:,3);
                                   tauyy => tau_bar(:,:,:,4); tauyz => tau_bar(:,:,:,5);
                                                              tauzz => tau_bar(:,:,:,6);
        ! x-direction     
        tmpsum = 0.D0
        !call this%gp_avg%transpose_y_to_x(tauxx,tmp2)
        !call this%der_avg%ddx(tmp2,tmp1, this%x_bc(1), this%x_bc(2))
        !call this%gp_avg%transpose_x_to_y(tmp1,tmp2)
        !tmpsum = tmpsum + tmp2 ! ddx tauxx 
        call this%der_avg%ddy(tauxy,tmp2, this%y_bc(1), this%y_bc(2))
        tmpsum = tmpsum + tmp2 ! ddy tauxy
        !call this%gp_avg%transpose_y_to_z(tauxz,tmp2)
        !call this%der_avg%ddz(tmp2,tmp1, this%z_bc(1), this%z_bc(2))
        !call this%gp_avg%transpose_z_to_y(tmp1,tmp2)
        tmpsum = tmpsum + tmp2 ! ddz tauxz
        diss_mass_flux = upp*tmpsum
        
        ! y-direction     
        !call this%gp_avg%transpose_y_to_x(tauxy,tmp2)
        !call this%der_avg%ddx(tmp2,tmp1, this%x_bc(1), this%x_bc(2))
        !call this%gp_avg%transpose_x_to_y(tmp1,tmp2)
        !tmpsum = tmpsum + tmp2 ! ddx tauyx 
        call this%der_avg%ddy(tauyy,tmp2, this%y_bc(1), this%y_bc(2))
        tmpsum = tmpsum + tmp2 ! ddy tauyy
        diss_mass_flux = vpp*tmpsum
        
        ! z-direction     
        !call this%gp_avg%transpose_y_to_x(tauxy,tmp2)
        !call this%der_avg%ddx(tmp2,tmp1, this%x_bc(1), this%x_bc(2))
        !call this%gp_avg%transpose_x_to_y(tmp1,tmp2)
        !tmpsum = tmpsum + tmp2 ! ddx tauyx 
        call this%der_avg%ddy(tauyz,tmp2, this%y_bc(1), this%y_bc(2))
        tmpsum = tmpsum + tmp2 ! ddy tauyy
        diss_mass_flux = wpp*tmpsum
    end subroutine

    subroutine get_baropycnal(this, grad_p_bar, u_pprime_bar, baropycnal )
        class(tkeBudget),                                                                         intent(inout) :: this
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),3), intent(in)    :: grad_p_bar
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),3), intent(in)    :: u_pprime_bar
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),   intent(out)   :: baropycnal 
        
        baropycnal = & 
              u_pprime_bar(:,:,:,1)*grad_p_bar(:,:,:,1) &
            + u_pprime_bar(:,:,:,2)*grad_p_bar(:,:,:,2) &
            + u_pprime_bar(:,:,:,3)*grad_p_bar(:,:,:,3)
        baropycnal = -baropycnal
    end subroutine

    subroutine tke_budget(this, rho, u, v, w, p, tauij, tke_old, tke_prefilter, tke_postfilter, tsim, dt)
        use RKCoeffs, only: RK45_steps
        class(tkeBudget),                                                                                  intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),                  intent(in)    :: rho, u, v, w, p
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3),      6),          intent(in)    :: tauij
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)),            intent(in)    :: tke_old
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),RK45_steps), intent(in)    :: tke_prefilter
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),RK45_steps), intent(in)    :: tke_postfilter
        real(rkind),                                                                                       intent(in)    :: tsim, dt

        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: ddt_tke
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: production
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: p_dil
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: baropycnal
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: dissipation
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: diss_mass_flux
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: diss_fluct
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: dissipation_num
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3)) :: transport 

        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3))   :: rho_bar, p_bar
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3))   :: u_tilde, v_tilde, w_tilde
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3))   :: tke
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),6) :: Rij
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),6) :: tau_bar
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),9) :: grad_u_tilde
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),3) :: grad_p_bar
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),3) :: u_pprime_bar

        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3))   :: u_pprime, v_pprime, w_pprime, p_prime, u_prime, v_prime, w_prime
        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3),9) :: grad_u_pprime
        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3),6) :: tau_prime
        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3))   :: tmp

        integer :: i
        character(len=clen) :: varname

        call this%reynolds_avg(rho, rho_bar)
        call this%favre_avg_and_fluct(rho, u, u_tilde, u_pprime)
        call this%favre_avg_and_fluct(rho, v, v_tilde, v_pprime)
        call this%favre_avg_and_fluct(rho, w, w_tilde, w_pprime)
        
        call this%reynolds_avg_and_fluct(u, tmp, u_prime)
        call this%reynolds_avg_and_fluct(v, tmp, v_prime)
        call this%reynolds_avg_and_fluct(w, tmp, w_prime)   
        call this%reynolds_avg_and_fluct(p, p_bar, p_prime)

        ! Get tke and rate of change
        tmp = half*rho*( u_pprime*u_pprime + v_pprime*v_pprime + w_pprime*w_pprime )
        call this%reynolds_avg(tmp, tke)
        ddt_tke = (tke - tke_old)/dt

        ! Get Reynolds stresses
        call this%get_reynolds_stress(rho, u_pprime, v_pprime, w_pprime, Rij)

        ! Get mean gradients
        call this%get_duidxj_avg(u_tilde, v_tilde, w_tilde, grad_u_tilde)
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

        ! Get mean and fluctuating shear stresses
        do i = 1, 6
            call this%reynolds_avg_and_fluct(tauij(:,:,:,i), tau_bar(:,:,:,i), tau_prime(:,:,:,i))
        end do
        
        !!!!!!!!!!!!!!!!!!!!!! TKE BUDGET TERMS !!!!!!!!!!!!!!!!!!!!!!!

        call this%get_production(rho_bar, Rij, grad_u_tilde, production)

        call this%get_dissipation(tau_prime, grad_u_pprime, dissipation)
        dissipation_num = sum(tke_prefilter - tke_postfilter, 4)/dt

        call this%get_transport(rho, p_prime, u_prime, v_prime, w_prime, &
                u_pprime, v_pprime, w_pprime, tau_prime, tmp, transport)

        call this%get_p_dil( p_prime, grad_u_pprime, p_dil )
       
        call this%get_mass_flux(tau_bar, u_pprime_bar, diss_mass_flux )
        
        call this%get_baropycnal(grad_p_bar, u_pprime_bar, baropycnal )
            
        !!!!!!!!!!!!!!!!!!!!!! WRITE !!!!!!!!!!!!!!!!!!!!!!!
        
        call this%tke_viz%start_viz(tsim)
        call this%tke_viz%write_attribute(1, [dt], 'dt', '/')

        call this%tke_viz%write_variable(rho_bar, 'rho_bar')
        call this%tke_viz%write_variable(u_tilde, 'u_tilde')
        call this%tke_viz%write_variable(v_tilde, 'v_tilde')
        call this%tke_viz%write_variable(w_tilde, 'w_tilde')
        call this%tke_viz%write_variable(p_bar,   'p_bar')

        call this%tke_viz%write_variable(tke, 'tke')
        call this%tke_viz%write_variable(tke_old, 'tke_old')
        do i=1,RK45_steps
            write(varname, '(A,I1.1)') 'tke_prefilter_', i
            call this%tke_viz%write_variable(tke_prefilter(:,:,:,i), trim(varname))
            write(varname, '(A,I1.1)') 'tke_postfilter_', i
            call this%tke_viz%write_variable(tke_postfilter(:,:,:,i), trim(varname))
        end do

        call this%tke_viz%write_variable(Rij(:,:,:,1), 'R_11')
        call this%tke_viz%write_variable(Rij(:,:,:,2), 'R_12')
        call this%tke_viz%write_variable(Rij(:,:,:,3), 'R_13')
        call this%tke_viz%write_variable(Rij(:,:,:,4), 'R_22')
        call this%tke_viz%write_variable(Rij(:,:,:,5), 'R_23')
        call this%tke_viz%write_variable(Rij(:,:,:,6), 'R_33')

        call this%tke_viz%write_variable(tau_bar(:,:,:,1), 'tau_bar_11')
        call this%tke_viz%write_variable(tau_bar(:,:,:,2), 'tau_bar_12')
        call this%tke_viz%write_variable(tau_bar(:,:,:,3), 'tau_bar_13')
        call this%tke_viz%write_variable(tau_bar(:,:,:,4), 'tau_bar_22')
        call this%tke_viz%write_variable(tau_bar(:,:,:,5), 'tau_bar_23')
        call this%tke_viz%write_variable(tau_bar(:,:,:,6), 'tau_bar_33')

        call this%tke_viz%write_variable(u_pprime_bar(:,:,:,1), 'u_pprime_bar')
        call this%tke_viz%write_variable(u_pprime_bar(:,:,:,2), 'v_pprime_bar')
        call this%tke_viz%write_variable(u_pprime_bar(:,:,:,3), 'w_pprime_bar')

        call this%tke_viz%write_variable(ddt_tke,         'TKE_rate')
        call this%tke_viz%write_variable(production,      'production')
        call this%tke_viz%write_variable(dissipation,     'dissipation')
        call this%tke_viz%write_variable(dissipation_num, 'dissipation_num')
        call this%tke_viz%write_variable(transport,       'transport')
        call this%tke_viz%write_variable(p_dil,           'p_dil')
        call this%tke_viz%write_variable(diss_mass_flux,  'diss_mass_flux')
        call this%tke_viz%write_variable(baropycnal,      'baropycnal')

        call this%tke_viz%end_viz()
    end subroutine


    subroutine get_rhoPsi_bar(this, rho, Ys, rhoPsi_bar)
        class(tkeBudget),                                                                                  intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),                  intent(in)    :: rho
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3),      this%ns),    intent(in)    :: Ys
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),this%ns),    intent(out)   :: rhoPsi_bar

        integer :: m

        do m = 1,this%ns
            ! Get rho*Psi average
            call this%reynolds_avg(rho*Ys(:,:,:,m)*(one-Ys(:,:,:,m)), rhoPsi_bar(:,:,:,m))
        end do
    end subroutine


    subroutine mixing_budget(this, rho, u, v, w, Ys, diff, rhoPsi_old, rhoPsi_prefilter, rhoPsi_postfilter, tsim, dt)
        use RKCoeffs, only: RK45_steps
        class(tkeBudget),                                                                                  intent(inout) :: this
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3)),                  intent(in)    :: rho, u, v, w
        real(rkind), dimension(this%avg%sz(1),      this%avg%sz(2),      this%avg%sz(3),      this%ns),    intent(in)    :: Ys, diff
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),this%ns),    intent(in)    :: rhoPsi_old
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),this%ns,RK45_steps), intent(in)    :: rhoPsi_prefilter
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),this%ns,RK45_steps), intent(in)    :: rhoPsi_postfilter
        real(rkind),                                                                                       intent(in)    :: tsim, dt

        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3))         :: rho_bar
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3))         :: u_tilde, v_tilde, w_tilde
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),this%ns) :: Psi_tilde

        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),this%ns) :: ddt_Psi
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),this%ns) :: scalar_dissipation
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),this%ns) :: num_scalar_dissipation
        real(rkind), dimension(this%avg%avg_size(1),this%avg%avg_size(2),this%avg%avg_size(3),this%ns) :: mix_turb_flux_x, mix_turb_flux_y, &
                                                                                                          mix_turb_flux_z

        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3))          :: u_pprime, v_pprime, w_pprime
        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3))          :: Psi_pprime
        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3), 3)       :: grad_Ys
        real(rkind), dimension(this%avg%sz(1), this%avg%sz(2), this%avg%sz(3))          :: tmp

        integer :: m, i
        character(len=clen) :: varname

        ! Get density average
        call this%reynolds_avg(rho, rho_bar)

        ! Get velocity favre average and fluctuations
        call this%favre_avg_and_fluct(rho, u, u_tilde, u_pprime)
        call this%favre_avg_and_fluct(rho, v, v_tilde, v_pprime)
        call this%favre_avg_and_fluct(rho, w, w_tilde, w_pprime)

        do m = 1,this%ns
            ! Get Psi favre average and fluctuations
            call this%favre_avg_and_fluct(rho, Ys(:,:,:,m)*(one-Ys(:,:,:,m)), Psi_tilde(:,:,:,m), Psi_pprime)

            ! Get Psi rate of change
            ddt_Psi(:,:,:,m) = (rho_bar*Psi_tilde(:,:,:,m) - rhoPsi_old(:,:,:,m))/dt

            ! Get massfraction gradients
            call gradient(this%gp, this%der, Ys(:,:,:,m), grad_Ys(:,:,:,1), grad_Ys(:,:,:,2), grad_Ys(:,:,:,3),&
                          this%x_bc, this%y_bc, this%z_bc)

            ! Get scalar_dissipation
            tmp = two * rho * diff(:,:,:,m) * (grad_Ys(:,:,:,1)**2 + grad_Ys(:,:,:,2)**2 + grad_Ys(:,:,:,3)**2)
            call this%reynolds_avg(tmp, scalar_dissipation(:,:,:,m))

            ! Get turbulent flux
            tmp = rho * Psi_pprime * u_pprime
            call this%reynolds_avg(tmp, mix_turb_flux_x(:,:,:,m))
            tmp = rho * Psi_pprime * v_pprime
            call this%reynolds_avg(tmp, mix_turb_flux_y(:,:,:,m))
            tmp = rho * Psi_pprime * w_pprime
            call this%reynolds_avg(tmp, mix_turb_flux_z(:,:,:,m))
        end do

        ! Get numerical dissipation (from filtering)
        num_scalar_dissipation = sum(rhoPsi_prefilter - rhoPsi_postfilter, 5)/dt

        ! Write out data to output file
        call this%mix_viz%start_viz(tsim)
        call this%mix_viz%write_attribute(1, [dt], 'dt', '/')

        call this%mix_viz%write_variable(rho_bar, 'rho_bar')

        do m = 1,this%ns
            write(varname, '(A,I2.2)') 'Psi_tilde_', m
            call this%mix_viz%write_variable(Psi_tilde(:,:,:,m), trim(varname))

            do i=1,RK45_steps
                write(varname, '(A,I2.2,A,I1.1)') 'rhoPsi_prefilter_', m, '_', i
                call this%tke_viz%write_variable(rhoPsi_prefilter(:,:,:,m,i), trim(varname))
                write(varname, '(A,I2.2,A,I1.1)') 'rhoPsi_postfilter_', m, '_', i
                call this%tke_viz%write_variable(rhoPsi_postfilter(:,:,:,m,i), trim(varname))
            end do

            write(varname, '(A,I2.2)') 'ddt_Psi_tilde_', m
            call this%mix_viz%write_variable(ddt_Psi(:,:,:,m), trim(varname))

            write(varname, '(A,I2.2)') 'scalar_dissipation_', m
            call this%mix_viz%write_variable(scalar_dissipation(:,:,:,m), trim(varname))

            write(varname, '(A,I2.2)') 'num_scalar_dissipation_', m
            call this%mix_viz%write_variable(num_scalar_dissipation(:,:,:,m), trim(varname))

            write(varname, '(A,I2.2)') 'mix_turb_flux_x_', m
            call this%mix_viz%write_variable(mix_turb_flux_x(:,:,:,m), trim(varname))

            write(varname, '(A,I2.2)') 'mix_turb_flux_y_', m
            call this%mix_viz%write_variable(mix_turb_flux_y(:,:,:,m), trim(varname))

            write(varname, '(A,I2.2)') 'mix_turb_flux_z_', m
            call this%mix_viz%write_variable(mix_turb_flux_z(:,:,:,m), trim(varname))

        end do

        call this%mix_viz%end_viz()
    end subroutine


end module
