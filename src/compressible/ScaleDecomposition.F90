module ScaleDecompositionMod

    use mpi
    use decomp_2d,        only: decomp_info, nrank
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero,half,one,two
    use exits,            only: GracefulExit, message
    use io_hdf5_stuff,    only: io_hdf5
    use DerivativesMod,   only: derivatives
    use FiltersMod,       only: filters

    implicit none

    type :: scaleDecomposition
        ! private

        type(decomp_info), pointer :: gp
        type(filters),     pointer :: gfil
        type(derivatives), pointer :: der
        
        type(io_hdf5)       :: tke_viz
        type(io_hdf5)       :: mix_viz

        integer :: ns
        integer :: species_id

        integer :: nfilters = 4

        real(rkind) :: dx, dy, dz
        integer, dimension(2) :: x_bc, y_bc, z_bc

    contains

        procedure          :: init
        procedure          :: filter
        procedure          :: favre_filter
        procedure          :: get_duidxj
        procedure          :: get_kinetic_energies
        procedure          :: get_kinetic_energies_turb_stress
        procedure          :: get_production
        procedure          :: get_mass_flux
        procedure          :: get_baropycnal
        procedure          :: get_pressure_dilatation
        procedure          :: get_dissipation
        procedure          :: tke_budget
        procedure          :: get_rhoPsi
        procedure          :: get_rhoPsi_turb_flux_generation
        procedure          :: get_diffusive_flux_generation
        procedure          :: get_transport_divergence
        procedure          :: mix_budget
        final              :: destructor

    end type

contains

    subroutine init(this, gp, der, gfil, mesh, dx, dy, dz, ns, x_bc, y_bc, z_bc, inputfile)
        class(scaleDecomposition)                   :: this
        class(decomp_info), target,      intent(in) :: gp
        class(derivatives), target,      intent(in) :: der
        class(filters),     target,      intent(in) :: gfil
        real(rkind), dimension(:,:,:,:), intent(in) :: mesh
        real(rkind),                     intent(in) :: dx, dy, dz
        integer,                         intent(in) :: ns
        integer, dimension(2),           intent(in) :: x_bc, y_bc, z_bc
        character(len=clen),             intent(in) :: inputfile

        logical             :: reduce_precision = .true.
        integer             :: nfilters = 4
        integer, dimension(3) :: subdomain_lo, subdomain_hi
        character(len=clen) :: outputdir, outputprefix
        integer             :: species_id = 1

        integer :: ioUnit

        namelist /SCALEDECOMP/  reduce_precision, nfilters, outputdir, subdomain_lo, subdomain_hi, species_id


        this%gp => gp
        this%der => der
        this%gfil => gfil

        this%dx = dx
        this%dy = dy
        this%dz = dz

        this%ns = ns
        this%species_id = species_id

        ! Set to full domain by default
        subdomain_lo = [1, 1, 1]
        subdomain_hi = [this%gp%xsz(1), this%gp%ysz(2), this%gp%zsz(3)]

        ! Read in parameters from the input file
        ioUnit = 121
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=SCALEDECOMP)
        close(ioUnit)

        this%x_bc = x_bc
        this%y_bc = y_bc
        this%z_bc = z_bc

        this%nfilters = nfilters

        write(outputprefix, '(A,I4.4)') 'TKE_scaledecomp_F', this%nfilters
        call message("Will save data with prefix "//trim(outputprefix))

        ! Initialize HDF5 output object (only z index of 1 since we're outputting a 1D field)
        call this%tke_viz%init(mpi_comm_world, this%gp, 'y', trim(outputdir), trim(outputprefix), reduce_precision=reduce_precision, &
                               read_only=.false., subdomain_lo=subdomain_lo, subdomain_hi=subdomain_hi, &
                               jump_to_last=.true.)

        write(outputprefix, '(A,I4.4)') 'Mix_scaledecomp_F', this%nfilters
        call message("Will save data with prefix "//trim(outputprefix))

        call this%mix_viz%init(mpi_comm_world, this%gp, 'y', trim(outputdir), trim(outputprefix), reduce_precision=reduce_precision, &
                               read_only=.false., subdomain_lo=subdomain_lo, subdomain_hi=subdomain_hi, &
                               jump_to_last=.true.)

        ! Write the coordinates of subdomain out
        call this%tke_viz%write_coords(mesh)
        call this%mix_viz%write_coords(mesh)

    end subroutine

    subroutine destructor(this)
        type(scaleDecomposition), intent(inout) :: this

        nullify(this%gp)
        nullify(this%der)
        nullify(this%gfil)

        call this%tke_viz%destroy()
        call this%mix_viz%destroy()

    end subroutine

    subroutine filter(this, var, var_f, x_bc, y_bc, z_bc)
        use operators, only: filter3D
        class(scaleDecomposition),                                                        intent(in)  :: this
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)),           intent(in)  :: var
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)),           intent(out) :: var_f
        integer,     dimension(2),                                              optional, intent(in)  :: x_bc, y_bc, z_bc

        integer,     dimension(2) :: x_bc_, y_bc_, z_bc_

        x_bc_ = this%x_bc; if (present(x_bc)) x_bc_ = x_bc
        y_bc_ = this%y_bc; if (present(y_bc)) y_bc_ = y_bc
        z_bc_ = this%z_bc; if (present(z_bc)) z_bc_ = z_bc

        var_f = var
        call filter3D( this%gp, this%gfil, var_f, this%nfilters, x_bc_, y_bc_, z_bc_ )
    end subroutine

    subroutine favre_filter(this, rho, rho_f, var, var_ff, x_bc, y_bc, z_bc)
        use operators, only: filter3D
        class(scaleDecomposition),                                                        intent(in)  :: this
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)),           intent(in)  :: rho, rho_f, var
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)),           intent(out) :: var_ff
        integer,     dimension(2),                                              optional, intent(in)  :: x_bc, y_bc, z_bc

        integer,     dimension(2) :: x_bc_, y_bc_, z_bc_

        x_bc_ = this%x_bc; if (present(x_bc)) x_bc_ = x_bc
        y_bc_ = this%y_bc; if (present(y_bc)) y_bc_ = y_bc
        z_bc_ = this%z_bc; if (present(z_bc)) z_bc_ = z_bc

        var_ff = rho*var
        call filter3D( this%gp, this%gfil, var_ff, this%nfilters, x_bc_, y_bc_, z_bc_ )
        var_ff = var_ff / rho_f
    end subroutine

    subroutine get_duidxj(this, u, v, w, duidxj)
        use operators, only: gradient
        class(scaleDecomposition),                                                         intent(in)  :: this
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)),            intent(in)  :: u, v, w
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3), 9), target, intent(out) :: duidxj

        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);

        call gradient(this%gp, this%der, u, dudx, dudy, dudz,-this%x_bc, this%y_bc, this%z_bc)
        call gradient(this%gp, this%der, v, dvdx, dvdy, dvdz, this%x_bc,-this%y_bc, this%z_bc)
        call gradient(this%gp, this%der, w, dwdx, dwdy, dwdz, this%x_bc, this%y_bc,-this%z_bc)

    end subroutine

    subroutine get_kinetic_energies(this, rho, u, v, w, KE_L, KE_S)
        class(scaleDecomposition),                                                 intent(in)  :: this
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(in)  :: rho, u, v, w
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(out) :: KE_L, KE_S

        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ) :: rho_f, u_ff, v_ff, w_ff, tmp

        ! Get filtered density, velocities, pressure
        call this%filter(rho, rho_f, this%x_bc, this%y_bc, this%z_bc)
        call this%favre_filter(rho, rho_f, u, u_ff,-this%x_bc, this%y_bc, this%z_bc)
        call this%favre_filter(rho, rho_f, v, v_ff, this%x_bc,-this%y_bc, this%z_bc)
        call this%favre_filter(rho, rho_f, w, w_ff, this%x_bc, this%y_bc,-this%z_bc)

        ! Large scale kinetic energy
        KE_L = half * rho_f * (u_ff*u_ff + v_ff*v_ff + w_ff*w_ff)

        ! Turbulent stresses: tilde(u_i*u_j) - tilde(u_i)*tilde(u_j)
        call this%favre_filter(rho, rho_f, u*u, tmp, this%x_bc, this%y_bc, this%z_bc)
        KE_S = (tmp - u_ff*u_ff)           ! ts_11 = tilde(u*u) - tilde(u)*tilde(u)

        call this%favre_filter(rho, rho_f, v*v, tmp, this%x_bc, this%y_bc, this%z_bc)
        KE_S = KE_S + (tmp - v_ff*v_ff)    ! ts_22 = tilde(v*v) - tilde(v)*tilde(v)

        call this%favre_filter(rho, rho_f, w*w, tmp, this%x_bc, this%y_bc, this%z_bc)
        KE_S = KE_S + (tmp - w_ff*w_ff)    ! ts_33 = tilde(w*w) - tilde(w)*tilde(w)

        ! Small scale kinetic energy
        KE_S = half * rho_f * KE_S

    end subroutine

    subroutine get_kinetic_energies_turb_stress(this, rho, rho_f, u, u_ff, v, v_ff, w, w_ff, KE_L, KE_S, turb_stress)
        class(scaleDecomposition),                                                 intent(in)  :: this
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(in)  :: rho, rho_f, u, u_ff, v, v_ff, w, w_ff
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(out) :: KE_L, KE_S
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3), 6), intent(out) :: turb_stress

        ! Large scale kinetic energy
        KE_L = half * rho_f * (u_ff*u_ff + v_ff*v_ff + w_ff*w_ff)

        ! Turbulent stresses: tilde(u_i*u_j) - tilde(u_i)*tilde(u_j)
        call this%favre_filter(rho, rho_f, u*u, turb_stress(:,:,:,1), this%x_bc, this%y_bc, this%z_bc)
        turb_stress(:,:,:,1) = turb_stress(:,:,:,1) - u_ff*u_ff    ! ts_11 = tilde(u*u) - tilde(u)*tilde(u)

        call this%favre_filter(rho, rho_f, u*v, turb_stress(:,:,:,2),-this%x_bc,-this%y_bc, this%z_bc)
        turb_stress(:,:,:,2) = turb_stress(:,:,:,2) - u_ff*v_ff    ! ts_12 = tilde(u*v) - tilde(u)*tilde(v)

        call this%favre_filter(rho, rho_f, u*w, turb_stress(:,:,:,3),-this%x_bc, this%y_bc,-this%z_bc)
        turb_stress(:,:,:,3) = turb_stress(:,:,:,3) - u_ff*w_ff    ! ts_13 = tilde(u*w) - tilde(u)*tilde(w)

        call this%favre_filter(rho, rho_f, v*v, turb_stress(:,:,:,4), this%x_bc, this%y_bc, this%z_bc)
        turb_stress(:,:,:,4) = turb_stress(:,:,:,4) - v_ff*v_ff    ! ts_22 = tilde(v*v) - tilde(v)*tilde(v)

        call this%favre_filter(rho, rho_f, v*w, turb_stress(:,:,:,5), this%x_bc,-this%y_bc,-this%z_bc)
        turb_stress(:,:,:,5) = turb_stress(:,:,:,5) - v_ff*w_ff    ! ts_23 = tilde(v*w) - tilde(v)*tilde(w)

        call this%favre_filter(rho, rho_f, w*w, turb_stress(:,:,:,6), this%x_bc, this%y_bc, this%z_bc)
        turb_stress(:,:,:,6) = turb_stress(:,:,:,6) - w_ff*w_ff    ! ts_33 = tilde(w*w) - tilde(w)*tilde(w)

        ! Small scale kinetic energy
        KE_S = half * rho_f * (turb_stress(:,:,:,1) + turb_stress(:,:,:,4) + turb_stress(:,:,:,6))

    end subroutine

    subroutine get_production(this, rho_f, u_ff, v_ff, w_ff, turb_stress, duidxj_ff, production)
        class(scaleDecomposition),                                                         intent(in)  :: this
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ),         intent(in)  :: rho_f, u_ff, v_ff, w_ff
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3), 6), target, intent(in)  :: turb_stress
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3), 9), target, intent(out) :: duidxj_ff
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ),         intent(out) :: production

        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: TSxx, TSxy, TSxz, TSyy, TSyz, TSzz

        dudx => duidxj_ff(:,:,:,1); dudy => duidxj_ff(:,:,:,2); dudz => duidxj_ff(:,:,:,3);
        dvdx => duidxj_ff(:,:,:,4); dvdy => duidxj_ff(:,:,:,5); dvdz => duidxj_ff(:,:,:,6);
        dwdx => duidxj_ff(:,:,:,7); dwdy => duidxj_ff(:,:,:,8); dwdz => duidxj_ff(:,:,:,9);

        TSxx => turb_stress(:,:,:,1); TSxy => turb_stress(:,:,:,2); TSxz => turb_stress(:,:,:,3);
                                      TSyy => turb_stress(:,:,:,4); TSyz => turb_stress(:,:,:,5);
                                                                    TSzz => turb_stress(:,:,:,6);

        call this%get_duidxj(u_ff,v_ff,w_ff,duidxj_ff)

        production = dudx*TSxx + dudy*TSxy + dudz*TSxz &
                   + dvdx*TSxy + dvdy*TSyy + dvdz*TSyz &
                   + dwdx*TSxz + dwdy*TSyz + dwdz*TSzz

        production = -rho_f * production

    end subroutine

    subroutine get_mass_flux(this, rho, rho_f, u, v, w, mass_flux)
        class(scaleDecomposition),                                                 intent(in)  :: this
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(in)  :: rho, rho_f, u, v, w
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3), 3), intent(out) :: mass_flux

        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)) :: vel_f

        call this%filter(rho*u, mass_flux(:,:,:,1), -this%x_bc, this%y_bc, this%z_bc)
        call this%filter(u, vel_f, -this%x_bc, this%y_bc, this%z_bc)
        mass_flux(:,:,:,1) = mass_flux(:,:,:,1) - rho_f*vel_f ! filter(rho*u) - filter(rho)*filter(u)
        
        call this%filter(rho*v, mass_flux(:,:,:,2), this%x_bc, -this%y_bc, this%z_bc)
        call this%filter(v, vel_f, this%x_bc, -this%y_bc, this%z_bc)
        mass_flux(:,:,:,2) = mass_flux(:,:,:,2) - rho_f*vel_f ! filter(rho*v) - filter(rho)*filter(v)
        
        call this%filter(rho*w, mass_flux(:,:,:,3), this%x_bc, this%y_bc, -this%z_bc)
        call this%filter(w, vel_f, this%x_bc, this%y_bc, -this%z_bc)
        mass_flux(:,:,:,3) = mass_flux(:,:,:,3) - rho_f*vel_f ! filter(rho*w) - filter(rho)*filter(w)
        
    end subroutine

    subroutine get_baropycnal(this, rho_f, p_f, mass_flux, baropycnal)
        use operators, only: gradient
        class(scaleDecomposition),                                                 intent(in)  :: this
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(in)  :: rho_f, p_f
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3), 3), intent(in)  :: mass_flux
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(out) :: baropycnal

        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)) :: dpdx, dpdy, dpdz

        call gradient(this%gp, this%der, p_f, dpdx, dpdy, dpdz, this%x_bc, this%y_bc, this%z_bc)

        baropycnal = ( dpdx*mass_flux(:,:,:,1) + dpdy*mass_flux(:,:,:,2) + dpdz*mass_flux(:,:,:,3) ) / rho_f
    end subroutine

    subroutine get_pressure_dilatation(this, p, p_f, u, v, w, pdil_L, pdil_S)
        use operators, only: divergence
        class(scaleDecomposition),                                              intent(in)  :: this
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)), intent(in)  :: p, p_f, u, v, w
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)), intent(out) :: pdil_L, pdil_S

        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)) :: u_f, v_f, w_f, divu

        call this%filter(u, u_f,-this%x_bc, this%y_bc, this%z_bc)
        call this%filter(v, v_f, this%x_bc,-this%y_bc, this%z_bc)
        call this%filter(w, w_f, this%x_bc, this%y_bc,-this%z_bc)

        ! call divergence( this%gp, this%der, u_f, v_f, w_f, divu, -this%x_bc, -this%y_bc, -this%z_bc )
        call divergence( this%gp, this%der, u_f, v_f, w_f, divu, this%x_bc, this%y_bc, this%z_bc )
        pdil_L = p_f*divu ! filter(p) * divergence( filter(velocity) )

        ! call divergence( this%gp, this%der, u, v, w, divu, -this%x_bc, -this%y_bc, -this%z_bc )
        call divergence( this%gp, this%der, u, v, w, divu, this%x_bc, this%y_bc, this%z_bc )
        divu = p*divu
        call this%filter(divu, pdil_S, this%x_bc, this%y_bc, this%z_bc)
        pdil_S = pdil_S - pdil_L    ! filter( p * divergence(velocity) ) - filter(p)*divergence( filter(velocity) )

    end subroutine

    subroutine get_dissipation(this, u, v, w, duidxj_ff, tauij, diss_L, diss_S)
        class(scaleDecomposition),                                                         intent(in)  :: this
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ),         intent(in)  :: u, v, w
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3), 9), target, intent(in)  :: duidxj_ff
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3), 6), target, intent(in)  :: tauij
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ),         intent(out) :: diss_L, diss_S

        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3), 9), target :: duidxj
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   )         :: tmp

        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: tauxx, tauxy, tauxz, tauyy, tauyz, tauzz

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        tauxx => tauij(:,:,:,1); tauxy => tauij(:,:,:,2); tauxz => tauij(:,:,:,3);
                                 tauyy => tauij(:,:,:,4); tauyz => tauij(:,:,:,5);
                                                          tauzz => tauij(:,:,:,6);

        call this%get_duidxj(u,v,w,duidxj)
        ! call get_tauij(duidxj,mir%mu,mir%bulk,tauij)

        ! Get full dissipation in diss_S
        diss_L = dudx*tauxx + dudy*tauxy + dudz*tauxz &
               + dvdx*tauxy + dvdy*tauyy + dvdz*tauyz &
               + dwdx*tauxz + dwdy*tauyz + dwdz*tauzz
        call this%filter(diss_L, diss_S, this%x_bc, this%y_bc, this%z_bc)

        dudx => duidxj_ff(:,:,:,1); dudy => duidxj_ff(:,:,:,2); dudz => duidxj_ff(:,:,:,3);
        dvdx => duidxj_ff(:,:,:,4); dvdy => duidxj_ff(:,:,:,5); dvdz => duidxj_ff(:,:,:,6);
        dwdx => duidxj_ff(:,:,:,7); dwdy => duidxj_ff(:,:,:,8); dwdz => duidxj_ff(:,:,:,9);
        
        ! Get large scale dissipation in diss_L
        call this%filter(tauxx, tmp, this%x_bc, this%y_bc, this%z_bc)
        diss_L = dudx*tmp
        call this%filter(tauxy, tmp,-this%x_bc,-this%y_bc, this%z_bc)
        diss_L = diss_L + (dudy+dvdx)*tmp
        call this%filter(tauxz, tmp,-this%x_bc, this%y_bc,-this%z_bc)
        diss_L = diss_L + (dudz+dwdx)*tmp
        call this%filter(tauyy, tmp, this%x_bc, this%y_bc, this%z_bc)
        diss_L = diss_L + dvdy*tmp
        call this%filter(tauyz, tmp, this%x_bc,-this%y_bc,-this%z_bc)
        diss_L = diss_L + (dvdz+dwdy)*tmp
        call this%filter(tauzz, tmp, this%x_bc, this%y_bc, this%z_bc)
        diss_L = diss_L + dwdz*tmp

        ! Get small scale dissipation
        diss_S = diss_S - diss_L

    end subroutine

    subroutine tke_budget(this, rho, u, v, w, p, tauij, KE_L_old, KE_S_old, KE_L_prefilter, KE_L_postfilter, KE_S_prefilter, KE_S_postfilter, tsim, dt)
        use RKCoeffs,   only: RK45_steps
        use reductions, only: P_SUM
        class(scaleDecomposition),                                                       intent(inout) :: this
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3)),            intent(in)    :: rho, u, v, w, p
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3),6),          intent(in)    :: tauij
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3)),            intent(in)    :: KE_L_old, KE_S_old
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3),RK45_steps), intent(in)    :: KE_L_prefilter,  KE_S_prefilter
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3),RK45_steps), intent(in)    :: KE_L_postfilter, KE_S_postfilter
        real(rkind),                                                                     intent(in)    :: tsim, dt

        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: rho_f, u_ff, v_ff, w_ff, p_f
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: KE_L, KE_S
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: ddt_KE_L, ddt_KE_S
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: production, baropycnal
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: pdil_L, pdil_S
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: diss_L, diss_S
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: diss_L_num, diss_S_num

        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3),3) :: mass_flux
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3),6) :: turb_stress
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3),9) :: duidxj_ff

        real(rkind) :: KE_L_int, KE_S_int, ddt_KE_L_int, ddt_KE_S_int
        real(rkind) :: production_int, baropycnal_int, pdil_L_int, pdil_S_int
        real(rkind) :: diss_L_int, diss_S_int, diss_L_num_int, diss_S_num_int

        ! Get filtered density, velocities, pressure
        call this%filter(rho, rho_f, this%x_bc, this%y_bc, this%z_bc)
        call this%favre_filter(rho, rho_f, u, u_ff,-this%x_bc, this%y_bc, this%z_bc)
        call this%favre_filter(rho, rho_f, v, v_ff, this%x_bc,-this%y_bc, this%z_bc)
        call this%favre_filter(rho, rho_f, w, w_ff, this%x_bc, this%y_bc,-this%z_bc)
        call this%filter(p, p_f, this%x_bc, this%y_bc, this%z_bc)

        ! Get large and small scale kinetic energies and turbulent stresses
        call this%get_kinetic_energies_turb_stress(rho, rho_f, u, u_ff, v, v_ff, w, w_ff, KE_L, KE_S, turb_stress)
        KE_L_int = P_SUM(KE_L)*this%dx*this%dy*this%dz
        KE_S_int = P_SUM(KE_S)*this%dx*this%dy*this%dz

        ! Get KE rate of change
        ddt_KE_L = (KE_L - KE_L_old)/dt
        ddt_KE_S = (KE_S - KE_S_old)/dt
        ddt_KE_L_int = P_SUM(ddt_KE_L)*this%dx*this%dy*this%dz
        ddt_KE_S_int = P_SUM(ddt_KE_S)*this%dx*this%dy*this%dz

        ! Get numerical dissipation (from filtering)
        diss_L_num = sum(KE_L_prefilter - KE_L_postfilter, 4)/dt
        diss_S_num = sum(KE_S_prefilter - KE_S_postfilter, 4)/dt
        diss_L_num_int = P_SUM(diss_L_num)*this%dx*this%dy*this%dz
        diss_S_num_int = P_SUM(diss_S_num)*this%dx*this%dy*this%dz

        ! Get production
        call this%get_production(rho_f, u_ff, v_ff, w_ff, turb_stress, duidxj_ff, production)
        production_int = P_SUM(production)*this%dx*this%dy*this%dz

        ! Get sub-filter mass flux
        call this%get_mass_flux(rho, rho_f, u, v, w, mass_flux)

        ! Get baropycnal work
        call this%get_baropycnal(rho_f, p_f, mass_flux, baropycnal)
        baropycnal_int = P_SUM(baropycnal)*this%dx*this%dy*this%dz

        ! Get large and small scale pressure-dilatation
        call this%get_pressure_dilatation(p, p_f, u, v, w, pdil_L, pdil_S)
        pdil_L_int = P_SUM(pdil_L)*this%dx*this%dy*this%dz
        pdil_S_int = P_SUM(pdil_S)*this%dx*this%dy*this%dz

        ! Get large and small scale dissipation
        call this%get_dissipation(u, v, w, duidxj_ff, tauij, diss_L, diss_S)
        diss_L_int = P_SUM(diss_L)*this%dx*this%dy*this%dz
        diss_S_int = P_SUM(diss_S)*this%dx*this%dy*this%dz

        ! Write out data to output file
        call this%tke_viz%start_viz(tsim)
        call this%tke_viz%write_attribute(1, [dt], 'dt', '/')

        call this%tke_viz%write_attribute(1, [KE_L_int], 'KE_L_int', '/')
        call this%tke_viz%write_attribute(1, [KE_S_int], 'KE_S_int', '/')

        call this%tke_viz%write_attribute(1, [ddt_KE_L_int], 'KE_L_rate_int', '/')
        call this%tke_viz%write_attribute(1, [ddt_KE_S_int], 'KE_S_rate_int', '/')

        call this%tke_viz%write_attribute(1, [production_int], 'production_int', '/')
        call this%tke_viz%write_attribute(1, [baropycnal_int], 'baropycnal_int', '/')

        call this%tke_viz%write_attribute(1, [pdil_L_int], 'pdil_L_int', '/')
        call this%tke_viz%write_attribute(1, [pdil_S_int], 'pdil_S_int', '/')

        call this%tke_viz%write_attribute(1, [diss_L_int],     'diss_L_int', '/')
        call this%tke_viz%write_attribute(1, [diss_S_int],     'diss_S_int', '/')
        call this%tke_viz%write_attribute(1, [diss_L_num_int], 'diss_L_num_int', '/')
        call this%tke_viz%write_attribute(1, [diss_S_num_int], 'diss_S_num_int', '/')

        call this%tke_viz%write_variable(rho_f,'rho_f')
        call this%tke_viz%write_variable(u_ff, 'u_ff')
        call this%tke_viz%write_variable(v_ff, 'v_ff')
        call this%tke_viz%write_variable(w_ff, 'w_ff')
        call this%tke_viz%write_variable(p_f,  'p_f')

        call this%tke_viz%write_variable(KE_L, 'KE_L')
        call this%tke_viz%write_variable(KE_S, 'KE_S')

        call this%tke_viz%write_variable(turb_stress(:,:,:,1),  'turb_stress_11')
        call this%tke_viz%write_variable(turb_stress(:,:,:,2),  'turb_stress_12')
        call this%tke_viz%write_variable(turb_stress(:,:,:,3),  'turb_stress_13')
        call this%tke_viz%write_variable(turb_stress(:,:,:,4),  'turb_stress_22')
        call this%tke_viz%write_variable(turb_stress(:,:,:,5),  'turb_stress_23')
        call this%tke_viz%write_variable(turb_stress(:,:,:,6),  'turb_stress_33')

        call this%tke_viz%write_variable(mass_flux(:,:,:,1), 'mass_flux_x')
        call this%tke_viz%write_variable(mass_flux(:,:,:,2), 'mass_flux_y')
        call this%tke_viz%write_variable(mass_flux(:,:,:,3), 'mass_flux_z')

        call this%tke_viz%write_variable(ddt_KE_L,        'KE_L_rate')
        call this%tke_viz%write_variable(ddt_KE_S,        'KE_S_rate')
        call this%tke_viz%write_variable(production,      'production')
        call this%tke_viz%write_variable(pdil_L,          'pdil_L')
        call this%tke_viz%write_variable(pdil_S,          'pdil_S')
        call this%tke_viz%write_variable(baropycnal,      'baropycnal')
        call this%tke_viz%write_variable(diss_L,          'diss_L')
        call this%tke_viz%write_variable(diss_S,          'diss_S')
        call this%tke_viz%write_variable(diss_L_num,      'diss_L_num')
        call this%tke_viz%write_variable(diss_S_num,      'diss_S_num')

        call this%tke_viz%end_viz()

    end subroutine

    subroutine get_rhoPsi(this, rho, Ys, rhoPsi)
        class(scaleDecomposition),                                                       intent(in)  :: this
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)         ), intent(in)  :: rho
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3), this%ns), intent(in)  :: Ys
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)         ), intent(out) :: rhoPsi

        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ) :: rho_f, Ys_ff

        ! Get filtered density, massfraction
        call this%filter(rho, rho_f, this%x_bc, this%y_bc, this%z_bc)
        call this%favre_filter(rho, rho_f, Ys(:,:,:,this%species_id), Ys_ff, this%x_bc, this%y_bc, this%z_bc)

        ! Large scale mix metric
        rhoPsi = rho_f * Ys_ff * (one - Ys_ff)

    end subroutine

    subroutine get_rhoPsi_turb_flux_generation(this, rho, rho_f, u, u_ff, v, v_ff, w, w_ff, Ys, Ys_ff, rhoPsi, turb_flux, grad_Ys_ff, Pi_Psi)
        use operators, only: gradient
        class(scaleDecomposition),                                                 intent(in)  :: this
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(in)  :: rho, u, v, w, Ys
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(in)  :: rho_f, u_ff, v_ff, w_ff, Ys_ff
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(out) :: rhoPsi, Pi_Psi
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3), 3), intent(out) :: turb_flux, grad_Ys_ff

        ! Large scale mix metric
        rhoPsi = rho_f * Ys_ff * (one - Ys_ff)

        ! Get turbulent flux
        call this%favre_filter(rho, rho_f, u*Ys, turb_flux(:,:,:,1),-this%x_bc, this%y_bc, this%z_bc)
        call this%favre_filter(rho, rho_f, v*Ys, turb_flux(:,:,:,2), this%x_bc,-this%y_bc, this%z_bc)
        call this%favre_filter(rho, rho_f, w*Ys, turb_flux(:,:,:,3), this%x_bc, this%y_bc,-this%z_bc)

        turb_flux(:,:,:,1) = turb_flux(:,:,:,1) - u_ff*Ys_ff
        turb_flux(:,:,:,2) = turb_flux(:,:,:,2) - v_ff*Ys_ff
        turb_flux(:,:,:,3) = turb_flux(:,:,:,3) - w_ff*Ys_ff

        ! Favre-filtered massfraction gradient
        call gradient(this%gp, this%der, Ys_ff, grad_Ys_ff(:,:,:,1), grad_Ys_ff(:,:,:,2), grad_Ys_ff(:,:,:,3), this%x_bc, this%y_bc, this%z_bc)

        ! Get turbulent generation term
        Pi_Psi = - two * rho_f * (turb_flux(:,:,:,1)*grad_Ys_ff(:,:,:,1) + turb_flux(:,:,:,2)*grad_Ys_ff(:,:,:,2) + turb_flux(:,:,:,3)*grad_Ys_ff(:,:,:,3))

    end subroutine

    subroutine get_diffusive_flux_generation(this, rho, rho_f, Ys, diff, grad_Ys_ff, J_ff, G_Psi)
        use operators, only: gradient
        class(scaleDecomposition),                                                 intent(in)  :: this
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(in)  :: rho, rho_f
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(in)  :: Ys, diff
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3), 3), intent(in)  :: grad_Ys_ff
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3), 3), intent(out) :: J_ff
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(out) :: G_Psi

        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)) :: dYsdx, dYsdy, dYsdz

        ! Full massfraction gradient
        call gradient(this%gp, this%der, Ys, dYsdx, dYsdy, dYsdz, this%x_bc, this%y_bc, this%z_bc)

        dYsdx = diff * dYsdx
        call this%favre_filter(rho, rho_f, dYsdx, J_ff(:,:,:,1),-this%x_bc, this%y_bc, this%z_bc)

        dYsdy = diff * dYsdy
        call this%favre_filter(rho, rho_f, dYsdy, J_ff(:,:,:,2), this%x_bc,-this%y_bc, this%z_bc)

        dYsdz = diff * dYsdz
        call this%favre_filter(rho, rho_f, dYsdz, J_ff(:,:,:,3), this%x_bc, this%y_bc,-this%z_bc)

        ! Large scale generation term
        G_psi = two * rho_f * (J_ff(:,:,:,1)*grad_Ys_ff(:,:,:,1) + J_ff(:,:,:,2)*grad_Ys_ff(:,:,:,2) + J_ff(:,:,:,3)*grad_Ys_ff(:,:,:,3))

    end subroutine

    subroutine get_transport_divergence(this, rho_f, u_ff, v_ff, w_ff, Ys_ff, rhoPsi, turb_flux, J_ff, &
                                        advection, turb_transport, diff_transport)
        use operators, only: divergence
        class(scaleDecomposition),                                                 intent(in)  :: this
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(in)  :: rho_f, u_ff, v_ff, w_ff
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(in)  :: Ys_ff, rhoPsi
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3), 3), intent(in)  :: J_ff, turb_flux
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(out) :: advection
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(out) :: turb_transport
        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)   ), intent(out) :: diff_transport

        real(rkind), dimension(this%gp%ysz(1), this%gp%ysz(2), this%gp%ysz(3)) :: T_x, T_y, T_z

        ! Large scale advection fluxes
        T_x = rhoPsi * u_ff
        T_y = rhoPsi * v_ff
        T_z = rhoPsi * w_ff

        ! Divergence of transport fluxes
        call divergence( this%gp, this%der, T_x, T_y, T_z, advection, this%x_bc, this%y_bc, this%z_bc )

        ! Turbulent transport fluxes
        T_x = rho_f * (one - two*Ys_ff) * (-turb_flux(:,:,:,1))
        T_y = rho_f * (one - two*Ys_ff) * (-turb_flux(:,:,:,2))
        T_z = rho_f * (one - two*Ys_ff) * (-turb_flux(:,:,:,3))

        ! Divergence of transport fluxes
        call divergence( this%gp, this%der, T_x, T_y, T_z, turb_transport, this%x_bc, this%y_bc, this%z_bc )

        ! Diffusive transport fluxes
        T_x = rho_f * (one - two*Ys_ff) * J_ff(:,:,:,1)
        T_y = rho_f * (one - two*Ys_ff) * J_ff(:,:,:,2)
        T_z = rho_f * (one - two*Ys_ff) * J_ff(:,:,:,3)

        ! Divergence of transport fluxes
        call divergence( this%gp, this%der, T_x, T_y, T_z, diff_transport, this%x_bc, this%y_bc, this%z_bc )

    end subroutine

    subroutine mix_budget(this, rho, u, v, w, Ys, diff, rhoPsi_old, rhoPsi_prefilter, rhoPsi_postfilter, tsim, dt)
        use RKCoeffs,   only: RK45_steps
        use reductions, only: P_SUM
        class(scaleDecomposition),                                                       intent(inout) :: this
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3)),            intent(in)    :: rho, u, v, w
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3),this%ns),    intent(in)    :: Ys, diff
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3)),            intent(in)    :: rhoPsi_old
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3),RK45_steps), intent(in)    :: rhoPsi_prefilter
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3),RK45_steps), intent(in)    :: rhoPsi_postfilter
        real(rkind),                                                                     intent(in)    :: tsim, dt

        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: rho_f, u_ff, v_ff, w_ff, Ys_ff
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: rhoPsi
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: ddt_rhoPsi
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: turb_generation
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: diff_generation
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: num_generation
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: advection
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: turb_transport
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3))   :: diff_transport

        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3),3) :: diff_flux
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3),3) :: turb_flux
        real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3),3) :: grad_Ys_ff

        real(rkind) :: rhoPsi_int, ddt_rhoPsi_int
        real(rkind) :: turb_generation_int, diff_generation_int
        real(rkind) :: num_generation_int

        ! Get filtered density, velocities, massfraction
        call this%filter(rho, rho_f, this%x_bc, this%y_bc, this%z_bc)
        call this%favre_filter(rho, rho_f, u, u_ff,-this%x_bc, this%y_bc, this%z_bc)
        call this%favre_filter(rho, rho_f, v, v_ff, this%x_bc,-this%y_bc, this%z_bc)
        call this%favre_filter(rho, rho_f, w, w_ff, this%x_bc, this%y_bc,-this%z_bc)
        call this%favre_filter(rho, rho_f, Ys(:,:,:,this%species_id), Ys_ff, this%x_bc, this%y_bc, this%z_bc)

        ! Get large scale mix, turbulent fluxes and turbulent generation
        call this%get_rhoPsi_turb_flux_generation(rho, rho_f, u, u_ff, v, v_ff, w, w_ff, &
                               Ys(:,:,:,this%species_id), Ys_ff, rhoPsi, turb_flux, grad_Ys_ff, turb_generation)
        rhoPsi_int = P_SUM(rhoPsi)*this%dx*this%dy*this%dz
        turb_generation_int = P_SUM(turb_generation)*this%dx*this%dy*this%dz

        ! Get rhoPsi rate of change
        ddt_rhoPsi = (rhoPsi - rhoPsi_old)/dt
        ddt_rhoPsi_int = P_SUM(ddt_rhoPsi)*this%dx*this%dy*this%dz

        ! Get numerical generation (from filtering)
        num_generation = sum(rhoPsi_prefilter - rhoPsi_postfilter, 4)/dt
        num_generation_int = P_SUM(num_generation)*this%dx*this%dy*this%dz

        ! Get large and small scale dissipation
        call this%get_diffusive_flux_generation(rho, rho_f, Ys(:,:,:,this%species_id), &
                                diff(:,:,:,this%species_id), grad_Ys_ff, diff_flux, diff_generation)
        diff_generation_int = P_SUM(diff_generation)*this%dx*this%dy*this%dz

        call this%get_transport_divergence(rho_f, u_ff, v_ff, w_ff, Ys_ff, rhoPsi, turb_flux, diff_flux, &
                                           advection, turb_transport, diff_transport)

        ! Write out data to output file
        call this%mix_viz%start_viz(tsim)
        call this%mix_viz%write_attribute(1, [dt], 'dt', '/')

        call this%mix_viz%write_attribute(1, [rhoPsi_int],          'rhoPsi_int',          '/')
        call this%mix_viz%write_attribute(1, [ddt_rhoPsi_int],      'rhoPsi_rate_int',     '/')
        call this%mix_viz%write_attribute(1, [turb_generation_int], 'turb_generation_int', '/')
        call this%mix_viz%write_attribute(1, [diff_generation_int], 'diff_generation_int', '/')
        call this%mix_viz%write_attribute(1, [num_generation_int],  'num_generation_int',  '/')

        call this%mix_viz%write_variable(rho_f, 'rho_f')
        call this%mix_viz%write_variable(Ys_ff, 'Ys_ff')

        call this%mix_viz%write_variable(rhoPsi, 'rhoPsi')

        call this%mix_viz%write_variable(turb_flux(:,:,:,1), 'turb_flux_x')
        call this%mix_viz%write_variable(turb_flux(:,:,:,2), 'turb_flux_y')
        call this%mix_viz%write_variable(turb_flux(:,:,:,3), 'turb_flux_z')

        call this%mix_viz%write_variable(diff_flux(:,:,:,1), 'diff_flux_x')
        call this%mix_viz%write_variable(diff_flux(:,:,:,2), 'diff_flux_y')
        call this%mix_viz%write_variable(diff_flux(:,:,:,3), 'diff_flux_z')

        call this%mix_viz%write_variable(ddt_rhoPsi,      'rhoPsi_rate')
        call this%mix_viz%write_variable(advection,       'advection')
        call this%mix_viz%write_variable(turb_transport,  'turb_transport')
        call this%mix_viz%write_variable(diff_transport,  'diff_transport')
        call this%mix_viz%write_variable(turb_generation, 'turb_generation')
        call this%mix_viz%write_variable(diff_generation, 'diff_generation')
        call this%mix_viz%write_variable(num_generation,  'num_generation')

        call this%mix_viz%end_viz()

    end subroutine

end module
