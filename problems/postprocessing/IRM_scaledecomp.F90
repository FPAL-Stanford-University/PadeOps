module IRM_scaledecomp_mod
    use mpi
    use kind_parameters,    only: rkind, clen, mpirkind, mpickind
    use constants,          only: zero, eps, half, one, three
    use miranda_reader_mod, only: miranda_reader
    use DerivativesMod,     only: derivatives
    use FiltersMod,         only: filters
    use decomp_2d,          only: nrank, nproc, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y, &
                                  decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize
    use exits,              only: message
    use io_hdf5_stuff,      only: io_hdf5

    implicit none

    type(miranda_reader)  :: mir
    type(io_hdf5)         :: viz
    type(derivatives)     :: der
    type(filters)         :: gfil

    character(len=clen)   :: inputdir, outputdir                                                 ! Input and output dirs
    integer               :: prow, pcol                                                          ! Procs in x and z
    logical               :: periodicx = .FALSE., periodicy = .FALSE., periodicz = .FALSE.       ! Periodic?
    character(len=4)      :: derivative_x = 'cd10', derivative_y = 'cd10', derivative_z = 'cd10' ! Derivative method to use

    integer               :: num_filter                                                          ! Number of times to filter (n^2 => k_max = 2*pi / (n*4*dx) )

    logical               :: writeviz                                                            ! Write out HDF5 viz files?
    integer, dimension(6) :: vizsteps                                                            ! Steps at which to write out viz files

    integer, dimension(2) :: x_bc, y_bc, z_bc                                                    ! Symmetric BC? (Not used if periodic)

    integer               :: XY_COMM
    integer               :: xyrank, xyproc

contains

    subroutine read_inputs(inputfile)
        character(len=clen), intent(in) :: inputfile
        integer :: iounit = 67

        integer :: x_bc1 = 0, x_bcn = 0
        integer :: y_bc1 = 0, y_bcn = 0
        integer :: z_bc1 = 0, z_bcn = 0

        namelist /INPUT/ inputdir, outputdir, &
                         periodicx, periodicy, periodicz, &
                         derivative_x, derivative_y, derivative_z, &
                         prow, pcol, writeviz, x_bc1, x_bcn, y_bc1, y_bcn, z_bc1, z_bcn, &
                         vizsteps, num_filter

        open(unit=iounit, file=trim(inputfile), form='FORMATTED')
        read(unit=iounit, NML=INPUT)
        close(iounit)

        x_bc = [x_bc1, x_bcn]
        y_bc = [y_bc1, y_bcn]
        z_bc = [z_bc1, z_bcn]

    end subroutine

    subroutine get_duidxj(u,v,w,duidxj)
        use operators, only: gradient
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)),            intent(in)  :: u, v, w
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 9), target, intent(out) :: duidxj

        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);

        call gradient(mir%gp,der,u,dudx,dudy,dudz,-x_bc, y_bc, z_bc)
        call gradient(mir%gp,der,v,dvdx,dvdy,dvdz, x_bc,-y_bc, z_bc)
        call gradient(mir%gp,der,w,dwdx,dwdy,dwdz, x_bc, y_bc,-z_bc)

    end subroutine

    subroutine get_tauij(duidxj,mu,beta,tauij)
        use constants, only: two, three, four
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 9), target, intent(in)  :: duidxj
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)   ),         intent(in)  :: mu, beta
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 6), target, intent(out) :: tauij
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: tauxx, tauxy, tauxz, tauyy, tauyz, tauzz

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        tauxx => tauij(:,:,:,1); tauxy => tauij(:,:,:,2); tauxz => tauij(:,:,:,3);
                                 tauyy => tauij(:,:,:,4); tauyz => tauij(:,:,:,5);
                                                          tauzz => tauij(:,:,:,6);

        tauxx = ( (four/three)*mu + beta )*dudx + (beta - (two/three)*mu)*(dvdy+dwdz)
        tauxy = mu * (dudy + dvdx)
        tauxz = mu * (dudz + dwdx)
        tauyy = ( (four/three)*mu + beta )*dvdy + (beta - (two/three)*mu)*(dudx+dwdz)
        tauyz = mu * (dvdz + dwdy)
        tauzz = ( (four/three)*mu + beta )*dwdz + (beta - (two/three)*mu)*(dudx+dvdy)
    end subroutine

    subroutine filter(fil, var, var_f, x_bc, y_bc, z_bc)
        use operators, only: filter3D
        type(filters),                                                       intent(in)  :: fil
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)), intent(in)  :: var
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)), intent(out) :: var_f
        integer,     dimension(2),                                           intent(in)  :: x_bc, y_bc, z_bc

        var_f = var
        call filter3D( mir%gp, fil, var_f, num_filter, x_bc, y_bc, z_bc )
    end subroutine

    subroutine favre_filter(fil, rho, rho_f, var, var_ff, x_bc, y_bc, z_bc)
        use operators, only: filter3D
        type(filters),                                                       intent(in)  :: fil
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)), intent(in)  :: rho, rho_f, var
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)), intent(out) :: var_ff
        integer,     dimension(2),                                           intent(in)  :: x_bc, y_bc, z_bc

        var_ff = rho*var
        call filter3D( mir%gp, fil, var_ff, num_filter, x_bc, y_bc, z_bc )
        var_ff = var_ff / rho_f
    end subroutine

    subroutine get_kinetic_energies_turb_stress(rho, rho_f, u, u_ff, v, v_ff, w, w_ff, KE_L, KE_S, turb_stress)
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)  ), intent(in)  :: rho, rho_f, u, u_ff, v, v_ff, w, w_ff
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)  ), intent(out) :: KE_L, KE_S
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3),6), intent(out) :: turb_stress

        ! Large scale kinetic energy
        KE_L = half * rho_f * (u_ff*u_ff + v_ff*v_ff + w_ff*w_ff)

        ! Turbulent stresses: tilde(u_i*u_j) - tilde(u_i)*tilde(u_j)
        call favre_filter(gfil, rho, rho_f, u*u, turb_stress(:,:,:,1), x_bc, y_bc, z_bc)
        turb_stress(:,:,:,1) = turb_stress(:,:,:,1) - u_ff*u_ff    ! ts_11 = tilde(u*u) - tilde(u)*tilde(u)

        call favre_filter(gfil, rho, rho_f, u*v, turb_stress(:,:,:,2),-x_bc,-y_bc, z_bc)
        turb_stress(:,:,:,2) = turb_stress(:,:,:,2) - u_ff*v_ff    ! ts_12 = tilde(u*v) - tilde(u)*tilde(v)

        call favre_filter(gfil, rho, rho_f, u*w, turb_stress(:,:,:,3),-x_bc, y_bc,-z_bc)
        turb_stress(:,:,:,3) = turb_stress(:,:,:,3) - u_ff*w_ff    ! ts_13 = tilde(u*w) - tilde(u)*tilde(w)

        call favre_filter(gfil, rho, rho_f, v*v, turb_stress(:,:,:,4), x_bc, y_bc, z_bc)
        turb_stress(:,:,:,4) = turb_stress(:,:,:,4) - v_ff*v_ff    ! ts_22 = tilde(v*v) - tilde(v)*tilde(v)

        call favre_filter(gfil, rho, rho_f, v*w, turb_stress(:,:,:,5), x_bc,-y_bc,-z_bc)
        turb_stress(:,:,:,5) = turb_stress(:,:,:,5) - v_ff*w_ff    ! ts_23 = tilde(v*w) - tilde(v)*tilde(w)

        call favre_filter(gfil, rho, rho_f, w*w, turb_stress(:,:,:,6), x_bc, y_bc, z_bc)
        turb_stress(:,:,:,6) = turb_stress(:,:,:,6) - w_ff*w_ff    ! ts_33 = tilde(w*w) - tilde(w)*tilde(w)

        ! Small scale kinetic energy
        KE_S = half * rho_f * (turb_stress(:,:,:,1) + turb_stress(:,:,:,4) + turb_stress(:,:,:,6))

    end subroutine

    subroutine get_production(rho_f, u_ff, v_ff, w_ff, turb_stress, duidxj_ff, production)
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)  ),         intent(in)  :: rho_f, u_ff, v_ff, w_ff
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3),6), target, intent(in)  :: turb_stress
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3),9), target, intent(out) :: duidxj_ff
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)  ),         intent(out) :: production

        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: TSxx, TSxy, TSxz, TSyy, TSyz, TSzz

        dudx => duidxj_ff(:,:,:,1); dudy => duidxj_ff(:,:,:,2); dudz => duidxj_ff(:,:,:,3);
        dvdx => duidxj_ff(:,:,:,4); dvdy => duidxj_ff(:,:,:,5); dvdz => duidxj_ff(:,:,:,6);
        dwdx => duidxj_ff(:,:,:,7); dwdy => duidxj_ff(:,:,:,8); dwdz => duidxj_ff(:,:,:,9);

        TSxx => turb_stress(:,:,:,1); TSxy => turb_stress(:,:,:,2); TSxz => turb_stress(:,:,:,3);
                                      TSyy => turb_stress(:,:,:,4); TSyz => turb_stress(:,:,:,5);
                                                                    TSzz => turb_stress(:,:,:,6);

        call get_duidxj(u_ff,v_ff,w_ff,duidxj_ff)

        production = dudx*TSxx + dudy*TSxy + dudz*TSxz &
                   + dvdx*TSxy + dvdy*TSyy + dvdz*TSyz &
                   + dwdx*TSxz + dwdy*TSyz + dwdz*TSzz

        production = -rho_f * production

    end subroutine

    subroutine get_mass_flux(rho, rho_f, u, v, w, mass_flux)
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)  ), intent(in)  :: rho, rho_f, u, v, w
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3),3), intent(out) :: mass_flux

        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)) :: vel_f

        call filter(gfil, rho*u, mass_flux(:,:,:,1),-x_bc, y_bc, z_bc)
        call filter(gfil, u, vel_f,-x_bc, y_bc, z_bc)
        mass_flux(:,:,:,1) = mass_flux(:,:,:,1) - rho_f*vel_f ! filter(rho*u) - filter(rho)*filter(u)
        
        call filter(gfil, rho*v, mass_flux(:,:,:,2), x_bc,-y_bc, z_bc)
        call filter(gfil, v, vel_f, x_bc,-y_bc, z_bc)
        mass_flux(:,:,:,2) = mass_flux(:,:,:,2) - rho_f*vel_f ! filter(rho*v) - filter(rho)*filter(v)
        
        call filter(gfil, rho*w, mass_flux(:,:,:,3), x_bc, y_bc,-z_bc)
        call filter(gfil, w, vel_f, x_bc, y_bc,-z_bc)
        mass_flux(:,:,:,3) = mass_flux(:,:,:,3) - rho_f*vel_f ! filter(rho*w) - filter(rho)*filter(w)
        
    end subroutine

    subroutine get_baropycnal(rho_f, p_f, mass_flux, baropycnal)
        use operators, only: gradient
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)  ), intent(in)  :: rho_f, p_f
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3),3), intent(in)  :: mass_flux
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)  ), intent(out) :: baropycnal

        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)) :: dpdx, dpdy, dpdz

        call gradient(mir%gp,der,p_f,dpdx,dpdy,dpdz,x_bc,y_bc,z_bc)

        baropycnal = ( dpdx*mass_flux(:,:,:,1) + dpdy*mass_flux(:,:,:,2) + dpdz*mass_flux(:,:,:,3) ) / rho_f

    end subroutine

    subroutine get_pressure_dilatation(p, p_f, u, v, w, pdil_L, pdil_S)
        use operators, only: divergence
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)  ), intent(in)  :: p, p_f, u, v, w
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)  ), intent(out) :: pdil_L, pdil_S

        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)) :: u_f, v_f, w_f, divu

        call filter(gfil, u, u_f,-x_bc, y_bc, z_bc)
        call filter(gfil, v, v_f, x_bc,-y_bc, z_bc)
        call filter(gfil, w, w_f, x_bc, y_bc,-z_bc)

        call divergence( mir%gp, der, u_f, v_f, w_f, divu,-x_bc,-y_bc,-z_bc )
        pdil_L = p_f*divu ! filter(p) * divergence( filter(velocity) )

        call divergence( mir%gp, der, u, v, w, divu,-x_bc,-y_bc,-z_bc )
        divu = p*divu
        call filter(gfil, divu, pdil_S, x_bc, y_bc, z_bc)
        pdil_S = pdil_S - pdil_L    ! filter( p * divergence(velocity) ) - filter(p)*divergence( filter(velocity) )

    end subroutine

    subroutine get_dissipation(u, v, w, duidxj_ff, tauij, diss_L, diss_S)
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)  ),         intent(in)    :: u, v, w
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3),9), target, intent(in)    :: duidxj_ff
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3),9), target, intent(inout) :: tauij
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)  ),         intent(out)   :: diss_L, diss_S

        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3),9), target :: duidxj
        real(rkind), dimension(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)  )         :: tmp

        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: tauxx, tauxy, tauxz, tauyy, tauyz, tauzz

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        tauxx => tauij(:,:,:,1); tauxy => tauij(:,:,:,2); tauxz => tauij(:,:,:,3);
                                 tauyy => tauij(:,:,:,4); tauyz => tauij(:,:,:,5);
                                                          tauzz => tauij(:,:,:,6);

        call get_duidxj(u,v,w,duidxj)
        call get_tauij(duidxj,mir%mu,mir%bulk,tauij)

        diss_L = dudx*tauxx + dudy*tauxy + dudz*tauxz &
               + dvdx*tauxy + dvdy*tauyy + dvdz*tauyz &
               + dwdx*tauxz + dwdy*tauyz + dwdz*tauzz
        call filter(gfil, diss_L, diss_S, x_bc, y_bc, z_bc)

        dudx => duidxj_ff(:,:,:,1); dudy => duidxj_ff(:,:,:,2); dudz => duidxj_ff(:,:,:,3);
        dvdx => duidxj_ff(:,:,:,4); dvdy => duidxj_ff(:,:,:,5); dvdz => duidxj_ff(:,:,:,6);
        dwdx => duidxj_ff(:,:,:,7); dwdy => duidxj_ff(:,:,:,8); dwdz => duidxj_ff(:,:,:,9);
        
        call filter(gfil, tauxx, tmp, x_bc, y_bc, z_bc)
        diss_L = dudx*tmp
        call filter(gfil, tauxy, tmp,-x_bc,-y_bc, z_bc)
        diss_L = diss_L + (dudy+dvdx)*tmp
        call filter(gfil, tauxz, tmp,-x_bc, y_bc,-z_bc)
        diss_L = diss_L + (dudz+dwdx)*tmp
        call filter(gfil, tauyy, tmp, x_bc, y_bc, z_bc)
        diss_L = diss_L + dvdy*tmp
        call filter(gfil, tauyz, tmp, x_bc,-y_bc,-z_bc)
        diss_L = diss_L + (dvdz+dwdy)*tmp
        call filter(gfil, tauzz, tmp, x_bc, y_bc, z_bc)
        diss_L = diss_L + dwdz*tmp

        diss_S = diss_S - diss_L

    end subroutine

end module

program IRM_scaledecomp
    use mpi
    use kind_parameters,    only: rkind, clen
    use constants,          only: zero
    use miranda_reader_mod, only: miranda_reader
    use DerivativesMod,     only: derivatives
    use FiltersMod,         only: filters
    use gridtools,          only: alloc_buffs
    use operators,          only: curl, gradient, divergence
    use reductions,         only: P_AVGZ, P_MAXVAL, P_MINVAL, P_SUM
    use timer,              only: tic, toc
    use exits,              only: message, GracefulExit

    use IRM_scaledecomp_mod

    implicit none

    character(len=clen) :: inputfile

    real(rkind), dimension(:,:,:,:), allocatable, target :: buffer

    real(rkind), dimension(:,:,:),   pointer :: rho_f, u_ff, v_ff, w_ff, p_f
    real(rkind), dimension(:,:,:),   pointer :: mu, ktc
    real(rkind), dimension(:,:,:),   pointer :: KE_L, production, baropycnal, pdil_L, diss_L
    real(rkind), dimension(:,:,:),   pointer :: KE_S, pdil_S, diss_S
    real(rkind), dimension(:,:,:,:), pointer :: Diff, turb_stress, mass_flux

    real(rkind), dimension(:,:,:,:), allocatable :: duidxj_ff, tauij

    real(rkind) :: KE_L_int, KE_S_int, production_int, baropycnal_int, pdil_L_int, pdil_S_int, diss_L_int, diss_S_int

    logical :: vizstep = .false.

    character(len=clen) :: time_message
    character(len=clen) :: outputfile
    integer :: iounit = 93

    integer :: ierr, istep, step, i

    !======================================================================================================!
    ! Note about usage: All 3D arrays are in Y decomposition.
    !                   All 2D arrays are also in Y decomposition (effectively 1D decomp)
    !                   All spectra are in Z, averaged over a region and are stored on every processor
    !                   It is assumed in some places that the number of species is 2
    !======================================================================================================!

    call MPI_Init(ierr)

    call tic()

    if( command_argument_count() .LT. 1 ) then
        call GracefulExit("Usage: "//NEW_LINE('A')//"    mpiexec -n 8 ./IRM_scaledecomp <input file>", 1729)
    end if

    call get_command_argument(1,inputfile)
    call read_inputs(inputfile)

    ! Initialize miranda_reader object
    call mir%init(inputdir, prow, pcol, periodicx, periodicy, periodicz)

    call message("Jobdir is " // adjustl(trim(inputdir)) )
    call message("Initialized the Miranda object" )
    
    call message("Reading in the grid" )
    call mir%read_grid()
    
    call message("Initializing the derivative object" )
    call der%init(                                      mir%gp, &
                          mir%dx,        mir%dy,        mir%dz, &
                   mir%periodicx, mir%periodicy, mir%periodicz, &
                    derivative_x,  derivative_y,  derivative_z  )

    ! Split MPI_COMM_WORLD into individual XY plane communicators
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, mir%gp%yst(3), nrank, XY_COMM, ierr)
    call MPI_COMM_RANK(XY_COMM, xyrank, ierr)
    call MPI_COMM_SIZE(XY_COMM, xyproc, ierr)

    call message("Initializing the gaussian filter object" )
    call gfil%init(                                     mir%gp, &
                   mir%periodicx, mir%periodicy, mir%periodicz, &
                      'gaussian',    'gaussian',    'gaussian'  )
    
    ! Allocate 3D buffer and associate pointers for convenience
    call alloc_buffs(buffer, 24+mir%ns, 'y', mir%gp); buffer = zero
    i = 1
    mu          => buffer(:,:,:,i);            i = i+1
    ktc         => buffer(:,:,:,i);            i = i+1
    Diff        => buffer(:,:,:,i:i+mir%ns-1); i = i+mir%ns
    rho_f       => buffer(:,:,:,i);            i = i+1
    u_ff        => buffer(:,:,:,i);            i = i+1
    v_ff        => buffer(:,:,:,i);            i = i+1
    w_ff        => buffer(:,:,:,i);            i = i+1
    p_f         => buffer(:,:,:,i);            i = i+1
    mass_flux   => buffer(:,:,:,i:i+2);        i = i+3
    turb_stress => buffer(:,:,:,i:i+5);        i = i+6
    KE_L        => buffer(:,:,:,i);            i = i+1
    KE_S        => buffer(:,:,:,i);            i = i+1
    production  => buffer(:,:,:,i);            i = i+1
    baropycnal  => buffer(:,:,:,i);            i = i+1
    pdil_L      => buffer(:,:,:,i);            i = i+1
    pdil_S      => buffer(:,:,:,i);            i = i+1
    diss_L      => buffer(:,:,:,i);            i = i+1
    diss_S      => buffer(:,:,:,i);            i = i+1

    allocate( duidxj_ff(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 9) )
    allocate( tauij    (mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 6) )

    ! Initialize visualization stuff
    if ( writeviz ) then
        call message("Initializing the HDF5 I/O object" )
        write(outputfile,'(A,I4.4)') 'IRM_scaledecomp_', num_filter
        call viz%init(mpi_comm_world, mir%gp, 'y', outputdir, outputfile, read_only=.false.)
        call viz%write_coords(mir%mesh) ! Write coordinates to file
    end if
    
    ! Open file for scalar outputs
    write(outputfile,'(2A,I4.4,A)') trim(outputdir), "/post_scalar_", num_filter,".dat"
    if(nrank == 0) then
        open(unit=iounit, file=trim(outputfile), form='FORMATTED', status='REPLACE')
        write(iounit,'(A,A25,8A26)') "#", "Time", "KE_L_int", "KE_S_int", "production_int", "baropycnal_int", "pdil_L_int", "pdil_S_int", "diss_L_int", "diss_S_int"
    end if

    call toc("Time to initialize everything: ")

    ! Loop through viz dumps
    do istep = 1,size(vizsteps,1)
        step = vizsteps(istep)

        if ( ( step < 0 ) .or. (step > mir%nsteps-1) ) then
            if (nrank == 0) then
                print '(A,I0,A)', "WARNING: Step ", step, " not available in the Miranda viz files. Check your input file for correctness."
            end if
        end if

        call tic()  ! Start timer

        call message(0,"Reading in data for step",step)

        ! Read in data for this viz dump
        call mir%read_data(step)
        
        ! Get filtered density, velocities, pressure
        call filter(gfil, mir%rho, rho_f, x_bc, y_bc, z_bc)
        call favre_filter(gfil, mir%rho, rho_f, mir%u, u_ff,-x_bc, y_bc, z_bc)
        call favre_filter(gfil, mir%rho, rho_f, mir%v, v_ff, x_bc,-y_bc, z_bc)
        call favre_filter(gfil, mir%rho, rho_f, mir%w, w_ff, x_bc, y_bc,-z_bc)
        call filter(gfil, mir%p, p_f, x_bc, y_bc, z_bc)

        ! Get large and small scale kinetic energies and turbulent stresses
        call get_kinetic_energies_turb_stress(mir%rho, rho_f, mir%u, u_ff, mir%v, v_ff, mir%w, w_ff, KE_L, KE_S, turb_stress)
        KE_L_int = P_SUM(KE_L)*mir%dx*mir%dy*mir%dz
        KE_S_int = P_SUM(KE_S)*mir%dx*mir%dy*mir%dz

        ! Get production
        call get_production(rho_f, u_ff, v_ff, w_ff, turb_stress, duidxj_ff, production)
        production_int = P_SUM(production)*mir%dx*mir%dy*mir%dz

        ! Get sub-filter mass flux
        call get_mass_flux(mir%rho, rho_f, mir%u, mir%v, mir%w, mass_flux)

        ! Get baropycnal work
        call get_baropycnal(rho_f, p_f, mass_flux, baropycnal)
        baropycnal_int = P_SUM(baropycnal)*mir%dx*mir%dy*mir%dz

        ! Get large and small scale pressure-dilatation
        call get_pressure_dilatation(mir%p, p_f, mir%u, mir%v, mir%w, pdil_L, pdil_S)
        pdil_L_int = P_SUM(pdil_L)*mir%dx*mir%dy*mir%dz
        pdil_S_int = P_SUM(pdil_S)*mir%dx*mir%dy*mir%dz

        ! Get large and small scale dissipation
        call get_dissipation(mir%u, mir%v, mir%w, duidxj_ff, tauij, diss_L, diss_S)
        diss_L_int = P_SUM(diss_L)*mir%dx*mir%dy*mir%dz
        diss_S_int = P_SUM(diss_S)*mir%dx*mir%dy*mir%dz

        call message(1,"Writing scalars to file")
        if(nrank == 0) then
            write(iounit,'(9ES26.16)') real(step*1.0d-4,rkind), KE_L_int, KE_S_int, production_int, baropycnal_int, pdil_L_int, pdil_S_int, diss_L_int, diss_S_int
        end if

        vizstep = .false.
        do i = 1,size(vizsteps,1)
            if (step == vizsteps(i)) vizstep = .true.
        end do
        ! Write out visualization files
        if ( writeviz .and. vizstep ) then
            call message(1,"Writing visualization files")
            call viz%start_viz(real(step*1.0d-4,rkind))

            call viz%write_variable(rho_f,      'rho_f')
            call viz%write_variable(u_ff,       'u_ff')
            call viz%write_variable(v_ff,       'v_ff')
            call viz%write_variable(w_ff,       'w_ff')
            call viz%write_variable(p_f,        'p_f')
            call viz%write_variable(KE_L,       'KE_L')
            call viz%write_variable(production, 'production')
            call viz%write_variable(baropycnal, 'baropycnal')
            call viz%write_variable(pdil_L,     'pdil_L')
            call viz%write_variable(pdil_S,     'pdil_S')
            call viz%write_variable(diss_L,     'diss_L')
            call viz%write_variable(diss_S,     'diss_S')

            call viz%write_variable(mass_flux(:,:,:,1), 'mass_flux_x')
            call viz%write_variable(mass_flux(:,:,:,2), 'mass_flux_y')
            call viz%write_variable(mass_flux(:,:,:,3), 'mass_flux_z')

            call viz%write_variable(turb_stress(:,:,:,1),  'turb_stress_11')
            call viz%write_variable(turb_stress(:,:,:,2),  'turb_stress_12')
            call viz%write_variable(turb_stress(:,:,:,3),  'turb_stress_13')
            call viz%write_variable(turb_stress(:,:,:,4),  'turb_stress_22')
            call viz%write_variable(turb_stress(:,:,:,5),  'turb_stress_23')
            call viz%write_variable(turb_stress(:,:,:,6),  'turb_stress_33')

            call viz%end_viz()
        end if

        write(time_message,'(A,I4,A)') "Time to postprocess step ", step, " :"
        call toc(trim(time_message))
    end do

    ! Close file for scalar inputs
    if(nrank == 0) close(iounit)

    ! Destroy all variables and exit cleanly
    if (allocated( buffer ))    deallocate( buffer   )
    if (allocated( duidxj_ff )) deallocate( duidxj_ff )
    if (allocated( tauij ))     deallocate( tauij )

    call der%destroy()
    call gfil%destroy()
    call mir%destroy()
    if ( writeviz ) then
        call viz%destroy()
    end if
    call MPI_Finalize(ierr)

end program 
