program test_TKEBudget
    use mpi
    use decomp_2d,           only: nrank, decomp_info
    use kind_parameters,     only: rkind, clen
    use constants,           only: zero
    use miranda_reader_mod,  only: miranda_reader
    use TKEBudgetMod,        only: tkeBudget
    use DerivativesMod,      only: derivatives
    use exits,               only: message, GracefulExit

    implicit none

    type(miranda_reader)  :: mir
    type(tkeBudget)       :: budget
    type(derivatives)     :: der

    character(len=clen) :: jobdir, outputdir
    integer :: prow = 0, pcol = 0

    integer, dimension(2) :: x_bc = [0, 1]
    integer, dimension(2) :: y_bc = [1, 1]
    integer, dimension(2) :: z_bc = [0, 0]

    logical :: periodicx = .false., periodicy = .false., periodicz = .true.
    logical, dimension(3) :: averaging_directions

    real(rkind) :: dt = real(1.0D-4, rkind)

    real(rkind), dimension(:,:,:),   allocatable :: tke_old, tke_prefilter
    real(rkind), dimension(:,:,:,:), allocatable :: tauij

    integer :: ierr, step

    call MPI_Init(ierr)

    if( iargc() .LT. 1 ) then
        call GracefulExit("Usage: "//NEW_LINE('A')//"    mpiexec -n 8 ./test_TKEBudget <jobdir>", 1729)
    end if


    call GETARG(1,jobdir)
    ! Initialize miranda_reader object
    call mir%init(jobdir, prow, pcol, periodicx, periodicy, periodicz)

    if (nrank == 0) then
        print *, "prow = ", prow, ", pcol = ", pcol
        print *, "Jobdir is "//trim(jobdir)
    end if

    ! Read in the grid
    call mir%read_grid()

    ! Initialize the derivative object
    call der%init(                            mir%gp, &
                        mir%dx,    mir%dy,    mir%dz, &
                     periodicx, periodicy, periodicz, &
                        "cd10",    "cd10",    "cd10" )

    ! Initialize budget object
    outputdir = "./outputs"
    budget = tkeBudget(mir%gp, der, mir%mesh, mir%dx, mir%dy, mir%dz, [periodicx, periodicy, periodicz], outputdir, x_bc, y_bc, z_bc, .true.)

    allocate( tke_old      (mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)) )
    allocate( tke_prefilter(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3)) )
    allocate(         tauij(mir%gp%ysz(1), mir%gp%ysz(2), mir%gp%ysz(3), 6) )

    tke_old = zero
    tke_prefilter = zero

    ! Read in data for time step 0
    do step = 0, mir%nsteps-1
        call message("Processing step ", step)
        call mir%read_data(step)

        call get_tauij(mir%gp, der, mir%u, mir%v, mir%w, mir%mu, mir%bulk, tauij, x_bc, y_bc, z_bc)
        call budget%tke_budget(mir%rho, mir%u, mir%v, mir%w, mir%p, tauij, tke_old, tke_prefilter, step*dt, dt)
    end do

    deallocate( tke_old )
    deallocate( tke_prefilter )
    deallocate(   tauij )

    call mir%destroy()
    call MPI_Finalize(ierr)

contains

    subroutine get_tauij(gp,der,u,v,w,mu,beta,tauij,x_bc,y_bc,z_bc)
        use constants, only: two, three, four
        use operators, only: gradient
        type(decomp_info),                                                  intent(in)  :: gp
        type(derivatives),                                                  intent(in)  :: der
        real(rkind), dimension(gp%ysz(1), gp%ysz(2), gp%ysz(3)   ),         intent(in)  :: u, v, w, mu, beta
        real(rkind), dimension(gp%ysz(1), gp%ysz(2), gp%ysz(3), 6), target, intent(out) :: tauij
        integer,     dimension(2),                                          intent(in)  :: x_bc, y_bc, z_bc

        real(rkind), dimension(gp%ysz(1), gp%ysz(2), gp%ysz(3), 9), target :: duidxj
        real(rkind), dimension(:,:,:), pointer :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
        real(rkind), dimension(:,:,:), pointer :: tauxx, tauxy, tauxz, tauyy, tauyz, tauzz

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3);
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6);
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9);
        
        tauxx => tauij(:,:,:,1); tauxy => tauij(:,:,:,2); tauxz => tauij(:,:,:,3);
                                 tauyy => tauij(:,:,:,4); tauyz => tauij(:,:,:,5);
                                                          tauzz => tauij(:,:,:,6);

        call gradient(gp,der,u,dudx,dudy,dudz,-x_bc, y_bc, z_bc)
        call gradient(gp,der,v,dvdx,dvdy,dvdz, x_bc,-y_bc, z_bc)
        call gradient(gp,der,w,dwdx,dwdy,dwdz, x_bc, y_bc,-z_bc)

        tauxx = ( (four/three)*mu + beta )*dudx + (beta - (two/three)*mu)*(dvdy+dwdz)
        tauxy = mu * (dudy + dvdx)
        tauxz = mu * (dudz + dwdx)
        tauyy = ( (four/three)*mu + beta )*dvdy + (beta - (two/three)*mu)*(dudx+dwdz)
        tauyz = mu * (dvdz + dwdy)
        tauzz = ( (four/three)*mu + beta )*dwdz + (beta - (two/three)*mu)*(dudx+dvdy)
    end subroutine


end program 
