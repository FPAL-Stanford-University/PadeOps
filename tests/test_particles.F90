program test_particles
    use kind_parameters, only: rkind, clen
    use constants, only: pi, two, one, imi, zero, half, kappa
    use reductions, only: p_maxval, p_sum, p_minval 
    use timer, only: tic, toc
    use decomp_2d
    use decomp_2d_io
    use particlesMod, only: particles
    use exits, only: message
    !use mpi

    implicit none 

    type(particles), allocatable :: LagrangianParticles
    integer, parameter :: nx = 192, ny = 144, nz = 120
    character(len=clen) :: inputfile  = "/home/aryam/runs/particles/prec_c3_19600/input.in"
    character(len=clen) :: outputfile
    real(rkind), dimension(:,:,:), allocatable :: xG, yG, zG
    real(rkind), dimension(:,:,:), allocatable :: u, v, w
    real(rkind), dimension(:),     allocatable :: xline, yline, zline
    real(rkind), parameter :: Lx = 36.0d0, Ly = 12.0d0, Lz = 10.0d0
    real(rkind) :: dx, dy, dz, current_time, dt
    type(decomp_info) :: gp
    integer :: idx, ix1, iy1, iz1, ixn, iyn, izn, i, j, k, ierr, prow = 0, pcol = 0
    integer :: itstep, num_time_steps = 20, RunID = 1, tviz_particles
    character(len=clen) :: field_dir = "/home/aryam/runs/min_turbn/c3r4n3/prec_data_new/"
    character(len=20) :: runIDstr
   ! Range of timesteps to process
   integer :: tStart = 53451, tEnd = 55801
   integer :: tIn
    real(rkind), parameter :: c = 1.0d0   
    real(rkind), parameter :: zn = 0.1d0   
    write(runIDstr,'("Run",I2.2)') RunId

    !call MPI_Init(ierr)
    !call decomp_2d_init(nx, ny, nz, prow, pcol)
    !call get_decomp_info(gp)

    !allocate(xG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(yG(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    !allocate(zG(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(u(gp%xsz(1),gp%xsz(2),gp%xsz(3)))
    !allocate(v(gp%xsz(1),gp%xsz(2),gp%xsz(3))); allocate(w(gp%xsz(1),gp%xsz(2),gp%xsz(3)))

    allocate(xG(nx,ny,nz), yG(nx,ny,nz), zG(nx,ny,nz))
    allocate( u(nx,ny,nz),  v(nx,ny,nz),  w(nx,ny,nz))
    allocate(xline(nx), yline(ny), zline(nz))

    dx = Lx/real(nx,rkind); dy = Ly/real(ny,rkind); dz = Lz/real(nz,rkind)
    !ix1 = gp%xst(1); iy1 = gp%xst(2); iz1 = gp%xst(3)
    !ixn = gp%xen(1); iyn = gp%xen(2); izn = gp%xen(3)
    ix1 = 1; iy1 = 1; iz1 = 1
    do k=1,size(xG,3)
        do j=1,size(xG,2)
            do i=1,size(xG,1)
                xG(i,j,k) = real( ix1 + i - 1, rkind ) * dx
                yG(i,j,k) = real( iy1 + j - 1, rkind ) * dy
                zG(i,j,k) = real( iz1 + k - 1, rkind ) * dz + dz/two
            end do
        end do
    end do
    xG = xG - dx; yG = yG - dy; zG = zG - dz 

    xline = xG(:,1,1)
    yline = yG(1,:,1)
    zline = zG(1,1,:)

    !u = one 
    !v = zero 
    !w = zero 

    allocate(LagrangianParticles)
    call LagrangianParticles%init(inputfile, Lx, Ly, Lz, nx, ny, nz, u, v, w, xline, yline, zline, .true., .true., .false.)

    current_time = zero
    dt = 0.001d0
    num_time_steps = 10000
    tviz_particles = 1


        
        itstep = 0
        
        do tIn = tStart, tEnd
            itstep = itstep + 1
            current_time = current_time + dt
        
            ! Read velocity fields for this timestep
            call read_field(field_dir, "uVel", tIn, u)
            call read_field(field_dir, "vVel", tIn, v)
            call read_field(field_dir, "wVel", tIn, w)
        
            ! Update particles
            call LagrangianParticles%update(dt, u, v, w)
        
            call message(1,"Finished Time Step:", itstep)
        
            ! Write particle data (one file per particle, append a row each timestep)
            call LagrangianParticles%write_viz(itstep, current_time, trim(runIDstr))
        end do

    call LagrangianParticles%destroy()
 
    deallocate(LagrangianParticles)
    deallocate(xG, yG, zG, u, v, w, xline, yline, zline)
    !call MPI_Finalize(ierr)
!end program

contains

subroutine read_field(field_dir, component, timestep, field)
        implicit none
        character(len=*), intent(in) :: field_dir     ! directory where files are stored
        character(len=*), intent(in) :: component     ! "uVel", "vVel", or "wVel"
        integer, intent(in)          :: timestep      ! timestep number, e.g. 53334
        real(rkind), dimension(:,:,:), intent(out) :: field

        character(len=clen) :: fname
        integer :: unit

        ! Construct filename, e.g. Run34_uVel_t053334.out
        write(fname,"(A,'Run34_',A,'_t',I6.6,'.out')") trim(field_dir), trim(component), timestep

        ! Open binary file and read into array
        open(newunit=unit, file=fname, status='old', action='read', form='unformatted', access = 'stream')
        read(unit) field
        close(unit)

    end subroutine read_field
end program
