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
    character(len=clen) :: inputfile  = "/home/aryam/runs/particles/second/input.in"
    character(len=clen) :: outputfile
    real(rkind), dimension(:,:,:), allocatable :: xG, yG, zG
    real(rkind), dimension(:,:,:), allocatable :: u, v, w
    real(rkind), dimension(:),     allocatable :: xline, yline, zline
    real(rkind), parameter :: Lx = 32.0d0, Ly = 12.0d0, Lz = 10.0d0
    real(rkind) :: dx, dy, dz, current_time, dt
    type(decomp_info) :: gp
    integer :: idx, ix1, iy1, iz1, ixn, iyn, izn, i, j, k, ierr, prow = 0, pcol = 0
    integer :: itstep, num_time_steps = 20, RunID = 1, tviz_particles

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

   ! u = one 
    v = zero 
    w = zero  
    
        ! --- Define log profile constants ---
    real(rkind), parameter :: c = 1.0d0   ! scaling constant
    real(rkind), parameter :: d = 0.1d0   ! roughness length (must be < min(z))

    ! --- Assign velocity fields ---
    do k=1,nz
        do j=1,ny
            do i=1,nx
                if (zG(i,j,k) > d) then
                    u(i,j,k) = c * log(zG(i,j,k)/d)
                else
                    u(i,j,k) = 0.0_rkind   ! avoid log(0) or negative
                end if
                !v(i,j,k) = zero
               ! w(i,j,k) = zero
            end do
        end do
    end do

    allocate(LagrangianParticles)
    call LagrangianParticles%init(inputfile, Lx, Ly, Lz, nx, ny, nz, u, v, w, xline, yline, zline, .true., .true., .false.)

    current_time = zero
    dt = 0.1d0
    num_time_steps = 20
    tviz_particles = 1

    do itstep = 1, num_time_steps

        current_time = current_time + dt
        call LagrangianParticles%update(dt, u, v, w)
 
        !call mpi_barrier(mpi_comm_world, ierr)
        !call tic()
        !do idx = 1,num_turbines
        !    call hawts_Rot(idx)%get_RHS(u, v, w, rhs1, rhsv, rhsw, inst_val)
        !end do 
        !call mpi_barrier(mpi_comm_world, ierr)
        !call toc()
        !call decomp_2d_write_one(1,rhs1,"temp_T1.bin", gp)

        call message(1,"Finished Time Step:", itstep)

        if(mod(itstep, tviz_particles) == 0) then
            write(outputfile,"(A3,I2.2,A12,I6.6,A4)") "Run",RunId,"_particles_t",itstep,".out"
            call LagrangianParticles%write_viz(current_time, outputfile)
        endif
        !print '(i5,1x,11(e19.12,1x))', iv, velocities(iv), inst_val, C_thrust, C_torque

    enddo

    call LagrangianParticles%destroy()
 
    deallocate(LagrangianParticles)
    deallocate(xG, yG, zG, u, v, w, xline, yline, zline)
    !call MPI_Finalize(ierr)
end program
