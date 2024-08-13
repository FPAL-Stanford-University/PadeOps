module Vortexring_data
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero, half, one, two, four, pi, imi
    use FiltersMod,       only: filters
    use decomp_2d,        only: decomp_info, nrank
    use basic_io,         only: read_2d_ascii 
    use reductions,       only: P_MAXVAL,P_MINVAL
    use exits,            only: message, message_min_max
    use mpi

    implicit none
    !!!! NOTE: Make sure to update this data according to the problem !!!!
    integer     :: ns     = 1
    real(rkind) :: Lx     = 0.768_rkind
    real(rkind) :: Ly     = 1.024_rkind
    real(rkind) :: Lz     = 1.024_rkind
    real(rkind) :: Dia    = 0.064_rkind
    real(rkind) :: Pr     = 0.7_rkind 
    real(rkind) :: Sc     = 1.0_rkind  
    real(rkind) :: gam    = 1.4_rkind
    real(rkind) :: Rgas   = 287.053_rkind
    real(rkind) :: mu_ref = 0.00001849_rkind 
    real(rkind) :: p_ref  = 101300.0_rkind
    real(rkind) :: rho_ref= 1.225_rkind
    real(rkind) :: x1, y1, z1
    real(rkind) :: xn, yn, zn
    logical     :: periodicx = .false., periodicy = .false., periodicz = .false. 
    logical     :: add_noise = .false.
    ! Gaussian filter for sponge
    type(filters) :: mygfil

contains

    subroutine scale_uvw_TI(u_pert, v_pert, w_pert, TI_reqd, U_scale)
    use kind_parameters,  only: rkind
    use constants,        only: two, zero
    use reductions,       only: P_MAXVAL,P_MINVAL, P_MEAN, P_SUM

    real(rkind), dimension(:,:),   intent(inout) :: u_pert, v_pert, w_pert
    real(rkind),                   intent(in)    :: TI_reqd, U_scale

    real(rkind) :: TI_current, TI_scale, sum_pert
    integer :: i, j, k

    sum_pert = zero
    do k = 1, size(u_pert,2)
       do j = 1, size(u_pert,1)
          sum_pert = sum_pert + (u_pert(j,k)**two + v_pert(j,k)**two + w_pert(j,k)**two)
       end do
    end do
    sum_pert = sum_pert/(size(u_pert,1)+size(u_pert,2))
    TI_current = sqrt(sum_pert/3)/U_scale
    
    !print*, ">>>>> Current TI =", TI_current
    TI_scale = TI_reqd/TI_current

    u_pert = TI_scale*u_pert
    v_pert = TI_scale*v_pert
    w_pert = TI_scale*w_pert

    end subroutine

    subroutine get_corr_field(u_pert, bx, by, bz, rx, ny, nz)
    use kind_parameters,  only: rkind

    real(rkind), dimension(:,:),   intent(inout):: u_pert
    real(rkind), dimension(:),     intent(in)   :: bx, by, bz
    real(rkind), dimension(:,:,:), intent(in)   :: rx
    integer,                        intent(in)   :: ny, nz
    
    real(rkind) :: b3d
    integer :: i, j, k, ii, jj, kk

    do k = 1, nz
       do j = 1, ny
           do kk = 1, size(bz)
               do jj = 1, size(by)
                   do ii = 1, size(bx)
                       b3d = bx(ii) * by(jj) * bz(kk)
                       u_pert(j,k) = u_pert(j,k) + b3d * rx(ii,jj+j,kk+k)
                   end do
               end do
           end do
       end do
    end do

    end subroutine

    subroutine get_filter_coeff(bx, nfx, lsx)
    use kind_parameters,  only: rkind
    use constants,        only: pi, zero, half, one, two, three, four, five, six, seven, eight

    real(rkind), dimension(:), intent(inout) :: bx
    integer,                   intent(in) :: nfx, lsx
    
    integer :: i, ii
    real(rkind) :: sum_bx, func

    sum_bx = zero
    
    ! Calculate sum of squares of the function values
    do i = 1, size(bx)
        ii = i - 1 - nfx
        func = exp((-pi*(ii**two))/(two*lsx**two))
        sum_bx = sum_bx + (func**two)
    end do
    
    sum_bx = sqrt(sum_bx)
    
    ! Normalize the function values
    do i = 1, size(bx)
        ii = i - 1 - nfx
        func = exp((-pi*(ii**two))/(two*lsx**two))
        bx(i) = func / sum_bx
    end do    

    end subroutine

    subroutine get_final_fields(y1d, z1d, dia, ny, nz, glob_npic, nfx, nfy, nfz, lsx, lsy, lsz, U_avg, &
                                                                 u2d_full, v2d_full, w2d_full, time_step, newTimestep)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one, two, three, four, five, six, seven, eight

    real(rkind), dimension(:),   intent(in)    :: y1d, z1d
    integer,                     intent(in)    :: ny, nz, glob_npic, nfx, nfy, nfz, lsx, lsy, lsz, time_step
    real(rkind),                 intent(in)    :: U_avg, dia
    real(rkind), dimension(:,:), intent(inout) :: u2d_full, v2d_full, w2d_full
    logical,                     intent(in)    :: newTimeStep

    integer :: i, j, k, idx,  ny_in, nz_in
    integer :: mpi_ierr, iseed, idx_2d_to_1d
    real(rkind) :: U_scale, TI_reqd, rad 
    integer, dimension(:), allocatable            :: fact
    real(rkind), dimension(:),        allocatable :: coeff_bx, coeff_by, coeff_bz, u1d, v1d, w1d
    real(rkind), dimension(:,:),      allocatable :: u_noise_2d, v_noise_2d, w_noise_2d
    real(rkind), dimension(:,:,:),    allocatable :: rand_x, rand_y, rand_z

    allocate(fact(1))
    fact = 1
    do i = 2, glob_npic
       if (mod(glob_npic, i) == 0) then
           fact = reshape((/ fact, i /), (/ size(fact) + 1 /))
       end if
    end do
    ny_in = fact((size(fact) / 2)+1)
    nz_in = glob_npic/ny_in

    if (time_step == 0 .and. newTimestep) then
         print*, ">>> Checking no.of points inside diameter = ", glob_npic
         print*, ">>> Checking U_avg within diameter = ", U_avg
         print*, ">>> Checking Factors of ", glob_npic," = " , fact
         print*, ">>> Checking noise matrix dimensions: ny_npic= ", ny_in, " nz_npic= ", nz_in
    end if

    allocate(u_noise_2d(ny_in, nz_in)); allocate(v_noise_2d(ny_in, nz_in)); allocate(w_noise_2d(ny_in, nz_in))
    allocate(rand_x(2*nfx+1, ny_in+2*nfy+1, nz_in+2*nfz+1))
    allocate(rand_y(2*nfx+1, ny_in+2*nfy+1, nz_in+2*nfz+1))
    allocate(rand_z(2*nfx+1, ny_in+2*nfy+1, nz_in+2*nfz+1))
    allocate(coeff_bx(2*nfx+1)); allocate(coeff_by(2*nfy+1)); allocate(coeff_bz(2*nfz+1))
    allocate(u1d(glob_npic));  allocate(v1d(glob_npic));  allocate(w1d(glob_npic))

    u_noise_2d = zero; v_noise_2d = zero; w_noise_2d = zero
    call random_seed(iseed)                 ! Initialize random number generator seed
    call random_number(rand_x)
    call random_number(rand_y)
    call random_number(rand_z)

    call get_filter_coeff(coeff_bx, nfx, lsx)
    call get_filter_coeff(coeff_by, nfy, lsy)
    call get_filter_coeff(coeff_bz, nfz, lsz)

    call get_corr_field(u_noise_2d, coeff_bx, coeff_by, coeff_bz, rand_x, ny_in, nz_in)
    call get_corr_field(v_noise_2d, coeff_bx, coeff_by, coeff_bz, rand_y, ny_in, nz_in)
    call get_corr_field(w_noise_2d, coeff_bx, coeff_by, coeff_bz, rand_z, ny_in, nz_in)
 
    ! Scale the correlated random velocity field to get correct TI
    TI_reqd = 0.05_rkind*U_avg; U_scale = 0.707_rkind
    call scale_uvw_TI(u_noise_2d, v_noise_2d, w_noise_2d, TI_reqd, U_scale)

    idx_2d_to_1d = 0
    do k = 1, nz_in
        do j = 1, ny_in
            idx_2d_to_1d = idx_2d_to_1d + 1
            u1d(idx_2d_to_1d) = u_noise_2d(j, k)
            v1d(idx_2d_to_1d) = v_noise_2d(j, k)
            w1d(idx_2d_to_1d) = w_noise_2d(j, k)
        end do
    end do

    idx = 0
    do k = 1, nz
       do j = 1, ny
          rad = sqrt(y1d(j)**2 + z1d(k)**2)
          if (rad .le. dia/2) then
             idx = idx + 1
             u2d_full(j,k)  = u1d(idx)
             v2d_full(j,k)  = v1d(idx)
             w2d_full(j,k)  = w1d(idx)
          endif
       end do
    end do

    deallocate(coeff_bx, coeff_by, coeff_bz)
    deallocate(rand_x, rand_y, rand_z)
    deallocate(u_noise_2d, v_noise_2d, w_noise_2d)
    deallocate(u1d, v1d, w1d)
    deallocate(fact)
    end subroutine

    subroutine find_noise(decomp, Lx, Ly, Lz, dia, x, y, z, u, v, w, u_noise, v_noise, w_noise, time_step, newTimeStep)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one, two, three, four, five, six, seven, eight
    use decomp_2d,        only: decomp_info, nrank, transpose_y_to_z, transpose_z_to_y
    use operators,        only: filter3D
    use reductions,       only: P_MAXVAL,P_MINVAL, P_MEAN, P_SUM

    type(decomp_info),               intent(in)    :: decomp
    real(rkind), dimension(:,:,:),   intent(in)    :: x, y, z
    real(rkind),                     intent(in)    :: Lx, Ly, Lz, dia
    real(rkind), dimension(:,:,:),   intent(in)    :: u,v,w
    real(rkind), dimension(:,:),     intent(inout) :: u_noise, v_noise, w_noise
    integer,                         intent(in)    :: time_step
    logical,                         intent(in)    :: newTimeStep
    
    integer :: i, j, k, ii, jj, kk, nx, ny, nz, nxl, nyl, nzl, lsx, lsy, lsz, nfx, nfy, nfz
    integer :: mpi_ierr, loc_npic, glob_npic, idx_rank, idx_tot
    real(rkind) :: loc_sum, glob_sum, rad, U_avg
    real(rkind), dimension(:),         allocatable :: y1d, z1d
    real(rkind), dimension(:,:,:),     allocatable :: ztmp

    ! Choose Length scale
    lsx = 1; lsy = 1; lsz = 1
    ! Choose filter width, should be atleast twice of length scale
    nfy = 2*lsy;  nfz = 2*lsz; nfx = nfz 
    ! Global size
    nx  = decomp%xsz(1);    ny  = decomp%ysz(2) ;  nz  = decomp%zsz(3)
    ! Local domain sizes
    nxl = decomp%ysz(1);    nyl = decomp%ysz(2);   nzl = decomp%ysz(3)
    ! Store y and z values to check points within circle
    allocate(y1d(ny)); allocate(z1d(nz)); allocate(ztmp(decomp%zsz(1), decomp%zsz(2),decomp%zsz(3)))
    y1d = y(1,:,1)
    call transpose_y_to_z(z,ztmp,decomp)   ! Decomposition in z
    z1d = ztmp(1,1,:)   

    ! Find number of points inside circle(glob_npic) and average vel(U_avg) within circle
    loc_sum = 0; loc_npic = 0
    do k = 1,nzl
       do j = 1,nyl
          rad = sqrt(y(1,j,1)**2 + z(1,1,k)**2)
          if (rad .le. dia/2) then
             loc_npic = loc_npic + 1
             loc_sum  = loc_sum + u(1,j,k)
          endif
       end do
    end do
    call mpi_allreduce(loc_sum, glob_sum, 1, mpirkind, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
    call mpi_allreduce(loc_npic, glob_npic, 1, mpi_integer, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)
    glob_npic = glob_npic/(nx/nxl)
    U_avg = glob_sum/glob_npic

    !print*, ">> Rank = ", nrank, ">> local sum = ",loc_sum, ">> Global sum = ", glob_sum
    !print*, ">> Rank = ", nrank, ">> local points = ",loc_npic, ">> Global points = ", glob_npic
    !print*, ">>> U_avg = ", U_avg


    ! Find correlated random velocity field within  circle
    if (nrank == 0) then
       call get_final_fields(y1d, z1d, dia, ny, nz, glob_npic, nfx, nfy, nfz, lsx, lsy, lsz, &
                                        U_avg, u_noise, v_noise, w_noise, time_step, newTimeStep) 
    end if

    call MPI_Bcast(u_noise, ny*nz, mpirkind, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_Bcast(v_noise, ny*nz, mpirkind, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_Bcast(w_noise, ny*nz, mpirkind, 0, MPI_COMM_WORLD, mpi_ierr)

    deallocate(y1d, z1d)
    deallocate(ztmp)
    end subroutine

    subroutine sponge_z(decomp, mygfil, z, Lz, u, v, w, p, rho, x_bc, y_bc, z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one, two, three, four, five, six, seven, eight
    use decomp_2d,        only: decomp_info, nrank
    use operators,        only: filter3D

    type(decomp_info),               intent(in)    :: decomp
    type(filters),                   intent(in)    :: mygfil
    real(rkind), dimension(:,:,:),   intent(in)    :: z
    real(rkind),                     intent(in)    :: Lz
    real(rkind), dimension(:,:,:),   intent(inout) :: u,v,w,p,rho
    integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc

    integer :: i, j, k
    real(rkind) :: dx, dy, dz, filpt, thickT
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: dumT, dumF
    character(len=clen) :: outputfile

    dz = Lz/real(decomp%zsz(3)-1,rkind)
    filpt = 0.03_rkind/dz 
    thickT = real(1.0D0, rkind)
        
    ! Gussian Filter for start
    do i=1,decomp%ysz(3)
       dumT(:,:,i)=half*(one-tanh( (real( decomp%yst(3) - 1 + i - 1, rkind)-filpt) / thickT ))
    end do
            
    dumF = u
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    u = u + dumT*(dumF-u) 

    dumF = v
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    v = v + dumT*(dumF-v)

    dumF = w
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    w = w + dumT*(dumF-w)

    dumF = p
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    p = p + dumT*(dumF-p)

    dumF = rho
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    rho = rho + dumT*(dumF-rho)

    ! Gaussian Filter for end
    do i=1,decomp%ysz(3)
        dumT(:,:,i)=half*(one-tanh( (real(decomp%zsz(3)- (decomp%yst(3) - 1 + i - 1), rkind)-filpt) / thickT ))
    end do

    !write(outputfile, '(a,i3.3,a)') 'dumT_', nrank, '.dat'
    !open(10,file=outputfile,status='unknown')
    !do i=1,decomp%ysz(3)
    !   write(10,'(2(e19.12),1x)') z(1,1,i), dumT(1,1,i)
    !end do
    !close(10)

    dumF = u
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    u = u + dumT*(dumF-u) 

    dumF = v
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    v = v + dumT*(dumF-v)

    dumF = w
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    w = w + dumT*(dumF-w)

    dumF = p
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    p = p + dumT*(dumF-p)

    dumF = rho
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    rho = rho + dumT*(dumF-rho)

    end subroutine


    subroutine sponge_y(decomp, mygfil, y, Ly, u, v, w, p, rho, x_bc, y_bc, z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one, two, three, four, five, six, seven, eight
    use decomp_2d,        only: decomp_info, nrank
    use operators,        only: filter3D

    type(decomp_info),               intent(in)    :: decomp
    type(filters),                   intent(in)    :: mygfil
    real(rkind), dimension(:,:,:),   intent(in)    :: y
    real(rkind),                     intent(in)    :: Ly
    real(rkind), dimension(:,:,:),   intent(inout) :: u,v,w,p,rho
    integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc

    integer :: i, j, k
    real(rkind) :: dx, dy, dz, filpt, thickT
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: dumT, dumF
    character(len=clen) :: outputfile

    dy = Ly/real(decomp%ysz(2)-1,rkind)
    filpt = 0.03_rkind/dy 
    thickT = real(1.0D0, rkind)
        
    ! Gussian Filter for bottom
    do i=1,decomp%ysz(2)
       dumT(:,i,:)=half*(one-tanh( (real( decomp%yst(2) - 1 + i - 1, rkind)-filpt) / thickT ))
    end do
            
    dumF = u
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    u = u + dumT*(dumF-u) 

    dumF = v
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    v = v + dumT*(dumF-v)

    dumF = w
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    w = w + dumT*(dumF-w)

    dumF = p
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    p = p + dumT*(dumF-p)

    dumF = rho
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    rho = rho + dumT*(dumF-rho)

    ! Gaussian Filter for top
    do i=1,decomp%ysz(2)
       dumT(:,i,:)=half*(one-tanh( (real(decomp%ysz(2)- (decomp%yst(2) - 1 + i - 1), rkind)-filpt) / thickT ))
    end do

    dumF = u
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    u = u + dumT*(dumF-u) 

    dumF = v
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    v = v + dumT*(dumF-v)

    dumF = w
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    w = w + dumT*(dumF-w)

    dumF = p
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    p = p + dumT*(dumF-p)

    dumF = rho
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    rho = rho + dumT*(dumF-rho)
   
    !if (nrank==0) then
    !write(outputfile, '(a,i3.3,a)') 'dump_y_', nrank, '.dat'
    !open(10,file=outputfile,status='unknown')
    !do i=1,decomp%ysz(2)
    !   write(10,'(2(e19.12),1x)') y(1,i,1), dumT(1,i,1)
    !end do
    !close(10)
    !endif

    end subroutine


    subroutine sponge_x(decomp, mygfil, x, Lx, u, v, w, p, rho, x_bc, y_bc, z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one, two, three, four, five, six, seven, eight
    use decomp_2d,        only: decomp_info, nrank
    use operators,        only: filter3D

    type(decomp_info),               intent(in)    :: decomp
    type(filters),                   intent(in)    :: mygfil
    real(rkind), dimension(:,:,:),   intent(in)    :: x
    real(rkind),                     intent(in)    :: Lx
    real(rkind), dimension(:,:,:),   intent(inout) :: u,v,w,p,rho
    integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc

    integer :: i, j, k
    real(rkind) :: dx, dy, dz, filpt, thickT
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: dumT, dumF
    character(len=clen) :: outputfile

    dx = Lx/real(decomp%xsz(1)-1,rkind)
    filpt = 0.06_rkind/dx 
    thickT = real(2.0D0, rkind)

    ! Gaussian Filter for right side of domain 
    do i=1,decomp%ysz(1)
       dumT(i,:,:)=half*(one-tanh( (real(decomp%xsz(1)- (decomp%yst(1) - 1 + i - 1), rkind)-filpt) / thickT ))
    end do

    !!! To check whether dumT is calculted correctly !!!
    !write(outputfile, '(a,i3.3,a)') 'dumT_', nrank, '.dat'
    !open(10,file=outputfile,status='unknown')
    !do i=1,decomp%ysz(1)
    !   write(10,'(2(e19.12),1x)') x(i,1,1), dumT(i,1,1)
    !end do
    !close(10)

    dumF = u
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    u = u + dumT*(dumF-u) 

    dumF = v
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    v = v + dumT*(dumF-v)

    dumF = w
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    w = w + dumT*(dumF-w)

    dumF = p
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    p = p + dumT*(dumF-p)

    dumF = rho
    call filter3D(decomp,mygfil,dumF,4,x_bc,y_bc,z_bc)
    rho = rho + dumT*(dumF-rho)

    end subroutine

    subroutine get_fitted_fields(tsim, P0, rho0, U0)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one, two
    use decomp_2d,        only: decomp_info, nrank
    use exits,            only: GracefulExit, message, nancheck

    real(rkind),                     intent(in)    :: tsim
    real(rkind),                     intent(out)   :: P0, rho0, U0
    real(rkind) :: a0, w0
    real(rkind) :: a1, a2, a3, a4, a5, a6, a7, a8
    real(rkind) :: b1, b2, b3, b4, b5, b6, b7, b8
    real(rkind) :: c1, c2, c3, c4, c5, c6, c7, c8
    
   if (tsim .le. 0.00138_rkind ) then
    !!!! Pressure   !!!! PR = 03; DRL = 165mm
        a1 =   97320000.0_rkind; b1 =  -0.0002813_rkind;  c1 =   0.0001025_rkind
        a2 =   109700.0_rkind;   b2 =   0.00003148_rkind; c2 =   0.0002566_rkind
        a3 =   87450.0_rkind;    b3 =   0.0007718_rkind;  c3 =   0.0005342_rkind
        a4 =   22370.0_rkind;    b4 =   0.0004977_rkind;  c4 =   0.0001356_rkind
        a5 =   4179.0_rkind;     b5 =   0.0007091_rkind;  c5 =   0.00009197_rkind
        a6 =   22190.0_rkind;    b6 =   0.0003435_rkind;  c6 =   0.00008325_rkind
        a7 =   81480.0_rkind;    b7 =   0.001516_rkind;   c7 =   0.0004892_rkind

        P0 = a1*exp(-((tsim-b1)/c1)**two) + a2*exp(-((tsim-b2)/c2)**two) +  &
             a3*exp(-((tsim-b3)/c3)**two) + a4*exp(-((tsim-b4)/c4)**two) + &
             a5*exp(-((tsim-b5)/c5)**two) + a6*exp(-((tsim-b6)/c6)**two) + &
             a7*exp(-((tsim-b7)/c7)**two)

    !!!! Density
        a1 = 5055.0_rkind;  b1 = -0.0003915_rkind; c1 = 0.0001282_rkind
        a2 = 0.9431_rkind;  b2 = 0.00003233_rkind; c2 = 0.0002381_rkind
        a3 = 1.144_rkind;   b3 = 0.0007858_rkind;  c3 = 0.0007648_rkind
        a4 = 0.1478_rkind;  b4 = 0.0004976_rkind;  c4 = 0.0001295_rkind
        a5 = 0.02531_rkind; b5 = 0.0007127_rkind;  c5 = 0.00007627_rkind
        a6 = 0.1846_rkind;  b6 = 0.0003433_rkind;  c6 = 0.00008253_rkind
        a7 = 0.7336_rkind;  b7 = 0.001618_rkind;   c7 = 0.0004885_rkind

        rho0 = a1*exp(-((tsim-b1)/c1)**two) + a2*exp(-((tsim-b2)/c2)**two) +  &
             a3*exp(-((tsim-b3)/c3)**two) + a4*exp(-((tsim-b4)/c4)**two) + &
             a5*exp(-((tsim-b5)/c5)**two) + a6*exp(-((tsim-b6)/c6)**two) + &
             a7*exp(-((tsim-b7)/c7)**two)

    !!!! Velocity
        a1 = -117.7_rkind;  b1 = 0.00007236_rkind;  c1 = 0.0001633_rkind
        a2 = -6.428_rkind;  b2 = 0.0002879_rkind;   c2 = 0.00002214_rkind
        a3 = -9.864_rkind;  b3 = 0.0003476_rkind;   c3 = 0.00005368_rkind
        a4 = -270000_rkind; b4 = -0.0009989_rkind;  c4 = 0.0003654_rkind
        a5 = zero;          b5 = -0.002467_rkind;   c5 = 0.00004653_rkind
        a6 = 2.411_rkind;   b6 = 0.0006228_rkind;   c6 = 0.00005292_rkind
        a7 = 653.9_rkind;   b7 = -0.001468_rkind;   c7 = 0.002015_rkind

        U0 = a1*exp(-((tsim-b1)/c1)**two) + a2*exp(-((tsim-b2)/c2)**two) +  &
             a3*exp(-((tsim-b3)/c3)**two) + a4*exp(-((tsim-b4)/c4)**two) + &
             a5*exp(-((tsim-b5)/c5)**two) + a6*exp(-((tsim-b6)/c6)**two) + &
             a7*exp(-((tsim-b7)/c7)**two)
    else if (tsim .gt. 0.00138_rkind .and. tsim .le. 0.0033828_rkind) then
    !!!! Pressure   !!!! PR = 10; DRL = 285mm
        a1 = 78580.0_rkind; b1 = 0.0033_rkind; c1 = 0.0009583_rkind
        a2 = 92070.0_rkind; b2 = 0.001702_rkind; c2 = 0.001258_rkind
        a3 = 27580.0_rkind; b3 = 0.0008947_rkind; c3 = 0.0005325_rkind
        a4 = 9775.0_rkind; b4 = 0.003481_rkind; c4 = 0.0001721_rkind
        a5 = 924.0_rkind; b5 = 0.002466_rkind; c5 = 0.000164_rkind

        P0 = a1*exp(-((tsim-b1)/c1)**two) + a2*exp(-((tsim-b2)/c2)**two) +  &
             a3*exp(-((tsim-b3)/c3)**two) + a4*exp(-((tsim-b4)/c4)**two) + &
             a5*exp(-((tsim-b5)/c5)**two)
!!!! Density
        a1 = 1.066_rkind; b1 = 0.003275_rkind; c1 = 0.001093_rkind
        a2 = 1.019_rkind; b2 = 0.001652_rkind; c2 = 0.001111_rkind
        a3 = 0.4825_rkind; b3 = 0.0007822_rkind; c3 = 0.0006249_rkind
        a4 = 0.09465_rkind; b4 = 0.003488_rkind; c4 = 0.0001782_rkind
        a5 = 0.00903_rkind; b5 = 0.002468_rkind; c5 = 0.0001766_rkind

        rho0 = a1*exp(-((tsim-b1)/c1)**two) + a2*exp(-((tsim-b2)/c2)**two) +  &
             a3*exp(-((tsim-b3)/c3)**two) + a4*exp(-((tsim-b4)/c4)**two) + &
             a5*exp(-((tsim-b5)/c5)**two)

        !!!! Velocity
        a1 = -19.99_rkind; b1 = 0.002963_rkind; c1 = 0.0006406_rkind
        a2 = 2.35_rkind; b2 = 0.001719_rkind; c2 = 0.0002337_rkind
        a3 = 133.0_rkind; b3 = 0.0008287_rkind; c3 = 0.0008501_rkind
        a4 = 0.0_rkind; b4 = 0.01234_rkind; c4 = 0.0006641_rkind
        a5 = -6.723_rkind; b5 = 0.002301_rkind; c5 = 0.0002398_rkind

        U0 = a1*exp(-((tsim-b1)/c1)**two) + a2*exp(-((tsim-b2)/c2)**two) +  &
             a3*exp(-((tsim-b3)/c3)**two) + a4*exp(-((tsim-b4)/c4)**two) + &
             a5*exp(-((tsim-b5)/c5)**two)
else if (tsim .gt. 0.0033828_rkind .and. tsim .le. 0.0050828_rkind ) then
        !!!! Pressure   !!!! PR = 10; DRL = 285mm
        a1 = 78460.0_rkind; b1 = 0.004608_rkind; c1 = 0.0007341_rkind
        a2 = 95270.0_rkind; b2 = 0.003383_rkind; c2 = 0.0009368_rkind
        a3 = 66360.0_rkind; b3 = 0.005373_rkind; c3 = 0.0004261_rkind
        a4 = 8981.0_rkind; b4 = 0.004871_rkind; c4 = 0.0002552_rkind
        a5 = -1203.0_rkind; b5 = 0.003868_rkind; c5 = 0.0001478_rkind

        P0 = a1*exp(-((tsim-b1)/c1)**two) + a2*exp(-((tsim-b2)/c2)**two) +  &
             a3*exp(-((tsim-b3)/c3)**two) + a4*exp(-((tsim-b4)/c4)**two) + &
             a5*exp(-((tsim-b5)/c5)**two)

        !!!! Density
        a1 = 1.23_rkind; b1 = 0.004301_rkind; c1 = 0.002246_rkind
        a2 = 0.3355_rkind; b2 = 0.00336_rkind; c2 = 0.0003216_rkind
        a3 = 0.2435_rkind; b3 = 0.005468_rkind; c3 = 0.0004594_rkind
        a4 = 0.02893_rkind; b4 = 0.004866_rkind; c4 = 0.0002126_rkind
        a5 = -0.1634_rkind; b5 = 0.003405_rkind; c5 = 0.0002323_rkind

        rho0 = a1*exp(-((tsim-b1)/c1)**two) + a2*exp(-((tsim-b2)/c2)**two) +  &
             a3*exp(-((tsim-b3)/c3)**two) + a4*exp(-((tsim-b4)/c4)**two) + &
             a5*exp(-((tsim-b5)/c5)**two)

        !!!! Velocity
        a1 = -9.215_rkind; b1 = 0.004314_rkind; c1 = 0.000329_rkind
        a2 = -18.74_rkind; b2 = 0.003636_rkind; c2 = 0.000485_rkind
        a3 = -30.64_rkind; b3 = 0.005678_rkind; c3 = 0.0003735_rkind
        a4 = -4.112_rkind; b4 = 0.004914_rkind; c4 = 0.0001964_rkind
        a5 = 0.4059_rkind; b5 = 0.003794_rkind; c5 = 0.000003343_rkind
        U0 = a1*exp(-((tsim-b1)/c1)**two) + a2*exp(-((tsim-b2)/c2)**two) +  &
             a3*exp(-((tsim-b3)/c3)**two) + a4*exp(-((tsim-b4)/c4)**two) + &
             a5*exp(-((tsim-b5)/c5)**two)
    end if    


    end subroutine

    subroutine stretched_coordinates(decomp, y, eta, ymetric, ymetric_flag, param1, param2, param3, param4)
    use constants,        only: zero, half, one
    use decomp_2d,        only: decomp_info, nrank
    use exits,            only: GracefulExit, message, nancheck

    type(decomp_info),               intent(in)    :: decomp
    real(rkind), dimension(:,:,:),   intent(inout) :: y
    real(rkind), dimension(:,:,:  ), intent(inout) :: eta
    logical,                         intent(in   ) :: ymetric
    integer,                         intent(in)    :: ymetric_flag
    real(rkind),                     intent(in)    :: param1, param2, param3, param4
    integer     :: i,j,k
    real(rkind) :: yfocus, ytau, ystart, yh, alpha, beta
    real(rkind) :: yfocus_adj, num, den, BB, yuniform_adj, BB2

    ! concentrate towards the center -- Pletcher, Tannehill, Anderson
    ! (Section 5.6, Transformation 3, pg. 332) 
    if(ymetric_flag==1) then
       ! Concentrate towards centre
       yfocus = param1; ytau = param2; ystart = param3; yh = param4
       yfocus_adj = yfocus - ystart
       num = one + (yfocus_adj/yh) * (exp( ytau) - one)
       den = one + (yfocus_adj/yh) * (exp(-ytau) - one)
       BB  = half/ytau*log(num/den)
       do k = 1,decomp%ysz(3)
          do j = 1,decomp%ysz(2)
             do i = 1,decomp%ysz(1)
                  yuniform_adj = (eta(i,j,k) - ystart) !/ yh
                  num = sinh(ytau*BB)
                  y(i,j,k) = yfocus_adj * (one + sinh(ytau * (yuniform_adj/yh-BB))/num) + ystart
             end do
          end do
       end do
    elseif(ymetric_flag==2) then
       ! concentrate towards the two ends for alpha=0.5, at left end (final point) for alpha=0
       alpha = param1; beta = param2; ystart = param3; yh = param4
       BB   = (beta + 1) / (beta - 1)
       do k = 1,decomp%ysz(3)
          do j = 1,decomp%ysz(2)
             do i = 1,decomp%ysz(1)
                  yuniform_adj = (eta(i,j,k) - ystart) / yh
                  BB2 = BB ** ( (yuniform_adj-alpha) / (1-alpha) )
                  num = ((beta+2*alpha)*BB2 - beta + 2*alpha ) * yh
                  y(i,j,k) = num/( (2*alpha+1)*(1+BB2) )   + ystart
             end do
          end do
       end do
    elseif(ymetric_flag==3) then
       ! concentrate at right end, towards initial point
       alpha = param1; beta = param2; ystart = param3; yh = param4 !alpha not required
       BB   = (beta + 1) / (beta - 1)
       do k = 1,decomp%ysz(3)
          do j = 1,decomp%ysz(2)
             do i = 1,decomp%ysz(1)
                  yuniform_adj = eta(i,j,k) - ystart
                  num= (beta+1) - (beta-1)*(BB**(1-yuniform_adj/yh))
                  den = (BB**(1-yuniform_adj/yh)) + 1
                  y(i,j,k) = yh*(num/den) + ystart
             end do
          end do
       end do

    elseif(ymetric_flag==10) then
       ! finite-difference evaluation of metrics (reduces order of accuracy)
       call GracefulExit("flag = 4 (finite-difference evaluation of metrics) is incomplete right now",21)
    endif
    end subroutine

end module


subroutine meshgen(decomp, dx, dy, dz, mesh, inputfile, xmetric, ymetric, zmetric, xi, eta, zeta, dxs, dys, dzs, xbuf, zbuf)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info, nrank, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use Vortexring_data

    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    character(len=*),                intent(in)    :: inputfile
    logical,                         intent(in   ) :: xmetric, ymetric, zmetric
    real(rkind), dimension(:,:,:  ), intent(inout) :: xi, eta, zeta
    real(rkind), dimension(:,:,:  ), intent(inout) :: dxs, dys, dzs
    real(rkind), dimension(:,:,:,:), target,intent(in):: xbuf, zbuf
    real(rkind), dimension(:,:,:), pointer :: xtmp1, xtmp2, ztmp1, ztmp2
    integer :: i,j,k,ioUnit, nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn
    real(rkind) :: xfocus, xtau, xh, xstart
    real(rkind) :: yfocus, ytau, yh, ystart
    real(rkind) :: zfocus, ztau, zh, zstart
    integer     ::  xmetric_flag,  ymetric_flag, zmetric_flag
    real(rkind), allocatable, dimension(:,:) :: metric_params
    character(len=clen) :: outputfile,str

    namelist /PROBINPUT/ ns, Lx, Ly, Lz, Dia, Pr, Sc, gam, Rgas, mu_ref, p_ref, rho_ref, add_noise
    namelist /METRICS/ xmetric_flag, ymetric_flag, zmetric_flag, metric_params

    ioUnit = 15
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)


    allocate(metric_params(3,5))    ! 3 :: (x,y,z); 5 :: max no of parameters
    metric_params = zero
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=METRICS)
    close(ioUnit)

    xtmp1 => xbuf(:,:,:,1); xtmp2 => xbuf(:,:,:,2) 
    ztmp1 => zbuf(:,:,:,1); ztmp2 => zbuf(:,:,:,2) 

    ! Global domain size 
    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)

    ! Need to set x, y and z as well as  dx, dy and dz
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        if (nrank == 0) then
            print *, "Domain size: ",Lx,Ly,Lz
        end if

        dx = Lx/real(nx-1,rkind)  
        dy = Ly/real(ny-1,rkind)
        dz = Lz/real(nz-1,rkind)  

        x1 = 0._rkind
        y1 = -Ly/2._rkind
        z1 = -Lz/2._rkind

        xn = Lx
        yn = Ly/2._rkind
        zn = Lz/2._rkind


        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = x1 + real( ix1 -1 + i - 1, rkind ) * dx
                    y(i,j,k) = y1 + real( iy1 -1 + j - 1, rkind ) * dy
                    z(i,j,k) = z1 + real( iz1 -1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    xfocus = metric_params(1,1);  xtau   = metric_params(1,2);  xstart = metric_params(1,3); xh = metric_params(1,4)
    yfocus = metric_params(2,1);  ytau   = metric_params(2,2);  ystart = metric_params(2,3); yh = metric_params(2,4)
    zfocus = metric_params(3,1);  ztau   = metric_params(3,2);  zstart = metric_params(3,3); zh = metric_params(3,4)

    if(nrank == zero) then
      print*, '>>xfocus=',xfocus, '>>xtau=',xtau, '>>xstart=',xstart, '>>xh=',xh, '>>xmetric=',xmetric, '>>xflag=',xmetric_flag
      print*, '>>yfocus=',yfocus, '>>ytau=',ytau, '>>ystart=',ystart, '>>yh=',yh, '>>ymetric=',ymetric, '>>yflag=',ymetric_flag
      print*, '>>zfocus=',zfocus, '>>ztau=',ztau, '>>zstart=',zstart, '>>zh=',zh, '>>zmetric=',zmetric, '>>zflag=',zmetric_flag
    endif
    
    if(xmetric) then
      xi = x
      call stretched_coordinates(decomp,x,xi,xmetric,xmetric_flag,metric_params(1,1),&
                                 metric_params(1,2),metric_params(1,3),metric_params(1,4))
    endif

    if(ymetric) then
      eta = y
      call stretched_coordinates(decomp,y,eta,ymetric,ymetric_flag,metric_params(2,1),&
                                 metric_params(2,2),metric_params(2,3),metric_params(2,4))
    endif

    if(zmetric) then
      zeta = z
      call stretched_coordinates(decomp,z,zeta,zmetric,zmetric_flag,metric_params(3,1),&
                                 metric_params(3,2),metric_params(3,3),metric_params(3,4))
    endif

    ! Grid width on stretched/uniform mesh
    call transpose_y_to_x(x,xtmp1,decomp)   ! Decomposition in x
    xtmp2(1,:,:) = xtmp1(2,:,:) - xtmp1(1,:,:)
    do i = 2, nx-1
       xtmp2(i,:,:) =  (xtmp1(i+1,:,:) - xtmp1(i-1,:,:))/2
    end do
    xtmp2(nx,:,:) = xtmp1(nx,:,:) - xtmp1(nx-1,:,:)
    call transpose_x_to_y(xtmp2,dxs,decomp)   ! Decomposition in x

    dys(:,1,:) = y(:,2,:) - y(:,1,:)       ! Base decomposition in Y
    do j=2, decomp%ysz(2)-1
       dys(:,j,:) =  (y(:,j+1,:) - y(:,j-1,:))/2
    end do
    dys(:,decomp%ysz(2),:) = y(:,decomp%ysz(2),:) - y(:,decomp%ysz(2)-1,:)       ! Base decomposition in Y

    call transpose_y_to_z(z,ztmp1,decomp)   ! Decomposition in z
    ztmp2(:,:,1) = ztmp1(:,:,2) - ztmp1(:,:,1)
    do k = 2, nz-1
       ztmp2(:,:,k) =  (ztmp1(:,:,k+1) - ztmp1(:,:,k-1))/2
    end do
    ztmp2(:,:,nz) =  ztmp1(:,:,nz) - ztmp1(:,:,nz-1)
    call transpose_z_to_y(ztmp2,dzs,decomp)   ! Decomposition in x

    !! Write grid width to a file
    !write(outputfile, '(a,i0,a)') 'grid_x_', nrank, '.dat'
    !open(11,file=outputfile,status='unknown')
    !do i=1,decomp%ysz(1)
    !   write(11,'(2(e19.12),1x)') x(i,1,1), dxs(i,1,1)
    !enddo
    !close(11)

    !if(nrank==0) then
    !  write(outputfile, '(a)') 'grid_y.dat'
    !  open(10,file=outputfile,status='unknown')
    !  do j=1,decomp%ysz(2)
    !     write(10,'(2(e19.12),1x)') y(1,j,1), dys(1,j,1)
    !  enddo
    !  close(10)
    !endif

    !write(outputfile, '(a,i0,a)') 'grid_z_', nrank, '.dat'
    !open(13,file=outputfile,status='unknown')
    !do k=1,decomp%ysz(3)
    !   write(13,'(2(e19.12),1x)') z(1,1,k), dzs(1,1,k)
    !enddo
    !close(13)

    end associate
    nullify(xtmp1)
    nullify(xtmp2)
    nullify(ztmp1)
    nullify(ztmp2)

end subroutine


subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tsim,tstop,dt,tviz)
    use kind_parameters,             only: rkind, clen
    use constants,                   only: zero,half,one,two,four,five,pi,eight
    use CompressibleGrid,            only: rho_index,u_index,v_index,w_index,&
                                           p_index,T_index,e_index,Ys_index
    use decomp_2d,                   only: decomp_info,nrank
    use MixtureEOSMod,               only: mixture
    use IdealGasEOS,                 only: idealgas
    use SutherLandViscosityMod,      only: sutherlandViscosity
    use ConstRatioBulkViscosityMod,  only: constRatioBulkViscosity
    use SutherLandConductivityMod,   only: SutherLandConductivity
    use ConstSchmidtDiffusivityMod,  only: constSchmidtDiffusivity
    use reductions,                  only: P_MAXVAL,P_MINVAL
    use exits,                       only: GracefulExit, message, nancheck
    use random,                      only: gaussian_random                  

    use Vortexring_data
    use mpi

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind),                     intent(inout) :: tsim, tstop, dt, tviz

    type(sutherlandViscosity) :: shearvisc
    type(constRatioBulkViscosity) :: bulkvisc
    type(sutherlandConductivity) :: thermcond
    real(rkind) :: S, Sk, T0
    integer :: i,j, k, iounit, nx, ny, nz, nxl, nyl, nzl
    character(len=clen) :: outputfile
    real(rkind), dimension(decomp%ysz(1)) :: x_new
    real(rkind), dimension(decomp%ysz(2)) :: y_new
    real(rkind), dimension(decomp%ysz(3)) :: z_new
    
    namelist /PROBINPUT/ ns, Lx, Ly, Lz, Dia, Pr, Sc, gam, Rgas, mu_ref, p_ref, rho_ref, add_noise

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    ! Global domain sizes
    nx = decomp%xsz(1);     ny = decomp%ysz(2);    nz = decomp%zsz(3)

    ! Local domain sizes
    nxl = decomp%ysz(1);    nyl = decomp%ysz(2);   nzl = decomp%ysz(3)
    
    !print *, "checking nx,ny,nz,n ",nx,ny,nz,nxl,nyl,nzl
    associate( rho => fields(:,:,:,rho_index), u  => fields(:,:,:,u_index),&
                 v => fields(:,:,:,  v_index), w  => fields(:,:,:,w_index),&
                 p => fields(:,:,:,  p_index), T  => fields(:,:,:,T_index),&
                 e => fields(:,:,:,  e_index), &
                Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if (mix%ns /= ns) call GracefulExit("Wrong number of species. Check your input file and make ns consistent with the problem file.",4562)

        T0     =  p_ref/(Rgas*rho_ref)        
        S      = 110.4_rkind
        Sk     = 194.0_rkind

        !print*, T0, mu_ref, Pr, S, Sk, Rgas
        !!!! Set each material's transport coefficient object
        shearvisc = sutherlandViscosity(mu_ref, T0, 1.5_rkind, S)
        bulkvisc  = constRatioBulkViscosity( zero )
        thermcond = sutherlandConductivity(Pr, T0, 1.5_rkind, Sk)
        call mix%set_material( 1, idealgas( gam, Rgas ), shearvisc = shearvisc, bulkvisc  = bulkvisc, thermcond = thermcond  )

        Ys(:,:,:,1)  = one 
        call mix%update(Ys)     
     
        ! Add base flow profiles
        u   = u + zero
        v   = v + zero
        w   = w + zero
        p   = p + p_ref
        rho = rho + rho_ref  
        T   = p/(rho*Rgas) 

        ! Initialize gaussian filter mygfil
        call mygfil%init(decomp, periodicx, periodicy, periodicz, "gaussian", "gaussian", "gaussian" )
    end associate
end subroutine


subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use MixtureEOSMod,    only: mixture
    use reductions,       only: P_MEAN
    use Vortexring_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    type(mixture),                   intent(in) :: mix
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    character(len=clen) :: outputfile,str
    integer :: i,outputunit=229

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        !write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/Vortexring_t.dat"

        !open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        !do i=1,decomp%ysz(1)
        !    write(outputunit,'(8ES26.16)') x(i,1,1), rho(i,1,1), u(i,1,1), e(i,1,1), p(i,1,1), &
        !                                   mu(i,1,1), bulk(i,1,1), kap(i,1,1)
        !
        !end do
        !close(outputunit)
    end associate
end subroutine


subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc,newTimeStep, time_step)
    use kind_parameters,  only: rkind
    use decomp_2d,        only: decomp_info, nrank, transpose_y_to_x, transpose_x_to_y
    use constants,        only: zero, half, one, two, three, four, five, six, seven, eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use MixtureEOSMod,    only: mixture
    use operators,        only: filter3D

    use Vortexring_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc
    logical,                         intent(in)    :: newTimeStep
    integer,                         intent(in)    :: time_step 

    integer :: i, j, k, nx, ny, nz, ix1_new, iy1_new, iz1_new, tidx
    real(rkind) :: dx, dy, dz,rad, filpt, thickT, U0, P0, rho0, T0
    real(rkind) :: umin, pmin, Tmin, rhomin, diff_u, diff_rho, diff_T, diff_p
    character(len=clen) :: outputfile
    real(rkind), dimension(:,:),       allocatable :: u_noise, v_noise, w_noise
    real(rkind), dimension(:,:,:),     allocatable :: u_xtmp, v_xtmp, w_xtmp

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

    !!!! NOTE: Make sure Vortexring_data, Lx, Ly, Lz, Rgas, filpt, thickT are mentioned correctly !!!!    
        call get_fitted_fields(tsim, P0, rho0, U0)

        umin = zero;             pmin = p_ref;              rhomin = rho_ref
        diff_u = half*(U0-umin); diff_p = half*(P0-pmin);   diff_rho = half*(rho0-rhomin)
          
        if(nrank == zero) then
        print*, ">>> tsim=",tsim, ">>> U0=",U0, ">>> P0=",P0, ">>> rho0=",rho0
        endif

        ! set Dirichlet BC for velocity, density, pressure
        if(decomp%yst(1) == 1) then 
            do k = 1, decomp%ysz(3) 
               do j = 1, decomp%ysz(2)
                   rad = sqrt(y(1,j,1)**2 + z(1,1,k)**2)
                   u(1,j,k)   =  umin + (one - tanh(1100.0_rkind*(rad-half*Dia))) * diff_u
                   p(1,j,k)   =  pmin + (one - tanh(1100.0_rkind*(rad-half*Dia))) * diff_p
                   rho(1,j,k) =  rhomin + (one - tanh(1100.0_rkind*(rad-half*Dia))) * diff_rho
                   T(1,j,k)   =  p(1,j,k) / (Rgas*rho(1,j,k))
               enddo
            enddo
            v(1,:,:)   = zero
            w(1,:,:)   = zero
         endif

        !!!!!!==========Add noise here=============!!!!!!!!!!!!!
        ! Global size
        nx  = decomp%xsz(1);    ny  = decomp%ysz(2) ;  nz  = decomp%zsz(3)
        ! Transpose u,v,w to x decomposition and add noise
        if (newTimeStep) then
          if (add_noise) then
            allocate(u_xtmp(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
            allocate(v_xtmp(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
            allocate(w_xtmp(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
            allocate(u_noise(ny, nz));      allocate(v_noise(ny, nz));         allocate(w_noise(ny, nz))
            call transpose_y_to_x(u,u_xtmp,decomp)  
            call transpose_y_to_x(v,v_xtmp,decomp)  
            call transpose_y_to_x(w,w_xtmp,decomp)  

            call find_noise(decomp, Lx, Ly, Lz, Dia, x, y, z, u, v, w, u_noise, v_noise, w_noise, time_step, newTimeStep)

            tidx = mod(time_step, nx) + 1 

            ! If base decomposition is in X
            ix1_new = decomp%xst(1); iy1_new = decomp%xst(2); iz1_new = decomp%xst(3) 
            ! print*, "I am rank=", nrank, "ix1=", ix1_new, "iy1", iy1_new ,"iz1", iz1_new
            do k = 1, decomp%xsz(3)       
               do j = 1, decomp%xsz(2)       
                  u_xtmp(1, j, k)  = u_xtmp(1, j, k) + u_noise(iy1_new + j - 1, iz1_new + k - 1)
                  v_xtmp(1, j, k)  = v_xtmp(1, j, k) + v_noise(iy1_new + j - 1, iz1_new + k - 1)
                  w_xtmp(1, j, k)  = w_xtmp(1, j, k) + w_noise(iy1_new + j - 1, iz1_new + k - 1)
               end do
            end do
            call transpose_x_to_y(u_xtmp,u,decomp)  
            call transpose_x_to_y(v_xtmp,v,decomp)  
            call transpose_x_to_y(w_xtmp,w,decomp)  

            deallocate(u_noise, v_noise, w_noise)
            deallocate(u_xtmp, v_xtmp, w_xtmp)
          end if
        end if

        !!!!! =============  Add Sponge+bulk for exit bc ==========!!!!!
        ! Gradually apply the exit boundary conditions
        ! Apply sponge in X-direction on right
        call  sponge_x(decomp, mygfil, x, Lx, u, v, w, p, rho, x_bc, y_bc, z_bc)

        ! Apply sponge in Y-direction on top and bottom
        call  sponge_y(decomp, mygfil, y, Ly, u, v, w, p, rho, x_bc, y_bc, z_bc)

        ! Apply sponge in Z-direction on front and back
        call  sponge_z(decomp, mygfil, z, Lz, u, v, w, p, rho, x_bc, y_bc, z_bc)
    end associate
end subroutine


subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim,sgsmodel)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,two
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info, nrank
    use MixtureEOSMod,    only: mixture
    use sgsmod_cgrid,     only: sgs_cgrid
    use exits,            only: message
    use reductions,       only: P_MAXVAL,P_MINVAL

    use Vortexring_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(mixture),                   intent(in) :: mix
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(sgs_cgrid), optional,       intent(in) :: sgsmodel

    real(rkind) :: dx, Ythick, oob
    integer :: ny  , j
    integer :: iounit = 229
    character(len=clen) :: outputfile
    !real(rkind), dimension(decomp%ysz(2)) :: cmodel_loc, cmodel_loc_Qjsgs, cmodel_loc_tke

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        call message(2,"Maximum u-velocity",P_MAXVAL(u))
        call message(2,"Maximum v-velocity",P_MAXVAL(v))
        call message(2,"Maximum pressure",P_MAXVAL(p))
        call message(2,"Maximum density",P_MAXVAL(rho))
        call message(2,"Maximum temperature",P_MAXVAL(T))
        call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))
        call message(2,"Maximum diffusivity",P_MAXVAL(diff))


    end associate

end subroutine

subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs)
    use kind_parameters, only: rkind
    use decomp_2d,       only: decomp_info
    use MixtureEOSMod,   only: mixture

    use Vortexring_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs

end subroutine
