module ShearLayer_data
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero, half, one, two, four, pi, imi
    use FiltersMod,       only: filters
    use decomp_2d,        only: decomp_info, nrank, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use basic_io,         only: read_2d_ascii 
    use reductions,       only: P_MAXVAL,P_MINVAL
    use exits,            only: message, message_min_max
    use mpi
    implicit none
    integer, parameter :: ns = 1

    ! Problem parameters
    real(rkind) :: Mc = 0.6_rkind           ! Convective Mach
    real(rkind) :: Re = 50.0_rkind          ! Reynolds number
    real(rkind) :: Sc = 1._rkind            ! Schmidt number
    real(rkind) :: p_ref = one              ! reference press
    real(rkind) :: T_ref = one              ! reference temp
    real(rkind) :: rho_ref = one            ! reference density
    real(rkind) :: rho_ratio = one          ! rho2/rho1
    real(rkind) :: dtheta0 = 1._rkind       ! Base profile thickness 
    real(rkind) :: noiseAmp = 1D-6          ! white noise amplitude
    character(len=clen) :: fname_prefix
    logical :: use_lstab = .false. 
    logical :: no_pert = .false. 
    
    ! Parameters for the 2 materials
    real(rkind):: gam=1.4_rkind
    real(rkind):: Pr=0.7_rkind 
    real(rkind):: Rgas=one

    ! Domain size data
    real(rkind) :: Ly = 59._rkind, Lx=68._rkind, Lz=68._rkind
    real(rkind) :: x1, y1, z1
    logical :: periodicx = .true., periodicy = .false., periodicz = .true.

    ! Gaussian filter for sponge
    type(filters) :: mygfil
contains
    subroutine stretched_coordinates(decomp, y, eta, ymetric, ymetric_flag, yfocus, ytau, ystart, yh)
    use constants,        only: zero, half, one
    use decomp_2d,        only: decomp_info, nrank
    use exits,            only: GracefulExit, message, nancheck

    type(decomp_info),               intent(in)    :: decomp
    real(rkind), dimension(:,:,:),   intent(inout) :: y
    real(rkind), dimension(:,:,:  ), intent(inout) :: eta
    logical,                         intent(in   ) :: ymetric
    integer,                         intent(in)    :: ymetric_flag
    real(rkind),                     intent(in)    :: yfocus, ytau, ystart, yh
    integer     :: i,j,k
    real(rkind) :: yfocus_adj, num, den, BB, yuniform_adj

    !print*, yfocus, ytau, ystart, yh, ymetric, ymetric_flag
    !print*, decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)

    if(ymetric_flag==1) then
       ! concentrate towards the center -- Pletcher, Tannehill, Anderson
       ! (Section 5.6, Transformation 3, pg. 332) 
       yfocus_adj = yfocus - ystart
       num = one + (yfocus_adj/yh) * (exp( ytau) - one)
       den = one + (yfocus_adj/yh) * (exp(-ytau) - one)
       BB = half/ytau*log(num/den)
       do k=1,decomp%ysz(3)
          do j=1,decomp%ysz(2)
             do i=1,decomp%ysz(1)
                ! adjust for starting point
                yuniform_adj = (eta(i,j,k) - ystart) !/ yh
                ! stretched location
                num = sinh(ytau*BB)
                y(i,j,k) = yfocus_adj * (one + sinh(ytau * (yuniform_adj/yh-BB))/num) + ystart
             end do
          end do
       end do
       elseif(ymetric_flag==2) then
          ! concentrate towards the two ends
          call GracefulExit("flag = 2 (concentrate towards two ends) is incomplete right now",21)
       elseif(ymetric_flag==3) then
          ! concentrate at arbitrary point
          call GracefulExit("flag = 3 (concentrate at arbitrary point) is incomplete right now",21)
       elseif(ymetric_flag==10) then
          ! finite-difference evaluation of metrics (reduces order of accuracy)
          call GracefulExit("flag = 4 (finite-difference evaluation of metrics) is incomplete right now",21)
       endif
    end subroutine

    subroutine perturb_potential_v2(gp,x,y,z,nx,ny,nz,nxl,nyl,nzl,Lx,Lz,u,v,w,p,rho,fname_prefix)
        use decomp_2d,        only: nrank
        use decomp_2d_io
        type(decomp_info), intent(in)               :: gp
        real(rkind), dimension(:,:,:), intent(in)   :: x,y,z
        real(rkind), dimension(:,:,:), intent(inout):: u,v,w,p,rho
        real(rkind), intent(in)                     :: Lx,Lz
        integer, intent(in)                         :: nx,ny,nz,nxl,nyl,nzl
        character(len=clen), intent(in) :: fname_prefix

        character(len=clen) :: tempname, fname
        real(rkind), dimension(:,:,:), allocatable :: uperturb, vperturb, wperturb,pperturb,rperturb
       
       ! real, dimension(:), allocatable :: A,e,x1d,y1d,z1d
       ! real :: kx, kz, phx, phz, kmin, Amax, du
       ! real :: sigma=10,pi = 3.14159265358
       ! integer :: ni,nj,nk,i, j, k, m, mx, mz, mode_min=4, nmodes_max=4, &
       !     nxmodes, nzmodes, mpi_ierr, ix1,iz1,ixn,izn
   
        ! Global size
       ! ni  = gp%xsz(1)
       ! nj  = gp%ysz(2)
       ! nk  = gp%zsz(3)

       ! ! If base decomposition is in Y
        !ix1 = gp%yst(1); 
        !iz1 = gp%yst(3)
        !ixn = gp%yen(1); 
        !izn = gp%yen(3)
       ! print '(a,5(i4,1x))', 'nrank = ', nrank, ix1, ixn, iz1, izn

        ! Store radius, theta, z 
       ! allocate(x1D(gp%ysz(1)),y1D(gp%ysz(2)),z1D(gp%ysz(3)))
       ! x1d = x(:,1,1)
       ! y1D = y(1,:,1) 
       ! z1d = z(1,1,:)

       ! ! Transverse mask 
        !allocate(e(nj),A(nj))
        !e = exp(-sigma*(y1D**2));
        !du = P_MAXVAL(u)-P_MINVAL(u)
        !Amax = 0.01*(du)
        !call mpi_bcast(Amax,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)
        !A = Amax*e

        call message(0,"Before adding perturbations")
        call message(2,"Maximum u", P_MAXVAL(abs(u)))
        call message(2,"Maximum v", P_MAXVAL(abs(v)))
        call message(2,"Maximum w", P_MAXVAL(abs(w)))
        call message(2,"Minimum u", P_MINVAL(abs(u)))
        call message(2,"Minimum v", P_MINVAL(abs(v)))
        call message(2,"Minimum w", P_MINVAL(abs(w)))
    
        !nxmodes = min(nmodes_max, ni/4)
        !nzmodes = min(nmodes_max, nk/4)
        
       ! do mx=mode_min,mode_min+nxmodes
       ! do mz=mode_min,mode_min+nzmodes
       !     kx = 2.D0*pi/Lx * mx
       !     kz = 2.D0*pi/Lz * mz
       !     if(nrank==0) then
       !         call random_number(phx)
       !         call random_number(phz)
       !     endif
       !     call mpi_bcast(phx,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)
       !     call mpi_bcast(phz,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)
       !     phx = phx*2.D0*pi
       !     phz = phz*2.D0*pi
       !     do i=1,gp%ysz(1) !ix1,ixn
       !     do k=1,gp%ysz(3) !iz1,izn
       !         !write(100+nrank,'(4(i4,1x),6(e19.12,1x))') i, k, size(x1D), size(z1D), phx, phz, x1D(i), kx*x1D(i), z1D(k), kz*z1D(k)
       !         !u(i,:,k) = u(i,:,k) + A*sin(kx*x1D(i)+phx)*sin(kz*z1D(k)+phz)
       !         !w(i,:,k) = w(i,:,k) + A*cos(kx*x1D(i)+phx)*cos(kz*z1D(k)+phz)
       !         !v(i,:,k) = v(i,:,k) + A/kz * sin(kx*x1D(i)+phx) * &
       !         !    ( (-2d0*sigma*y1d)*cos(kz*z1D(k)+phz) + kx*sin(kz*z1D(k)+phz) )
       !         write(100+nrank,*) i, k
       !     enddo
       !     enddo
       ! enddo
       ! enddo
       ! print '(a,6(i4,1x))', 'nrank = ', nrank, size(A), size(e), size(x1d), size(y1d), size(z1d)
       ! !deallocate(A,e,x1d,y1d,z1d)

    
        allocate(uperturb(nxl, nyl, nzl))
        allocate(vperturb(nxl, nyl, nzl))
        allocate(wperturb(nxl, nyl, nzl))
        allocate(pperturb(nxl, nyl, nzl))
        allocate(rperturb(nxl, nyl, nzl))
        
        write(tempname,"(A7,A3,I4.4,A,I4.4,A,I4.4,A4)") "perturb", "_u_",nx, "_",ny,"_",nz,".dat"
        fname = fname_prefix(:len_trim(fname_prefix))//"/"//trim(tempname)
        call decomp_2d_read_one(2,uperturb,fname, gp)
        
        write(tempname,"(A7,A3,I4.4,A,I4.4,A,I4.4,A4)") "perturb", "_v_",nx, "_",ny,"_",nz,".dat"
        fname = fname_prefix(:len_trim(fname_prefix))//"/"//trim(tempname)
        call decomp_2d_read_one(2,vperturb,fname, gp)
      
        write(tempname,"(A7,A3,I4.4,A,I4.4,A,I4.4,A4)") "perturb", "_w_",nx, "_",ny,"_",nz,".dat"
        fname = fname_prefix(:len_trim(fname_prefix))//"/"//trim(tempname)
        call decomp_2d_read_one(2,wperturb,fname, gp)

        write(tempname,"(A7,A3,I4.4,A,I4.4,A,I4.4,A4)") "perturb", "_p_",nx, "_",ny,"_",nz,".dat"
        fname = fname_prefix(:len_trim(fname_prefix))//"/"//trim(tempname)
        call decomp_2d_read_one(2,pperturb,fname, gp)

        write(tempname,"(A7,A3,I4.4,A,I4.4,A,I4.4,A4)") "perturb", "_r_",nx, "_",ny,"_",nz,".dat"
        fname = fname_prefix(:len_trim(fname_prefix))//"/"//trim(tempname)
        call decomp_2d_read_one(2,rperturb,fname, gp)

        u = u + uperturb * exp(-(y**2.0_rkind)/8.0_rkind)
        v = v + vperturb * exp(-(y**2.0_rkind)/8.0_rkind)
        w = w + wperturb * exp(-(y**2.0_rkind)/8.0_rkind)
        p = p + pperturb * exp(-(y**2.0_rkind)/8.0_rkind)
        rho = rho + rperturb* exp(-(y**2.0_rkind)/8.0_rkind)

        deallocate(uperturb)
        deallocate(vperturb)
        deallocate(wperturb)
        deallocate(pperturb)
        deallocate(rperturb)

        call message(0,"After adding perturbations")
        call message(2,"Maximum u", P_MAXVAL(abs(u)))
        call message(2,"Maximum v", P_MAXVAL(abs(v)))
        call message(2,"Maximum w", P_MAXVAL(abs(w)))
        call message(2,"Minimum u", P_MINVAL(abs(u)))
        call message(2,"Minimum v", P_MINVAL(abs(v)))
        call message(2,"Minimum w", P_MINVAL(abs(w)))
    end subroutine 

    ! KVM 2021
    subroutine my_ddy(decomp,dy,f,dfdy)
        type(decomp_info), intent(in)               :: decomp
        real(rkind),intent(in)                      :: dy
        real(rkind), dimension(:,:,:), intent(in)   :: f
        real(rkind), dimension(:,:,:), intent(inout):: dfdy

        integer :: j,ni,nj,nk, iy1,iyn,mpi_ierr 

        nj = decomp%ysz(2)
        ! If base decomposition is in Y
        iy1 = decomp%yst(2)
        iyn = decomp%yen(2)
         
        ! df/dy: O(h^4) centered diff
        !http://www.dam.brown.edu/people/alcyew/handouts/numdiff.pdf 
        dfdy = 0.d0
        do j=3+iy1,iyn-3
            dfdy(:,j,:) = ( -f(:,j+2,:) + f(:,j-2,:) &
                       + 8.d0*f(:,j+1,:) - 8.d0*f(:,j-1,:)) / (12.d0*dy)
        enddo
            
    end subroutine 

end module


subroutine meshgen(decomp, dx, dy, dz, mesh, inputfile, xmetric, ymetric, zmetric, xi, eta, zeta, dxs, dys, dzs, xbuf, zbuf)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use decomp_2d,        only: decomp_info, nrank, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use exits,            only: GracefulExit, message, nancheck
    use ShearLayer_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    character(len=*),                intent(in)    :: inputfile
    logical,                         intent(in   ) :: xmetric, ymetric, zmetric
    real(rkind), dimension(:,:,:  ), intent(inout) :: xi, eta, zeta
    real(rkind), dimension(:      ), intent(inout) :: dxs, dys, dzs
    real(rkind), dimension(:,:,:,:), target,intent(in):: xbuf, zbuf
    integer :: i,j,k,ioUnit, nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn
    real(rkind) :: xfocus, xtau, xh, xstart
    real(rkind) :: yfocus, ytau, yh, ystart
    real(rkind) :: zfocus, ztau, zh, zstart
    integer     ::  xmetric_flag,  ymetric_flag, zmetric_flag
    real(rkind), allocatable, dimension(:,:) :: metric_params
    character(len=clen) :: outputfile,str
    real(rkind), dimension(:,:,:), pointer :: xtmp1, xtmp2, ztmp1, ztmp2

    namelist /PROBINPUT/ Lx, Ly, Lz,Mc, Re, Pr, Sc,&
                        T_ref, p_ref, rho_ref, rho_ratio,&
                        noiseAmp, fname_prefix, use_lstab, no_pert

    namelist /METRICS/ xmetric_flag, ymetric_flag, zmetric_flag, metric_params

    ioUnit = 11
    open(unit=ioUnit, file=inputfile, form='FORMATTED')
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
            print *, "Domain size: ",Lx,Ly,Lz ! ,ix1,iy1,iz1,ixn,iyn,izn
        end if
        dx = Lx/real(nx,rkind)
        dy = Ly/real(ny,rkind)
        dz = Lz/real(nz,rkind)
        x1 = 0._rkind! -Lx/2._rkind
        y1 = -Ly/2._rkind
        z1 = 0._rkind!-Lz/2._rkind
        
        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = x1 + real( ix1 + i - 1, rkind ) * dx
                    y(i,j,k) = y1 + real( iy1 + j - 1, rkind ) * dy
                    z(i,j,k) = z1 + real( iz1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    xfocus = metric_params(1,1);  xtau   = metric_params(1,2);  xstart = metric_params(1,3); xh = metric_params(1,4)
    yfocus = metric_params(2,1);  ytau   = metric_params(2,2);  ystart = metric_params(2,3); yh = metric_params(2,4)
    zfocus = metric_params(3,1);  ztau   = metric_params(3,2);  zstart = metric_params(3,3); zh = metric_params(3,4)

    !print*, xfocus, xtau, xstart, xh, xmetric, xmetric_flag
    !print*, yfocus, ytau, ystart, yh, ymetric, ymetric_flag
    !print*, zfocus, ztau, zstart, zh, zmetric, zmetric_flag
    
    if(xmetric) then
      xi = x
      call stretched_coordinates(decomp, x, xi, xmetric, xmetric_flag, xfocus, xtau, xstart, xh)
    endif

    if(ymetric) then
      eta = y
      call stretched_coordinates(decomp, y, eta, ymetric, ymetric_flag, yfocus, ytau, ystart, yh)
    endif

    if(zmetric) then
      zeta = z
      call stretched_coordinates(decomp, z, zeta, zmetric, zmetric_flag, zfocus, ztau, zstart, zh)
    endif


    !print*, "nrank=",nrank, "ix1=",ix1, "ixn=",ixn
    !print*, "nrank=",nrank, "iy1=",iy1, "iyn=",iyn
    !print*, "nrank=",nrank, "iz1=",iz1, "izn=",izn

    ! Grid width on stretched/uniform mesh
    call transpose_y_to_x(x,xtmp1,decomp)   ! Decomposition in x
    xtmp2(1,:,:) = xtmp1(1,:,:) - x1
    do i = 2, nx
       xtmp2(i,:,:) =  xtmp1(i,:,:) - xtmp1(i-1,:,:)
    end do
   
    do i = 1, decomp%ysz(1)
       dxs(i) = xtmp2(ix1 + i -1, 1, 1)
    end do

    dys(1) = y(1,1,1) - y1  ! Base decomposition in Y
    do j=2, decomp%ysz(2)
       dys(j) =  y(1,j,1) - y(1,j-1,1)
    end do

    call transpose_y_to_z(z,ztmp1,decomp)   ! Decomposition in z
    ztmp2(:,:,1) = ztmp1(:,:,1) - z1
    do k = 2, nz
       ztmp2(:,:,k) =  ztmp1(:,:,k) - ztmp1(:,:,k-1)
    end do
   
    do k = 1, decomp%ysz(3)
       dzs(k) = ztmp2(1, 1, iz1 + k -1)
    end do

    !!! Write grid width to a file
    write(outputfile, '(a,i0,a)') 'grid_x_', nrank, '.dat'
    open(11,file=outputfile,status='unknown')
    do i=1,decomp%ysz(1)
       write(11,'(2(e19.12),1x)') x(i,1,1), dxs(i)
    enddo
    close(11)

    if(nrank==0) then
      write(outputfile, '(a)') 'grid_y.dat'
      open(10,file=outputfile,status='unknown')
      do j=1,decomp%ysz(2)
         write(10,'(2(e19.12),1x)') y(1,j,1), dys(j)
      enddo
      close(10)
    endif

    write(outputfile, '(a,i0,a)') 'grid_z_', nrank, '.dat'
    open(13,file=outputfile,status='unknown')
    do k=1,decomp%ysz(3)
       write(13,'(2(e19.12),1x)') z(1,1,k), dzs(k)
    enddo
    close(13)
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
    use SutherLandConductivityMod,   only:SutherLandConductivity
    use ConstSchmidtDiffusivityMod,  only: constSchmidtDiffusivity
    use reductions,                  only: P_MAXVAL,P_MINVAL
    use exits,                       only: GracefulExit, message, nancheck
    use random,                      only: gaussian_random                  

    use ShearLayer_data
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
    !real(rkind), dimension(:,:), allocatable :: tmp2D
    real(rkind) :: mu_ref, c1,c2,du, Rgas1, Rgas2,lambda,kx, kz, ph, S, Sk, T0
    integer :: i,iounit, nx, ny, nz, nxl, nyl, nzl
    
    namelist /PROBINPUT/ Lx, Ly, Lz,Mc, Re, Pr, Sc,&
                        T_ref, p_ref, rho_ratio, rho_ref, &
                        noiseAmp, fname_prefix, use_lstab, no_pert 
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
        
        ! First set the materials for rho_ratio = rho2/rho1
        Rgas1 = (one+rho_ratio)/two
        Rgas2 = (one+rho_ratio)/(two*rho_ratio)
        Rgas = Rgas1
        c1 = sqrt(gam*p_ref/(rho_ref/Rgas1))
        c2 = sqrt(gam*p_ref/(rho_ref/Rgas2))
        du = Mc*(c1+c2)
        mu_ref = one/Re
  
        T0 =  101325.0_rkind/(287.0_rkind*1.2_rkind)
        
        S = 110.4_rkind/T0
        Sk = 194.0_rkind/T0

        ! Set each material's transport coefficient object
        shearvisc = sutherlandViscosity(mu_ref, T_ref, 1.5_rkind, S)
        bulkvisc  = constRatioBulkViscosity( zero )
        thermcond = sutherlandConductivity(Pr, T_ref, 1.5_rkind, Sk)
        call mix%set_material( 1, idealgas( gam, Rgas ),&
             shearvisc = shearvisc, & 
             bulkvisc  = bulkvisc, &
             thermcond = thermcond  )

        ! Set mass diffusivity object (Ensure that all units are consistent)
        lambda = (rho_ratio-1)/(rho_ratio+1)
        !call mix%set_massdiffusivity( constSchmidtDiffusivity( mu_ref,rho_ref,Sc))
        Ys(:,:,:,1)  = one !- half*(one+lambda*tanh(y/(two*dtheta0)))
        !Ys(:,:,:,2)  = one - Ys(:,:,:,1)
        call mix%update(Ys)

        ! Add base flow profiles
        lambda = (rho_ratio-1)/(rho_ratio+1)
        u = u + du*half*tanh(y/(two*dtheta0))
        v = v + zero
        w = w + zero
        p = p + p_ref
        rho = rho + rho_ref*(1 + lambda*tanh(y/(two*dtheta0)))
        T = p/(rho*Rgas) 

        ! Perturbations: this must be specific for each problem.
        if (use_lstab) then
            call GracefulExit("LSTAB init is deprecated.",4562)
            !call lstab_pert(decomp,x,z,fname_prefix,u,2)
            !call lstab_pert(decomp,x,z,fname_prefix,v,3)
            !call lstab_pert(decomp,x,z,fname_prefix,w,4)
            !call lstab_pert(decomp,x,z,fname_prefix,rho,1)
        else if (no_pert) then
            call message(0,"No perturbations")
        else
            call perturb_potential_v2(decomp,x,y,z,nx,ny,nz,nxl,nyl,nzl,Lx,Lz,u,v,w,p,rho,fname_prefix)
        endif

        ! Initialize gaussian filter mygfil
        call mygfil%init(decomp, periodicx, periodicy, periodicz, &
                        "gaussian", "gaussian", "gaussian" )
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
    use ShearLayer_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    type(mixture),                   intent(in) :: mix
    real(rkind),                     intent(in) :: dx,dy,dz,tsim
    integer,                         intent(in) :: vizcount
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) ::&
    tke,ubase
    real(rkind) :: tke_mean,tke0
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

        write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/ShearLayer1Mat_t.dat"!, &
            !vizcount, ".dat"

        ! Get TKE
        ubase = 0.5*(1+tanh(y))
        tke = half*rho*((u-ubase)**2 + v*v + w*w)

        !write(str,'(I4.4)') decomp%ysz(2)
        !write(outputfile,'(2A)') trim(outputdir),"/ShearLayer_"//trim(str)//".dat"

        if (vizcount == 0) then
            ! On the first step, get initial disturbance energy and
            ! write the header for the output file.
            tke0 = P_MEAN( tke )
            if (nrank == 0) then
                open(unit=outputunit, file=trim(outputfile), &
                    form='FORMATTED', status='REPLACE')
                write(outputunit,'(3A26)') "Time", "TKE"
            end if
        else
            ! Open the previously started file
            if (nrank == 0) then
                open(unit=outputunit, file=trim(outputfile), &
                    form='FORMATTED', position='APPEND', status='OLD')
            end if
        end if

        tke_mean = P_MEAN(tke)
        if (nrank == 0) then
            write(outputunit,'(3ES26.16)') tsim, tke_mean
            close(outputunit)
        end if

        !!!! write out cmodel_local and cmodel_local_Qjsgs if using
        !!!! DynamicProcedureType==1 so that we can analyze if the
        !!!! LD implementation is correct

    end associate
end subroutine


subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
    use kind_parameters,  only: rkind
    use constants,        only: zero, half, one
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use operators,        only: filter3D

    use ShearLayer_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc

    integer :: i
    real(rkind) :: dy, filpt, thickT
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: dumT, dumF

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        ! Sponge+bulk for exit bc
        ! Gradually apply the exit boundary conditions
        dy = Ly/real(decomp%ysz(2)-1,rkind)
        filpt = 3.0_rkind/dy 
        thickT = real(5.0D0, rkind)
        
        ! Gussian Filter for 
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

        ! Gussian Filter for the top
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

    use ShearLayer_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(mixture),                   intent(in) :: mix
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(sgs_cgrid),optional,        intent(in) :: sgsmodel

    real(rkind) :: dx, Ythick, oob
    integer :: ny  , j
    integer :: iounit = 229
    character(len=clen) :: outputfile
    real(rkind), dimension(decomp%ysz(2)) :: cmodel_loc, cmodel_loc_Qjsgs, cmodel_loc_tke

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        call message(2,"Maximum v-velocity",P_MAXVAL(v))
        call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))
        call message(2,"Maximum diffusivity",P_MAXVAL(diff))

        !if(useSGS)
          !if(sgsmodel%DynamicProcedureType==1) then
             call message_min_max(2,"Bounds for LD-Coeff-tke  : ",     &
                    sgsmodel%get_Max_LocalDynamicProcedure_Coeff_tke(),   &
                    sgsmodel%get_Min_LocalDynamicProcedure_Coeff_tke())
             call message_min_max(2,"Bounds for LD-Coeff      : ",     &
                    sgsmodel%get_Max_LocalDynamicProcedure_Coeff(),   &
                    sgsmodel%get_Min_LocalDynamicProcedure_Coeff())
             call message_min_max(2,"Bounds for LD-Coeff-Qjsgs: ",     &
                    sgsmodel%get_Max_LocalDynamicProcedure_Coeff_Qjsgs(),   &
                    sgsmodel%get_Min_LocalDynamicProcedure_Coeff_Qjsgs())
          !endif

           if(mod(step,100)==0) then
              cmodel_loc = sgsmodel%get_LocalDynamicProcedure_Coeff()
              cmodel_loc_Qjsgs = sgsmodel%get_LocalDynamicProcedure_Coeff_Qjsgs()
              cmodel_loc_tke   = sgsmodel%get_LocalDynamicProcedure_Coeff_tke()
              if(nrank==0) then
                  write(outputfile, '(a,i5.5,a)') 'cmodel_', step, '.dat'
                  open(10,file=outputfile,status='unknown')
                  do j=1,decomp%ysz(2)
                     write(10,'(4(e19.12),1x)') y(1,j,1), cmodel_loc(j), cmodel_loc_Qjsgs(j), cmodel_loc_tke(j)
                     !write(10,'(3(e19.12),1x)') y(1,j,1), cmodel_loc(j), cmodel_loc_Qjsgs(j)
                  enddo
                  close(10)
              endif
           endif

       !! endif

    end associate

end subroutine


subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs,Wcnsrv,tkeb,tsim_0,dtheta_0)
    use CompressibleGrid,   only: rho_index,u_index,v_index,w_index,&
                                  p_index,T_index,e_index,Ys_index
    use kind_parameters,    only: rkind
    use decomp_2d,          only: decomp_info,nrank
    use MixtureEOSMod,      only: mixture
    use reductions,         only: P_MAXVAL,P_MINVAL
    use TKEBudgetMod,       only: tkeBudget
    use ShearLayer_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs

    real(rkind), dimension(:,:,:,:), optional,intent(in)    :: Wcnsrv
    type(tkeBudget), optional,       intent(inout) :: tkeb
    real(rkind), optional,           intent(in)    :: tsim_0 ! the previous time 
    real(rkind), optional,           intent(in)    :: dtheta_0 ! the previous L99 
        
    real(rkind), dimension(:,:,:), allocatable :: &
        rbar,qtilde,vtilde,Rij,src,rhs_mn
    real(rkind), dimension(:,:,:), allocatable :: qpp,vpp
    real(rkind) :: dy
    integer :: i,j,nx,ny,nz,iy1,iyn,mpi_ierr
    integer :: mass_index, mom_index, TE_index
    logical :: forcing
    
    if (present(tkeb)) then
        forcing=.true.
    else 
        forcing = .false.
    endif

    if (forcing) then   
    
    ! If base decomposition is in Y
    iy1 = decomp%yst(2)
    iyn = decomp%yen(2)

    nx = tkeb%avg%avg_size(1)
    ny = tkeb%avg%avg_size(2)
    nz = tkeb%avg%avg_size(3)
    allocate(rbar   (nx,ny,nz))
    allocate(qtilde (nx,ny,nz))
    allocate(vtilde (nx,ny,nz))
    allocate(Rij    (nx,ny,nz))
    allocate(src    (nx,ny,nz))
    allocate(rhs_mn (nx,ny,nz))
    nx = tkeb%avg%sz(1)
    ny = tkeb%avg%sz(2)
    nz = tkeb%avg%sz(3)
    allocate(qpp(nx,ny,nz))
    allocate(vpp(nx,ny,nz))
    
    associate( rho => fields(:,:,:,rho_index), u  => fields(:,:,:,u_index),&
                 v => fields(:,:,:,  v_index), w  => fields(:,:,:,w_index),&
                 p => fields(:,:,:,  p_index), T  => fields(:,:,:,T_index),&
                 e => fields(:,:,:,  e_index), &
                Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1), &
                 y => mesh(:,:,:,2) )
     
    !mpi_bcast(factor,1,mpirkind,0,MPI_COMM_WORLD,mpi_ierr)
    !reynolds_avg_and_fluct(this, f, f_bar, f_prime)
    !favre_avg_and_fluct(this, rho, f, f_tilde, f_pprime)
    !my_ddy(decomp,y,f,dfdy)
       
    ! 1: Y1
    ! 2: Y2
    ! 3: rhou
    ! 4: rhov
    ! 5: rhow
    ! 6: TE
    mass_index = ns
    mom_index = ns+1
    TE_index = mom_index+3

    dy = abs(y(1,iy1,1)-y(1,iy1+1,1))
    call tkeb%reynolds_avg(rho, rbar)
    call tkeb%favre_avg_and_fluct(rho, v, vtilde, vpp)

    ! Mass (fraction) source:
    if (ns.gt.1) then
        do i = 1,ns
            call tkeb%favre_avg_and_fluct(rho, Ys(:,:,:,i), qtilde, qpp) 
            call tkeb%favre_avg(rho,qpp*vpp, Rij)
            call my_ddy(decomp,dy,rbar*(qtilde*vtilde+Rij),src)
            do j=iy1,iyn 
                rhs(:,j,:,i) = rhs(:,j,:,i) + src(1,j,1)
            enddo
        enddo 
    else
        call my_ddy(decomp,dy,rbar*vtilde,src)
        do j=iy1,iyn 
            rhs(:,j,:,mass_index) = rhs(:,j,:,mass_index) + src(1,j,1)
        enddo
    endif
    !call message(4,"Maximum src rho",P_MAXVAL(src))
   
    !call tkeb%reynolds_avg(rhs(:,:,:,2),rhs_mn) 
    !if (nrank==0) then
    !    open(1,file='./src_Y2.txt')
    !    write (1,*) qtilde, src,rhs_mn
    !    close(1)
    !endif
    
    ! X momentum source:
    call tkeb%favre_avg_and_fluct(rho, u, qtilde, qpp)
    call tkeb%favre_avg(rho,qpp*vpp, Rij)
    call my_ddy(decomp,dy,rbar*(qtilde*vtilde+Rij),src)
    do j=iy1,iyn 
        rhs(:,j,:,mom_index) = rhs(:,j,:,mom_index) + src(1,j,1)
    enddo
    !call message(4,"Maximum src momx",P_MAXVAL(src))
    !call tkeb%reynolds_avg(rhs(:,:,:,mom_index),rhs_mn) 
    !if (nrank==0) then
    !    open(1,file='./src_ru.txt')
    !    write (1,*) qtilde, src,rhs_mn
    !    close(1)
    !endif
    
    ! Y momentum source:
    call tkeb%favre_avg(rho,vpp*vpp, Rij)
    call my_ddy(decomp,dy,rbar*(vtilde*vtilde+Rij),src)
    do j=iy1,iyn 
        rhs(:,j,:,mom_index+1) = rhs(:,j,:,mom_index+1) + src(1,j,1)
    enddo
    !call message(4,"Maximum src momy",P_MAXVAL(src))
    !call tkeb%reynolds_avg(rhs(:,:,:,mom_index+1),rhs_mn) 
    !if (nrank==0) then
    !    open(1,file='./src_rv.txt')
    !    write (1,*) vtilde, src,rhs_mn
    !    close(1)
    !endif
    
    
    ! Z momentum source:
    !call tkeb%favre_avg_and_fluct(rho, w, qtilde, qpp)
    !call tkeb%favre_avg(rho,qpp*vpp, Rij)
    !call my_ddy(decomp,dy,rbar*(qtilde*vtilde+Rij),src)
    !do j=iy1,iyn 
    !    !print *, j,rhs(1,j,1,mom_index+2),src(1,j,1)
    !    rhs(:,j,:,mom_index+2) = rhs(:,j,:,mom_index+2) + src(1,j,1)
    !enddo
    !call tkeb%reynolds_avg(rhs(:,:,:,mom_index+2),rhs_mn) 
    !call message(4,"Maximum src momz",P_MAXVAL(src))
    !if (nrank==0) then
    !    open(1,file='./src_rw.txt')
    !    write (1,*) qtilde, src,rhs_mn

    !    close(1)
    !endif

    
    ! Energy source: e
    call tkeb%reynolds_avg(rhs(:,:,:,TE_index),rhs_mn) 
    call tkeb%favre_avg_and_fluct(rho, e, qtilde, qpp)
    call tkeb%favre_avg(rho,qpp*vpp, Rij)
    call my_ddy(decomp,dy,rbar*(qtilde*vtilde+Rij),src)
    do j=iy1,iyn 
        rhs(:,j,:,TE_index) = rhs(:,j,:,TE_index) + src(1,j,1)
    enddo
    !call message(4,"Maximum src TE(e)",P_MAXVAL(src))
    !if (nrank==0) then
    !    open(1,file='./src_re.txt')
    !    write (1,*) qtilde, src, rhs_mn
    !    close(1)
    !endif
    ! Energy source: p
    call tkeb%reynolds_avg_and_fluct(p, qtilde, qpp)
    call tkeb%reynolds_avg(qpp*vpp, Rij)
    call my_ddy(decomp,dy,qtilde*vtilde+Rij,src)
    do j=iy1,iyn 
        rhs(:,j,:,TE_index) = rhs(:,j,:,TE_index) + src(1,j,1)
    enddo

    !call tkeb%reynolds_avg(rhs(:,:,:,mom_index+2),rhs_mn) 
    !call message(4,"Maximum src TE(p)",P_MAXVAL(src))
    !if (nrank==0) then
    !    open(1,file='./src_rp.txt')
    !    write (1,*) qtilde, src, rhs_mn
    !    close(1)
    !endif

    deallocate(rbar,qtilde,vtilde,Rij,src,rhs_mn)
    deallocate(qpp,vpp)

    end associate
    endif !if forcing
end subroutine
