module Cavity_data
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero, half, one, two, four, pi, imi, three
    use FiltersMod,       only: filters
    use MultiBlockTopologyMod, only: multiblocktopol
    use decomp_2d,        only: decomp_info, nrank
    use basic_io,         only: read_2d_ascii 
    use reductions,       only: P_MAXVAL,P_MINVAL
    use exits,            only: message, message_min_max
    use mpi

    implicit none
    !!!! NOTE: Make sure to update this data according to the problem !!!!
    integer     :: ns     = 1, ybc1 = 0, ybcn = 0
    real(rkind) :: Lx     = four*pi
    real(rkind) :: Ly     = two
    real(rkind) :: Lz     = four*pi/three
    real(rkind) :: Pr     = 0.71_rkind 
    real(rkind) :: Sc     = 1.0_rkind  
    real(rkind) :: gam    = 1.4_rkind
    real(rkind) :: rho_ref= 1.0_rkind
    real(rkind) :: Tw     = 1.0_rkind
    real(rkind) :: Rgas   = 1.0_rkind
    real(rkind) :: Re     = 3000.0_rkind
    real(rkind) :: Mc     = 1.5_rkind
    real(rkind) :: x1, y1, z1
    real(rkind) :: xn, yn, zn
    logical     :: periodicx = .true., periodicy = .false., periodicz = .true. 
    logical     :: add_pert = .true., xplbc_recycle = .true.
    character(len=clen) :: fname_prefix
    ! Gaussian filter for sponge
    type(filters) :: mygfil
    integer :: n_prof
    real(rkind) :: yloc,val,uq=zero
    real(rkind), allocatable :: y_prof(:), u_prof(:), inp_prfl(:)

contains

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
    integer :: ntf  ! ntf is the parameter represents number of times the filter is applied
    real(rkind) :: dx, dy, dz, filpt, thickT
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: dumT, dumF
    character(len=clen) :: outputfile

    dx = Lx/real(decomp%xsz(1)-1,rkind)
    filpt = 2.00_rkind/dx 
    thickT = real(0.3D0, rkind)
    ntf = 4

    ! Gaussian Filter for right side of domain 
    do i=1,decomp%ysz(1)
       dumT(i,:,:)=half*(one-tanh( (real(decomp%xsz(1)- (decomp%yst(1) - 1 + i - 1), rkind)-filpt) / thickT ))
    end do

   ! To check whether dumT is calculted correctly !!!
    write(outputfile, '(a,i3.3,a)') 'dumT_', nrank, '.dat'
    open(10,file=outputfile,status='unknown')
    do i=1,decomp%ysz(1)
       write(10,'(2(e19.12),1x)') x(i,1,1), dumT(i,1,1)
    end do
    close(10)

    dumF = u
    call filter3D(decomp,mygfil,dumF,ntf,x_bc,y_bc,z_bc)
    u = u + dumT*(dumF-u) 

    dumF = v
    call filter3D(decomp,mygfil,dumF,ntf,x_bc,y_bc,z_bc)
    v = v + dumT*(dumF-v)

    dumF = w
    call filter3D(decomp,mygfil,dumF,ntf,x_bc,y_bc,z_bc)
    w = w + dumT*(dumF-w)

    dumF = p
    call filter3D(decomp,mygfil,dumF,ntf,x_bc,y_bc,z_bc)
    p = p + dumT*(dumF-p)

    dumF = rho
    call filter3D(decomp,mygfil,dumF,ntf,x_bc,y_bc,z_bc)
    rho = rho + dumT*(dumF-rho)

  end subroutine

  subroutine sponge_y(decomp, mygfil, y, Ly, u, v, w, p, rho, x_bc, y_bc,z_bc)
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
    integer :: ntf  ! ntf is the parameter represents number of times the filter is applied
    real(rkind) :: dx, dy, dz, filpt, thickT
    real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: dumT, dumF
    character(len=clen) :: outputfile
    real(rkind) :: y_start, thickness, y_top

    dy = Ly/real(decomp%ysz(2)-1,rkind)
    filpt = 0.08_rkind/dy
    thickT = real(0.9D0, rkind)
    ntf = 4
    
    !y_top = maxval(y)
    !y_start = 2.8_rkind       ! <-- your desired start
    !thickness = 0.05_rkind     ! adjust smoothness
    
    !do i=1,decomp%ysz(2)
    !   dumT(:,i,:) = half * ( one - tanh( ( (y_top - y(:,i,:)) - (y_top -y_start) ) /thickness ) )
    !end do
    ! Gaussian Filter for top
    do i=1,decomp%ysz(2)
       dumT(:,i,:)=half*(one-tanh( (real(decomp%ysz(2)- (decomp%yst(2) - 1 + i - 1), rkind)-filpt) / thickT ))
    end do

   ! write(outputfile, '(a,i3.3,a)') 'dumT_', nrank, '.dat'
   ! open(10,file=outputfile,status='unknown')
   ! do i=1,decomp%ysz(2)
   !    write(10,'(2(e19.12),1x)') y(1,i,1), dumT(1,i,1)
   ! end do
   ! close(10)

    dumF = u
    call filter3D(decomp,mygfil,dumF,ntf,x_bc,y_bc,z_bc)
    u = u + dumT*(dumF-u)

    dumF = v
    call filter3D(decomp,mygfil,dumF,ntf,x_bc,y_bc,z_bc)
    v = v + dumT*(dumF-v)

    dumF = w
    call filter3D(decomp,mygfil,dumF,ntf,x_bc,y_bc,z_bc)
    w = w + dumT*(dumF-w)

    dumF = p
    call filter3D(decomp,mygfil,dumF,ntf,x_bc,y_bc,z_bc)
    p = p + dumT*(dumF-p)

    dumF = rho
    call filter3D(decomp,mygfil,dumF,ntf,x_bc,y_bc,z_bc)
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
       ! concentrate towards the two ends
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
       ! concentrate at arbitrary point
       call GracefulExit("flag = 3 (concentrate at arbitrary point) is incomplete right now",21)
    elseif(ymetric_flag==10) then
       ! finite-difference evaluation of metrics (reduces order of accuracy)
       call GracefulExit("flag = 4 (finite-difference evaluation of metrics) is incomplete right now",21)
    endif
    end subroutine
         
    subroutine perturb_potential_v2(gp,x,y,z,nx,ny,nz,nxl,nyl,nzl,Lx,Ly,Lz,u,v,w,p,rho,fname_prefix,Re)
        use decomp_2d,        only: nrank
        use decomp_2d_io
        use constants,        only: half,one, two, four, pi
        type(decomp_info), intent(in)               :: gp
        real(rkind), dimension(:,:,:), intent(in)   :: x,y,z
        real(rkind), dimension(:,:,:), intent(inout):: u,v,w,p,rho
        real(rkind), intent(in)                     :: Lx,Ly,Lz,Re
        integer, intent(in)                         :: nx,ny,nz,nxl,nyl,nzl
        character(len=clen), intent(in) :: fname_prefix

        character(len=clen) :: tempname, fname
        real(rkind), dimension(:,:,:), allocatable :: uperturb, vperturb, wperturb,pperturb,rperturb
        real(rkind) :: amp
        real(rkind), dimension(:,:,:), allocatable :: rand_u      
 
        call message(0,"Before adding perturbations")
        call message(2,"Maximum u", P_MAXVAL(u))
        call message(2,"Maximum v", P_MAXVAL(v))
        call message(2,"Maximum w", P_MAXVAL(w))
        call message(2,"Minimum u", P_MINVAL(u))
        call message(2,"Minimum v", P_MINVAL(v))
        call message(2,"Minimum w", P_MINVAL(w))
      
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

        u = u + uperturb 
        v = v + vperturb 
        w = w + wperturb
        !p = p + pperturb 
        !rho = rho + rperturb

        !deallocate(uperturb)
        !deallocate(vperturb)
        !deallocate(wperturb)
        !deallocate(pperturb)
        !deallocate(rperturb)

        !allocate(rand_u(nxl, nyl, nzl))
        !call random_number(rand_u)
        !rand_u = two*rand_u - one
        !u = u + u*0.5*rand_u 
        !deallocate(rand_u)

        call message(0,"After adding perturbations")
        call message(2,"Maximum u", P_MAXVAL(u))
        call message(2,"Maximum v", P_MAXVAL(v))
        call message(2,"Maximum w", P_MAXVAL(w))
        call message(2,"Minimum u", P_MINVAL(u))
        call message(2,"Minimum v", P_MINVAL(v))
        call message(2,"Minimum w", P_MINVAL(w))
    end subroutine
       
    subroutine interp_profile(y_prof, u_prof, n, yq, uq)
        implicit none
    
        integer, intent(in) :: n
        real(rkind), intent(in) :: y_prof(n), u_prof(n), yq
        real(rkind), intent(out) :: uq
        integer :: i
    
        if (yq <= y_prof(1)) then
            uq = 0.0_rkind
        endif
    
        if (yq >= y_prof(n)) then
            uq = 1.0_rkind
        endif
    
        do i = 1, n-1
            if (yq >= y_prof(i) .and. yq <= y_prof(i+1)) then
                uq = u_prof(i) + (u_prof(i+1)-u_prof(i)) * (yq - y_prof(i)) / (y_prof(i+1)-y_prof(i))
            endif
        end do
    
    end subroutine

end module


subroutine meshgen(decomp, dx, dy, dz, mesh, inputfile, xmetric, ymetric, zmetric, xi, eta, zeta, dxs, dys, dzs, xbuf, zbuf)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info, nrank, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use Cavity_data

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
    integer :: i,j,k,ioUnit, nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn, nxl,nyl,nzl
    integer     ::  xmetric_flag,  ymetric_flag, zmetric_flag
    real(rkind), allocatable, dimension(:,:) :: metric_params
    character(len=clen) :: outputfile,str

    namelist /PROBINPUT/ ns, Lx, Ly, Lz, y1, Pr, Sc, gam, rho_ref, Tw, Re, Mc, add_pert, fname_prefix, xplbc_recycle
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

    ! Local domain sizes
    nxl = decomp%ysz(1);    nyl = decomp%ysz(2);   nzl = decomp%ysz(3)

    ! If base decomposition is in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)

    ! Need to set x, y and z as well as  dx, dy and dz
    ! For now no stretching
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        if (nrank == 0) then
            print *, "Domain size: ",Lx,Ly,Lz
        end if

        dx = Lx/real(nx-1,rkind)  !not periodic
        dy = Ly/real(ny-1,rkind)  !not periodic
        dz = Lz/real(nz-0,rkind)  !periodic    
        x1 = 0._rkind;                               z1 = 0._rkind
        xn = Lx;              yn =  y1 + Ly;         zn = Lz

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k) = x1 + real( ix1 - 1 + i - 1, rkind ) * dx
                    y(i,j,k) = y1 + real( iy1 - 1 + j - 1, rkind ) * dy
                    z(i,j,k) = z1 + real( iz1 - 1 + k - 1, rkind ) * dz
                end do
            end do
        end do

    
    if(nrank == 0) then
      print*, '>>alpha=',metric_params(2,1), '>>beta=',metric_params(2,2), '>>ystart=',metric_params(2,3), '>>yh=', metric_params(2,4)
      print*, '>>ymetric=',ymetric, '>>yflag=',ymetric_flag
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
    do i = 1, decomp%xsz(1)-1
       xtmp2(i,:,:) =  xtmp1(i+1,:,:) - xtmp1(i,:,:)
    end do
    !xtmp2(nx,:,:) = xn - xtmp1(nx,:,:)                                    ! would have worked if x was periodic
    xtmp2(nx,:,:) = xtmp1(decomp%xsz(1),:,:) - xtmp1(decomp%xsz(1)-1,:,:)  ! not periodic; ensures dxs is not zero anywhere
    call transpose_x_to_y(xtmp2,dxs,decomp)   ! Decomposition in x

    dys(:,1,:) = y(:,2,:) - y(:,1,:)       ! Base decomposition in Y
    do j=2, decomp%ysz(2)-1
       dys(:,j,:) =  (y(:,j+1,:) - y(:,j-1,:))/2
    end do
    dys(:,decomp%ysz(2),:) = y(:,decomp%ysz(2),:) - y(:,decomp%ysz(2)-1,:)       ! Base decomposition in Y

    call transpose_y_to_z(z,ztmp1,decomp)   ! Decomposition in z
    do k = 1, nz-1
       ztmp2(:,:,k) =  ztmp1(:,:,k+1) - ztmp1(:,:,k)
    end do
    ztmp2(:,:,nz) = zn - ztmp1(:,:,nz)
    call transpose_z_to_y(ztmp2,dzs,decomp)   ! Decomposition in x

   ! !! Write grid width to a file
   ! write(outputfile, '(a,i0,a)') 'grid_x_', nrank, '.dat'
   ! open(11,file=outputfile,status='unknown')
   ! do i=1,decomp%ysz(1)
   !    write(11,'(2(e19.12),1x)') x(i,1,1), dxs(i,1,1)
   ! enddo
   ! close(11)

   ! if(nrank==0) then
   !   write(outputfile, '(a)') 'grid_y.dat'
   !   open(10,file=outputfile,status='unknown')
   !   do j=1,decomp%ysz(2)
   !      write(10,'(2(e19.12),1x)') y(1,j,1), dys(1,j,1)
   !   enddo
   !   close(10)
   ! endif

   ! write(outputfile, '(a,i0,a)') 'grid_z_', nrank, '.dat'
   ! open(13,file=outputfile,status='unknown')
   ! do k=1,decomp%ysz(3)
   !    write(13,'(2(e19.12),1x)') z(1,1,k), dzs(1,1,k)
   ! enddo
   ! close(13)
    end associate
    nullify(xtmp1)
    nullify(xtmp2)
    nullify(ztmp1)
    nullify(ztmp2)

end subroutine


subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tsim,tstop,dt,tviz)
    use kind_parameters,             only: rkind, clen
    use constants,                   only: zero,half,one,two,four,five,pi,eight, three
    use CompressibleGrid,            only: rho_index,u_index,v_index,w_index,&
                                           p_index,T_index,e_index,Ys_index
    use decomp_2d,                   only: decomp_info,nrank
    use MixtureEOSMod,               only: mixture
    use IdealGasEOS,                 only: idealgas
    use PowerLawViscosityMod,        only: powerLawViscosity
    use ConstRatioBulkViscosityMod,  only: constRatioBulkViscosity
    use ConstPrandtlConductivityMod, only: constPrandtlConductivity
    use ConstSchmidtDiffusivityMod,  only: constSchmidtDiffusivity
    use reductions,                  only: P_MAXVAL,P_MINVAL
    use exits,                       only: GracefulExit, message, nancheck
    use random,                      only: gaussian_random                  

    use Cavity_data
    use mpi

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind),                     intent(inout) :: tsim, tstop, dt, tviz

    type(powerLawViscosity) :: shearvisc
    type(constRatioBulkViscosity) :: bulkvisc
    type(constPrandtlConductivity) :: thermcond
    real(rkind) :: S, Sk, T0, var, mu_ref, umax, eta
    integer :: i,j, k, iounit, nx, ny, nz, nxl, nyl, nzl
    character(len=clen) :: outputfile
    real(rkind), dimension(decomp%ysz(1)) :: x_new
    real(rkind), dimension(decomp%ysz(2)) :: y_new
    real(rkind), dimension(decomp%ysz(3)) :: z_new

    
    namelist /PROBINPUT/ ns, Lx, Ly, Lz, y1, Pr, Sc, gam, rho_ref, Tw, Re, Mc, add_pert, fname_prefix, xplbc_recycle

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

    n_prof = 30
    allocate(y_prof(n_prof), u_prof(n_prof))
    
    open(10, file="input_profile.dat", status="old")
    
    do i = 1, n_prof
        read(10, *) y_prof(i), u_prof(i)
    end do
    
    close(10) 
    ! Global domain sizes
    nx = decomp%xsz(1);     ny = decomp%ysz(2);    nz = decomp%zsz(3)

    ! Local domain sizes
    nxl = decomp%ysz(1);    nyl = decomp%ysz(2);   nzl = decomp%ysz(3)
    
    associate( rho => fields(:,:,:,rho_index), u  => fields(:,:,:,u_index),&
                 v => fields(:,:,:,  v_index), w  => fields(:,:,:,w_index),&
                 p => fields(:,:,:,  p_index), T  => fields(:,:,:,T_index),&
                 e => fields(:,:,:,  e_index), &
                Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1), &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        if (mix%ns /= ns) call GracefulExit("Wrong number of species. Check your input file and make ns consistent with the problem file.",4562)
        Rgas = one/(gam*(Mc**two))
        mu_ref = one/Re
        !print*,'Re=',Re
        !!!! Set each material's transport coefficient object
        shearvisc = powerLawViscosity( mu_ref, Tw, 0.7_rkind)
        bulkvisc  = constRatioBulkViscosity( zero )
        thermcond = constPrandtlConductivity( Pr )
        call mix%set_material( 1, idealgas( gam, Rgas ), shearvisc = shearvisc, bulkvisc  = bulkvisc, thermcond = thermcond  )

        Ys(:,:,:,1)  = one 
        call mix%update(Ys)     
        allocate(inp_prfl(nyl))    
        ! Add base flow profiles
        !u = 1.0d0 !(1-y**2)
        do k = 1, nzl
          do j = 1, nyl
        
            yloc = y(1, j, k)
        
            call interp_profile(y_prof, u_prof, n_prof, yloc, uq)
            do i = 1, nxl
              u(i,j,k) = uq
            end do
            inp_prfl(j)=uq 
          end do
        end do

!        do k=1,nzl
!          do j=1,nyl
!            do i=1,nxl
!              !! laminar inflow
!              eta = y(i,j,k)/1.4_rkind
!              if(y(i,j,k) < zero) then
!                u(i,j,k) = zero
!              elseif(eta<=1.0_rkind) then
!                !u(i,j,k) = 1.8*eta - 1.2*eta**2 + 0.4*eta**3
!                
!                u(i,j,k) = 
!              else
!                u(i,j,k) = 1.0_rkind
!              endif
!              !!T(i,j,k) = (1.5_rkind*(1-y(i,j,k)**4)*(gam-1)*Pr*(Mc**2)/three) +  Tw
!
!              ! turbulent inflow
!              !var = (1-abs(y(i,j,k)))*Re
!              !if (var .lt. 10) then
!              !    u(i,j,k) = var
!              !else
!              !    u(i,j,k) = 2.5_rkind*log(var) + 5.5_rkind
!              !endif
!            end do
!          end do
!        end do

        ! turbulent inflow
        !umax = p_maxval(u)
        !u = u/umax * 1.8d0

        v   = zero
        w   = zero
        rho = rho_ref  
        T   = Tw 
        p   = rho*Rgas*T

        if (add_pert) then
            call perturb_potential_v2(decomp,x,y,z,nx,ny,nz,nxl,nyl,nzl,Lx,Ly,Lz,u,v,w,p,rho,fname_prefix,Re)
        endif

        T = p/(rho*Rgas)     
        call message(2,"Maximum u-velocity",P_MAXVAL(u))

        ! Initialize gaussian filter mygfil
        call mygfil%init(decomp, periodicx, periodicy, periodicz, "gaussian", "gaussian", "gaussian" )
    end associate
    deallocate(y_prof,u_prof)
end subroutine

subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use MixtureEOSMod,    only: mixture
    use reductions,       only: P_MEAN
    use Cavity_data

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

        !write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/Channel_t.dat"

        !open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        !do
        !close(outputunit)
    end associate
end subroutine


subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc,newTimeStep, time_step, xplbc, xplbcInflow, numtbc, tbcIn, xplbcin_type, useMultiBlock, mbtopology)
    use kind_parameters,  only: rkind
    use decomp_2d,        only: decomp_info, nrank, transpose_y_to_x, transpose_x_to_y
    use constants,        only: zero, half, one, two, three, four, five, six, seven, eight
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use MixtureEOSMod,    only: mixture
    use operators,        only: filter3D
    use Cavity_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc
    logical,                         intent(in)    :: newTimeStep
    integer,                         intent(in)    :: time_step 
    logical,                         intent(in)    :: xplbc
    real(rkind), dimension(:,:,:,:), intent(in)    :: xplbcInflow
    integer,                         intent(in)    :: numtbc
    real(rkind), dimension(:),       intent(in)    :: tbcIn
    integer,                         intent(in)    :: xplbcin_type
    logical, optional,               intent(in)    :: useMultiBlock
    type(multiblocktopol), optional, intent(in)    :: mbtopology

    integer :: i, j, k, nx, ny, nz, ix1_new, iy1_new, iz1_new, tidx, ncycles
    integer :: ist, ien, jlo, jst, jen, kst, ken, imb, i_intbd, ttind
    real(rkind) :: dx, dy, dz,rad, filpt, thickT, U0, P0, rho0, T0, Rgas_Tw, alpf, tbcmax
    real(rkind) :: umin, pmin, Tmin, rhomin, diff_u, diff_rho, diff_T, diff_p, onemalpf
    real(rkind) :: umax, vmax, wmax, rmax, pmax
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

        Rgas_Tw = mix%material(1)%mat%Rgas * Tw


        ! set Dirichlet BC at the inlet
        if(decomp%yst(1) == 1) then 
          if(xplbc) then
             ! use x-plane boundary condition from a previous simulation
             ! temporal interpolation index and factor
             if(xplbc_recycle) then
                 tbcmax = maxval(tbcIn)
                 ncycles = floor(tsim/tbcmax)
                 ttind = minloc(abs(tbcIn-(tsim-ncycles*tbcmax)), 1)
             else
                 ttind = minloc(abs(tbcIn-tsim),1)
                 if(tbcIn(ttind) > tsim) ttind = ttind-1
             endif
             if(ttind==0) then
                 ! tsim is before smallest time where inflow is read in
                 ttind = ttind+1
                 alpf = zero
             elseif(ttind==numtbc) then
                 ! tsim is after largest time where inflow is read in
                 ttind = ttind-1
                 alpf = one
             else
                 ! tsim is within the range of times where inflow is read in
                 alpf = (tsim - tbcIn(ttind)) / (tbcIn(ttind+1) - tbcIn(ttind))
             endif
             onemalpf = one - alpf

             if(xplbcin_type==1) then
               do k = 1, decomp%ysz(3) 
                 do j = 1, decomp%ysz(2)
                   u(1,j,k)   =  onemalpf * xplbcInflow(j,k,ttind,1) + alpf * xplbcInflow(j,k,ttind+1,1)
                   v(1,j,k)   =  onemalpf * xplbcInflow(j,k,ttind,2) + alpf * xplbcInflow(j,k,ttind+1,2)
                   w(1,j,k)   =  onemalpf * xplbcInflow(j,k,ttind,3) + alpf * xplbcInflow(j,k,ttind+1,3)
                   p(1,j,k)   =  onemalpf * xplbcInflow(j,k,ttind,4) + alpf * xplbcInflow(j,k,ttind+1,4)
                   rho(1,j,k) =  onemalpf * xplbcInflow(j,k,ttind,5) + alpf * xplbcInflow(j,k,ttind+1,5)
                 enddo
               enddo
             elseif(xplbcin_type==2) then
               do k = 1, decomp%ysz(3) 
                 do j = 1, decomp%ysz(2)
                   u(1,j,k)   =  u(1,j,k)   + onemalpf * xplbcInflow(j,k,ttind,1) + alpf * xplbcInflow(j,k,ttind+1,1)
                   v(1,j,k)   =  v(1,j,k)   + onemalpf * xplbcInflow(j,k,ttind,2) + alpf * xplbcInflow(j,k,ttind+1,2)
                   w(1,j,k)   =  w(1,j,k)   + onemalpf * xplbcInflow(j,k,ttind,3) + alpf * xplbcInflow(j,k,ttind+1,3)
                   p(1,j,k)   =  p(1,j,k)   + onemalpf * xplbcInflow(j,k,ttind,4) + alpf * xplbcInflow(j,k,ttind+1,4)
                   rho(1,j,k) =  rho(1,j,k) + onemalpf * xplbcInflow(j,k,ttind,5) + alpf * xplbcInflow(j,k,ttind+1,5)
                 enddo
               enddo
             endif
 
             umax = maxval(abs(u)); vmax = maxval(abs(v));  wmax = maxval(abs(w)); 
             pmax = maxval(abs(p)); rmax = maxval(abs(rho)); 
             if(nrank==0) print '(a,e19.12,1x,a,i6.6,1x,a,i4.4,1x,a,5(e19.12,1x))', 'x-inflow bc tsim= ', tsim, ' ttind=', ttind, 'rank=', nrank, "uvwpr=", umax, vmax, wmax, pmax, rmax
          else
               do k = 1, decomp%ysz(3) 
                 do j = 1, decomp%ysz(2)
                   u(1,j,k)   =  inp_prfl(j)
                   rho(1,j,k) =  rho_ref
                   !T(1,j,k)   =  Tw
                   v(1,j,k)   = zero
                   w(1,j,k)   = zero
                   p(1,j,k)   =  rho_ref * Rgas_Tw
                 enddo
               enddo
          endif
        endif
        ! set Dirichlet BC at top and bottom
        do k = 1,decomp%ysz(3) 
           u(:,1,k) = zero;                  !u(:,decomp%ysz(2),k) = zero
           v(:,1,k) = zero;                  !v(:,decomp%ysz(2),k) = zero
           w(:,1,k) = zero;                  !w(:,decomp%ysz(2),k) = zero
           !T(:,1,k) = Tw;                   T(:,decomp%ysz(2),k) = Tw
           p(:,1,k) = rho(:,1,k)*Rgas_Tw;    !p(:,decomp%ysz(2),k) = rho(:,decomp%ysz(2),k)*Rgas_Tw
        end do
        if(present(useMultiBlock)) then
         if(useMultiBlock) then
            ! set Dirichlet BC at bottom block of multiblock 
            do imb = 1, mbtopology%y_num_blocks
              !jlo = mbtopology%yst(2, imb)
              !ist = mbtopology%yst(1, imb);   ien = mbtopology%yen(1, imb)
              !kst = mbtopology%yst(3, imb);   ken = mbtopology%yen(3, imb)
              !do k = kst, ken
              !    u(ist:ien, jlo, k) = zero
              !    v(ist:ien, jlo, k) = zero
              !    w(ist:ien, jlo, k) = zero
              !    T(ist:ien, jlo, k) = Tw
              !enddo

              ! lower boundary (jst, jen should be the same)
              ist = mbtopology%ybclo_st(1,imb); jst = mbtopology%ybclo_st(2,imb); kst = mbtopology%ybclo_st(3,imb)
              ien = mbtopology%ybclo_en(1,imb); jen = mbtopology%ybclo_en(2,imb); ken = mbtopology%ybclo_en(3,imb)
              do k = kst, ken
                  u(ist:ien, jst, k) = zero
                  v(ist:ien, jst, k) = zero
                  w(ist:ien, jst, k) = zero
                  !T(ist:ien, jst, k) = Tw
                  p(ist:ien, jst, k) = rho(ist:ien, jst, k) * Rgas_Tw
              enddo

              ! upper boundary (jst, jen should be the same)
              !ist = mbtopology%ybchi_st(1,imb); jst = mbtopology%ybchi_st(2,imb); kst = mbtopology%ybchi_st(3,imb)
              !ien = mbtopology%ybchi_en(1,imb); jen = mbtopology%ybchi_en(2,imb); ken = mbtopology%ybchi_en(3,imb)
!              do k = kst, ken
!                  u(ist:ien, jst, k) = zero
!                  v(ist:ien, jst, k) = zero
!                  w(ist:ien, jst, k) = zero
!                  !T(ist:ien, jst, k) = Tw
!                  p(ist:ien, jst, k) = rho(ist:ien, jst, k) * Rgas_Tw
!              enddo
!
              !print *, 'Num-internal-boundaries-left: nrank=', nrank, 'num_blocks=', mbtopology%y_num_blocks, 'num_int_bdries=',mbtopology%y_num_intbd_left
              ! left internal boundary (ist, ien should be the same)
              do i_intbd = 1, mbtopology%y_num_intbd_left(imb)
                ist = mbtopology%y_intbd_left_st(1,i_intbd) 
                jst = mbtopology%y_intbd_left_st(2,i_intbd)
                kst = mbtopology%y_intbd_left_st(3,i_intbd)
                ien = mbtopology%y_intbd_left_en(1,i_intbd)
                jen = mbtopology%y_intbd_left_en(2,i_intbd)
                ken = mbtopology%y_intbd_left_en(3,i_intbd)
                do k = kst, ken
                 do j = jst, jen
                    u(ist, j, k) = zero
                    v(ist, j, k) = zero
                    w(ist, j, k) = zero
                    !T(ist, j, k) = Tw
                    p(ist, j, k) = rho(ist, j, k) * Rgas_Tw
                 enddo
                enddo
              enddo

              ! right internal boundary (ist, ien should be the same)
              do i_intbd = 1, mbtopology%y_num_intbd_rght(imb)
                ist = mbtopology%y_intbd_rght_st(1,i_intbd) 
                jst = mbtopology%y_intbd_rght_st(2,i_intbd)
                kst = mbtopology%y_intbd_rght_st(3,i_intbd)
                ien = mbtopology%y_intbd_rght_en(1,i_intbd)
                jen = mbtopology%y_intbd_rght_en(2,i_intbd)
                ken = mbtopology%y_intbd_rght_en(3,i_intbd)
                do k = kst, ken
                 do j = jst, jen
                    u(ist, j, k) = zero
                    v(ist, j, k) = zero
                    w(ist, j, k) = zero
                    !T(ist, j, k) = Tw
                    p(ist, j, k) = rho(ist, j, k) * Rgas_Tw
                 enddo
                enddo
              enddo
            enddo
         endif
        endif
        !u(1,:,:)=one
        !p   = rho*Rgas*T

        !!!!! =============  Add Sponge+bulk for exit bc ==========!!!!!
        ! Gradually apply the exit boundary conditions
        ! Apply sponge in X-direction on right
        call  sponge_x(decomp, mygfil, x, Lx, u, v, w, p, rho, x_bc, y_bc, z_bc)
        call  sponge_y(decomp, mygfil, y, Ly, u, v, w, p, rho, x_bc, y_bc, z_bc)

    end associate
end subroutine


subroutine hook_timestep(decomp,der,dx,dy,dz,mesh,fields,mix,step,tsim,outputdir,sgsmodel)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,two
    use CompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info, nrank
    use DerivativesMod,     only: derivatives
    use MixtureEOSMod,    only: mixture
    use sgsmod_cgrid,     only: sgs_cgrid
    use exits,            only: message
    use reductions,       only: P_MAXVAL,P_MINVAL

    use Cavity_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    real(rkind),                     intent(in) :: dx,dy,dz
    type(mixture),                   intent(in) :: mix
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    character(len=*),                intent(in) :: outputdir
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(sgs_cgrid), optional,       intent(in) :: sgsmodel

    integer :: ny  , j, my_step = 0
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
        call message(2,"Maximum u-velocity",P_MAXVAL(u))
        call message(2,"Maximum v-velocity",P_MAXVAL(v))
        call message(2,"Maximum pressure",P_MAXVAL(p))
        call message(2,"Maximum density",P_MAXVAL(rho))
        call message(2,"Maximum temperature",P_MAXVAL(T))
        call message(2,"Maximum shear viscosity",P_MAXVAL(mu))
        call message(2,"Maximum bulk viscosity",P_MAXVAL(bulk))
        call message(2,"Maximum conductivity",P_MAXVAL(kap))
        call message(2,"Maximum diffusivity",P_MAXVAL(diff))


        !!!if(useSGS)
         ! !!if(sgsmodel%DynamicProcedureType==1) then
         !    call message_min_max(2,"Bounds for LD-Coeff-tke  : ",     &
         !           sgsmodel%get_Max_LocalDynamicProcedure_Coeff_tke(),   &
         !           sgsmodel%get_Min_LocalDynamicProcedure_Coeff_tke())
         !    call message_min_max(2,"Bounds for LD-Coeff      : ",     &
         !           sgsmodel%get_Max_LocalDynamicProcedure_Coeff(),   &
         !           sgsmodel%get_Min_LocalDynamicProcedure_Coeff())
         !    call message_min_max(2,"Bounds for LD-Coeff-Qjsgs: ",     &
         !           sgsmodel%get_Max_LocalDynamicProcedure_Coeff_Qjsgs(),   &
         !           sgsmodel%get_Min_LocalDynamicProcedure_Coeff_Qjsgs())
         ! !endif

         !  if(mod(step,1000)==0) then
         !     my_step = my_step + 1
         !     cmodel_loc = sgsmodel%get_LocalDynamicProcedure_Coeff()
         !     cmodel_loc_Qjsgs = sgsmodel%get_LocalDynamicProcedure_Coeff_Qjsgs()
         !     cmodel_loc_tke   = sgsmodel%get_LocalDynamicProcedure_Coeff_tke()
         !     if(nrank==0) then
         !         write(outputfile, '(a,i7.7,a)') 'cmodel_', step, '.dat'
         !         open(10,file=outputfile,status='unknown')
         !         do j=1,decomp%ysz(2)
         !            write(10,'(4(e19.12),1x)') y(1,j,1), cmodel_loc(j), cmodel_loc_Qjsgs(j), cmodel_loc_tke(j)
         !            !write(10,'(3(e19.12),1x)') y(1,j,1), cmodel_loc(j), cmodel_loc_Qjsgs(j)
         !         enddo
         !         close(10)
         !     endif
         !  endif

       !! endif

    end associate

end subroutine

subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs,der,dt,step,dys)
    use CompressibleGrid,   only: rho_index,u_index,v_index,w_index,&
                                  p_index,T_index,e_index,Ys_index,mu_index
    use kind_parameters,    only: rkind
    use AveragingMod,       only: averaging
    use DerivativesMod,     only: derivatives
    use constants,          only: one, zero, two
    use decomp_2d,          only: decomp_info,nrank
    use MixtureEOSMod,      only: mixture
    use reductions,         only: P_MAXVAL,P_MINVAL
    use Cavity_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(derivatives),               intent(in)    :: der
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim,dt
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    integer,                         intent(in)    :: step
    real(rkind), dimension(:,:,:),       intent(in)    :: dys

    !integer :: mass_index, mom_index, TE_index, i, j, k, ioUnit, nxl, nyl, nzl, nx, ny, nz, mpi_ierr, ierr
    !real(rkind) :: f_src = 0._rkind, q0_flux = 2.0_rkind, mu_bar, mu_locsum, mu_globsum, alpha, beta, q_flux_old, q_flux_new,u_bulk
    !real(rkind), allocatable, dimension(:)     :: u_locsum, u_globsum, rhou_locsum, rhou_globsum, rho_locsum, rho_globsum,dys_1d
    !real(rkind), allocatable, dimension(:,:,:) :: du_bardy, u_bar
    !integer, dimension(3) :: st
    !integer :: my_step=0
    !character(len=clen) :: outputfile
    !real(rkind) :: src


    !associate( rho => fields(:,:,:,rho_index), u  => fields(:,:,:,u_index),&
    !             v => fields(:,:,:,  v_index), w  => fields(:,:,:,w_index),&
    !             p => fields(:,:,:,  p_index), T  => fields(:,:,:,T_index),&
    !             e => fields(:,:,:,  e_index), mu => fields(:,:,:,mu_index ),&
    !            Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1), &
    !             y => mesh(:,:,:,2) )
    !    ! Set mass, momentum and energy indices in Wcnsrv
    !    mass_index = 1
    !    mom_index  = mass_index + ns
    !    TE_index   = mom_index + 3
    !  
    !    ! Global domain sizes
    !    nx = decomp%xsz(1);     ny = decomp%ysz(2);    nz = decomp%zsz(3)
    !    ! Local domain sizes
    !    nxl = decomp%ysz(1);    nyl = decomp%ysz(2);   nzl = decomp%ysz(3)
    !   
    !    allocate(u_locsum(nyl));      allocate(u_globsum(nyl)) ;  allocate(dys_1d(nyl))  
    !    allocate(rhou_locsum(nyl));   allocate(rhou_globsum(nyl)) ;  
    !    allocate(rho_locsum(nyl));    allocate(rho_globsum(nyl)) ;  
    !    allocate(u_bar(nxl,nyl,nzl)); allocate(du_bardy(nxl,nyl,nzl))
 
    !    dys_1d = dys(1,:,1)
    !    
    !    mu_locsum = 0
    !    do k=1,nzl
    !       do i=1,nxl
    !          mu_locsum = mu_locsum + mu(i,1,k)
    !       end do
    !    end do
    !    call mpi_allreduce(mu_locsum, mu_globsum, 1, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
    !    mu_bar  = mu_globsum/(nx*nz) 
    !    !print*, nrank, mu(1,1,1), mu_locsum, mu_globsum, mu_bar       

    !    u_locsum = 0; rhou_locsum = 0; rho_locsum = 0
    !    do k=1,nzl
    !      do j=1,nyl
    !       do i=1,nxl
    !            u_locsum(j)    = u_locsum(j)    + u(i,j,k)
    !            rhou_locsum(j) = rhou_locsum(j) + u(i,j,k)*rho(i,j,k)
    !            rho_locsum(j)  = rho_locsum(j)  + rho(i,j,k)
    !       end do
    !      end do
    !    end do
    !    call mpi_allreduce(u_locsum, u_globsum, nyl, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
    !    call mpi_allreduce(rhou_locsum, rhou_globsum, nyl, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
    !    call mpi_allreduce(rho_locsum, rho_globsum, nyl, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
    !    do k=1,nzl
    !      do j=1,nyl
    !       do i=1,nxl
    !            u_bar(i,j,k) = u_globsum(j)/(nx*nz)
    !       end do
    !      end do
    !    end do
    !    
    !    call der%ddy(u_bar,du_bardy,ybc1,ybcn) 
    !    
    !    q_flux_old = sum( (rhou_globsum/(nx*nz)) * dys_1d)
    !    q_flux_new = q_flux_old - dt*( Ly*Lz*f_src + two*Lz*mu_bar*du_bardy(1,1,1))
    !    alpha  = two/dt;  beta = -0.2_rkind/dt
    !    f_src  = f_src + (dt/(Ly*Lz))*(alpha*(q_flux_new-q0_flux) + beta*(q_flux_old-q0_flux))
    !    u_bulk = sum( (rhou_globsum/(nx*nz)) * dys_1d)/sum( (rho_globsum/(nx*nz))* dys_1d )

    !    !if(step==0) then
    !    !   my_step = my_step + 1
    !    !   if (nrank==0) then
    !    !       write(outputfile, '(a,i3.3,a)') 'dump_uavg_', my_step, '.dat'
    !    !       open(10,file=outputfile,status='unknown')
    !    !       do i=1,decomp%ysz(2)
    !    !           write(10,'(3(e19.12),1x)') y(1,i,1), u_bar(1,i,1), du_bardy(1,i,1)
    !    !       end do
    !    !       close(10)
    !    !   endif
    !    !endif
    !    
    !    src = -mu_bar*du_bardy(1,1,1)
    !    if (nrank==0)then
    !       print*, '>> Mass flux=',0.5*q_flux_new, '>> Bulk Vel=',u_bulk , '>> f_src=',f_src
    !    endif
    !    ! X momentum source:
    !    rhs(:,:,:,mom_index) = rhs(:,:,:,mom_index) - f_src
    !    !rhs(:,:,:,mom_index) = rhs(:,:,:,mom_index)  - src

    !    ! Energy source: e
    !    !rhs(:,:,:,TE_index) = rhs(:,:,:,TE_index) - f_src*u_bulk
    !    rhs(:,:,:,TE_index) = rhs(:,:,:,TE_index) - f_src*u(:,:,:)
    !    !rhs(:,:,:,TE_index) = rhs(:,:,:,TE_index)  -  src*u(:,:,:)
    !    
    !    deallocate(u_locsum);  deallocate(rhou_locsum); deallocate(dys_1d)
    !    deallocate(u_globsum); deallocate(rhou_globsum)
    !end associate
end subroutine
