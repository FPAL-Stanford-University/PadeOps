module Channel_cvl_data
    use kind_parameters,  only: rkind, mpirkind, clen
    use constants,        only: zero, half, one, two, four, pi, imi, three
    use FiltersMod,       only: filters
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
    real(rkind) :: Re     = 3000.0_rkind
    real(rkind) :: Mc     = 1.5_rkind
    real(rkind) :: x1, y1, z1
    real(rkind) :: xn, yn, zn
    logical     :: periodicx = .true., periodicy = .false., periodicz = .true. 
    logical     :: add_pert = .true.
    character(len=clen) :: fname_prefix
    ! Gaussian filter for sponge
    type(filters) :: mygfil

contains

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
    filpt = 0.08_rkind/dy 
    thickT = real(0.9D0, rkind)
        
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

    subroutine stretched_coordinates(decomp, y, eta, ymetric_flag, param1, param2, param3, param4)
    use constants,        only: zero, half, one
    use decomp_2d,        only: decomp_info, nrank
    use exits,            only: GracefulExit, message, nancheck

    type(decomp_info),               intent(in)    :: decomp
    real(rkind), dimension(:,:,:),   intent(inout) :: y
    real(rkind), dimension(:,:,:  ), intent(inout) :: eta
    integer,                         intent(in)    :: ymetric_flag
    real(rkind),                     intent(in)    :: param1, param2, param3, param4
    integer     :: i,j,k
    real(rkind) :: yfocus, ytau, ystart, yh, alpha, beta
    real(rkind) :: yfocus_adj, num, den, BB, yuniform_adj, BB2

    ! concentrate towards the center -- Pletcher, Tannehill, Anderson
    ! (Section 5.6, Transformation 3, pg. 332) 
    if(ymetric_flag==0) then
       y = eta
    elseif(ymetric_flag==1) then
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
       ! concentrate towards the two ends for alpha=0.5, at right end (final point) for alpha=0
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
       ! concentrate at leftt end, towards initial point
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

    elseif(ymetric_flag==4) then
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
       
end module


subroutine meshgen(decomp, dxi, deta, dzeta, mesh, inputfile, meshcvl, dxs, dys, dzs, xbuf, zbuf)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info, nrank, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use Channel_cvl_data

    implicit none

    type(decomp_info),                                          intent(in)    :: decomp
    real(rkind),                                                intent(inout) :: dxi,deta,dzeta
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh,meshcvl
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:  ), intent(inout) :: dxs, dys, dzs
    real(rkind), dimension(:,:,:,:), target,intent(in):: xbuf, zbuf
    real(rkind), dimension(:,:,:), pointer :: xtmp1, xtmp2, ztmp1, ztmp2
    integer :: i,j,k,ioUnit, nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn, nxl,nyl,nzl
    integer     ::  xmetric_flag=0,  ymetric_flag=0, zmetric_flag=0, curvil_flag = 0
    real(rkind), allocatable, dimension(:,:) :: metric_params
    real(rkind) :: xfocus, xtau, xh, xstart
    real(rkind) :: yfocus, ytau, yh, ystart
    real(rkind) :: zfocus, ztau, zh, zstart
    character(len=clen) :: outputfile,str

    namelist /PROBINPUT/ ns, Lx, Ly, Lz, Pr, Sc, gam, rho_ref, Tw, Re, Mc, add_pert, fname_prefix
    namelist /METRICS/ xmetric_flag, ymetric_flag, zmetric_flag, metric_params, curvil_flag
    
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
    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3), & 
               xi => meshcvl(:,:,:,1), eta => meshcvl(:,:,:,2), zeta => meshcvl(:,:,:,3) )
        if (nrank == 0) then
            print *, "Domain size: ",Lx,Ly,Lz
        end if

        dxi   = Lx/real(nx-0,rkind)  !periodic
        deta  = Ly/real(ny-1,rkind)
        dzeta = Lz/real(nz-0,rkind)  !periodic    

        x1 = 0._rkind;        y1 = -Ly/2._rkind;     z1 = 0._rkind
        xn = Lx;              yn =  Ly/2._rkind;     zn = Lz

        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    xi(i,j,k)   = x1 + real( ix1 - 1 + i - 1, rkind ) * dxi
                    eta(i,j,k)  = y1 + real( iy1 - 1 + j - 1, rkind ) * deta
                    zeta(i,j,k) = z1 + real( iz1 - 1 + k - 1, rkind ) * dzeta
                end do
            end do
        end do

    xfocus = metric_params(1,1);  xtau   = metric_params(1,2);  xstart = metric_params(1,3); xh = metric_params(1,4)
    yfocus = metric_params(2,1);  ytau   = metric_params(2,2);  ystart = metric_params(2,3); yh = metric_params(2,4)
    zfocus = metric_params(3,1);  ztau   = metric_params(3,2);  zstart = metric_params(3,3); zh = metric_params(3,4)
   
    if(nrank == 0) then
      print*, '>>yfocus=',yfocus, '>>ytau=',ytau, '>>ystart=',ystart, '>>yh=',yh, '>>yflag=',ymetric_flag
    endif

    call stretched_coordinates(decomp,x,xi,xmetric_flag,metric_params(1,1),&
                               metric_params(1,2),metric_params(1,3),metric_params(1,4))

    call stretched_coordinates(decomp,y,eta,ymetric_flag,metric_params(2,1),&
                               metric_params(2,2),metric_params(2,3),metric_params(2,4))
    
    call stretched_coordinates(decomp,z,zeta,zmetric_flag,metric_params(3,1),&
                               metric_params(3,2),metric_params(3,3),metric_params(3,4))


    ! Grid width on stretched/uniform mesh
    call transpose_y_to_x(x,xtmp1,decomp)   ! Decomposition in x
    do i = 1, nx-1
       xtmp2(i,:,:) =  xtmp1(i+1,:,:) - xtmp1(i,:,:)
    end do
    xtmp2(nx,:,:) = xn - xtmp1(nx,:,:)
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
    deallocate(metric_params)   

end subroutine


subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tsim,tstop,dt,tviz)
    use kind_parameters,             only: rkind, clen
    use constants,                   only: zero,half,one,two,four,five,pi,eight, three
    use CurvilCompressibleGrid,      only: rho_index,u_index,v_index,w_index,&
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

    use Channel_cvl_data
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
    real(rkind) :: S, Sk, T0, Rgas, var, mu_ref
    integer :: i,j, k, iounit, nx, ny, nz, nxl, nyl, nzl
    character(len=clen) :: outputfile
    real(rkind), dimension(decomp%ysz(1)) :: x_new
    real(rkind), dimension(decomp%ysz(2)) :: y_new
    real(rkind), dimension(decomp%ysz(3)) :: z_new
    
    namelist /PROBINPUT/ ns, Lx, Ly, Lz, Pr, Sc, gam, rho_ref, Tw, Re, Mc, add_pert, fname_prefix

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)

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

        !!!! Set each material's transport coefficient object
        shearvisc = powerLawViscosity( mu_ref, Tw, 0.7_rkind)
        bulkvisc  = constRatioBulkViscosity( zero )
        thermcond = constPrandtlConductivity( Pr )
        call mix%set_material( 1, idealgas( gam, Rgas ), shearvisc = shearvisc, bulkvisc  = bulkvisc, thermcond = thermcond  )

        Ys(:,:,:,1)  = one 
        call mix%update(Ys)     
     
        ! Add base flow profiles
        u = (1-y**2)*1.5

        !do k=1,nzl
        !    do j=1,nyl
        !        do i=1,nxl
        !           u(i,j,k) = 1.5_rkind*(1-y(i,j,k)**2)
        !           !T(i,j,k) = (1.5_rkind*(1-y(i,j,k)**4)*(gam-1)*Pr*(Mc**2)/three) +  Tw
        !           !var = (1-abs(y(i,j,k)))*Re
        !           !if (var .lt. 10) then
        !           !    u(i,j,k) = var
        !           !else
        !           !    u(i,j,k) = 2.5_rkind*log(var) + 5.5_rkind
        !           !endif

        !        end do
        !    end do
        !end do

        v   = zero
        w   = zero
        rho = rho_ref  
        T   = Tw 
        p   = rho*Rgas*T

        if (add_pert) then
            call perturb_potential_v2(decomp,x,y,z,nx,ny,nz,nxl,nyl,nzl,Lx,Ly,Lz,u,v,w,p,rho,fname_prefix,Re)
        endif

        T = p/(rho*Rgas)     

        ! Initialize gaussian filter mygfil
        call mygfil%init(decomp, periodicx, periodicy, periodicz, "gaussian", "gaussian", "gaussian" )
    end associate
end subroutine


subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use CurvilCompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use MixtureEOSMod,    only: mixture
    use reductions,       only: P_MEAN
    use Channel_cvl_data

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

        !write(outputfile,'(2A,I4.4,A)') trim(outputdir),"/Channel_cvl_t.dat"

        !open(unit=outputunit, file=trim(outputfile), form='FORMATTED')
        !do
        !close(outputunit)
    end associate
end subroutine


subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc,newTimeStep, time_step)
    use kind_parameters,  only: rkind
    use decomp_2d,        only: decomp_info, nrank, transpose_y_to_x, transpose_x_to_y
    use constants,        only: zero, half, one, two, three, four, five, six, seven, eight
    use CurvilCompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use MixtureEOSMod,    only: mixture
    use operators,        only: filter3D

    use Channel_cvl_data

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

        !!!!! =============  Add Sponge+bulk for exit bc ==========!!!!!
        ! Apply sponge in Y-direction on top and bottom
        !call  sponge_y(decomp, mygfil, y, Ly, u, v, w, p, rho, x_bc, y_bc, z_bc)

        ! set Dirichlet BC for velocity, Temperature
        do k = 1,decomp%ysz(3) 
           u(:,1,k) = zero;          u(:,decomp%ysz(2),k) = zero
           v(:,1,k) = zero;          v(:,decomp%ysz(2),k) = zero
           w(:,1,k) = zero;          w(:,decomp%ysz(2),k) = zero
           T(:,1,k) = Tw;            T(:,decomp%ysz(2),k) = Tw
        end do
         


    end associate
end subroutine


subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim,sgsmodel)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,two
    use CurvilCompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info, nrank
    use MixtureEOSMod,    only: mixture
    use sgsmod_cgrid,     only: sgs_cgrid
    use exits,            only: message
    use reductions,       only: P_MAXVAL,P_MINVAL

    use Channel_cvl_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(mixture),                   intent(in) :: mix
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(sgs_cgrid), optional,       intent(in) :: sgsmodel

    real(rkind) :: dx, Ythick, oob
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

subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs,der,dt,step,dys,detady)
    use CurvilCompressibleGrid,   only: rho_index,u_index,v_index,w_index,&
                                  p_index,T_index,e_index,Ys_index,mu_index
    use kind_parameters,    only: rkind
    use AveragingMod,       only: averaging
    use DerivativesMod,     only: derivatives
    use constants,          only: one, zero, two
    use decomp_2d,          only: decomp_info,nrank
    use MixtureEOSMod,      only: mixture
    use reductions,         only: P_MAXVAL,P_MINVAL
    use channel_cvl_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(derivatives),               intent(in)    :: der
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim,dt
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
    integer,                         intent(in)    :: step
    real(rkind), dimension(:,:,:),   intent(in)    :: dys,detady

    integer :: mass_index, mom_index, TE_index, i, j, k, ioUnit, nxl, nyl, nzl, nx, ny, nz, mpi_ierr, ierr
    real(rkind) :: f_src = 0._rkind, q0_flux = 2.0_rkind, mu_bar, mu_locsum, mu_globsum, alpha, beta, q_flux_old, q_flux_new,u_bulk
    real(rkind), allocatable, dimension(:)     :: u_locsum, u_globsum, rhou_locsum, rhou_globsum, rho_locsum, rho_globsum,dys_1d
    real(rkind), allocatable, dimension(:,:,:) :: du_bardeta, u_bar, du_bardy
    integer :: my_step=0
    character(len=clen) :: outputfile
    real(rkind) :: src


    associate( rho => fields(:,:,:,rho_index), u  => fields(:,:,:,u_index),&
                 v => fields(:,:,:,  v_index), w  => fields(:,:,:,w_index),&
                 p => fields(:,:,:,  p_index), T  => fields(:,:,:,T_index),&
                 e => fields(:,:,:,  e_index), mu => fields(:,:,:,mu_index ),&
                Ys => fields(:,:,:,Ys_index:Ys_index+mix%ns-1), &
                 y => mesh(:,:,:,2) )
        ! Set mass, momentum and energy indices in Wcnsrv
        mass_index = 1
        mom_index  = mass_index + ns
        TE_index   = mom_index + 3
      
        ! Global domain sizes
        nx = decomp%xsz(1);     ny = decomp%ysz(2);    nz = decomp%zsz(3)
        ! Local domain sizes
        nxl = decomp%ysz(1);    nyl = decomp%ysz(2);   nzl = decomp%ysz(3)
       
        allocate(u_locsum(nyl));      allocate(u_globsum(nyl)) ;  allocate(dys_1d(nyl))  
        allocate(rhou_locsum(nyl));   allocate(rhou_globsum(nyl)) ;  
        allocate(rho_locsum(nyl));    allocate(rho_globsum(nyl)) ;  
        allocate(u_bar(nxl,nyl,nzl)); allocate(du_bardy(nxl,nyl,nzl)); allocate(du_bardeta(nxl,nyl,nzl))
 
        dys_1d = dys(1,:,1)
        
        mu_locsum = 0
        do k=1,nzl
           do i=1,nxl
              mu_locsum = mu_locsum + mu(i,1,k)
           end do
        end do
        call mpi_allreduce(mu_locsum, mu_globsum, 1, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
        mu_bar  = mu_globsum/(nx*nz) 
        !print*, nrank, mu(1,1,1), mu_locsum, mu_globsum, mu_bar       

        u_locsum = 0; rhou_locsum = 0; rho_locsum = 0
        do k=1,nzl
          do j=1,nyl
           do i=1,nxl
                u_locsum(j)    = u_locsum(j)    + u(i,j,k)
                rhou_locsum(j) = rhou_locsum(j) + u(i,j,k)*rho(i,j,k)
                rho_locsum(j)  = rho_locsum(j)  + rho(i,j,k)
           end do
          end do
        end do
        call mpi_allreduce(u_locsum, u_globsum, nyl, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
        call mpi_allreduce(rhou_locsum, rhou_globsum, nyl, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
        call mpi_allreduce(rho_locsum, rho_globsum, nyl, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)
        do k=1,nzl
          do j=1,nyl
           do i=1,nxl
                u_bar(i,j,k) = u_globsum(j)/(nx*nz)
           end do
          end do
        end do
        
        call der%ddeta(u_bar,du_bardeta,ybc1,ybcn)    !du_bardxi = 0; du_bardzeta = 0
        !dfdy = dxidy*dfdxi + detady*dfdeta + dzetady*dfdzeta
        du_bardy = detady * du_bardeta 
        
        q_flux_old = sum( (rhou_globsum/(nx*nz)) * dys_1d)
        q_flux_new = q_flux_old - dt*( Ly*Lz*f_src + two*Lz*mu_bar*du_bardy(1,1,1))
        alpha  = two/dt;  beta = -0.2_rkind/dt
        f_src  = f_src + (dt/(Ly*Lz))*(alpha*(q_flux_new-q0_flux) + beta*(q_flux_old-q0_flux))
        u_bulk = sum( (rhou_globsum/(nx*nz)) * dys_1d)/sum( (rho_globsum/(nx*nz))* dys_1d )

        !if(step==0) then
        !   my_step = my_step + 1
        !   if (nrank==0) then
        !       write(outputfile, '(a,i3.3,a)') 'dump_uavg_', my_step, '.dat'
        !       open(10,file=outputfile,status='unknown')
        !       do i=1,decomp%ysz(2)
        !           write(10,'(3(e19.12),1x)') y(1,i,1), u_bar(1,i,1), du_bardy(1,i,1)
        !       end do
        !       close(10)
        !   endif
        !endif
        
        src = -mu_bar*du_bardy(1,1,1)
        if (nrank==0)then
           print*, '>> Mass flux=',0.5*q_flux_new, '>> Bulk Vel=',u_bulk , '>> f_src=',f_src
        endif
        ! X momentum source:
        rhs(:,:,:,mom_index) = rhs(:,:,:,mom_index) - f_src
        !rhs(:,:,:,mom_index) = rhs(:,:,:,mom_index)  - src

        ! Energy source: e
        !rhs(:,:,:,TE_index) = rhs(:,:,:,TE_index) - f_src*u_bulk
        rhs(:,:,:,TE_index) = rhs(:,:,:,TE_index) - f_src*u(:,:,:)
        !rhs(:,:,:,TE_index) = rhs(:,:,:,TE_index)  -  src*u(:,:,:)
        
        deallocate(u_locsum);    deallocate(u_globsum) ;    deallocate(dys_1d)  
        deallocate(rhou_locsum); deallocate(rhou_globsum) ;  
        deallocate(rho_locsum);  deallocate(rho_globsum) ;  
        deallocate(u_bar);       deallocate(du_bardy);      deallocate(du_bardeta)
    end associate
end subroutine
