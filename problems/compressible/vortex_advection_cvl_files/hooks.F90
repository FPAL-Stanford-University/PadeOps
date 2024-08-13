module Vortex_advection_cvl_data
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
    real(rkind) :: Lx     = 15.0_rkind
    real(rkind) :: Ly     = 15.0_rkind
    real(rkind) :: p_inf  = 1.0_rkind
    real(rkind) :: rho_inf= 1.0_rkind
    real(rkind) :: M_inf  = 0.1_rkind
    real(rkind) :: gam    = 1.4_rkind
    real(rkind) :: Rgas   = 1.0_rkind
    real(rkind) :: R_v    = 1.0_rkind
    ! Gaussian filter for sponge
    !type(filters) :: mygfil

contains

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

end module

subroutine meshgen(decomp, dxi, deta, dzeta, mesh, inputfile, meshcvl, dxs, dys, dzs, xbuf, zbuf)
    use kind_parameters,  only: rkind
    use constants,        only: half,one
    use decomp_2d,        only: decomp_info, nrank, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use Vortex_advection_cvl_data

    implicit none

    type(decomp_info),               intent(in)    :: decomp
    real(rkind),                     intent(inout) :: dxi,deta,dzeta
    real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
    character(len=*),                intent(in)    :: inputfile
    real(rkind), dimension(:,:,:,:), intent(inout) :: meshcvl
    real(rkind), dimension(:,:,:  ), intent(inout) :: dxs, dys, dzs
    real(rkind), dimension(:,:,:,:), target,intent(in):: xbuf, zbuf
    real(rkind), dimension(:,:,:), pointer :: xtmp1, xtmp2, ztmp1, ztmp2
    integer :: i,j,k,ioUnit, nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn
    real(rkind) :: x1,xn,y1,yn,z1,zn
    real(rkind) :: Ampx = 1.0_rkind, Ampy = 2.0_rkind, num = 6, phi = 0.00_rkind, ome_t = 0.25_rkind
    character(len=clen) :: outputfile,str
    integer     ::  xmetric_flag=0,  ymetric_flag=0, zmetric_flag=0, curvil_flag = 0
    real(rkind) :: xfocus, xtau, xh, xstart
    real(rkind) :: yfocus, ytau, yh, ystart
    real(rkind) :: zfocus, ztau, zh, zstart
    real(rkind), allocatable, dimension(:,:) :: metric_params

    namelist /PROBINPUT/ ns, Lx, Ly, p_inf, rho_inf, M_inf, gam, Rgas, R_v 
    namelist /METRICS/ xmetric_flag, ymetric_flag, zmetric_flag, curvil_flag, metric_params

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
    associate( x  => mesh(:,:,:,1),    y   => mesh(:,:,:,2),    z    => mesh(:,:,:,3), &
               xi => meshcvl(:,:,:,1), eta => meshcvl(:,:,:,2), zeta => meshcvl(:,:,:,3) )
        if (nrank == 0) then
            print *, "Domain size: ", Lx, Ly
        end if

        dxi   = Lx/real(nx-1,rkind)   !periodic
        deta  = Ly/real(ny-1,rkind)   !periodic
        dzeta = dxi  !2D 

        x1 = -Lx/2._rkind; xn = Lx/2._rkind
        y1 = -Ly/2._rkind; yn = Ly/2._rkind
        z1 = zero;         zn = zero


        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    xi(i,j,k)   = x1 + real( ix1 -1 + i - 1, rkind ) * dxi
                    eta(i,j,k)  = y1 + real( iy1 -1 + j - 1, rkind ) * deta
                    zeta(i,j,k) = z1 + real( iz1 -1 + k - 1, rkind ) * dzeta
                end do
            end do
        end do
  
       xfocus = metric_params(1,1);  xtau   = metric_params(1,2);  xstart = metric_params(1,3); xh = metric_params(1,4)
       yfocus = metric_params(2,1);  ytau   = metric_params(2,2);  ystart = metric_params(2,3); yh = metric_params(2,4)
       if(nrank == 0) then
          print*, '>>xfocus=',xfocus, '>>xtau=',xtau, '>>xstart=',xstart, '>>xh=',xh, '>>xflag=',xmetric_flag
          print*, '>>yfocus=',yfocus, '>>ytau=',ytau, '>>ystart=',ystart, '>>yh=',yh, '>>yflag=',ymetric_flag
       endif
        !z = zeta
        !call stretched_coordinates(decomp,x,xi,xmetric_flag,metric_params(1,1),&
        !                            metric_params(1,2),metric_params(1,3),metric_params(1,4))
        !call stretched_coordinates(decomp,y,eta,ymetric_flag,metric_params(2,1),&
        !                            metric_params(2,2),metric_params(2,3),metric_params(2,4))
        do k=1,size(mesh,3)
            do j=1,size(mesh,2)
                do i=1,size(mesh,1)
                    x(i,j,k)   = x1 + (real( ix1 -1 + i - 1, rkind ) + Ampx*sin(num*pi*real(iy1 -1 + j - 1, rkind)*deta/Ly + real(iy1 -1 + i, rkind)*phi/(decomp%xsz(1)-1)) )* dxi
                    y(i,j,k)   = y1 + (real( iy1 -1 + j - 1, rkind ) + Ampy*sin(num*pi*real( ix1 -1 + i - 1, rkind )*dxi/Lx + real(iy1 -1 + j, rkind)*phi/(decomp%ysz(2)-1)) )* deta
                    z(i,j,k)   = z1 + real( iz1 -1 + k - 1, rkind ) * dzeta
                end do
            end do
        end do

    ! Grid width on stretched/uniform mesh
    call transpose_y_to_x(x,xtmp1,decomp)                   ! Decomposition in x
    do i = 1, nx-1
       xtmp2(i,:,:) =  xtmp1(i+1,:,:) - xtmp1(i,:,:)
    end do
    !xtmp2(nx,:,:) = xn - xtmp1(nx,:,:)
    xtmp2(nx,:,:) = xtmp2(nx-1,:,:)
    call transpose_x_to_y(xtmp2,dxs,decomp)   

    do j = 1, decomp%ysz(2)-1
       dys(:,j,:) =  y(:,j+1,:) - y(:,j,:)
    end do
    !dys(:,decomp%ysz(2),:) = yn - y(:,decomp%ysz(2),:)      ! Base decomposition in Y
    dys(:,decomp%ysz(2),:) = dys(:,decomp%ysz(2)-1,:)     ! Base decomposition in Y

    dzs = dzeta
    !call transpose_y_to_z(z,ztmp1,decomp)                   ! Decomposition in z
    !do k = 1, nz-1
    !   ztmp2(:,:,k) =  ztmp1(:,:,k+1) - ztmp1(:,:,k)
    !end do
    !ztmp2(:,:,nz) = zn - ztmp1(:,:,nz)
    !call transpose_z_to_y(ztmp2,dzs,decomp) 

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
    deallocate(metric_params)
end subroutine


subroutine initfields(decomp,dxi,deta,dzeta,inputfile,mesh,fields,mix,tsim,tstop,dt,tviz)
    use kind_parameters,  only: rkind
    use constants,        only: zero,half,one,two,pi,eight
    use CurvilCompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index
    use decomp_2d,        only: decomp_info
    use MixtureEOSMod,    only: mixture
    use Vortex_advection_cvl_data
    

    implicit none
    character(len=*),                intent(in)    :: inputfile
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(inout) :: mix
    real(rkind),                     intent(in)    :: dxi,deta,dzeta
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind),                     intent(inout) :: tsim,tstop,dt,tviz

    real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)) :: tmp
    real(rkind) :: xc = 0.0_rkind, yc = 0.0_rkind, rad_sq, C_v, u_inf    
    integer     :: i,j,k,ioUnit

    namelist /PROBINPUT/ ns, Lx, Ly, p_inf, rho_inf, M_inf, gam, Rgas, R_v 
    ioUnit = 15
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=PROBINPUT)
    close(ioUnit)
    
    associate( rho => fields(:,:,:,rho_index), u => fields(:,:,:,u_index), &
                 v => fields(:,:,:,  v_index), w => fields(:,:,:,w_index), &
                 p => fields(:,:,:,  p_index), T => fields(:,:,:,T_index), &
                 e => fields(:,:,:,  e_index),                             &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        u_inf = M_inf*sqrt(gam) 
        C_v = 0.02*R_v*u_inf
        rho = rho_inf
        do k=1,decomp%ysz(3)
            do j=1,decomp%ysz(2)
                do i=1,decomp%ysz(1)
                   rad_sq   = ( (x(i,j,k)-xc)**two + (y(i,j,k)-yc)**two ) / (R_v**two)
                   u(i,j,k) = u_inf - C_v*(y(i,j,k)-yc)*exp(- rad_sq/two) / (R_v**two)
                   v(i,j,k) = C_v*(x(i,j,k)-xc)*exp(- rad_sq/two) / (R_v**two)
                   p(i,j,k) = p_inf - C_v**two * rho(i,j,k) *exp(- rad_sq) / (two*R_v**two)
                end do
            end do
        end do

    end associate
end subroutine

subroutine hook_output(decomp,der,dxi,deta,dzeta,outputdir,mesh,fields,mix,tsim,vizcount)
    use kind_parameters,  only: rkind,clen
    use constants,        only: zero,half,one,two,pi,eight
    use CurvilCompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use decomp_2d,        only: decomp_info
    use DerivativesMod,   only: derivatives
    use MixtureEOSMod,    only: mixture
    use reductions,       only: P_MEAN
    use Vortex_advection_cvl_data

    implicit none
    character(len=*),                intent(in) :: outputdir
    type(decomp_info),               intent(in) :: decomp
    type(derivatives),               intent(in) :: der
    type(mixture),                   intent(in) :: mix
    real(rkind),                     intent(in) :: dxi,deta,dzeta,tsim
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

    end associate
end subroutine

subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc,newTimeStep, time_step)
    use kind_parameters,  only: rkind
    use decomp_2d,        only: decomp_info, nrank, transpose_y_to_x, transpose_x_to_y
    use constants,        only: zero, half, one, two, three, four, five, six, seven, eight
    use CurvilCompressibleGrid, only: rho_index,u_index,v_index,w_index,p_index,T_index,e_index,mu_index,bulk_index,kap_index,Ys_index
    use MixtureEOSMod,    only: mixture
    use operators,        only: filter3D

    use Vortex_advection_cvl_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc
    logical,                         intent(in)    :: newTimeStep
    integer,                         intent(in)    :: time_step 

    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )


        !!!!! =============  Add Sponge+bulk for exit bc ==========!!!!!
        ! Gradually apply the exit boundary conditions
        ! Apply sponge in X-direction on right
        !call  sponge_x(decomp, mygfil, x, Lx, u, v, w, p, rho, x_bc, y_bc, z_bc)

        ! Apply sponge in Y-direction on top and bottom
        !call  sponge_y(decomp, mygfil, y, Ly, u, v, w, p, rho, x_bc, y_bc, z_bc)

        ! Apply sponge in Z-direction on front and back
        !call  sponge_z(decomp, mygfil, z, Lz, u, v, w, p, rho, x_bc, y_bc, z_bc)
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

    use Vortex_advection_cvl_data

    implicit none
    type(decomp_info),               intent(in) :: decomp
    type(mixture),                   intent(in) :: mix
    integer,                         intent(in) :: step
    real(rkind),                     intent(in) :: tsim
    real(rkind), dimension(:,:,:,:), intent(in) :: mesh
    real(rkind), dimension(:,:,:,:), intent(in) :: fields
    type(sgs_cgrid), optional,       intent(in) :: sgsmodel


    associate( rho    => fields(:,:,:, rho_index), u   => fields(:,:,:,  u_index), &
                 v    => fields(:,:,:,   v_index), w   => fields(:,:,:,  w_index), &
                 p    => fields(:,:,:,   p_index), T   => fields(:,:,:,  T_index), &
                 e    => fields(:,:,:,   e_index), mu  => fields(:,:,:, mu_index), &
                 bulk => fields(:,:,:,bulk_index), kap => fields(:,:,:,kap_index), &
                 Ys   => fields(:,:,:,Ys_index:Ys_index+mix%ns-1),                 &
                 diff => fields(:,:,:,Ys_index+mix%ns:Ys_index+2*mix%ns-1),        &
                 x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )
        call message(2,"Maximum u-velocity",P_MAXVAL(u))
        call message(2,"Minimum u-velocity",P_MINVAL(u))
        call message(2,"Maximum v-velocity",P_MAXVAL(v))
        call message(2,"Minimum v-velocity",P_MINVAL(v))
        call message(2,"Maximum w-velocity",P_MAXVAL(w))
        call message(2,"Minimum w-velocity",P_MINVAL(w))
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

    use Vortex_advection_cvl_data

    implicit none
    type(decomp_info),               intent(in)    :: decomp
    type(mixture),                   intent(in)    :: mix
    real(rkind),                     intent(in)    :: tsim
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(in)    :: fields
    real(rkind), dimension(:,:,:,:), intent(inout) :: rhs

end subroutine
