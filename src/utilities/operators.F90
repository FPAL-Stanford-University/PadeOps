module operators
    use kind_parameters, only: rkind, clen
    use FiltersMod, only: filters
    use gridtools, only: alloc_buffs, destroy_buffs
    use decomp_2d, only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                    transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
    use DerivativesMod,  only: derivatives
    use exits,           only: GracefulExit
   
    implicit none

    interface crossprod
        module procedure crossprod_components, crossprod_arrays, crossprod_vec_tens
    end interface

    interface dotprod
        module procedure dotprod_components, dotprod_symmtens_vec, dotprod_tens_vec
    end interface

    interface curl
        module procedure curl_vec, curl_tens
    end interface

    interface gradient
        module procedure gradient_vec, gradient_scalar
    end interface

    interface divergence
        module procedure divergence_vec, divergence_symmtens, divergence_tens
    end interface

contains

    subroutine gradient_scalar(decomp, der, x, y, z, f, dfdx, dfdy, dfdz, coordsys, x_bc_, y_bc_, z_bc_)
        type(decomp_info), intent(in) :: decomp
        type(derivatives), intent(in) :: der
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(in)  :: x,y,z,f
        real(rkind), dimension(size(f,1), size(f,2), size(f,3)),             intent(out) :: dfdx
        real(rkind), dimension(size(f,1), size(f,2), size(f,3)),             intent(out) :: dfdy
        real(rkind), dimension(size(f,1), size(f,2), size(f,3)),             intent(out) :: dfdz
        integer,                         intent(in) :: coordsys
        integer, dimension(2), optional, intent(in) :: x_bc_, y_bc_, z_bc_
        integer, dimension(2) :: x_bc, y_bc, z_bc

        !real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xtmp,xdum
        !real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: ztmp,zdum
        
        call check_size(decomp,dfdx,41);   call check_size(decomp,dfdy,42);   call check_size(decomp,dfdz,43)
        call check_size(decomp,f,   44);

        x_bc = 0; if (present(x_bc_)) x_bc = x_bc_
        y_bc = 0; if (present(y_bc_)) y_bc = y_bc_
        z_bc = 0; if (present(z_bc_)) z_bc = z_bc_

        call der_x(decomp,der,f,dfdx,x_bc(1),x_bc(2))
        call der_y(decomp,der,f,dfdy,y_bc(1),y_bc(2))
        call der_z(decomp,der,f,dfdz,z_bc(1),z_bc(2))

        select case(coordsys)
        case (1) ! Cylindrical
            dfdy = dfdy/x
        case (2) ! Spherical
            call GracefulExit("Gradient of a scalar in spherical coordinates not implemented yet", 234)
        end select

    end subroutine 

    subroutine gradient_vec(decomp,der,x,y,z,fx,fy,fz,dfxdx,dfxdy,dfxdz,dfydx,dfydy,dfydz,dfzdx,dfzdy,dfzdz,coordsys,x_bc_,y_bc_,z_bc_)
        type(decomp_info), intent(in) :: decomp
        type(derivatives), intent(in) :: der
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(in)   :: x,y,z,fx,fy,fz
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(out)  :: dfxdx,dfxdy,dfxdz
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(out)  :: dfydx,dfydy,dfydz
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(out)  :: dfzdx,dfzdy,dfzdz
        integer, dimension(2), optional, intent(in) :: x_bc_, y_bc_, z_bc_
        integer,               intent(in) :: coordsys

        integer, dimension(2) :: x_bc, y_bc, z_bc

        !real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xtmp,xdum
        !real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: ztmp,zdum
       
        call check_size(decomp,fx,    1);   call check_size(decomp,fy,    2);   call check_size(decomp,fz,    3)
        call check_size(decomp,dfxdx, 4);   call check_size(decomp,dfxdy, 5);   call check_size(decomp,dfxdz, 6)
        call check_size(decomp,dfydx, 7);   call check_size(decomp,dfydy, 8);   call check_size(decomp,dfydz, 9)
        call check_size(decomp,dfzdx,10);   call check_size(decomp,dfzdy,11);   call check_size(decomp,dfzdz,12)

        x_bc = 0; if (present(x_bc_)) x_bc = x_bc_
        y_bc = 0; if (present(y_bc_)) y_bc = y_bc_
        z_bc = 0; if (present(z_bc_)) z_bc = z_bc_

        ! Get Y derivatives
        call der_y(decomp,der,fx,dfxdy,y_bc(1),y_bc(2))
        call der_y(decomp,der,fy,dfydy,y_bc(1),y_bc(2))
        call der_y(decomp,der,fz,dfzdy,y_bc(1),y_bc(2))

        ! Get X derivatives
        call der_x(decomp,der,fx,dfxdx,x_bc(1),x_bc(2))
        call der_x(decomp,der,fy,dfydx,x_bc(1),x_bc(2))
        call der_x(decomp,der,fz,dfzdx,x_bc(1),x_bc(2))

        ! Get Z derivatives
        call der_z(decomp,der,fx,dfxdz,z_bc(1),z_bc(2))
        call der_z(decomp,der,fy,dfydz,z_bc(1),z_bc(2))
        call der_z(decomp,der,fz,dfzdz,z_bc(1),z_bc(2))

        select case(coordsys)
        case(1) ! Cylindrical
            dfxdy = (dfxdy - fy)/x
            dfydy = (dfydy + fx)/x
            dfzdy =  dfzdy/x
        case(2) ! Spherical
            call GracefulExit("Gradient of a vector in spherical coordinates not implemented yet", 234)
        end select

    end subroutine 

    subroutine curl_vec(decomp, der, x, y, z, u, v, w, curlu, coordsys, x_bc_, y_bc_, z_bc_)
        type(decomp_info), intent(in) :: decomp
        type(derivatives), intent(in) :: der
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(in)  :: x, y, z, u, v, w
        real(rkind), dimension(size(u,1), size(u,2), size(u,3),3),           intent(out) :: curlu
        integer,                         intent(in) :: coordsys
        integer, dimension(2), optional, intent(in) :: x_bc_, y_bc_, z_bc_
        integer, dimension(2) :: x_bc, y_bc, z_bc

        !real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xtmp,xdum
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: ytmp
        !real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: ztmp,zdum
        
        call check_size(decomp,u,             13);   call check_size(decomp,v,              14);   
        call check_size(decomp,w,             15);   call check_size(decomp,curlu(:,:,:,1), 16);
        call check_size(decomp,curlu(:,:,:,2),17);   call check_size(decomp,curlu(:,:,:,3), 18)

        x_bc = 0; if (present(x_bc_)) x_bc = x_bc_
        y_bc = 0; if (present(y_bc_)) y_bc = y_bc_
        z_bc = 0; if (present(z_bc_)) z_bc = z_bc_

        select case(coordsys)
        case (0) ! Cartesian
            ! --------- Step 1 :: Get dw/dy - dvdz ----------
            call der_y(decomp, der, w, curlu(:,:,:,1), y_bc(1), y_bc(2))
            call der_z(decomp, der, v, ytmp,           z_bc(1), z_bc(2))
            curlu(:,:,:,1) = curlu(:,:,:,1) - ytmp

            ! ---------- Step 2 :: Get du/dz - dwdx ----------
            call der_z(decomp, der, u, curlu(:,:,:,2), z_bc(1), z_bc(2))
            call der_x(decomp, der, w, ytmp,           x_bc(1), x_bc(2))
            curlu(:,:,:,2) = curlu(:,:,:,2) - ytmp

            ! ---------- Step 3 :: Get dv/dx - dudy ----------
            call der_x(decomp, der, v, curlu(:,:,:,3), x_bc(1), x_bc(2))
            call der_y(decomp, der, u, ytmp,           y_bc(1), y_bc(2))
            curlu(:,:,:,3) = curlu(:,:,:,3) - ytmp
        case (1) ! Cylindrical
            ! --------- Step 1 :: Get dw/dy - dvdz ----------
            call der_y(decomp, der, w, curlu(:,:,:,1), y_bc(1), y_bc(2))
            call der_z(decomp, der, v, ytmp,           z_bc(1), z_bc(2))
            curlu(:,:,:,1) = curlu(:,:,:,1)/x - ytmp

            ! ---------- Step 2 :: Get du/dz - dwdx ----------
            call der_z(decomp, der, u, curlu(:,:,:,2), z_bc(1), z_bc(2))
            call der_x(decomp, der, w, ytmp,           x_bc(1), x_bc(2))
            curlu(:,:,:,2) = curlu(:,:,:,2) - ytmp

            ! ---------- Step 3 :: Get dv/dx - dudy ----------
            call der_x(decomp, der, v, curlu(:,:,:,3), x_bc(1), x_bc(2))
            call der_y(decomp, der, u, ytmp,           y_bc(1), y_bc(2))
            curlu(:,:,:,3) = curlu(:,:,:,3) - ytmp/x + v/x
        case (2) ! Spherical
            call GracefulExit("Curl of a vector in spherical coordinates not implemented yet", 234)
        end select

    end subroutine 

    subroutine curl_tens(decomp,der,x,y,z,fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz,ttxx,ttxy,ttxz,ttyx,ttyy,ttyz,ttzx,ttzy,ttzz,coordsys,x_bc_,y_bc_,z_bc_)
        type(decomp_info), intent(in) :: decomp
        type(derivatives), intent(in) :: der
        real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), intent(in)  :: x,y,z,fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz
        real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), intent(out) :: ttxx,ttxy,ttxz,ttyx,ttyy,ttyz,ttzx,ttzy,ttzz
        integer, dimension(2), optional, intent(in) :: x_bc_, y_bc_, z_bc_
        integer,               intent(in) :: coordsys

        integer, dimension(2) :: x_bc, y_bc, z_bc
        !real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xtmp,xdum
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: ytmp
        !real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: ztmp,zdum
        
        call check_size(decomp,fxx,19);   call check_size(decomp,fxy,20);   call check_size(decomp,fxz,21);
        call check_size(decomp,fyx,22);   call check_size(decomp,fyy,23);   call check_size(decomp,fyz,24);
        call check_size(decomp,fzx,25);   call check_size(decomp,fzy,26);   call check_size(decomp,fzz,27);
        call check_size(decomp,ttxx,28);   call check_size(decomp,ttxy,29);   call check_size(decomp,ttxz,30);
        call check_size(decomp,ttyx,31);   call check_size(decomp,ttyy,32);   call check_size(decomp,ttyz,33);
        call check_size(decomp,ttzx,34);   call check_size(decomp,ttzy,35);   call check_size(decomp,ttzz,36);

        x_bc = 0; if (present(x_bc_)) x_bc = x_bc_
        y_bc = 0; if (present(y_bc_)) y_bc = y_bc_
        z_bc = 0; if (present(z_bc_)) z_bc = z_bc_

        select case(coordsys)
        case(0) ! Cartesian
            !-------------- Step 1 :: ttxx --------------
            call der_y(decomp, der, fzx, ttxx, y_bc(1), y_bc(2))
            call der_z(decomp, der, fyx, ytmp, z_bc(1), z_bc(2))
            ttxx = ttxx - ytmp ! dfzx/dy - dfyx/dz

            !-------------- Step 2 :: ttxy --------------
            call der_y(decomp, der, fzy, ttxy, -y_bc(1), -y_bc(2))
            call der_z(decomp, der, fyy, ytmp,  z_bc(1),  z_bc(2))
            ttxy = ttxy - ytmp ! dfzy/dy - dfyy/dz

            !-------------- Step 3 :: ttxz --------------
            call der_y(decomp, der, fzz, ttxy,  y_bc(1),  y_bc(2))
            call der_z(decomp, der, fyz, ytmp, -z_bc(1), -z_bc(2))
            ttxz = ttxz - ytmp ! dfzz/dy - dfyz/dz

            !-------------- Step 4 :: ttyx --------------
            call der_z(decomp, der, fxx, ttyx,  z_bc(1),  z_bc(2))
            call der_x(decomp, der, fzx, ytmp, -x_bc(1), -x_bc(2))
            ttyx = ttyx - ytmp ! dfxx/dz - dfzx/dx

            !-------------- Step 5 :: ttyy --------------
            call der_z(decomp, der, fxy, ttyy, z_bc(1), z_bc(2))
            call der_x(decomp, der, fzy, ytmp, x_bc(1), x_bc(2))
            ttyy = ttyy - ytmp ! dfxy/dz - dfzy/dx

            !-------------- Step 6 :: ttyz --------------
            call der_z(decomp, der, fxz, ttyz, -z_bc(1), -z_bc(2))
            call der_x(decomp, der, fzz, ytmp,  x_bc(1),  x_bc(2))
            ttyz = ttyz - ytmp ! dfxz/dz - dfzz/dx

            !-------------- Step 7 :: ttzx --------------
            call der_x(decomp, der, fyx, ttzx, -x_bc(1), -x_bc(2))
            call der_y(decomp, der, fxx, ytmp,  y_bc(1),  y_bc(2))
            ttzx = ttzx - ytmp ! dfyx/dx - dfxx/dy

            !-------------- Step 8 :: ttzy --------------
            call der_x(decomp, der, fyy, ttzy,  x_bc(1),  x_bc(2))
            call der_y(decomp, der, fxy, ytmp, -y_bc(1), -y_bc(2))
            ttzy = ttzy - ytmp ! dfyy/dx - dfxy/dy

            !-------------- Step 9 :: ttzz --------------
            call der_x(decomp, der, fyz, ttzz, x_bc(1), x_bc(2))
            call der_y(decomp, der, fxz, ytmp, y_bc(1), y_bc(2))
            ttzz = ttzz - ytmp ! dfyz/dx - dfxz/dy
        case(1) ! Cylindrical
            !-------------- Step 1 :: ttxx --------------
            call der_y(decomp, der, fzx, ttxx, y_bc(1), y_bc(2))
            call der_z(decomp, der, fyx, ytmp, z_bc(1), z_bc(2))
            ttxx = (ttxx - fzy)/x - ytmp ! dfzx/dy - dfyx/dz

            !-------------- Step 2 :: ttxy --------------
            call der_y(decomp, der, fzy, ttxy, -y_bc(1), -y_bc(2))
            call der_z(decomp, der, fyy, ytmp,  z_bc(1),  z_bc(2))
            ttxy = (ttxy + fzx)/x - ytmp ! dfzy/dy - dfyy/dz

            !-------------- Step 3 :: ttxz --------------
            call der_y(decomp, der, fzz, ttxy,  y_bc(1),  y_bc(2))
            call der_z(decomp, der, fyz, ytmp, -z_bc(1), -z_bc(2))
            ttxz = ttxz/x - ytmp ! dfzz/dy - dfyz/dz

            !-------------- Step 4 :: ttyx --------------
            call der_z(decomp, der, fxx, ttyx,  z_bc(1),  z_bc(2))
            call der_x(decomp, der, fzx, ytmp, -x_bc(1), -x_bc(2))
            ttyx = ttyx - ytmp ! dfxx/dz - dfzx/dx

            !-------------- Step 5 :: ttyy --------------
            call der_z(decomp, der, fxy, ttyy, z_bc(1), z_bc(2))
            call der_x(decomp, der, fzy, ytmp, x_bc(1), x_bc(2))
            ttyy = ttyy - ytmp ! dfxy/dz - dfzy/dx

            !-------------- Step 6 :: ttyz --------------
            call der_z(decomp, der, fxz, ttyz, -z_bc(1), -z_bc(2))
            call der_x(decomp, der, fzz, ytmp,  x_bc(1),  x_bc(2))
            ttyz = ttyz - ytmp ! dfxz/dz - dfzz/dx

            !-------------- Step 7 :: ttzx --------------
            call der_x(decomp, der, fyx, ttzx, -x_bc(1), -x_bc(2))
            call der_y(decomp, der, fxx, ytmp,  y_bc(1),  y_bc(2))
            ttzx = ttzx - (ytmp - fxy - fyx)/x ! dfyx/dx - dfxx/dy

            !-------------- Step 8 :: ttzy --------------
            call der_x(decomp, der, fyy, ttzy,  x_bc(1),  x_bc(2))
            call der_y(decomp, der, fxy, ytmp, -y_bc(1), -y_bc(2))
            ttzy = ttzy - (ytmp + fxx - fyy)/x ! dfyy/dx - dfxy/dy

            !-------------- Step 9 :: ttzz --------------
            call der_x(decomp, der, fyz, ttzz, x_bc(1), x_bc(2))
            call der_y(decomp, der, fxz, ytmp, y_bc(1), y_bc(2))
            ttzz = ttzz - (ytmp - fyz)/x ! dfyz/dx - dfxz/dy
        case(2) ! Spherical
            call GracefulExit("Curl of a tensor in spherical coordinates not implemented yet", 234)
        end select

    end subroutine 

    subroutine divergence_vec(decomp, der, x, y, z, u, v, w, div, coordsys, x_bc_, y_bc_, z_bc_)
        type(decomp_info), intent(in) :: decomp
        type(derivatives), intent(in) :: der
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(in)  :: x, y, z, u, v, w
        real(rkind), dimension(size(u,1), size(u,2), size(u,3)),             intent(out) :: div
        integer, dimension(2), optional, intent(in) :: x_bc_, y_bc_, z_bc_
        integer,                         intent(in) :: coordsys
        integer, dimension(2) :: x_bc, y_bc, z_bc

        !real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xtmp,xdum
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: ytmp1,ytmp2
        !real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: ztmp,zdum
        
        call check_size(decomp,u,37);   call check_size(decomp,v,38);   call check_size(decomp,w,39);
        call check_size(decomp,div,40)

        x_bc = 0; if (present(x_bc_)) x_bc = x_bc_
        y_bc = 0; if (present(y_bc_)) y_bc = y_bc_
        z_bc = 0; if (present(z_bc_)) z_bc = z_bc_

        call der_x(decomp, der, u,div  ,x_bc(1),x_bc(2))
        call der_y(decomp, der, v,ytmp1,y_bc(1),y_bc(2))
        call der_z(decomp, der, w,ytmp2,z_bc(1),z_bc(2))

        select case(coordsys)
        case(0) ! Cartesian
            div = div + ytmp1 + ytmp2
        case(1) ! Cylindrical
            div = div + ytmp2 + (ytmp1 + u)/x
        case(2) ! Spherical
            call GracefulExit("Divergence of a vector in spherical coordinates not implemented yet", 234)
        end select

    end subroutine

    subroutine divergence_symmtens(decomp, der, x, y, z, fxx, fxy, fxz, fyy, fyz, fzz, divx, divy, divz, coordsys, x_bc_, y_bc_, z_bc_)
        type(decomp_info), intent(in) :: decomp
        type(derivatives), intent(in) :: der
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(in)  :: x,y,z,fxx,fxy,fxz,fyy,fyz,fzz
        real(rkind), dimension(size(x,1), size(x,2), size(x,3)),             intent(out) :: divx,divy,divz
        integer,                         intent(in) :: coordsys
        integer, dimension(2), optional, intent(in) :: x_bc_, y_bc_, z_bc_
        integer, dimension(2) :: x_bc, y_bc, z_bc

        !real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xtmp,xdum
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: ytmp1,ytmp2
        !real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: ztmp,zdum
        
        call check_size(decomp,fxx, 45);   call check_size(decomp,fxy, 46);  call check_size(decomp,fxz, 47);
        call check_size(decomp,fyy, 48);   call check_size(decomp,fyz, 49);  call check_size(decomp,fzz, 50);
        call check_size(decomp,divx,51);   call check_size(decomp,divy,51);  call check_size(decomp,divz,51);

        x_bc = 0; if (present(x_bc_)) x_bc = x_bc_
        y_bc = 0; if (present(y_bc_)) y_bc = y_bc_
        z_bc = 0; if (present(z_bc_)) z_bc = z_bc_

        select case(coordsys)
        case(0) ! Cartesian
            ! -------------------Step 1 :: divx------------------------
            call der_x(decomp, der, fxx, divx,  x_bc(1),x_bc(2))
            call der_y(decomp, der, fxy, ytmp1, y_bc(1),y_bc(2))
            call der_z(decomp, der, fxz, ytmp2, z_bc(1),z_bc(2))
            divx = divx + ytmp1 + ytmp2
            ! -------------------Step 2 :: divy------------------------
            call der_x(decomp, der, fxy, divy,  x_bc(1),x_bc(2))
            call der_y(decomp, der, fyy, ytmp1, y_bc(1),y_bc(2))
            call der_z(decomp, der, fyz, ytmp2, z_bc(1),z_bc(2))
            divy = divy + ytmp1 + ytmp2
            ! -------------------Step 3 :: divz------------------------
            call der_x(decomp, der, fxz, divz,  x_bc(1),x_bc(2))
            call der_y(decomp, der, fyz, ytmp1, y_bc(1),y_bc(2))
            call der_z(decomp, der, fzz, ytmp2, z_bc(1),z_bc(2))
            divz = divz + ytmp1 + ytmp2
            ! ---------------------------------------------------------
        case(1) ! Cylindrical
            ! -------------------Step 1 :: divx------------------------
            call der_x(decomp, der, fxx, divx,  x_bc(1),x_bc(2))
            call der_y(decomp, der, fxy, ytmp1, y_bc(1),y_bc(2))
            call der_z(decomp, der, fxz, ytmp2, z_bc(1),z_bc(2))
            divx = divx + (ytmp1 + fxx - fyy)/x + ytmp2
            ! -------------------Step 2 :: divy------------------------
            call der_x(decomp, der, fxy, divy,  x_bc(1),x_bc(2))
            call der_y(decomp, der, fyy, ytmp1, y_bc(1),y_bc(2))
            call der_z(decomp, der, fyz, ytmp2, z_bc(1),z_bc(2))
            divy = divy + (ytmp1 + fxy + fxy)/x + ytmp2
            ! -------------------Step 3 :: divz------------------------
            call der_x(decomp, der, fxz, divz,  x_bc(1),x_bc(2))
            call der_y(decomp, der, fyz, ytmp1, y_bc(1),y_bc(2))
            call der_z(decomp, der, fzz, ytmp2, z_bc(1),z_bc(2))
            divz = divz + (ytmp1 + fxz)/x + ytmp2
            ! ---------------------------------------------------------
        case(2) ! Spherical
            call GracefulExit("Divergence of a tensor in spherical coordinates not implemented yet", 234)
        end select

    end subroutine

    subroutine divergence_tens(decomp, der, x, y, z, fxx, fxy, fxz, fyx, fyy, fyz, fzx, fzy, fzz, divx, divy, divz, coordsys, x_bc_, y_bc_, z_bc_)
        type(decomp_info), intent(in) :: decomp
        type(derivatives), intent(in) :: der
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(in)  :: x,y,z,fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz
        real(rkind), dimension(size(x,1), size(x,2), size(x,3)),             intent(out) :: divx,divy,divz
        integer,                         intent(in) :: coordsys
        integer, dimension(2), optional, intent(in) :: x_bc_, y_bc_, z_bc_
        integer, dimension(2) :: x_bc, y_bc, z_bc

        !real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xtmp,xdum
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: ytmp1,ytmp2
        !real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: ztmp,zdum
        
        call check_size(decomp,fxx, 45);   call check_size(decomp,fxy, 46);  call check_size(decomp,fxz, 47);
        call check_size(decomp,fyy, 48);   call check_size(decomp,fyz, 49);  call check_size(decomp,fzz, 50);
        call check_size(decomp,divx,51);   call check_size(decomp,divy,51);  call check_size(decomp,divz,51);

        x_bc = 0; if (present(x_bc_)) x_bc = x_bc_
        y_bc = 0; if (present(y_bc_)) y_bc = y_bc_
        z_bc = 0; if (present(z_bc_)) z_bc = z_bc_


        select case(coordsys)
        case(0) ! Cartesian
            ! -------------------Step 1 :: divx------------------------
            call der_x(decomp, der, fxx, divx,  x_bc(1),x_bc(2))
            call der_y(decomp, der, fyx, ytmp1, y_bc(1),y_bc(2))
            call der_z(decomp, der, fzx, ytmp2, z_bc(1),z_bc(2))
            divx = divx + ytmp1 + ytmp2
            ! -------------------Step 2 :: divy------------------------
            call der_x(decomp, der, fxy, divy,  x_bc(1),x_bc(2))
            call der_y(decomp, der, fyy, ytmp1, y_bc(1),y_bc(2))
            call der_z(decomp, der, fzy, ytmp2, z_bc(1),z_bc(2))
            divy = divy + ytmp1 + ytmp2
            ! -------------------Step 3 :: divz------------------------
            call der_x(decomp, der, fzx, divz,  x_bc(1),x_bc(2))
            call der_y(decomp, der, fzy, ytmp1, y_bc(1),y_bc(2))
            call der_z(decomp, der, fzz, ytmp2, z_bc(1),z_bc(2))
            divz = divz + ytmp1 + ytmp2
            ! ---------------------------------------------------------
        case(1) ! Cylindrical
            ! -------------------Step 1 :: divx------------------------
            call der_x(decomp, der, fxx, divx,  x_bc(1),x_bc(2))
            call der_y(decomp, der, fyx, ytmp1, y_bc(1),y_bc(2))
            call der_z(decomp, der, fzx, ytmp2, z_bc(1),z_bc(2))
            divx = divx + (ytmp1 + fxx - fyy)/x + ytmp2
            ! -------------------Step 2 :: divy------------------------
            call der_x(decomp, der, fxy, divy,  x_bc(1),x_bc(2))
            call der_y(decomp, der, fyy, ytmp1, y_bc(1),y_bc(2))
            call der_z(decomp, der, fzy, ytmp2, z_bc(1),z_bc(2))
            divy = divy + (ytmp1 + fxy + fyx)/x + ytmp2
            ! -------------------Step 3 :: divz------------------------
            call der_x(decomp, der, fxz, divz,  x_bc(1),x_bc(2))
            call der_y(decomp, der, fyz, ytmp1, y_bc(1),y_bc(2))
            call der_z(decomp, der, fzz, ytmp2, z_bc(1),z_bc(2))
            divz = divz + (ytmp1 + fxz)/x + ytmp2
            ! ---------------------------------------------------------
        case(2) ! Spherical
            call GracefulExit("Divergence of a tensor in spherical coordinates not implemented yet", 234)
        end select

    end subroutine

    subroutine crossprod_vec_tens(ux,uy,uz,fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz,ttxx,ttxy,ttxz,ttyx,ttyy,ttyz,ttzx,ttzy,ttzz)
        real(kind=rkind), dimension(:,:,:), intent(in) :: ux,uy,uz,fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz
        real(kind=rkind), dimension(:,:,:), intent(out) :: ttxx,ttxy,ttxz,ttyx,ttyy,ttyz,ttzx,ttzy,ttzz

        ttxx = uy*fzx - uz*fyx;        ttxy = uy*fzy - uz*fyy;        ttxz = uy*fzz - uz*fyz;
        ttyx = uz*fxx - ux*fzx;        ttyy = uz*fxy - ux*fzy;        ttyz = uz*fxz - ux*fzz;
        ttzx = ux*fyx - uy*fxx;        ttzy = ux*fyy - uy*fxy;        ttzz = ux*fyz - uy*fxz;
    
    end subroutine crossprod_vec_tens
    
    subroutine dotprod_components(ux,uy,uz,vx,vy,vz,r)
        real(kind=rkind), dimension(:,:,:), intent(in) :: ux,uy,uz,vx,vy,vz
        real(kind=rkind), dimension(:,:,:), intent(out) :: r
    
        r = ux*ux + uy*uy + uz*uz
    
    end subroutine dotprod_components

    subroutine dotprod_tens_vec(fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz,ux,uy,uz,rx,ry,rz)
        real(kind=rkind), dimension(:,:,:), intent(in)  :: fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz,ux,uy,uz
        real(kind=rkind), dimension(:,:,:), intent(out) :: rx,ry,rz

        rx = fxx*ux + fxy*uy + fxz*uz 
        ry = fyx*ux + fyy*uy + fyz*uz 
        rx = fzx*ux + fzy*uy + fzz*uz 

    end subroutine dotprod_tens_vec
    
    subroutine dotprod_symmtens_vec(fxx,fxy,fxz,fyy,fyz,fzz,ux,uy,uz,rx,ry,rz)
        real(kind=rkind), dimension(:,:,:), intent(in)  :: fxx,fxy,fxz,fyy,fyz,fzz,ux,uy,uz
        real(kind=rkind), dimension(:,:,:), intent(out) :: rx,ry,rz

        rx = fxx*ux + fxy*uy + fxz*uz 
        ry = fxy*ux + fyy*uy + fyz*uz 
        rx = fxz*ux + fyz*uy + fzz*uz 

    end subroutine dotprod_symmtens_vec

    subroutine crossprod_components(rx,ry,rz,ux,uy,uz,vx,vy,vz)
        real(kind=rkind), dimension(:,:,:), intent(in) :: ux,uy,uz,vx,vy,vz
        real(kind=rkind), dimension(:,:,:), intent(out) :: rx,ry,rz
    
        rx = uy*vz - uz*vy
        ry = uz*vx - ux*vz
        rz = ux*vy - uy*vx
    
    end subroutine crossprod_components
    
    subroutine crossprod_arrays(r,u,v)
        real(kind=rkind), dimension(:,:,:,:), intent(in) :: u,v
        real(kind=rkind), dimension(:,:,:,:), intent(out) :: r
    
        r(:,:,:,1) = u(:,:,:,2)*v(:,:,:,3) - u(:,:,:,3)*v(:,:,:,2)
        r(:,:,:,2) = u(:,:,:,3)*v(:,:,:,1) - u(:,:,:,1)*v(:,:,:,3)
        r(:,:,:,3) = u(:,:,:,1)*v(:,:,:,2) - u(:,:,:,2)*v(:,:,:,1)
    
    end subroutine crossprod_arrays

    subroutine filter3D(decomp,fil,arr,numtimes, x_bc_, y_bc_, z_bc_)
        type(decomp_info), intent(in) :: decomp
        type(filters),     intent(in) :: fil
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(inout)  :: arr
        integer, optional, intent(in) :: numtimes
        integer, dimension(2), optional, intent(in) :: x_bc_, y_bc_, z_bc_
        integer, dimension(2) :: x_bc, y_bc, z_bc

        integer :: times2fil, idx
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: tmp_in_y
        real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: tmp1_in_x, tmp2_in_x
        real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: tmp1_in_z, tmp2_in_z


        if (present(numtimes)) then
            times2fil = numtimes
        else
            times2fil = 1
        end if

        x_bc = 0; if (present(x_bc_)) x_bc = x_bc_
        y_bc = 0; if (present(y_bc_)) y_bc = y_bc_
        z_bc = 0; if (present(z_bc_)) z_bc = z_bc_

        ! First filter in y
        call fil%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
        ! Subsequent refilters 
        do idx = 1,times2fil-1
            arr = tmp_in_y
            call fil%filtery(arr,tmp_in_y,y_bc(1),y_bc(2))
        end do

        ! Then transpose to x
        call transpose_y_to_x(tmp_in_y,tmp1_in_x,decomp)

        ! First filter in x
        call fil%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
        ! Subsequent refilters
        do idx = 1,times2fil-1
            tmp1_in_x = tmp2_in_x
            call fil%filterx(tmp1_in_x,tmp2_in_x,x_bc(1),x_bc(2))
        end do

        ! Now transpose back to y
        call transpose_x_to_y(tmp2_in_x,tmp_in_y,decomp)

        ! Now transpose to z
        call transpose_y_to_z(tmp_in_y,tmp1_in_z,decomp)

        !First filter in z
        call fil%filterz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
        ! Subsequent refilters
        do idx = 1,times2fil-1
            tmp1_in_z = tmp2_in_z
            call fil%filterz(tmp1_in_z,tmp2_in_z,z_bc(1),z_bc(2))
        end do

        ! Now transpose back to y
        call transpose_z_to_y(tmp2_in_z,arr,decomp)

        ! Finished
    end subroutine

    subroutine der_x(decomp, der, fin, fout, bc1, bc2)
        type(decomp_info), intent(in) :: decomp
        type(derivatives), intent(in) :: der
        real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), intent(in)  :: fin
        real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), intent(out) :: fout
        integer,           intent(in) :: bc1, bc2

        real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xtmp,xdum

        call transpose_y_to_x( fin, xtmp, decomp)
        call der%ddx(xtmp,xdum, bc1, bc2)
        call transpose_x_to_y(xdum, fout, decomp )

    end subroutine der_x

    subroutine der_y(decomp, der, fin, fout, bc1, bc2)
        type(decomp_info), intent(in) :: decomp
        type(derivatives), intent(in) :: der
        real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), intent(in)  :: fin
        real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), intent(out) :: fout
        integer,           intent(in) :: bc1, bc2

        call der%ddy( fin, fout, bc1, bc2 )

    end subroutine der_y

    subroutine der_z(decomp, der, fin, fout, bc1, bc2)
        type(decomp_info), intent(in) :: decomp
        type(derivatives), intent(in) :: der
        real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), intent(in)  :: fin
        real(rkind), dimension(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), intent(out) :: fout
        integer,           intent(in) :: bc1, bc2

        real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: ztmp,zdum

        call transpose_y_to_z( fin, ztmp, decomp)
        call der%ddz( ztmp, zdum, bc1, bc2 )
        call transpose_z_to_y(zdum,fout, decomp)

    end subroutine der_z

    subroutine check_size(decomp, f, flag)
        type(decomp_info), intent(in) :: decomp
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(in)   :: f
        integer, intent(in) :: flag
        character(len=clen) :: msg

        write(msg,'(a,i4)') "Size of array inconsistent with Y decomp in operators. Check entry number ", flag
        if ( (size(f,1) .NE. decomp%ysz(1)) .AND. (size(f,2) .NE. decomp%ysz(2)) .AND. (size(f,3) .NE. decomp%ysz(3)) ) then
            call GracefulExit(trim(msg), 234)
        end if
    end subroutine check_size

end module
