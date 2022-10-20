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
        module procedure crossprod_components, crossprod_arrays
    end interface

contains

    subroutine gradient(decomp, der, f, dfdx, dfdy, dfdz, x_bc_, y_bc_, z_bc_)
        type(decomp_info), intent(in) :: decomp
        type(derivatives), intent(in) :: der
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(in)  :: f
        real(rkind), dimension(size(f,1), size(f,2), size(f,3)),             intent(out) :: dfdx
        real(rkind), dimension(size(f,1), size(f,2), size(f,3)),             intent(out) :: dfdy
        real(rkind), dimension(size(f,1), size(f,2), size(f,3)),             intent(out) :: dfdz
        integer, dimension(2), optional, intent(in) :: x_bc_, y_bc_, z_bc_
        integer, dimension(2) :: x_bc, y_bc, z_bc

        real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xtmp,xdum
        real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: ztmp,zdum
        
        if ( (size(f,1) .NE. decomp%ysz(1)) .AND. (size(f,2) .NE. decomp%ysz(2)) .AND. (size(f,3) .NE. decomp%ysz(3)) ) then
            call GracefulExit(&
              "Either size of input array to gradient operator is "//&
              "inconsistent with decomp or not in Y decomp. Other "//&
              "decomps have yet to be implemented",234)
        end if

        x_bc = 0; if (present(x_bc_)) x_bc = x_bc_
        y_bc = 0; if (present(y_bc_)) y_bc = y_bc_
        z_bc = 0; if (present(z_bc_)) z_bc = z_bc_

        ! Get Y derivatives
        call der%ddy(f,dfdy,y_bc(1),y_bc(2))

        ! Get X derivatives
        call transpose_y_to_x(f,xtmp,decomp)
        call der%ddx(xtmp,xdum,x_bc(1),x_bc(2))
        call transpose_x_to_y(xdum,dfdx)

        ! Get Z derivatives
        call transpose_y_to_z(f,ztmp,decomp)
        call der%ddz(ztmp,zdum,z_bc(1),z_bc(2))
        call transpose_z_to_y(zdum,dfdz)

    end subroutine 

    subroutine curl(decomp, der, u, v, w, curlu, x_bc_, y_bc_, z_bc_)
        type(decomp_info), intent(in) :: decomp
        type(derivatives), intent(in) :: der
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(in)  :: u, v, w
        real(rkind), dimension(size(u,1), size(u,2), size(u,3),3),           intent(out) :: curlu
        integer, dimension(2), optional, intent(in) :: x_bc_, y_bc_, z_bc_
        integer, dimension(2) :: x_bc, y_bc, z_bc

        real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xtmp,xdum
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: ytmp
        real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: ztmp,zdum
        
        if (    (size(u,1) .NE. decomp%ysz(1)) .AND. (size(u,2) .NE. decomp%ysz(2)) .AND. (size(u,3) .NE. decomp%ysz(3)) ) then
            call GracefulExit(&
              "Either size of input array 'u' to the curl operator is "//&
              "inconsistent with decomp or not in Y decomp. Other "//&
              "decomps have yet to be implemented",234)
        end if
        if (    (size(v,1) .NE. decomp%ysz(1)) .AND. (size(v,2) .NE. decomp%ysz(2)) .AND. (size(v,3) .NE. decomp%ysz(3)) ) then
            call GracefulExit(&
              "Either size of input array 'v' to the curl operator is "//&
              "inconsistent with decomp or not in Y decomp. Other&
                             & decomps have yet to be implemented",234)
        end if
        if (    (size(w,1) .NE. decomp%ysz(1)) .AND. (size(w,2) .NE. decomp%ysz(2)) .AND. (size(w,3) .NE. decomp%ysz(3)) ) then
            call GracefulExit(&
              "Either size of input array 'w' to the curl operator is "//&
              "inconsistent with decomp or not in Y decomp. Other&
                             & decomps have yet to be implemented",234)
        end if

        x_bc = 0; if (present(x_bc_)) x_bc = x_bc_
        y_bc = 0; if (present(y_bc_)) y_bc = y_bc_
        z_bc = 0; if (present(z_bc_)) z_bc = z_bc_

        ! Get dw/dy
        call der%ddy( w, curlu(:,:,:,1), y_bc(1), y_bc(2) )

        ! Get dv/dz
        call transpose_y_to_z( v, ztmp, decomp)
        call der%ddz(ztmp,zdum,z_bc(1),z_bc(2))
        call transpose_z_to_y(zdum,ytmp)

        curlu(:,:,:,1) = curlu(:,:,:,1) - ytmp ! dw/dy - dv/dz

        ! Get du/dz
        call transpose_y_to_z( u, ztmp, decomp)
        call der%ddz(ztmp,zdum,z_bc(1),z_bc(2))
        call transpose_z_to_y(zdum, curlu(:,:,:,2) )

        ! Get dw/dx
        call transpose_y_to_x( w, xtmp, decomp)
        call der%ddx(xtmp,xdum,x_bc(1),x_bc(2))
        call transpose_x_to_y(xdum, ytmp )

        curlu(:,:,:,2) = curlu(:,:,:,2) - ytmp ! du/dz - dw/dx

        ! Get dv/dx
        call transpose_y_to_x( v, xtmp, decomp)
        call der%ddx(xtmp,xdum,x_bc(1),x_bc(2))
        call transpose_x_to_y(xdum, curlu(:,:,:,3) )

        ! Get du/dy
        call der%ddy( u, ytmp, y_bc(1), y_bc(2) )

        curlu(:,:,:,3) = curlu(:,:,:,3) - ytmp ! dv/dx - du/dy

    end subroutine 

    subroutine divergence(decomp, der, u, v, w, div, x_bc_, y_bc_, z_bc_)
        type(decomp_info), intent(in) :: decomp
        type(derivatives) :: der
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), intent(in)  :: u, v, w
        real(rkind), dimension(size(u,1), size(u,2), size(u,3)),             intent(out) :: div
        integer, dimension(2), optional, intent(in) :: x_bc_, y_bc_, z_bc_
        integer, dimension(2) :: x_bc, y_bc, z_bc

        real(rkind), dimension(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)) :: xtmp,xdum
        real(rkind), dimension(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)) :: ytmp
        real(rkind), dimension(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)) :: ztmp,zdum
        
        if ( (size(u,1) .NE. decomp%ysz(1)) .AND. (size(u,2) .NE. decomp%ysz(2)) .AND. (size(u,3) .NE. decomp%ysz(3)) ) then
            call GracefulExit(&
              "Either size of input array to divergence operator is "//&
              "inconsistent with decomp or not in Y decomp. Other&
                             & decomps have yet to be implemented",234)
        end if

        x_bc = 0; if (present(x_bc_)) x_bc = x_bc_
        y_bc = 0; if (present(y_bc_)) y_bc = y_bc_
        z_bc = 0; if (present(z_bc_)) z_bc = z_bc_

        ! Get Y derivatives
        call der%ddy(v,div,-y_bc(1),-y_bc(2))

        ! Get X derivatives
        call transpose_y_to_x(u,xtmp,decomp)
        call der%ddx(xtmp,xdum,-x_bc(1),-x_bc(2))
        call transpose_x_to_y(xdum,ytmp)

        div = div + ytmp

        ! Get Z derivatives
        call transpose_y_to_z(w,ztmp,decomp)
        call der%ddz(ztmp,zdum,-z_bc(1),-z_bc(2))
        call transpose_z_to_y(zdum,ytmp)

        div = div + ytmp

    end subroutine

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

        if ( (size(arr,1) .NE. decomp%ysz(1)) .AND. (size(arr,2) .NE. decomp%ysz(2)) .AND. (size(arr,3) .NE. decomp%ysz(3)) ) then
            call GracefulExit(&
              "Either size of input array to divergence operator is "//&
              "inconsistent with decomp or not in Y decomp. Other&
                             & decomps have yet to be implemented",234)
        end if

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

end module
