module hitCD10_timestep

    use kind_parameters, only: rkind
    use constants,       only: zero
    use DerivativesMod,  only: derivatives
    use decomp_2d,       only: decomp_info,transpose_x_to_y, transpose_y_to_x, &
                                           transpose_y_to_z, transpose_z_to_y
    use hitCD10_RHS,     only: calc_xRHS, calc_yRHS, calc_zRHS

    implicit none

contains

    subroutine timestep_x_to_z(ux,vx,wx,uy,vy,wy,uz,vz,wz,xRHS,yRHS,zRHS,Re,xder,yder,zder,gp)
        real(rkind), dimension(:,:,:),   intent(inout) :: ux,vx,wx
        real(rkind), dimension(:,:,:),   intent(inout) :: uy,vy,wy
        real(rkind), dimension(:,:,:),   intent(inout) :: uz,vz,wz
        real(rkind), dimension(:,:,:,:), intent(inout) :: xRHS,yRHS,zRHS
        real(rkind),                     intent(in)    :: Re
        class( decomp_info ),            intent(in)    :: gp
        class( derivatives ),            intent(in)    :: xder,yder,zder
        
        real(rkind), dimension( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) :: xwrk1
        real(rkind), dimension( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) :: ywrk1
        real(rkind), dimension( gp%zsz(1), gp%zsz(2), gp%zsz(3) ) :: zwrk1
        
        integer :: i

        xRHS = zero
    
        ! Get X derivatives 
        call calc_xRHS(ux,vx,wx,xRHS,Re,xder,xwrk1)
    
        ! Transpose data to Y decomp
        call transpose_x_to_y(ux,uy,gp)
        call transpose_x_to_y(vx,vy,gp)
        call transpose_x_to_y(wx,wy,gp)
        do i=1,3
            call transpose_x_to_y(xRHS(:,:,:,i),yRHS(:,:,:,i),gp)
        end do
    
        ! Get Y derivatives 
        call calc_yRHS(uy,vy,wy,yRHS,Re,yder,ywrk1)
    
        ! Transpose data to Z decomp
        call transpose_y_to_z(uy,uz,gp)
        call transpose_y_to_z(vy,vz,gp)
        call transpose_y_to_z(wy,wz,gp)
        do i=1,3
            call transpose_y_to_z(yRHS(:,:,:,i),zRHS(:,:,:,i),gp)
        end do
    
        ! Get Z derivatives 
        call calc_zRHS(uz,vz,wz,zRHS,Re,zder,zwrk1)
    
    end subroutine
    
    subroutine timestep_z_to_x(ux,vx,wx,uy,vy,wy,uz,vz,wz,xRHS,yRHS,zRHS,Re,xder,yder,zder,gp)
        real(rkind), dimension(:,:,:),   intent(inout) :: ux,vx,wx
        real(rkind), dimension(:,:,:),   intent(inout) :: uy,vy,wy
        real(rkind), dimension(:,:,:),   intent(inout) :: uz,vz,wz
        real(rkind), dimension(:,:,:,:), intent(inout) :: xRHS,yRHS,zRHS
        real(rkind),                     intent(in)    :: Re
        class( decomp_info ),            intent(in)    :: gp
        class( derivatives ),            intent(in)    :: xder,yder,zder
        
        real(rkind), dimension( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) :: xwrk1
        real(rkind), dimension( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) :: ywrk1
        real(rkind), dimension( gp%zsz(1), gp%zsz(2), gp%zsz(3) ) :: zwrk1
        
        integer :: i

        zRHS = zero
    
        ! Get Z derivatives 
        call calc_zRHS(uz,vz,wz,zRHS,Re,zder,zwrk1)
    
        ! Transpose data to Y decomp
        call transpose_z_to_y(uz,uy,gp)
        call transpose_z_to_y(vz,vy,gp)
        call transpose_z_to_y(wz,wy,gp)
        do i=1,3
            call transpose_z_to_y(zRHS(:,:,:,i),yRHS(:,:,:,i),gp)
        end do
    
        ! Get Y derivatives 
        call calc_yRHS(uy,vy,wy,yRHS,Re,yder,ywrk1)
    
        ! Transpose data to X decomp
        call transpose_y_to_x(uy,ux,gp)
        call transpose_y_to_x(vy,vx,gp)
        call transpose_y_to_x(wy,wx,gp)
        do i=1,3
            call transpose_y_to_x(yRHS(:,:,:,i),xRHS(:,:,:,i),gp)
        end do
    
        ! Get X derivatives 
        call calc_xRHS(ux,vx,wx,xRHS,Re,xder,xwrk1)
    
    end subroutine

    subroutine get_divergence_x(ux,vx,wx,vy,wy,wz,xder,yder,zder,gp,divz)
        real(rkind), dimension(:,:,:),   intent(in)    :: ux,vx,wx
        real(rkind), dimension(:,:,:),   intent(inout) :: vy,wy
        real(rkind), dimension(:,:,:),   intent(inout) :: wz
        class( decomp_info ),            intent(in)    :: gp
        class( derivatives ),            intent(in)    :: xder,yder,zder
        real(rkind), dimension( gp%zsz(1), gp%zsz(2), gp%zsz(3) ), intent(inout) :: divz
        
        real(rkind), dimension( gp%xsz(1), gp%xsz(2), gp%xsz(3) ) :: divx
        real(rkind), dimension( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) :: divy,ywrk
        real(rkind), dimension( gp%zsz(1), gp%zsz(2), gp%zsz(3) ) :: zwrk
        
        ! Get X derivatives 
        call xder%ddx(ux,divx)
   
        ! Transpose data to Y decomp
        call transpose_x_to_y(vx,vy,gp)
        call transpose_x_to_y(wx,wy,gp)
        call transpose_x_to_y(divx,divy,gp)
   
        ! Get Y derivatives 
        call yder%ddy(vy,ywrk)
        divy = divy + ywrk
    
        ! Transpose data to Z decomp
        call transpose_y_to_z(wy,wz,gp)
        call transpose_y_to_z(divy,divz,gp)
    
        ! Get Z derivatives 
        call zder%ddz(wz,zwrk)
        divz = divz + zwrk
        
    end subroutine

end module
