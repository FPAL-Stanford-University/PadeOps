module hitCD10_RHS

    use kind_parameters, only: rkind
    use DerivativesMod,  only: derivatives
    implicit none

contains

    subroutine calc_xRHS(ux,vx,wx,xRHS,Re,xder,xwrk1)
        class( derivatives ),            intent(in)    :: xder
        real(rkind), dimension(:,:,:),   intent(in)    :: ux,vx,wx
        real(rkind), dimension(:,:,:,:), intent(inout) :: xRHS
        real(rkind),                     intent(in)    :: Re
        real(rkind), dimension(:,:,:),   intent(inout) :: xwrk1
       
        ! Get all the X derivatives
        associate( RHSx => xRHS(:,:,:,1), RHSy => xRHS(:,:,:,2), RHSz => xRHS(:,:,:,3) ) 
    
            call xder%ddx(ux*ux,xwrk1)
            RHSx = RHSx - xwrk1
            call xder%d2dx2(ux,xwrk1)
            RHSx = RHSx + xwrk1/Re
    
            call xder%ddx(ux*vx,xwrk1)
            RHSy = RHSy - xwrk1
            call xder%d2dx2(vx,xwrk1)
            RHSy = RHSy + xwrk1/Re
    
            call xder%ddx(ux*wx,xwrk1)
            RHSz = RHSz - xwrk1
            call xder%d2dx2(wx,xwrk1)
            RHSz = RHSz + xwrk1/Re
    
        end associate
    end subroutine
    
    subroutine calc_yRHS(uy,vy,wy,yRHS,Re,yder,ywrk1)
        class( derivatives ),            intent(in)    :: yder
        real(rkind), dimension(:,:,:),   intent(in)    :: uy,vy,wy
        real(rkind), dimension(:,:,:,:), intent(inout) :: yRHS
        real(rkind),                     intent(in)    :: Re
        real(rkind), dimension(:,:,:),   intent(inout) :: ywrk1
        
        ! Get all the Y derivatives
        associate( RHSx => yRHS(:,:,:,1), RHSy => yRHS(:,:,:,2), RHSz => yRHS(:,:,:,3) ) 
    
            call yder%ddy(uy*vy,ywrk1)
            RHSx = RHSx - ywrk1
            call yder%d2dy2(uy,ywrk1)
            RHSx = RHSx + ywrk1/Re
    
            call yder%ddy(vy*vy,ywrk1)
            RHSy = RHSy - ywrk1 
            call yder%d2dy2(vy,ywrk1)
            RHSy = RHSy + ywrk1/Re
    
            call yder%ddy(vy*wy,ywrk1)
            RHSz = RHSz - ywrk1
            call yder%d2dy2(wy,ywrk1)
            RHSz = RHSz + ywrk1/Re
    
        end associate
    end subroutine
    
    subroutine calc_zRHS(uz,vz,wz,zRHS,Re,zder,zwrk1)
        class( derivatives ),            intent(in)    :: zder
        real(rkind), dimension(:,:,:),   intent(in)    :: uz,vz,wz
        real(rkind), dimension(:,:,:,:), intent(inout) :: zRHS
        real(rkind),                     intent(in)    :: Re
        real(rkind), dimension(:,:,:),   intent(inout) :: zwrk1
        
        ! Get all the Z derivatives
        associate( RHSx => zRHS(:,:,:,1), RHSy => zRHS(:,:,:,2), RHSz => zRHS(:,:,:,3) ) 
    
            call zder%ddz(uz*wz,zwrk1)
            RHSx = RHSx - zwrk1
            call zder%d2dz2(uz,zwrk1)
            RHSx = RHSx + zwrk1/Re
    
            call zder%ddz(wz*vz,zwrk1)
            RHSy = RHSy - zwrk1 
            call zder%d2dz2(vz,zwrk1)
            RHSy = RHSy + zwrk1/Re
    
            call zder%ddz(wz*wz,zwrk1)
            RHSz = RHSz - zwrk1
            call zder%d2dz2(wz,zwrk1)
            RHSz = RHSz + zwrk1/Re
    
        end associate
    end subroutine

end module
