subroutine init_mgm(this, dx, dy, dz)
   use constants, only: two, eight
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), intent(in) :: dx, dy, dz
   real(rkind)  :: deltaLES, ceps_2 = 1.0_rkind, cepsT_2 = 1.0_rkind 

   this%isEddyViscosityModel = .false. 
   
   if (.not. this%isPeriodic) then
      deltaLES = (1.5d0*dx*1.5d0*dy*dz)**(1.d0/3.d0)
   else
      deltaLES =  (1.5d0*dx*1.5d0*dy*1.5d0*dz)**(1.d0/3.d0)
   end if 

   this%cmgm_x = dx**two / 12.0d0
   this%cmgm_y = dy**two / 12.0d0 
   this%cmgm_z = dz**two / 12.0d0
   this%PrCpfac= this%Cp / this%Pr
   this%c1     = eight * (deltaLES**two) * ceps_2 
   this%c2     = two * (deltaLES**two) * cepsT_2
   call message(1,"MGM model initialized") 
end subroutine

subroutine destroy_mgm(this)
   class(sgs_cgrid), intent(inout) :: this
   this%isEddyViscosityModel = .false. 
end subroutine

subroutine get_tausgs_mgm(cmgm_x, cmgm_y, cmgm_z, c1, rho, duidxj, Sij, nxL, nyL, nzL, tausgs)
   use contants, only: eps,two, zero
   integer, intent(in) :: nxL, nyL, nzL
   real(rkind), intent(in) :: camd_x, camd_y, camd_z, c1
   real(rkind), dimension(:,:,:,:), intent(out) :: tausgs
   real(rkind), dimension(:,:,:,:), intent(in)  :: Sij , duidxj
   real(rkind), dimension(:,:,:) , intent(in)   :: rho
   real(rkind)  :: G11, G12, G13, G22, G23, G33 , Gmm , ckl
   integer :: i,j,k

   do k = 1,nzL
      do j = 1,nyL 
         do i = 1,nxL
            G11=cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,1)+cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,2)+cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,3)
            G12=cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,4)+cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,5)+cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,6)
            G13=cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,7)+cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,8)+cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,9)
            G22=cmgm_x*duidxj(i,j,k,4)*duidxj(i,j,k,1)+cmgm_y*duidxj(i,j,k,5)*duidxj(i,j,k,2)+cmgm_z*duidxj(i,j,k,6)*duidxj(i,j,k,3)
            G23=cmgm_x*duidxj(i,j,k,4)*duidxj(i,j,k,4)+cmgm_y*duidxj(i,j,k,5)*duidxj(i,j,k,5)+cmgm_z*duidxj(i,j,k,6)*duidxj(i,j,k,6)
            G33=cmgm_x*duidxj(i,j,k,4)*duidxj(i,j,k,7)+cmgm_y*duidxj(i,j,k,5)*duidxj(i,j,k,8)+cmgm_z*duidxj(i,j,k,6)*duidxj(i,j,k,9)
            Gmm=G11+G22+G33+eps
            
            ckl=-(G11*Sij(i,j,k,1)+G22*Sij(i,j,k,4)+G33*Sij(i,j,k,6)+2*(G12*Sij(i,j,k,2)+G13*Sij(i,j,k,3)+G23*Sij(i,j,k,5)))/Gmm
            
            tausgs(i,j,k,1) = c1 * (ckl**two) * (G11/Gmm) * rho(i,j,k)
            tausgs(i,j,k,2) = c1 * (ckl**two) * (G12/Gmm) * rho(i,j,k)
            tausgs(i,j,k,3) = c1 * (ckl**two) * (G13/Gmm) * rho(i,j,k)
            tausgs(i,j,k,4) = c1 * (ckl**two) * (G22/Gmm) * rho(i,j,k)
            tausgs(i,j,k,5) = c1 * (ckl**two) * (G23/Gmm) * rho(i,j,k)
            tausgs(i,j,k,6) = c1 * (ckl**two) * (G33/Gmm) * rho(i,j,k)
 
            if (ckl .GT. zero) then
            ! tausgs(i,j,k,:)=tausgs(i,j,k,:)
            else
            tausgs(i,j,k,:) = zero
            end if

         end do
      end do
   end do
end subroutine


subroutine get_Qjsgs_mgm(cmgm_x, cmgm_y, cmgm_z, c2, PrCpfac, rho, duidxj,gradT, Sij, nxL, nyL, nzL, tausgs, Qsgs)
   use contants, only: eps,two, zero
   integer, intent(in) :: nxL, nyL, nzL
   real(rkind), intent(in) :: camd_x, camd_y, camd_z, c2, PrCpfac
   real(rkind), dimension(:,:,:,:), intent(out) :: Qsgs
   real(rkind), dimension(:,:,:,:), intent(in)  :: Sij , gradT,duidxj, tausgs
   real(rkind), dimension(:,:,:) , intent(in)   :: rho
   real(rkind)  :: G11, G12, G13, G22, G23, G33 , Gmm , ckl, G1T, G2T, GT3, modG
   integer :: i,j,k

   do k = 1,nzL
      do j = 1,nyL 
         do i = 1,nxL
            G11=cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,1)+cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,2)+cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,3)
            G12=cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,4)+cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,5)+cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,6)
            G13=cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,7)+cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,8)+cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,9)
            G22=cmgm_x*duidxj(i,j,k,4)*duidxj(i,j,k,1)+cmgm_y*duidxj(i,j,k,5)*duidxj(i,j,k,2)+cmgm_z*duidxj(i,j,k,6)*duidxj(i,j,k,3)
            G23=cmgm_x*duidxj(i,j,k,4)*duidxj(i,j,k,4)+cmgm_y*duidxj(i,j,k,5)*duidxj(i,j,k,5)+cmgm_z*duidxj(i,j,k,6)*duidxj(i,j,k,6)
            G33=cmgm_x*duidxj(i,j,k,4)*duidxj(i,j,k,7)+cmgm_y*duidxj(i,j,k,5)*duidxj(i,j,k,8)+cmgm_z*duidxj(i,j,k,6)*duidxj(i,j,k,9)
            Gmm=G11+G22+G33+eps
            G1T=cmgm_x*duidxj(i,j,k,1)*gradT(i,j,k,1)+cmgm_y*duidxj(i,j,k,2)*gradT(i,j,k,2)+cmgm_z*duidxj(i,j,k,3)*gradT(i,j,k,3)
            G2T=cmgm_x*duidxj(i,j,k,4)*gradT(i,j,k,1)+cmgm_y*duidxj(i,j,k,5)*gradT(i,j,k,2)+cmgm_z*duidxj(i,j,k,6)*gradT(i,j,k,3)
            G3T=cmgm_x*duidxj(i,j,k,7)*gradT(i,j,k,1)+cmgm_y*duidxj(i,j,k,8)*gradT(i,j,k,2)+cmgm_z*duidxj(i,j,k,9)*gradT(i,j,k,3)
            modG=sqrt(G1T**two + G2T**two + G3T**two)          
            
            ckl=-(G11*Sij(i,j,k,1)+G22*Sij(i,j,k,4)+G33*Sij(i,j,k,6)+2*(G12*Sij(i,j,k,2)+G13*Sij(i,j,k,3)+G23*Sij(i,j,k,5)))/Gmm
            
            Qsgs(i,j,k,1) = c2 * (ckl**two) * (G11/Gmm) * rho(i,j,k)
            Qsgs(i,j,k,2) = c2 * (ckl**two) * (G12/Gmm) * rho(i,j,k)
            tausgs(i,j,k,3) = c2 * (ckl**two) * (G13/Gmm) * rho(i,j,k)
 
            if (ckl .GT. zero) then
            ! tausgs(i,j,k,:)=tausgs(i,j,k,:)
            else
            tausgs(i,j,k,:) = zero
            end if

         end do
      end do
   end do
