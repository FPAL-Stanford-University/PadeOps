subroutine init_mgm(this)
   use constants, only: two, eight
   class(sgs_cgrid), intent(inout) :: this
   real(rkind)  :: deltaLES, ceps_2 = 1.0_rkind, cepsT_2 = 1.0_rkind 

   this%isEddyViscosityModel = .false. 
   this%isEddyDiffModel = .false. 
   
   this%cmgm_x = this%dx**2 / 12.0d0
   this%cmgm_y = this%dy**2 / 12.0d0 
   this%cmgm_z = this%dz**2 / 12.0d0
   this%PrCpfac= this%Cp / this%Pr
   this%c1_mgm = eight * (this%deltaLES**two) * ceps_2 
   this%c2_mgm = two * (this%deltaLES**two) * cepsT_2

   call message(1,"MGM model initialized") 

end subroutine

subroutine destroy_mgm(this)

   class(sgs_cgrid), intent(inout) :: this
   this%isEddyViscosityModel = .false. 

end subroutine

subroutine get_tausgs_mgm(rho, duidxj, gradT, tausgs, Qsgs)
   use contants, only: eps,two, zero
   real(rkind), dimension(this%nxL,this%nyL,this%nzL  ), intent(in)  :: rho
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,9), intent(in)  :: duidxj
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,6), intent(out) :: tausgs

   real(rkind)  :: G11, G12, G13, G22, G23, G33 , Gmm , ckl, multfactor
   integer :: i, j, k

   do k = 1, this%nzL
      do j = 1, this%nyL 
         do i = 1, this%nxL
            G11 = this%cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,1) + this%cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,2) + &
                  this%cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,3)

            G12 = this%cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,4) + this%cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,5) + &
                  this%cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,6)

            G13 = this%cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,7) + this%cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,8) + &
                  this%cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,9)

            G22 = this%cmgm_x*duidxj(i,j,k,4)*duidxj(i,j,k,1) + this%cmgm_y*duidxj(i,j,k,5)*duidxj(i,j,k,2) + &
                  this%cmgm_z*duidxj(i,j,k,6)*duidxj(i,j,k,3)

            G23 = this%cmgm_x*duidxj(i,j,k,4)*duidxj(i,j,k,4) + this%cmgm_y*duidxj(i,j,k,5)*duidxj(i,j,k,5) + &
                  this%cmgm_z*duidxj(i,j,k,6)*duidxj(i,j,k,6)
            G33 = this%cmgm_x*duidxj(i,j,k,4)*duidxj(i,j,k,7) + this%cmgm_y*duidxj(i,j,k,5)*duidxj(i,j,k,8) + &
                  this%cmgm_z*duidxj(i,j,k,6)*duidxj(i,j,k,9)

            Gmm = G11 + G22 + G33 + eps
            
            ckl = - (  G11 * this%S_ij(i,j,k,1) + G22 * this%S_ij(i,j,k,4) + G33 * this%S_ij(i,j,k,6) + &
                two * (G12 * this%S_ij(i,j,k,2) + G13 * this%S_ij(i,j,k,3) + G23 * this%S_ij(i,j,k,5) ) ) / Gmm
           
            multfactor = this%c1_mgm * ckl * ckl * rho(i,j,k) / Gmm
 
            tausgs(i,j,k,1) = multfactor * G11
            tausgs(i,j,k,2) = multfactor * G12
            tausgs(i,j,k,3) = multfactor * G13
            tausgs(i,j,k,4) = multfactor * G22
            tausgs(i,j,k,5) = multfactor * G23
            tausgs(i,j,k,6) = multfactor * G33
 
            if (ckl .GT. zero) then
            ! tausgs(i,j,k,:)=tausgs(i,j,k,:)
            else
            tausgs(i,j,k,:) = zero
            end if

         end do
      end do
   end do
end subroutine


subroutine get_Qjsgs_mgm(rho, duidxj, gradT, Qjsgs)
   use contants, only: eps,two, zero
   real(rkind), dimension(this%nxL,this%nyL,this%nzL  ), intent(in)  :: rho
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,9), intent(in)  :: duidxj
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,3), intent(in)  :: gradT
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,3), intent(out) :: Qjsgs

   real(rkind)  :: G11, G12, G13, G22, G23, G33, Gmm, ckl, cnT, G1T, G2T, GT3, modG, multfactor
   integer :: i, j, k

   do k = 1, this%nzL
      do j = 1, this%nyL 
         do i = 1, this%nxL
            G11 = this%cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,1) + this%cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,2) + &
                  this%cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,3)

            G12 = this%cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,4) + this%cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,5) + &
                  this%cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,6)

            G13 = this%cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,7) + this%cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,8) + &
                  this%cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,9)

            G22 = this%cmgm_x*duidxj(i,j,k,4)*duidxj(i,j,k,1) + this%cmgm_y*duidxj(i,j,k,5)*duidxj(i,j,k,2) + &
                  this%cmgm_z*duidxj(i,j,k,6)*duidxj(i,j,k,3)

            G23 = this%cmgm_x*duidxj(i,j,k,4)*duidxj(i,j,k,4) + this%cmgm_y*duidxj(i,j,k,5)*duidxj(i,j,k,5) + &
                  this%cmgm_z*duidxj(i,j,k,6)*duidxj(i,j,k,6)

            G33 = this%cmgm_x*duidxj(i,j,k,4)*duidxj(i,j,k,7) + this%cmgm_y*duidxj(i,j,k,5)*duidxj(i,j,k,8) + &
                  this%cmgm_z*duidxj(i,j,k,6)*duidxj(i,j,k,9)

            Gmm = G11 + G22 + G33 + eps

            G1T = this%cmgm_x*duidxj(i,j,k,1)*gradT(i,j,k,1) + this%cmgm_y*duidxj(i,j,k,2)*gradT(i,j,k,2) + &
                  this%cmgm_z*duidxj(i,j,k,3)*gradT(i,j,k,3)

            G2T = this%cmgm_x*duidxj(i,j,k,4)*gradT(i,j,k,1) + this%cmgm_y*duidxj(i,j,k,5)*gradT(i,j,k,2) + &
                  this%cmgm_z*duidxj(i,j,k,6)*gradT(i,j,k,3)

            G3T = this%cmgm_x*duidxj(i,j,k,7)*gradT(i,j,k,1) + this%cmgm_y*duidxj(i,j,k,8)*gradT(i,j,k,2) + &
                  this%cmgm_z*duidxj(i,j,k,9)*gradT(i,j,k,3)

            modG = sqrt(G1T**2 + G2T**2 + G3T**2) + eps
            
            ckl = -( G11 * this%S_ij(i,j,k,1) + G22 * this%S_ij(i,j,k,4) + G33 * this%S_ij(i,j,k,6) + &
                two*(G12 * this%S_ij(i,j,k,2) + G13 * this%S_ij(i,j,k,3) + G23 * this%S_ij(i,j,k,5) ) ) / Gmm

            cnT = -( G1T * gradT(i,j,k,1) + G2T * gradT(i,j,k,2) + G3T * gradT(i,j,k,3) )/ modG
            
            multfactor = this%c2_mgm * ckl * cnT * rho(i,j,k) * this%PrCpfac / modG
            
            Qjsgs(i,j,k,1) = multfactor * G1T
            Qjsgs(i,j,k,2) = multfactor * G2T
            Qjsgs(i,j,k,3) = multfactor * G3T
 
            if ((ckl .GT. zero) .AND. (cnT .GT. zero)) then
                ! Qjsgs(i,j,k,:)=Qjsgs(i,j,k,:)
            else
                Qjsgs(i,j,k,:) = zero
            end if

         end do
      end do
   end do
end subroutine
