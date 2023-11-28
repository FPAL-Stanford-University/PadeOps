subroutine init_mgm(this)
   class(sgs_cgrid), intent(inout) :: this
  ! real(rkind)  :: deltaLES, ceps_2 = 1.0_rkind, cepsT_2 = 1.0_rkind 
   integer :: i, j, k

   this%isEddyViscosityModel = .false. 
   this%isEddyDiffModel = .false. 
   
   do k = 1,this%nzL
      do j = 1,this%nyL
         do i = 1,this%nxL
            this%cmgm_x(i) = this%dxs(i)**2 / 12.0d0
            this%cmgm_y(j) = this%dys(j)**2 / 12.0d0 
            this%cmgm_z(k) = this%dzs(k)**2 / 12.0d0
            this%c1_mgm(i,j,k) = eight * (this%deltaLES(i,j,k)**2) !* ceps_2 
            this%c2_mgm(i,j,k) = two * (this%deltaLES(i,j,k)**2) !* cepsT_2
         end do
      end do
   end do

   this%PrCpfac= this%Cp / this%Pr
   this%cmodel_global =  this%Csgs
   this%cmodel_global_Qjsgs =  this%Csgs/this%Prsgs

   call message(1,"MGM model initialized") 

end subroutine

subroutine destroy_mgm(this)

   class(sgs_cgrid), intent(inout) :: this
   this%isEddyViscosityModel = .false. 

end subroutine

!subroutine get_tausgs_mgm(this, rho, duidxj, tausgs)
subroutine get_mgm_kernel(this, rho, duidxj, Sij, gradT, tausgs, Qjsgs)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL,this%nyL,this%nzL  ), intent(in)  :: rho
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,9), intent(in)  :: duidxj
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,6), intent(in)  :: Sij
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,3), intent(in)  :: gradT
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,6), intent(out) :: tausgs
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,3), intent(out) :: Qjsgs

   real(rkind)  :: G11, G12, G13, G22, G23, G33 , Gmm , ckl, multfactor, maxrho
   real(rkind)  :: G1T, G2T, G3T, cnT, modG
   integer :: i, j, k

   do k = 1, this%nzL
      do j = 1, this%nyL 
         do i = 1, this%nxL
            !! For tausgs
            G11 = this%cmgm_x(i)*duidxj(i,j,k,1)*duidxj(i,j,k,1) + this%cmgm_y(j)*duidxj(i,j,k,2)*duidxj(i,j,k,2) + &
                  this%cmgm_z(k)*duidxj(i,j,k,3)*duidxj(i,j,k,3)

            G12 = this%cmgm_x(i)*duidxj(i,j,k,1)*duidxj(i,j,k,4) + this%cmgm_y(j)*duidxj(i,j,k,2)*duidxj(i,j,k,5) + &
                  this%cmgm_z(k)*duidxj(i,j,k,3)*duidxj(i,j,k,6)

            G13 = this%cmgm_x(i)*duidxj(i,j,k,1)*duidxj(i,j,k,7) + this%cmgm_y(j)*duidxj(i,j,k,2)*duidxj(i,j,k,8) + &
                  this%cmgm_z(k)*duidxj(i,j,k,3)*duidxj(i,j,k,9)

            G22 = this%cmgm_x(i)*duidxj(i,j,k,4)*duidxj(i,j,k,4) + this%cmgm_y(j)*duidxj(i,j,k,5)*duidxj(i,j,k,5) + &
                  this%cmgm_z(k)*duidxj(i,j,k,6)*duidxj(i,j,k,6)

            G23 = this%cmgm_x(i)*duidxj(i,j,k,4)*duidxj(i,j,k,7) + this%cmgm_y(j)*duidxj(i,j,k,5)*duidxj(i,j,k,8) + &
                  this%cmgm_z(k)*duidxj(i,j,k,6)*duidxj(i,j,k,9)

            G33 = this%cmgm_x(i)*duidxj(i,j,k,7)*duidxj(i,j,k,7) + this%cmgm_y(j)*duidxj(i,j,k,8)*duidxj(i,j,k,8) + &
                  this%cmgm_z(k)*duidxj(i,j,k,9)*duidxj(i,j,k,9)

            Gmm = G11 + G22 + G33 + eps
            
            ckl = - (  G11 * Sij(i,j,k,1) + G22 * Sij(i,j,k,4) + G33 * Sij(i,j,k,6) + &
                two * (G12 * Sij(i,j,k,2) + G13 * Sij(i,j,k,3) + G23 * Sij(i,j,k,5) ) ) / Gmm
           
            multfactor = this%c1_mgm(i,j,k) * ckl * ckl * rho(i,j,k) / Gmm
 
            if (ckl .GE. zero) then
                ! tausgs(i,j,k,:)=tausgs(i,j,k,:)
                tausgs(i,j,k,1) = multfactor * G11
                tausgs(i,j,k,2) = multfactor * G12
                tausgs(i,j,k,3) = multfactor * G13
                tausgs(i,j,k,4) = multfactor * G22
                tausgs(i,j,k,5) = multfactor * G23
                tausgs(i,j,k,6) = multfactor * G33
            else
                tausgs(i,j,k,:) = zero
            end if

            !! For Qjsgs
            G1T = this%cmgm_x(i)*duidxj(i,j,k,1)*gradT(i,j,k,1) + this%cmgm_y(j)*duidxj(i,j,k,2)*gradT(i,j,k,2) + &
                  this%cmgm_z(k)*duidxj(i,j,k,3)*gradT(i,j,k,3)

            G2T = this%cmgm_x(i)*duidxj(i,j,k,4)*gradT(i,j,k,1) + this%cmgm_y(j)*duidxj(i,j,k,5)*gradT(i,j,k,2) + &
                  this%cmgm_z(k)*duidxj(i,j,k,6)*gradT(i,j,k,3)

            G3T = this%cmgm_x(i)*duidxj(i,j,k,7)*gradT(i,j,k,1) + this%cmgm_y(j)*duidxj(i,j,k,8)*gradT(i,j,k,2) + &
                  this%cmgm_z(k)*duidxj(i,j,k,9)*gradT(i,j,k,3)

            modG = sqrt(G1T**2 + G2T**2 + G3T**2) + eps
            
            cnT = -( G1T * gradT(i,j,k,1) + G2T * gradT(i,j,k,2) + G3T * gradT(i,j,k,3) )/ modG
            multfactor = this%c2_mgm(i,j,k) * ckl * cnT * rho(i,j,k) * this%PrCpfac / modG

            if ((ckl .GE. zero) .AND. (cnT .GE. zero)) then
                ! Qjsgs(i,j,k,:)=Qjsgs(i,j,k,:)
                Qjsgs(i,j,k,1) = multfactor * G1T
                Qjsgs(i,j,k,2) = multfactor * G2T
                Qjsgs(i,j,k,3) = multfactor * G3T
            else
                Qjsgs(i,j,k,:) = zero
            end if

         end do
      end do
   end do

end subroutine

subroutine multiply_by_model_coefficient_mgm(this, tausgs, Qjsgs)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,6), intent(inout)  :: tausgs
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,3), intent(inout)  :: Qjsgs

   integer  :: ii, j, k

   if( (this%DynamicProcedureType==0) .or. (this%DynamicProcedureType==2)) then
      !! constant coefficient or Global-Dynamic Procedure
      do ii = 1, 6
          tausgs(:,:,:,ii) = this%cmodel_global * tausgs(:,:,:,ii)
      enddo 

      do ii = 1, 3
          Qjsgs(:,:,:,ii) = this%cmodel_global_Qjsgs * Qjsgs(:,:,:,ii)
      enddo 
   elseif(this%DynamicProcedureType==1) then
      !! Local-Dynamic Procedure - averaged in (x,z); function of (y)
      do ii = 1, 6
        do k = 1, this%nzL
          do j = 1, this%nyL
            tausgs(:,j,k,ii) = this%cmodel_local(j) * tausgs(:,j,k,ii)
          end do
        end do
      end do

      do ii = 1, 3
        do k = 1, this%nzL
          do j = 1, this%nyL
            Qjsgs(:,j,k,ii) = this%cmodel_local_Qjsgs(j) * Qjsgs(:,j,k,ii)
          end do
        end do
      end do
   endif

end subroutine

subroutine multiply_by_model_coefficient_mgm_Qjsgs(this, Qjsgs)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,3), intent(inout)  :: Qjsgs

!!   integer  :: ii, j, k
!!
!!   if( (this%DynamicProcedureType==0) .or. (this%DynamicProcedureType==2)) then
!!      !! constant coefficient or Global-Dynamic Procedure
!!      do ii = 1, 3
!!          Qjsgs(:,:,:,ii) = this%cmodel_global_Qjsgs * Qjsgs(:,:,:,ii)
!!      enddo 
!!   elseif(this%DynamicProcedureType==1) then
!!      !! Local-Dynamic Procedure - averaged in (x,z); function of (y)
!!      do ii = 1, 3
!!        do k = 1, this%nzL
!!          do j = 1, this%nyL
!!            Qjsgs(:,j,k,ii) = this%cmodel_local_Qjsgs(j) * Qjsgs(:,j,k,ii)
!!          end do
!!        end do
!!      end do
!!   endif

end subroutine

!subroutine get_Qjsgs_mgm(this, rho, duidxj, gradT, Qjsgs)
subroutine get_Qjsgs_mgm_kernel(this, rho, duidxj, Sij, gradT, Qjsgs)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL,this%nyL,this%nzL  ), intent(in)  :: rho
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,9), intent(in)  :: duidxj
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,6), intent(in)  :: Sij
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,3), intent(in)  :: gradT
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,3), intent(out) :: Qjsgs

!!   real(rkind)  :: G11, G12, G13, G22, G23, G33, Gmm, ckl, cnT, G1T, G2T, G3T, modG, multfactor
!!   integer :: i, j, k
!!
!!   do k = 1, this%nzL
!!      do j = 1, this%nyL 
!!         do i = 1, this%nxL
!!            G11 = this%cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,1) + this%cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,2) + &
!!                  this%cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,3)
!!
!!            G12 = this%cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,4) + this%cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,5) + &
!!                  this%cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,6)
!!
!!            G13 = this%cmgm_x*duidxj(i,j,k,1)*duidxj(i,j,k,7) + this%cmgm_y*duidxj(i,j,k,2)*duidxj(i,j,k,8) + &
!!                  this%cmgm_z*duidxj(i,j,k,3)*duidxj(i,j,k,9)
!!
!!            G22 = this%cmgm_x*duidxj(i,j,k,4)*duidxj(i,j,k,4) + this%cmgm_y*duidxj(i,j,k,5)*duidxj(i,j,k,5) + &
!!                  this%cmgm_z*duidxj(i,j,k,6)*duidxj(i,j,k,6)
!!
!!            G23 = this%cmgm_x*duidxj(i,j,k,4)*duidxj(i,j,k,7) + this%cmgm_y*duidxj(i,j,k,5)*duidxj(i,j,k,8) + &
!!                  this%cmgm_z*duidxj(i,j,k,6)*duidxj(i,j,k,9)
!!
!!            G33 = this%cmgm_x*duidxj(i,j,k,7)*duidxj(i,j,k,7) + this%cmgm_y*duidxj(i,j,k,8)*duidxj(i,j,k,8) + &
!!                  this%cmgm_z*duidxj(i,j,k,9)*duidxj(i,j,k,9)
!!
!!            Gmm = G11 + G22 + G33 + eps
!!
!!            G1T = this%cmgm_x*duidxj(i,j,k,1)*gradT(i,j,k,1) + this%cmgm_y*duidxj(i,j,k,2)*gradT(i,j,k,2) + &
!!                  this%cmgm_z*duidxj(i,j,k,3)*gradT(i,j,k,3)
!!
!!            G2T = this%cmgm_x*duidxj(i,j,k,4)*gradT(i,j,k,1) + this%cmgm_y*duidxj(i,j,k,5)*gradT(i,j,k,2) + &
!!                  this%cmgm_z*duidxj(i,j,k,6)*gradT(i,j,k,3)
!!
!!            G3T = this%cmgm_x*duidxj(i,j,k,7)*gradT(i,j,k,1) + this%cmgm_y*duidxj(i,j,k,8)*gradT(i,j,k,2) + &
!!                  this%cmgm_z*duidxj(i,j,k,9)*gradT(i,j,k,3)
!!
!!            modG = sqrt(G1T**2 + G2T**2 + G3T**2) + eps
!!            
!!            ckl = -( G11 * Sij(i,j,k,1) + G22 * Sij(i,j,k,4) + G33 * Sij(i,j,k,6) + &
!!                two*(G12 * Sij(i,j,k,2) + G13 * Sij(i,j,k,3) + G23 * Sij(i,j,k,5) ) ) / Gmm
!!
!!            cnT = -( G1T * gradT(i,j,k,1) + G2T * gradT(i,j,k,2) + G3T * gradT(i,j,k,3) )/ modG
!!            
!!            multfactor = this%c2_mgm * ckl * cnT * rho(i,j,k) * this%PrCpfac / modG
!!            
!!            Qjsgs(i,j,k,1) = multfactor * G1T
!!            Qjsgs(i,j,k,2) = multfactor * G2T
!!            Qjsgs(i,j,k,3) = multfactor * G3T
!! 
!!            
!!
!!            if ((ckl .GE. zero) .AND. (cnT .GE. zero)) then
!!                ! Qjsgs(i,j,k,:)=Qjsgs(i,j,k,:)
!!            else
!!                Qjsgs(i,j,k,:) = zero
!!            end if
!!
!!         end do
!!      end do
!!   end do
!!   !print*, sum(abs(Qjsgs(:,:,:,1))), sum(abs(Qjsgs(:,:,:,2))), sum(abs(Qjsgs(:,:,:,3)))
end subroutine
