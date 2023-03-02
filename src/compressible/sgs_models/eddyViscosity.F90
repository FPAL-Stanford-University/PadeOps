subroutine get_Sij_from_duidxj(this, duidxj, Sij)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL, this%nyL, this%nzL, 9), intent(in)  :: duidxj
   real(rkind), dimension(this%nxL, this%nyL, this%nzL, 6), intent(out) :: Sij
   integer :: i, j, k

   do k = 1, this%nzL
      do j = 1, this%nyL
         !$omp simd
         do i = 1, this%nxL
            Sij(i,j,k,1) = duidxj(i,j,k,1) ! S11 = dudx
            Sij(i,j,k,4) = duidxj(i,j,k,5) ! S22 = dvdy
            Sij(i,j,k,6) = duidxj(i,j,k,9) ! S33 = dwdz

            Sij(i,j,k,2) = 0.5d0*(duidxj(i,j,k,2) + duidxj(i,j,k,4)) ! S12 = 0.5*(dudy + dvdx)
            Sij(i,j,k,3) = 0.5d0*(duidxj(i,j,k,3) + duidxj(i,j,k,7)) ! S13 = 0.5*(dudz + dwdx)
            Sij(i,j,k,5) = 0.5d0*(duidxj(i,j,k,6) + duidxj(i,j,k,8)) ! S23 = 0.5*(dvdz + dwdy)
         end do 
      end do 
   end do 
end subroutine

!subroutine get_modS_Sii(this, modS_sq, Sii)
subroutine get_modS_Sii(this,Sij, modS_sq, Sii)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,6), intent(in)  :: Sij
   real(rkind), dimension(this%nxL, this%nyL, this%nzL),    intent(out) :: modS_sq, Sii

   integer :: i, j, k

   modS_sq =                 Sij(:,:,:,1) * Sij(:,:,:,1)  ! S11*S11
   modS_sq = modS_sq + two * Sij(:,:,:,2) * Sij(:,:,:,2)  ! S12*S12 + S21*S21
   modS_sq = modS_sq + two * Sij(:,:,:,3) * Sij(:,:,:,3)  ! S13*S13 + S31*S31
   modS_sq = modS_sq +       Sij(:,:,:,4) * Sij(:,:,:,4)  ! S22*S22
   modS_sq = modS_sq + two * Sij(:,:,:,5) * Sij(:,:,:,5)  ! S23*S23 + S32*S32
   modS_sq = modS_sq +       Sij(:,:,:,6) * Sij(:,:,:,6)  ! S33*S33
   modS_sq = two * modS_sq                                            ! (2 * S_ij * S_ij)
   
   Sii = third * (Sij(:,:,:,1) + Sij(:,:,:,4) + Sij(:,:,:,6))

end subroutine


subroutine get_qtke(this, rho, modS_sq, qtke)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL,this%nyL,this%nzL), intent(in) :: rho, modS_sq
   real(rkind), dimension(this%nxL,this%nyL,this%nzL), intent(out):: qtke

   qtke = twothird * (this%deltaLES**2) * rho * modS_sq

end subroutine

!subroutine get_SGS_kernel(this,duidxj, modS_sq)
subroutine get_SGS_kernel(this,duidxj, Sij, modS_sq, nusgs)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,9), intent(in) :: duidxj
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,6), intent(in)  :: Sij
   real(rkind), dimension(this%nxL,this%nyL,this%nzL),   intent(in) :: modS_sq
   real(rkind), dimension(this%nxL,this%nyL,this%nzL),   intent(out) :: nusgs

   select case(this%SGSmodelID) 
   case (0)
      ! Smagorinsky
      call this%get_smagorinsky_kernel(modS_sq, nusgs)
   case (1)
      ! Sigma
      call this%get_sigma_kernel(duidxj, nusgs)
   case (2)
      ! AMD  
     !call this%get_amd_kernel(duidxj)
      call this%get_amd_kernel(duidxj, Sij, nusgs)
   end select

end subroutine

subroutine get_Qjsgs_eddy_kernel(this, rho, nusgs, gradT, Qjsgs)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL,this%nyL,this%nzL),   intent(in)  :: rho, nusgs
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,3), intent(in)  :: gradT
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,3), intent(out) :: Qjsgs
   integer :: i,j,k
  
  ! if (this%isTurbPrandtlconst) then
  !    this%kapsgs = nusgs ! /this%PrSGS
  ! else
  !   ! call this%get_amd_Dkernel(duidxj, gradT, this%kapsgs)
  ! endif
  !    
   do k = 1, 3
      !Qjsgs(:,:,:,k) = - rho * this%Cp * this%kapsgs * gradT(:,:,:,k)
      Qjsgs(:,:,:,k) = - rho * this%Cp * nusgs * gradT(:,:,:,k)
   end do

end subroutine


subroutine multiply_by_model_constant(this,qtke)
   class(sgs_cgrid), intent(inout) :: this 
   real(rkind), dimension(this%nxL,this%nyL,this%nzL), intent(inout)  :: qtke
 
   integer  :: j, k

   if( (this%DynamicProcedureType==0) .or. (this%DynamicProcedureType==2)) then
      !! constant coefficient or Global-Dynamic Procedure
      this%nusgs = this%cmodel_global * this%nusgs
      qtke       = this%Ctke  * qtke
   elseif(this%DynamicProcedureType==1) then
      !! Local-Dynamic Procedure - averaged in (x,z); function of (y)
      do k = 1, this%nzL
         do j = 1, this%nyL
            this%nusgs(:,j,k) = this%cmodel_local(j)     * this%nusgs(:,j,k)
            qtke(:,j,k)       = this%cmodel_local_tke(j) * qtke(:,j,k)
         end do
      end do
   endif
     
end subroutine

subroutine multiply_by_model_coefficient_eddy_Qjsgs(this, Qjsgs)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,3), intent(inout)  :: Qjsgs

   integer  :: ii, j, k

   if( (this%DynamicProcedureType==0) .or. (this%DynamicProcedureType==2)) then
      !! constant coefficient or Global-Dynamic Procedure
      do ii = 1, 3
          Qjsgs(:,:,:,ii) = this%cmodel_global_Qjsgs * Qjsgs(:,:,:,ii)
      enddo 
   elseif(this%DynamicProcedureType==1) then
      !! Local-Dynamic Procedure - averaged in (x,z); function of (y)
      do ii = 1, 3
        do k = 1, this%nzL
          do j = 1, this%nyL
            Qjsgs(:,j,k,ii) = this%cmodel_local_Qjsgs(j) * Qjsgs(:,j,k,ii)
          end do
        end do
      end do
   endif

end subroutine
