subroutine get_Sij_from_duidxj(this, duidxj)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL, this%nyL, this%nzL, 9), intent(in) :: duidxj
   integer :: i, j, k

   do k = 1, this%nzL
      do j = 1, this%nyL
         !$omp simd
         do i = 1, this%nxL
            this%S_ij(i,j,k,1) = duidxj(i,j,k,1) ! S11 = dudx
            this%S_ij(i,j,k,4) = duidxj(i,j,k,5) ! S22 = dvdy
            this%S_ij(i,j,k,6) = duidxj(i,j,k,9) ! S33 = dwdz

            this%S_ij(i,j,k,2) = 0.5d0*(duidxj(i,j,k,2) + duidxj(i,j,k,4)) ! S12 = 0.5*(dudy + dvdx)
            this%S_ij(i,j,k,3) = 0.5d0*(duidxj(i,j,k,3) + duidxj(i,j,k,7)) ! S13 = 0.5*(dudz + dwdx)
            this%S_ij(i,j,k,5) = 0.5d0*(duidxj(i,j,k,6) + duidxj(i,j,k,8)) ! S23 = 0.5*(dvdz + dwdy)
         end do 
      end do 
   end do 
end subroutine

subroutine get_modS_Sii(this, modS_sq, Sii)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL, this%nyL, this%nzL, 9), intent(in)  :: duidxj
   real(rkind), dimension(this%nxL, this%nyL, this%nzL),    intent(out) :: modS_sq, Sii
   integer :: i, j, k

   modS_sq =              this%S_ij(:,:,:,1) * this%S_ij(:,:,:,1)  ! S11*S11
   modS_sq = modS + two * this%S_ij(:,:,:,2) * this%S_ij(:,:,:,2)  ! S12*S12 + S21*S21
   modS_sq = modS + two * this%S_ij(:,:,:,3) * this%S_ij(:,:,:,3)  ! S13*S13 + S31*S31
   modS_sq = modS +       this%S_ij(:,:,:,4) * this%S_ij(:,:,:,4)  ! S22*S22
   modS_sq = modS + two * this%S_ij(:,:,:,5) * this%S_ij(:,:,:,5)  ! S23*S23 + S32*S32
   modS_sq = modS +       this%S_ij(:,:,:,6) * this%S_ij(:,:,:,6)  ! S33*S33
   modS_sq = two * modS                                            ! (2 * S_ij * S_ij)
   
   Sii = third * (this%S_ij(:,:,:,1) + this%S_ij(:,:,:,4) + this%S_ij(:,:,:,6))

end subroutine

subroutine get_SGS_kernel(this,duidxj, modS_sq)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,9), intent(in) :: duidxj
   real(rkind), dimension(this%nxL,this%nyL,this%nzL),   intent(in) :: modS_sq

   select case(this%SGSmodelID) 
   case (0)
      ! Smagorinsky
      call get_smagorinsky_kernel(modS_sq)
   case (1)
      ! Sigma
      call get_sigma_kernel(duidxj)
   case (2)
      ! AMD  
      call get_amd_kernel(duidxj)
   end select

end subroutine

subroutine multiply_by_model_constant(this)
  class(sgs_cgrid), intent(inout) :: this 
 
  this%nusgs = this%cmodel_global*this%nusgs
     
end subroutine
