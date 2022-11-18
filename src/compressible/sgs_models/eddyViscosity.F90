subroutine get_Sij_from_duidxj(duidxj, Sij, nxL, nyL, nzL)
   integer, intent(in) :: nxL, nyL, nzL
   real(rkind), dimension(nxL,nyL,nzL,9), intent(in) :: duidxj
   real(rkind), dimension(nxL,nyL,nzL,6), intent(out) :: Sij
   integer :: i, j, k

   do k = 1,nzL
      do j = 1,nyL
         !$omp simd
         do i = 1,nxL
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


subroutine get_SGS_kernel(this,duidxj)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3),9), intent(in) :: duidxj

   select case(this%mid) 
   case (0)
      ! Smagorinsky
      call get_smagorinsky_kernel(this%S_ij,this%nusgs,this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3))
   case (1)
      ! Sigma
      call get_sigma_kernel(this%nusgs, duidxj,this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3))
   case (2)
      ! AMD  
      call get_amd_kernel(this%nusgs, this%camd_x, this%camd_y, this%camd_z, duidxj, this%S_ij,this%decomp%ysz(1),this%decomp%ysz(2),this%decomp%ysz(3))
   end select

end subroutine

subroutine multiply_by_model_constant(this)
  class(sgs_cgrid), intent(inout) :: this 
 
  this%nu_sgs = this%cmodel_global*this%nusgs
     
end subroutine
