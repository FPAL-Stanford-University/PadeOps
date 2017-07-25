subroutine init_amd(this, dx, dy, dz, Csgs)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), intent(in) :: dx, dy, dz, Csgs

   this%useCglobal = .true. 
   this%isEddyViscosityModel = .true. 
   
   this%camd_x = Csgs*dx*sqrt(1.d0/12.d0)
   this%camd_y = Csgs*dy*sqrt(1.d0/12.d0)
   this%camd_z = Csgs*dz*this%PadeDer%getApproxPoincareConstant()
   this%cmodel_global = one  ! Anisotropic model constants 
   call message(1,"AMD model initialized")

end subroutine

subroutine destroy_amd(this)
   class(sgs_igrid), intent(inout) :: this
   this%isEddyViscosityModel = .true. 
end subroutine



subroutine get_amd_kernel(nu_sgs, camd_x, camd_y, camd_z, duidxj, Sij, nxL, nyL, nzL)
   integer, intent(in) :: nxL, nyL, nzL
   real(rkind), intent(in) :: camd_x, camd_y, camd_z
   real(rkind), intent(in), dimension(nxL,nyL,nzL,9) :: duidxj
   real(rkind), intent(in), dimension(nxL,nyL,nzL,6) :: Sij
   real(rkind), intent(out), dimension(nxL,nyL,nzL)  :: nu_sgs
 
   real(rkind) :: num, den
   integer :: i, j, k

   !dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3)
   !dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6)
   !dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9)

   !S11 => Sij(:,:,:,1); S12 => Sij(:,:,:,2); S13 => Sij(:,:,:,3)
   !S22 => Sij(:,:,:,4); S23 => Sij(:,:,:,5); S33 => Sij(:,:,:,6)

   do k = 1,nzL
      do j = 1,nyL
         !$omp simd 
         do i = 1,nxL
            
            ! First compute the numerator
            num =      ((duidxj(i,j,k,1)*camd_x)*(duidxj(i,j,k,1)*camd_x) + (duidxj(i,j,k,2)*camd_y)*(duidxj(i,j,k,2)*camd_y) &
                        + (duidxj(i,j,k,3)*camd_z)*(duidxj(i,j,k,3)*camd_z))*Sij(i,j,k,1)

            num = num + ((duidxj(i,j,k,4)*camd_x)*(duidxj(i,j,k,4)*camd_x) + (duidxj(i,j,k,5)*camd_y)*(duidxj(i,j,k,5)*camd_y) &
                        + (duidxj(i,j,k,6)*camd_z)*(duidxj(i,j,k,6)*camd_z))*Sij(i,j,k,4)

            num = num + ((duidxj(i,j,k,7)*camd_x)*(duidxj(i,j,k,7)*camd_x) + (duidxj(i,j,k,8)*camd_y)*(duidxj(i,j,k,8)*camd_y) &
                        + (duidxj(i,j,k,9)*camd_z)*(duidxj(i,j,k,9)*camd_z))*Sij(i,j,k,6)

            num = num + two*((duidxj(i,j,k,1)*camd_x)*(duidxj(i,j,k,4)*camd_x) + (duidxj(i,j,k,2)*camd_y)*(duidxj(i,j,k,5)*camd_y) &
                        + (duidxj(i,j,k,3)*camd_z)*(duidxj(i,j,k,6)*camd_z))*Sij(i,j,k,2)
            
            num = num + two*((duidxj(i,j,k,1)*camd_x)*(duidxj(i,j,k,7)*camd_x) + (duidxj(i,j,k,2)*camd_y)*(duidxj(i,j,k,8)*camd_y) &
                        + (duidxj(i,j,k,3)*camd_z)*(duidxj(i,j,k,9)*camd_z))*Sij(i,j,k,3)

            num = num + two*((duidxj(i,j,k,4)*camd_x)*(duidxj(i,j,k,7)*camd_x) + (duidxj(i,j,k,5)*camd_y)*(duidxj(i,j,k,8)*camd_y) &
                        + (duidxj(i,j,k,6)*camd_z)*(duidxj(i,j,k,9)*camd_z))*Sij(i,j,k,5)

            ! Now compute the denominator
            den =  duidxj(i,j,k,1)*duidxj(i,j,k,1) +  duidxj(i,j,k,2)*duidxj(i,j,k,2) + duidxj(i,j,k,3)*duidxj(i,j,k,3) &
                 + duidxj(i,j,k,4)*duidxj(i,j,k,4) +  duidxj(i,j,k,5)*duidxj(i,j,k,5) + duidxj(i,j,k,6)*duidxj(i,j,k,6) &
                 + duidxj(i,j,k,7)*duidxj(i,j,k,7) +  duidxj(i,j,k,8)*duidxj(i,j,k,8) + duidxj(i,j,k,9)*duidxj(i,j,k,9) 
                     
            nu_sgs(i,j,k) = max(-num/(den + 1.d-32),zero) 
         end do 
      end do
   end do 
end subroutine 
