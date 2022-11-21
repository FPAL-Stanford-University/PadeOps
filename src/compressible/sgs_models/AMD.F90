subroutine init_amd(this)
   use constants, only: one
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), intent(in) :: dx, dy, dz, Csgs

   this%useCglobal = .true. 
   this%isEddyViscosityModel = .true. 

   if(this%isEddyDiffModel) then
       call message(1, "!!!! WARNING !!!! AMD model is being used as an Eddy Diffusivity Model")
       call message(1, "!!!! WARNING !!!! PrSGS = ", this%PrSGS)
   endif
   
   this%camd_x = this%Csgs*1.5d0*this%dx*sqrt(1.d0/12.0d0)
   this%camd_y = this%Csgs*1.5d0*this%dy*sqrt(1.d0/12.0d0)
   this%camd_z = this%Csgs*1.5d0*this%dz*sqrt(1.d0/12.0d0)          !this%PadeDer%getApproxPoincareConstant()
   this%cmodel_global = one                                    ! Anisotropic model constants 
   call message(1,"AMD model initialized")

end subroutine

subroutine destroy_amd(this)
   class(sgs_cgrid), intent(inout) :: this

   this%isEddyViscosityModel = .true. 

end subroutine

subroutine get_amd_kernel(duidxj)
   real(rkind), intent(in), dimension(this%nxL,this%nyL,this%nzL,9) :: duidxj
   real(rkind) :: num, den
   integer :: i, j, k

   !dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3)
   !dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6)
   !dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9)

   !S11 => Sij(:,:,:,1); S12 => Sij(:,:,:,2); S13 => Sij(:,:,:,3)
   !S22 => Sij(:,:,:,4); S23 => Sij(:,:,:,5); S33 => Sij(:,:,:,6)


   do k = 1, this%nzL
      do j = 1, this%nyL
         !$omp simd 
         do i = 1, this%nxL
            
            ! First compute the numerator
            num =           ((duidxj(i,j,k,1)*this%camd_x)*(duidxj(i,j,k,1)*this%camd_x) &
                      +      (duidxj(i,j,k,2)*this%camd_y)*(duidxj(i,j,k,2)*this%camd_y) &
                      +      (duidxj(i,j,k,3)*this%camd_z)*(duidxj(i,j,k,3)*this%camd_z)) * Sij(i,j,k,1)

            num = num +     ((duidxj(i,j,k,4)*this%camd_x)*(duidxj(i,j,k,4)*this%camd_x) &
                      +      (duidxj(i,j,k,5)*this%camd_y)*(duidxj(i,j,k,5)*this%camd_y) &
                      +      (duidxj(i,j,k,6)*this%camd_z)*(duidxj(i,j,k,6)*this%camd_z)) * Sij(i,j,k,4)

            num = num +     ((duidxj(i,j,k,7)*this%camd_x)*(duidxj(i,j,k,7)*this%camd_x) &
                      +      (duidxj(i,j,k,8)*this%camd_y)*(duidxj(i,j,k,8)*this%camd_y) &
                      +      (duidxj(i,j,k,9)*this%camd_z)*(duidxj(i,j,k,9)*this%camd_z)) * Sij(i,j,k,6)

            num = num + two*((duidxj(i,j,k,1)*this%camd_x)*(duidxj(i,j,k,4)*this%camd_x) &
                           + (duidxj(i,j,k,2)*this%camd_y)*(duidxj(i,j,k,5)*this%camd_y) &
                           + (duidxj(i,j,k,3)*this%camd_z)*(duidxj(i,j,k,6)*this%camd_z)) * Sij(i,j,k,2)
            
            num = num + two*((duidxj(i,j,k,1)*this%camd_x)*(duidxj(i,j,k,7)*this%camd_x) &
                           + (duidxj(i,j,k,2)*this%camd_y)*(duidxj(i,j,k,8)*this%camd_y) &
                           + (duidxj(i,j,k,3)*this%camd_z)*(duidxj(i,j,k,9)*this%camd_z)) * Sij(i,j,k,3)

            num = num + two*((duidxj(i,j,k,4)*this%camd_x)*(duidxj(i,j,k,7)*this%camd_x) &
                           + (duidxj(i,j,k,5)*this%camd_y)*(duidxj(i,j,k,8)*this%camd_y) &
                           + (duidxj(i,j,k,6)*this%camd_z)*(duidxj(i,j,k,9)*this%camd_z)) * Sij(i,j,k,5)

            ! Now compute the denominator
            den =  duidxj(i,j,k,1)*duidxj(i,j,k,1) +  duidxj(i,j,k,2)*duidxj(i,j,k,2) + duidxj(i,j,k,3)*duidxj(i,j,k,3) &
                 + duidxj(i,j,k,4)*duidxj(i,j,k,4) +  duidxj(i,j,k,5)*duidxj(i,j,k,5) + duidxj(i,j,k,6)*duidxj(i,j,k,6) &
                 + duidxj(i,j,k,7)*duidxj(i,j,k,7) +  duidxj(i,j,k,8)*duidxj(i,j,k,8) + duidxj(i,j,k,9)*duidxj(i,j,k,9) 
                     
            this%nusgs(i,j,k) = max(-num/(den + 1.d-15),zero) 
         end do 
      end do
   end do 
end subroutine 



subroutine get_amd_Dkernel(duidxj, gradT)
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,9), intent(in) :: duidxj
   real(rkind), dimension(this%nxL,this%nyL,this%nzL,3), intent(in) :: gradT
   real(rkind) :: num, den
   integer :: i, j, k


   do k = 1,this%nzL
     do j = 1,this%nyL
        !$omp simd 
        do i = 1,this%nxL
           ! First compute the numerator
           num =       ((duidxj(i,j,k,1)*this%camd_x)*(gradT(i,j,k,1)*this%camd_x) &
                     +  (duidxj(i,j,k,2)*this%camd_y)*(gradT(i,j,k,2)*this%camd_y) &
                     +  (duidxj(i,j,k,3)*this%camd_z)*(gradT(i,j,k,3)*this%camd_z)) * gradT(i,j,k,1)

           num = num + ((duidxj(i,j,k,4)*this%camd_x)*(gradT(i,j,k,1)*this%camd_x) &
                     +  (duidxj(i,j,k,5)*this%camd_y)*(gradT(i,j,k,2)*this%camd_y) &
                     +  (duidxj(i,j,k,6)*this%camd_z)*(gradT(i,j,k,3)*this%camd_z)) * gradT(i,j,k,2)

           num = num + ((duidxj(i,j,k,7)*this%camd_x)*(gradT(i,j,k,1)*this%camd_x) &
                      + (duidxj(i,j,k,8)*this%camd_y)*(gradT(i,j,k,2)*this%camd_y) &
                      + (duidxj(i,j,k,9)*this%camd_z)*(gradT(i,j,k,3)*this%camd_z)) * gradT(i,j,k,3)

           ! Now get the denominator
           den = gradT(i,j,k,1)*gradT(i,j,k,1) + gradT(i,j,k,2)*gradT(i,j,k,2) +  gradT(i,j,k,3)*gradT(i,j,k,3)  

           this%kappasgs(i,j,k) = max(-num/(den + 1.d-15),zero)
        end do 
      end do
   end do 

end subroutine 
