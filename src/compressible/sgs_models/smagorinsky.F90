subroutine init_smagorinsky(this)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind) :: deltaLES

   ! Set the type of mnodel constant (default is wall function). 
   ! Can be reset to true via dynamic procedure initialization, 
   ! in case the global dynamic procedure is used. 
      
   !this%useCglobal = .true.
   !this%cmodel_global = (this%Csgs*this%deltaLES)**2
   this%isEddyViscosityModel = .true. 
   this%cmodel_global = this%Csgs !* this%deltaLES**2 ! one       
   this%cmodel_global_Qjsgs =  this%Csgs/this%Prsgs

   call message(1,"Smagorinsky model initialized")
end subroutine

subroutine destroy_smagorinsky(this)
   class(sgs_cgrid), intent(inout) :: this

   this%isEddyViscosityModel = .false. 

end subroutine

subroutine get_smagorinsky_kernel(this, modS_sq, nusgs)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), dimension(this%nxL,this%nyL,this%nzL), intent(in)  :: modS_sq
   real(rkind), dimension(this%nxL,this%nyL,this%nzL), intent(out) :: nusgs
   !real(rkind), dimension(nxL, nyL) :: S
   
   !!do k = 1,nzL
   !!   do j = 1,nyL
   !!      !$omp simd 
   !!      do i = 1,nxL
   !!         S = Sij(i,j,k,1)*Sij(i,j,k,1) ! S11*S11
   !!         S = S + 2.d0*(Sij(i,j,k,2)*Sij(i,j,k,2)) ! S12*S12 + S21*S21
   !!         S = S + 2.d0*(Sij(i,j,k,3)*Sij(i,j,k,3)) ! S13*S13 + S31*S31
   !!         S = S + (Sij(i,j,k,4)*Sij(i,j,k,4)) ! S22*S22
   !!         S = S + 2.d0*(Sij(i,j,k,5)*Sij(i,j,k,5)) ! S23*S23 + S32*S32
   !!         S = S + (Sij(i,j,k,6)*Sij(i,j,k,6)) ! S33*S33

   !!         ! Now do modS = sqrt(2* S_ij*S_ij)
   !!         S = 2.d0*S
   !!         nusgs(i,j,k) = sqrt(S)
   !!      end do 
   !!   end do 
   !!end do

   nusgs = sqrt(modS_sq) 
   nusgs = sqrt(modS_sq) * (this%deltaLES**2)



end subroutine
