subroutine DoStandardDynamicProcedure(this, uE, vE, wE, uhatE, vhatE, whatE, duidxjEhat)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)), intent(in) :: uE, vE, wE
   complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3))  , intent(in) :: uhatE, vhatE, whatE
   complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3),9), intent(in) :: duidxjEhat
   integer :: idx
   
   ! STEP 1: Test filter velocities (required for Lij)
   call this%TestFilter_Cmplx_to_Real(uhatE, this%ui_Filt(:,:,:,1))   
   call this%TestFilter_Cmplx_to_Real(vhatE, this%ui_Filt(:,:,:,2))   
   call this%TestFilter_Cmplx_to_Real(whatE, this%ui_Filt(:,:,:,3))   

   ! STEP 2: Compute Lij
   this%buff1 = uE*uE
   call this%TestFilter_real_to_real(this%buff1, this%buff2)
   this%Lij(:,:,:,1) = -this%buff2 + this%ui_Filt(:,:,:,1)*this%ui_Filt(:,:,:,1)
   
   this%buff1 = uE*vE
   call this%TestFilter_real_to_real(this%buff1, this%buff2)
   this%Lij(:,:,:,2) = -this%buff2 + this%ui_Filt(:,:,:,1)*this%ui_Filt(:,:,:,2)

   this%buff1 = uE*wE
   call this%TestFilter_real_to_real(this%buff1, this%buff2)
   this%Lij(:,:,:,3) = -this%buff2 + this%ui_Filt(:,:,:,1)*this%ui_Filt(:,:,:,3)
   
   this%buff1 = vE*vE
   call this%TestFilter_real_to_real(this%buff1, this%buff2)
   this%Lij(:,:,:,4) = -this%buff2 + this%ui_Filt(:,:,:,2)*this%ui_Filt(:,:,:,2)
   
   this%buff1 = vE*wE
   call this%TestFilter_real_to_real(this%buff1, this%buff2)
   this%Lij(:,:,:,5) = -this%buff2 + this%ui_Filt(:,:,:,2)*this%ui_Filt(:,:,:,3)

   this%buff1 = wE*wE
   call this%TestFilter_real_to_real(this%buff1, this%buff2)
   this%Lij(:,:,:,6) = -this%buff2 + this%ui_Filt(:,:,:,3)*this%ui_Filt(:,:,:,3)

   
   ! STEP 3: Compute M_ij
   ! Part a: Compute \tilde{duidxj}
   do idx = 1,9
      call this%TestFilter_cmplx_to_real(duidxjEhat(:,:,:,idx),this%alphaij_Filt(:,:,:,idx))
   end do
   ! Part b: Compute \tilde{Sij}, NOTE: Mij is used to store filtered Sij
   call get_Sij_from_duidxj(this%alphaij_Filt, this%Mij, size(this%Mij,1), size(this%Mij,2), size(this%Mij,3))
   ! Part c: Compute \tilde{D_SGS}
   select case (this%mid)
   case (0) ! smagorinsky
      call get_smagorinsky_kernel(this%Mij,this%Dsgs_filt, &
                              this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3))
   case (1) ! sigma
      call get_sigma_kernel(this%Dsgs_filt, this%alphaij_Filt, &
                              this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3))
   end select
   ! Part d: Compute the rest of it
   do idx = 1,6
      this%buff1 = this%Dsgs*this%S_ij_E(:,:,:,idx)
      call this%TestFilter_real_to_real(this%buff1,this%buff2)
      this%buff1 = (this%deltaRat*this%deltaRat)*this%Dsgs_filt*this%Mij(:,:,:,idx)
      this%Mij(:,:,:,idx) = this%buff1 - this%buff2
   end do 


   ! STEP 4: Compute the numerator
   this%buff1 = this%Lij(:,:,:,1)*this%Mij(:,:,:,1)
   do idx = 2,6
      this%buff1 = this%buff1 + this%Lij(:,:,:,idx)*this%Mij(:,:,:,idx)
   end do 

   ! STEP 5: Compute the denominator
   this%buff2 = this%Mij(:,:,:,1)*this%Mij(:,:,:,1)
   do idx = 2,6
      this%buff2 = this%buff2 + this%Mij(:,:,:,idx)*this%Mij(:,:,:,idx)
   end do
   this%buff2 = 2.d0 * this%buff2

   ! STEP 6: Get the planar averages
   call this%planarAverage(this%buff1)
   call this%planarAverage(this%buff2)

   ! STEP 7: get the fraction and clip
   this%buff1 = this%buff1/(this%buff2 + 1.d-15)
   where (this%buff1 < 0.d0) 
      this%buff1 = 0.d0
   end where

   ! STEP 8: Interpolate to cells 
end subroutine




