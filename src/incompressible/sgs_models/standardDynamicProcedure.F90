subroutine DoStandardDynamicProcedureScalar(this, u, v, w, T, That, dTdx, dTdy, dTdz)
   class(sgs_igrid), intent(inout), target :: this
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)),   intent(in) :: u, v, w, T
   complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3))  , intent(in) :: That
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: dTdx, dTdy, dTdz
  
   real(rkind) :: deltaRat_sq
   real(rkind), dimension(:,:,:),   pointer :: buff1, T_Filt, numerator, denominator
   real(rkind), dimension(:,:,:,:), pointer :: Ri, Pi, dTdxi_filt

   ! =========================== SETUP =================================
   Ri => this%Lij(:,:,:,1:3)
   Pi => this%Lij(:,:,:,4:6)
   dTdxi_filt => this%Sij_filt(:,:,:,1:3)
   T_Filt => this%Sij_filt(:,:,:,4)
   numerator   => this%rbuffxC(:,:,:,1)
   denominator => this%rbuffxC(:,:,:,2)
   buff1 => this%rbuffxC(:,:,:,3)


   ! =========================== CORE  =================================
   ! Step 1: Test filter dTdxi
   call this%TestFilter_Real_to_Real(dTdx,dTdxi_filt(:,:,:,1))
   call this%TestFilter_Real_to_Real(dTdy,dTdxi_filt(:,:,:,2))
   call this%TestFilter_Real_to_Real(dTdz,dTdxi_filt(:,:,:,3))

   ! Step 2: Do nothing... already have the smagorinsky kernels in
   ! this%Dsgs_filt and this%Dsgs

   ! Step 3: Compute Ri
   deltaRat_sq = this%DeltaRat*this%DeltaRat
   buff1 = this%Dsgs*dTdx
   call this%TestFilter_Real_To_Real_inplace(buff1)
   Ri(:,:,:,1) = buff1 - deltaRat_sq*this%Dsgs_filt*dTdxi_filt(:,:,:,1)
   
   buff1 = this%Dsgs*dTdy
   call this%TestFilter_Real_To_Real_inplace(buff1)
   Ri(:,:,:,2) = buff1 - deltaRat_sq*this%Dsgs_filt*dTdxi_filt(:,:,:,2)
   
   buff1 = this%Dsgs*dTdz
   call this%TestFilter_Real_To_Real_inplace(buff1)
   Ri(:,:,:,3) = buff1 - deltaRat_sq*this%Dsgs_filt*dTdxi_filt(:,:,:,3)

   ! Step 4: Compute T_filt (already have ui_filt in this%ui_filt)
   call this%TestFilter_Cmplx_to_Real(That, T_Filt)   

   ! Step 5: Compute Pi
   buff1 = u*T
   call this%TestFilter_Real_to_Real_inplace(buff1)
   Pi(:,:,:,1) = buff1 - this%ui_Filt(:,:,:,1)*T_Filt

   buff1 = v*T
   call this%TestFilter_Real_to_Real_inplace(buff1)
   Pi(:,:,:,2) = buff1 - this%ui_Filt(:,:,:,2)*T_Filt

   buff1 = w*T
   call this%TestFilter_Real_to_Real_inplace(buff1)
   Pi(:,:,:,3) = buff1 - this%ui_Filt(:,:,:,3)*T_Filt

   ! Step 6: Compute the numerator
   numerator = Ri(:,:,:,1)*Pi(:,:,:,1) + Ri(:,:,:,2)*Pi(:,:,:,2) + Ri(:,:,:,3)*Pi(:,:,:,3)

   ! Step 7: Compute the denominator 
   denominator = Ri(:,:,:,1)*Ri(:,:,:,1) + Ri(:,:,:,2)*Ri(:,:,:,2) + Ri(:,:,:,3)*Ri(:,:,:,3) 


   ! =========================  WRAPUP  ================================
   call this%planarAverage_and_TakeRatio(numerator,denominator,this%BetaDynProc_C, this%Beta_1d)


end subroutine 


subroutine DoStandardDynamicProcedure(this, u, v, w, uhat, vhat, what, Sij)
   class(sgs_igrid), intent(inout), target :: this
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)),   intent(in) :: u, v, w
   complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3))  , intent(in) :: uhat, vhat, what
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),6), intent(in) :: Sij
   
   real(rkind) :: deltaRat_sq
   real(rkind), dimension(:,:,:),   pointer ::  numerator, denominator, buff1
   real(rkind), dimension(:,:,:,:), pointer :: Mij
   integer :: idx 

   ! =========================== SETUP =================================
   Mij => this%tau_ij
   numerator   => this%rbuffxC(:,:,:,1)
   denominator => this%rbuffxC(:,:,:,2)
   buff1 => this%rbuffxC(:,:,:,3)

   ! =========================== CORE  =================================
   
   ! Step 1: Test filter Sij
   do idx = 1,6
      call this%TestFilter_Real_To_Real(Sij(:,:,:,idx), this%Sij_filt(:,:,:,idx))
   end do 

   ! Step 2: Compute Smagorinsky kernel using test filtered Sij
   call get_smagorinsky_kernel(this%Sij_filt,this%Dsgs_filt,this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3))


   ! Step 3: Compute Mij
   deltaRat_sq = this%DeltaRat*this%DeltaRat
   do idx = 1,6
      buff1 = this%Dsgs*Sij(:,:,:,idx)
      call this%TestFilter_Real_To_Real_inplace(buff1)
      Mij(:,:,:,idx) = buff1 - deltaRat_sq*this%Dsgs_filt*this%Sij_Filt(:,:,:,idx)
   end do 

   ! Step 4: Compute ui_filt
   call this%TestFilter_Cmplx_to_Real(uhat, this%ui_Filt(:,:,:,1))   
   call this%TestFilter_Cmplx_to_Real(vhat, this%ui_Filt(:,:,:,2))   
   call this%TestFilter_Cmplx_to_Real(what, this%ui_Filt(:,:,:,3))   

   ! Step 5: Compute Lij
   ! L11
   buff1 = u*u
   call this%TestFilter_Real_to_Real_inplace(buff1)
   this%Lij(:,:,:,1) = buff1 - this%ui_Filt(:,:,:,1)*this%ui_Filt(:,:,:,1)
   
   !L12
   buff1 = u*v
   call this%TestFilter_Real_to_Real_inplace(buff1)
   this%Lij(:,:,:,2) = buff1 - this%ui_Filt(:,:,:,1)*this%ui_Filt(:,:,:,2)

   !L13
   buff1 = u*w
   call this%TestFilter_Real_to_Real_inplace(buff1)
   this%Lij(:,:,:,3) = buff1 - this%ui_Filt(:,:,:,1)*this%ui_Filt(:,:,:,3)

   !L22
   buff1 = v*v
   call this%TestFilter_Real_to_Real_inplace(buff1)
   this%Lij(:,:,:,4) = buff1 - this%ui_Filt(:,:,:,2)*this%ui_Filt(:,:,:,2)

   !L23
   buff1 = v*w
   call this%TestFilter_Real_to_Real_inplace(buff1)
   this%Lij(:,:,:,5) = buff1 - this%ui_Filt(:,:,:,2)*this%ui_Filt(:,:,:,3)

   !L33
   buff1 = w*w
   call this%TestFilter_Real_to_Real_inplace(buff1)
   this%Lij(:,:,:,6) = buff1 - this%ui_Filt(:,:,:,3)*this%ui_Filt(:,:,:,3)


   ! Step 6: Compute the numerator LijMij
   numerator = this%Lij(:,:,:,1)*Mij(:,:,:,1)
   do idx = 2,6
      numerator = numerator + this%Lij(:,:,:,idx)*Mij(:,:,:,idx)
   end do 
   numerator = numerator + this%Lij(:,:,:,2)*Mij(:,:,:,2) ! Off diagonal M21*L21 = M12*L12 needs to be counted twice 
   numerator = numerator + this%Lij(:,:,:,3)*Mij(:,:,:,3) ! Off diagonal M31*L31 = M13*L13 needs to be counted twice
   numerator = numerator + this%Lij(:,:,:,5)*Mij(:,:,:,5) ! Off diagonal M32*L32 = M23*L23 needs to be counted twice

   ! Step 7: Compute the denominator MijMij
   denominator = Mij(:,:,:,1)*Mij(:,:,:,1)
   do idx = 2,6
      denominator = denominator + Mij(:,:,:,idx)*Mij(:,:,:,idx)
   end do 
   denominator = denominator + Mij(:,:,:,2)*Mij(:,:,:,2)  ! Off diagonal M21*M21 = M12*M12 needs to be counted twice
   denominator = denominator + Mij(:,:,:,3)*Mij(:,:,:,3)  ! Off diagonal M31*M31 = M13*M13 needs to be counted twice
   denominator = denominator + Mij(:,:,:,5)*Mij(:,:,:,5)  ! Off diagonal M32*M32 = M23*M23 needs to be counted twice
   
   denominator = 2.d0*denominator
   
   ! =========================  WRAPUP  ================================
   call this%planarAverage_and_TakeRatio(numerator,denominator,this%LambdaDynProc_C, this%Lambda_1d)


   ! What happens after this?
   ! this%nSGS is computed as this%nuSGS = this%LambdaDynProc_C*this%this%Dsgs
end subroutine 

subroutine TestFilter_Real_to_Real(this, f,ffil)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)), intent(out)  :: ffil
   
   call this%spectC%fft(f, this%Tfilhat)
   call this%spectC%testFilter_ip(this%Tfilhat)
   if (this%useVerticalTfilter) then
      call transpose_y_to_z(this%Tfilhat, this%Tfilhatz1, this%sp_gpC)
      call this%gaussianTestFilterZ%filter3(this%Tfilhatz1, this%Tfilhatz2, this%sp_gpC%zsz(1), this%sp_gpC%zsz(2))
      call transpose_z_to_y(this%Tfilhatz2, this%Tfilhat, this%sp_gpC)
   end if
   call this%spectC%ifft(this%Tfilhat, ffil)

end subroutine

subroutine TestFilter_Real_to_Real_inplace(this, f)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)), intent(inout)  :: f
   
   call this%spectC%fft(f, this%Tfilhat)
   call this%spectC%testFilter_ip(this%Tfilhat)
   if (this%useVerticalTfilter) then
      call transpose_y_to_z(this%Tfilhat, this%Tfilhatz1, this%sp_gpC)
      call this%gaussianTestFilterZ%filter3(this%Tfilhatz1, this%Tfilhatz2, this%sp_gpC%zsz(1), this%sp_gpC%zsz(2))
      call transpose_z_to_y(this%Tfilhatz2, this%Tfilhat, this%sp_gpC)
   end if
   call this%spectC%ifft(this%Tfilhat, f)

end subroutine

subroutine TestFilter_Cmplx_to_Real(this, fhat, f)
   class(sgs_igrid), intent(inout) :: this
   complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(1),this%sp_gpC%ysz(1)), intent(in) :: fhat
   real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)), intent(out)  :: f

   call this%spectC%TestFilter_oop(fhat, this%Tfilhat)
   if (this%useVerticalTfilter) then
      call transpose_y_to_z(this%Tfilhat, this%Tfilhatz1, this%sp_gpC)
      call this%gaussianTestFilterZ%filter3(this%Tfilhatz1, this%Tfilhatz2, this%sp_gpC%zsz(1), this%sp_gpC%zsz(2))
      call transpose_z_to_y(this%Tfilhatz2, this%Tfilhat, this%sp_gpC)
   end if
   call this%spectC%ifft(this%Tfilhat, f)

end subroutine

subroutine planarAverage_and_TakeRatio(this, num, den, ratC, rat1d)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)), intent(in)  :: num, den
   real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)), intent(out) :: ratC
   real(rkind), dimension(this%gpC%zsz(3)), intent(out) :: rat1d
   integer :: idx
   
   call transpose_x_to_y(num, this%rbuffyC(:,:,:,1), this%gpC)
   call transpose_y_to_z(this%rbuffyC(:,:,:,1), this%rbuffzC(:,:,:,1), this%gpC)
   do idx = 1,this%gpC%zsz(3)
      this%rbuffzC(:,:,idx,1) = max(p_sum(sum(this%rbuffzC(:,:,idx,1)))*this%meanfact, zero)
   end do 
   
   call transpose_x_to_y(den, this%rbuffyC(:,:,:,1), this%gpC)
   call transpose_y_to_z(this%rbuffyC(:,:,:,1), this%rbuffzC(:,:,:,2), this%gpC)
   do idx = 1,this%gpC%zsz(3)
      this%rbuffzC(:,:,idx,2) = this%rbuffzC(1,1,idx,1)/(p_sum(sum(this%rbuffzC(:,:,idx,2)))*this%meanfact + 1.d-18)
   end do

   rat1d = this%rbuffzC(1,1,:,2)
   call transpose_z_to_y(this%rbuffzC(:,:,:,2), this%rbuffyC(:,:,:,1), this%gpC)
   call transpose_y_to_x(this%rbuffyC(:,:,:,1), ratC, this%gpC)

end subroutine









