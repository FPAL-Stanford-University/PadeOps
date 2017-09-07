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
   call this%get_SGS_kernel_E(this%alphaij_Filt, this%Mij, this%Dsgs_filt)
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

   ! STEP 6: Get the planar average and interpolate
   select case(this%AverageType)
   case(1) ! Planar Average
     !call this%planarAverage(this%buff1)
     !call this%planarAverage(this%buff2)
     call this%planarAverageAndInterpolateToCells(this%buff1, this%buff2, this%rbuffxC(:,:,:,1))
     this%cmodelE = this%buff1(1,1,:)
     this%cmodelC = this%rbuffxC(1,1,:,1)
   case(2) ! Gaussian Filter
     where(this%buff1 < zero)
       this%buff1 = zero
     end where
     
     call this%gaussFilter3D(this%buf1)
     call this%gaussFilter3D(this%buf2)
     this%cmodelE_local = this%buf1/(this%buf2 + 1.0D-14)

     call transpose_x_to_y(this%cmodelE_local, this%rbuffyE(:,:,:,1), this%gpE)
     call transpose_y_to_z(this%rbuffyE(:,:,:,1), this%rbuffzE(:,:,:,1), this%gpE)
     this%rbuffzC(:,:,1:this%gpC%zsz(3),1) = 0.5d0*(this%rbuffzE(:,:,1:this%gpC%zsz(3),1)+this%rbuffzE(:,:,2:this%gpC%zsz(3)+1,1))
     call transpose_z_to_y(this%rbuffzC(:,:,:,1), this%rbuffyC(:,:,:,1), this%gpC)
     call transpose_y_to_x(this%rbuffyC(:,:,:,1), cmodelC_local, this%gpC)

   end select
end subroutine


subroutine gaussFilter3D(this, fin)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)), intent(inout) :: fin

   call this%gfiltx%filter1(fin,                   this%rbuffxE(:,:,:,1), this%gpE%xsz(2), this%gpE%xsz(3), 0, 0)
   call transpose_x_to_y   (this%rbuffxE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
   call this%gfilty%filter2(this%rbuffyE(:,:,:,1), this%rbuffyE(:,:,:,2), this%gpE%ysz(1), this%gpE%ysz(3), 0, 0)
   call transpose_y_to_z   (this%rbuffyE(:,:,:,2), this%rbuffzE(:,:,:,1), this%gpE)
   call this%gfiltz%filter3(this%rbuffzE(:,:,:,1), this%rbuffzE(:,:,:,2), this%gpE%zsz(1), this%gpE%zsz(2), 0, 0)

   call transpose_z_to_y   (this%rbuffzE(:,:,:,2), this%rbuffyE(:,:,:,1), this%gpE)
   call transpose_y_to_x   (this%rbuffyE(:,:,:,1), fin,                   this%gpE)
end subroutine 


subroutine planarAverageAndInterpolateToCells(this, numE, denE, ratC)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)), intent(inout) :: numE
   real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)), intent(in)    :: denE
   real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)), intent(out)   :: ratC
   integer :: idx
   
   call transpose_x_to_y(numE, this%rbuffyE(:,:,:,1), this%gpE)
   call transpose_y_to_z(this%rbuffyE(:,:,:,1), this%rbuffzE(:,:,:,1), this%gpE)
   do idx = 1,this%gpE%zsz(3)
      this%rbuffzE(:,:,idx,1) = max(p_sum(sum(this%rbuffzE(:,:,idx,1)))*this%meanfact, zero)
   end do 
   
   call transpose_x_to_y(denE, this%rbuffyE(:,:,:,1), this%gpE)
   call transpose_y_to_z(this%rbuffyE(:,:,:,1), this%rbuffzE(:,:,:,2), this%gpE)
   do idx = 1,this%gpE%zsz(3)
      this%rbuffzE(:,:,idx,2) = p_sum(sum(this%rbuffzE(:,:,idx,2)))*this%meanfact
   end do
   this%rbuffzE(:,:,:,1) = this%rbuffzE(:,:,:,1)/(this%rbuffzE(:,:,:,2) + 1.d-14)
   this%cmodel_allZ = this%rbuffzE(1,1,:,1)

   this%rbuffzC(:,:,1:this%gpC%zsz(3),1) = 0.5d0*(this%rbuffzE(:,:,1:this%gpC%zsz(3),1)+this%rbuffzE(:,:,2:this%gpC%zsz(3)+1,1))
   call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
   call transpose_y_to_x(this%rbuffyE(:,:,:,1), numE, this%gpE)

   call transpose_z_to_y(this%rbuffzC(:,:,:,1), this%rbuffyC(:,:,:,1), this%gpC)
   call transpose_y_to_x(this%rbuffyC(:,:,:,1), ratC, this%gpC)
end subroutine


