subroutine allocateMemory_DynamicProcedure(this, computeFbody)
   class(sgs_igrid), intent(inout), target :: this
   logical, intent(out) :: computeFbody
   integer :: ierr

   this%useDynamicProcedure = .true.
   select case( this%dynamicProcedureType) 
   case (1)
      allocate(this%alphaij_Filt(this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),9))
      allocate(this%Mij         (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),6))
      allocate(this%Lij         (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),6))
      allocate(this%ui_Filt     (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),3))
      allocate(this%Sij_Filt    (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),6))
      allocate(this%Dsgs_Filt   (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3)))
      allocate(this%buff1       (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3)))
      allocate(this%buff2       (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3)))
   case (2)
      allocate(this%Sij_Filt    (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),6))
      allocate(this%alphaij_Filt(this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),9))
      allocate(this%tauijWM_Filt(this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),2))
      allocate(this%fi_Filt     (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),3))
      allocate(this%ui_Filt     (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),3))
      allocate(this%Dsgs_Filt   (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3)))
      allocate(this%buff1       (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3)))
      allocate(this%buff2       (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3)))
      allocate(this%Tfilhat     (this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)))
      allocate(this%Tfilhatz1   (this%sp_gpE%zsz(1),this%sp_gpE%zsz(2),this%sp_gpE%zsz(3)))
      allocate(this%Tfilhatz2   (this%sp_gpE%zsz(1),this%sp_gpE%zsz(2),this%sp_gpE%zsz(3)))
   end select

   this%Dsgs => this%nu_SGS_E

   if((.not. associated(this%fi)) .and. (this%DynamicProcedureType==2)) then
      call GracefulExit("Global dynamic procedure needs body force term",1111)
   endif

   this%deltaRat = 2.d0

   if (this%DynamicProcedureType == 2) then
      computeFbody = .true. ! Need body forces for global dynamic procedure
   else
      computeFbody = .false. 
   end if
   ierr = this%gaussianTestFilterZ%init(this%gpE%zsz(3), .false.)
end subroutine

subroutine destroyMemory_DynamicProcedure(this)
   class(sgs_igrid), intent(inout) :: this
     
   nullify(this%Dsgs)
   if (allocated(this%Lij         )) deallocate(this%Lij         )
   if (allocated(this%Mij         )) deallocate(this%Mij         )
   if (allocated(this%Sij_Filt    )) deallocate(this%Sij_Filt    )
   if (allocated(this%alphaij_Filt)) deallocate(this%alphaij_Filt)
   if (allocated(this%tauijWM_Filt)) deallocate(this%tauijWM_Filt)
   if (allocated(this%fi_Filt     )) deallocate(this%fi_Filt     )
   if (allocated(this%ui_Filt     )) deallocate(this%ui_Filt     )
   if (allocated(this%Dsgs_Filt   )) deallocate(this%Dsgs_Filt   )
   if (allocated(this%buff1       )) deallocate(this%buff1       )
   if (allocated(this%buff2       )) deallocate(this%buff2       )
end subroutine


subroutine applyDynamicProcedure(this, uE, vE, wE, uhatE, vhatE, whatE, duidxjE, duidxjEhat)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)), intent(in) :: uE, vE, wE
   complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3))  , intent(in) :: uhatE, vhatE, whatE
   complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3),9), intent(in) :: duidxjEhat
   real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),9), intent(in) :: duidxjE

   select case(this%dynamicProcedureType) 
   case (1) ! Standard (planar averaged) dynamic procedure
      call this%DoStandardDynamicProcedure(uE, vE, wE, uhatE, vhatE, whatE, duidxjEhat)
   case (2) ! Global Dynamic Procedure
      call this%DoGlobalDynamicProcedure(uhatE, vhatE, whatE, uE, vE, wE, duidxjEhat, duidxjE) ! Pass in the relevant stuff, finish the procedure implementation
   end select 

end subroutine

subroutine planarAverage(this, f)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)), intent(inout) :: f
   integer :: idx

   call transpose_x_to_y(f, this%rbuffyE(:,:,:,1), this%gpE)
   call transpose_y_to_z(this%rbuffyE(:,:,:,1), this%rbuffzE(:,:,:,1), this%gpE)
   do idx = 1,this%gpE%zsz(3)
      this%rbuffzE(:,:,idx,1) = p_sum(sum(this%rbuffzE(:,:,idx,1)))*this%meanfact
   end do 
   call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
   call transpose_y_to_x(this%rbuffyE(:,:,:,1), f, this%gpE)

end subroutine

subroutine TestFilter_Real_to_Real(this, f,ffil)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)), intent(out)  :: ffil
   
   call this%spectE%fft(f, this%Tfilhat)
   call this%spectE%testFilter_ip(this%Tfilhat)
   if (this%useVerticalTfilter) then
      call transpose_y_to_z(this%Tfilhat, this%Tfilhatz1, this%sp_gpE)
      call this%gaussianTestFilterZ%filter3(this%Tfilhatz1, this%Tfilhatz2, this%sp_gpE%zsz(1), this%sp_gpE%zsz(2))
      call transpose_z_to_y(this%Tfilhatz2, this%Tfilhat, this%sp_gpE)
   end if
   call this%spectE%ifft(this%Tfilhat, ffil)

end subroutine

subroutine TestFilter_Real_to_Real_ip(this, f)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)), intent(inout)  :: f
   
   call this%spectE%fft(f, this%Tfilhat)
   call this%spectE%testFilter_ip(this%Tfilhat)
   if (this%useVerticalTfilter) then
      call transpose_y_to_z(this%Tfilhat, this%Tfilhatz1, this%sp_gpE)
      call this%gaussianTestFilterZ%filter3(this%Tfilhatz1, this%Tfilhatz2, this%sp_gpE%zsz(1), this%sp_gpE%zsz(2))
      call transpose_z_to_y(this%Tfilhatz2, this%Tfilhat, this%sp_gpE)
   end if
   call this%spectE%ifft(this%Tfilhat, f)

end subroutine

subroutine TestFilter_Cmplx_to_Real(this, fhat, f)
   class(sgs_igrid), intent(inout) :: this
   complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(1),this%sp_gpE%ysz(1)), intent(in) :: fhat
   real(rkind), dimension(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)), intent(out)  :: f

   call this%spectE%TestFilter_oop(fhat, this%Tfilhat)
   if (this%useVerticalTfilter) then
      call transpose_y_to_z(this%Tfilhat, this%Tfilhatz1, this%sp_gpE)
      call this%gaussianTestFilterZ%filter3(this%Tfilhatz1, this%Tfilhatz2, this%sp_gpE%zsz(1), this%sp_gpE%zsz(2))
      call transpose_z_to_y(this%Tfilhatz2, this%Tfilhat, this%sp_gpE)
   end if
   call this%spectE%ifft(this%Tfilhat, f)

end subroutine
