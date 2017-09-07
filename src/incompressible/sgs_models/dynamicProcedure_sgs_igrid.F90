subroutine allocateMemory_DynamicProcedure(this, computeFbody, dx, dy, dz)
   class(sgs_igrid), intent(inout), target :: this
   logical, intent(out) :: computeFbody
   real(rkind), intent(in) :: dx, dy, dz
   integer :: ierr

   this%useDynamicProcedure = .true.
   select case( this%dynamicProcedureType) 
   case (1)
      !this%useCglobal = .false.
      allocate(this%alphaij_Filt(this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),9))
      allocate(this%Mij         (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),6))
      allocate(this%Lij         (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),6))
      allocate(this%ui_Filt     (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),3))
      allocate(this%Sij_Filt    (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),6))
      allocate(this%Dsgs_Filt   (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3)))
      allocate(this%buff1       (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3)))
      allocate(this%buff2       (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3)))

      select case(this%AverageType)
      case(0) ! Volume 
         this%modelConstType = 0
         if(this%isPeriodic) then
           call this%spectC%init_TestFilter(this%deltaRat, dx, dy, dz)
         else
           call this%spectE%init_TestFilter(this%deltaRat, dx, dy, dz)
         endif
      case(1) ! Planar
         if(this%isPeriodic) then
           call GracefulExit("You cannot use planar averaging if the problem is periodic in Z",12)
         endif
         this%modelConstType = 1
         allocate(this%cmodel_allZ(this%gpE%zsz(3)))
         allocate(this%cmodelC(this%gpC%xsz(3)))
         allocate(this%cmodelE(this%gpE%xsz(3)))
         call this%spectE%init_TestFilter(this%deltaRat, dx, dy, dz)
      case(2)
         this%modelConstType = 2
         allocate(this%cmodelC_local(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
         allocate(this%cmodelE_local(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)))
         call this%gfiltx%init(this%gpE%xsz(1), .true.)
         call this%gfilty%init(this%gpE%ysz(2), .true.)
         call this%gfiltz%init(this%gpE%zsz(3), this%isPeriodic)
         if(this%isPeriodic) then
           call this%spectC%init_TestFilter(this%deltaRat, dx, dy, dz)
         else
           call this%spectE%init_TestFilter(this%deltaRat, dx, dy, dz)
         endif
      end select
   case (2)
      !this%useCglobal = .true.
      this%modelConstType = 0
      allocate(this%Sij_Filt    (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),6))
      allocate(this%alphaij_Filt(this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),9))
      allocate(this%tauijWM_Filt(this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),2))
      allocate(this%fiE         (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),3))
      allocate(this%fi_Filt     (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),3))
      allocate(this%ui_Filt     (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3),3))
      allocate(this%Dsgs_Filt   (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3)))
      allocate(this%buff1       (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3)))
      allocate(this%buff2       (this%gpE%xsz(1)   ,this%gpE%xsz(2)   ,this%gpE%xsz(3)))
   end select

   this%Dsgs => this%nu_SGS_E
   allocate(this%Tfilhat     (this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)))
   allocate(this%Tfilhatz1   (this%sp_gpE%zsz(1),this%sp_gpE%zsz(2),this%sp_gpE%zsz(3)))
   allocate(this%Tfilhatz2   (this%sp_gpE%zsz(1),this%sp_gpE%zsz(2),this%sp_gpE%zsz(3)))

   if((.not. associated(this%fxC)) .and. (this%DynamicProcedureType==2)) then
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
   if (allocated(this%fiE         )) deallocate(this%fiE         )
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
   if(this%isPeriodic) then
        ! transpose to z
        call transpose_y_to_z(this%Tfilhat, this%cbuffzE(:,:,:,1), this%sp_gpE)
        
        ! use spectC to filter :: the filter call should use spectC and not spectE
        call this%spectC%testFilter_ip_periodic(this%cbuffzE(:,:,:,1))

        ! transpose back to y
        call transpose_z_to_y(this%cbuffzE(:,:,:,1), this%Tfilhat, this%sp_gpE)
   else
       call this%spectE%testFilter_ip(this%Tfilhat)
   endif
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
   if(this%isPeriodic) then
        ! transpose to z
        call transpose_y_to_z(this%Tfilhat, this%cbuffzE(:,:,:,1), this%sp_gpE)
        
        ! use spectC to filter :: the filter call should use spectC and not spectE
        call this%spectC%testFilter_ip_periodic(this%cbuffzE(:,:,:,1))

        ! transpose back to y
        call transpose_z_to_y(this%cbuffzE(:,:,:,1), this%Tfilhat, this%sp_gpE)
   else
       call this%spectE%testFilter_ip(this%Tfilhat)
   endif
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

   if(this%isPeriodic) then
        ! transpose to z
        call transpose_y_to_z(fhat, this%cbuffzE(:,:,:,1), this%sp_gpE)
        
        ! use spectC to filter :: the filter call should use spectC and not spectE
        call this%spectC%testFilter_ip_periodic(this%cbuffzE(:,:,:,1))

        ! transpose back to y
        call transpose_z_to_y(this%cbuffzE(:,:,:,1), this%Tfilhat, this%sp_gpE)
   else
       call this%spectE%TestFilter_oop(fhat, this%Tfilhat)
   endif
   if (this%useVerticalTfilter) then
      call transpose_y_to_z(this%Tfilhat, this%Tfilhatz1, this%sp_gpE)
      call this%gaussianTestFilterZ%filter3(this%Tfilhatz1, this%Tfilhatz2, this%sp_gpE%zsz(1), this%sp_gpE%zsz(2))
      call transpose_z_to_y(this%Tfilhatz2, this%Tfilhat, this%sp_gpE)
   end if
   call this%spectE%ifft(this%Tfilhat, f)

end subroutine
