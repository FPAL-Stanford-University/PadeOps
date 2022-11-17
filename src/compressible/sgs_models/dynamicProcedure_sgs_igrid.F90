subroutine allocateMemory_DynamicProcedure(this, computeFbody, deltaRat)
   class(sgs_igrid), intent(inout), target :: this
   logical, intent(out) :: computeFbody
   real(rkind), intent(in) :: deltaRat 
   integer :: ierr

   if (this%mid .ne. 0) then
      call gracefulExit("Dynamic procedure only supported for smagorinsky",312)
   end if

   this%usingDynamicProcedureMomentum = .true.
   allocate(this%LambdaDynProc_C(this%gpC%xsz(1) ,this%gpC%xsz(2) ,this%gpC%xsz(3)))
   
   if (this%UseDynamicProcedureScalar) then
      allocate(this%BetaDynProc_C(this%gpC%xsz(1) ,this%gpC%xsz(2) ,this%gpC%xsz(3)))
   end if 

   allocate(this%Beta_1d(this%gpC%zsz(3)))
   allocate(this%Lambda_1d(this%gpC%zsz(3)))
   ierr = this%gaussianTestFilterZ%init(this%gpC%zsz(3), .false.)
   
   this%deltaRat = deltaRat
   call this%spectC%InitTestFilter(this%deltaRat)
   allocate(this%Dsgs     (this%gpC%xsz(1)   ,this%gpC%xsz(2)   ,this%gpC%xsz(3)))
   allocate(this%Dsgs_Filt(this%gpC%xsz(1)   ,this%gpC%xsz(2)   ,this%gpC%xsz(3)))
   allocate(this%Lij      (this%gpC%xsz(1)   ,this%gpC%xsz(2)   ,this%gpC%xsz(3),6))
   allocate(this%Sij_Filt (this%gpC%xsz(1)   ,this%gpC%xsz(2)   ,this%gpC%xsz(3),6))
   allocate(this%ui_Filt  (this%gpC%xsz(1)   ,this%gpC%xsz(2)   ,this%gpC%xsz(3),3))
   
   allocate(this%Tfilhat(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))

   computeFbody = .false. 
end subroutine

subroutine destroyMemory_DynamicProcedure(this)
   class(sgs_igrid), intent(inout) :: this
   
   if (allocated(this%Tfilhat)) deallocate(This%Tfilhat)

   if (allocated(this%Lij            )) deallocate(this%Lij         )
   if (allocated(this%Sij_Filt       )) deallocate(this%Sij_Filt    )
   if (allocated(this%ui_Filt        )) deallocate(this%ui_Filt     )
   if (allocated(this%Dsgs_Filt      )) deallocate(this%Dsgs_Filt   )
   if (allocated(this%Dsgs           )) deallocate(this%Dsgs        )
   if (allocated(this%LambdaDynProc_C)) deallocate(this%LambdaDynProc_C   )
   if (allocated(this%BetaDynProc_C  )) deallocate(this%BetaDynProc_C   )
end subroutine


!subroutine applyDynamicProcedure(this, uE, vE, wE, uhatE, vhatE, whatE, duidxjE, duidxjEhat)
!   class(sgs_igrid), intent(inout) :: this
!   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: uE, vE, wE
!   complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3))  , intent(in) :: uhatE, vhatE, whatE
!   complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3),9), intent(in) :: duidxjEhat
!   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),9), intent(in) :: duidxjE
!
!   select case(this%dynamicProcedureType) 
!   case (1) ! Standard (planar averaged) dynamic procedure
!      call this%DoStandardDynamicProcedure(uE, vE, wE, uhatE, vhatE, whatE, duidxjEhat)
!   case (2) ! Global Dynamic Procedure
!      call this%DoGlobalDynamicProcedure(uhatE, vhatE, whatE, uE, vE, wE, duidxjEhat, duidxjE) ! Pass in the relevant stuff, finish the procedure implementation
!   end select 
!
!end subroutine

!subroutine planarAverage(this, f)
!   class(sgs_igrid), intent(inout) :: this
!   real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)), intent(inout) :: f
!   integer :: idx
!
!   call transpose_x_to_y(f, this%rbuffyE(:,:,:,1), this%gpC)
!   call transpose_y_to_z(this%rbuffyE(:,:,:,1), this%rbuffzE(:,:,:,1), this%gpC)
!   do idx = 1,this%gpC%zsz(3)
!      this%rbuffzE(:,:,idx,1) = p_sum(sum(this%rbuffzE(:,:,idx,1)))*this%meanfact
!   end do 
!   call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpC)
!   call transpose_y_to_x(this%rbuffyE(:,:,:,1), f, this%gpC)
!
!end subroutine

