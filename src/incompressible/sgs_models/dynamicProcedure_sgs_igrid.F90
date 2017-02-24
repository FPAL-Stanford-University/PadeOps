subroutine allocateMemory_DynamicProcedure(this)
  class(sgs_igrid), intent(inout), target :: this

  allocate(this%Sij_Filt    (this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),6))
  allocate(this%alphaij_Filt(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),9))
  allocate(this%tauijWM_Filt(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),2))
  allocate(this%fi_Filt     (this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),3))
  allocate(this%ui_Filt     (this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),3))
  allocate(this%Dsgs_Filt   (this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)))
  allocate(this%buff1       (this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)))
  allocate(this%buff2       (this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)))

  this%Dsgs => this%nu_SGS_E

  if((.not. associated(this%fi)) .and. (this%DynamicProcedureType==2)) then
     call GracefulExit("Global dynamic procedure needs body force term",1111)
  endif

  !select case(this%DynamicProcedureType)
  !case (1)
  !   allocate(this%rbuff_DynProc(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),1))
  !   allocate(this%cbuff_DynProc(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3),1))
  !case (2)
  !   print*, "Global Dynamic Procedure is Incomplete."
  !end select
end subroutine

subroutine destroyMemory_DynamicProcedure(this)
  class(sgs_igrid), intent(inout) :: this
    
  nullify(this%Dsgs) 
  deallocate(this%Sij_Filt)
  deallocate(this%alphaij_Filt)
  deallocate(this%tauijWM_Filt)
  deallocate(this%fi_Filt)
  deallocate(this%ui_Filt)
  deallocate(this%Dsgs_Filt)
  deallocate(this%buff1)
  deallocate(this%buff2)
  !deallocate(this%rbuff_DynProc, this%cbuff_DynProc)
end subroutine
