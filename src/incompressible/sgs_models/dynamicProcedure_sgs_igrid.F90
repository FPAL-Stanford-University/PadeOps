subroutine allocateMemory_DynamicProcedure(this)
  class(sgs_igrid), intent(inout) :: this

  select case(this%DynamicProcedureType)
  case (1)
     allocate(this%rbuff_DynProc(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),1))
     allocate(this%cbuff_DynProc(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3),1))
  case (2)
     print*, "Global Dynamic Procedure is Incomplete."
  end select
end subroutine

subroutine destroyMemory_DynamicProcedure(this)
  class(sgs_igrid), intent(inout) :: this
     
  deallocate(this%rbuff_DynProc, this%cbuff_DynProc)
end subroutine
