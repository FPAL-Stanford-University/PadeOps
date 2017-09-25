pure function get_GlobalConstant(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
   
   val = this%cmodel_global
end function

pure function get_Ustar(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
  
   val = this%ustar
end function

pure function get_InvObLength(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
   
   val = this%InvObLength
end function

pure function get_umean(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
  
   val = this%umn

end function

pure function get_vmean(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
  
   val = this%vmn

end function

pure function get_uspeedmean(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
  
   val = this%uspmn

end function

pure function get_wTh_surf(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val

   val = this%wTh_surf 

end function

pure function get_dynamicProcedureType(this) result(val)
   class(sgs_igrid), intent(in) :: this
   integer                     :: val
  
   val = this%DynamicProcedureType 

end function
