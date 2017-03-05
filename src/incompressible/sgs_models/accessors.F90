pure function get_GlobalConstant(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
   
   val = this%cmodel_global
end function

pure function get_Ustar(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
  
   if (this%useWallModel) then
      val = this%ustar
   else
      val = 1.d0
   end if
end function

pure function get_InvObLength(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
   
   val = this%InvObLength
end function

pure function get_umean(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
  
   if (this%useWallModel) then
      val = this%umn
   else
      val = 0.d0
   end if

end function

pure function get_vmean(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
  
   if (this%useWallModel) then
      val = this%vmn
   else
      val = 0.d0
   end if

end function

pure function get_uspeedmean(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
  
   if (this%useWallModel) then
      val = this%uspmn
   else
      val = 0.d0
   end if

end function
