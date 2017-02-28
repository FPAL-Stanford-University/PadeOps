pure function getGlobalConstant(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
   
   val = this%cmodel_global
end function

pure function getUstar(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
   
   val = this%ustar
end function

pure function getInvObLength(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
   
   val = this%InvObLength
end function
