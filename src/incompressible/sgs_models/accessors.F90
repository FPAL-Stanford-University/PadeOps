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


pure function getMax_DynSmagConst(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind) :: val

   if (this%usingDynamicProcedureMomentum) then
      val = maxval(this%lambda_1d) 
   else
      val = 0.d0
   end if

end function


pure function getMax_DynPrandtl(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind) :: val

   if (this%useDynamicProcedureScalar) then
      val = maxval(this%beta_1d)
      !val = maxval(this%BetaDynProc_C)
   else
      val = 0.d0
   end if

end function

subroutine set_buoyancyFactor(this, buoyancyFact)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), intent(in) :: buoyancyFact

   this%buoyancyFact = buoyancyFact 

end subroutine 

pure function usingDynProc(this) result(val)
   class(sgs_igrid), intent(in) :: this
   logical :: val

   val = this%usingDynamicProcedureMomentum

end function

pure function get_dynamicProcedureType(this) result(val)
   class(sgs_igrid), intent(in) :: this
   integer                     :: val
  
   val = this%DynamicProcedureType 

end function

subroutine populate_tauij_E_to_C(this)
    class(sgs_igrid), intent(inout) :: this

    ! This subroutine interpolates tau_13 and tau_23 to tau_13C and tau_23C (for
    ! post-processing) 

    call transpose_x_to_y(this%tau_13,this%rbuffyE(:,:,:,1),this%gpE)
    call transpose_y_to_z(this%rbuffyE(:,:,:,1),this%rbuffzE(:,:,:,1),this%gpE)
    call this%PadeDer%interpz_E2C(this%rbuffzE(:,:,:,1),this%rbuffzC(:,:,:,1),0,0)
    call transpose_z_to_y(this%rbuffzC(:,:,:,1),this%rbuffyC(:,:,:,1),this%gpC)
    call transpose_y_to_x(this%rbuffyC(:,:,:,1),this%tau_13C,this%gpC)

    call transpose_x_to_y(this%tau_23,this%rbuffyE(:,:,:,1),this%gpE)
    call transpose_y_to_z(this%rbuffyE(:,:,:,1),this%rbuffzE(:,:,:,1),this%gpE)
    call this%PadeDer%interpz_E2C(this%rbuffzE(:,:,:,1),this%rbuffzC(:,:,:,1),0,0)
    call transpose_z_to_y(this%rbuffzC(:,:,:,1),this%rbuffyC(:,:,:,1),this%gpC)
    call transpose_y_to_x(this%rbuffyC(:,:,:,1),this%tau_23C,this%gpC)
end subroutine 
