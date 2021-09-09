pure function get_GlobalConstant(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
   
   val = this%cmodel_global
end function

pure function get_z0(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
   
   val = this%z0
end function

pure function get_Ustar(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
  
   val = this%ustar
end function

pure function get_Ustar_upstream(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
  
   val = this%ustar_upstream
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

pure function get_Tmean(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
  
   val = this%Tmn

end function

pure function get_Tsurf(this) result(val)
   class(sgs_igrid), intent(in) :: this
   real(rkind)                     :: val
  
   val = this%Tsurf

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

function get_uw_surf(this) result(val)
   class(sgs_igrid), intent(inout) :: this
   real(rkind)                     :: val
  
   call transpose_x_to_y(this%tau_13,this%rbuffyE(:,:,:,1),this%gpE)
   call transpose_y_to_z(this%rbuffyE(:,:,:,1),this%rbuffzE(:,:,:,1),this%gpE)
   call this%PadeDer%interpz_E2C(this%rbuffzE(:,:,:,1),this%rbuffzC(:,:,:,1),0,0)
   call transpose_z_to_y(this%rbuffzC(:,:,:,1),this%rbuffyC(:,:,:,1),this%gpC)
   call transpose_y_to_x(this%rbuffyC(:,:,:,1),this%tau_13C,this%gpC)
   call this%spectC%fft(this%tau_13C,this%cbuffyC(:,:,:,1))
   this%uw_surf = real(this%cbuffyC(1,1,1,1),rkind)*this%meanFact

   val = this%uw_surf
end function

function get_vw_surf(this) result(val)
   class(sgs_igrid), intent(inout) :: this
   real(rkind)                     :: val
  
   call transpose_x_to_y(this%tau_23,this%rbuffyE(:,:,:,1),this%gpE)
   call transpose_y_to_z(this%rbuffyE(:,:,:,1),this%rbuffzE(:,:,:,1),this%gpE)
   call this%PadeDer%interpz_E2C(this%rbuffzE(:,:,:,1),this%rbuffzC(:,:,:,1),0,0)
   call transpose_z_to_y(this%rbuffzC(:,:,:,1),this%rbuffyC(:,:,:,1),this%gpC)
   call transpose_y_to_x(this%rbuffyC(:,:,:,1),this%tau_23C,this%gpC)
   call this%spectC%fft(this%tau_23C,this%cbuffyC(:,:,:,1))
   this%vw_surf = real(this%cbuffyC(1,1,1,1),rkind)*this%meanFact

   val = this%vw_surf
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
    use mpi
    use kind_parameters, only: mpirkind
    class(sgs_igrid), intent(inout) :: this

    real(rkind) :: uwloc(2), uwglob(2)
    integer :: ierr

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

    ! set uw_surf and vw_surf
    uwloc(1:2) = zero
    if(this%gpC%xst(3)==1) then
        uwloc(1) = sum(this%tau_13C(:,:,1))
        uwloc(2) = sum(this%tau_23C(:,:,1))
    endif
    call mpi_reduce(uwloc, uwglob, 2, mpirkind, MPI_SUM, 0, mpi_comm_world, ierr)
    this%uw_surf = uwglob(1)*this%meanfact
    this%vw_surf = uwglob(2)*this%meanfact
end subroutine 

subroutine get_z0varstats(this, fout)
    use mpi
    use kind_parameters, only: mpirkind
    class(sgs_igrid), intent(in) :: this
    real(rkind), dimension(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)) :: fout

    fout(:,:,1) = this%ustarsqvar
    fout(:,:,2) = this%Uxvar
    fout(:,:,3) = this%Uyvar
    fout(:,:,4) = this%filteredSpeedSq(:,:,1)
    if(this%WallModel == 4) then
        fout(:,:,5) = this%nlptype
    endif

end subroutine 

pure function get_is_z0_varying(this) result(val)
   class(sgs_igrid), intent(in) :: this
   logical                      :: val
   
   val = this%is_z0_varying
end function

