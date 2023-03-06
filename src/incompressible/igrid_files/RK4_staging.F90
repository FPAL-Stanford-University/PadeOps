subroutine advance_RK4_all_Stages(this,dtForced)
  class(igrid), intent(inout), target :: this
  real(rkind), intent(in), optional :: dtforced
  if(present(dtforced)) then
    this%dt = dtforced
  else
    ! Step 0: Compute TimeStep 
    call this%compute_deltaT
  endif

  call this%advance_RK4_Stage(stage=1,dtforced=this%dt)
  call this%advance_RK4_Stage(stage=2)
  call this%advance_RK4_Stage(stage=3)
  call this%advance_RK4_Stage(stage=4)
  
  ! Wrap up this time step 
  call this%wrapup_timestep() 

end subroutine

subroutine advance_RK4_Stage(this, stage, dtforced)
  class(igrid), intent(inout), target :: this
  integer, intent(in) :: stage
  real(rkind), intent(in), optional :: dtforced
  real(rkind), dimension(4) :: a, b
  if (stage == 1) then
    if(present(dtforced)) then
      this%dt = dtforced
    else
      ! Step 0: Compute TimeStep 
      call this%compute_deltaT
    endif

    this%du = im0
    this%dv = im0
    this%dw = im0
  end if
  
  a = [0.d0, 0.5d0, 0.5d0, 1.d0]
  b = [1.d0, 2.d0, 2.d0, 1.d0]/6.d0

  this%ustar = this%uhat + a(stage)*this%dt*this%u_rhs
  this%vstar = this%vhat + a(stage)*this%dt*this%v_rhs
  this%wstar = this%what + a(stage)*this%dt*this%w_rhs
  
  this%uhat => this%ustar
  this%vhat => this%vstar
  this%what => this%wstar
  
  call this%project_and_prep(.false.)
  call this%populate_rhs()
  call this%reset_pointers()

  this%du = this%du + b(stage)*this%dt*this%u_rhs
  this%dv = this%dv + b(stage)*this%dt*this%v_rhs
  this%dw = this%dw + b(stage)*this%dt*this%w_rhs

  if (stage == 4) then
    this%uhat = this%uhat + this%du 
    this%vhat = this%vhat + this%dv 
    this%what = this%what + this%dw 
    
    call this%project_and_prep(.false.)
  end if
end subroutine
