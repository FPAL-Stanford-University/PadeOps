
   subroutine advance_SSP_RK45_all_stages(this, dtforced)
       class(igrid), intent(inout), target :: this
       real(rkind), intent(in), optional :: dtforced
       if(present(dtforced)) then
         this%dt = dtforced
       else
         ! Step 0: Compute TimeStep 
         call this%compute_deltaT
       endif
       ! Advance the first step
       call this%advance_SSP_RK45_Stage_1()
       call this%advance_SSP_RK45_Stage_2()
       call this%advance_SSP_RK45_Stage_3()
       call this%advance_SSP_RK45_Stage_4()
       call this%advance_SSP_RK45_Stage_5()

       ! Wrap up this time step
       call this%wrapup_timestep()

   end subroutine


   subroutine advance_SSP_RK45_Stage_1(this)
       use RK_coefficients
       class(igrid), intent(inout), target :: this
       integer :: idx
        
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! STAGE 1 !!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! First stage - everything is where it's supposed to be
       if (this%AlreadyHaveRHS) then
           this%AlreadyHaveRHS = .false.
       else
           call this%populate_rhs()
       end if
       this%newTimeStep = .false.
       ! Set the pointers 
       this%uhat1 => this%uExtra(:,:,:,1); this%vhat1 => this%vExtra(:,:,:,1);
       this%what1 => this%wExtra(:,:,:,1); this%That1 => this%TExtra(:,:,:,1);
       ! Do the time step
       this%uhat1 = this%uhat + b01*this%dt*this%u_rhs 
       this%vhat1 = this%vhat + b01*this%dt*this%v_rhs 
       this%what1 = this%what + b01*this%dt*this%w_rhs 
       if (this%isStratified .or. this%initspinup) this%That1 = this%That + b01*this%dt*this%T_rhs
       if (this%useScalars) then
          do idx = 1,this%n_scalars
             this%scalars(idx)%Fhat1 = this%scalars(idx)%Fhat + b01*this%dt*this%scalars(idx)%rhs
          end do 
       end if
       
       ! Now set pointers so that things operate on uhat1, vhat1, etc.
       this%uhat => this%uExtra(:,:,:,1); this%vhat => this%vExtra(:,:,:,1)
       this%what => this%wExtra(:,:,:,1); this%That => this%TExtra(:,:,:,1)
       if (this%useScalars) then
          do idx = 1,this%n_scalars
              this%scalars(idx)%Fhat => this%scalars(idx)%Sfields(:,:,:,2) ! Fhat1 is the second index of Sfields
          end do 
       end if

       ! Now perform the projection and prep for next stage
       
       call this%project_and_prep(this%fastCalcPressure)

   end subroutine

   subroutine advance_SSP_RK45_Stage_2(this)
       use RK_coefficients
       class(igrid), intent(inout), target :: this
       integer :: idx

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! STAGE 2 !!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call this%populate_rhs()
       ! reset u, v, w pointers
       call this%reset_pointers()
       ! Set the new pointers
       this%uhat2 => this%uExtra(:,:,:,2); this%vhat2 => this%vExtra(:,:,:,2)
       this%what2 => this%wExtra(:,:,:,2); this%That2 => this%TExtra(:,:,:,2)
       ! Do the time step
       this%uhat2 = a20*this%uhat + a21*this%uhat1 + b12*this%dt*this%u_rhs
       this%vhat2 = a20*this%vhat + a21*this%vhat1 + b12*this%dt*this%v_rhs
       this%what2 = a20*this%what + a21*this%what1 + b12*this%dt*this%w_rhs
       if (this%isStratified .or. this%initspinup) this%That2 = a20*this%That + a21*this%That1 + b12*this%dt*this%T_rhs
       if (this%useScalars) then
          do idx = 1,this%n_scalars
             this%scalars(idx)%Fhat2 = a20*this%scalars(idx)%Fhat + a21*this%scalars(idx)%Fhat1 + b12*this%dt*this%scalars(idx)%rhs
          end do 
       end if
       ! now set the u, v, w, pointers to u2, v2, w2
       this%uhat => this%uExtra(:,:,:,2); this%vhat => this%vExtra(:,:,:,2)
       this%what => this%wExtra(:,:,:,2); this%That => this%TExtra(:,:,:,2)
       if (this%useScalars) then
          do idx = 1,this%n_scalars
           this%scalars(idx)%Fhat => this%scalars(idx)%Sfields(:,:,:,3) ! Fhat2 is the third index of Sfields
          end do 
       end if
       ! Now perform the projection and prep for next stage
       call this%project_and_prep(.false.)

   end subroutine


   subroutine advance_SSP_RK45_Stage_3(this)
       use RK_coefficients
       class(igrid), intent(inout), target :: this
       integer :: idx


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! STAGE 3 !!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call this%populate_rhs()
       ! reset u, v, w pointers
       call this%reset_pointers()
       ! Set the pointers 
       this%uhat3 => this%uExtra(:,:,:,3); this%vhat3 => this%vExtra(:,:,:,3)
       this%what3 => this%wExtra(:,:,:,3); this%That3 => this%TExtra(:,:,:,3)
       ! Do the time step 
       this%uhat3 = a30*this%uhat + a32*this%uhat2 + b23*this%dt*this%u_rhs
       this%vhat3 = a30*this%vhat + a32*this%vhat2 + b23*this%dt*this%v_rhs
       this%what3 = a30*this%what + a32*this%what2 + b23*this%dt*this%w_rhs
       if (this%isStratified .or. this%initspinup) this%That3 = a30*this%That + a32*this%That2 + b23*this%dt*this%T_rhs
       if (this%useScalars) then
          do idx = 1,this%n_scalars
             this%scalars(idx)%Fhat3 = a30*this%scalars(idx)%Fhat + a32*this%scalars(idx)%Fhat2 + b23*this%dt*this%scalars(idx)%rhs
          end do 
       end if
       ! now set u, v, w pointers to point to u3, v3, w3
       this%uhat => this%uExtra(:,:,:,3); this%vhat => this%vExtra(:,:,:,3)
       this%what => this%wExtra(:,:,:,3); this%That => this%TExtra(:,:,:,3)
       if (this%useScalars) then
          do idx = 1,this%n_scalars
           this%scalars(idx)%Fhat => this%scalars(idx)%Sfields(:,:,:,4) ! Fhat3 is the third index of Sfields
          end do 
       end if
       ! Now perform the projection and prep for next time step
       call this%project_and_prep(.false.)

   end subroutine

   
   subroutine advance_SSP_RK45_Stage_4(this)
       use RK_coefficients
       class(igrid), intent(inout), target :: this
       integer :: idx

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! STAGE 4 !!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call this%populate_rhs()
       ! reset u, v, w pointers
       call this%reset_pointers()
       ! Set the pointers 
       !this%uhat4 => this%uhat; this%vhat4 => this%vhat
       !this%what4 => this%what; this%That4 => this%That
       ! Do the time step 
       this%uhat = a40*this%uhat + a43*this%uhat3 + b34*this%dt*this%u_rhs
       this%vhat = a40*this%vhat + a43*this%vhat3 + b34*this%dt*this%v_rhs
       this%what = a40*this%what + a43*this%what3 + b34*this%dt*this%w_rhs
       if (this%isStratified .or. this%initspinup) this%That = a40*this%That + a43*this%That3 + b34*this%dt*this%T_rhs
       if (this%useScalars) then
          do idx = 1,this%n_scalars
             this%scalars(idx)%Fhat = a40*this%scalars(idx)%Fhat + a43*this%scalars(idx)%Fhat3 + b34*this%dt*this%scalars(idx)%rhs
          end do 
       end if
       ! < IMPORTANT: no need to do anything here since uhat4 is really just uhat > 
       ! Now perform the projection and prep for next time step
       call this%project_and_prep(.false.)

   end subroutine


   subroutine advance_SSP_RK45_Stage_5(this)
       use RK_coefficients
       class(igrid), intent(inout), target :: this
       integer :: idx


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!! STAGE 5 !!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! First reset urhs pointers to point to spare buffers
       this%u_rhs => this%uRHSExtra(:,:,:,1); this%v_rhs => this%vRHSExtra(:,:,:,1);
       this%w_rhs => this%wRHSExtra(:,:,:,1); this%T_rhs => this%TRHSExtra(:,:,:,1);
       if (this%useScalars) then
          do idx = 1,this%n_scalars
           this%scalars(idx)%rhs => this%scalars(idx)%rhs_storage(:,:,:,2) 
          end do 
       end if
       call this%populate_rhs()
       ! Reset pointers
       call this%reset_pointers(resetRHS=.true.)
       ! Do the time step 
       this%uhat = a52*this%uhat2 + a53*this%uhat3 + b35*this%dt*this%u_rhs + a54*this%uhat + b45*this%dt*this%uRHSExtra(:,:,:,1)
       this%vhat = a52*this%vhat2 + a53*this%vhat3 + b35*this%dt*this%v_rhs + a54*this%vhat + b45*this%dt*this%vRHSExtra(:,:,:,1)
       this%what = a52*this%what2 + a53*this%what3 + b35*this%dt*this%w_rhs + a54*this%what + b45*this%dt*this%wRHSExtra(:,:,:,1)
       if (this%isStratified .or. this%initspinup) this%That = a52*this%That2 + a53*this%That3 + b35*this%dt*this%T_rhs + a54*this%That + b45*this%dt*this%TRHSExtra(:,:,:,1) 
       if (this%useScalars) then
          do idx = 1,this%n_scalars
             this%scalars(idx)%Fhat = a52*this%scalars(idx)%Fhat2 + a53*this%scalars(idx)%Fhat3  &
                                  & + b35*this%dt*this%scalars(idx)%rhs + a54*this%scalars(idx)%Fhat &
                                  & + b45*this%dt*this%scalars(idx)%rhs_storage(:,:,:,2)
          end do 
       end if
       ! Now perform the projection and prep for next time step
       call this%project_and_prep(.false.)
       
   end subroutine
