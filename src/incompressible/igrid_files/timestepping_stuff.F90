   subroutine reset_pointers(this, resetRHS)
       class(igrid), intent(inout), target :: this
       logical, intent(in), optional :: resetRHS 
       integer :: idx 

       this%uhat => this%SfieldsC(:,:,:,1); 
       this%vhat => this%SfieldsC(:,:,:,2); 
       this%what => this%SfieldsE(:,:,:,1)
       if ((this%isStratified) .or. (this%initspinup)) then
           this%That => this%SfieldsC(:,:,:,4)
       end if

       if (present(resetRHS)) then
           if (resetRHS) then
              this%u_rhs => this%rhsC(:,:,:,1); 
              this%v_rhs => this%rhsC(:,:,:,2); 
              this%w_rhs => this%rhsE(:,:,:,1)
              if ((this%isStratified) .or. (this%initspinup)) then
                  this%T_rhs =>  this%rhsC(:,:,:,3)
              end if
           end if 
       end if 

       if (this%useScalars) then
           do idx = 1,this%n_scalars
              if (present(resetRHS)) then
                 call this%scalars(idx)%reset_pointers(resetRHS)
              else
                 call this%scalars(idx)%reset_pointers()
              end if  
           end do 
       end if 

   end subroutine

   
   subroutine TVD_RK3(this, dtforced)
       class(igrid), intent(inout), target :: this
       real(rkind), intent(in), optional :: dtforced
       integer :: idx 

       if(present(dtforced)) then
         this%dt = dtforced
       else
         ! Step 0: Compute TimeStep 
         call this%compute_deltaT
       endif
      
       !!! STAGE 1
       ! First stage - everything is where it's supposed to be
       if (this%AlreadyHaveRHS) then
           this%AlreadyHaveRHS = .false.
       else
           call this%populate_rhs()
       end if
       !print*, sum(abs(this%u_rhs))
       !print*, sum(abs(this%v_rhs))
       !print*, sum(abs(this%w_rhs))
       !print*, sum(abs(this%T_rhs))
     
       this%newTimeStep = .false. 
       this%uhat1 = this%uhat + this%dt*this%u_rhs 
       this%vhat1 = this%vhat + this%dt*this%v_rhs 
       this%what1 = this%what + this%dt*this%w_rhs 
       if (this%isStratified .or. this%initspinup) this%That1 = this%That + this%dt*this%T_rhs
       if (this%useScalars) then
          do idx = 1,this%n_scalars
             this%scalars(idx)%Fhat1 = this%scalars(idx)%Fhat + this%dt*this%scalars(idx)%rhs
          end do 
       end if
       
       ! Now set pointers so that things operate on uhat1, vhat1, etc.
       this%uhat => this%SfieldsC2(:,:,:,1); this%vhat => this%SfieldsC2(:,:,:,2); this%what => this%SfieldsE2(:,:,:,1); 
       if (this%isStratified .or. this%initspinup) this%That => this%SfieldsC2(:,:,:,3)
       if (this%useScalars) then
          do idx = 1,this%n_scalars
              this%scalars(idx)%Fhat => this%scalars(idx)%Sfields(:,:,:,2) 
          end do 
       end if

       ! Now perform the projection and prep for next stage
       call this%project_and_prep(this%fastCalcPressure)
        

       !!! STAGE 2
       ! Second stage - u, v, w are really pointing to u1, v1, w1 (which is
       ! what we want. 
       call this%populate_rhs()

       ! reset u, v, w pointers
       call this%reset_pointers()
       this%uhat1 = (3.d0/4.d0)*this%uhat + (1.d0/4.d0)*this%uhat1 + (1.d0/4.d0)*this%dt*this%u_rhs
       this%vhat1 = (3.d0/4.d0)*this%vhat + (1.d0/4.d0)*this%vhat1 + (1.d0/4.d0)*this%dt*this%v_rhs
       this%what1 = (3.d0/4.d0)*this%what + (1.d0/4.d0)*this%what1 + (1.d0/4.d0)*this%dt*this%w_rhs
       if (this%isStratified .or. this%initspinup) this%That1 = (3.d0/4.d0)*this%That + (1.d0/4.d0)*this%That1 + (1.d0/4.d0)*this%dt*this%T_rhs
       if (this%useScalars) then
          do idx = 1,this%n_scalars
             this%scalars(idx)%Fhat1 = (3.d0/4.d0)*this%scalars(idx)%Fhat + (1.d0/4.d0)*this%scalars(idx)%Fhat1 &
                                   & + (1.d0/4.d0)*this%dt*this%scalars(idx)%rhs
          end do 
       end if
       ! now set the u, v, w, pointers to u1, v1, w1
       this%uhat => this%SfieldsC2(:,:,:,1); this%vhat => this%SfieldsC2(:,:,:,2); this%what => this%SfieldsE2(:,:,:,1); 
       if (this%isStratified .or. this%initspinup) this%That => this%SfieldsC2(:,:,:,3)
       if (this%useScalars) then
          do idx = 1,this%n_scalars
           this%scalars(idx)%Fhat => this%scalars(idx)%Sfields(:,:,:,2) 
          end do 
       end if

       ! Now perform the projection and prep for next stage
       call this%project_and_prep(.false.)


       !!! STAGE 3 (Final Stage)
       ! Third stage - u, v, w are really pointing to u2, v2, w2 (which is what
       ! we really want. 
       call this%populate_rhs()
       
       ! reset u, v, w pointers
       call this%reset_pointers()
       this%uhat = (1.d0/3.d0)*this%uhat + (2.d0/3.d0)*this%uhat1 + (2.d0/3.d0)*this%dt*this%u_rhs
       this%vhat = (1.d0/3.d0)*this%vhat + (2.d0/3.d0)*this%vhat1 + (2.d0/3.d0)*this%dt*this%v_rhs
       this%what = (1.d0/3.d0)*this%what + (2.d0/3.d0)*this%what1 + (2.d0/3.d0)*this%dt*this%w_rhs
       if (this%isStratified .or. this%initspinup) this%That = (1.d0/3.d0)*this%That + (2.d0/3.d0)*this%That1 + (2.d0/3.d0)*this%dt*this%T_rhs
       if (this%useScalars) then
          do idx = 1,this%n_scalars
             this%scalars(idx)%Fhat  = (1.d0/3.d0)*this%scalars(idx)%Fhat + (2.d0/3.d0)*this%scalars(idx)%Fhat1 &
                                   & + (2.d0/3.d0)*this%dt*this%scalars(idx)%rhs
          end do 
       end if
       
       ! Now perform the projection and prep for next time step
       call this%project_and_prep(.false.)


       ! Wrap up this time step 
       call this%wrapup_timestep() 
   
   end subroutine


   subroutine SSP_RK45(this, dtforced)
       class(igrid), intent(inout), target :: this
       real(rkind), intent(in), optional :: dtforced
       real(rkind), parameter :: b01 = 0.39175222657189d0 , b12 = 0.368410593050371d0, b23 = 0.25189177427169d0, b34 = 0.54497475022852d0
       real(rkind), parameter :: b35 = 0.06369246866629d0 , b45 = 0.22600748323690d0
       
       real(rkind), parameter :: a20 = 0.444370493651235d0, a21 = 0.555629506348765d0
       real(rkind), parameter :: a30 = 0.620101851488403d0, a32 = 0.379898148511597d0
       real(rkind), parameter :: a40 = 0.17807995439313d0 , a43 = 0.821920045606868d0
       real(rkind), parameter :: a52 = 0.517231671970585d0, a53 = 0.096059710526147d0, a54 = 0.386708617503269d0

       integer :: idx 

       if(present(dtforced)) then
         this%dt = dtforced
       else
         ! Step 0: Compute TimeStep 
         call this%compute_deltaT
       endif
       
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

       ! Wrap up this time step 
       call this%wrapup_timestep() 

   end subroutine

   subroutine populate_rhs(this, CopyForDNSpress, CopyForFringePress, copyForTurbinePress)
      class(igrid), intent(inout) :: this
      !integer,           intent(in)    :: RKstage
      logical, intent(in), optional :: CopyForDNSpress, CopyForFringePress, copyForTurbinePress
      logical :: copyFringeRHS, copyTurbRHS
      integer :: idx 
       

      if (present(copyForTurbinePress)) then
         copyTurbRHS = copyForTurbinePress
      else
         copyTurbRHS = .false. 
      end if

      ! Step 1: Non Linear Term 
      if (useSkewSymm) then
          call this%addNonLinearTerm_skewSymm()
      else
          call this%AddNonLinearTerm_Rot()
      end if

      ! Step 2: Coriolis Term
      if (this%useCoriolis) then
          call this%AddCoriolisTerm()
      end if 
      
      
      ! Step 3a: Extra Forcing 
      if (this%useExtraForcing) then
          call this%addExtraForcingTerm()
      end if

      ! Step 3b: Wind Turbines
      !if (this%useWindTurbines .and. (RKstage==1)) then
      if (this%useWindTurbines) then
         if (copyTurbRHS) then
          call this%WindTurbineArr%getForceRHS(this%dt, this%u, this%v, this%wC, this%u_rhs, this%v_rhs, this%w_rhs, &
                & this%newTimestep, this%inst_horz_avg_turb, uturb=this%urhs_turbine, vturb=this%vrhs_turbine, wturb=this%wrhs_turbine) 
         else
          !if (allocated(this%inst_horz_avg_turb)) then
              call this%WindTurbineArr%getForceRHS(this%dt, this%u, this%v, this%wC,&
                                   this%u_rhs, this%v_rhs, this%w_rhs, this%newTimestep, this%inst_horz_avg_turb)
          !else
          !    call this%WindTurbineArr%getForceRHS(this%dt, this%u, this%v, this%wC,&
          !                         this%u_rhs, this%v_rhs, this%w_rhs, this%newTimestep)
          !end if
         end if
      end if 
     
      ! Step 4: Buoyance + Sponge (inside Buoyancy)
      if (this%isStratified .or. this%initspinup) then
          call this%addBuoyancyTerm()
      end if 
      
       if (present(CopyForDNSpress)) then
           if (CopyForDNSpress) then
              if (copyTurbRHS) then
                 this%urhs_dns = this%u_rhs - this%urhs_turbine 
                 this%vrhs_dns = this%v_rhs - this%vrhs_turbine
                 this%wrhs_dns = this%w_rhs - this%wrhs_turbine
              else
                 this%urhs_dns = this%u_rhs
                 this%vrhs_dns = this%v_rhs
                 this%wrhs_dns = this%w_rhs
              end if
           end if   
       end if

      ! Step 6: SGS and Viscous Stress Terms
      if (this%useSGS) then
          call this%sgsmodel%getRHS_SGS(this%u_rhs, this%v_rhs, this%w_rhs,      this%duidxjC, this%duidxjE, &
                                        this%uhat,  this%vhat,  this%whatC,      this%That,    this%u,       &
                                        this%v,     this%wC,    this%newTimeStep,this%dTdxC,   this%dTdyC,   & 
                                        this%dTdzC, this%dTdxE, this%dTdyE, this%dTdzE)

          if (this%isStratified .or. this%initspinup) then
             call this%sgsmodel%getRHS_SGS_Scalar(this%T_rhs, this%dTdxC, this%dTdyC, this%dTdzC, this%dTdzE, &
                                        this%u, this%v, this%wC, this%T, this%That, this%duidxjC, this%turbPr)
          end if
   
      else
           if (.not. this%isInviscid) then
               call this%addViscousTerm()
           end if

      end if

       ! Step 7: Fringe source term if fringe is being used (non-periodic)
       if (present(copyForFringePress)) then
           copyFringeRHS = copyForFringePress
       else
           copyFringeRHS = .false. 
       end if
       if (this%usedoublefringex) then
           if (copyFringeRHS) then
               call this%fringe_x1%addFringeRHS(this%dt, this%u_rhs, this%v_rhs, this%w_rhs, this%u, this%v, this%w, &
                                   this%urhs_fringe, this%vrhs_fringe, this%wrhs_fringe, addF=.false. )
               call this%fringe_x2%addFringeRHS(this%dt, this%u_rhs, this%v_rhs, this%w_rhs, this%u, this%v, this%w, &
                                   this%urhs_fringe, this%vrhs_fringe, this%wrhs_fringe, addF=.true. )
           else
               call this%fringe_x1%addFringeRHS(this%dt, this%u_rhs, this%v_rhs, this%w_rhs, this%u, this%v, this%w)
               call this%fringe_x2%addFringeRHS(this%dt, this%u_rhs, this%v_rhs, this%w_rhs, this%u, this%v, this%w)
           end if 
       else
           if (this%useFringe) then
              if (copyFringeRHS) then
                  call this%fringe_x%addFringeRHS(this%dt, this%u_rhs, this%v_rhs, this%w_rhs, this%u, this%v, this%w, & 
                                   this%urhs_fringe, this%vrhs_fringe, this%wrhs_fringe, addF=.false. )
              else
                  call this%fringe_x%addFringeRHS(this%dt, this%u_rhs, this%v_rhs, this%w_rhs, this%u, this%v, this%w)
              end if 
           end if 
       end if 

       ! Step 8: HIT forcing source term
       if (this%useHITForcing) then
           call this%hitforce%getRHS_HITForcing(this%u_rhs, this%v_rhs, this%w_rhs, this%uhat, this%vhat, this%what, this%newTimeStep)
       end if

       ! Step 9: Frame rotatio PI controller to fix yaw angle at a given height
       if (this%useControl .AND. abs(180.d0/pi*this%angleHubHeight)>0.0d0) then
           call this%angCont_yaw%update_RHS_control(this%dt, this%u_rhs, this%v_rhs, &
                         this%w_rhs, this%u, this%v, this%newTimeStep, this%angleHubHeight, this%wFilt, this%deltaGalpha, this%zHubIndex, this%angleTrigger)
           this%totalAngle = this%totalAngle + this%angleHubHeight
           !if (this%newTimeStep)  then
               !this%G_alpha = this%G_alpha - deltaGalpha
               !this%frameAngle = this%frameAngle + deltaGalpha 
               !this%coriolis_omegaX = cos(this%latitude*pi/180.d0)*sin(this%totalAngle)
               !this%coriolis_omegaZ = sin(this%latitude*pi/180.d0)
               !this%coriolis_omegaY = cos(this%latitude*pi/180.d0)*cos(this%totalAngle)
           !end if
                       !this%angleHubHeight = 0.d0
           !do j = 1, this%gpC%xsz(1)
           !    do i = 1, this%gpC%xsz(2)
           !          this%angleHubHeight = this%angleHubHeight + this%rbuffxC(j,i,zHubIndex,1)
           !    enddo
           !enddo
           !this%angleHubHeight = this%angleHubHeight / (float(this%gpC%xsz(1)) * float(this%gpC%xsz(1)))
       end if 

       ! Step 10: Add sponge
       if (this%useSponge) then
           call this%addSponge()
       end if
       
       ! Step 11: Populate RHS for scalars
       !if (allocated(this%scalars)) then
       if (this%useScalars) then
           do idx = 1,this%n_scalars
              call this%scalars(idx)%populateRHS(this%dt)
           end do 
       end if
   end subroutine

   subroutine project_and_prep(this, AlreadyProjected)
       class(igrid), intent(inout) :: this
       logical, intent(in) :: AlreadyProjected
       integer :: idx 

       ! Step 1: Dealias
       call this%dealiasFields()

       ! Step 2: Pressure projection
       if (.not. AlreadyProjected) then
           call this%padepoiss%PressureProjection(this%uhat,this%vhat,this%what)
           if (mod(this%step,this%t_DivergenceCheck) == 0) then
               call this%padepoiss%DivergenceCheck(this%uhat, this%vhat, this%what, this%divergence,.true.)
           end if 
       end if 

       ! Step 3: Take it back to physical fields
       call this%spectC%ifft(this%uhat,this%u)
       call this%spectC%ifft(this%vhat,this%v)
       call this%spectE%ifft(this%what,this%w)
       if (this%isStratified .or. this%initspinup) call this%spectC%ifft(this%That,this%T)
   
       ! STEP 4: Interpolate the cell center values of w
       !call this%compute_and_bcast_surface_Mn()
       call this%interp_PrimitiveVars()

       ! STEP 5: Compute duidxjC 
       call this%compute_duidxj()
       if (this%isStratified .or. this%initspinup) call this%compute_dTdxi() 
        
       ! STEP 6: Prep scalar fields if being used. 
       if (this%useScalars) then
          do idx = 1,this%n_scalars
             call this%scalars(idx)%prep_scalar()
          end do 
       end if
   end subroutine
   
   subroutine wrapup_timestep(this)
       class(igrid), intent(inout) :: this

       logical :: forceWrite, exitStat, forceDumpPressure, restartWrite, forceDumpProbes 
       integer :: ierr = -1, ierr2

       ! STEP 1: Update Time, BCs and record probe data
       this%step = this%step + 1; this%tsim = this%tsim + this%dt
       this%newTimeStep = .true. 
       if (this%useControl .AND. abs(180.d0/pi*this%angleHubHeight) > this%angleTrigger) then
           this%G_alpha = this%G_alpha - this%deltaGalpha
           this%frameAngle = this%frameAngle + this%deltaGalpha 
       end if
       
       !if (this%isStratified) then
       !    if (this%botBC_Temp == 0) then
       !        this%Tsurf = this%Tsurf0 + this%dTsurf_dt*this%tsim
       !    end if
       !end if  
       if (this%PreprocessForKS) this%KSupdated = .false. 
       if (this%useProbes) call this%updateProbes()
       if (this%computevorticity) call this%compute_vorticity

       ierr = -1; forceWrite = .FALSE.; exitstat = .FALSE.; forceDumpPressure = .FALSE.; 
       forceDumpProbes = .false.; restartWrite = .FALSE. 

       if(this%tsim > this%tstop) then
         forceWrite = .TRUE.
         restartWrite = .TRUE.
         if (this%useProbes) forceDumpProbes = .TRUE.
         call message(0,"The simulation has ended.")
         call message(1,"Dumping a restart file.")
       endif

       if (this%useSystemInteractions) then
           if (mod(this%step,this%tSystemInteractions) == 0) then
               exitStat = check_exit(this%controlDir)
               if (exitStat) forceWrite = .true.
               open(777,file=trim(this%controlDir)//"/dumppdo",status='old',iostat=ierr)
               if(ierr==0) then
                   forceWrite = .TRUE.
                   call message(1, "Forced Dump because found file dumppdo")
                   call message(2, "Current Time Step is:", this%step)
                   if (nrank .ne. 0) close(777)
                   call mpi_barrier(mpi_comm_world, ierr2)
                   !if(nrank==0) close(777, status='delete')
                   if (this%deleteInstructions) then
                      if(nrank==0) close(777, status='delete')
                   else
                      if(nrank==0) close(777)
                   end if
               else
                   close(777)
               endif
          

               if (this%useProbes) then
                   open(777,file=trim(this%controlDir)//"/dumpprobes",status='old',iostat=ierr)
                   if (ierr==0) then
                       forceDumpProbes = .true. 
                       call message(1, "Forced Dump for PROBES because found file dumpprobes")
                       call message(2, "Current Time Step is:", this%step)
                       if (nrank .ne. 0) close(777)
                       call mpi_barrier(mpi_comm_world, ierr2)
                       !if(nrank==0) close(777, status='delete')
                       if (this%deleteInstructions) then
                          if(nrank==0) close(777, status='delete')
                       else
                          if(nrank==0) close(777)
                       end if
                   else
                       close(777)
                   end if
               end if

               if (this%storePressure) then 
                   open(777,file=trim(this%controlDir)//"/prsspdo",status='old',iostat=ierr)
                   if (ierr == 0) then
                       forceDumpPressure = .true.
                       call message(1,"Force Dump for pressure reqested; file prsspdo found.")
                       if (nrank .ne. 0) close(777)
                       call mpi_barrier(mpi_comm_world, ierr2)
                       !if (nrank==0) close(777, status='delete')
                       if (this%deleteInstructions) then
                          if(nrank==0) close(777, status='delete')
                       else
                          if(nrank==0) close(777)
                       end if
                   else
                       close(777)
                   end if
               end if

               open(777,file=trim(this%controlDir)//"/dumprestart",status='old',iostat=ierr)
               if(ierr==0) then
                   restartWrite = .TRUE.
                   call message(1, "Restart Dump because found file dumprestart")
                   call message(2, "Current Time Step is:", this%step)
                   if (nrank .ne. 0) close(777)
                   call mpi_barrier(mpi_comm_world, ierr2)
                   !if(nrank==0) close(777, status='delete')
                   if (this%deleteInstructions) then
                      if(nrank==0) close(777, status='delete')
                   else
                      if(nrank==0) close(777)
                   end if
               else
                   close(777)
               endif

           end if 
       end if 

       ! STEP 2: Do logistical stuff
       if (this%fastCalcPressure) then
           call this%computePressure()
       else
           if ((this%storePressure)) then
               if ((mod(this%step,this%P_compFreq)==0) .or. (forceDumpPressure)) then
                   call this%computePressure()
               end if 
               if ( (mod(this%step,this%P_dumpFreq) == 0).or. (forceDumpPressure)) then
                   call this%dumpFullField(this%pressure,"prss")
                   if (this%computeDNSpressure) call this%dumpFullField(this%pressure_dns,'pdns')
                   if (this%computeturbinepressure) call this%dumpFullField(this%pressure_turbine,'ptrb')
                   if (this%computefringepressure) call this%dumpFullField(this%pressure_fringe,'pfrn')
               end if 
           end if 
       end if

      if ((forceWrite.or.(mod(this%step,this%tid_compStats)==0)).and.(this%tsim > this%tSimStartStats) ) then
          if (this%timeAvgFullFields) then
              ! rhs needs to be evaluated before computing statistics to ensure
              ! that tauSGS are consistent with the velocities and velocity
              ! derivatives, which is needed for correct SGS dissipation
              if (.not. this%AlreadyHaveRHS) then
                  call this%populate_rhs()
                  this%AlreadyHaveRHS = .true. 
              end if 

              call this%compute_stats3D()
          else
          !    call this%compute_stats()
          end if 
      end if 

      if ((forceWrite.or.(mod(this%step,this%tid_statsDump)==0)).and.(this%tsim > this%tSimStartStats) ) then
         if (this%timeAvgFullFields) then
              call this%dump_stats3D()
              call mpi_barrier(mpi_comm_world, ierr)
              !stop
          else
              !call this%compute_stats()
              !call this%dump_stats()
          end if
      end if 

      if ( restartWrite .or. (mod(this%step,this%t_restartDump) == 0) ) then
          call this%dumpRestartfile()
          call message(0,"Scheduled restart file dumped.")
      end if
      
      if ( (forceWrite .or. ((mod(this%step,this%t_planeDump) == 0) .and. &
               (this%step .ge. this%t_start_planeDump) .and. (this%step .le. this%t_stop_planeDump))) .and. (this%dumpPlanes)) then
          if (this%PreprocessForKS) then
              if (.not. this%KSupdated) then
                  call this%LES2KS%applyFilterForKS(this%u, this%v, this%wC)
                  this%KSupdated = .true. 
              end if
          end if
          call this%dump_planes()
      end if 


      if (this%useProbes) then
          if (forceDumpProbes) then
              call this%dumpProbes()    
              call message(0,"Performed a forced dump for probes.")
          end if
      end if

      if (this%vizDump_Schedule == 1) then
           if (this%DumpThisStep) then
              call message(2,"Performing a fixed timed visualization dump at time:", this%tsim)
              call message(2,"This time step used a deltaT:",this%dt)
              call this%dump_visualization_files()
           end if 
      else
           if (mod(this%step,this%t_dataDump) == 0) then
              call message(0,"Scheduled visualization dump.")
              call this%dump_visualization_files()
           end if
      end if 


      if (forceWrite) then
         call message(2,"Performing a forced visualization dump.")
         call this%dump_visualization_files()
      end if

      if (this%initspinup) then
       if (this%tsim > this%Tstop_initspinup) then
           this%initspinup = .false.
           call message(0,"Initialization spin up turned off. No active scalar in the problem.")
       end if 
      end if 
      ! Exit if the exitpdo file was detected earlier at the beginning of this
      ! subroutine
      if(exitstat) call GracefulExit("Found exitpdo file in control directory",1234)

   end subroutine
   
   subroutine AdamsBashforth(this)
       class(igrid), intent(inout) :: this
       real(rkind) :: abf1, abf2

       ! Step 0: Compute TimeStep 
       call this%compute_deltaT
       this%dtRat = this%dt/this%dtOld

       ! Step 1: Get the RHS
       if (this%AlreadyHaveRHS) then
           this%AlreadyHaveRHS = .false.
       else
           call this%populate_rhs()
       end if 

       ! Step 2: Time Advance
       if (this%step == 0) then
           this%uhat = this%uhat + this%dt*this%u_rhs 
           this%vhat = this%vhat + this%dt*this%v_rhs 
           this%what = this%what + this%dt*this%w_rhs 
           if (this%isStratified .or. this%initspinup) then
               this%That = this%That + this%dt*this%T_rhs
           end if
       else
           abf1 = (one + half*this%dtRat)*this%dt
           abf2 = -half*this%dtRat*this%dt
           this%uhat = this%uhat + abf1*this%u_rhs + abf2*this%u_Orhs
           this%vhat = this%vhat + abf1*this%v_rhs + abf2*this%v_Orhs
           this%what = this%what + abf1*this%w_rhs + abf2*this%w_Orhs
           if (this%isStratified .or. this%initspinup) then
               this%That = this%That + abf1*this%T_rhs + abf2*this%T_Orhs
           end if 
       end if 

       ! Step 3: Pressure Project and prep for the next step
       call this%project_and_prep(.false.)

       ! Step 4: Store the RHS values for the next use
       this%u_Orhs = this%u_rhs; this%v_Orhs = this%v_rhs; this%w_Orhs = this%w_rhs
       if (this%isStratified .or. this%initspinup) this%T_Orhs = this%T_rhs
       this%dtOld = this%dt

       ! Step 5: Do end of time step operations (I/O, stats, etc.)
       call this%wrapup_timestep()
   end subroutine
   
   function get_dt(this, recompute) result(val)
       class(igrid), intent(inout) :: this 
       logical, intent(in), optional :: recompute
       real(rkind) :: val

       if (present(recompute)) then
          if (recompute) call this%compute_deltaT
       end if


       val = this%dt

   end function
   
   subroutine compute_deltaT(this)
       use reductions, only: p_maxval
       class(igrid), intent(inout), target :: this
       real(rkind) :: TSmax, Tsim_next 
       real(rkind), dimension(:,:,:), pointer :: rb1, rb2
       real(rkind), dimension(5) :: dtmin
       integer :: idx

       rb1 => this%rbuffxC(:,:,:,1)
       rb2 => this%rbuffxC(:,:,:,2)

       if (this%useCFL) then
           rb1 = (one/this%dx)*this%u
           rb2 = (one/this%dy)*this%v 
           rb1 = abs(rb1) 
           rb2 = abs(rb2)
           rb1 = rb1 + rb2
           rb2 = (one/this%dz)*this%wC
           rb2 = abs(rb2)
           rb1 = rb1 + rb2
           TSmax = p_maxval(rb1)
           dtmin(1)= this%CFL/TSmax
              
           if (.not. this%isInviscid) then
              dtmin(2) = this%CviscDT*0.5d0*this%Re*(min(this%dx,this%dy,this%dz)**2)
              if (this%isStratified .and. (this%PrandtlFluid > 1.d0)) then 
                 dtmin(2) = dtmin(2) / this%PrandtlFluid  
              end if 
           else
              dtmin(2) = 1.d15
           end if
           
           if (associated(this%nu_SGS)) then
              dtmin(3) = this%CviscDT*0.25d0*(min(this%dx,this%dy,this%dz)**2)/(p_maxval(this%nu_SGS) + 1.d-18)
           else
              dtmin(3) = 1.d15
           end if 

           if (associated(this%kappaSGS)) then
              dtmin(4) = this%CviscDT*0.5d0*(min(this%dx,this%dy,this%dz)**2)/(p_maxval(this%kappaSGS) + 1.d-18)
           else
              dtmin(4) = 1.d15
           end if 

           if (associated(this%kappa_bounding)) then
              dtmin(5) = this%CviscDT*0.5d0*(min(this%dx,this%dy,this%dz)**2)/(p_maxval(this%kappa_bounding) + 1.d-18)
           else
              dtmin(5) = 1.d15
           end if 

           this%dt = minval(dtmin)
           idx = minloc(dtmin, DIM=1)
           select case(idx)
           case (1) 
              this%dtlimit = " Convective CFL"
           case(2)
              this%dtlimit = " Viscous"
           case(3)
              this%dtlimit = " SGS viscosity"
           case(4)
              this%dtlimit = " SGS scalar diffusivity"
           case(5)
              this%dtlimit = " Scalar bounding diffusivity"
           end select 
        else
           this%dtlimit = "invalid/fixed dt"
       end if 

       if (this%vizDump_Schedule == 1) then
           this%DumpThisStep = .false. 
           Tsim_next = this%tsim + this%dt
           if (Tsim_next > this%t_NextDump) then
               this%dt = this%t_nextDump - this%tsim
               this%t_NextDump = this%t_NextDump + this%deltaT_dump 
               this%DumpThisStep = .true.
           end if
       end if 

   end subroutine
