   subroutine populate_rhs(this, CopyForDNSpress, CopyForFringePress, copyForTurbinePress)
      class(igrid), intent(inout) :: this
      !integer,           intent(in)    :: RKstage
      logical, intent(in), optional :: CopyForDNSpress, CopyForFringePress, copyForTurbinePress
      logical :: copyFringeRHS, copyTurbRHS
       

      if (present(copyForTurbinePress)) then
         copyTurbRHS = copyForTurbinePress
      else
         copyTurbRHS = .false. 
      end if

       if (present(copyForFringePress)) then
           copyFringeRHS = copyForFringePress
       else
           copyFringeRHS = .false. 
       end if
      
       ! Step 1: Non Linear Term 
       if (useSkewSymm) then
           call this%addNonLinearTerm_skewSymm(this%u_rhs, this%v_rhs, this%w_rhs)
       else
           call this%AddNonLinearTerm_Rot(this%u_rhs, this%v_rhs, this%w_rhs)
       end if

       ! Step 2: Coriolis Term
       if (this%useCoriolis) call this%AddCoriolisTerm(this%u_rhs, this%v_rhs, this%w_rhs)
       !
       ! Step 3b: Wind Turbines
       if (this%useWindTurbines) then          
          if (copyTurbRHS) then
           call this%WindTurbineArr%getForceRHS(this%dt, this%u, this%v, this%wC, this%u_rhs, this%v_rhs, this%w_rhs, &
                 this%newTimestep, this%inst_horz_avg_turb, uturb=this%urhs_turbine, vturb=this%vrhs_turbine, wturb=this%wrhs_turbine)
          else
               call this%WindTurbineArr%getForceRHS(this%dt, this%u, this%v, this%wC,&
                                    this%u_rhs, this%v_rhs, this%w_rhs, this%newTimestep, this%inst_horz_avg_turb)
          end if
       end if 
     
       ! Step 4: Buoyance + Sponge (inside Buoyancy)
       if (this%isStratified .or. this%initspinup) then
           call this%addBuoyancyTerm(this%u_rhs, this%v_rhs, this%w_rhs)
       end if 
       
       if (present(CopyForDNSpress)) then
           if (CopyForDNSpress) then
              if (copyTurbRHS) then
                 this%urhs_dns = this%u_rhs - this%urhs_turbine; this%vrhs_dns = this%v_rhs - this%vrhs_turbine; this%wrhs_dns = this%w_rhs - this%wrhs_turbine
              else
                 this%urhs_dns = this%u_rhs; this%vrhs_dns = this%v_rhs; this%wrhs_dns = this%w_rhs
              end if
           end if   
       end if

       ! Step 6a: SGS and viscous stress Terms
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
           ! IMPORTANT: If SGS model is used, the viscous term is evaluated
           ! as part of the SGS stress tensor. 
           if (.not. this%isInviscid) then
               call this%addViscousTerm(this%u_rhs, this%v_rhs, this%w_rhs)
           end if
       end if


       call this%populate_RHS_extraTerms(copyFringeRHS, .false.)

   end subroutine

   
   subroutine populate_rhs_for_budgets(this, CopyForDNSpress, CopyForFringePress, copyForTurbinePress)
      class(igrid), intent(inout) :: this
      !integer,           intent(in)    :: RKstage
      logical, intent(in), optional :: CopyForDNSpress, CopyForFringePress, copyForTurbinePress
      logical :: copyFringeRHS, copyTurbRHS
      
      call this%set_budget_rhs_to_zero()

      if (present(copyForTurbinePress)) then
         copyTurbRHS = copyForTurbinePress
      else
         copyTurbRHS = .false. 
      end if

       if (present(copyForFringePress)) then
           copyFringeRHS = copyForFringePress
       else
           copyFringeRHS = .false. 
       end if
      
       ! Step 1: Non Linear Term 
       if (useSkewSymm) then
           !call this%addNonLinearTerm_skewSymm(this%ucon, this%vcon, this%wcon)
           call this%addNonLinearTerm_skewSymm(this%u_rhs, this%v_rhs, this%w_rhs)
       else
           !call this%AddNonLinearTerm_Rot(this%ucon, this%vcon, this%wcon)
           call this%AddNonLinearTerm_Rot(this%u_rhs, this%v_rhs, this%w_rhs)
       end if
       if(associated(this%ucon)) then
           this%ucon = this%u_rhs 
           this%vcon = this%v_rhs 
           this%wcon = this%w_rhs
       endif

       ! Step 2: Coriolis Term
       if ((this%useCoriolis) .and. associated(this%ucor)) then 
           call this%AddCoriolisTerm(this%ucor, this%vcor, this%wcor)
           this%u_rhs = this%u_rhs + this%ucor
           this%v_rhs = this%v_rhs + this%vcor
           this%w_rhs = this%w_rhs + this%wcor
       end if 

       ! Step 3b: Wind Turbines
       if ((this%useWindTurbines) .and. associated(this%uturb)) then
          if (copyTurbRHS) then
           call this%WindTurbineArr%getForceRHS(this%dt, this%u, this%v, this%wC, this%uturb, this%vturb, this%wturb, &
                 & this%newTimestep, this%inst_horz_avg_turb, uturb=this%urhs_turbine, vturb=this%vrhs_turbine, wturb=this%wrhs_turbine) 
          else
               call this%WindTurbineArr%getForceRHS(this%dt, this%u, this%v, this%wC,&
                                    !this%uturb, this%v_rhs, this%w_rhs, this%newTimestep, this%inst_horz_avg_turb)
                                    this%uturb, this%vturb, this%wturb, this%newTimestep, this%inst_horz_avg_turb)
          end if
          this%u_rhs = this%u_rhs + this%uturb
          this%v_rhs = this%v_rhs + this%vturb
          this%w_rhs = this%w_rhs + this%wturb
       end if 
     
       ! Step 4: Buoyance + Sponge (inside Buoyancy)
       if ((this%isStratified .or. this%initspinup) .and. (associated(this%wb))) then
           call this%addBuoyancyTerm(this%u_rhs, this%v_rhs, this%wb)
           this%w_rhs = this%w_rhs + this%wb
       end if 
       
       if (present(CopyForDNSpress)) then
           if (CopyForDNSpress) then
              if (copyTurbRHS) then
                 this%urhs_dns = this%u_rhs - this%urhs_turbine; this%vrhs_dns = this%v_rhs - this%vrhs_turbine; this%wrhs_dns = this%w_rhs - this%wrhs_turbine
              else
                 this%urhs_dns = this%u_rhs; this%vrhs_dns = this%v_rhs; this%wrhs_dns = this%w_rhs
              end if
           end if   
       end if

       ! Step 6a: SGS  Stress Terms
       if (this%useSGS) then
         if(associated(this%usgs)) then
           call this%sgsmodel%getRHS_SGS(this%usgs, this%vsgs, this%wsgs,      this%duidxjC, this%duidxjE,    &
                                         this%uhat,  this%vhat,  this%whatC,      this%That,    this%u,       &
                                         this%v,     this%wC,    this%newTimeStep,this%dTdxC,   this%dTdyC,   & 
                                         this%dTdzC, this%dTdxE, this%dTdyE, this%dTdzE)
         else
           call this%sgsmodel%getRHS_SGS(this%u_rhs, this%v_rhs, this%w_rhs,      this%duidxjC, this%duidxjE, &
                                         this%uhat,  this%vhat,  this%whatC,      this%That,    this%u,       &
                                         this%v,     this%wC,    this%newTimeStep,this%dTdxC,   this%dTdyC,   & 
                                         this%dTdzC, this%dTdxE, this%dTdyE, this%dTdzE)
         end if

           if (this%isStratified .or. this%initspinup) then
              call this%sgsmodel%getRHS_SGS_Scalar(this%T_rhs, this%dTdxC, this%dTdyC, this%dTdzC, this%dTdzE, &
                                         this%u, this%v, this%wC, this%T, this%That, this%duidxjC, this%turbPr)
           end if

         if(associated(this%usgs)) then
           this%u_rhs = this%u_rhs + this%usgs
           this%v_rhs = this%v_rhs + this%vsgs
           this%w_rhs = this%w_rhs + this%wsgs
         end if
           
           ! viscous term evaluate separately  
           if ((.not. this%isInviscid) .and. associated(this%uvisc)) then
               call this%addViscousTerm(this%uvisc, this%vvisc, this%wvisc)
               
               ! Now correct the SGS term 
               this%usgs = this%usgs - this%uvisc 
               this%vsgs = this%vsgs - this%vvisc 
               this%wsgs = this%wsgs - this%wvisc 
           end if  
          
       else
           ! IMPORTANT: If SGS model is used, the viscous term is evaluated
           ! as part of the SGS stress tensor. 
           if ((.not. this%isInviscid)) then
               call this%addViscousTerm(this%u_rhs, this%v_rhs, this%w_rhs)
           end if
       end if


       call this%populate_RHS_extraTerms(copyFringeRHS, .true.)

   end subroutine

   subroutine populate_RHS_extraTerms(this, copyFringeRHS, storeForBudget)
       class(igrid), intent(inout) :: this
       logical, intent(in) :: copyFringeRHS, storeForBudget
       integer :: idx 
       
       ! Step 7a: Extra Forcing 
       if (this%useExtraForcing) then
           call this%addExtraForcingTerm()
       end if
       
       ! Step 7b: HIT forcing source term
       if (this%useHITForcing) then
           if(storeForBudget) then
               call this%hitforce%getRHS_HITForcing(this%HITforcing_x, this%HITforcing_y, this%HITforcing_z, this%uhat, this%vhat, this%what, this%newTimeStep)
               !print '(a,i4.4,3(1x,e19.12))', '--------', nrank, maxval(abs(this%HITforcing_x)), maxval(abs(this%HITforcing_y)), maxval(abs(this%HITforcing_z))
               this%u_rhs = this%u_rhs + this%HITforcing_x
               this%v_rhs = this%v_rhs + this%HITforcing_y
               this%w_rhs = this%w_rhs + this%HITforcing_z
           else
               call this%hitforce%getRHS_HITForcing(this%u_rhs, this%v_rhs, this%w_rhs, this%uhat, this%vhat, this%what, this%newTimeStep)
           end if
       end if

       if (this%useHITRealSpaceLinearForcing) then
           this%rbuffxC(:,:,:,1) = this%u/this%HITForceTimeScale
           call this%spectC%fft(this%rbuffxC(:,:,:,1),this%cbuffyC(:,:,:,1))
           this%u_rhs = this%u_rhs + this%cbuffyC(:,:,:,1)
           if(storeForBudget) then
               this%HITforcing_x = this%HITforcing_x + this%cbuffyC(:,:,:,1)
           end if

           this%rbuffxC(:,:,:,1) = this%v/this%HITForceTimeScale
           call this%spectC%fft(this%rbuffxC(:,:,:,1),this%cbuffyC(:,:,:,1))
           this%v_rhs = this%v_rhs + this%cbuffyC(:,:,:,1)
           if(storeForBudget) then
               this%HITforcing_y = this%HITforcing_y + this%cbuffyC(:,:,:,1)
           end if
           
           this%rbuffxE(:,:,:,1) = this%w/this%HITForceTimeScale
           call this%spectE%fft(this%rbuffxE(:,:,:,1),this%cbuffyE(:,:,:,1))
           this%w_rhs = this%w_rhs + this%cbuffyE(:,:,:,1)
           if(storeForBudget) then
               this%HITforcing_z = this%HITforcing_z + this%cbuffyE(:,:,:,1)
           end if
       end if 
       
       ! Step 8: Fringe and sponge source terms
       if (this%useSponge) then
           call this%addSponge()
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
           if (this%isStratified .or. this%initspinup) then
               call this%fringe_x1%addFringeRHS_scalar(this%dt, this%T_rhs, this%T)
               call this%fringe_x2%addFringeRHS_scalar(this%dt, this%T_rhs, this%T)
           end if
       else
           if (this%useFringe) then
              if (copyFringeRHS) then
                  call this%fringe_x%addFringeRHS(this%dt, this%u_rhs, this%v_rhs, this%w_rhs, this%u, this%v, this%w, & 
                                   this%urhs_fringe, this%vrhs_fringe, this%wrhs_fringe, addF=.false. )
              else
                  call this%fringe_x%addFringeRHS(this%dt, this%u_rhs, this%v_rhs, this%w_rhs, this%u, this%v, this%w)
              end if 
              if (this%isStratified .or. this%initspinup) then
                  call this%fringe_x%addFringeRHS_scalar(this%dt, this%T_rhs, this%T)
              end if
           end if 
       end if 

       ! Step 9: Frame rotatio PI controller to fix yaw angle at a given height
       if (this%useControl .AND. abs(180.d0/pi*this%angleHubHeight)>0.0d0) then
           call this%angCont_yaw%update_RHS_control(this%dt, this%u_rhs, this%v_rhs, &
                         this%w_rhs, this%u, this%v, this%newTimeStep, this%angleHubHeight, this%wFilt, this%deltaGalpha, this%zHubIndex, this%angleTrigger)
           this%totalAngle = this%totalAngle + this%angleHubHeight
       end if 

       ! Step 10: Populate RHS for scalars
       !if (allocated(this%scalars)) then
       if (this%useScalars) then
           do idx = 1,this%n_scalars
              call this%scalars(idx)%populateRHS(this%dt)
           end do 
       end if

       ! Step 11: Add statified forcing 
       if (this%isStratified .and. this%useforcedStratification) then
            call this%addForcedStratification()
       end if 

   end subroutine 
