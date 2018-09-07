   !! STATISTICS !!

   !--------------------------------Beginning 3D Statistics----------------------------------------------
   subroutine init_stats3D(this)
       use exits, only: message
       class(igrid), intent(inout), target :: this
       integer :: nstatsvar, nhorzavgvars, nstv, nhzv


       nstatsvar = 12; nhorzavgvars = 17
       if(this%fastCalcPressure .or. this%storePressure) then
          nstatsvar = nstatsvar + 4
          nhorzavgvars = nhorzavgvars + 3
       endif
       if(.not. this%isInviscid) then
          nstatsvar = nstatsvar + 4
          nhorzavgvars = nhorzavgvars + 4
       endif
       if(this%useSGS) then
          nstatsvar = nstatsvar + 10
          nhorzavgvars = nhorzavgvars + 10
       endif
       if(this%useWindTurbines) then
          nstatsvar = nstatsvar + 4
          nhorzavgvars = nhorzavgvars + 5
       endif
       if(this%isStratified) then
          nstatsvar = nstatsvar + 5
          nhorzavgvars = nhorzavgvars + 5
          if(this%useSGS) then
             nstatsvar = nstatsvar + 3
             nhorzavgvars = nhorzavgvars + 3
          endif
       endif

       allocate(this%debugavg(5),this%debuginst(5))
       allocate(this%inst_horz_avg(5)) ! [ustar, uw, vw, Linv, wT]
       allocate(this%runningSum_sc(5))
       allocate(this%stats3D(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),nstatsvar))
       allocate(this%horzavgstats(this%nz,nhorzavgvars))

       if(this%useWindTurbines) then
           allocate(this%runningSum_sc_turb(8*this%WindTurbineArr%nTurbines))
           allocate(this%runningSum_turb   (8*this%WindTurbineArr%nTurbines))
       endif

       if (this%computeSpectra) then
           allocate(this%xspectra_mean(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))   ! ensure that number of variables for which spectrum is to be computed is smaller than nyg
       end if 

       ! mean velocities
       this%u_mean3D => this%stats3D(:,:,:,1);  this%v_mean3D  => this%stats3D(:,:,:,2);  this%w_mean3D => this%stats3D(:,:,:,3) 
       ! mean squared velocities
       this%uu_mean3D => this%stats3D(:,:,:,4); this%uv_mean3D => this%stats3D(:,:,:,5); this%uw_mean3D => this%stats3D(:,:,:,6)
                                                this%vv_mean3D => this%stats3D(:,:,:,7); this%vw_mean3D => this%stats3D(:,:,:,8) 
                                                                                         this%ww_mean3D => this%stats3D(:,:,:,9)

       ! triple product of velocities
       this%tketurbtranspx_mean3D => this%stats3D(:,:,:,10)
       this%tketurbtranspy_mean3D => this%stats3D(:,:,:,11)
       this%tketurbtranspz_mean3D => this%stats3D(:,:,:,12)
       nstv = 3 + 6 + 3

       if(this%fastCalcPressure .or. this%storePressure) then
          this%p_mean3D  => this%stats3D(:,:,:,nstv+1)
          this%pu_mean3D => this%stats3D(:,:,:,nstv+2)
          this%pv_mean3D => this%stats3D(:,:,:,nstv+3)
          this%pw_mean3D => this%stats3D(:,:,:,nstv+4)
          nstv = nstv + 4
       endif

       if(.not. this%isInviscid) then
          this%viscdisp_mean3D => this%stats3D(:,:,:,nstv+1)
          this%Siju1_mean3D    => this%stats3D(:,:,:,nstv+2)
          this%Siju2_mean3D    => this%stats3D(:,:,:,nstv+3)
          this%Siju3_mean3D    => this%stats3D(:,:,:,nstv+4)
          nstv = nstv + 4
       endif

       if(this%useSGS)  then
          ! SGS stresses
          this%tau11_mean3D => this%stats3D(:,:,:,nstv+1); this%tau12_mean3D => this%stats3D(:,:,:,nstv+2); this%tau13_mean3D => this%stats3D(:,:,:,nstv+3)
                                                           this%tau22_mean3D => this%stats3D(:,:,:,nstv+4); this%tau23_mean3D => this%stats3D(:,:,:,nstv+5) 
                                                                                                            this%tau33_mean3D => this%stats3D(:,:,:,nstv+6)
          ! SGS dissipation
          this%sgsdissp_mean3D => this%stats3D(:,:,:,nstv+7)

          this%tauu1_mean3D    => this%stats3D(:,:,:,nstv+8)
          this%tauu2_mean3D    => this%stats3D(:,:,:,nstv+9)
          this%tauu3_mean3D    => this%stats3D(:,:,:,nstv+10)

          nstv = nstv + 10

       endif

       if(this%useWindTurbines)  then
          this%turbfx_mean3D => this%stats3D(:,:,:,nstv+1)
          this%turbfy_mean3D => this%stats3D(:,:,:,nstv+2)
          this%turbfz_mean3D => this%stats3D(:,:,:,nstv+3)
          this%uturbf_mean3D => this%stats3D(:,:,:,nstv+4)
          nstv = nstv + 4
       endif

       if (this%isStratified) then
          this%TT_mean3D => this%stats3D(:,:,:,nstv+1);  this%wT_mean3D => this%stats3D(:,:,:,nstv+2);  this%vT_mean3D => this%stats3D(:,:,:,nstv+3)
          this%uT_mean3D => this%stats3D(:,:,:,nstv+4);  this%T_mean3D  => this%stats3D(:,:,:,nstv+5);
          nstv = nstv + 5
          if(this%useSGS) then
             this%q1_mean3D => this%stats3D(:,:,:,nstv+1);  this%q2_mean3D => this%stats3D(:,:,:,nstv+2);  this%q3_mean3D => this%stats3D(:,:,:,nstv+3)
             nstv = nstv + 3
          endif
       end if 

       ! horizontal averages
       ! mean velocities
       this%u_mean   => this%horzavgstats(:,1);  this%v_mean    => this%horzavgstats(:,2);  this%w_mean   => this%horzavgstats(:,3) 
       ! mean squared velocities
       this%uu_mean => this%horzavgstats(:,4); this%uv_mean => this%horzavgstats(:,5); this%uw_mean => this%horzavgstats(:,6)
                                               this%vv_mean => this%horzavgstats(:,7); this%vw_mean => this%horzavgstats(:,8) 
                                                                                       this%ww_mean => this%horzavgstats(:,9)
       ! Dispersive stresses
       this%disperuw_mean => this%horzavgstats(:,10); this%dispervw_mean => this%horzavgstats(:,11)

       ! Mean Kinetic Energy Equation: advection, turbulent transport, dissipation
       this%mkeadv_mean => this%horzavgstats(:,12); this%mkett_mean => this%horzavgstats(:,13); this%mkedisp_mean => this%horzavgstats(:,14)
       ! Turbulent Kinetic Energy Equation: advection, turbulent transport, shear production
       this%tkeadv_mean => this%horzavgstats(:,15); this%tkett_mean => this%horzavgstats(:,16); this%tkeprod_mean => this%horzavgstats(:,17)
       nhzv = 17

       if(this%fastCalcPressure .or. this%storePressure) then
          this%p_mean  => this%horzavgstats(:,nhzv+1)
          this%mkept_mean => this%horzavgstats(:,nhzv+2)
          this%tkept_mean => this%horzavgstats(:,nhzv+3)
          nhzv = nhzv + 3
       endif

       if(.not. this%isInviscid) then
          this%mkevdif_mean => this%horzavgstats(:,nhzv+1)
          this%mkevdsp_mean => this%horzavgstats(:,nhzv+2)
          this%tkevdif_mean => this%horzavgstats(:,nhzv+3)
          this%tkevdsp_mean => this%horzavgstats(:,nhzv+4)
          nhzv = nhzv + 4
       endif

       if(this%useSGS) then
          ! SGS stresses
          this%tau11_mean => this%horzavgstats(:,nhzv+1); this%tau12_mean => this%horzavgstats(:,nhzv+2); this%tau13_mean => this%horzavgstats(:,nhzv+3)
                                                          this%tau22_mean => this%horzavgstats(:,nhzv+4); this%tau23_mean => this%horzavgstats(:,nhzv+5) 
                                                                                                          this%tau33_mean => this%horzavgstats(:,nhzv+6)
          ! SGS dissipation
          this%mkesgst_mean => this%horzavgstats(:,nhzv+7)
          this%mkesgsd_mean => this%horzavgstats(:,nhzv+8)
          this%tkesgst_mean => this%horzavgstats(:,nhzv+9)
          this%tkesgsd_mean => this%horzavgstats(:,nhzv+10)
          nhzv = nhzv + 10

       endif

       if(this%useWindTurbines) then
          this%turbfx_mean => this%horzavgstats(:,nhzv+1)
          this%turbfy_mean => this%horzavgstats(:,nhzv+2)
          this%turbfz_mean => this%horzavgstats(:,nhzv+3)
          this%mketurbf_mean => this%horzavgstats(:,nhzv+4)
          this%tketurbf_mean => this%horzavgstats(:,nhzv+5)
          nhzv = nhzv + 5
       endif

       if (this%isStratified) then
           this%T_mean  => this%horzavgstats(:,nhzv+1);  this%uT_mean => this%horzavgstats(:,nhzv+2);  this%vT_mean => this%horzavgstats(:,nhzv+3)
           this%wT_mean => this%horzavgstats(:,nhzv+4);  this%TT_mean => this%horzavgstats(:,nhzv+5)
           nhzv = nhzv + 5 
          if(this%useSGS) then
             this%q1_mean => this%horzavgstats(:,nhzv+1); this%q2_mean => this%horzavgstats(:,nhzv+2);  this%q3_mean => this%horzavgstats(:,nhzv+3)
             nhzv = nhzv + 3
          endif
       end if 

       this%tidSUM = 0
       this%tprev2 = -1; this%tprev1 = -1;
       this%stats3D = zero
       this%horzavgstats = zero
       this%debugavg = zero
       this%debuginst = zero
       this%inst_horz_avg = zero
       this%runningSum_sc = zero
       if(this%useWindTurbines) then
           this%inst_horz_avg_turb = zero
           this%runningSum_sc_turb = zero
           this%runningSum_turb    = zero
       endif
       this%xspectra_mean = zero

       if((nhzv .ne. nhorzavgvars) .or. (nstv .ne. nstatsvar)) then
           call message(0,"Error in init_stats3D")
           write(*,*) 'nhzv = ', nhzv, nhorzavgvars
           write(*,*) 'nstv = ', nstv, nstatsvar
       endif

       call message(0,"Done init_stats3D")
   end subroutine

   subroutine compute_stats3D(this)
       use kind_parameters, only: mpirkind
       class(igrid), intent(inout), target :: this
       real(rkind), dimension(:,:,:), pointer :: rbuff0, rbuff1, rbuff2, rbuff2E, rbuff3E, rbuff3, rbuff1E, rbuff4
       integer :: j, k, jindx, ierr
       real(rkind),    dimension(:,:,:), pointer :: dudx, dudy
       real(rkind),    dimension(:,:,:), pointer :: dvdx, dvdy
       real(rkind),    dimension(:,:,:), pointer :: dwdz
       real(rkind),    dimension(:,:,:), pointer :: dvdzC, dudzC
       real(rkind),    dimension(:,:,:), pointer :: dwdxC, dwdyC

       rbuff0  => this%rbuffxC(:,:,:,1); rbuff1  => this%rbuffxC(:,:,:,2);
       rbuff2  => this%rbuffyC(:,:,:,1);
       rbuff2E => this%rbuffyE(:,:,:,1); rbuff3E => this%rbuffzE(:,:,:,1);
       rbuff3 => this%rbuffzC(:,:,:,1);  rbuff1E => this%rbuffxE(:,:,:,1)
       rbuff4 => this%rbuffzC(:,:,:,2);

       dudx  => this%duidxjC(:,:,:,1); dudy  => this%duidxjC(:,:,:,2); dudzC => this%duidxjC(:,:,:,3); 
       dvdx  => this%duidxjC(:,:,:,4); dvdy  => this%duidxjC(:,:,:,5); dvdzC => this%duidxjC(:,:,:,6); 
       dwdxC => this%duidxjC(:,:,:,7); dwdyC => this%duidxjC(:,:,:,8); dwdz  => this%duidxjC(:,:,:,9); 

       this%tidSUM = this%tidSUM + 1

       ! compute u*w on E, interpolate to C
       rbuff1E = this%uE * this%w
       call transpose_x_to_y(rbuff1E,rbuff2E,this%gpE)
       call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
       call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
       call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
       call transpose_y_to_x(rbuff2,rbuff1,this%gpC)

       ! compute v*w on E, interpolate to C
       rbuff1E = this%vE * this%w
       call transpose_x_to_y(rbuff1E,rbuff2E,this%gpE)
       call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
       call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
       call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
       call transpose_y_to_x(rbuff2,rbuff0,this%gpC)

       ! Compute u,v,wC - mean
       if(this%normByUstar) then
           this%u_mean3D = this%u_mean3D + this%u/this%sgsmodel%get_ustar()
           this%v_mean3D = this%v_mean3D + this%v/this%sgsmodel%get_ustar()
           this%w_mean3D = this%w_mean3D + this%wC/this%sgsmodel%get_ustar()

           this%uu_mean3D = this%uu_mean3D + this%u * this%u /this%sgsmodel%get_ustar()**2
           this%uv_mean3D = this%uv_mean3D + this%u * this%v /this%sgsmodel%get_ustar()**2
           this%uw_mean3D = this%uw_mean3D + rbuff1          /this%sgsmodel%get_ustar()**2
           this%vv_mean3D = this%vv_mean3D + this%v * this%v /this%sgsmodel%get_ustar()**2
           this%vw_mean3D = this%vw_mean3D + rbuff0          /this%sgsmodel%get_ustar()**2
           this%ww_mean3D = this%ww_mean3D + this%wC* this%wC/this%sgsmodel%get_ustar()**2
       else
           this%u_mean3D = this%u_mean3D + this%u
           this%v_mean3D = this%v_mean3D + this%v
           this%w_mean3D = this%w_mean3D + this%wC

           this%uu_mean3D = this%uu_mean3D + this%u * this%u
           this%uv_mean3D = this%uv_mean3D + this%u * this%v
           this%uw_mean3D = this%uw_mean3D + rbuff1
           this%vv_mean3D = this%vv_mean3D + this%v * this%v
           this%vw_mean3D = this%vw_mean3D + rbuff0
           this%ww_mean3D = this%ww_mean3D + this%wC* this%wC
       endif

       ! triple correlation for transport term in TKE budget -- triple product should be dealiased
       rbuff0 = this%u*this%u + this%v*this%v + this%wC*this%wC
       this%tketurbtranspx_mean3D = this%tketurbtranspx_mean3D + this%u *rbuff0
       this%tketurbtranspy_mean3D = this%tketurbtranspy_mean3D + this%v *rbuff0
       this%tketurbtranspz_mean3D = this%tketurbtranspz_mean3D + this%wC*rbuff0

       if(this%fastCalcPressure .or. this%storePressure) then
           if(this%normByUstar) then
               this%p_mean3D  = this%p_mean3D  + this%pressure          /this%sgsmodel%get_ustar()**2
               this%pu_mean3D = this%pu_mean3D + this%pressure * this%u /this%sgsmodel%get_ustar()**3
               this%pv_mean3D = this%pv_mean3D + this%pressure * this%v /this%sgsmodel%get_ustar()**3
               this%pw_mean3D = this%pw_mean3D + this%pressure * this%wC/this%sgsmodel%get_ustar()**3
           else
               this%p_mean3D  = this%p_mean3D  + this%pressure
               this%pu_mean3D = this%pu_mean3D + this%pressure * this%u
               this%pv_mean3D = this%pv_mean3D + this%pressure * this%v
               this%pw_mean3D = this%pw_mean3D + this%pressure * this%wC
           endif
       endif

       if(.not. this%isInviscid) then
           ! for viscous dissipation
           rbuff0 = dudx*dudx + dvdy*dvdy + dwdz*dwdz &
                  + half*( (dudy + dvdx)**2 + (dudzC + dwdxC)**2 + (dvdzC + dwdyC)**2) ! half here is two/four

           if(this%normByUstar) then
               ! viscous dissipation
               this%viscdisp_mean3D = this%viscdisp_mean3D + rbuff0/(this%sgsmodel%get_ustar()**2)

               ! viscous diffusion
               this%Siju1_mean3D = this%Siju1_mean3D + (dudx*this%u + &
                                                        half*(dudy  + dvdx )*this%v + &
                                                        half*(dudzC + dwdxC)*this%wC)/(this%sgsmodel%get_ustar()**2)

               this%Siju2_mean3D = this%Siju2_mean3D + (half*(dudy  + dvdx )*this%u + &
                                                        dvdy*this%v + &
                                                        half*(dvdzC + dwdyC)*this%wC)/(this%sgsmodel%get_ustar()**2)

               this%Siju3_mean3D = this%Siju3_mean3D + (half*(dudzC + dwdxC)*this%u + &
                                                        half*(dvdzC + dwdyC)*this%v + &
                                                        dwdz*this%wC)/(this%sgsmodel%get_ustar()**2)
           else
               ! viscous dissipation
               this%viscdisp_mean3D = this%viscdisp_mean3D + rbuff0

               ! viscous diffusion
               this%Siju1_mean3D = this%Siju1_mean3D + (dudx*this%u + &
                                                        half*(dudy  + dvdx )*this%v + &
                                                        half*(dudzC + dwdxC)*this%wC)

               this%Siju2_mean3D = this%Siju2_mean3D + (half*(dudy  + dvdx )*this%u + &
                                                        dvdy*this%v + &
                                                        half*(dvdzC + dwdyC)*this%wC)

               this%Siju3_mean3D = this%Siju3_mean3D + (half*(dudzC + dwdxC)*this%u + &
                                                        half*(dvdzC + dwdyC)*this%v + &
                                                        dwdz*this%wC)
           endif
       endif

       if(this%useSGS) then

           !write(300+nrank,'(i4,3(e19.12,1x))') 1, this%inst_horz_avg(2:3)
           call mpi_bcast(this%inst_horz_avg(2:3),2,mpirkind,0,mpi_comm_world,ierr)
           !write(300+nrank,'(i4,3(e19.12,1x))') 2, this%inst_horz_avg(2:3)


           ! interpolate tau13 from E to C
           call transpose_x_to_y(this%tau13,rbuff2E,this%gpE)
           call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE) 
           !if(nrank==0) then
           !    write(*,*) '-----------------'
           !    write(*,*) rbuff3E(1,1,1:2)
           !endif
           !write(200+nrank,*) 1, this%tsim, this%inst_horz_avg(2)
           !rbuff3E(:,:,1) = this%sgsmodel%tauijWMhat_in_Z(:,:,1,1)    !this%inst_horz_avg(2)     !--- =-(this%sgsmodel%get_ustar()**2) is correct only for Moeng's Wall Model, not for Bou-Zeid's model
           !if(nrank==0) then
           !    write(*,*) rbuff3E(1,1,1:2)
           !endif
           !this%debuginst(1) = p_sum(sum(rbuff3E(:,:,1)))*this%meanFact
           !write(200+nrank,*) 2, this%tsim, this%debuginst(1)
           !this%debuginst(2) = p_sum(sum(rbuff3E(:,:,2)))*this%meanFact
           !this%debuginst(3) = p_sum(sum(rbuff3E(:,:,3)))*this%meanFact
           call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
           !this%debuginst(4) = p_sum(sum(rbuff3(:,:,1)))*this%meanFact
           !this%debuginst(5) = p_sum(sum(rbuff3(:,:,2)))*this%meanFact
           !this%debugavg(:) = this%debugavg(:) + this%debuginst(:)
           !if(nrank==0) then
           !    write(nrank+100,'(11(e19.12,1x))') this%tsim, this%debuginst, this%debugavg/real(this%tidSUM, rkind)
           !    !write(*,*) '-----------------'
           !endif
           call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
           call transpose_y_to_x(rbuff2,this%tauSGS_ij(:,:,:,3),this%gpC)

           ! interpolate tau23 from E to C
           call transpose_x_to_y(this%tau23,rbuff2E,this%gpE)
           call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
           !rbuff3E(:,:,1) = this%sgsmodel%tauijWMhat_in_Z(:,:,1,2)    !this%inst_horz_avg(3)     !--- =-(this%sgsmodel%get_ustar()**2) is correct only for Moeng's Wall Model, not for Bou-Zeid's model
           call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
           call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
           call transpose_y_to_x(rbuff2,this%tauSGS_ij(:,:,:,5),this%gpC)

           ! sgs dissipation
           !rbuff1 = this%tauSGS_ij(:,:,:,1)*this%tauSGS_ij(:,:,:,1) + &
           !         this%tauSGS_ij(:,:,:,4)*this%tauSGS_ij(:,:,:,4) + &
           !         this%tauSGS_ij(:,:,:,6)*this%tauSGS_ij(:,:,:,6)
           !rbuff1 = rbuff1 + two*(this%tauSGS_ij(:,:,:,2)*this%tauSGS_ij(:,:,:,2) + &
           !                       this%tauSGS_ij(:,:,:,3)*this%tauSGS_ij(:,:,:,3) + &
           !                       this%tauSGS_ij(:,:,:,5)*this%tauSGS_ij(:,:,:,5) )
           !rbuff1 = rbuff1/(this%nu_SGS + 1.0d-14)         ! note: factor of half is in dump_stats
           rbuff1 = this%tauSGS_ij(:,:,:,1)*dudx + this%tauSGS_ij(:,:,:,4)*dvdy + this%tauSGS_ij(:,:,:,6)*dwdz &
                  + this%tauSGS_ij(:,:,:,2)*(dudy + dvdx) + this%tauSGS_ij(:,:,:,3)*(dudzC + dwdxC) &
                  + this%tauSGS_ij(:,:,:,5)*(dvdzC + dwdyC)

           if(this%normByUstar) then
               this%tau11_mean3D = this%tau11_mean3D + this%tauSGS_ij(:,:,:,1)/(this%sgsmodel%get_ustar()**2)
               this%tau12_mean3D = this%tau12_mean3D + this%tauSGS_ij(:,:,:,2)/(this%sgsmodel%get_ustar()**2)
               this%tau13_mean3D = this%tau13_mean3D + this%tauSGS_ij(:,:,:,3)/(this%sgsmodel%get_ustar()**2)
               this%tau22_mean3D = this%tau22_mean3D + this%tauSGS_ij(:,:,:,4)/(this%sgsmodel%get_ustar()**2)
               this%tau23_mean3D = this%tau23_mean3D + this%tauSGS_ij(:,:,:,5)/(this%sgsmodel%get_ustar()**2)
               this%tau33_mean3D = this%tau33_mean3D + this%tauSGS_ij(:,:,:,6)/(this%sgsmodel%get_ustar()**2)

               ! factor of H in normalization is missing from all statistics below
               this%sgsdissp_mean3D = this%sgsdissp_mean3D + rbuff1/(this%sgsmodel%get_ustar()**3)

               this%tauu1_mean3D = this%tauu1_mean3D + (this%tauSGS_ij(:,:,:,1)*this%u + &
                                                        this%tauSGS_ij(:,:,:,2)*this%v + &
                                                        this%tauSGS_ij(:,:,:,3)*this%wC)/(this%sgsmodel%get_ustar()**3)

               this%tauu2_mean3D = this%tauu2_mean3D + (this%tauSGS_ij(:,:,:,2)*this%u + &
                                                        this%tauSGS_ij(:,:,:,4)*this%v + &
                                                        this%tauSGS_ij(:,:,:,5)*this%wC)/(this%sgsmodel%get_ustar()**3)

               this%tauu3_mean3D = this%tauu3_mean3D + (this%tauSGS_ij(:,:,:,3)*this%u + &
                                                        this%tauSGS_ij(:,:,:,5)*this%v + &
                                                        this%tauSGS_ij(:,:,:,6)*this%wC)/(this%sgsmodel%get_ustar()**3)
           else
               this%tau11_mean3D = this%tau11_mean3D + this%tauSGS_ij(:,:,:,1)
               this%tau12_mean3D = this%tau12_mean3D + this%tauSGS_ij(:,:,:,2)
               this%tau13_mean3D = this%tau13_mean3D + this%tauSGS_ij(:,:,:,3)
               this%tau22_mean3D = this%tau22_mean3D + this%tauSGS_ij(:,:,:,4)
               this%tau23_mean3D = this%tau23_mean3D + this%tauSGS_ij(:,:,:,5)
               this%tau33_mean3D = this%tau33_mean3D + this%tauSGS_ij(:,:,:,6)

               this%sgsdissp_mean3D = this%sgsdissp_mean3D + rbuff1

               this%tauu1_mean3D = this%tauu1_mean3D + (this%tauSGS_ij(:,:,:,1)*this%u + &
                                                        this%tauSGS_ij(:,:,:,2)*this%v + &
                                                        this%tauSGS_ij(:,:,:,3)*this%wC)

               this%tauu2_mean3D = this%tauu2_mean3D + (this%tauSGS_ij(:,:,:,2)*this%u + &
                                                        this%tauSGS_ij(:,:,:,4)*this%v + &
                                                        this%tauSGS_ij(:,:,:,5)*this%wC)

               this%tauu3_mean3D = this%tauu3_mean3D + (this%tauSGS_ij(:,:,:,3)*this%u + &
                                                        this%tauSGS_ij(:,:,:,5)*this%v + &
                                                        this%tauSGS_ij(:,:,:,6)*this%wC)
           endif

       endif

       if(this%useWindTurbines) then
          if(this%normByUstar) then
             this%turbfx_mean3D = this%turbfx_mean3D + this%WindTurbineArr%fx/(this%sgsmodel%get_ustar()**3)
             this%turbfy_mean3D = this%turbfy_mean3D + this%WindTurbineArr%fy/(this%sgsmodel%get_ustar()**3)
             this%turbfz_mean3D = this%turbfz_mean3D + this%WindTurbineArr%fz/(this%sgsmodel%get_ustar()**3)
             this%uturbf_mean3D = this%uturbf_mean3D + (this%u *this%WindTurbineArr%fx + &
                                                        this%v *this%WindTurbineArr%fy + &
                                                        this%wC*this%WindTurbineArr%fz)/(this%sgsmodel%get_ustar()**3)
          else
             this%turbfx_mean3D = this%turbfx_mean3D + this%WindTurbineArr%fx
             this%turbfy_mean3D = this%turbfy_mean3D + this%WindTurbineArr%fy
             this%turbfz_mean3D = this%turbfz_mean3D + this%WindTurbineArr%fz
             this%uturbf_mean3D = this%uturbf_mean3D + this%u *this%WindTurbineArr%fx + &
                                                       this%v *this%WindTurbineArr%fy + &
                                                       this%wC*this%WindTurbineArr%fz
          endif
       endif

       if(this%isStratified) then
           ! compute T*w on E, interpolate to C
           rbuff1E = this%TE * this%w
           call transpose_x_to_y(rbuff1E,rbuff2E,this%gpE)
           call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
           call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
           call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
           call transpose_y_to_x(rbuff2,rbuff0,this%gpC)

           if(this%normByUstar) then
               this%T_mean3D = this%T_mean3D + this%T*this%sgsmodel%get_ustar()/this%wTh_surf
               this%uT_mean3D = this%uT_mean3D + this%T * this%u /this%wTh_surf
               this%vT_mean3D = this%vT_mean3D + this%T * this%v /this%wTh_surf
               this%wT_mean3D = this%wT_mean3D + rbuff0          /this%wTh_surf
               this%TT_mean3D = this%TT_mean3D + this%T * this%T*(this%sgsmodel%get_ustar()/this%wTh_surf)**2
           else
               this%T_mean3D = this%T_mean3D + this%T
               this%uT_mean3D = this%uT_mean3D + this%T * this%u
               this%vT_mean3D = this%vT_mean3D + this%T * this%v
               this%wT_mean3D = this%wT_mean3D + rbuff0
               this%TT_mean3D = this%TT_mean3D + this%T * this%T
           endif
           if(this%useSGS) then
               ! interpolate q3 from E to C
               call transpose_x_to_y(this%q3,rbuff2E,this%gpE)
               call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
               rbuff3E(:,:,1) = this%wTh_surf
               call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
               call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
               call transpose_y_to_x(rbuff2,rbuff1,this%gpC)

               if(this%normByUstar) then
                   this%q1_mean3D = this%q1_mean3D + this%q1/this%wTh_surf
                   this%q2_mean3D = this%q2_mean3D + this%q2/this%wTh_surf
                   this%q3_mean3D = this%q3_mean3D + rbuff1/this%wTh_surf
               else
                   this%q1_mean3D = this%q1_mean3D + this%q1
                   this%q2_mean3D = this%q2_mean3D + this%q2
                   this%q3_mean3D = this%q3_mean3D + rbuff1
               endif
           endif

       endif

       if (this%computeSpectra) then
           ! compute 1D spectra ---- make sure that number of variables for which spectra are computed is smaller than nyg
           ! For each variable, at each y, z, location, x-spectrum is computed first, and then averaged over time and y-direction
           jindx = 1    ! u
           call this%spectC%fft1_x2y(this%u,this%cbuffyC(:,:,:,1))
           do k = 1, size(this%cbuffyC, 3)
             do j = 1, size(this%cbuffyC, 2)
               this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%cbuffyC(:,j,k,1))
             end do
           end do

           jindx = 2    ! v
           call this%spectC%fft1_x2y(this%v,this%cbuffyC(:,:,:,1))
           do k = 1, size(this%cbuffyC, 3)
             do j = 1, size(this%cbuffyC, 2)
               this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%cbuffyC(:,j,k,1))
             end do
           end do
           
           jindx = 3    ! w
           call this%spectC%fft1_x2y(this%wC,this%cbuffyC(:,:,:,1))
           do k = 1, size(this%cbuffyC, 3)
             do j = 1, size(this%cbuffyC, 2)
               this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%cbuffyC(:,j,k,1))
             end do
           end do
           
           jindx = 4    ! KE
           call this%spectC%fft1_x2y(half*(this%u**2+this%v**2+this%wC**2),this%cbuffyC(:,:,:,1))
           do k = 1, size(this%cbuffyC, 3)
             do j = 1, size(this%cbuffyC, 2)
               this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%cbuffyC(:,j,k,1))
             end do
           end do
           
           if(this%isStratified) then
               jindx = jindx + 1    ! T
               call this%spectC%fft1_x2y(this%T,this%cbuffyC(:,:,:,1))
               do k = 1, size(this%cbuffyC, 3)
                 do j = 1, size(this%cbuffyC, 2)
                   this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%cbuffyC(:,j,k,1))
                 end do
               end do
           endif
           
           if(this%fastCalcPressure .or. this%storePressure) then
               jindx = jindx + 1    ! p
               call this%spectC%fft1_x2y(this%pressure,this%cbuffyC(:,:,:,1))
               do k = 1, size(this%cbuffyC, 3)
                 do j = 1, size(this%cbuffyC, 2)
                   this%xspectra_mean(:,jindx,k) = this%xspectra_mean(:,jindx,k) + abs(this%cbuffyC(:,j,k,1))
                 end do
               end do
           endif
       end if 

       ! instantaneous horizontal averages of some quantities
       this%inst_horz_avg(1) = this%sgsmodel%get_ustar()
       ! this%inst_horz(2) and (3) are computed on nrank==0 proc in getRHS_SGS_WallM
       ! broadcast to all other procs above in this subroutine
       ! do nothing about inst_horz_avg(2:3) herre
       if(this%isStratified) then
           this%inst_horz_avg(4) = this%invObLength
           this%inst_horz_avg(5) = this%wTh_surf
       endif
       ! this%inst_horz_avg_turb(1:5*this%WindTurbineArr%nTurbines) is computed in this%WindTurbineArr%getForceRHS
       this%runningSum_sc = this%runningSum_sc + this%inst_horz_avg
           !write(200+nrank,*) 3, this%tsim, this%runningSum_sc(2)
       if(this%useWindTurbines) this%runningSum_sc_turb = this%runningSum_sc_turb + this%inst_horz_avg_turb

       nullify(rbuff0,rbuff1,rbuff2,rbuff3,rbuff2E,rbuff3E,rbuff4,rbuff1E)
       nullify(dudx, dudy, dudzC, dvdx, dvdy, dvdzC, dwdxC, dwdyC, dwdz)

   end subroutine

   subroutine Delete_file_if_present(this, tempname)
       class(igrid), intent(inout) :: this
       character(len=clen), intent(in) :: tempname
       character(len=clen) :: fname

       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       open(unit=11, file=fname, status='old',iostat=ierr); if(ierr == 0) close(11, status='delete')

   end subroutine

   subroutine DeletePrevStats3DFiles(this)
       class(igrid), intent(inout) :: this
       character(len=clen) :: tempname

       if(nrank==0) then
         ! delete stats files corresponding to tprev2
         write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_um_t",this%tprev2,".3Dstt";  call this%Delete_file_if_present(tempname)
         write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_vm_t",this%tprev2,".3Dstt";  call this%Delete_file_if_present(tempname)
         write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_wm_t",this%tprev2,".3Dstt";  call this%Delete_file_if_present(tempname)
         write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uum_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
         write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uvm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
         write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uwm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
         write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_vvm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
         write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_vwm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
         write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_wwm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
         write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tktt_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
         write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tkma_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
         write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tkpr_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
         write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_mktt_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
         write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_mkma_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
         if(this%fastCalcPressure .or. this%storePressure) then
             write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_pm_t",this%tprev2,".3Dstt";   call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_mkpt_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tkpt_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
         endif
         if(this%useSGS) then
             write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t11_t",this%tprev2,".3Dstt";    call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t12_t",this%tprev2,".3Dstt";    call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t13_t",this%tprev2,".3Dstt";    call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t22_t",this%tprev2,".3Dstt";    call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t23_t",this%tprev2,".3Dstt";    call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t33_t",this%tprev2,".3Dstt";    call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_mksgsd_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_tksgsd_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_mksgst_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_tksgst_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
         endif

         if(this%useWindTurbines) then
             write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_trbx_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_trby_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_trbz_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_mktrbf_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_tktrbf_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
         endif

         if(this%isStratified) then
             write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_Tm_t",this%tprev2,".3Dstt";  call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uTm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_vTm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_wTm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
             write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_TTm_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
             if(this%useSGS) then
               write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_q1_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
               write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_q2_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
               write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_q3_t",this%tprev2,".3Dstt"; call this%Delete_file_if_present(tempname)
             endif
         endif
       endif
   end subroutine

   subroutine dump_stats3D(this)
       use basic_io, only: write_2d_ascii
       use decomp_2d_io
       use kind_parameters, only: mpirkind
       class(igrid), intent(inout), target :: this
     ! compute horizontal averages and dump .stt files
     ! overwrite previously written out 3D stats dump
       real(rkind), dimension(:,:,:), pointer :: rbuff0, rbuff1, rbuff2, rbuff3, rbuff4, rbuff5, rbuff6, rbuff3E, rbuff4E
       real(rkind), dimension(:,:,:), pointer :: S11_mean3D, S12_mean3D, S13_mean3D, S22_mean3D, S23_mean3D, S33_mean3D
       complex(rkind), dimension(:,:,:), pointer :: cbuffy1, cbuffy2
       real(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)) :: tmpvar
       real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),6), target :: Stmp
       character(len=clen) :: tempname, fname
       real(rkind) :: tidSUMreal, normfac
       integer :: tid, dirid, decompdir, jindx, nspectra

       tid = this%step

       ! Ensure only two sets of 3Dstats files are kept
       if(this%tprev2 > 0) then
         call this%DeletePrevStats3DFiles()
       endif
       ! Now update tprev2 and tprev1 
       this%tprev2 = this%tprev1
       this%tprev1 = this%step

       rbuff0 => this%rbuffxC(:,:,:,2);
       rbuff1 => this%rbuffxC(:,:,:,1);  rbuff2 => this%rbuffyC(:,:,:,1)
       rbuff3 => this%rbuffzC(:,:,:,1);  rbuff4 => this%rbuffzC(:,:,:,2)
       rbuff5 => this%rbuffzC(:,:,:,3);  rbuff6 => this%rbuffzC(:,:,:,4)

       cbuffy1 => this%cbuffyC(:,:,:,1); cbuffy2 => this%cbuffyC(:,:,:,2)
       rbuff3E => this%rbuffzE(:,:,:,1); rbuff4E => this%rbuffzE(:,:,:,2);

       tidSUMreal = real(this%tidSUM, rkind)

       ! u_avg
       call transpose_x_to_y(this%u_mean3D/tidSUMreal, rbuff2, this%gpC)
       call transpose_y_to_z(rbuff2,                   rbuff3, this%gpC)
       call this%compute_z_mean(rbuff3, this%u_mean)
       write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_um_t",this%step,".3Dstt"
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(3, rbuff3, fname)

       ! v_avg
       call transpose_x_to_y(this%v_mean3D/tidSUMreal, rbuff2, this%gpC)
       call transpose_y_to_z(rbuff2,                   rbuff4, this%gpC)
       call this%compute_z_mean(rbuff4, this%v_mean)
       write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_vm_t",this%step,".3Dstt"
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(3, rbuff4, fname)
     
       ! w_avg
       call transpose_x_to_y(this%w_mean3D/tidSUMreal, rbuff2, this%gpC)
       call transpose_y_to_z(rbuff2,                   rbuff5, this%gpC)
       call this%compute_z_mean(rbuff5, this%w_mean)
       write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_wm_t",this%step,".3Dstt"
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(3, rbuff5, fname)
     
       ! uu_avg - u_avg*u_avg - Total (1D) and Reynolds (3D)
       call transpose_x_to_y(this%uu_mean3D/tidSUMreal, rbuff2, this%gpC)
       call transpose_y_to_z(rbuff2,                    rbuff6, this%gpC)
       call this%compute_z_mean(rbuff6, this%uu_mean)
       this%uu_mean = this%uu_mean - this%u_mean*this%u_mean    !--Total
       rbuff6 = rbuff6 - rbuff3 * rbuff3                        !--Reynolds
       write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uum_t",this%step,".3Dstt"
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(3, rbuff6, fname)
     
       ! uv_avg - u_avg*v_avg - Total (1D) and Reynolds (3D)
       call transpose_x_to_y(this%uv_mean3D/tidSUMreal, rbuff2, this%gpC)
       call transpose_y_to_z(rbuff2,                    rbuff6, this%gpC)
       call this%compute_z_mean(rbuff6, this%uv_mean)
       this%uv_mean = this%uv_mean - this%u_mean*this%v_mean    !--Total
       rbuff6 = rbuff6 - rbuff3 * rbuff4                        !--Reynolds
       write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uvm_t",this%step,".3Dstt"
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(3, rbuff6, fname)
     
       ! uw_avg - u_avg*w_avg - Reynolds
       call transpose_x_to_y(this%uw_mean3D/tidSUMreal, rbuff2, this%gpC)
       call transpose_y_to_z(rbuff2,                    rbuff6, this%gpC)
       rbuff6 = rbuff6 - rbuff3*rbuff5
       call this%compute_z_mean(rbuff6, this%uw_mean)
       write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uwm_t",this%step,".3Dstt"
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(3, rbuff6, fname)
     
       ! uw_avg - u_avg*w_avg - Dispersive
       rbuff6 = rbuff3*rbuff5
       call this%compute_z_mean(rbuff6, this%disperuw_mean)
       this%disperuw_mean = this%disperuw_mean - this%u_mean*this%w_mean
     
       ! vv_avg - v_avg*v_avg - Total (1D) and Reynolds (3D)
       call transpose_x_to_y(this%vv_mean3D/tidSUMreal, rbuff2, this%gpC)
       call transpose_y_to_z(rbuff2,                    rbuff6, this%gpC)
       call this%compute_z_mean(rbuff6, this%vv_mean)
       this%vv_mean = this%vv_mean - this%v_mean*this%v_mean    !--Total
       rbuff6 = rbuff6 - rbuff4 * rbuff4                        !--Reynolds
       write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_vvm_t",this%step,".3Dstt"
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(3, rbuff6, fname)
     
       ! vw_avg - v_avg*w_avg - Reynolds
       call transpose_x_to_y(this%vw_mean3D/tidSUMreal, rbuff2, this%gpC)
       call transpose_y_to_z(rbuff2,                    rbuff6, this%gpC)
       rbuff6 = rbuff6 - rbuff4*rbuff5
       call this%compute_z_mean(rbuff6, this%vw_mean)
       write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_vwm_t",this%step,".3Dstt"
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(3, rbuff6, fname)
     
       ! vw_avg - v_avg*w_avg - Dispersive
       rbuff6 = rbuff4*rbuff5
       call this%compute_z_mean(rbuff6, this%dispervw_mean)
       this%dispervw_mean = this%dispervw_mean - this%v_mean*this%w_mean
     
       ! ww_avg - w_avg*w_avg - Total (1D) and Reynolds (3D)
       call transpose_x_to_y(this%ww_mean3D/tidSUMreal, rbuff2, this%gpC)
       call transpose_y_to_z(rbuff2,                    rbuff6, this%gpC)
       call this%compute_z_mean(rbuff6, this%ww_mean)
       this%ww_mean = this%ww_mean - this%w_mean*this%w_mean    !--Total
       rbuff6 = rbuff6 - rbuff5 * rbuff5                        !--Reynolds
       write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_wwm_t",this%step,".3Dstt"
       fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
       call decomp_2d_write_one(3, rbuff6, fname)
     
       !-----Turbulent transport of TKE budget-------
         ! triple correlation for transport term in TKE budget
         ! x term in xdecomp
         !rbuff1 = this%u_mean3D*this%u_mean3D + this%v_mean3D*this%v_mean3D + this%w_mean3D*this%w_mean3D
         rbuff1 = this%tketurbtranspx_mean3D/tidSUMreal - two*(this%u_mean3D*this%uu_mean3D + this%v_mean3D*this%uv_mean3D + this%w_mean3D*this%uw_mean3D)/(tidSumreal**2) &
                + ( two*(this%u_mean3D*this%u_mean3D + this%v_mean3D*this%v_mean3D + this%w_mean3D*this%w_mean3D)/tidSUMreal**2 - (this%uu_mean3D + this%vv_mean3D + this%ww_mean3D)/tidSUMreal )*(this%u_mean3D/tidSUMreal)
         call this%spectC%fft(rbuff1,cbuffy1)   
         call this%spectC%mTimes_ik1_ip(cbuffy1)
         call this%spectC%ifft(cbuffy1,rbuff0)

         ! add y term in x decomp
         rbuff1 = this%tketurbtranspy_mean3D/tidSUMreal - two*(this%u_mean3D*this%uv_mean3D + this%v_mean3D*this%vv_mean3D + this%w_mean3D*this%vw_mean3D)/(tidSumreal**2) &
                + ( two*(this%u_mean3D*this%u_mean3D + this%v_mean3D*this%v_mean3D + this%w_mean3D*this%w_mean3D)/tidSUMreal**2 - (this%uu_mean3D + this%vv_mean3D + this%ww_mean3D)/tidSUMreal )*(this%v_mean3D/tidSUMreal)
         call this%spectC%fft(rbuff1,cbuffy1)   
         call this%spectC%mTimes_ik2_ip(cbuffy1)
         call this%spectC%ifft(cbuffy1,rbuff1)
         rbuff0 = rbuff0 + rbuff1

         ! take sum of x and y terms to z decomp
         call transpose_x_to_y(rbuff0, rbuff2, this%gpC)
         call transpose_y_to_z(rbuff2, rbuff4, this%gpC)

         ! compute z term in z decomp
         rbuff1 = this%tketurbtranspz_mean3D/tidSUMreal - two*(this%u_mean3D*this%uw_mean3D + this%v_mean3D*this%vw_mean3D + this%w_mean3D*this%ww_mean3D)/(tidSumreal**2) &
                + ( two*(this%u_mean3D*this%u_mean3D + this%v_mean3D*this%v_mean3D + this%w_mean3D*this%w_mean3D)/tidSUMreal**2 - (this%uu_mean3D + this%vv_mean3D + this%ww_mean3D)/tidSUMreal )*(this%w_mean3D/tidSUMreal)
         call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
         call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
         ! interpolate rbuff3 from C to E
         call this%Pade6opZ%interpz_C2E(rbuff3, rbuff3E, 0,0)
         call this%Pade6opZ%ddz_E2C(rbuff3E,rbuff3,0,0)

         ! add x and y terms to z term
         rbuff3 = rbuff3 + rbuff4
         rbuff3 = -half*rbuff3

         call this%compute_z_mean(rbuff3, this%tkett_mean)
         write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tktt_t",this%step,".3Dstt"
         fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
         call decomp_2d_write_one(3, rbuff3, fname)
       !-----Done turbulent transport of TKE budget-------

       !-----Mean advection term of TKE budget------
         ! compute tke first
         rbuff1 = this%u_mean3D*this%u_mean3D + this%v_mean3D*this%v_mean3D + this%w_mean3D*this%w_mean3D
         rbuff1 = -rbuff1/tidSUMreal**2
         rbuff1 = rbuff1 + (this%uu_mean3D + this%vv_mean3D + this%ww_mean3D)/tidSUMreal

         ! transpose to z and take z derivative
         call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
         call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
         call this%Pade6opZ%ddz_C2E(rbuff3,rbuff3E,0,0)
         
         ! transpose w_mean3D to z, interpolate to E and multiply
         call transpose_x_to_y(this%w_mean3D/tidSUMreal, rbuff2, this%gpC)
         call transpose_y_to_z(rbuff2,                   rbuff4, this%gpC)
         call this%Pade6opZ%interpz_C2E(rbuff4,rbuff4E,0,0)
         rbuff4E = rbuff4E * rbuff3E

         ! interpolate E to C
         call this%Pade6opZ%interpz_E2C(rbuff4E, rbuff4, 0,0)


         ! x derivative
         call this%spectC%fft(rbuff1,cbuffy1)   
         call this%spectC%mTimes_ik1_oop(cbuffy1, cbuffy2)
         call this%spectC%ifft(cbuffy2,rbuff0)

         ! y derivative
         call this%spectC%fft(rbuff1,cbuffy1)   
         call this%spectC%mTimes_ik2_oop(cbuffy1, cbuffy2)
         call this%spectC%ifft(cbuffy2,rbuff1)
         rbuff0 = this%u_mean3D*rbuff0 + this%v_mean3D*rbuff1

         ! transpose sum of x and y parts to z and add z part
         call transpose_x_to_y(rbuff0/tidSUMreal, rbuff2, this%gpC)
         call transpose_y_to_z(rbuff2,            rbuff3, this%gpC)
         rbuff3 = rbuff3 + rbuff4
         rbuff3 = -half*rbuff3

         ! write outputs
         call this%compute_z_mean(rbuff3, this%tkeadv_mean)
         write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tkma_t",this%step,".3Dstt"
         fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
         call decomp_2d_write_one(3, rbuff3, fname)
       !-----Done mean advection term of TKE budget------

       S11_mean3D => Stmp(:,:,:,1);   S12_mean3D => Stmp(:,:,:,2);    S13_mean3D => Stmp(:,:,:,3) 
                                      S22_mean3D => Stmp(:,:,:,4);    S23_mean3D => Stmp(:,:,:,5) 
                                                                      S33_mean3D => Stmp(:,:,:,6) 
       call this%compute_Sijmean(Stmp)

       ! ---- Shear production in TKE budget equation-------
         rbuff1 = - (this%uu_mean3D/tidSUMreal - this%u_mean3D*this%u_mean3D/tidSUMreal**2)*S11_mean3D &
                  - (this%vv_mean3D/tidSUMreal - this%v_mean3D*this%v_mean3D/tidSUMreal**2)*S22_mean3D &
                  - (this%ww_mean3D/tidSUMreal - this%w_mean3D*this%w_mean3D/tidSUMreal**2)*S33_mean3D
         rbuff1 = rbuff1 - two*( &
                  + (this%uv_mean3D/tidSUMreal - this%u_mean3D*this%v_mean3D/tidSUMreal**2)*S12_mean3D &
                  + (this%uw_mean3D/tidSUMreal - this%u_mean3D*this%w_mean3D/tidSUMreal**2)*S13_mean3D &
                  + (this%vw_mean3D/tidSUMreal - this%v_mean3D*this%w_mean3D/tidSUMreal**2)*S23_mean3D )

         call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
         call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
         call this%compute_z_mean(rbuff3, this%tkeprod_mean)
         write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tkpr_t",this%step,".3Dstt"
         fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
         call decomp_2d_write_one(3, rbuff3, fname)
       ! ---- Done Shear production in TKE budget equation-------
     
       ! ---- Turbulent dissipation in MKE budget equation-------
         ! --- Identical to negative of shear production in TKE equation---
         this%mkedisp_mean = -this%tkeprod_mean
       ! ---- Done Turbulent dissipation in MKE budget equation-------
       
       !-----Transport of turbulent stresses in MKE budget-------
         ! x term in xdecomp
         rbuff1 = - (this%uu_mean3D/tidSUMreal - this%u_mean3D*this%u_mean3D/tidSUMreal**2)*this%u_mean3D/tidSUMreal &
                  - (this%uv_mean3D/tidSUMreal - this%u_mean3D*this%v_mean3D/tidSUMreal**2)*this%v_mean3D/tidSUMreal &
                  - (this%uw_mean3D/tidSUMreal - this%u_mean3D*this%w_mean3D/tidSUMreal**2)*this%w_mean3D/tidSUMreal
         call this%spectC%fft(rbuff1,cbuffy1)   
         call this%spectC%mTimes_ik1_ip(cbuffy1)
         call this%spectC%ifft(cbuffy1,rbuff0)

         ! y term in xdecomp
         rbuff1 = - (this%uv_mean3D/tidSUMreal - this%u_mean3D*this%v_mean3D/tidSUMreal**2)*this%u_mean3D/tidSUMreal &
                  - (this%vv_mean3D/tidSUMreal - this%v_mean3D*this%v_mean3D/tidSUMreal**2)*this%v_mean3D/tidSUMreal &
                  - (this%vw_mean3D/tidSUMreal - this%v_mean3D*this%w_mean3D/tidSUMreal**2)*this%w_mean3D/tidSUMreal
         call this%spectC%fft(rbuff1,cbuffy1)   
         call this%spectC%mTimes_ik2_oop(cbuffy1, cbuffy2)
         call this%spectC%ifft(cbuffy2,rbuff1)

         ! transpose sum of x and y parts to z
         call transpose_x_to_y(rbuff0+rbuff1, rbuff2, this%gpC)
         call transpose_y_to_z(rbuff2,        rbuff4, this%gpC)

         ! compute z term in z decomp
         rbuff1 = - (this%uw_mean3D/tidSUMreal - this%u_mean3D*this%w_mean3D/tidSUMreal**2)*this%u_mean3D/tidSUMreal &
                  - (this%vw_mean3D/tidSUMreal - this%v_mean3D*this%w_mean3D/tidSUMreal**2)*this%v_mean3D/tidSUMreal &
                  - (this%ww_mean3D/tidSUMreal - this%w_mean3D*this%w_mean3D/tidSUMreal**2)*this%w_mean3D/tidSUMreal
         call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
         call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
         ! interpolate rbuff3 from C to E
         call this%Pade6opZ%interpz_C2E(rbuff3, rbuff3E, 0,0)
         call this%Pade6opZ%ddz_E2C(rbuff3E,rbuff3,0,0)

         ! add x and y terms to z term
         rbuff3 = rbuff3 + rbuff4

         call this%compute_z_mean(rbuff3, this%mkett_mean)
         write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_mktt_t",this%step,".3Dstt"
         fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
         call decomp_2d_write_one(3, rbuff3, fname)
       !-----Done Transport of turbulent stresses in MKE budget-------

       !-----Mean advection term of MKE budget------
         ! compute mke first
         rbuff1 = this%u_mean3D*this%u_mean3D + this%v_mean3D*this%v_mean3D + this%w_mean3D*this%w_mean3D
         rbuff1 = rbuff1/tidSUMreal**2

         ! transpose to z and take z derivative
         call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
         call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
         call this%Pade6opZ%ddz_C2E(rbuff3,rbuff3E,0,0)
         
         ! transpose w_mean3D to z, interpolate to E and multiply
         call transpose_x_to_y(this%w_mean3D/tidSUMreal, rbuff2, this%gpC)
         call transpose_y_to_z(rbuff2,                   rbuff4, this%gpC)
         call this%Pade6opZ%interpz_C2E(rbuff4,rbuff4E,0,0)
         rbuff4E = rbuff4E * rbuff3E

         ! interpolate E to C
         call this%Pade6opZ%interpz_E2C(rbuff4E, rbuff4, 0,0)


         ! x derivative
         call this%spectC%fft(rbuff1,cbuffy1)   
         call this%spectC%mTimes_ik1_oop(cbuffy1, cbuffy2)
         call this%spectC%ifft(cbuffy2,rbuff0)

         ! y derivative
         call this%spectC%fft(rbuff1,cbuffy1)   
         call this%spectC%mTimes_ik2_oop(cbuffy1, cbuffy2)
         call this%spectC%ifft(cbuffy2,rbuff1)
         rbuff0 = this%u_mean3D*rbuff0 + this%v_mean3D*rbuff1

         ! transpose sum of x and y parts to z and add z part
         call transpose_x_to_y(rbuff0/tidSUMreal, rbuff2, this%gpC)
         call transpose_y_to_z(rbuff2,            rbuff3, this%gpC)
         rbuff3 = rbuff3 + rbuff4
         rbuff3 = -half*rbuff3

         ! write outputs
         call this%compute_z_mean(rbuff3, this%mkeadv_mean)
         write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_mkma_t",this%step,".3Dstt"
         fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
         call decomp_2d_write_one(3, rbuff3, fname)
       !-----Done mean advection term of MKE budget------

       if(this%fastCalcPressure .or. this%storePressure) then
           ! p_avg
           call transpose_x_to_y(this%p_mean3D/tidSUMreal, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2,                   rbuff6, this%gpC)
           call this%compute_z_mean(rbuff6, this%p_mean)
           write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_pm_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff6, fname)

           ! compute ddxj(p_avg * uavg_j) -- Pressure transport term in mean KE eqn
             ! ddx(p_avg*uavg_1) in x decomp
             rbuff1 = this%p_mean3D*this%u_mean3D/(tidSUMreal**2)
             call this%spectC%fft(rbuff1,cbuffy1)   
             call this%spectC%mTimes_ik1_oop(cbuffy1, cbuffy2)
             call this%spectC%ifft(cbuffy2,rbuff0)

             ! ddy(p_avg*uavg_2) in x decomp
             rbuff1 = this%p_mean3D*this%v_mean3D/(tidSUMreal**2)
             call this%spectC%fft(rbuff1,cbuffy1)   
             call this%spectC%mTimes_ik2_oop(cbuffy1, cbuffy2)
             call this%spectC%ifft(cbuffy2,rbuff1)
             rbuff0 = rbuff0 + rbuff1

             ! take sum of x and y terms to z decomp
             call transpose_x_to_y(rbuff0, rbuff2, this%gpC)
             call transpose_y_to_z(rbuff2, rbuff4, this%gpC)

             ! ddz(p_avg*uavg_3) in z decomp
             rbuff1 = this%p_mean3D*this%w_mean3D/(tidSUMreal**2)
             ! transpose to z and take z derivative
             call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
             call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
             call this%Pade6opZ%ddz_C2E(rbuff3,rbuff3E,0,0)
             call this%Pade6opZ%interpz_E2C(rbuff3E, rbuff3, 0,0)       

             rbuff3 = -(rbuff3 + rbuff4)

           call this%compute_z_mean(rbuff3, this%mkept_mean)       ! Pressure transport mean (term in mean KE eqn)
           write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_mkpt_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
           ! Done computing ddxj(p_avg * uavg_j) -- Pressure transport term in mean KE eqn

     
           ! compute ddxj(p' * u'_j) -- Pressure transport term in turbulent KE eqn
             ! pu_avg - u_avg*p_avg
             rbuff1 = this%pu_mean3D/tidSUMreal - this%u_mean3D * this%p_mean3D / tidSUMreal**2
             call this%spectC%fft(rbuff1,cbuffy1)   
             call this%spectC%mTimes_ik1_oop(cbuffy1, cbuffy2)
             call this%spectC%ifft(cbuffy2,rbuff0)

             ! pv_avg - v_avg*p_avg - Reynolds only
             rbuff1 = this%pv_mean3D/tidSUMreal - this%v_mean3D * this%p_mean3D / tidSUMreal**2
             call this%spectC%fft(rbuff1,cbuffy1)   
             call this%spectC%mTimes_ik2_oop(cbuffy1, cbuffy2)
             call this%spectC%ifft(cbuffy2,rbuff1)
             rbuff0 = rbuff0 + rbuff1

             ! take sum of x and y terms to z decomp
             call transpose_x_to_y(rbuff0, rbuff2, this%gpC)
             call transpose_y_to_z(rbuff2, rbuff4, this%gpC)

             ! pw_avg - w_avg*p_avg - Reynolds only
             rbuff1 = this%pw_mean3D/tidSUMreal - this%w_mean3D * this%p_mean3D / tidSUMreal**2
             call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
             call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
             call this%Pade6opZ%ddz_C2E(rbuff3,rbuff3E,0,0)
             call this%Pade6opZ%interpz_E2C(rbuff3E, rbuff3, 0,0)       

             rbuff3 = -(rbuff3 + rbuff4)

           call this%compute_z_mean(rbuff3, this%tkept_mean)       ! Pressure transport turbulent (term in turbulent KE eqn)
           write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_tkpt_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
           ! Done computing ddxj(p' * u'_j) -- Pressure transport term in turbulent KE eqn

       endif

       if(.not. this%isInviscid) then
          ! mkevdif_mean; mkevdsp_mean; tkevdif_mean; tkevdsp_mean to be written
       endif

       if(this%useSGS) then
           ! tau11SGS_avg
           call transpose_x_to_y(this%tau11_mean3D/tidSUMreal, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2,                       rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%tau11_mean)
           write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t11_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
     
           ! tau12SGS_avg
           call transpose_x_to_y(this%tau12_mean3D/tidSUMreal, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2,                       rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%tau12_mean)
           write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t12_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
     
           ! tau13SGS_avg
           call transpose_x_to_y(this%tau13_mean3D/tidSUMreal, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2,                       rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%tau13_mean)
           write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t13_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
     
           ! tau22SGS_avg
           call transpose_x_to_y(this%tau22_mean3D/tidSUMreal, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2,                       rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%tau22_mean)
           write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t22_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
     
           ! tau23SGS_avg
           call transpose_x_to_y(this%tau23_mean3D/tidSUMreal, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2,                       rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%tau23_mean)
           write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t23_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
     
           ! tau33SGS_avg
           call transpose_x_to_y(this%tau33_mean3D/tidSUMreal, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2,                       rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%tau33_mean)
           write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_t33_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
     
           ! SGS dissipation in MKE equation
           rbuff1 = (    this%tau11_mean3D*S11_mean3D + this%tau22_mean3D*S22_mean3D + this%tau33_mean3D*S33_mean3D + & 
                    two*(this%tau12_mean3D*S12_mean3D + this%tau13_mean3D*S13_mean3D + this%tau23_mean3D*S23_mean3D)  )/tidSUMreal
           call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%mkesgsd_mean)
           write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_mksgsd_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)

           ! SGS dissipation in TKE equation
           rbuff1 = this%sgsdissp_mean3D/tidSUMreal - rbuff1
           call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%tkesgsd_mean)
           write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_tksgsd_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
     
           ! Mean transport of SGS stress in MKE equation
             ! x term in xdecomp
             rbuff1 = (this%u_mean3D*this%tau11_mean3D + this%v_mean3D*this%tau12_mean3D + this%w_mean3D*this%tau13_mean3D)/tidSUMreal**2
             call this%spectC%fft(rbuff1,cbuffy1)   
             call this%spectC%mTimes_ik1_oop(cbuffy1, cbuffy2)
             call this%spectC%ifft(cbuffy2,rbuff0)

             ! y term in xdecomp
             rbuff1 = (this%u_mean3D*this%tau12_mean3D + this%v_mean3D*this%tau22_mean3D + this%w_mean3D*this%tau23_mean3D)/tidSUMreal**2
             call this%spectC%fft(rbuff1,cbuffy1)   
             call this%spectC%mTimes_ik2_oop(cbuffy1, cbuffy2)
             call this%spectC%ifft(cbuffy2,rbuff1)
             rbuff0 = rbuff0 + rbuff1

             ! take sum of x and y terms to z decomp
             call transpose_x_to_y(rbuff0, rbuff2, this%gpC)
             call transpose_y_to_z(rbuff2, rbuff4, this%gpC)

             ! z term in zdecomp
             rbuff1 = (this%u_mean3D*this%tau13_mean3D + this%v_mean3D*this%tau23_mean3D + this%w_mean3D*this%tau33_mean3D)/tidSUMreal**2
             call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
             call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
             call this%Pade6opZ%ddz_C2E(rbuff3,rbuff3E,0,0)
             call this%Pade6opZ%interpz_E2C(rbuff3E, rbuff3, 0,0)       

             rbuff3 = -(rbuff3 + rbuff4)

             call this%compute_z_mean(rbuff3, this%mkesgst_mean)       ! Pressure transport turbulent (term in turbulent KE eqn)
             write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_mksgst_t",this%step,".3Dstt"
             fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
             call decomp_2d_write_one(3, rbuff3, fname)
           ! Done mean transport of SGS stress in MKE equation

           ! Turbulent transport of SGS stress in TKE equation
             ! x term in xdecomp
             rbuff1 = this%tauu1_mean3D/tidSUMreal - (this%u_mean3D*this%tau11_mean3D + this%v_mean3D*this%tau12_mean3D + this%w_mean3D*this%tau13_mean3D)/tidSUMreal**2
             call this%spectC%fft(rbuff1,cbuffy1)   
             call this%spectC%mTimes_ik1_oop(cbuffy1, cbuffy2)
             call this%spectC%ifft(cbuffy2,rbuff0)

             ! y term in xdecomp
             rbuff1 = this%tauu2_mean3D/tidSUMreal - (this%u_mean3D*this%tau12_mean3D + this%v_mean3D*this%tau22_mean3D + this%w_mean3D*this%tau23_mean3D)/tidSUMreal**2
             call this%spectC%fft(rbuff1,cbuffy1)   
             call this%spectC%mTimes_ik2_oop(cbuffy1, cbuffy2)
             call this%spectC%ifft(cbuffy2,rbuff1)
             rbuff0 = rbuff0 + rbuff1

             ! take sum of x and y terms to z decomp
             call transpose_x_to_y(rbuff0, rbuff2, this%gpC)
             call transpose_y_to_z(rbuff2, rbuff4, this%gpC)

             ! z term in zdecomp
             rbuff1 = this%tauu3_mean3D/tidSUMreal - (this%u_mean3D*this%tau13_mean3D + this%v_mean3D*this%tau23_mean3D + this%w_mean3D*this%tau33_mean3D)/tidSUMreal**2
             call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
             call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
             call this%Pade6opZ%ddz_C2E(rbuff3,rbuff3E,0,0)
             call this%Pade6opZ%interpz_E2C(rbuff3E, rbuff3, 0,0)       

             rbuff3 = -(rbuff3 + rbuff4)

             call this%compute_z_mean(rbuff3, this%tkesgst_mean)       ! Pressure transport turbulent (term in turbulent KE eqn)
             write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_tksgst_t",this%step,".3Dstt"
             fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
             call decomp_2d_write_one(3, rbuff3, fname)
           ! Done turbulent transport of SGS stress in TKE equation
       endif

       if(this%useWindTurbines) then
          !----Turbine work term in MKE budget-------
          rbuff1 = ( this%u_mean3D*this%turbfx_mean3D + this%v_mean3D*this%turbfy_mean3D + this%w_mean3D*this%turbfy_mean3D ) / tidSUMreal**2
          call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
          call this%compute_z_mean(rbuff3, this%mketurbf_mean)
          write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_mktrbf_t",this%step,".3Dstt"
          fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
          call decomp_2d_write_one(3, rbuff3, fname)
          !----Done turbine work term in MKE budget-------

          !----Turbine work term in TKE budget-------
          rbuff1 = this%uturbf_mean3D/tidSUMreal -  ( this%u_mean3D*this%turbfx_mean3D + &
                   this%v_mean3D*this%turbfy_mean3D + this%w_mean3D*this%turbfy_mean3D ) / tidSUMreal**2
          call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
          call this%compute_z_mean(rbuff3, this%tketurbf_mean)
          write(tempname,"(A3,I2.2,A9,I6.6,A6)") "Run",this%runID, "_tktrbf_t",this%step,".3Dstt"
          fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
          call decomp_2d_write_one(3, rbuff3, fname)
          !----Done turbine work term in TKE budget-------

          ! write out turbfx_mean3D
          call transpose_x_to_y(this%turbfx_mean3D/tidSUMreal, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2,                        rbuff3, this%gpC)
          call this%compute_z_mean(rbuff3, this%turbfx_mean)
          write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_trbx_t",this%step,".3Dstt"
          fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
          call decomp_2d_write_one(3, rbuff3, fname)

          ! write out turbfy_mean3D
          call transpose_x_to_y(this%turbfy_mean3D/tidSUMreal, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2,                        rbuff3, this%gpC)
          call this%compute_z_mean(rbuff3, this%turbfy_mean)
          write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_trby_t",this%step,".3Dstt"
          fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
          call decomp_2d_write_one(3, rbuff3, fname)

          ! write out turbfz_mean3D
          call transpose_x_to_y(this%turbfz_mean3D/tidSUMreal, rbuff2, this%gpC)
          call transpose_y_to_z(rbuff2,                        rbuff3, this%gpC)
          call this%compute_z_mean(rbuff3, this%turbfz_mean)
          write(tempname,"(A3,I2.2,A7,I6.6,A6)") "Run",this%runID, "_trbz_t",this%step,".3Dstt"
          fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
          call decomp_2d_write_one(3, rbuff3, fname)
       endif

       if(this%isStratified) then
           ! T_avg
           call transpose_x_to_y(this%T_mean3D/tidSUMreal, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2,                   rbuff6, this%gpC)
           call this%compute_z_mean(rbuff6, this%T_mean)
           write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_Tm_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff6, fname)
     
           ! uT_avg - u_avg*T_avg - Reynolds only
           rbuff1 = this%uT_mean3D/tidSUMreal - this%u_mean3D * this%T_mean3D / tidSUMreal**2
           call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%uT_mean)       !--Reynolds
           write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_uTm_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
     
           ! vT_avg - v_avg*T_avg - Reynolds only
           rbuff1 = this%vT_mean3D/tidSUMreal - this%v_mean3D * this%T_mean3D / tidSUMreal**2
           call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%vT_mean)
           write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_vTm_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
     
           ! wT_avg - w_avg*T_avg - Reynolds only
           rbuff1 = this%wT_mean3D/tidSUMreal - this%w_mean3D * this%T_mean3D / tidSUMreal**2
           call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%wT_mean)
           write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_wTm_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)
     
           ! TT_avg - T_avg*T_avg - Reynolds only
           rbuff1 = this%TT_mean3D/tidSUMreal - this%T_mean3D * this%T_mean3D / tidSUMreal**2
           call transpose_x_to_y(rbuff1, rbuff2, this%gpC)
           call transpose_y_to_z(rbuff2, rbuff3, this%gpC)
           call this%compute_z_mean(rbuff3, this%TT_mean)
           write(tempname,"(A3,I2.2,A6,I6.6,A6)") "Run",this%runID, "_TTm_t",this%step,".3Dstt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_one(3, rbuff3, fname)

           if(this%isStratified) then
               ! q1_mean
               call transpose_x_to_y(this%q1_mean3D/tidSUMreal, rbuff2, this%gpC)
               call transpose_y_to_z(rbuff2,                    rbuff3, this%gpC)
               call this%compute_z_mean(rbuff3, this%q1_mean)
               write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_q1_t",this%step,".3Dstt"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_one(3, rbuff3, fname)
    
               ! q2_mean
               call transpose_x_to_y(this%q2_mean3D/tidSUMreal, rbuff2, this%gpC)
               call transpose_y_to_z(rbuff2,                    rbuff3, this%gpC)
               call this%compute_z_mean(rbuff3, this%q2_mean)
               write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_q2_t",this%step,".3Dstt"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_one(3, rbuff3, fname)
    
               ! q3_mean
               call transpose_x_to_y(this%q3_mean3D/tidSUMreal, rbuff2, this%gpC)
               call transpose_y_to_z(rbuff2,                    rbuff3, this%gpC)
               call this%compute_z_mean(rbuff3, this%q3_mean)
               write(tempname,"(A3,I2.2,A5,I6.6,A6)") "Run",this%runID, "_q3_t",this%step,".3Dstt"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_one(3, rbuff3, fname)
           endif
       endif

       call message(1, "Dumped 3D stats files")

       nullify(rbuff1,rbuff2,rbuff3,rbuff4,rbuff5,rbuff6, rbuff0, cbuffy1, cbuffy2, rbuff3E, rbuff4E)

       ! dump horizontal averages
       if(this%useWindTurbines) then
           this%runningSum_turb = zero
           call MPI_reduce(this%runningSum_sc_turb, this%runningSum_turb, 8*this%WindTurbineArr%nTurbines, mpirkind, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       endif
       if (nrank == 0) then
           write(tempname,"(A3,I2.2,A2,I6.6,A4)") "Run", this%RunID,"_t",tid,".stt"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call write_2d_ascii(this%horzavgstats,fname)

           write(tempname,"(A3,I2.2,A2,I6.6,A4)") "Run", this%RunID,"_t",tid,".sth"   ! time and horz averages of scalars
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           open(unit=771,file=fname,status='unknown')
           if(this%useWindTurbines) then
               write(771,'(e19.12,1x,i7,1x,8008(e19.12,1x))') this%tsim, this%tidSUM, this%runningSum_sc/tidSUMreal, this%runningSum_turb/tidSUMreal ! change if using more than 1000 turbines
           else
               write(771,'(e19.12,1x,i7,1x,5(e19.12,1x))') this%tsim, this%tidSUM, this%runningSum_sc/tidSUMreal
           endif
           close(771)
       end if
       !write(200+nrank,*) 4, this%tsim, this%runningSum_sc(2)/tidSUMreal
       call message(1, "Just dumped a .stt file")
       call message(2, "Number ot tsteps averaged:",this%tidSUM)

       if (this%computeSpectra) then
           ! Dump horizontally averaged x-spectra
           dirid = 2; decompdir = 2

           ! --- only 4, 5 or 6 planes of xspextra_mean in y-direction are being used
           nspectra = 4
           if(this%isStratified)                             nspectra = nspectra + 1
           if(this%fastCalcPressure .or. this%storePressure) nspectra = nspectra+1

           ! --- for k1 = 1, multiplication factor is 1.0,
           ! --- for k1 = 2:Nx/2+1, multiplication factor is 2.0
           normfac = two/real(size(this%cbuffyC(:,:,:,1),2),rkind)/tidSUMreal
           tmpvar(1:this%sp_gpC%ysz(1),1:nspectra,:) = normfac*this%xspectra_mean(1:this%sp_gpC%ysz(1),1:nspectra,:)
           if(this%sp_gpC%yst(1)==1) then
               tmpvar(1,1:nspectra,:) = half*tmpvar(1,1:nspectra,:)
           endif

           jindx = 1 ! u
           write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%RunID,"_specu_t",tid,".out"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%sp_gpC)

           jindx = 2 ! v
           write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%RunID,"_specv_t",tid,".out"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%sp_gpC)

           jindx = 3 ! w
           write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%RunID,"_specw_t",tid,".out"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%sp_gpC)

           jindx = 4 ! KE
           write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%RunID,"_speck_t",tid,".out"
           fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
           call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%sp_gpC)

           if(this%isStratified) then
               jindx = jindx + 1 ! T
               write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%RunID,"_specT_t",tid,".out"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%sp_gpC)
           endif

           if(this%fastCalcPressure .or. this%storePressure) then
               jindx = jindx + 1 ! p
               write(tempname,"(A3,I2.2,A8,I6.6,A4)") "Run", this%RunID,"_specp_t",tid,".out"
               fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
               call decomp_2d_write_plane(decompdir, tmpvar, dirid, jindx, fname, this%sp_gpC)
           endif
       end if 

       nullify(S11_mean3D, S12_mean3D, S13_mean3D, S22_mean3D, S23_mean3D, S33_mean3D)

   end subroutine

   subroutine finalize_stats3D(this)
       class(igrid), intent(inout) :: this

       nullify(this%u_mean, this%v_mean, this%w_mean, this%uu_mean, this%uv_mean, this%uw_mean, this%vv_mean, this%vw_mean, this%ww_mean, this%disperuw_mean, this%dispervw_mean)
       nullify(this%mkeadv_mean, this%mkett_mean, this%mkedisp_mean, this%tkeadv_mean, this%tkett_mean, this%tkeprod_mean)
       nullify(this%u_mean3D, this%v_mean3D, this%w_mean3D, this%uu_mean3D, this%uv_mean3D, this%uw_mean3D, this%vv_mean3D, this%vw_mean3D, this%ww_mean3D, this%tketurbtranspx_mean3D, this%tketurbtranspy_mean3D, this%tketurbtranspz_mean3D)

       if(this%fastCalcPressure .or. this%storePressure) then
           nullify(this%p_mean,   this%mkept_mean,  this%tkept_mean)
           nullify(this%p_mean3D, this%pu_mean3D, this%pv_mean3D, this%pw_mean3D)
       endif

       if(.not. this%isInviscid) then
           nullify(this%mkevdif_mean, this%mkevdsp_mean, this%tkevdif_mean, this%tkevdsp_mean)
           nullify(this%viscdisp_mean3D, this%Siju1_mean3D, this%Siju2_mean3D, this%Siju3_mean3D)
       endif

       if(this%useSGS) then
          nullify(this%mkesgst_mean, this%mkesgsd_mean, this%tkesgst_mean, this%tkesgsd_mean)
          nullify(this%tau11_mean, this%tau12_mean, this%tau13_mean, this%tau22_mean, this%tau23_mean, this%tau33_mean)
          nullify(this%tau11_mean3D, this%tau12_mean3D, this%tau13_mean3D, this%tau22_mean3D, this%tau23_mean3D, this%tau33_mean3D)
          nullify(this%sgsdissp_mean3D, this%tauu1_mean3D, this%tauu2_mean3D, this%tauu3_mean3D)
       endif

       if(this%useWindTurbines) then
           nullify(this%turbfx_mean, this%turbfy_mean, this%turbfz_mean, this%tketurbf_mean, this%mketurbf_mean)
           nullify(this%turbfx_mean3D, this%turbfy_mean3D, this%turbfz_mean3D, this%uturbf_mean3D)
       endif

       if(this%isStratified) then
           nullify(this%T_mean, this%uT_mean, this%vT_mean, this%wT_mean, this%TT_mean)
           nullify(this%T_mean3D, this%uT_mean3D, this%vT_mean3D, this%wT_mean3D, this%TT_mean3D)
          if(this%useSGS) then
             nullify(this%q1_mean, this%q2_mean, this%q3_mean)
             nullify(this%q1_mean3D, this%q2_mean3D, this%q3_mean3D)
          endif
       endif

       !deallocate(this%stats3D, this%horzavgstats, this%inst_horz_avg, this%runningSum_sc, this%xspectra_mean)
       if(this%useWindTurbines) deallocate(this%inst_horz_avg_turb, this%runningSum_sc_turb, this%runningSum_turb)
   end subroutine


   !--------------------------------Done 3D Statistics----------------------------------------------

   !!--------------------------------Beginning 1D Statistics----------------------------------------------
   !subroutine init_stats( this)
   !    class(igrid), intent(inout), target :: this
   !    type(decomp_info), pointer  :: gpC

   !    gpC => this%gpC
   !    this%tidSUM = 0

   !    if (this%isStratified) then
   !        allocate(this%zStats2dump(this%nz,33))
   !        allocate(this%runningSum(this%nz,33))
   !        allocate(this%TemporalMnNOW(this%nz,33))
   !        allocate(this%runningSum_sc(5))
   !        allocate(this%inst_horz_avg(5))
   !    else
   !        allocate(this%zStats2dump(this%nz,25))
   !        allocate(this%runningSum(this%nz,25))
   !        allocate(this%TemporalMnNOW(this%nz,25))
   !        allocate(this%runningSum_sc(3))
   !        allocate(this%inst_horz_avg(3))
   !    end if 

   !    if(this%useWindTurbines) then
   !        allocate(this%inst_horz_avg_turb(8*this%WindTurbineArr%nTurbines))
   !        allocate(this%runningSum_sc_turb(8*this%WindTurbineArr%nTurbines))
   !        allocate(this%runningSum_turb   (8*this%WindTurbineArr%nTurbines))
   !    endif

   !    ! mean velocities
   !    this%u_mean => this%zStats2dump(:,1);  this%v_mean  => this%zStats2dump(:,2);  this%w_mean => this%zStats2dump(:,3) 

   !    ! mean squared velocities
   !    this%uu_mean => this%zStats2dump(:,4); this%uv_mean => this%zStats2dump(:,5); this%uw_mean => this%zStats2dump(:,6)
   !                                           this%vv_mean => this%zStats2dump(:,7); this%vw_mean => this%zStats2dump(:,8) 
   !                                                                                  this%ww_mean => this%zStats2dump(:,9)

   !    ! SGS stresses
   !    this%tau11_mean => this%zStats2dump(:,10); this%tau12_mean => this%zStats2dump(:,11); this%tau13_mean => this%zStats2dump(:,12)
   !                                               this%tau22_mean => this%zStats2dump(:,13); this%tau23_mean => this%zStats2dump(:,14) 
   !                                                                                          this%tau33_mean => this%zStats2dump(:,15)

   !    ! SGS dissipation
   !    this%sgsdissp_mean => this%zStats2dump(:,16)

   !    ! velocity derivative products - for viscous dissipation
   !    this%viscdisp_mean => this%zStats2dump(:,17)

   !    ! means of velocity derivatives
   !    this%S11_mean => this%zStats2dump(:,18); this%S12_mean => this%zStats2dump(:,19); this%S13_mean => this%zStats2dump(:,20)
   !                                             this%S22_mean => this%zStats2dump(:,21); this%S23_mean => this%zStats2dump(:,22)
   !                                                                                      this%S33_mean => this%zStats2dump(:,23)

   !    ! SGS model coefficient
   !    this%sgscoeff_mean => this%zStats2dump(:,24)

   !    this%PhiM => this%zStats2dump(:,25)
   !    
   !    if (this%isStratified) then
   !        this%TT_mean => this%zStats2dump(:,30);  this%wT_mean => this%zStats2Dump(:,29);  this%vT_mean => this%zStats2Dump(:,28)
   !        this%uT_mean => this%zStats2dump(:,27);  this%T_mean => this%zStats2Dump(:,26); this%q1_mean => this%zStats2Dump(:,31)
   !        this%q2_mean => this%zStats2dump(:,32);  this%q3_mean => this%zStats2Dump(:,33)
   !    end if 
   !    this%runningSum_sc = zero
   !    this%runningSum = zero
   !    this%TemporalMnNOW = zero
   !    this%zStats2dump = zero
   !    this%inst_horz_avg = zero
   !    if(this%useWindTurbines) then
   !        this%inst_horz_avg_turb = zero
   !        this%runningSum_sc_turb = zero
   !        this%runningSum_turb    = zero
   !    endif
   !    nullify(gpC)
   !end subroutine

   !subroutine compute_stats(this)
   !    class(igrid), intent(inout), target :: this
   !    type(decomp_info), pointer :: gpC
   !    real(rkind), dimension(:,:,:), pointer :: rbuff1, rbuff2, rbuff3E, rbuff2E, rbuff3, rbuff4, rbuff5, rbuff5E, rbuff4E, rbuff6E, rbuff6!, rbuff7

   !    rbuff1  => this%rbuffxC(:,:,:,1); rbuff2  => this%rbuffyC(:,:,:,1);
   !    rbuff2E => this%rbuffyE(:,:,:,1); rbuff3E => this%rbuffzE(:,:,:,1);
   !    rbuff3 => this%rbuffzC(:,:,:,1); rbuff4E => this%rbuffzE(:,:,:,2);
   !    rbuff4 => this%rbuffzC(:,:,:,2); rbuff5E => this%rbuffzE(:,:,:,3)
   !    rbuff5 => this%rbuffzC(:,:,:,3); rbuff6E => this%rbuffzE(:,:,:,4)
   !    rbuff6 => this%rbuffzC(:,:,:,4); !rbuff7 => this%rbuffzC(:,:,:,5); 
   !    gpC => this%gpC

   !    this%tidSUM = this%tidSUM + 1


   !    ! Compute u - mean 
   !    call transpose_x_to_y(this%u,rbuff2,this%gpC)
   !    call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !    call this%compute_z_mean(rbuff3, this%u_mean)
   !    !if (this%normByustar) this%u_mean = this%u_mean/this%sgsmodel%get_ustar()

   !    ! Compute v - mean 
   !    call transpose_x_to_y(this%v,rbuff2,this%gpC)
   !    call transpose_y_to_z(rbuff2,rbuff4,this%gpC)
   !    call this%compute_z_mean(rbuff4, this%v_mean)
   !    !if (this%normByustar)this%v_mean = this%v_mean/this%sgsmodel%get_ustar()

   !    ! Compute wC - mean 
   !    call transpose_x_to_y(this%wC,rbuff2,this%gpC)
   !    call transpose_y_to_z(rbuff2,rbuff5,this%gpC)
   !    call this%compute_z_mean(rbuff5, this%w_mean)
   !    !if (this%normByustar)this%w_mean = this%w_mean/this%sgsmodel%get_ustar()

   !    ! take w from x -> z decomp
   !    call transpose_x_to_y(this%w,rbuff2E,this%gpE)
   !    call transpose_y_to_z(rbuff2E,rbuff5E,this%gpE)

   !    ! take uE from x -> z decomp
   !    call transpose_x_to_y(this%uE,rbuff2E,this%gpE)
   !    call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)

   !    ! take vE from x -> z decomp
   !    call transpose_x_to_y(this%vE,rbuff2E,this%gpE)
   !    call transpose_y_to_z(rbuff2E,rbuff4E,this%gpE)

   !    ! uu mean
   !    rbuff6 = rbuff3*rbuff3
   !    call this%compute_z_mean(rbuff6, this%uu_mean)
   !    !if (this%normByustar)this%uu_mean = this%uu_mean/(this%sgsmodel%get_ustar()**2)

   !    ! uv mean
   !    rbuff6 = rbuff3*rbuff4
   !    call this%compute_z_mean(rbuff6, this%uv_mean)
   !    !if (this%normByustar)this%uv_mean = this%uv_mean/(this%sgsmodel%get_ustar()**2)

   !    ! uw mean
   !    rbuff6E = rbuff3E*rbuff5E
   !    call this%OpsPP%InterpZ_Edge2Cell(rbuff6E,rbuff6)
   !    call this%compute_z_mean(rbuff6, this%uw_mean)
   !    !if (this%normByustar)this%uw_mean = this%uw_mean/(this%sgsmodel%get_ustar()**2)

   !    ! vv mean 
   !    rbuff6 = rbuff4*rbuff4
   !    call this%compute_z_mean(rbuff6, this%vv_mean)
   !    !if (this%normByustar)this%vv_mean = this%vv_mean/(this%sgsmodel%get_ustar()**2)

   !    ! vw mean 
   !    rbuff6E = rbuff4E*rbuff5E
   !    call this%OpsPP%InterpZ_Edge2Cell(rbuff6E,rbuff6)
   !    call this%compute_z_mean(rbuff6, this%vw_mean)
   !    !if (this%normByustar)this%vw_mean = this%vw_mean/(this%sgsmodel%get_ustar()**2)

   !    ! ww mean 
   !    rbuff6 = rbuff5*rbuff5
   !    call this%compute_z_mean(rbuff6, this%ww_mean)
   !    !if (this%normByustar)this%ww_mean = this%ww_mean/(this%sgsmodel%get_ustar()**2)

   !    ! Statified Stuff
   !    if (this%isStratified) then
   !        ! T mean
   !        call transpose_x_to_y(this%T,rbuff2,this%gpC)
   !        call transpose_y_to_z(rbuff2,rbuff6,this%gpC)
   !        call this%compute_z_mean(rbuff6, this%T_mean)

   !        ! uT mean
   !        rbuff3 = rbuff3*rbuff6
   !        call this%compute_z_mean(rbuff3, this%uT_mean)
   !        
   !        ! vT mean
   !        rbuff4 = rbuff4*rbuff6
   !        call this%compute_z_mean(rbuff4, this%vT_mean)
   !        
   !        ! wT mean
   !        rbuff5 = rbuff5*rbuff6
   !        call this%compute_z_mean(rbuff5, this%wT_mean)
   !        
   !        ! TT mean
   !        rbuff6 = rbuff6*rbuff6
   !        call this%compute_z_mean(rbuff6, this%TT_mean)
   !    end if 


   !    if (this%useSGS) then
   !        ! tau_11
   !        call transpose_x_to_y(this%tauSGS_ij(:,:,:,1),rbuff2,this%gpC)
   !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !        call this%compute_z_mean(rbuff3, this%tau11_mean)
   !        if (this%normByustar)this%tau11_mean = this%tau11_mean/(this%sgsmodel%get_ustar()**2)

   !        ! tau_12
   !        call transpose_x_to_y(this%tauSGS_ij(:,:,:,2),rbuff2,this%gpC)
   !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !        call this%compute_z_mean(rbuff3, this%tau12_mean)
   !        if (this%normByustar)this%tau12_mean = this%tau12_mean/(this%sgsmodel%get_ustar()**2)

   !        ! tau_13
   !        call transpose_x_to_y(this%tau13,rbuff2E,this%gpE)
   !        call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
   !        rbuff3E(:,:,1) = -(this%sgsmodel%get_ustar()**2)
   !        call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
   !        call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
   !        call transpose_y_to_x(rbuff2,this%tauSGS_ij(:,:,:,3),this%gpC)
   !        call this%compute_z_mean(rbuff3, this%tau13_mean)
   !        if (this%normByustar)this%tau13_mean = this%tau13_mean/(this%sgsmodel%get_ustar()**2)

   !        ! tau_22
   !        call transpose_x_to_y(this%tauSGS_ij(:,:,:,4),rbuff2,this%gpC)
   !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !        call this%compute_z_mean(rbuff3, this%tau22_mean)
   !        if (this%normByustar)this%tau22_mean = this%tau22_mean/(this%sgsmodel%get_ustar()**2)

   !        ! tau_23
   !        call transpose_x_to_y(this%tau23,rbuff2E,this%gpE)
   !        call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
   !        rbuff3E(:,:,1) = -(this%sgsmodel%get_ustar()**2)*this%sgsmodel%get_vmean()/this%sgsmodel%get_umean()
   !        call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
   !        call transpose_z_to_y(rbuff3,rbuff2,this%gpC)
   !        call transpose_y_to_x(rbuff2,this%tauSGS_ij(:,:,:,5),this%gpC)
   !        call this%compute_z_mean(rbuff3, this%tau23_mean)
   !        if (this%normByustar)this%tau23_mean = this%tau23_mean/(this%sgsmodel%get_ustar()**2)

   !        ! tau_33
   !        call transpose_x_to_y(this%tauSGS_ij(:,:,:,6),rbuff2,this%gpC)
   !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !        call this%compute_z_mean(rbuff3, this%tau33_mean)
   !        if (this%normByustar)this%tau33_mean = this%tau33_mean/(this%sgsmodel%get_ustar()**2)


   !        ! sgs dissipation
   !        rbuff1 = this%tauSGS_ij(:,:,:,1)*this%tauSGS_ij(:,:,:,1) + &
   !                 this%tauSGS_ij(:,:,:,2)*this%tauSGS_ij(:,:,:,2) + &
   !                 this%tauSGS_ij(:,:,:,3)*this%tauSGS_ij(:,:,:,3)
   !        rbuff1 = rbuff1 + two*(this%tauSGS_ij(:,:,:,4)*this%tauSGS_ij(:,:,:,4) + &
   !                               this%tauSGS_ij(:,:,:,5)*this%tauSGS_ij(:,:,:,5) + &
   !                               this%tauSGS_ij(:,:,:,6)*this%tauSGS_ij(:,:,:,6) )
   !        rbuff1 = rbuff1/(this%nu_SGS + 1.0d-14)         ! note: factor of half is in dump_stats

   !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
   !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !        call this%compute_z_mean(rbuff3, this%sgsdissp_mean)

   !        ! viscous dissipation- *****????? Is rbuff1 contaminated after transpose_x_to_y? *****?????
   !        rbuff1 = rbuff1/(this%nu_SGS + 1.0d-14)        ! note: factor of fourth is in dump_stats
   !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
   !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !        call this%compute_z_mean(rbuff3, this%viscdisp_mean)

   !        ! note: factor of half in all S_** is in dump_stats
   !        ! S_11
   !        rbuff1 = this%tauSGS_ij(:,:,:,1)/(this%nu_SGS + 1.0d-14)
   !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
   !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !        call this%compute_z_mean(rbuff3, this%S11_mean)

   !        ! S_12
   !        rbuff1 = this%tauSGS_ij(:,:,:,2)/(this%nu_SGS + 1.0d-14)
   !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
   !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !        call this%compute_z_mean(rbuff3, this%S12_mean)

   !        ! S_13
   !        rbuff1 = this%tauSGS_ij(:,:,:,3)/(this%nu_SGS + 1.0d-14)
   !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
   !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !        call this%compute_z_mean(rbuff3, this%S13_mean)

   !        ! S_22
   !        rbuff1 = this%tauSGS_ij(:,:,:,4)/(this%nu_SGS + 1.0d-14)
   !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
   !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !        call this%compute_z_mean(rbuff3, this%S22_mean)

   !        ! S_23
   !        rbuff1 = this%tauSGS_ij(:,:,:,5)/(this%nu_SGS + 1.0d-14)
   !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
   !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !        call this%compute_z_mean(rbuff3, this%S23_mean)

   !        ! S_33
   !        rbuff1 = this%tauSGS_ij(:,:,:,6)/(this%nu_SGS + 1.0d-14)
   !        call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
   !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !        call this%compute_z_mean(rbuff3, this%S33_mean)

   !        ! sgs coefficient
   !        call transpose_x_to_y(this%c_SGS,rbuff2,this%gpC)
   !        call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !        !call this%compute_z_mean(rbuff3, this%sgscoeff_mean)    ! -- averaging not needed
   !        this%sgscoeff_mean(:) = rbuff3(1,1,:)
   !   
   !        
   !        if (this%isStratified) then
   !            ! q1
   !            call transpose_x_to_y(this%q1,rbuff2,this%gpC)
   !            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !            call this%compute_z_mean(rbuff3, this%q1_mean)    

   !            ! q2
   !            call transpose_x_to_y(this%q2,rbuff2,this%gpC)
   !            call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !            call this%compute_z_mean(rbuff3, this%q2_mean)    

   !            ! q3
   !            call transpose_x_to_y(this%q3,rbuff2E,this%gpE)
   !            call transpose_y_to_z(rbuff2E,rbuff3E,this%gpE)
   !            rbuff3E(:,:,1) = this%wTh_surf
   !            call this%OpsPP%InterpZ_Edge2Cell(rbuff3E,rbuff3)
   !            call this%compute_z_mean(rbuff3, this%q3_mean)    
   !        end if 
   !    end if

   !    rbuff1 = this%duidxjC(:,:,:,3)*this%mesh(:,:,:,3)
   !    call transpose_x_to_y(rbuff1,rbuff2,this%gpC)
   !    call transpose_y_to_z(rbuff2,rbuff3,this%gpC)
   !    call this%compute_z_mean(rbuff3, this%PhiM)
   !    if (this%useSGS) then
   !     this%PhiM = this%PhiM*kappa/this%sgsmodel%get_ustar()
   !    end if

   !    this%runningSum = this%runningSum + this%zStats2dump

   !    !write(*,*) 'In stats'
   !    !write(*,*) 'umean', maxval(this%u_mean), minval(this%u_mean)
   !    !write(*,*) 'vmean', maxval(this%v_mean), minval(this%v_mean)
   !    !write(*,*) 'wmean', maxval(this%w_mean), minval(this%w_mean)

   !    ! instantaneous horizontal averages of some quantities
   !    if (this%useSGS) this%inst_horz_avg(1) = this%sgsmodel%get_ustar()
   !    ! this%inst_horz(2) and (3) are computed in getRHS_SGS_WallM
   !    if(this%isStratified) then
   !        this%inst_horz_avg(4) = this%invObLength
   !        this%inst_horz_avg(5) = this%wTh_surf
   !    endif
   !    ! this%inst_horz_avg_turb(1:5*this%WindTurbineArr%nTurbines) is computed in this%WindTurbineArr%getForceRHS
   !    this%runningSum_sc = this%runningSum_sc + this%inst_horz_avg
   !    if(this%useWindTurbines) this%runningSum_sc_turb = this%runningSum_sc_turb + this%inst_horz_avg_turb

   !end subroutine 

   !subroutine dump_stats(this)
   !    use basic_io, only: write_2d_ascii, write_2D_binary
   !    use exits, only: message
   !    use kind_parameters, only: clen, mpirkind
   !    use mpi
   !    class(igrid), intent(inout), target :: this
   !    character(len=clen) :: fname
   !    character(len=clen) :: tempname
   !    integer :: tid, ierr

   !    this%TemporalMnNOW = this%runningSum/real(this%tidSUM,rkind)
   !    tid = this%step

   !    ! compute (u_i'u_j')
   !    this%TemporalMnNOW(:,4) = this%TemporalMnNOW(:,4) - this%TemporalMnNOW(:,1)*this%TemporalMnNOW(:,1)
   !    this%TemporalMnNOW(:,5) = this%TemporalMnNOW(:,5) - this%TemporalMnNOW(:,1)*this%TemporalMnNOW(:,2)
   !    this%TemporalMnNOW(:,6) = this%TemporalMnNOW(:,6) - this%TemporalMnNOW(:,1)*this%TemporalMnNOW(:,3)
   !    this%TemporalMnNOW(:,7) = this%TemporalMnNOW(:,7) - this%TemporalMnNOW(:,2)*this%TemporalMnNOW(:,2)
   !    this%TemporalMnNOW(:,8) = this%TemporalMnNOW(:,8) - this%TemporalMnNOW(:,2)*this%TemporalMnNOW(:,3)
   !    this%TemporalMnNOW(:,9) = this%TemporalMnNOW(:,9) - this%TemporalMnNOW(:,3)*this%TemporalMnNOW(:,3)

   !    if (this%isStratified) then
   !        this%TemporalMnNOW(:,27) = this%TemporalMnNOW(:,27) - this%TemporalMnNOW(:,1)*this%TemporalMnNOW(:,26)
   !        this%TemporalMnNOW(:,28) = this%TemporalMnNOW(:,28) - this%TemporalMnNOW(:,2)*this%TemporalMnNOW(:,26)
   !        this%TemporalMnNOW(:,29) = this%TemporalMnNOW(:,29) - this%TemporalMnNOW(:,3)*this%TemporalMnNOW(:,26)
   !        this%TemporalMnNOW(:,30) = this%TemporalMnNOW(:,30) - this%TemporalMnNOW(:,26)*this%TemporalMnNOW(:,26)
   !    end if 

   !    ! compute sgs dissipation
   !    this%TemporalMnNOW(:,16) = half*this%TemporalMnNOW(:,16)

   !    ! compute viscous dissipation
   !    this%TemporalMnNOW(:,17) = this%TemporalMnNOW(:,17) - (                        &
   !                               this%TemporalMnNOW(:,18)*this%TemporalMnNOW(:,18) + &
   !                               this%TemporalMnNOW(:,21)*this%TemporalMnNOW(:,21) + &
   !                               this%TemporalMnNOW(:,23)*this%TemporalMnNOW(:,23) + &
   !                          two*(this%TemporalMnNOW(:,19)*this%TemporalMnNOW(:,19) + &
   !                               this%TemporalMnNOW(:,20)*this%TemporalMnNOW(:,20) + & 
   !                               this%TemporalMnNOW(:,22)*this%TemporalMnNOW(:,22)))
   !    this%TemporalMnNOW(:,17) = half*this%TemporalMnNOW(:,17)/this%Re     ! note: this is actually 2/Re*(..)/4

   !    if(this%useWindTurbines) then
   !        this%runningSum_turb = zero
   !        call MPI_reduce(this%runningSum_sc_turb, this%runningSum_turb, 8*this%WindTurbineArr%nTurbines, mpirkind, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
   !    endif
   !    if (nrank == 0) then
   !        write(tempname,"(A3,I2.2,A2,I6.6,A4)") "Run", this%RunID,"_t",tid,".stt"
   !        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
   !        call write_2d_ascii(this%TemporalMnNOW,fname)
   !        !call write_2D_binary(TemporalMnNOW,fname)

   !        write(tempname,"(A3,I2.2,A2,I6.6,A4)") "Run", this%RunID,"_t",tid,".sth"   ! time and horz averages of scalars
   !        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
   !        open(unit=771,file=fname,status='unknown')
   !        if(this%useWindTurbines) then
   !            write(771,'(e19.12,1x,i7,1x,8008(e19.12,1x))') this%tsim, this%tidSUM, this%runningSum_sc/real(this%tidSUM,rkind), this%runningSum_turb/real(this%tidSUM,rkind) ! change if using more than 1000 turbines
   !        else
   !            write(771,'(e19.12,1x,i7,1x,5(e19.12,1x))') this%tsim, this%tidSUM, this%runningSum_sc/real(this%tidSUM,rkind)
   !        endif
   !        close(771)
   !    end if
   !    call message(1, "Just dumped a .stt file")
   !    call message(2, "Number ot tsteps averaged:",this%tidSUM)

   !end subroutine

   !subroutine compute_z_fluct(this,fin)
   !    use reductions, only: P_SUM
   !    class(igrid), intent(in), target :: this
   !    real(rkind), dimension(:,:,:), intent(inout) :: fin
   !    integer :: k
   !    real(rkind) :: fmean

   !    do k = 1,size(fin,3)
   !        fmean = P_SUM(sum(fin(:,:,k)))/(real(this%nx,rkind)*real(this%ny,rkind))
   !        fin(:,:,k) = fin(:,:,k) - fmean
   !    end do 

   !end subroutine

   !subroutine compute_y_mean(this, arr_in, arr_out)
   !    use reductions, only: P_SUM
   !    class(igrid), intent(in), target :: this
   !    real(rkind), dimension(:,:,:), intent(in) :: arr_in
   !    real(rkind), dimension(:,:), intent(out) :: arr_out
   !    integer :: k, i

   !    do k = 1,size(arr_in,3)
   !      do i = 1,size(arr_in,1)
   !        arr_out(i,k) = P_SUM(sum(arr_in(:,:,k)))/(real(this%nx,rkind)*real(this%ny,rkind))
   !    end do 

   !end subroutine

   subroutine compute_z_mean(this, arr_in, vec_out)
       use reductions, only: P_SUM
       class(igrid), intent(in), target :: this
       real(rkind), dimension(:,:,:), intent(in) :: arr_in
       real(rkind), dimension(:), intent(out) :: vec_out
       integer :: k

       do k = 1,size(arr_in,3)
           vec_out(k) = P_SUM(sum(arr_in(:,:,k)))/(real(this%nx,rkind)*real(this%ny,rkind))
       end do 

   end subroutine

   !subroutine finalize_stats(this)
   !    class(igrid), intent(inout) :: this

   !    nullify(this%u_mean, this%v_mean, this%w_mean, this%uu_mean, this%uv_mean, this%uw_mean, this%vv_mean, this%vw_mean, this%ww_mean)
   !    nullify(this%tau11_mean, this%tau12_mean, this%tau13_mean, this%tau22_mean, this%tau23_mean, this%tau33_mean)
   !    nullify(this%S11_mean, this%S12_mean, this%S13_mean, this%S22_mean, this%S23_mean, this%S33_mean)
   !    nullify(this%sgsdissp_mean, this%viscdisp_mean, this%sgscoeff_mean)
   !    if (allocated(this%zStats2dump)) deallocate(this%zStats2dump)
   !    if (allocated(this%runningSum)) deallocate(this%runningSum)
   !    if (allocated(this%TemporalMnNow)) deallocate(this%TemporalMnNOW)
   !    if (allocated(this%runningSum_sc)) deallocate(this%runningSum_sc)
   !    if(this%useWindTurbines) deallocate(this%inst_horz_avg_turb, this%runningSum_sc_turb, this%runningSum_turb)
   !end subroutine 

   !!--------------------------------Done 1D Statistics----------------------------------------------
