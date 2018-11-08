    pure subroutine get_geostrophic_forcing(this, Fg_x, Fg_y)
        class(igrid), intent(in) :: this
        real(rkind), intent(out) :: Fg_x, Fg_y
        real(rkind) :: gx, gy 

        gx = this%G_GEOSTROPHIC*cos(this%G_ALPHA*pi/180.d0)
        gy = this%G_GEOSTROPHIC*sin(this%G_ALPHA*pi/180.d0)

        Fg_x = -this%coriolis_omegaZ*(two/this%Ro)*gy
        Fg_y =  this%coriolis_omegaZ*(two/this%Ro)*gx


    end subroutine 


    subroutine addCoriolisTerm(this, urhs, vrhs, wrhs)
       class(igrid), intent(inout), target :: this
       complex(rkind), dimension(:,:,:), pointer :: ybuffE, ybuffC1, ybuffC2, zbuffC, zbuffE
       complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout) :: urhs, vrhs
       complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: wrhs
       !real(rkind) :: frameAngle!, latitude

       ybuffE => this%cbuffyE(:,:,:,1)
       ybuffC1 => this%cbuffyC(:,:,:,1)
       ybuffC2 => this%cbuffyC(:,:,:,2)
       zbuffE => this%cbuffzE(:,:,:,1)
       zbuffC => this%cbuffzC(:,:,:,1)


       if (this%newTimestep) then
           if (this%assume_fplane) then
               this%coriolis_omegaZ   = sin(this%latitude*pi/180.d0)
               this%coriolis_omegaY = 0.d0
               this%coriolis_omegaX = 0.d0
           else
               this%coriolis_omegaX = cos(this%latitude*pi/180.d0)*sin(this%frameAngle*pi/180.d0)
               this%coriolis_omegaZ = sin(this%latitude*pi/180.d0)
               this%coriolis_omegaY = cos(this%latitude*pi/180.d0)*cos(this%frameAngle*pi/180.d0)
           end if
           this%rbuffxC(:,:,:,1) = this%G_GEOSTROPHIC*cos(this%G_ALPHA*pi/180.d0)
           call this%spectC%fft(this%rbuffxC(:,:,:,1),this%Gxhat)
           this%rbuffxE(:,:,:,1) = this%G_GEOSTROPHIC*cos(this%G_ALPHA*pi/180.d0)
           call this%spectE%fft(this%rbuffxE(:,:,:,1),this%Gxhat_Edge)
           this%rbuffxC(:,:,:,1) = this%G_GEOSTROPHIC*sin(this%G_ALPHA*pi/180.d0)
           call this%spectC%fft(this%rbuffxC(:,:,:,1),this%Gyhat)
           this%rbuffxE(:,:,:,1) = this%G_GEOSTROPHIC*sin(this%G_ALPHA*pi/180.d0)
           call this%spectE%fft(this%rbuffxE(:,:,:,1),this%Gyhat_Edge)
       end if
       ! u equation 
       ybuffC1    = (two/this%Ro)*(-this%coriolis_omegaY*this%whatC - this%coriolis_omegaZ*(this%GyHat - this%vhat))
       !this%u_rhs = this%u_rhs  + ybuffC1
       urhs = urhs + ybuffC1

       ! v equation
       !ybuffC2    = (two/this%Ro)*(this%coriolis_omegaZ*(this%GxHat - this%uhat))
       ybuffC2    = (two/this%Ro)*(this%coriolis_omegaZ*(this%GxHat - this%uhat) + this%coriolis_omegaX*this%whatC)
       !this%v_rhs = this%v_rhs + ybuffC2 
       vrhs = vrhs + ybuffC2

       ! w equation 
       ! The real equation is given as:
       ! this%w_rhs = this%w_rhs - this%coriolis_omegaY*(this%GxHat - this%uhat)/this%Ro
       ! But we evaluate this term as:
     
       ybuffE = (two/this%Ro)*(-this%coriolis_omegaY*(this%Gxhat_Edge - this%uEhat) + this%coriolis_omegaX*(this%Gyhat_Edge - this%vEhat)) 
       if (this%spectE%CarryingZeroK) then
           ybuffE(1,1,:) = cmplx(zero,zero,rkind)
       end if 
       !this%w_rhs = this%w_rhs + ybuffE 
       wrhs = wrhs + ybuffE
       
       ! The residual quantity (Gx - <u>)*cos(alpha)/Ro is accomodated in
       ! pressure

       if (this%storeFbody) then
           call this%spectC%ifft(ybuffC1,this%fbody_x)
           call this%spectC%ifft(ybuffC2,this%fbody_y)
           call this%spectE%ifft(ybuffE ,this%fbody_z)
       end if
   end subroutine  

   subroutine addSponge(this)
       class(igrid), intent(inout), target :: this
       complex(rkind), dimension(:,:,:), pointer :: deviationC

       deviationC => this%cbuffyC(:,:,:,1)
       
       deviationC = this%uhat
       if (this%spectC%carryingZeroK) then
           deviationC(1,1,:) = cmplx(zero,zero,rkind)
       end if 
       this%u_rhs = this%u_rhs - (this%RdampC/this%dt)*deviationC

       deviationC = this%vhat
       if (this%spectC%carryingZeroK) then
           deviationC(1,1,:) = cmplx(zero,zero,rkind)
       end if 
       this%v_rhs = this%v_rhs - (this%RdampC/this%dt)*deviationC
       
       this%w_rhs = this%w_rhs - (this%RdampE/this%dt)*this%what ! Mean of w is always zero 

       deviationC = this%That 
       if (this%spectC%carryingZeroK) then
           deviationC(1,1,:) = cmplx(zero,zero,rkind)
       end if 
       this%T_rhs = this%T_rhs - (this%RdampC/this%dt)*deviationC
   end subroutine  
   

   subroutine addNonLinearTerm_Rot(this, u_rhs, v_rhs, w_rhs)
       class(igrid), intent(inout), target :: this
       real(rkind),    dimension(:,:,:), pointer :: dudy, dudz, dudx
       real(rkind),    dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
       real(rkind),    dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
       real(rkind),    dimension(:,:,:), pointer :: dvdzC, dudzC
       real(rkind),    dimension(:,:,:), pointer :: dwdxC, dwdyC
       real(rkind),    dimension(:,:,:), pointer :: T1C, T2C, T1E, T2E 
       complex(rkind), dimension(:,:,:), pointer :: fT1C, fT2C, fT1E, fT2E 
       complex(rkind), dimension(:,:,:), pointer :: tzC, tzE
       complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2), this%sp_gpC%ysz(3)), intent(inout) :: u_rhs, v_rhs
       complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2), this%sp_gpE%ysz(3)), intent(inout) :: w_rhs

       dudx  => this%duidxjC(:,:,:,1); dudy  => this%duidxjC(:,:,:,2); dudzC => this%duidxjC(:,:,:,3); 
       dvdx  => this%duidxjC(:,:,:,4); dvdy  => this%duidxjC(:,:,:,5); dvdzC => this%duidxjC(:,:,:,6); 
       dwdxC => this%duidxjC(:,:,:,7); dwdyC => this%duidxjC(:,:,:,8); dwdz  => this%duidxjC(:,:,:,9); 

       !dwdx => this%duidxjE(:,:,:,1); dwdy => this%duidxjE(:,:,:,2);
       !dudz => this%duidxjE(:,:,:,3); dvdz => this%duidxjE(:,:,:,4);

       dwdx => this%duidxjE(:,:,:,7); dwdy => this%duidxjE(:,:,:,8);
       dudz => this%duidxjE(:,:,:,3); dvdz => this%duidxjE(:,:,:,6);

       T1C => this%rbuffxC(:,:,:,1); T2C => this%rbuffxC(:,:,:,2)
       T1E => this%rbuffxE(:,:,:,1); T2E => this%rbuffxE(:,:,:,2)
      
       fT1C => this%cbuffyC(:,:,:,1); fT2C => this%cbuffyC(:,:,:,2)
       fT1E => this%cbuffyE(:,:,:,1); fT2E => this%cbuffyE(:,:,:,2)
       
       tzC => this%cbuffzC(:,:,:,1); tzE => this%cbuffzE(:,:,:,1)


       T1C = dvdx - dudy
       T1C = T1c*this%v
       call this%spectC%fft(T1C,fT1C)
       T2E = dwdx - dudz
       T2E = T2E*this%w
       call this%spectE%fft(T2E,fT2E)
       call transpose_y_to_z(fT2E,tzE, this%sp_gpE)
       call this%Pade6opZ%interpz_E2C(tzE,tzC,0,0)
       !call transpose_z_to_y(tzC,this%u_rhs, this%sp_gpC)
       call transpose_z_to_y(tzC,u_rhs, this%sp_gpC)
       !this%u_rhs = this%u_rhs + fT1C
       u_rhs = u_rhs + fT1C


       T1C = dudy - dvdx
       T1C = T1C*this%u
       call this%spectC%fft(T1C,fT1C)
       T2E = dwdy - dvdz
       T2E = T2E*this%w
       call this%spectE%fft(T2E,fT2E)
       call transpose_y_to_z(fT2E,tzE, this%sp_gpE)
       call this%Pade6opZ%interpz_E2C(tzE,tzC,0,0)
       !call transpose_z_to_y(tzC,this%v_rhs, this%sp_gpC)
       call transpose_z_to_y(tzC,v_rhs, this%sp_gpC)
       !this%v_rhs = this%v_rhs + fT1C
       v_rhs = v_rhs + fT1C

       T1E = dudz - dwdx
       T1E = T1E*this%uE
       T2E = dvdz - dwdy
       T2E = T2E*this%vE
       T1E = T1E + T2E
       !call this%spectE%fft(T1E,this%w_rhs)
       call this%spectE%fft(T1E,w_rhs)

       if (this%isStratified .or. this%initspinup) then
           T1C = -this%u*this%dTdxC 
           T2C = -this%v*this%dTdyC
           T1C = T1C + T2C
           call this%spectC%fft(T1c,fT1C)
           T1E = -this%w*this%dTdzE
           call this%spectE%fft(T1E,fT1E)
           call transpose_y_to_z(fT2E,tzE, this%sp_gpE)
           call this%Pade6opZ%interpz_E2C(tzE,tzC,WdTdzBC_bottom,WdTdzBC_top)
           call transpose_z_to_y(tzC,this%T_rhs, this%sp_gpC)
           this%T_rhs = this%T_rhs + fT1C
       end if

   end subroutine

   subroutine addNonLinearTerm_skewSymm(this, urhs, vrhs, wrhs)
       class(igrid), intent(inout), target :: this
       real(rkind),    dimension(:,:,:), pointer :: dudy, dudz, dudx
       real(rkind),    dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
       real(rkind),    dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
       real(rkind),    dimension(:,:,:), pointer :: dvdzC, dudzC
       real(rkind),    dimension(:,:,:), pointer :: dwdxC, dwdyC
       real(rkind),    dimension(:,:,:), pointer :: T1C, T2C, T1E, T2E 
       complex(rkind), dimension(:,:,:), pointer :: fT1C, fT2C, fT1E, fT2E 
       complex(rkind), dimension(:,:,:), pointer :: tzC, tzE
       complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout) :: urhs, vrhs
       complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: wrhs

       dudx  => this%duidxjC(:,:,:,1); dudy  => this%duidxjC(:,:,:,2); dudzC => this%duidxjC(:,:,:,3); 
       dvdx  => this%duidxjC(:,:,:,4); dvdy  => this%duidxjC(:,:,:,5); dvdzC => this%duidxjC(:,:,:,6); 
       dwdxC => this%duidxjC(:,:,:,7); dwdyC => this%duidxjC(:,:,:,8); dwdz  => this%duidxjC(:,:,:,9); 

       dwdx => this%duidxjE(:,:,:,7); dwdy => this%duidxjE(:,:,:,8);
       dudz => this%duidxjE(:,:,:,3); dvdz => this%duidxjE(:,:,:,6);

       T1C => this%rbuffxC(:,:,:,1); T2C => this%rbuffxC(:,:,:,2)
       T1E => this%rbuffxE(:,:,:,1); T2E => this%rbuffxE(:,:,:,2)
      
       fT1C => this%cbuffyC(:,:,:,1); fT2C => this%cbuffyC(:,:,:,2)
       fT1E => this%cbuffyE(:,:,:,1); fT2E => this%cbuffyE(:,:,:,2)
       
       tzC => this%cbuffzC(:,:,:,1); tzE => this%cbuffzE(:,:,:,1)


       T1C = dudx*this%u
       T2C = dudy*this%v
       T1C = T1C + T2C
       T1E = dudz*this%w
       call this%spectC%fft(T1C,fT1C)
       call this%spectE%fft(T1E,fT1E)
       call transpose_y_to_z(fT1E,tzE, this%sp_gpE)
       call this%Pade6opZ%interpz_E2C(tzE,tzC,WdUdzBC_bottom,WdUdzBC_top)
       !call transpose_z_to_y(tzC,this%u_rhs, this%sp_gpC)
       call transpose_z_to_y(tzC,urhs, this%sp_gpC)
       !this%u_rhs = this%u_rhs + fT1C
       urhs = urhs + fT1C
       
       T1C = dvdx*this%u
       T2C = dvdy*this%v
       T1C = T1C + T2C
       T1E = dvdz*this%w
       call this%spectC%fft(T1C,fT1C)
       call this%spectE%fft(T1E,fT1E)
       call transpose_y_to_z(fT1E,tzE, this%sp_gpE)
       call this%Pade6opZ%interpz_E2C(tzE,tzC,WdVdzBC_bottom,WdVdzBC_top)
       !call transpose_z_to_y(tzC,this%v_rhs, this%sp_gpC)
       call transpose_z_to_y(tzC,vrhs, this%sp_gpC)
       !this%v_rhs = this%v_rhs + fT1C
       vrhs = vrhs + fT1C
       
       T1E = dwdx*this%uE
       T2E = dwdy*this%vE
       T2E = T1E + T2E
       call this%spectE%fft(T2E,fT2E)
       T1C = dwdz*this%wC
       call this%spectC%fft(T1C,fT1C)
       call transpose_y_to_z(fT1C,tzC, this%sp_gpC)
       call this%Pade6opZ%interpz_C2E(tzC,tzE,WdWdzBC_bottom,WdWdzBC_top)
       !call transpose_z_to_y(tzE,this%w_rhs, this%sp_gpE)
       call transpose_z_to_y(tzE,wrhs, this%sp_gpE)
       !this%w_rhs = this%w_rhs + fT2E
       wrhs = wrhs + fT2E

       T1C = this%u*this%u
       call this%spectC%fft(T1C,fT1C)
       call this%spectC%mtimes_ik1_ip(fT1C)
       !this%u_rhs = this%u_rhs + fT1C
       urhs = urhs + fT1C

       T1C = this%v*this%v
       call this%spectC%fft(T1C,fT1C)
       call this%spectC%mtimes_ik2_ip(fT1C)
       !this%v_rhs = this%v_rhs + fT1C
       vrhs = vrhs + fT1C

       T1C = this%wC*this%wC
       call this%spectC%fft(T1C,fT1C)
       call transpose_y_to_z(fT1C,tzC,this%sp_gpC)
       call this%Pade6opZ%ddz_C2E(tzC,tzE,WWBC_bottom,WWBC_top)
       call transpose_z_to_y(tzE,fT1E,this%sp_gpE)
       !this%w_rhs = this%w_rhs + fT1E
       wrhs = wrhs + fT1E

       T1C = this%u*this%v
       call this%spectC%fft(T1C,fT1C)
       call this%spectC%mtimes_ik2_oop(fT1C,fT2C)
       !this%u_rhs = this%u_rhs + fT2C
       urhs = urhs + fT2C
       call this%spectC%mtimes_ik1_ip(fT1C)
       !this%v_rhs = this%v_rhs + fT1C
       vrhs = vrhs + fT1C

       T1E = this%uE*this%w
       call this%spectE%fft(T1E,fT1E)
       call transpose_y_to_z(fT1E,TzE,this%sp_gpE)
       call this%Pade6opZ%ddz_E2C(tzE,tzC,UWBC_bottom,UWBC_top)
       call transpose_z_to_y(tzC,fT1C,this%sp_gpC)
       ! this%u_rhs = this%u_rhs + fT1C
       urhs = urhs + fT1C
       
       call this%spectE%mtimes_ik1_ip(fT1E)
       !this%w_rhs = this%w_rhs + fT1E
       wrhs = wrhs + fT1E


       T1E = this%vE*this%w
       call this%spectE%fft(T1E,fT1E)
       call transpose_y_to_z(fT1E,TzE,this%sp_gpE)
       call this%Pade6opZ%ddz_E2C(tzE,tzC,VWBC_bottom,VWBC_top)
       call transpose_z_to_y(tzC,fT1C,this%sp_gpC)
       !this%v_rhs = this%v_rhs + fT1C
       vrhs = vrhs + fT1C

       call this%spectE%mtimes_ik2_ip(fT1E)
       !this%w_rhs = this%w_rhs + fT1E
       wrhs = wrhs + fT1E

       !this%u_rhs = -half*this%u_rhs
       !this%v_rhs = -half*this%v_rhs
       !this%w_rhs = -half*this%w_rhs
       urhs = -half*urhs
       vrhs = -half*vrhs
       wrhs = -half*wrhs


       if (this%isStratified .or. this%initspinup) then
           T1C = -this%u*this%T
           call this%spectC%fft(T1C,this%T_rhs)
           call this%spectC%mtimes_ik1_ip(this%T_rhs)
           T1C = -this%v*this%T
           call this%spectC%fft(T1C,fT1C)
           call this%spectC%mtimes_ik2_ip(fT1C)
           this%T_rhs = this%T_rhs + fT1C
           T1E = -this%w * this%TE
           call this%spectE%fft(T1E,fT1E)
           call transpose_y_to_z(fT1E,TzE,this%sp_gpE)
           call this%Pade6opZ%ddz_E2C(tzE,tzC,WTBC_bottom,WTBC_top)
           call transpose_z_to_y(tzC,fT1C,this%sp_gpC)
           this%T_rhs = this%T_rhs + fT1C
       end if 
   end subroutine

   subroutine addExtraForcingTerm(this)
       class(igrid), intent(inout) :: this
       this%u_rhs = this%u_rhs + this%dpF_dxhat
       if (this%storeFbody) then
           this%fbody_x = this%fbody_x + this%dpFdx 
       end if
   end subroutine

   subroutine addBuoyancyTerm(this, urhs, vrhs, wrhs)
       class(igrid), intent(inout), target :: this
       complex(rkind), dimension(:,:,:), pointer :: fT1E, fT1C
       real(rkind), dimension(:,:,:), pointer :: rbuffE
       complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2), this%sp_gpE%ysz(3)), intent(inout) :: wrhs
       complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2), this%sp_gpC%ysz(3)), intent(inout) :: urhs, vrhs
       integer :: mind

       fT1E => this%cbuffyE(:,:,:,1)
       fT1C => this%cbuffyC(:,:,:,1)
       rbuffE => this%rbuffxE(:,:,:,1)

       !fT1E = (this%TEhat)/(this%ThetaRef*this%Fr*this%Fr)
       if(this%useMoisture) then
           mind = this%moistureIndex
           call transpose_y_to_z(this%scalars(mind)%fhat, this%cbuffzC(:,:,:,1), this%sp_gpC)
           call this%Pade6opZ%interpz_C2E(this%cbuffzC(:,:,:,1), this%cbuffzE(:,:,:,1), this%scalars(mind)%BC_bottom, this%scalars(mind)%BC_top)
           call transpose_z_to_y(this%cbuffzE(:,:,:,1), this%cbuffyE(:,:,:,1), this%sp_gpC)
           fT1E = (this%TEhat + this%moistureFactor*this%cbuffyE(:,:,:,1))*this%BuoyancyFact ! See definition of buoyancy factor in init 
       !else
           !fT1E = (this%TEhat)*this%BuoyancyFact ! See definition of buoyancy factor in init 
           wrhs = wrhs + fT1E
           return  ! MOISTURE ADD WILL ONLY WORK IN Z EQUATION
       endif

       !this%w_rhs = this%w_rhs + fT1E 
       select case (this%BuoyancyDirection)
       case(1)
           fT1C = (this%That)*this%BuoyancyFact ! See definition of buoyancy factor in init 
           urhs = urhs + fT1C
       case(2)
           fT1C = (this%That)*this%BuoyancyFact ! See definition of buoyancy factor in init 
           vrhs = vrhs + fT1C
       case(3)
           fT1E = (this%TEhat)*this%BuoyancyFact ! See definition of buoyancy factor in init 
           if (this%spectE%carryingZeroK) then
               fT1E(1,1,:) = cmplx(zero,zero,rkind)
           end if 
           wrhs = wrhs + fT1E
           if (this%storeFbody) then
               call this%spectE%ifft(fT1E, rbuffE)
               this%fbody_z = this%fbody_z + rbuffE
           end if
       end select 


   end subroutine

   
   subroutine addForcedStratification(this)
       class(igrid), intent(inout) :: this 
       real(rkind) :: molecularDiff
       integer :: i, j, k

       molecularDiff = one/(this%Re*this%PrandtlFluid)
       do k = 1,size(this%T_rhs,3)
          do j = 1,size(this%T_rhs,2)
             !$omp simd
             do i = 1,size(this%T_rhs,1)
                this%T_rhs(i,j,k) = this%T_rhs(i,j,k) - molecularDiff*this%uhat(i,j,k)
             end do 
          end do 
       end do 


   end subroutine 

   subroutine addViscousTerm(this, u_rhs, v_rhs, w_rhs)
       class(igrid), intent(inout) :: this
       integer :: i, j, k
       real(rkind) :: oneByRe, molecularDiff
       complex(rkind) :: tmp1, tmp2
       complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2), this%sp_gpC%ysz(3)), intent(inout) :: u_rhs, v_rhs
       complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2), this%sp_gpE%ysz(3)), intent(inout) :: w_rhs

       oneByRe = one/this%Re

       do k = 1,size(this%u_rhs,3)
          do j = 1,size(this%u_rhs,2)
             !$omp simd
             do i = 1,size(this%u_rhs,1)
                 tmp1 = -this%spectC%kabs_sq(i,j,k)*this%uhat(i,j,k) + this%d2udz2hatC(i,j,k)
                 tmp2 = -this%spectC%kabs_sq(i,j,k)*this%vhat(i,j,k) + this%d2vdz2hatC(i,j,k)
                 !this%u_rhs(i,j,k) = this%u_rhs(i,j,k) + oneByRe*tmp1
                 !this%v_rhs(i,j,k) = this%v_rhs(i,j,k) + oneByRe*tmp2
                 u_rhs(i,j,k) = u_rhs(i,j,k) + oneByRe*tmp1
                 v_rhs(i,j,k) = v_rhs(i,j,k) + oneByRe*tmp2
              end do
           end do
        end do

       do k = 1,size(this%w_rhs,3)
          do j = 1,size(this%w_rhs,2)
             !$omp simd
             do i = 1,size(this%w_rhs,1)
                 tmp1 = -this%spectE%kabs_sq(i,j,k)*this%what(i,j,k) + this%d2wdz2hatE(i,j,k)
                 !this%w_rhs(i,j,k) = this%w_rhs(i,j,k) + oneByRe*tmp1
                 w_rhs(i,j,k) = w_rhs(i,j,k) + oneByRe*tmp1
              end do
           end do
        end do

        if (this%isStratified) then
           molecularDiff = one/(this%Re*this%PrandtlFluid)
           do k = 1,size(this%T_rhs,3)
              do j = 1,size(this%T_rhs,2)
                 !$omp simd
                 do i = 1,size(this%T_rhs,1)
                    tmp1 = -this%spectC%kabs_sq(i,j,k)*this%That(i,j,k) + this%d2Tdz2hatC(i,j,k) 
                    this%T_rhs(i,j,k) = this%T_rhs(i,j,k) + molecularDiff*tmp1
                 end do 
              end do 
           end do 
           
        end if

   end subroutine

   pure subroutine S_sponge_smooth(x, output)
      real(rkind), dimension(:,:,:), intent(in)    :: x
      real(rkind), dimension(:,:,:), intent(out)   :: output
      integer :: i, j, k
      real(rkind) :: exparg

      do i = 1,size(x,3)
        do j = 1,size(x,2)
           do k = 1,size(x,1)
               if (x(k,j,i) .le. 0.d0) then
               output(k,j,i) = 0.d0
               else if (x(k,j,i) .ge. 1.d0) then
                   output(k,j,i) = 1.d0
               else
                   exparg = 1.d0/(x(k,j,i) - 1.d0 + 1.0D-32) + 1.d0/(x(k,j,i) + 1.0D-32)
                   exparg = min(exparg,708.0d0) ! overflows if exparg > 709. need a better fix for this
                   output(k,j,i) = 1.d0/(1.d0 + exp(exparg))
               end if
           end do 
        end do 
      end do

   end subroutine
