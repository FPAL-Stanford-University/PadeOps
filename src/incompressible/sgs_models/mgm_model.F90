    subroutine get_MGM_Op(this, tau11, tau12, tau13C, tau22, tau23C, tau33,                      &
                                tau13, tau23,                                                    &
                                 dudx,  dudy,  dudzC,  dvdx,  dvdy,  dvdzC, dwdxC, dwdyC, dwdz,  &
                                 dudz,  dvdz,  dwdx,   dwdy,                                     &
                                 dxsq, dysq, dzsq, mgm_option, maxDissp)
        use constants, only : eps
        class(sgs), intent(inout), target :: this
        real(rkind), dimension(:,:,:), pointer, intent(out) :: tau11, tau12, tau13C, tau22, tau23C, tau33                     ! Cell quantities
        real(rkind), dimension(:,:,:), pointer, intent(out) :: tau13, tau23                                                   ! Edge quantities
        real(rkind), dimension(:,:,:), pointer, intent(in)  :: dudx, dudy, dudzC, dvdx, dvdy, dvdzC, dwdxC, dwdyC, dwdz       ! Cell quantities
        real(rkind), dimension(:,:,:), pointer, intent(in)  :: dudz, dvdz, dwdx, dwdy                                         ! Edge quantities
        real(rkind), intent(in) :: dxsq, dysq, dzsq
        integer, intent(in) :: mgm_option
        real(rkind), intent(out) :: maxDissp

        real(rkind), dimension(:,:,:), pointer :: trG, dissp    ! Cell quantities
        real(rkind), dimension(:,:,:), pointer :: G11E, G22E, G33E, G12E, G13E, G23E, dudxE, dvdyE, dwdzE ! Edge quantities

        integer :: nz
        real(rkind) :: rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,rdum7,rdum8,rdum9,rdum10,rdum11,rdum12

        ! --needed memory MGMbuffs(:,:,:,1:2)
        trG=>this%MGMbuffs(:,:,:,1); dissp=>this%MGMbuffs(:,:,:,2)

        ! G_ij is stored in tau_ij
        tau11  = dxsq * dudx * dudx  + dysq * dudy * dudy  + dzsq * dudzC * dudzC
        tau12  = dxsq * dudx * dvdx  + dysq * dudy * dvdy  + dzsq * dudzC * dvdzC 
        tau13C = dxsq * dudx * dwdxC + dysq * dudy * dwdyC + dzsq * dudzC * dwdz
        tau22  = dxsq * dvdy * dvdx  + dysq * dvdy * dvdy  + dzsq * dvdzC * dvdzC 
        tau23C = dxsq * dvdy * dwdxC + dysq * dvdy * dwdyC + dzsq * dvdzC * dwdz 
        tau33  = dxsq * dwdz * dwdxC + dysq * dwdyC* dwdyC + dzsq * dwdz  * dwdz

        print*, "--------------------------------------"
        print*, dudx(1,1,1), dudy(1,1,1), dudz(1,1,1)
        print*, dvdx(1,1,1), dvdy(1,1,1), dvdz(1,1,1)
        print*, dwdx(1,1,1), dwdy(1,1,1), dwdz(1,1,1)
        print*, "--------------------------------------"

        print*, tau11(1,1,1), tau12(1,1,1), tau13C(1,1,1)
        print*, tau22(1,1,1), tau23C(1,1,1)
        print*, tau33(1,1,1)

        print*, "--------------------------------------"
        !rdum1 = p_maxval(tau11); rdum2 = p_minval(tau11)
        !rdum3 = p_maxval(tau12); rdum4 = p_minval(tau12)
        !rdum5 = p_maxval(tau22); rdum6 = p_minval(tau22)
        !rdum7 = p_maxval(tau33); rdum8 = p_minval(tau33)
        !rdum9 = p_maxval(tau13C); rdum10 = p_minval(tau13C)
        !rdum11 = p_maxval(tau23C); rdum12 = p_minval(tau23C)
        !if(nrank==0) then
        !  write(*,*) '------------Beginnning MGM---------------'
        !  write(*,*) '--dx, dy, dz--', dxsq, dysq, dzsq
        !  write(*,*) '--G11--', rdum1, rdum2
        !  write(*,*) '--G12--', rdum3, rdum4
        !  write(*,*) '--G22--', rdum5, rdum6
        !  write(*,*) '--G33--', rdum7, rdum8
        !  write(*,*) '--G13--', rdum9, rdum10
        !  write(*,*) '--G23--', rdum11, rdum12
        !endif

        !rdum1 = p_maxval((dudx));  rdum2 = p_maxval((dudy));  rdum3 = p_maxval((dudzC))
        !rdum4 = p_maxval((dvdx));  rdum5 = p_maxval((dvdy));  rdum6 = p_maxval((dvdzC))
        !rdum7 = p_maxval((dwdxC)); rdum8 = p_maxval((dwdyC)); rdum9 = p_maxval((dwdz))
        !if(nrank==0) then
        !  write(*,'(a,3(1x,e19.12))') '--dudx, dudy, dudz--', rdum1, rdum2, rdum3
        !  write(*,'(a,3(1x,e19.12))') '--dvdx, dvdy, dvdz--', rdum4, rdum5, rdum6
        !  write(*,'(a,3(1x,e19.12))') '--dwdx, dwdy, dwdz--', rdum7, rdum8, rdum9
        !  write(*,*) '------------Beginnning MGM---------------'
        !endif

        trG = tau11 + tau22 + tau33         ! actually this is trace of G
        trG = max(trG, eps)                 ! to guard against division by zero

        ! compute G_ij * S_ij
        dissp = -tau13C * (dudzC + dwdx)
        dissp = dissp -(tau12 * (dudy + dvdx) + tau23C * (dvdzC + dwdxC))
        dissp = dissp - (tau11*dudx + tau22*dvdy + tau33*dwdz)
        dissp = max(dissp, zero)   ! to preclude backscatter
        rdum3 = p_maxval(dissp)!; rdum4 = p_minval(dissp)
        maxDissp = rdum3!  max(rdum3, rdum4)
        !if(nrank==0) then
        !  write(*,*) '------------Dissp---------------'
        !  write(*,*) '--Bef backscat--', rdum1, rdum2
        !  write(*,*) '--Aft backscat--', rdum3, rdum4
        !  !write(*,*) '--mconst--',  this%mconst
        !  write(*,*) '------------Dissp---------------'
        !endif

        dissp = dissp/trG     ! note tau13C is trace of G
        dissp = dissp*dissp
        dissp = dissp/trG

        tau11 = this%mconst * dissp * tau11
        tau12 = this%mconst * dissp * tau12
        tau22 = this%mconst * dissp * tau22
        tau33 = this%mconst * dissp * tau33

            print*, this%mconst
        if(mgm_option == 1) then
            ! ------ option 1 -------!
            ! compute tau13C and tau23C and then interpolate to edges
            tau13C = this%mconst * dissp * tau13C 
            tau23C = this%mconst * dissp * tau23C
            

            ! interpolate(tau13C, tau13)
              call transpose_x_to_y(tau13C,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              nz = size(this%rtmpZ,3)
              this%rtmpzE(:,:,nz+1) = two*this%rtmpZ(:,:,nz) - this%rtmpzE(:,:,nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,tau13,this%gpE)
    
            ! interpolate(tau23C, tau23)
              call transpose_x_to_y(tau23C,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              nz = size(this%rtmpZ,3)
              this%rtmpzE(:,:,nz+1) = two*this%rtmpZ(:,:,nz) - this%rtmpzE(:,:,nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,tau23,this%gpE)
            ! ------ option 1 -------!
  
        else if(mgm_option == 2) then
            ! ------ option 2 -------!
            ! alternatively, interpolate dissp, trG, G13 and G23 to edges and compute tau13
            ! --needed memory: this%MGMbuffsE(:,:,:,1)
    
            ! interpolate(trG -> tau13) ! contains trGE
              call transpose_x_to_y(trG,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              nz = size(this%rtmpZ,3)
              this%rtmpzE(:,:,nz+1) = two*this%rtmpZ(:,:,nz) - this%rtmpzE(:,:,nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,tau13,this%gpE)
    
            dissp = dissp*trG
            ! interpolate(dissp -> tau23)  ! contains disspE
              call transpose_x_to_y(dissp,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              nz = size(this%rtmpZ,3)
              this%rtmpzE(:,:,nz+1) = two*this%rtmpZ(:,:,nz) - this%rtmpzE(:,:,nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,tau23,this%gpE)
              tau23 = max(tau23, zero)   ! to preclude backscatter
              rdum3 = p_maxval(tau23)!; rdum4 = p_minval(tau23)
              maxDissp = max(maxDissp, rdum3)
              !if(nrank==0) then
              !  write(*,*) '------------DisspE Opt2---------'
              !  write(*,*) '--Bef backscat--', rdum1, rdum2
              !  write(*,*) '--Aft backscat--', rdum3, rdum4
              !  write(*,*) '------------DisspE Opt2---------'
              !endif

            tau23 = tau23 / (tau13 + eps)        ! contains disspE/trGE; guard against division by zero
    
            ! interpolate(G13 -> this%MGMbuffsE(:,:,:,1))     ! tau13C contains G13; this%MGMbuffsE(:,:,:,1) contains G13E
              call transpose_x_to_y(tau13C,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              nz = size(this%rtmpZ,3)
              this%rtmpzE(:,:,nz+1) = two*this%rtmpZ(:,:,nz) - this%rtmpzE(:,:,nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,this%MGMbuffsE(:,:,:,1),this%gpE)
            tau13 = this%mconst * tau23 * this%MGMbuffsE(:,:,:,1)
    
            ! interpolate(G23 -> this%MGMbuffsE(:,:,:,1))     ! tau23C contains G23; this%MGMbuffsE(:,:,:,1) contains G23E
              call transpose_x_to_y(tau23C,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              nz = size(this%rtmpZ,3)
              this%rtmpzE(:,:,nz+1) = two*this%rtmpZ(:,:,nz) - this%rtmpzE(:,:,nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,this%MGMbuffsE(:,:,:,1),this%gpE)
            tau23 = this%mconst * tau23 * this%MGMbuffsE(:,:,:,1)

            ! ------ option 2 -------!
  
        else if(mgm_option == 3) then
            ! ------ option 3 -------!
            ! alternatively, interpolate everything needed to edges and compute tau13, tau23
            ! --needed memory: this%MGMbuffsE(:,:,:,1:9)
            dudxE => this%MGMbuffsE(:,:,:,1);  dvdyE => this%MGMbuffsE(:,:,:,2);  dwdzE => this%MGMbuffsE(:,:,:,3);
            G11E  => this%MGMbuffsE(:,:,:,4);  G12E  => this%MGMbuffsE(:,:,:,5);  G13E  => this%MGMbuffsE(:,:,:,6)
            G22E  => this%MGMbuffsE(:,:,:,7);  G23E  => this%MGMbuffsE(:,:,:,8);  G33E  => this%MGMbuffsE(:,:,:,9)
  
            ! interpolate(dudx -> dudxE)
              call transpose_x_to_y(dudx,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              nz = size(this%rtmpZ,3)
              this%rtmpzE(:,:,nz+1) = two*this%rtmpZ(:,:,nz) - this%rtmpzE(:,:,nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,dudxE,this%gpE)
  
            ! interpolate(dudy -> tau23)    ! tau23 is dudyE here
              call transpose_x_to_y(dudy,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              nz = size(this%rtmpZ,3)
              this%rtmpzE(:,:,nz+1) = two*this%rtmpZ(:,:,nz) - this%rtmpzE(:,:,nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,tau23,this%gpE)
  
            ! interpolate(dvdx -> tau13)    ! tau13 is dvdxE here
              call transpose_x_to_y(dvdx,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              nz = size(this%rtmpZ,3)
              this%rtmpzE(:,:,nz+1) = two*this%rtmpZ(:,:,nz) - this%rtmpzE(:,:,nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,tau13,this%gpE)
  
            ! interpolate(dvdy -> dvdyE)
              call transpose_x_to_y(dvdy,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              nz = size(this%rtmpZ,3)
              this%rtmpzE(:,:,nz+1) = two*this%rtmpZ(:,:,nz) - this%rtmpzE(:,:,nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,dvdyE,this%gpE)
  
            ! interpolate(dwdz -> dwdzE)
              call transpose_x_to_y(dwdz,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              nz = size(this%rtmpZ,3)
              this%rtmpzE(:,:,nz+1) = two*this%rtmpZ(:,:,nz) - this%rtmpzE(:,:,nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,dwdzE,this%gpE)
  
            G11E = dxsq * dudxE * dudxE + dysq * tau23 * tau23 + dzsq * dudz  * dudz
            G12E = dxsq * dudxE * tau13 + dysq * tau23 * dvdyE + dzsq * dudz  * dvdz 
            G13E = dxsq * dudxE * dwdx  + dysq * tau23 * dwdy  + dzsq * dudz  * dwdzE
            G22E = dxsq * dvdyE * tau13 + dysq * dvdyE * dvdyE + dzsq * dvdz  * dvdz 
            G23E = dxsq * dvdyE * dwdx  + dysq * dvdyE * dwdy  + dzsq * dvdz  * dwdzE 
            G33E = dxsq * dwdzE * dwdx  + dysq * dwdy  * dwdy  + dzsq * dwdzE * dwdzE
  
            ! tau23 actually contains disspE now; tau23 is overwritten, no longer have dudyE
            tau23 = -G12E*(tau23 + tau13)
            tau23 = tau23 - (G13E*(dudz + dwdx) + G23E*(dvdz + dwdy))
            tau23 = tau23 - (G11E*dudxE + G22E*dvdyE + G33E*dwdzE)
            tau23 = max(tau23, zero)   ! to preclude backscatter
            rdum3 = p_maxval(tau23)!; rdum4 = p_minval(tau23)
            maxDissp = max(maxDissp, rdum3)
            !if(nrank==0) then
            !  write(*,*) '------------DisspE Opt3---------'
            !  write(*,*) '--Bef backscat--', rdum1, rdum2
            !  write(*,*) '--Aft backscat--', rdum3, rdum4
            !  write(*,*) '------------DisspE Opt3---------'
            !endif

            ! tau13 actually contains trace of G now; tau13 is overwritten, no longer have dvdxE
            tau13 = G11E + G22E + G33E
            tau13 = max(tau13, eps)        ! to guard against division by zero
  
            ! note tau13 is trace of G; tau23 is disspE
            tau23 = tau23/tau13
            tau23 = tau23*tau23
            tau23 = tau23/tau13

            tau13 = this%mconst * tau23 * G13E     ! tau23 contains disspE/trGE
            tau23 = this%mconst * tau23 * G23E
  
            ! Done
            nullify(G11E,G12E,G13E,G22E,G23E,G33E,dudxE,dvdyE,dwdzE)
            ! ------ option 3 -------!

        endif

        nullify(trG,dissp)

        !rdum1 = p_maxval(tau11); rdum2 = p_minval(tau11)
        !rdum3 = p_maxval(tau12); rdum4 = p_minval(tau12)
        !rdum5 = p_maxval(tau22); rdum6 = p_minval(tau22)
        !rdum7 = p_maxval(tau33); rdum8 = p_minval(tau33)
        !rdum9 = p_maxval(tau13); rdum10 = p_minval(tau13)
        !rdum11 = p_maxval(tau23); rdum12 = p_minval(tau23)
        !if(nrank==0) then
        !  write(*,*) '--In MGM--', mgm_option
        !  write(*,*) '--tau11--', rdum1, rdum2
        !  write(*,*) '--tau12--', rdum3, rdum4
        !  write(*,*) '--tau22--', rdum5, rdum6
        !  write(*,*) '--tau33--', rdum7, rdum8
        !  write(*,*) '--tau13--', rdum9, rdum10
        !  write(*,*) '--tau23--', rdum11, rdum12
        !endif

    end subroutine 
