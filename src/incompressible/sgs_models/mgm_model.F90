    subroutine get_MGM_Op(this, duidxjC, duidxjE)!, maxDissp)
        use constants, only : eps
        class(sgs), intent(inout), target :: this
        real(rkind)   , dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),9), intent(in), target :: duidxjC
        real(rkind)   , dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),4), intent(in), target :: duidxjE
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz
        real(rkind), dimension(:,:,:), pointer :: dvdx, dvdy, dvdz
        real(rkind), dimension(:,:,:), pointer :: dwdx, dwdy, dwdz
        real(rkind), dimension(:,:,:), pointer :: dvdzC, dudzC
        real(rkind), dimension(:,:,:), pointer :: dwdxC, dwdyC
        real(rkind), dimension(:,:,:), pointer :: tau11, tau12, tau13, tau13C, tau22, tau23, tau23C, tau33
        !real(rkind), intent(out) :: maxDissp

        real(rkind), dimension(:,:,:), pointer :: trG, dissp, disspE    ! Cell quantities
        real(rkind), dimension(:,:,:), pointer :: G13E, G23E, dudxE, dvdyE, dwdzE ! Edge quantities

        real(rkind) :: rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,rdum7,rdum8,rdum9,rdum10,rdum11,rdum12

        dudx  => duidxjC(:,:,:,1); dudy  => duidxjC(:,:,:,2); dudzC => duidxjC(:,:,:,3);
        dvdx  => duidxjC(:,:,:,4); dvdy  => duidxjC(:,:,:,5); dvdzC => duidxjC(:,:,:,6);
        dwdxC => duidxjC(:,:,:,7); dwdyC => duidxjC(:,:,:,8); dwdz  => duidxjC(:,:,:,9); 
        dwdx => duidxjE(:,:,:,1); dwdy => duidxjE(:,:,:,2);
        dudz => duidxjE(:,:,:,3); dvdz => duidxjE(:,:,:,4);


        tau11 => this%rbuff(:,:,:,1); tau12 => this%rbuff(:,:,:,2) ; tau13C => this%rbuff(:,:,:,3)
        tau22 => this%rbuff(:,:,:,4); tau23C => this%rbuff(:,:,:,5); tau33 => this%rbuff(:,:,:,6)
        tau13 => this%rbuffE(:,:,:,1); tau23 => this%rbuffE(:,:,:,2) 



        trG=>this%rbuff(:,:,:,7); dissp=>this%rbuff(:,:,:,8)
        disspE=>this%rbuffE(:,:,:,3)

        dudxE => this%MGMbuffsE(:,:,:,1);  dvdyE => this%MGMbuffsE(:,:,:,2);  dwdzE => this%MGMbuffsE(:,:,:,3);
        G13E  => this%MGMbuffsE(:,:,:,4); G23E  => this%MGMbuffsE(:,:,:,5)

        ! G_ij is stored in tau_ij
        tau11  = this%dxsq * dudx * dudx  + this%dysq * dudy * dudy  + this%dzsq * dudzC * dudzC
        tau12  = this%dxsq * dudx * dvdx  + this%dysq * dudy * dvdy  + this%dzsq * dudzC * dvdzC 
        tau13C = this%dxsq * dudx * dwdxC + this%dysq * dudy * dwdyC + this%dzsq * dudzC * dwdz
        tau22  = this%dxsq * dvdx * dvdx  + this%dysq * dvdy * dvdy  + this%dzsq * dvdzC * dvdzC 
        tau23C = this%dxsq * dvdx * dwdxC + this%dysq * dvdy * dwdyC + this%dzsq * dvdzC * dwdz 
        tau33  = this%dxsq * dwdxC * dwdxC + this%dysq * dwdyC* dwdyC + this%dzsq * dwdz  * dwdz

        !print *, "G11", tau11(3,2,3)
        !print *, "G22", tau22(3,2,3)
        !print *, "G33", tau33(3,2,3)
        !print *, "G12", tau12(3,2,3)
        !print *, "G13", tau13C(3,2,3)
        !print *, "G23", tau23C(3,2,3)
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
        dissp = -tau13C * (dudzC + dwdxC)
        dissp = dissp -(tau12 * (dudy + dvdx) + tau23C * (dvdzC + dwdyC))
        dissp = dissp - (tau11*dudx + tau22*dvdy + tau33*dwdz)
        dissp = max(dissp, zero)   ! to preclude backscatter
        rdum3 = maxval(dissp)
        !maxDissp = p_maxval(rdum3)!; rdum4 = p_minval(dissp)
        !if(nrank==0) then
        !  write(*,*) '------------Dissp---------------'
        !  write(*,*) '--Bef backscat--', rdum1, rdum2
        !  write(*,*) '--Aft backscat--', rdum3, rdum4
        !  !write(*,*) '--mconst--',  this%mconst
        !  write(*,*) '------------Dissp---------------'
        !endif

        dissp = dissp**2     ! note tau13C is trace of G
        trG = trG**3
        dissp = dissp/trG     ! note tau13C is trace of G
        
        tau11 = this%mconst * dissp * tau11
        tau12 = this%mconst * dissp * tau12
        tau22 = this%mconst * dissp * tau22
        tau33 = this%mconst * dissp * tau33


            ! ------ option 3 -------!
            ! alternatively, interpolate everything needed to edges and compute tau13, tau23
            ! --needed memory: this%MGMbuffsE(:,:,:,1:9)
  
            ! interpolate(dudx -> dudxE)
              call transpose_x_to_y(dudx,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              this%rtmpzE(:,:,this%nz+1) = two*this%rtmpZ(:,:,this%nz) - this%rtmpzE(:,:,this%nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,dudxE,this%gpE)
  
            ! interpolate(dudy -> tau23)    ! tau23 is dudyE here
              call transpose_x_to_y(dudy,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              this%rtmpzE(:,:,this%nz+1) = two*this%rtmpZ(:,:,this%nz) - this%rtmpzE(:,:,this%nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,tau23,this%gpE)
  
            ! interpolate(dvdx -> tau13)    ! tau13 is dvdxE here
              call transpose_x_to_y(dvdx,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              this%rtmpzE(:,:,this%nz+1) = two*this%rtmpZ(:,:,this%nz) - this%rtmpzE(:,:,this%nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,tau13,this%gpE)
  
            ! interpolate(dvdy -> dvdyE)
              call transpose_x_to_y(dvdy,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              this%rtmpzE(:,:,this%nz+1) = two*this%rtmpZ(:,:,this%nz) - this%rtmpzE(:,:,this%nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,dvdyE,this%gpE)
  
            ! interpolate(dwdz -> dwdzE)
              call transpose_x_to_y(dwdz,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              this%rtmpzE(:,:,this%nz+1) = two*this%rtmpZ(:,:,this%nz) - this%rtmpzE(:,:,this%nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,dwdzE,this%gpE)
  
            G13E = this%dxsq * dudxE * dwdx  + this%dysq * tau23 * dwdy  + this%dzsq * dudz  * dwdzE
            G23E = this%dxsq * tau13 * dwdx  + this%dysq * dvdyE * dwdy  + this%dzsq * dvdz  * dwdzE 
 
            ! interpolate(dissp -> disspE)
              call transpose_x_to_y(dissp,this%rtmpY,this%gpC)
              call transpose_y_to_z(this%rtmpY,this%rtmpZ,this%gpC)
              call this%Ops2ndOrder%InterpZ_Cell2Edge(this%rtmpZ,this%rtmpZE,zero,zero)
              this%rtmpzE(:,:,this%nz+1) = two*this%rtmpZ(:,:,this%nz) - this%rtmpzE(:,:,this%nz)
              call transpose_z_to_y(this%rtmpzE,this%rtmpYE,this%gpE)
              call transpose_y_to_x(this%rtmpYE,disspE,this%gpE)
  

            tau13 = this%mconst * disspE * G13E
            tau23 = this%mconst * disspE * G23E
 
            ! Done
            nullify(G13E,G23E,dudxE,dvdyE,dwdzE)
            ! ------ option 3 -------!

        nullify(trG,dissp)


    end subroutine 
