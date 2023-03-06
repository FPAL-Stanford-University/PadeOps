!subroutine init_dynamic_procedure(this, testfilter_x, testfilter_y, testfilter_z, ybuf, DeltaRatio)
subroutine init_dynamic_procedure(this, testfilter_x, testfilter_y, testfilter_z, DeltaRatio)
  class(sgs_cgrid), target, intent(inout) :: this
  character(len=clen), intent(in) :: testfilter_x, testfilter_y, testfilter_z
  !real(rkind), dimension(:,:,:,:), intent(in), target :: ybuf
  real(rkind), intent(in) :: DeltaRatio

   if(this%DynamicProcedureType==1) then
      allocate(this%cmodel_local(this%nyL))
      allocate(this%cmodel_local_Qjsgs(this%nyL))
      allocate(this%cmodel_local_tke(this%nyL))
      allocate(this%ybuflocal(this%nxL, this%nyL, this%nzL, 34))

      this%uFil => this%ybuf(:,:,:,1)
      this%vFil => this%ybuf(:,:,:,2)
      this%wFil => this%ybuf(:,:,:,3)
      this%rhoFil  => this%ybuf(:,:,:,4)
      this%numer   => this%ybuf(:,:,:,5)
      this%denom   => this%ybuf(:,:,:,6)

      this%duiFildxj => this%ybuflocal(:,:,:,1:9)
      this%duFildx   => this%ybuflocal(:,:,:,1)
      this%duFildy   => this%ybuflocal(:,:,:,2)
      this%duFildz   => this%ybuflocal(:,:,:,3)
      this%dvFildx   => this%ybuflocal(:,:,:,4)
      this%dvFildy   => this%ybuflocal(:,:,:,5)
      this%dvFildz   => this%ybuflocal(:,:,:,6)
      this%dwFildx   => this%ybuflocal(:,:,:,7)
      this%dwFildy   => this%ybuflocal(:,:,:,8)
      this%dwFildz   => this%ybuflocal(:,:,:,9)
      this%SFil_ij   => this%ybuflocal(:,:,:,10:15)
      this%tausgsFil => this%ybuflocal(:,:,:,16:21)
      this%Lij       => this%ybuflocal(:,:,:,22)
      this%Mij       => this%ybuflocal(:,:,:,23)
      this%TFil      => this%ybuflocal(:,:,:,24)

      this%gradTFil  => this%ybuflocal(:,:,:,25:27)
      this%dTFildx   => this%ybuflocal(:,:,:,25)
      this%dTFildy   => this%ybuflocal(:,:,:,26)
      this%dTFildz   => this%ybuflocal(:,:,:,27)
      this%QjsgsFil  => this%ybuflocal(:,:,:,28:30)

      this%nusgsFil   => this%ybuflocal(:,:,:,31)
      this%modSFil_sq => this%ybuflocal(:,:,:,32) 
      this%SiiFil     => this%ybuflocal(:,:,:,33)
      this%qtkeFil     => this%ybuflocal(:,:,:,34)
      this%DeltaRatioSq = DeltaRatio * DeltaRatio

   endif

   ! Allocate testFil
   if ( allocated(this%testfil) ) deallocate(this%testfil)
   allocate(this%testfil)
   call this%testfil%init(this%decomp, this%periodicx, this%periodicy, this%periodicz, testfilter_x, testfilter_y, testfilter_z) 

end subroutine

subroutine destroy_dynamic_procedure(this)
  class(sgs_cgrid), intent(inout) :: this

   if ( allocated(this%testfil) ) deallocate(this%testfil)

   if(this%DynamicProcedureType==1) then
      nullify(this%QjsgsFil, this%dTFildz, this%dTFildy, this%dTFildx)
      nullify(this%gradTFil, this%TFil, this%Mij, this%Lij, this%tausgsFil)
      nullify(this%SFil_ij, this%dwFildz, this%dwFildy, this%dwFildx)
      nullify(this%dvFildz, this%dvFildy, this%dvFildx)
      nullify(this%duFildz, this%duFildy, this%duFildx)
      nullify(this%duiFildxj, this%denom, this%numer)
      nullify(this%rhoFil, this%wFil, this%vFil, this%uFil)
      nullify(this%nusgsFil, this%SiiFil, this%modSFil_sq, this%qtkeFil)

      if(allocated(this%ybuflocal))          deallocate(this%ybuflocal)
      if(allocated(this%cmodel_local_Qjsgs)) deallocate(this%cmodel_local_Qjsgs)
      if(allocated(this%cmodel_local_tke))   deallocate(this%cmodel_local_tke)
      if(allocated(this%cmodel_local))       deallocate(this%cmodel_local)
   endif

end subroutine

subroutine gradient(this, f, dfdx, dfdy, dfdz, x_bc, y_bc, z_bc)
    class(sgs_cgrid),target, intent(inout) :: this
    real(rkind), intent(in ), dimension(this%nxL, this%nyL, this%nzL) :: f
    real(rkind), intent(out), dimension(this%nxL, this%nyL, this%nzL) :: dfdx
    real(rkind), intent(out), dimension(this%nxL, this%nyL, this%nzL) :: dfdy
    real(rkind), intent(out), dimension(this%nxL, this%nyL, this%nzL) :: dfdz
    integer, dimension(2), optional, intent(in) :: x_bc, y_bc, z_bc

    real(rkind), dimension(:,:,:), pointer :: tmp1_in_x, tmp2_in_x, tmp1_in_z, tmp2_in_z
    integer :: lastx, lastz
    
    ! Allocate pointers for the needed buffers 
    ! 2 buffers in x and z are needed
    lastx = size(this%xbuf,4)
    tmp1_in_x => this%xbuf(:,:,:,lastx)
    tmp2_in_x => this%xbuf(:,:,:,lastx-1)

    lastz = size(this%zbuf,4)
    tmp1_in_z => this%zbuf(:,:,:,lastz)
    tmp2_in_z => this%zbuf(:,:,:,lastz-1)
 
    ! Get Y derivatives
    call this%der%ddy(f,dfdy,y_bc(1),y_bc(2))

    ! Get X derivatives
    call transpose_y_to_x(f, tmp1_in_x, this%decomp)
    call this%der%ddx(tmp1_in_x, tmp2_in_x, x_bc(1), x_bc(2))
    call transpose_x_to_y(tmp2_in_x, dfdx)

    ! Get Z derivatives
    call transpose_y_to_z(f, tmp1_in_z, this%decomp)
    call this%der%ddz(tmp1_in_z, tmp2_in_z, z_bc(1), z_bc(2))
    call transpose_z_to_y(tmp2_in_z, dfdz)

end subroutine 

subroutine doLocalDynamicProcedure_mgm(this, rho, u, v, w, tausgs)
  class(sgs_cgrid), intent(inout) :: this
  real(rkind), dimension(this%nxL, this%nyL, this%nzL),    intent(in)  :: rho, u, v, w
  real(rkind), dimension(this%nxL, this%nyL, this%nzL, 6), intent(in)  :: tausgs

  integer :: j

  ! space needed: rhoFil, uFil, vFil, wFil, duiFildxj (9), SFil_ij (6), 
 
  ! compute filtered velocities -- this is needed to compute \hat{\bar{S}}_ij
  call this%filter_xyz(u, this%uFil);
  call this%filter_xyz(v, this%vFil);
  call this%filter_xyz(w, this%wFil);
  call this%filter_xyz(rho, this%rhoFil);

  ! compute gradient; Sij; model kernel
  call this%gradient(this%uFil,this%duFildx,this%duFildy,this%duFildz,-this%x_bc, this%y_bc, this%z_bc)
  call this%gradient(this%vFil,this%dvFildx,this%dvFildy,this%dvFildz, this%x_bc,-this%y_bc, this%z_bc)
  call this%gradient(this%wFil,this%dwFildx,this%dwFildy,this%dwFildz, this%x_bc, this%y_bc,-this%z_bc)
  call this%get_Sij_from_duidxj(this%duiFildxj, this%SFil_ij)
  call this%get_mgm_kernel(this%rhoFil, this%duiFildxj, this%SFil_ij, this%tausgsFil)
  this%tausgsFil = this%deltaRatioSq * this%tausgsFil
 
  ! compute filtered density-velocity products \hat{rho*u} -- this is needed to compute L_ij
  call this%filter_xyz(rho*u, this%uFil);
  call this%filter_xyz(rho*v, this%vFil);
  call this%filter_xyz(rho*w, this%wFil);
 
  ! ---11 component--
  call this%filter_xyz(rho*u*u, this%Lij);             this%Lij = this%Lij - this%uFil*this%uFil/this%rhoFil
  call this%filter_xyz(tausgs(:,:,:,1), this%Mij);     this%Mij = this%Mij - this%tausgsFil(:,:,:,1)
  this%numer = this%Lij*this%Mij;                      this%denom = this%Mij*this%Mij

  ! ---12 component--
  call this%filter_xyz(rho*u*v, this%Lij);             this%Lij = this%Lij - this%uFil*this%vFil/this%rhoFil
  call this%filter_xyz(tausgs(:,:,:,2), this%Mij);     this%Mij = this%Mij - this%tausgsFil(:,:,:,2)
  this%numer   = this%numer + two*this%Lij*this%Mij;   this%denom = this%denom + two*this%Mij*this%Mij

  ! ---13 component--
  call this%filter_xyz(rho*u*w, this%Lij);             this%Lij = this%Lij - this%uFil*this%wFil/this%rhoFil
  call this%filter_xyz(tausgs(:,:,:,3), this%Mij);     this%Mij = this%Mij - this%tausgsFil(:,:,:,3)
  this%numer   = this%numer + two*this%Lij*this%Mij;   this%denom = this%denom + two*this%Mij*this%Mij

  ! ---22 component--
  call this%filter_xyz(rho*v*v, this%Lij);             this%Lij = this%Lij - this%vFil*this%vFil/this%rhoFil
  call this%filter_xyz(tausgs(:,:,:,4), this%Mij);     this%Mij = this%Mij - this%tausgsFil(:,:,:,4)
  this%numer   = this%numer + this%Lij*this%Mij;       this%denom = this%denom + this%Mij*this%Mij

  ! ---23 component--
  call this%filter_xyz(rho*v*w, this%Lij);             this%Lij = this%Lij - this%vFil*this%wFil/this%rhoFil
  call this%filter_xyz(tausgs(:,:,:,5), this%Mij);     this%Mij = this%Mij - this%tausgsFil(:,:,:,5)
  this%numer   = this%numer + two*this%Lij*this%Mij;   this%denom = this%denom + two*this%Mij*this%Mij

  ! ---33 component--
  call this%filter_xyz(rho*w*w, this%Lij);             this%Lij = this%Lij - this%wFil*this%wFil/this%rhoFil
  call this%filter_xyz(tausgs(:,:,:,6), this%Mij);     this%Mij = this%Mij - this%tausgsFil(:,:,:,6)
  this%numer   = this%numer + this%Lij*this%Mij;       this%denom = this%denom + this%Mij*this%Mij

  do j = 1, this%nyL
      this%cmodel_local(j) = max(-p_sum(sum(this%numer(:,j,:))), zero) / (p_sum(sum(this%denom(:,j,:))) + 1.0d-18)
  enddo

  !! allows reuse of several filtered fields for calcualtion of Qjsgs coefficient
  this%preComputed_SFil_duFil = .true.

end subroutine

subroutine doLocalDynamicProcedure_mgm_Qjsgs(this, rho, u, v, w, T, Qjsgs)
  class(sgs_cgrid), intent(inout) :: this
  real(rkind), dimension(this%nxL, this%nyL, this%nzL),    intent(in)  :: rho, u, v, w, T
  real(rkind), dimension(this%nxL, this%nyL, this%nzL, 3), intent(in)  :: Qjsgs

  integer :: j

  if(.not. this%preComputed_SFil_duFil) then
      ! compute filtered velocities -- this is needed to compute \hat{\bar{S}}_ij
      call this%filter_xyz(u, this%uFil);
      call this%filter_xyz(v, this%vFil);
      call this%filter_xyz(w, this%wFil);
      call this%filter_xyz(rho, this%rhoFil);

      ! compute gradient; Sij; tausgs model kernel
      call this%gradient(this%uFil,this%duFildx,this%duFildy,this%duFildz,-this%x_bc, this%y_bc, this%z_bc)
      call this%gradient(this%vFil,this%dvFildx,this%dvFildy,this%dvFildz, this%x_bc,-this%y_bc, this%z_bc)
      call this%gradient(this%wFil,this%dwFildx,this%dwFildy,this%dwFildz, this%x_bc, this%y_bc,-this%z_bc)
      call this%get_Sij_from_duidxj(this%duiFildxj, this%SFil_ij)
  endif

  ! space needed: TFil, gradTFil (3), QjsgsFil (3), 
  ! compute filtered Temperature - this is needed to compute \hat{\gradT}_j
  call this%filter_xyz(T, this%TFil);
  call this%gradient(this%TFil,this%dTFildx,this%dTFildy,this%dTFildz, this%x_bc, this%y_bc, this%z_bc)
  call this%get_Qjsgs_mgm_kernel(this%rhoFil, this%duiFildxj, this%SFil_ij, this%gradTFil, this%QjsgsFil)
  this%QjsgsFil = this%deltaRatioSq * this%QjsgsFil
 
  ! compute filtered density-velocity products \hat{rho*u} -- this is needed to compute L_ij
  call this%filter_xyz(rho*u, this%uFil);
  call this%filter_xyz(rho*v, this%vFil);
  call this%filter_xyz(rho*w, this%wFil);
  call this%filter_xyz(rho*T, this%TFil);
 
  ! ---1 component--
  call this%filter_xyz(rho*u*T, this%Lij);            this%Lij = this%Lij - this%uFil*this%TFil/this%rhoFil
  call this%filter_xyz(Qjsgs(:,:,:,1), this%Mij);     this%Mij = this%Mij - this%QjsgsFil(:,:,:,1)
  this%numer = this%Lij*this%Mij;                     this%denom = this%Mij*this%Mij

  ! ---2 component--
  call this%filter_xyz(rho*v*T, this%Lij);            this%Lij = this%Lij - this%vFil*this%TFil/this%rhoFil
  call this%filter_xyz(Qjsgs(:,:,:,2), this%Mij);     this%Mij = this%Mij - this%QjsgsFil(:,:,:,2)
  this%numer = this%numer + this%Lij*this%Mij;        this%denom = this%denom + this%Mij*this%Mij

  ! ---3 component--
  call this%filter_xyz(rho*w*T, this%Lij);            this%Lij = this%Lij - this%wFil*this%TFil/this%rhoFil
  call this%filter_xyz(Qjsgs(:,:,:,3), this%Mij);     this%Mij = this%Mij - this%QjsgsFil(:,:,:,3)
  this%numer = this%numer + this%Lij*this%Mij;        this%denom = this%denom + this%Mij*this%Mij
  this%numer = this%numer* this%Cp

  !M_ij consists of Cp. So, include Cp in numer
  do j = 1, this%nyL
      this%cmodel_local_Qjsgs(j) = max(-p_sum(sum(this%numer(:,j,:))), zero) / (p_sum(sum(this%denom(:,j,:))) + 1.0d-18)
  enddo

  this%preComputed_SFil_duFil = .false.

end subroutine

subroutine doLocalDynamicProcedure_eddy(this, rho, u, v, w, nusgs, Sij, Sii, modS_sq, qtke)
  class(sgs_cgrid), intent(inout) :: this
  real(rkind), dimension(this%nxL, this%nyL, this%nzL),    intent(in)  :: rho, u, v, w
  real(rkind), dimension(this%nxL, this%nyL, this%nzL),    intent(in)  :: nusgs
  real(rkind), dimension(this%nxL, this%nyL, this%nzL, 6), intent(in)  :: Sij
  real(rkind), dimension(this%nxL, this%nyL, this%nzL),    intent(in)  :: Sii, modS_sq, qtke

  real(rkind), dimension(this%nxL, this%nyL, this%nzL)  ::  mij_part1, mij_part2
  integer :: j

  ! space needed: rhoFil, uFil, vFil, wFil, duiFildxj (9), SFil_ij (6), 
 
  ! compute filtered velocities -- this is needed to compute \hat{\bar{S}}_ij
  call this%filter_xyz(u, this%uFil);
  call this%filter_xyz(v, this%vFil);
  call this%filter_xyz(w, this%wFil);
  call this%filter_xyz(rho, this%rhoFil);

  ! compute gradient; Sij; model kernel
  call this%gradient(this%uFil,this%duFildx,this%duFildy,this%duFildz,-this%x_bc, this%y_bc, this%z_bc)
  call this%gradient(this%vFil,this%dvFildx,this%dvFildy,this%dvFildz, this%x_bc,-this%y_bc, this%z_bc)
  call this%gradient(this%wFil,this%dwFildx,this%dwFildy,this%dwFildz, this%x_bc, this%y_bc,-this%z_bc)

  call this%get_Sij_from_duidxj(this%duiFildxj, this%SFil_ij)
  call this%get_modS_Sii(this%SFil_ij,this%modSFil_sq, this%SiiFil)
  call this%get_qtke(this%rhoFil, this%modSFil_sq, this%qtkeFil)
  call this%get_SGS_kernel(this%duiFildxj, this%SFil_ij, this%modSFil_sq, this%nusgsFil)
  this%qtkeFil   = this%deltaRatioSq * this%qtkeFil
  this%nusgsFil = this%deltaRatioSq * this%nusgsFil
 
  ! compute filtered density-velocity products \hat{rho*u} -- this is needed to compute L_ij
  call this%filter_xyz(rho*u, this%uFil);
  call this%filter_xyz(rho*v, this%vFil);
  call this%filter_xyz(rho*w, this%wFil);
  
  !!!!!!!--------- To find C_I---------------------!!!!!!!! 
  ! ---11 component--
  call this%filter_xyz(rho*u*u, this%Lij);      this%Lij = this%Lij - this%uFil*this%uFil/this%rhoFil
  call this%filter_xyz(qtke, this%Mij);         this%Mij = this%Mij - this%qtkeFil
  this%numer = this%Lij;                        this%denom = -this%Mij

  ! ---22 component--
  call this%filter_xyz(rho*v*v, this%Lij);      this%Lij = this%Lij - this%vFil*this%vFil/this%rhoFil
  this%numer = this%numer + this%Lij;  
  
  ! ---33 component--
  call this%filter_xyz(rho*w*w, this%Lij);      this%Lij = this%Lij - this%wFil*this%wFil/this%rhoFil
  this%numer = this%numer + this%Lij;  
  
  do j = 1, this%nyL
      !this%cmodel_local_tke(j) = max(p_sum(sum(this%numer(:,j,:))), zero) / (p_sum(sum(this%denom(:,j,:))) + 1.0d-18)
      this%cmodel_local_tke(j) = max(p_sum(sum(this%numer(:,j,:)))/ (p_sum(sum(this%denom(:,j,:))) + 1.0d-18),zero)
  enddo

  !!!!!!!--------- To find C-----------------------!!!!!!!! 
  ! ---11 component--
  call this%filter_xyz(rho*u*u, this%Lij);            this%Lij = this%Lij - this%uFil*this%uFil/this%rhoFil
  mij_part1 = rho*nusgs*(Sij(:,:,:,1)-Sii);           mij_part2= this%rhoFil*this%nusgsFil*(this%SFil_ij(:,:,:,1)-this%SiiFil)          
  call this%filter_xyz(mij_part1, this%Mij);          this%Mij = this%Mij - mij_part2
  this%numer = this%Lij*this%Mij;                     this%denom = this%Mij*this%Mij

  ! ---12 component--
  call this%filter_xyz(rho*u*v, this%Lij);            this%Lij = this%Lij - this%uFil*this%vFil/this%rhoFil
  mij_part1 = rho*nusgs*Sij(:,:,:,2);                 mij_part2= this%rhoFil*this%nusgsFil*this%SFil_ij(:,:,:,2)
  call this%filter_xyz(mij_part1, this%Mij);          this%Mij = this%Mij - mij_part2
  this%numer = this%numer + two*this%Lij*this%Mij;    this%denom = this%denom + two*this%Mij*this%Mij

  ! ---13 component--
  call this%filter_xyz(rho*u*w, this%Lij);            this%Lij = this%Lij - this%uFil*this%wFil/this%rhoFil
  mij_part1 = rho*nusgs*Sij(:,:,:,3);                 mij_part2= this%rhoFil*this%nusgsFil*this%SFil_ij(:,:,:,3)      
  call this%filter_xyz(mij_part1, this%Mij);          this%Mij = this%Mij - mij_part2
  this%numer = this%numer + two*this%Lij*this%Mij;    this%denom = this%denom + two*this%Mij*this%Mij

  ! ---22 component--
  call this%filter_xyz(rho*v*v, this%Lij);            this%Lij = this%Lij - this%vFil*this%vFil/this%rhoFil
  mij_part1 = rho*nusgs*(Sij(:,:,:,4)-Sii);           mij_part2= this%rhoFil*this%nusgsFil*(this%SFil_ij(:,:,:,4)-this%SiiFil)          
  call this%filter_xyz(mij_part1, this%Mij);          this%Mij = this%Mij - mij_part2
  this%numer = this%numer + this%Lij*this%Mij;        this%denom = this%denom + this%Mij*this%Mij

  ! ---23 component--
  call this%filter_xyz(rho*v*w, this%Lij);            this%Lij = this%Lij - this%vFil*this%wFil/this%rhoFil
  mij_part1 = rho*nusgs*Sij(:,:,:,5);                 mij_part2= this%rhoFil*this%nusgsFil*this%SFil_ij(:,:,:,5)        
  call this%filter_xyz(mij_part1, this%Mij);          this%Mij = this%Mij - mij_part2
  this%numer = this%numer + two*this%Lij*this%Mij;    this%denom = this%denom + two*this%Mij*this%Mij

  ! ---33 component--
  call this%filter_xyz(rho*w*w, this%Lij);            this%Lij = this%Lij - this%wFil*this%wFil/this%rhoFil
  mij_part1 = rho*nusgs*(Sij(:,:,:,6)-Sii);           mij_part2= this%rhoFil*this%nusgsFil*(this%SFil_ij(:,:,:,6)-this%SiiFil)          
  call this%filter_xyz(mij_part1, this%Mij);          this%Mij = this%Mij - mij_part2
  this%numer = this%numer + this%Lij*this%Mij;        this%denom = this%denom + this%Mij*this%Mij
  this%denom = two*this%denom;    !-ve sign ??           

  do j = 1, this%nyL
      this%cmodel_local(j) = max(p_sum(sum(this%numer(:,j,:))), zero) / (p_sum(sum(this%denom(:,j,:))) + 1.0d-18)
  enddo

  if(nrank==0) then
    do j=1,this%nyL
    write(100,'(i5,1x,2(e19.12,1x))') j, this%cmodel_local_tke(j), this%cmodel_local(j)
    end do
  endif

  !! allows reuse of several filtered fields for calcualtion of Qjsgs coefficient
  this%preComputed_SFil_duFil = .true.

end subroutine

subroutine doLocalDynamicProcedure_eddy_Qjsgs(this, rho, u, v, w, T, Qjsgs)
  class(sgs_cgrid), intent(inout) :: this
  real(rkind), dimension(this%nxL, this%nyL, this%nzL),    intent(in)  :: rho, u, v, w, T
  real(rkind), dimension(this%nxL, this%nyL, this%nzL, 3), intent(in)  :: Qjsgs

  integer :: j,k

  if(.not. this%preComputed_SFil_duFil) then
      ! compute filtered velocities -- this is needed to compute \hat{\bar{S}}_ij
      call this%filter_xyz(u, this%uFil);
      call this%filter_xyz(v, this%vFil);
      call this%filter_xyz(w, this%wFil);
      call this%filter_xyz(rho, this%rhoFil);
      ! compute gradient; Sij; tausgs model kernel
      call this%gradient(this%uFil,this%duFildx,this%duFildy,this%duFildz,-this%x_bc, this%y_bc, this%z_bc)
      call this%gradient(this%vFil,this%dvFildx,this%dvFildy,this%dvFildz, this%x_bc,-this%y_bc, this%z_bc)
      call this%gradient(this%wFil,this%dwFildx,this%dwFildy,this%dwFildz, this%x_bc, this%y_bc,-this%z_bc)
      call this%get_Sij_from_duidxj(this%duiFildxj, this%SFil_ij)
      call this%get_modS_Sii(this%SFil_ij,this%modSFil_sq, this%SiiFil)
      call this%get_SGS_kernel(this%duiFildxj, this%SFil_ij, this%modSFil_sq, this%nusgsFil)
      this%nusgsFil = this%deltaRatioSq * this%nusgsFil
  endif

  ! space needed: TFil, gradTFil (3), QjsgsFil (3), 
  ! compute filtered Temperature - this is needed to compute \hat{\gradT}_j
  call this%filter_xyz(T, this%TFil);
  call this%gradient(this%TFil,this%dTFildx,this%dTFildy,this%dTFildz, this%x_bc, this%y_bc, this%z_bc)
  do k = 1, this%nzL
     do j = 1, this%nyL
        this%nusgsFil(:,j,k) = this%cmodel_local(j) * this%nusgsFil(:,j,k)
     end do
  end do
  call this%get_Qjsgs_eddy_kernel(this%rhoFil, this%nusgsFil, this%duiFildxj, this%gradTFil, this%QjsgsFil, this%DeltaRatioSq)
  !this%QjsgsFil = this%QjsgsFil !*this%deltaRatioSq   !!nusgsFil already contains deltaRatioSq 

  ! compute filtered density-velocity products \hat{rho*u} -- this is needed to compute L_ij
  call this%filter_xyz(rho*u, this%uFil);
  call this%filter_xyz(rho*v, this%vFil);
  call this%filter_xyz(rho*w, this%wFil);
  call this%filter_xyz(rho*T, this%TFil);
 
  ! ---1 component--    
  call this%filter_xyz(rho*u*T, this%Lij);            this%Lij = this%Lij -this%uFil*this%TFil/this%rhoFil
  call this%filter_xyz(Qjsgs(:,:,:,1), this%Mij);     this%Mij = this%Mij - this%QjsgsFil(:,:,:,1)
  this%numer = this%Lij*this%Mij;                     this%denom = this%Mij*this%Mij

  ! ---2 component--
  call this%filter_xyz(rho*v*T, this%Lij);            this%Lij = this%Lij - this%vFil*this%TFil/this%rhoFil
  call this%filter_xyz(Qjsgs(:,:,:,2), this%Mij);     this%Mij = this%Mij - this%QjsgsFil(:,:,:,2)
  this%numer = this%numer + this%Lij*this%Mij;        this%denom = this%denom + this%Mij*this%Mij

  ! ---3 component--
  call this%filter_xyz(rho*w*T, this%Lij);            this%Lij = this%Lij - this%wFil*this%TFil/this%rhoFil
  call this%filter_xyz(Qjsgs(:,:,:,3), this%Mij);     this%Mij = this%Mij - this%QjsgsFil(:,:,:,3)
  this%numer = this%numer + this%Lij*this%Mij;        this%denom = this%denom + this%Mij*this%Mij
  this%numer = this%numer* this%Cp      !M_ij consists of Cp. So, include Cp in numer
  
  do j = 1, this%nyL     
      this%cmodel_local_Qjsgs(j) = max(p_sum(sum(this%numer(:,j,:))), zero) / (p_sum(sum(this%denom(:,j,:))) + 1.0d-18)
  enddo

  this%preComputed_SFil_duFil = .false.


end subroutine


subroutine filter_xyz(this, arrIn, arrOut)
  class(sgs_cgrid), target                            , intent(inout) :: this
  real(rkind), dimension(this%nxL, this%nyL, this%nzL), intent(in ) :: arrIn
  real(rkind), dimension(this%nxL, this%nyL, this%nzL), intent(out) :: arrOut
  
  real(rkind), dimension(:,:,:), pointer :: tmp1_in_x, tmp1_in_z, tmp2_in_x, tmp2_in_z
  integer :: lastx, lastz


  ! Allocate pointers for the needed buffers 
  ! 2 buffers in x and z are needed
  lastx = size(this%xbuf,4)
  tmp1_in_x => this%xbuf(:,:,:,lastx)
  tmp2_in_x => this%xbuf(:,:,:,lastx-1)

  lastz = size(this%zbuf,4)
  tmp1_in_z => this%zbuf(:,:,:,lastz)
  tmp2_in_z => this%zbuf(:,:,:,lastz-1)
 
 
  if(this%filter_in_y) then 
    ! filter in y
    call this%testfil%filtery(arrIn, arrOut, this%y_bc(1), this%y_bc(2))
  endif
  
  if(this%filter_in_x) then
    ! Then transpose to x
    call transpose_y_to_x(arrOut, tmp1_in_x, this%decomp)

    ! filter in x
    call this%testfil%filterx(tmp1_in_x, tmp2_in_x, this%x_bc(1), this%x_bc(2))

    ! Now transpose back to y
    call transpose_x_to_y(tmp2_in_x, arrOut, this%decomp)
  endif

  if(this%filter_in_z) then
    ! Now transpose to z
    call transpose_y_to_z(arrOut, tmp1_in_z, this%decomp)

    ! filter in z
    call this%testfil%filterz(tmp1_in_z, tmp2_in_z, this%z_bc(1), this%z_bc(2))

    ! Now transpose back to y
    call transpose_z_to_y(tmp2_in_z, arrOut, this%decomp)
  endif

  ! Finished
end subroutine
