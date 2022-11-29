subroutine initLargeScales(this, tid, rid)
    class(igrid), intent(inout) :: this
    integer, intent(in) :: tid, rid

    ! Fill u, v and w 
    call this%readLargeScales(tid, rid)
    !call this%readRestartFile(tid, rid)

    ! Compute uhat, vhat and what 
    call this%spectC%fft(this%u,this%uhat)   
    call this%spectC%fft(this%v,this%vhat)   
    call this%spectE%fft(this%w,this%what) 

    ! Interpolate primitive variables 
    call this%interp_PrimitiveVars()

    ! Compute gradients 
    call this%compute_duidxj()

    ! Now u, v, wC and duidxjC are all accessible
end subroutine


subroutine HaloUpdateField(this, u, uh)
    use procgrid_mod, only: num_pad 
    class(igrid), intent(inout) :: this 
    real(rkind), dimension(1:this%gpC%xsz(1), 1:this%gpC%xsz(2), 1:this%gpC%xsz(3)), intent(in) :: u 
    real(rkind), dimension(-num_pad+1:, -num_pad+1:, -num_pad+1:), intent(out) :: uh 

    uh(1:this%gpC%xsz(1), 1:this%gpC%xsz(2), 1:this%gpC%xsz(3)) = u

    call this%pg%halo_exchange(uh)
end subroutine 

subroutine HaloUpdateVelocities(this, uh, vh, wh, duidxj_h)
    use procgrid_mod, only: num_pad 
    class(igrid), intent(inout) :: this 
    real(rkind), dimension(-num_pad+1:, -num_pad+1:, -num_pad+1:), intent(out) :: uh, vh, wh
    real(rkind), dimension(-num_pad+1:, -num_pad+1:, -num_pad+1:,:), intent(out) :: duidxj_h
    integer :: idx

    call this%HaloUpdateField(this%u, uh)
    call this%HaloUpdateField(this%v, vh)
    call this%HaloUpdateField(this%w, wh)

    do idx = 1,size(duidxj_h,4)
        call this%HaloUpdateField(this%duidxjC(:,:,:,idx), duidxj_h(:,:,:,idx))
    end do 

end subroutine


subroutine ProjectToFixBC(this)
    class(igrid), intent(inout) :: this 

    ! Assume we have u, v and wC 

    ! Get w from wC
    call this%spectC%fft(this%u ,this%uhat )   
    call this%spectC%fft(this%v ,this%vhat )   
    call this%spectC%fft(this%wC,this%whatC)   

    call transpose_y_to_z(this%whatC, this%cbuffzC(:,:,:,1),this%sp_gpC)
    call this%Pade6opZ%interpz_C2E(this%cbuffzC(:,:,:,1),this%cbuffzE(:,:,:,1),0,0)
    call transpose_z_to_y(this%cbuffzE(:,:,:,1),this%what,this%sp_gpE)

    call this%padepoiss%PressureProjection(this%uhat,this%vhat,this%what)

    call transpose_y_to_z(this%what, this%cbuffzE(:,:,:,1),this%sp_gpE)
    call this%Pade6opZ%interpz_E2C(this%cbuffzE(:,:,:,1),this%cbuffzC(:,:,:,1),0,0)
    call transpose_z_to_y(this%cbuffzC(:,:,:,1), this%whatC,this%sp_gpC)

    call this%spectC%ifft(this%uhat ,this%u )   
    call this%spectC%ifft(this%vhat ,this%v )   
    call this%spectC%ifft(this%whatC,this%wC)   

end subroutine

subroutine getPressure(this)
    class(igrid), intent(inout) :: this

end subroutine 


subroutine readLargeScales(this,tid, rid)
    use decomp_2d_io
    class(igrid), intent(inout) :: this
    integer, intent(in) :: rid, tid
    character(len=clen) :: tempname, fname
    
   call readField3D(rid, tid, this%inputDir, "uVel", this%u , this%gpC)
   call readField3D(rid, tid, this%inputDir, "vVel", this%v , this%gpC)
   call readField3D(rid, tid, this%inputDir, "wVel", this%wC, this%gpC)

   call this%spectC%fft(this%wC,this%whatC)   
   call transpose_y_to_z(this%whatC, this%cbuffzC(:,:,:,1),this%sp_gpC)
   call this%Pade6opZ%interpz_C2E(this%cbuffzC(:,:,:,1),this%cbuffzE(:,:,:,1),0,0)
   call transpose_z_to_y(this%cbuffzE(:,:,:,1),this%what,this%sp_gpE)
   call this%spectE%ifft(this%what, this%w)
end subroutine

subroutine fixGradientsForPeriodicity(this,periodicBCs)
  class(igrid), intent(inout) :: this
  logical, dimension(3), intent(in) :: periodicBCs

  if ((.not. periodicBCs(1)) .or. (.not. periodicBCs(2))) then
    if (.not. allocated(this%cd6opX)) then
      allocate(this%cd6opX)
      ierr = this%cd6opX%init(this%nx, this%dx, periodicBCs(1), 0, 0)
    end if
    if (.not. allocated(this%cd6opY)) then
      allocate(this%cd6opY)
      ierr = this%cd6opY%init(this%ny, this%dy, periodicBCs(2), 0, 0)
    end if
    if (.not. allocated(this%cd6opZ)) then
      allocate(this%cd6opZ)
      ierr = this%cd6opZ%init(this%nz, this%dz, periodicBCs(3), 0, 0)
    end if

    ! x - derivatives
  
    call this%cd6opX%dd1(this%u,this%duidxjC(:,:,:,1),this%gpC%xsz(2),this%gpC%xsz(3))
    call this%dfdyC2CinX(this%u,this%duidxjC(:,:,:,2))
    call this%dfdzC2CinX(this%u,this%duidxjC(:,:,:,3))
    
    call this%cd6opX%dd1(this%v,this%duidxjC(:,:,:,4),this%gpC%xsz(2),this%gpC%xsz(3))
    call this%dfdyC2CinX(this%v,this%duidxjC(:,:,:,5))
    call this%dfdzC2CinX(this%v,this%duidxjC(:,:,:,6))
    
    call this%cd6opX%dd1(this%w,this%duidxjC(:,:,:,7),this%gpC%xsz(2),this%gpC%xsz(3))
    call this%dfdyC2CinX(this%w,this%duidxjC(:,:,:,8))
    call this%dfdzC2CinX(this%w,this%duidxjC(:,:,:,9))
    
  end if


end subroutine

subroutine dfdyC2CinX(this, f, dfdy)
  class(igrid), intent(inout) :: this
  real(rkind), dimension(:,:,:), intent(in) :: f
  real(rkind), dimension(:,:,:), intent(out) :: dfdy

  call transpose_x_to_y(f, this%rbuffyC(:,:,:,1), this%gpC)
  call this%cd6opY%dd2(this%rbuffyC(:,:,:,1), this%rbuffyC(:,:,:,2),&
    this%gpC%ysz(1),this%gpC%ysz(3))
  call transpose_y_to_x(this%rbuffyC(:,:,:,2), dfdy, this%gpC)

end subroutine

subroutine dfdzC2CinX(this, f, dfdz)
  class(igrid), intent(inout) :: this
  real(rkind), dimension(:,:,:), intent(in) :: f
  real(rkind), dimension(:,:,:), intent(out) :: dfdz

  call transpose_x_to_y(f, this%rbuffyC(:,:,:,1), this%gpC)
  call transpose_y_to_z(this%rbuffyC(:,:,:,1), this%rbuffzC(:,:,:,1), this%gpC)
  call this%cd6opZ%dd3(this%rbuffzC(:,:,:,1), this%rbuffzC(:,:,:,2),&
    this%gpC%zsz(1),this%gpC%zsz(2))
  call transpose_z_to_y(this%rbuffzC(:,:,:,2), this%rbuffyC(:,:,:,2), this%gpC) 
  call transpose_y_to_x(this%rbuffyC(:,:,:,2), dfdz, this%gpC)

end subroutine
