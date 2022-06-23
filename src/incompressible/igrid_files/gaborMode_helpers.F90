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


subroutine HaloUpdateVelocities(this, uh, vh, wh, duidxj_h)
    use decomp_2d, only: update_halo
    class(igrid), intent(inout) :: this 
    real(rkind), dimension(:,:,:), allocatable, intent(out) :: uh, vh, wh
    real(rkind), dimension(:,:,:,:), allocatable, intent(out) :: duidxj_h
    real(rkind), dimension(:,:,:), allocatable :: buff
    integer :: idx

    if (allocated(uh)) deallocate(uh) 
    if (allocated(vh)) deallocate(vh) 
    if (allocated(wh)) deallocate(wh) 
    if (allocated(duidxj_h)) deallocate(duidxj_h)

    call update_halo(this%u, uh, 1, this%gpC)
    call update_halo(this%v, vh, 1, this%gpC)
    call update_halo(this%w, wh, 1, this%gpC)

    allocate(duidxj_h(size(uh,1),size(uh,2),size(uh,3),9))
    do idx = 1,9
      call update_halo(this%duidxjC(:,:,:,idx), buff, 1, this%gpC)
      duidxj_h(:,:,:,idx) = buff
    end do
    deallocate(buff)
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
