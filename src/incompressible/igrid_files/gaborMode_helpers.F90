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
