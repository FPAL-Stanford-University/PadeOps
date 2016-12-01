subroutine get_AMD_Op(this, nuSGS, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, S11, S12, S13, S22, S23, S33)
    class(sgs), intent(inout), target :: this
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: dudx, dudy, dudz 
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: dvdx, dvdy, dvdz 
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: dwdx, dwdy, dwdz 
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: S11, S12, S13 
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: S22, S23, S33
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(out) :: nuSGS

     this%AMD_den = dudx*dudx + dvdx*dvdx + dwdx*dwdx + dudy*dudy + dvdy*dvdy + &
                 & dwdy*dwdy  + dudz*dudz + dvdz*dvdz + dwdz*dwdz  
   
     ! k = 1            
     this%AMD_num = (this%cx**2)*(dudx*dudx*S11 + dvdx*dvdx*S22 + dwdx*dwdx*S33 + 2.d0*dudx*dvdx*S12 + 2.d0*dudx*dwdx*S13 + 2.d0*dvdx*dwdx*S23)
     
     ! k = 2
     this%AMD_num = this%AMD_num + (this%cy**2)*(dudy*dudy*S11 + dvdy*dvdy*S22 + dwdy*dwdy*S33 + 2.d0*dudy*dvdy*S12 + 2.d0*dudy*dwdy*S13 + 2.d0*dvdy*dwdy*S23)

     ! k = 3
     this%AMD_num = this%AMD_num + (this%cz**2)*(dudz*dudz*S11 + dvdz*dvdz*S22 + dwdz*dwdz*S33 + 2.d0*dudz*dvdz*S12 + 2.d0*dudz*dwdz*S13 + 2.d0*dvdz*dwdz*S23)
    
     nuSGS = this%AMD_num/this%AMD_den
     nuSGS = -nuSGS
     nuSGS = max(zero,nuSGS)

end subroutine

subroutine get_AMD_Op_strat(this, nuSGS, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, S11, S12, S13, S22, S23, S33, dTdx, dTdy)
    class(sgs), intent(inout), target :: this
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: dudx, dudy, dudz 
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: dvdx, dvdy, dvdz 
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: dwdx, dwdy, dwdz 
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: S11, S12, S13 
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: S22, S23, S33
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: dTdx, dTdy
    real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(out) :: nuSGS

     this%AMD_den = dudx*dudx + dvdx*dvdx + dwdx*dwdx + dudy*dudy + dvdy*dvdy + &
                 & dwdy*dwdy  + dudz*dudz + dvdz*dvdz + dwdz*dwdz  
   
     ! k = 1            
     this%AMD_num = (this%cx**2)*(dudx*dudx*S11 + dvdx*dvdx*S22 + dwdx*dwdx*S33 + 2.d0*dudx*dvdx*S12 + 2.d0*dudx*dwdx*S13 + 2.d0*dvdx*dwdx*S23)
     
     ! k = 2
     this%AMD_num = this%AMD_num + (this%cy**2)*(dudy*dudy*S11 + dvdy*dvdy*S22 + dwdy*dwdy*S33 + 2.d0*dudy*dvdy*S12 + 2.d0*dudy*dwdy*S13 + 2.d0*dvdy*dwdy*S23)

     ! k = 3
     this%AMD_num = this%AMD_num + (this%cz**2)*(dudz*dudz*S11 + dvdz*dvdz*S22 + dwdz*dwdz*S33 + 2.d0*dudz*dvdz*S12 + 2.d0*dudz*dwdz*S13 + 2.d0*dvdz*dwdz*S23)
   
     ! Buoyancy term
     this%AMD_num = this%AMD_num + (1.d0/(this%Theta0*this%Fr*this%Fr))*(dwdx*dTdx + dwdy*dTdy  + dwdz*this%dTdzC_diff)

     nuSGS = this%AMD_num/this%AMD_den
     nuSGS = -nuSGS
     nuSGS = max(zero,nuSGS)

end subroutine
