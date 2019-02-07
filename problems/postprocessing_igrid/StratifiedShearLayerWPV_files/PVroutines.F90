module PVroutines
  use kind_parameters, only: rkind, clen
  use igrid_Operators, only: igrid_ops
  implicit none
   contains


  subroutine ComputeF(ops, wantPi, ufluct, vfluct, w, T, Fx, Fy, Fz, buff1, buff2, buff3, buff4)
     class(igrid_ops), intent(inout) :: ops
     real(rkind), dimension(ops%gp%xsz(1),ops%gp%xsz(2),ops%gp%xsz(3)), intent(in)  :: ufluct, vfluct, w, T
     real(rkind), dimension(ops%gp%xsz(1),ops%gp%xsz(2),ops%gp%xsz(3)), intent(out) :: Fx, Fy, Fz 
     real(rkind), dimension(ops%gp%xsz(1),ops%gp%xsz(2),ops%gp%xsz(3)), intent(inout) :: buff1, buff2, buff3, buff4
     logical, intent(in) :: wantPi

     !call ops%initfilter(ops%gp%xsz(1)/4,ops%gp%ysz(2)/4,4)
     !call ops%filterfield(T)

     call ops%getCurl(ufluct,vfluct,w, buff1,buff2,buff3,0,0,0,0) 
     call ops%getGradient(T, Fx,Fy,Fz,0,0)
     
     buff4 = Fx*Fx + Fy*Fy + Fz*Fz
     buff4 = sqrt(buff4) + 1.d-14 ! |gradT|
     !call ops%dealias(buff4)
     
     Fx = Fx / buff4
     Fy = Fy / buff4
     Fz = Fz / buff4  ! gradT/|gradT|
     !call ops%dealias(Fx)
     !call ops%dealias(Fy)
     !call ops%dealias(Fz)
     !call ops%filterfield_inplace(Fx)
     !call ops%filterfield_inplace(Fy)
     !call ops%filterfield_inplace(Fz)
  

     buff4 = buff1 * Fx
     buff4 = buff4 + buff2 * Fy
     buff4 = buff4 + buff3 * Fz ! omega \cdot gradT/|gradT|
     !call ops%dealias(buff4)

     if (wantPi) then ! Calculate -omega_Pi Components
       buff1 = -Fx * buff4
       buff2 = -Fy * buff4
       buff3 = -Fz * buff4
     else           ! Solving rather for -omega_Xi
       buff1 = -buff1 + Fx * buff4
       buff2 = -buff2 + Fy * buff4
       buff3 = -buff3 + Fz * buff4
     end if
     !call ops%softdealias(buff1)
     !call ops%softdealias(buff2)
     !call ops%softdealias(buff3)

     call ops%getCurl(buff1,buff2,buff3, Fx,Fy,Fz,0,0,0,0)

     !call ops%getFluct_from_MeanZ(buff1,Fx)
     !call ops%getFluct_from_MeanZ(buff2,Fy)
     !call ops%getFluct_from_MeanZ(buff3,Fz)

     buff1 = -buff1 ! To use the omega_Pi/Xi for later!
     buff2 = -buff2
     buff3 = -buff3
  
  end subroutine

end module



