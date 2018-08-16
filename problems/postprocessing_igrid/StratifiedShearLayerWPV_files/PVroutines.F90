module PVroutines
  use kind_parameters, only: rkind, clen
  use igrid_Operators, only: igrid_ops
  implicit none
   contains


  subroutine ComputeF(ops, ufluct, vfluct, w, T, Fx, Fy, Fz, buff1, buff2, buff3, buff4)
     class(igrid_ops), intent(inout) :: ops
     real(rkind), dimension(ops%gp%xsz(1),ops%gp%xsz(2),ops%gp%xsz(3)), intent(in)  :: ufluct, vfluct, w, T
     real(rkind), dimension(ops%gp%xsz(1),ops%gp%xsz(2),ops%gp%xsz(3)), intent(out) :: Fx, Fy, Fz 
     real(rkind), dimension(ops%gp%xsz(1),ops%gp%xsz(2),ops%gp%xsz(3)), intent(inout) :: buff1, buff2, buff3, buff4

      
     call ops%getCurl(ufluct,vfluct,w, buff1,buff2,buff3,1,1,1,1) 
     call ops%getGradient(T, Fx,Fy,Fz,1,1)
     buff4 = Fx*Fx + Fy*Fy + Fz*Fz
     buff4 = buff4**(1.d0/2.d0) ! |gradT|
     Fx = Fx / buff4
     Fy = Fy / buff4
     Fz = Fz / buff4 ! gradT/|gradT|
  
     buff4 = buff1 * Fx
     buff4 = buff4 + buff2 * Fy
     buff4 = buff4 + buff3 * Fz ! omega \cdot gradT/|gradT|
  
     buff1 = -Fx * buff4
     buff2 = -Fy * buff4
     buff3 = -Fz * buff4
  
     call ops%dealias(buff1)
     call ops%dealias(buff2)
     call ops%dealias(buff3)
     call ops%getCurl(buff1,buff2,buff3, Fx,Fy,Fz,0,0,0,0)
  
  
  end subroutine

end module



