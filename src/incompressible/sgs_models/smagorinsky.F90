subroutine init_smagorinsky(this,dx,dy,dz,Cs,ncWall,z0,useWallDamping, zMeshC, zMeshE)
   class(sgs_igrid), intent(inout) :: this
   logical, intent(in) :: useWallDamping 
   real(rkind), intent(in) :: z0, ncWall, Cs, dx, dy, dz
   real(rkind), intent(in), dimension(:) :: zMeshC, zMeshE
   real(rkind) :: deltaLES
   real(rkind), parameter :: kappa = 0.4d0

   deltaLES = (1.5d0*dx*1.5d0*dy*dz)**(1.d0/3.d0)
   
   allocate(this%cmodelC(size(zMeshC)))
   allocate(this%cmodelE(size(zMeshE)))
   if (useWallDamping) then
      this%cmodelC = ( Cs**(-real(ncWall,rkind)) + (kappa*(zMeshC/deltaLES + &
          & z0/deltaLES))**(-real(ncWall,rkind))  )**(-one/real(ncWall,rkind))
      this%cmodelE = ( Cs**(-real(ncWall,rkind)) + (kappa*(zMeshE/deltaLES + &
          & z0/deltaLES))**(-real(ncWall,rkind))  )**(-one/real(ncWall,rkind))
      this%cmodelC = (deltaLES*this%cmodelC)**2    
      this%cmodelE = (deltaLES*this%cmodelE)**2    
   else
      this%cmodelC = (Cs*deltaLES)**2
      this%cmodelE = (Cs*deltaLES)**2
   end if
   this%isEddyViscosityModel = .true. 

   call message(1,"Smagorinsky model initialized")
end subroutine

pure subroutine get_nu_smagorinsky(Sij, nuSGS, nxL, nyL, nzL)
   integer, intent(in) :: nxL, nyL, nzL
   real(rkind), dimension(nxL,nyL,nzL), intent(out) :: nuSGS
   real(rkind), dimension(nxL,nyL,nzL,6), intent(in)  :: Sij
   real(rkind), dimension(nxL, nyL) :: S
   integer :: k

   
   do k = 1,nzL
      ! S12 
      S = Sij(:,:,k,1)*Sij(:,:,k,1) ! S11*S11
      S = S + 2.d0*(Sij(:,:,k,2)*Sij(:,:,k,2)) ! S12*S12 + S21*S21
      S = S + 2.d0*(Sij(:,:,k,3)*Sij(:,:,k,3)) ! S13*S13 + S31*S31
      S = S + (Sij(:,:,k,4)*Sij(:,:,k,4)) ! S22*S22
      S = S + 2.d0*(Sij(:,:,k,5)*Sij(:,:,k,5)) ! S23*S23 + S32*S32
      S = S + (Sij(:,:,k,6)*Sij(:,:,k,6)) ! S33*S33

      ! Now do modS = sqrt(2* S_ij*S_ij)
      S = 2.d0*S
      nuSGS(:,:,k) = sqrt(S)
   end do


end subroutine


subroutine destroy_smagorinsky(this)
   class(sgs_igrid), intent(inout) :: this
   deallocate(this%cmodelC, this%cmodelE)
   this%isEddyViscosityModel = .false. 

end subroutine
