subroutine init_smagorinsky(this,dx,dy,dz,Cs,ncWall,z0,useWallDamping, zMeshC, zMeshE)
   class(sgs_igrid), intent(inout) :: this
   logical, intent(in) :: useWallDamping 
   real(rkind), intent(in) :: z0, ncWall, Cs, dx, dy, dz
   real(rkind), intent(in), dimension(:) :: zMeshC, zMeshE
   real(rkind) :: deltaLES
   real(rkind), parameter :: kappa = 0.4d0

   ! Set the type of mnodel constant (default is wall function). 
   ! Can be reset to true via dynamic procedure initialization, 
   ! in case the global dynamic procedure is used. 
   
   
   if (.not. this%isPeriodic) then
      deltaLES = (1.5d0*dx*1.5d0*dy*dz)**(1.d0/3.d0)
   else
      deltaLES =  (1.5d0*dx*1.5d0*dy*1.5d0*dz)**(1.d0/3.d0)
   end if 
   
   if (useWallDamping) then
      this%useCglobal = .false. 
      this%cmodelC = ( Cs**(-real(ncWall,rkind)) + (kappa*(zMeshC/deltaLES + &
          & z0/deltaLES))**(-real(ncWall,rkind))  )**(-one/real(ncWall,rkind))
      this%cmodelE = ( Cs**(-real(ncWall,rkind)) + (kappa*(zMeshE/deltaLES + &
          & z0/deltaLES))**(-real(ncWall,rkind))  )**(-one/real(ncWall,rkind))
      this%cmodelC = (deltaLES*this%cmodelC)**2    
      this%cmodelE = (deltaLES*this%cmodelE)**2    
   else
      this%useCglobal = .true. 
      this%cmodelC = (Cs*deltaLES)**2
      this%cmodelE = (Cs*deltaLES)**2
      this%cmodel_global = (Cs*deltaLES)**2
   end if
   this%isEddyViscosityModel = .true. 

   call message(1,"Smagorinsky model initialized")
end subroutine

subroutine get_smagorinsky_kernel(Sij, nuSGS, nxL, nyL, nzL)
   integer, intent(in) :: nxL, nyL, nzL
   real(rkind), dimension(nxL,nyL,nzL), intent(out) :: nuSGS
   real(rkind), dimension(nxL,nyL,nzL,6), intent(in)  :: Sij
   !real(rkind), dimension(nxL, nyL) :: S
   real(rkind) :: S
   integer :: i,j,k

   
   do k = 1,nzL
      do j = 1,nyL
         !$omp simd 
         do i = 1,nxL
            S = Sij(i,j,k,1)*Sij(i,j,k,1) ! S11*S11
            S = S + 2.d0*(Sij(i,j,k,2)*Sij(i,j,k,2)) ! S12*S12 + S21*S21
            S = S + 2.d0*(Sij(i,j,k,3)*Sij(i,j,k,3)) ! S13*S13 + S31*S31
            S = S + (Sij(i,j,k,4)*Sij(i,j,k,4)) ! S22*S22
            S = S + 2.d0*(Sij(i,j,k,5)*Sij(i,j,k,5)) ! S23*S23 + S32*S32
            S = S + (Sij(i,j,k,6)*Sij(i,j,k,6)) ! S33*S33

            ! Now do modS = sqrt(2* S_ij*S_ij)
            S = 2.d0*S
            nuSGS(i,j,k) = sqrt(S)
         end do 
      end do 
   end do


end subroutine


subroutine destroy_smagorinsky(this)
   class(sgs_igrid), intent(inout) :: this
   deallocate(this%cmodelC, this%cmodelE)
   this%isEddyViscosityModel = .false. 

end subroutine
