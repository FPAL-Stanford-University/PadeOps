subroutine get_Sij_from_duidxj(duidxj, Sij, nxL, nyL, nzL)
   integer, intent(in) :: nxL, nyL, nzL
   real(rkind), dimension(nxL,nyL,nzL,9), intent(in) :: duidxj
   real(rkind), dimension(nxL,nyL,nzL,6), intent(out) :: Sij
   integer :: i, j, k

   do k = 1,nzL
      do j = 1,nyL
         !$omp simd
         do i = 1,nxL
            Sij(i,j,k,1) = duidxj(i,j,k,1) ! S11 = dudx
            Sij(i,j,k,4) = duidxj(i,j,k,5) ! S22 = dvdy
            Sij(i,j,k,6) = duidxj(i,j,k,9) ! S33 = dwdz

            Sij(i,j,k,2) = 0.5d0*(duidxj(i,j,k,2) + duidxj(i,j,k,4)) ! S12 = 0.5*(dudy + dvdx)
            Sij(i,j,k,3) = 0.5d0*(duidxj(i,j,k,3) + duidxj(i,j,k,7)) ! S13 = 0.5*(dudz + dwdx)
            Sij(i,j,k,5) = 0.5d0*(duidxj(i,j,k,6) + duidxj(i,j,k,8)) ! S23 = 0.5*(dvdz + dwdy)
         end do 
      end do 
   end do 
end subroutine

subroutine allocateMemory_EddyViscosity(this)
  class(sgs_igrid), intent(inout) :: this
  
  allocate(this%nu_sgs_C(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
  allocate(this%nu_sgs_E(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)))
  allocate(this%S_ij_C(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),6))
  allocate(this%S_ij_E(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),6))

  if (this%isStratified .or. this%initspinup) then
      allocate(this%kappa_sgs_C(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)))
      allocate(this%kappa_sgs_E(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)))
  end if
end subroutine

subroutine get_SGS_kernel(this,duidxjC, duidxjE)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3),9), intent(in) :: duidxjC
   real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3),9), intent(in) :: duidxjE

   select case(this%mid) 
   case (0)
      ! Smagorinsky
      call get_smagorinsky_kernel(this%S_ij_C,this%nu_sgs_C, &
                                 this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3))
      
      if (this%explicitCalcEdgeEddyViscosity) then
         call get_smagorinsky_kernel(this%S_ij_E,this%nu_sgs_E, &
                              this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3))
      else
         call this%interpolate_eddy_viscosity()                     
      end if
   case (1)
      ! Sigma
      call get_sigma_kernel(this%nu_sgs_C, duidxjC, this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3))
      
      if (this%explicitCalcEdgeEddyViscosity) then
         call get_sigma_kernel(this%nu_sgs_E, duidxjE, this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3))
      else
         call this%interpolate_eddy_viscosity()                     
      end if
   case (2)
      ! AMD 
      call get_amd_kernel(this%nu_sgs_C, this%camd_x, this%camd_y, this%camd_z, duidxjC, this%S_ij_C, &
                                 this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3))
      
      if (this%explicitCalcEdgeEddyViscosity) then
         call get_amd_kernel(this%nu_sgs_E, this%camd_x, this%camd_y, this%camd_z, duidxjE, this%S_ij_E, &
                                 this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3))
      else
         call this%interpolate_eddy_viscosity()                     
      end if
   end select

end subroutine

subroutine multiply_by_model_constant(this)
  class(sgs_igrid), intent(inout) :: this 
  integer :: k

  if (this%useCglobal) then
      this%nu_sgs_C = this%cmodel_global*this%nu_sgs_C
      this%nu_sgs_E = this%cmodel_global*this%nu_sgs_E
  else
      do k = 1,size(this%nu_sgs_C,3)
         this%nu_sgs_C(:,:,k) = this%cmodelC(k)*this%nu_sgs_C(:,:,k)
      end do 
      do k = 1,size(this%nu_sgs_E,3)
         this%nu_sgs_E(:,:,k) = this%cmodelE(k)*this%nu_sgs_E(:,:,k)
      end do 
  end if


end subroutine

subroutine interpolate_eddy_viscosity(this)
  class(sgs_igrid), intent(inout) :: this 

  call transpose_x_to_y(this%nu_sgs_C,this%rbuffyC(:,:,:,1), this%gpC)
  call transpose_y_to_z(this%rbuffyC(:,:,:,1), this%rbuffzC(:,:,:,1), this%gpC)
  call this%PadeDer%interpz_C2E(this%rbuffzC(:,:,:,1), this%rbuffzE(:,:,:,1),0,0)
  call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
  call transpose_y_to_x(this%rbuffyE(:,:,:,1), this%nu_sgs_E, this%gpE)

  where (this%nu_sgs_E < 0.d0) 
   this%nu_sgs_E = 0.d0
  end where
end subroutine

subroutine destroyMemory_EddyViscosity(this)
  class(sgs_igrid), intent(inout) :: this
 
  deallocate(this%S_ij_C, this%S_ij_E)
  deallocate(this%nu_sgs_C, this%nu_sgs_E)
end subroutine
