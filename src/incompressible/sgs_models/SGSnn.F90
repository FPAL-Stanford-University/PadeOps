subroutine init_SGSnn(this)
   use fortran_assert, only: assert
   class(sgs_igrid), intent(inout) :: this

   this%isEddyViscosityModel = .false.

   call message(1,"NN model initialized")
end subroutine

subroutine getSijRijforNNmod(duidxj, Sij, Rij, nxL, nyL, nzL)
   integer, intent(in) :: nxL, nyL, nzL
   real(rkind), dimension(nxL,nyL,nzL,9), intent(in) :: duidxj
   real(rkind), dimension(nxL,nyL,nzL,9), intent(out) :: Sij, Rij
   integer :: i, j, k

   do k = 1,nzL
      do j = 1,nyL
         !$omp simd
         do i = 1,nxL
            Sij(i,j,k,1) = duidxj(i,j,k,1) ! S11 = dudx
            Sij(i,j,k,2) = 0.5d0*(duidxj(i,j,k,2) + duidxj(i,j,k,4)) ! S12 = 0.5*(dudy + dvdx)
            Sij(i,j,k,3) = 0.5d0*(duidxj(i,j,k,3) + duidxj(i,j,k,7)) ! S13 = 0.5*(dudz + dwdx)
           
            Sij(i,j,k,5) = duidxj(i,j,k,5) ! S22 = dvdy
            Sij(i,j,k,6) = 0.5d0*(duidxj(i,j,k,6) + duidxj(i,j,k,8)) ! S23 = 0.5*(dvdz + dwdy)
            
            Sij(i,j,k,9) = duidxj(i,j,k,9) ! S33 = dwdz

            Rij(i,j,k,2) = 0.5d0*(duidxj(i,j,k,2) - duidxj(i,j,k,4)) ! R12 = 0.5*(dudy - dvdx)
            Rij(i,j,k,3) = 0.5d0*(duidxj(i,j,k,3) - duidxj(i,j,k,7)) ! R13 = 0.5*(dudz - dwdx)
            Rij(i,j,k,6) = 0.5d0*(duidxj(i,j,k,6) - duidxj(i,j,k,8)) ! R23 = 0.5*(dvdz - dwdy)
         end do 
      end do 
   end do 
   Sij(:,:,:,4) = Sij(:,:,:,2) ! S21 = S12
   Sij(:,:,:,7) = Sij(:,:,:,3) ! S31 = S13
   Sij(:,:,:,8) = Sij(:,:,:,6) ! S32 = S23

   Rij(:,:,:,1) = 0.d0
   Rij(:,:,:,5) = 0.d0
   Rij(:,:,:,9) = 0.d0

   Rij(:,:,:,4) = -Rij(:,:,:,2) ! R21 = -R12
   Rij(:,:,:,7) = -Rij(:,:,:,3) ! R31 = -R13
   Rij(:,:,:,8) = -Rij(:,:,:,6) ! R32 = -R23
end subroutine

subroutine compute_tauij_NN(this)
  class(sgs_igrid), intent(inout) :: this

  ! TODO: Andy, compute your invariants, call your NN, compute tauij
end subroutine
