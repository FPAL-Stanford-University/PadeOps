subroutine init_sigma(this, dx, dy, dz, Cs)
   class(sgs_cgrid), intent(inout) :: this
   real(rkind), intent(in) :: dx, dy, dz, Cs
   real(rkind) :: deltaLES

   ! Set the type of mnodel constant (default is constant global). 
   ! Can be reset to false via dynamic procedure initialization, 
   ! in case the dynamic procedure is planar averages
   this%useCglobal = .true. 
   
  
   if (.not. this%isPeriodic) then
      deltaLES = (1.5d0*dx*1.5d0*dy*dz)**(1.d0/3.d0)
   else
      deltaLES =  (1.5d0*dx*1.5d0*dy*1.5d0*dz)**(1.d0/3.d0)
   end if 
   this%cmodel_global = (Cs*deltaLES)**2
  
   this%isEddyViscosityModel = .true. 
   call message(1,"Sigma model initialized")
end subroutine

subroutine destroy_sigma(this)
   class(sgs_cgrid), intent(inout) :: this

   this%isEddyViscosityModel = .false. 

end subroutine

subroutine get_sigma_kernel(nu_sgs, duidxj, nxL, nyL, nzL)
   integer, intent(in) :: nxL, nyL, nzL
   real(rkind), intent(in), dimension(nxL,nyL,nzL,9):: duidxj
   real(rkind), intent(out), dimension(nxL,nyL,nzL) :: nu_sgs
   real(rkind) :: G11, G12, G13, G22, G23, G33
   real(rkind) :: I1, I2, I3, I1sq, I1cu
   real(rkind) :: alpha1, alpha2, alpha1sqrt
   real(rkind) :: alpha1tmp, alpha3
   real(rkind) :: sigma1, sigma2, sigma3, sigma1sq
   integer :: i,j,k

   do k = 1,nzL
      do j = 1,nyL 
         !$omp simd 
         do i = 1,nxL
            !print '(2(i5,1x),9(e19.12,1x))', i, j, duidxj(i,j,k,:) 
            G11=duidxj(i,j,k,1)*duidxj(i,j,k,1)+duidxj(i,j,k,4)*duidxj(i,j,k,4)+duidxj(i,j,k,7)*duidxj(i,j,k,7)
            G12=duidxj(i,j,k,1)*duidxj(i,j,k,2)+duidxj(i,j,k,4)*duidxj(i,j,k,5)+duidxj(i,j,k,7)*duidxj(i,j,k,8)
            G13=duidxj(i,j,k,1)*duidxj(i,j,k,3)+duidxj(i,j,k,4)*duidxj(i,j,k,6)+duidxj(i,j,k,7)*duidxj(i,j,k,9)
            G22=duidxj(i,j,k,2)*duidxj(i,j,k,2)+duidxj(i,j,k,5)*duidxj(i,j,k,5)+duidxj(i,j,k,8)*duidxj(i,j,k,8)
            G23=duidxj(i,j,k,2)*duidxj(i,j,k,3)+duidxj(i,j,k,5)*duidxj(i,j,k,6)+duidxj(i,j,k,8)*duidxj(i,j,k,9)
            G33=duidxj(i,j,k,3)*duidxj(i,j,k,3)+duidxj(i,j,k,6)*duidxj(i,j,k,6)+duidxj(i,j,k,9)*duidxj(i,j,k,9)
            
            I1   = G11 + G22 + G33
            I1sq = I1*I1
            I1cu = I1sq*I1

            I2 = -G11*G11 - G22*G22 - G33*G33
            I2 = I2 - two*G12*G12 - two*G13*G13
            I2 = I2 - two*G23*G23
            I2 = I2 + I1sq
            I2 = half*I2

            I3 = G11*(G22*G33 - G23*G23)
            I3 = I3 + G12*(G13*G23 - G12*G33)
            I3 = I3 + G13*(G12*G23 - G22*G13)

            alpha1 = I1sq/nine - I2/three
            alpha1 = max(alpha1,zero)

            alpha2 = I1cu/27._rkind - I1*I2/six + I3/two
            alpha1sqrt = sqrt(alpha1)
            alpha1tmp = alpha1*alpha1sqrt
            alpha1tmp = alpha2/(alpha1tmp + 1.d-13)
            alpha1tmp = min(alpha1tmp,one)
            alpha1tmp = max(alpha1tmp,-one)
            alpha1tmp = acos(alpha1tmp)
            alpha3 = (one/three)*(alpha1tmp)

            sigma1sq = I1/three + two*alpha1sqrt*cos(alpha3)
            sigma1sq = max(sigma1sq,zero)
            sigma1 = sqrt(sigma1sq)

            sigma2 = pi/three + alpha3
            sigma2 = (-two)*alpha1sqrt*cos(sigma2)
            sigma2 = sigma2 + I1/three
            sigma2 = max(sigma2,zero)
            sigma2 = sqrt(sigma2)
            
            sigma3 = pi/three - alpha3
            sigma3 = (-two)*alpha1sqrt*cos(sigma3)
            sigma3 = sigma3 + I1/three
            sigma3 = max(sigma3,zero)
            sigma3 = sqrt(sigma3)

            nu_sgs(i,j,k) = sigma3*(sigma1 - sigma2)*(sigma2 - sigma3)/(sigma1sq + 1.d-15)

         end do 
      end do 
   end do

   !print*, duidxj(5,6,3,:)
   !print*, sum(abs(nu_sgs(:,:,:)))
   !print*, sum(abs(duidxj(:,:,:,1))), sum(abs(duidxj(:,:,:,2))), sum(abs(duidxj(:,:,:,3)))
   !print*, sum(abs(duidxj(:,:,:,4))), sum(abs(duidxj(:,:,:,5))), sum(abs(duidxj(:,:,:,6)))
   !print*, sum(abs(duidxj(:,:,:,7))), sum(abs(duidxj(:,:,:,8))), sum(abs(duidxj(:,:,:,9)))

end subroutine
