subroutine compute_Tscale(this, u, v, w) 
   class(sgs_igrid), intent(inout), target :: this
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: u, v, w 
   real(rkind), dimension(:,:,:), pointer :: rb1, rb2
  
   rb1 => this%rbuffxC(:,:,:,1)
   rb2 => this%rbuffxC(:,:,:,2)
   rb1 = (one/this%dx)*u
   rb2 = (one/this%dy)*v 
   rb1 = abs(rb1) 
   rb2 = abs(rb2)
   rb1 = rb1 + rb2
   rb2 = (one/this%dz)*w
   rb2 = abs(rb2)
   rb1 = rb1 + rb2
   this%Tscale = one/p_maxval(rb1)

end subroutine 


subroutine compute_kappa_bounding(this, T, dTdx, dTdy, dTdz, u, v, w)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: T, dTdx, dTdy, dTdz, u, v, w
   integer :: i, j, k
   real(rkind) :: num, den, delta, lowbound, highbound, e


   select case (this%kappa_bounding_scheme)
   case (Cook04)
       do k = 1,size(dTdx,3)
          do j = 1,size(dTdx,2)
             !$omp simd
             do i = 1,size(dTdx,1)
                den = sqrt(dTdx(i,j,k)*dTdx(i,j,k) + dTdy(i,j,k)*dTdy(i,j,k) +  dTdz(i,j,k)*dTdz(i,j,k))
                if (den < 1.d-18) then 
                   delta = sqrt(this%dx*this%dx + this%dy*this%dy + this%dz*this%dz)
                else
                   num = this%dx*abs(dTdx(i,j,k)) + this%dy*abs(dTdy(i,j,k)) + this%dz*abs(dTdz(i,j,k))
                   delta = num/(den + 1.d-20) 
                end if 
                !this%eta_boundingC(i,j,k)   = abs(T(i,j,k) - this%lowbound) + this%lowbound - this%highbound + abs(this%highbound - T(i,j,k))
                !this%kappa_boundingC(i,j,k) = this%Cy*(delta*delta/this%Tscale)*(this%eta_boundingC(i,j,k))
                this%kappa_boundingC(i,j,k) = this%Cy*(delta*delta/this%Tscale)*(&
                  abs(T(i,j,k) - this%lowbound) + this%lowbound - this%highbound + abs(this%highbound - T(i,j,k)))
             end do
          end do 
       end do 
       !call this%filter_and_interpolate_kappa_bounding(this%eta_boundingC  ,this%eta_boundingE  )
       call this%filter_and_interpolate_kappa_bounding(this%kappa_boundingC,this%kappa_boundingE)
   case (CookMod1)
       lowbound  = this%lowbound  + this%shrink_scalar_bounds_fact 
       highbound = this%highbound - this%shrink_scalar_bounds_fact 
       this%eta_boundingC   = abs(T - lowbound) + lowbound - highbound + abs(highbound - T )
       call this%filter_and_interpolate_kappa_bounding(this%eta_boundingC,this%eta_boundingE)
       this%kappa_boundingC = this%Cy*(abs(u)*this%dx + abs(v)*this%dy + abs(w)*this%dz)*this%eta_boundingC
       call this%interpolate_C2E(this%kappa_boundingC,this%kappa_boundingE)
   case (CookMod2)
       do k = 1,size(dTdx,3)
           do j = 1,size(dTdx,2)
               !$omp simd
               do i = 1,size(dTdx,1)
                   e = min(this%shrink_scalar_bounds_fact,this%dx*abs(dTdx(i,j,k))+this%dy*abs(dTdy(i,j,k))+this%dz*abs(dTdz(i,j,k)))
                   lowbound  = this%lowbound  + e 
                   highbound = this%highbound - e 
                   this%eta_boundingC(i,j,k)   = abs(T(i,j,k) - lowbound) + lowbound - highbound + abs(highbound - T(i,j,k))
                   this%kappa_boundingC(i,j,k) = this%Cy*(abs(u(i,j,k))*this%dx + abs(v(i,j,k))*this%dy + abs(w(i,j,k))*this%dz)
               end do
           end do
       end do
       call this%filter_and_interpolate_kappa_bounding(this%eta_boundingC,this%eta_boundingE)
       this%kappa_boundingC = this%kappa_boundingC*this%eta_boundingC
       call this%interpolate_C2E(this%kappa_boundingC,this%kappa_boundingE)
   end select
end subroutine

subroutine interpolate_C2E(this,fieldC,fieldE)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(:,:,:), intent(in) :: fieldC
   real(rkind), dimension(:,:,:), intent(out) :: fieldE
   call transpose_x_to_y(fieldC               , this%rbuffyC(:,:,:,1), this%gpC)
   call transpose_y_to_z(this%rbuffyC(:,:,:,1), this%rbuffzC(:,:,:,1), this%gpC)
   call this%PadeDer%interpz_C2E(this%rbuffzC(:,:,:,1), this%rbuffzE(:,:,:,1),0,0)
   call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
   call transpose_y_to_x(this%rbuffyE(:,:,:,1), fieldE               , this%gpE)
end subroutine

subroutine filter_and_interpolate_kappa_bounding(this,fieldC,fieldE)
   class(sgs_igrid), intent(inout) :: this
   real(rkind), dimension(:,:,:), intent(inout) :: fieldC, fieldE

   call this%gaussianX%filter1(fieldC, this%rbuffxC(:,:,:,1), this%gpC%xsz(2),this%gpC%xsz(3))
   call transpose_x_to_y(this%rbuffxC(:,:,:,1),this%rbuffyC(:,:,:,1), this%gpC)
   call this%gaussianY%filter2(this%rbuffyC(:,:,:,1), this%rbuffyC(:,:,:,2), this%gpC%ysz(1),this%gpC%ysz(3))
   call transpose_y_to_z(this%rbuffyC(:,:,:,2), this%rbuffzC(:,:,:,1), this%gpC)
   
   call this%gaussianZ%filter3(this%rbuffzC(:,:,:,1), this%rbuffzC(:,:,:,2), this%gpC%zsz(1),this%gpC%zsz(2))
   call this%PadeDer%interpz_C2E(this%rbuffzC(:,:,:,2), this%rbuffzE(:,:,:,1),0,0)
   
   call transpose_z_to_y(this%rbuffzC(:,:,:,2), this%rbuffyC(:,:,:,1), this%gpC)
   call transpose_y_to_x(this%rbuffyC(:,:,:,1), fieldC, this%gpC)

   call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
   call transpose_y_to_x(this%rbuffyE(:,:,:,1), fieldE, this%gpE)
end subroutine 
