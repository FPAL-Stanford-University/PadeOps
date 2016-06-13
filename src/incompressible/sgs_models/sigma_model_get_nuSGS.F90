    subroutine get_SIGMA_Op(this, nuSGS, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
        class(sgs), intent(inout), target :: this
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: dudx, dudy, dudz 
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: dvdx, dvdy, dvdz 
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: dwdx, dwdy, dwdz 
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(out) :: nuSGS

        real(rkind), dimension(:,:,:), pointer :: G11, G12, G13, G22, G23, G33
        real(rkind), dimension(:,:,:), pointer :: I1, I2, I3, I1sq, I1cu
        real(rkind), dimension(:,:,:), pointer :: alpha1, alpha2, alpha3 
        real(rkind), dimension(:,:,:), pointer :: sigma1, sigma2, sigma3 
        real(rkind), dimension(:,:,:), pointer :: alpha1sqrt, sigma1sq, alpha1tmp 

        G11=>this%SIGMAbuffs(:,:,:,1); G12=>this%SIGMAbuffs(:,:,:,2); G13=>this%SIGMAbuffs(:,:,:,3)
        G22=>this%SIGMAbuffs(:,:,:,4); G23=>this%SIGMAbuffs(:,:,:,5); G33=>this%SIGMAbuffs(:,:,:,6)
        
        I1=>this%SIGMAbuffs(:,:,:,7); I2=>this%SIGMAbuffs(:,:,:,8); I3=>this%SIGMAbuffs(:,:,:,9)
        I1sq=>this%SIGMAbuffs(:,:,:,10); I1cu=>this%SIGMAbuffs(:,:,:,11)

        alpha1=>this%SIGMAbuffs(:,:,:,12); alpha2=>this%SIGMAbuffs(:,:,:,13); alpha3=>this%SIGMAbuffs(:,:,:,14)

        G11 = dudx*dudx + dvdx*dvdx + dwdx*dwdx 
        G12 = dudx*dudy + dvdx*dvdy + dwdx*dwdy 
        G13 = dudx*dudz + dvdx*dvdz + dwdx*dwdz
        G22 = dudy*dudy + dvdy*dvdy + dwdy*dwdy 
        G23 = dudy*dudz + dvdy*dvdz + dwdy*dwdz 
        G33 = dudz*dudz + dvdz*dvdz + dwdz*dwdz

        I1 = G11 + G22 + G33
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
       
        alpha1sqrt => this%SIGMAbuffs(:,:,:,4)
        alpha1sqrt = sqrt(alpha1)    
        alpha1tmp => this%SIGMAbuffs(:,:,:,5) 
        alpha1tmp = alpha1*alpha1sqrt
        alpha1tmp = alpha2/(alpha1tmp)
        alpha1tmp = min(alpha1tmp,one)
        alpha1tmp = max(alpha1tmp,-one)
        alpha1tmp = acos(alpha1tmp)
        alpha3 = (one/three)*(alpha1tmp)
  
        sigma1 => this%SIGMAbuffs(:,:,:,9); sigma2 => this%SIGMAbuffs(:,:,:,10); sigma3 => this%SIGMAbuffs(:,:,:,11)
        sigma1sq => this%SIGMAbuffs(:,:,:,12)

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
            
        nuSGS = sigma3*(sigma1 - sigma2)*(sigma2 - sigma3)/(sigma1sq + 1.d-15)
    end subroutine 
