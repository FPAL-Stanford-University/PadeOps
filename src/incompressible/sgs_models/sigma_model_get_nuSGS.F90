        G11 => this%rbuff(:,:,:,1); G12 => this%rbuff(:,:,:,2); G13 => this%rbuff(:,:,:,3)
        G22 => this%rbuff(:,:,:,4); G23 => this%rbuff(:,:,:,5); G33 => this%rbuff(:,:,:,6)
        
        I1 => this%rbuff(:,:,:,8); I2 => this%rbuff(:,:,:,9); I3 => this%rbuff(:,:,:,10)
        I1sq => this%rbuff(:,:,:,11); I1cu => this%rbuff(:,:,:,12)

        alpha1 => this%rbuff(:,:,:,1); alpha2 => this%rbuff(:,:,:,2); alpha3 => this%rbuff(:,:,:,3)

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
       
        alpha1sqrt => this%rbuff(:,:,:,4)
        alpha1sqrt = sqrt(alpha1)    
        alpha1tmp => this%rbuff(:,:,:,5) 
        alpha1tmp = alpha1*alpha1sqrt
        alpha1tmp = alpha2/(alpha1tmp)
        alpha1tmp = min(alpha1tmp,one)
        alpha1tmp = max(alpha1tmp,-one)
        alpha1tmp = acos(alpha1tmp)
        alpha3 = (one/three)*(alpha1tmp)
          
  
        sigma1 => this%rbuff(:,:,:,9); sigma2 => this%rbuff(:,:,:,10); sigma3 => this%rbuff(:,:,:,11)
        sigma1sq => this%rbuff(:,:,:,12)

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

        if (present(nuSGSfil)) then
            nuSGSfil = sigma3*(sigma1 - sigma2)*(sigma2 - sigma3)/(sigma1sq + 1.d-15)
            call this%spect%fft(nuSGSfil,this%nuSGShat)
            call this%spect%dealias(this%nuSGShat)
            call this%spect%ifft(this%nuSGShat,nuSGSfil)  
        else
            this%nuSGS = sigma3*(sigma1 - sigma2)*(sigma2 - sigma3)/(sigma1sq + 1.d-15)
            call this%spect%fft(this%nuSGS,this%nuSGShat)
            call this%spect%dealias(this%nuSGShat)
            call this%spect%ifft(this%nuSGShat,this%nuSGS)  
        end if 
