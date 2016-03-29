            if(.not. present(nuSGSfil)) then ! This is the standard nuSGS
                S11 => this%rbuff(:,:,:,13); S12 => this%rbuff(:,:,:,14); S13 => this%rbuff(:,:,:,15)
                S22 => this%rbuff(:,:,:,16); S23 => this%rbuff(:,:,:,17); S33 => this%rbuff(:,:,:,18)
                Snorm => this%nuSGS
            else ! This is the filtered usual nuSGS
                S11 => this%rbuff(:,:,:,1); S12 => this%rbuff(:,:,:,2); S13 => this%rbuff(:,:,:,3)
                S22 => this%rbuff(:,:,:,4); S23 => this%rbuff(:,:,:,5); S33 => this%rbuff(:,:,:,6)
                Snorm => nuSGSfil
            end if  
           
            Snorm = S12*S12 + S13*S13 + S23*S23
            Snorm = two*Snorm
            Snorm = Snorm + S11*S11 + S22*S22 + S33*S33
            Snorm = two*Snorm 
            Snorm = sqrt(Snorm)
