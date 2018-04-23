            if (present(nuSGSfil)) then ! This is the standard nuSGS
                S11 => this%rbuff(:,:,:,1); S12 => this%rbuff(:,:,:,2); S13 => this%rbuff(:,:,:,3)
                S22 => this%rbuff(:,:,:,4); S23 => this%rbuff(:,:,:,5); S33 => this%rbuff(:,:,:,6)
                nuSGSfil = S12*S12 + S13*S13 + S23*S23
                nuSGSfil = two*nuSGSfil
                nuSGSfil = nuSGSfil + S11*S11 + S22*S22 + S33*S33
                nuSGSfil = two*nuSGSfil
                nuSGSfil = sqrt(nuSGSfil)
            else ! This is the filtered usual nuSGS
                S11 => this%rbuff(:,:,:,13); S12 => this%rbuff(:,:,:,14); S13 => this%rbuff(:,:,:,15)
                S22 => this%rbuff(:,:,:,16); S23 => this%rbuff(:,:,:,17); S33 => this%rbuff(:,:,:,18)
                this%nuSGS = S12*S12 + S13*S13 + S23*S23
                this%nuSGS = two*this%nuSGS
                this%nuSGS = this%nuSGS + S11*S11 + S22*S22 + S33*S33
                this%nuSGS = two*this%nuSGS
                this%nuSGS = sqrt(this%nuSGS)
            end if  
            nullify(S11,S22,S33,S12,S13,S23)

