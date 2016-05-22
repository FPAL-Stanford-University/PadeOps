        real(rkind), parameter :: b = (one/ten)/two , a = (three/two)/two

        if (this%isbotSided) then
            rhs(:,:,1) = 0.5d0*(fE(:,:,1) + fE(:,:,2)) 
        else 
            if (this%isBotEven) then
                rhs(:,:,1) = (b)*(fE(:,:,3) + fE(:,:,2)) + (a)*(fE(:,:,2) + fE(:,:,1))
            else
                rhs(:,:,1) = (b)*(fE(:,:,3) - fE(:,:,2)) + (a)*(fE(:,:,2) + fE(:,:,1))
            end if 
        end if 
        rhs(:,:,2:this%n-1) = (b)*(fE(:,:,4:this%nE) + fE(:,:,1:this%nE-3)) + (a)*(fE(:,:,3:this%nE-1) + fE(:,:,2:this%nE-2))
        
        if (this%isTopEven) then
            rhs(:,:,this%n) = (b)*(fE(:,:,this%nE-1) + fE(:,:,this%nE-2)) + (a)*(fE(:,:,this%nE) + fE(:,:,this%nE-1))
        else
            rhs(:,:,this%n) = (b)*(-fE(:,:,this%nE-1) + fE(:,:,this%nE-2)) + (a)*(fE(:,:,this%nE) + fE(:,:,this%nE-1))
        end if 
