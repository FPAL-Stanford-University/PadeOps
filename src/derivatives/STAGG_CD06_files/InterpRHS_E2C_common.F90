        real(rkind), parameter :: b = (one/ten)/two , a = (three/two)/two
        real(rkind), parameter :: alpha0 = 1.d0
        real(rkind), parameter :: a0 = (1.d0/16.d0)*(5.d0 - alpha0), b0 = (1.d0/16.d0)*(9.d0*alpha0 + 15.d0)
        real(rkind), parameter :: c0 = (1.d0/16.d0)*(9.d0*alpha0 - 5.d0), d0 = (1.d0/16.d0)*(1.d0 - alpha0)

        if (this%isbotSided) then
            rhs(:,:,1) = a0*fE(:,:,1) + b0*fE(:,:,2) + c0*fE(:,:,3) + d0*fE(:,:,4)
        else 
            if (this%isBotEven) then
                rhs(:,:,1) = (b)*(fE(:,:,3) + fE(:,:,2)) + (a)*(fE(:,:,2) + fE(:,:,1))
            else
                rhs(:,:,1) = (b)*(fE(:,:,3) - fE(:,:,2)) + (a)*(fE(:,:,2) + fE(:,:,1))
            end if 
        end if 
        rhs(:,:,2:this%n-1) = (b)*(fE(:,:,4:this%nE) + fE(:,:,1:this%nE-3)) + (a)*(fE(:,:,3:this%nE-1) + fE(:,:,2:this%nE-2))
       
        if (this%isTopSided) then 
            rhs(:,:,this%n) = a0*fE(:,:,this%nE  ) + b0*fE(:,:,this%nE-1) & 
                            + c0*fE(:,:,this%nE-2) + d0*fE(:,:,this%nE-3)
        else
            if (this%isTopEven) then
                rhs(:,:,this%n) = (b)*(fE(:,:,this%nE-1) + fE(:,:,this%nE-2)) + (a)*(fE(:,:,this%nE) + fE(:,:,this%nE-1))
            else
                rhs(:,:,this%n) = (b)*(-fE(:,:,this%nE-1) + fE(:,:,this%nE-2)) + (a)*(fE(:,:,this%nE) + fE(:,:,this%nE-1))
            end if
        end if 
