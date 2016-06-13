        real(rkind), parameter :: a = (63._rkind/62._rkind)/one , b = (17._rkind/62._rkind)/three
        real(rkind) :: a06, b06
        a06 = a * this%onebydx; b06 = b * this%onebydx

        if (this%isBotSided) then
            rhs(:,:,1) = this%onebydx*(fE(:,:,2) - fE(:,:,1))
        else 
            if (this%isBotEven) then
                rhs(:,:,1)  = b06*(fE(:,:,3) - fE(:,:,2)) &
                            + a06*(fE(:,:,2) - fE(:,:,1)) 
            else
                rhs(:,:,1)  = b06*(fE(:,:,3) + fE(:,:,2)) &
                            + a06*(fE(:,:,2) - fE(:,:,1)) 
            end if 
        end if 

        rhs(:,:,2:this%n-1) = b06*(fE(:,:,4:this%nE)   - fE(:,:,1:this%nE-3)) &
                        + a06*(fE(:,:,3:this%nE-1) - fE(:,:,2:this%nE-2)) 
        
        if (this%isTopEven) then
            rhs(:,:,this%n) = b06*(fE(:,:,this%nE-1) - fE(:,:,this%nE-2))  &
                        + a06*(fE(:,:,this%nE  ) - fE(:,:,this%nE-1)) 
        else
            rhs(:,:,this%n) = -b06*(fE(:,:,this%nE-1) + fE(:,:,this%nE-2)) &
                        +  a06*(fE(:,:,this%nE  ) - fE(:,:,this%nE-1))
        end if 
