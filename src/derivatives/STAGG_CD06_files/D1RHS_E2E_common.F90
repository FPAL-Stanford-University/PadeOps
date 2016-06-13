        real(rkind), parameter :: a = (14._rkind/9._rkind)/two , b = (1._rkind/9._rkind)/four
        real(rkind) :: a06, b06
        a06 = a * this%onebydx; b06 = b * this%onebydx

        rhs(:,:,2:this%nE-1) = a06*(fE(:,:,3:this%nE) - fE(:,:,1:this%nE-2))
        rhs(:,:,3:this%nE-2) = rhs(:,:,3:this%nE-2) + b06*(fE(:,:,5:this%nE) - fE(:,:,1:this%nE-4))
 
        if (this%isBotEven) then
            rhs(:,:,1) = 0
            rhs(:,:,2) = rhs(:,:,2) + b06*(fE(:,:,4) - fE(:,:,2))
        else
            rhs(:,:,1) = b06*(fE(:,:,3) + fE(:,:,3)) + a06*(fE(:,:,2) + fE(:,:,2))
            rhs(:,:,2) = rhs(:,:,2) + b06*(fE(:,:,4) + fE(:,:,2)) 
        end if 

        if (this%isTopEven) then
            rhs(:,:,this%nE) = 0
            rhs(:,:,this%nE-1) = rhs(:,:,this%nE-1) + b06*(fE(:,:,this%nE-1) - fE(:,:,this%nE-3))
        else
            rhs(:,:,this%nE)   = -b06*(fE(:,:,this%nE-2) + fE(:,:,this%nE-2)) - a06*(fE(:,:,this%nE-1) + fE(:,:,this%nE-1))
            rhs(:,:,this%nE-1) = rhs(:,:,this%nE-1) - b06*(fE(:,:,this%nE-1) + fE(:,:,this%nE-3))
        end if 
