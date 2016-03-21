        real(rkind), parameter :: a = 12._rkind/11._rkind , b = (3._rkind/11._rkind)/four
        real(rkind) :: a06, b06

        a06 = a * this%onebydx2; b06 = b * this%onebydx2

        rhs(:,:,2:this%nE-1) = a06*(fE(:,:,3:this%nE) + fE(:,:,1:this%nE-2))
        rhs(:,:,3:this%nE-2) = rhs(:,:,3:this%nE-2) + b06*(fE(:,:,5:this%nE) + fE(:,:,1:this%nE-4)) 
 
        if (this%isBotEven) then
            rhs(:,:,1) = b06*(fE(:,:,3) + fE(:,:,3)) + a06*(fE(:,:,2) + fE(:,:,2)) 
            rhs(:,:,2) = rhs(:,:,2) + b06*(fE(:,:,4) + fE(:,:,2)) 
        else
            rhs(:,:,1) = zero 
            rhs(:,:,2) = rhs(:,:,2) + b06*(fE(:,:,4) - fE(:,:,2)) 
        end if 

        if (this%isTopEven) then
            rhs(:,:,this%nE) = two*b06*fE(:,:,this%nE-2) + two*a06*fE(:,:,this%nE-1)
            rhs(:,:,this%nE-1) = rhs(:,:,this%nE-1) + b06*(fE(:,:,this%nE-1) + fE(:,:,this%nE-3))
        else
            rhs(:,:,this%nE) = zero
            rhs(:,:,this%nE-1)= rhs(:,:,this%nE-1) + b06*(-fE(:,:,this%nE-1) + fE(:,:,this%nE-3))
        end if

        rhs = rhs - two*(b06 + a06)*fE 
