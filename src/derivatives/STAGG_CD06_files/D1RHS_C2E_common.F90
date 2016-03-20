        real(rkind), parameter :: a = (63._rkind/62._rkind)/one , b = (17._rkind/62._rkind)/three
        real(rkind) :: a06, b06
        a06 = a * this%onebydx; b06 = b * this%onebydx

        rhs(:,:,2:this%nE-1) = a06*(fC(:,:,2:this%n) - fC(:,:,1:this%n-1))
        rhs(:,:,3:this%nE-2) = rhs(:,:,3:this%nE-2) + b06*(fC(:,:,4:this%n) - fC(:,:,1:this%n-3))

        if (this%isBotEven) then
            rhs(:,:,1) = zero
            rhs(:,:,2) = rhs(:,:,2) + b06*(fC(:,:,3) - fC(:,:,1))
        else
            rhs(:,:,1) = two*b06*fC(:,:,2) + two*a06*fC(:,:,1)
            rhs(:,:,2) = rhs(:,:,2) + b06*(fC(:,:,3) + fC(:,:,1))
        end if 

        if (this%isTopEven) then
            rhs(:,:,this%nE)    = zero 
            rhs(:,:,this%nE-1) = rhs(:,:,this%nE-1) + b06*(fC(:,:,this%n) - fC(:,:,this%n-2))
        else
            rhs(:,:,this%nE)    = -two*b06*(fC(:,:,this%n-1)) - two*a06*(fC(:,:,this%n))
            rhs(:,:,this%nE-1) = rhs(:,:,this%nE-1) - b06*(fC(:,:,this%n) + fC(:,:,this%n-2))
        end if 
