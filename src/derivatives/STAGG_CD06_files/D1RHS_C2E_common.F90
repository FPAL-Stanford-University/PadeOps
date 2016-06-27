        real(rkind), parameter :: a = (63._rkind/62._rkind)/one , b = (17._rkind/62._rkind)/three
        real(rkind), parameter :: a0 = -71.d0/24.d0, b0 = 47.d0/ 8.d0
        real(rkind), parameter :: c0 = -31.d0/ 8.d0, d0 = 23.d0/24.d0
        real(rkind), parameter :: a1 = 12.d0/11.d0
        real(rkind) :: a06, b06
        a06 = a * this%onebydx; b06 = b * this%onebydx

        rhs(:,:,2:this%nE-1) = a06*(fC(:,:,2:this%n) - fC(:,:,1:this%n-1))
        rhs(:,:,3:this%nE-2) = rhs(:,:,3:this%nE-2) + b06*(fC(:,:,4:this%n) - fC(:,:,1:this%n-3))

        if (this%isBotSided) then
            rhs(:,:,1) = (a0*fC(:,:,1) + b0*fC(:,:,2) + c0*fC(:,:,3) + d0*fC(:,:,4))*this%onebydx 
            rhs(:,:,2) = (fC(:,:,2) - fC(:,:,1))*(a1*this%onebydx)
        else
            if (this%isBotEven) then
                rhs(:,:,1) = zero
                rhs(:,:,2) = rhs(:,:,2) + b06*(fC(:,:,3) - fC(:,:,1))
            else
                rhs(:,:,1) = two*b06*fC(:,:,2) + two*a06*fC(:,:,1)
                rhs(:,:,2) = rhs(:,:,2) + b06*(fC(:,:,3) + fC(:,:,1))
            end if 
        end if 

        if (this%isTopSided) then
            rhs(:,:,this%nE-1) = (fC(:,:,this%n) - fC(:,:,this%n-1))*(a1*this%onebydx)
            rhs(:,:,this%nE) = (a0*fC(:,:,this%n) + b0*fC(:,:,this%n-1) + c0*fC(:,:,this%n-2) &
                            & + d0*fC(:,:,this%n-3))*(-this%onebydx)
        else
            if (this%isTopEven) then
                rhs(:,:,this%nE-1) = rhs(:,:,this%nE-1) + b06*(fC(:,:,this%n) - fC(:,:,this%n-2))
                rhs(:,:,this%nE)    = zero 
            else
                rhs(:,:,this%nE-1) = rhs(:,:,this%nE-1) - b06*(fC(:,:,this%n) + fC(:,:,this%n-2))
                rhs(:,:,this%nE)    = -two*b06*(fC(:,:,this%n-1)) - two*a06*(fC(:,:,this%n))
            end if 
        end if 
