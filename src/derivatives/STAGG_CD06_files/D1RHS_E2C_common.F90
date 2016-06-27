        real(rkind), parameter :: a = (63._rkind/62._rkind)/one , b = (17._rkind/62._rkind)/three
        real(rkind), parameter :: alpha1 = 37.d0/183.d0, alpha0 = -1.d0
        real(rkind), parameter :: a0 = (1.d0/24.d0)*(alpha0 - 23.d0), b0 = (1.d0/8.d0)*(-9.d0*alpha0 + 7.d0)
        real(rkind), parameter :: c0 = (1.d0/8.d0 )*(9.d0*alpha0 + 1.d0),d0 = -(1.d0/24.d0)*(alpha0 + 1.d0)
        real(rkind), parameter :: a1 = (3.d0/8.d0)*(3.d0 - 2.d0*alpha1), b1 = (1.d0/8.d0)*(-1.d0 + 22.d0*alpha1)
        real(rkind) :: a06, b06
        a06 = a * this%onebydx; b06 = b * this%onebydx

        rhs(:,:,2:this%n-1) = b06*(fE(:,:,4:this%nE)   - fE(:,:,1:this%nE-3)) &
                        + a06*(fE(:,:,3:this%nE-1) - fE(:,:,2:this%nE-2)) 

        if (this%isBotSided) then
            rhs(:,:,1) = w0s*this%onebydx*(a0*fE(:,:,1) + b0*fE(:,:,2) + &
                         & c0*fE(:,:,3) + d0*fE(:,:,4))
            rhs(:,:,2) = w1s*this%onebydx*((-b1/3.d0)*fE(:,:,1) + (-a1)*fE(:,:,2) + &
                         & (a1)*fE(:,:,3) + (b1/3.d0)*fE(:,:,4))
        else 
            if (this%isBotEven) then
                rhs(:,:,1)  = b06*(fE(:,:,3) - fE(:,:,2)) &
                            + a06*(fE(:,:,2) - fE(:,:,1)) 
            else
                rhs(:,:,1)  = b06*(fE(:,:,3) + fE(:,:,2)) &
                            + a06*(fE(:,:,2) - fE(:,:,1)) 
            end if 
        end if 

        if (this%isTopSided) then
            rhs(:,:,this%n  ) = -w0s*this%onebydx*(a0*fE(:,:,this%nE) + b0*fE(:,:,this%nE-1) + &
                         & c0*fE(:,:,this%nE-2) + d0*fE(:,:,this%nE-3))
            rhs(:,:,this%n-1) = w1s*this%onebydx*((-b1/3.d0)*fE(:,:,this%nE-3) + (-a1)*fE(:,:,this%nE-2) + &
                         & (a1)*fE(:,:,this%nE-1) + (b1/3.d0)*fE(:,:,this%nE))
        else 
            if (this%isTopEven) then
                rhs(:,:,this%n) = b06*(fE(:,:,this%nE-1) - fE(:,:,this%nE-2))  &
                            + a06*(fE(:,:,this%nE  ) - fE(:,:,this%nE-1)) 
            else
                rhs(:,:,this%n) = -b06*(fE(:,:,this%nE-1) + fE(:,:,this%nE-2)) &
                            +  a06*(fE(:,:,this%nE  ) - fE(:,:,this%nE-1))
            end if 
        end if 
