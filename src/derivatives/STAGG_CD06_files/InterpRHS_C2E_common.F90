        real(rkind), parameter :: b = (one/ten)/two , a = (three/two)/two
        real(rkind), parameter :: a0 = 15.d0/8.d0, b0 = -5.d0/4.d0, c0 = 3.d0/8.d0
        real(rkind), parameter :: alpha1 = 1.d0/6.d0, a1 = (1.d0/8.d0)*(9.d0 + 10.d0*alpha1)

        rhs(:,:,2:this%nE-1) = (a)*(fC(:,:,2:this%n) + fC(:,:,1:this%n-1))
        rhs(:,:,3:this%nE-2) = rhs(:,:,3:this%nE-2) + (b)*(fC(:,:,4:this%n) + fC(:,:,1:this%n-3)) 

        if (this%isBotSided) then
            rhs(:,:,1) = a0*fC(:,:,1) + b0*fC(:,:,2) + c0*fC(:,:,3) 
            rhs(:,:,2) = (a1/2.d0)*(fC(:,:,2) + fC(:,:,1))
        else
            if (this%isBotEven) then
                rhs(:,:,1) = two*b*fC(:,:,2) + two*a*fC(:,:,1)
                rhs(:,:,2) = rhs(:,:,2) + (b)*(fC(:,:,3) + fC(:,:,1))
            else
                rhs(:,:,1) = zero
                rhs(:,:,2) = rhs(:,:,2) + (b)*(fC(:,:,3) - fC(:,:,1))
            end if 
        end if 

        if (this%isTopSided) then
            rhs(:,:,this%nE) = a0*fC(:,:,this%n) + b0*fC(:,:,this%n-1) + c0*fC(:,:,this%n-2) 
            rhs(:,:,this%nE-1) = (a1/2.d0)*(fC(:,:,this%n) + fC(:,:,this%n-1))
        else
            if (this%isTopEven) then
                rhs(:,:,this%nE) =   two*b*fC(:,:,this%n-1) + two*a*fC(:,:,this%n)
                rhs(:,:,this%nE-1) = rhs(:,:,this%nE-1) + (b)*(fC(:,:,this%n) + fC(:,:,this%n-2))
            else
                rhs(:,:,this%nE) = zero
                rhs(:,:,this%nE-1) = rhs(:,:,this%nE-1) + (b)*(-fC(:,:,this%n) + fC(:,:,this%n-2)) 
            end if
        end if  
