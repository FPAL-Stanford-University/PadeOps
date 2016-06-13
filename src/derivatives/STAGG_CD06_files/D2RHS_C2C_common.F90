        real(rkind), parameter :: a = 12._rkind/11._rkind , b = (3._rkind/11._rkind)/four
        real(rkind) :: a06, b06

        a06 = a * this%onebydx2; b06 = b * this%onebydx2

        rhs(:,:,2:this%n-1) = a06*(fC(:,:,3:this%n) + fC(:,:,1:this%n-2))
        rhs(:,:,3:this%n-2) = rhs(:,:,3:this%n-2) + b06*(fC(:,:,5:this%n) + fC(:,:,1:this%n-4))

 
        if (this%isBotEven) then
            rhs(:,:,1) = b06*(fC(:,:,3) + fC(:,:,2)) + a06*(fC(:,:,2) + fC(:,:,1)) 
            rhs(:,:,2) = rhs(:,:,2) + b06*(fC(:,:,4) + fC(:,:,1))
        else
            rhs(:,:,1) = b06*(fC(:,:,3) - fC(:,:,2)) + a06*(fC(:,:,2) - fC(:,:,1))
            rhs(:,:,2) = rhs(:,:,2) + b06*(fC(:,:,4) - fC(:,:,1))
        end if 

        if (this%isTopEven) then
            rhs(:,:,this%n)  =b06*(fC(:,:,this%n-1) + fC(:,:,this%n-2)) + a06*(fC(:,:,this%n) + fC(:,:,this%n-1))
            rhs(:,:,this%n-1)=rhs(:,:,this%n-1)                     + b06*(fC(:,:,this%n) + fC(:,:,this%n-3)) 
        else
            rhs(:,:,this%n) =b06*(-fC(:,:,this%n-1) + fC(:,:,this%n-2)) + a06*(-fC(:,:,this%n) + fC(:,:,this%n-1))
            rhs(:,:,this%n-1)=rhs(:,:,this%n-1)                     + b06*(-fC(:,:,this%n) + fC(:,:,this%n-3)) 
        end if

        rhs = rhs - two*(b06 + a06)*fC 
