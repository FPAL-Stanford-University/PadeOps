        ! Step 1
        y(:,:,2) = y(:,:,2) - this%penta(2,8)*y(:,:,1)
        do k = 3,this%n
            y(:,:,k) = y(:,:,k) - this%penta(k,9)*y(:,:,k-2) - this%penta(k,8)*y(:,:,k-1)
        end do 

        ! Step 2
        y(:,:,this%n) = y(:,:,this%n)*this%penta(this%n,7)
        
        y(:,:,this%n-1) = y(:,:,this%n-1)*this%penta(this%n-1,7) - this%penta(this%n-1,10)*y(:,:,this%n)
        do k = this%n-2,1,-1
            y(:,:,k) = y(:,:,k)*this%penta(k,7) - y(:,:,k+2)*this%penta(k,5)*this%penta(k,7) - y(:,:,k+1)*this%penta(k,10)
        end do 
