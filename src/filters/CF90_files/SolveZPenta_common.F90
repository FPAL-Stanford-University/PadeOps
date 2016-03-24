        ! Step 1
        y(:,:,2) = y(:,:,2) - penta(2,8)*y(:,:,1)
        do k = 3,this%n
            y(:,:,k) = y(:,:,k) - penta(k,9)*y(:,:,k-2) - penta(k,8)*y(:,:,k-1)
        end do 

        ! Step 2
        y(:,:,this%n) = y(:,:,this%n)*penta(this%n,7)
        
        y(:,:,this%n-1) = y(:,:,this%n-1)*penta(this%n-1,7) - penta(this%n-1,10)*y(:,:,this%n)
        do k = this%n-2,1,-1
            y(:,:,k) = y(:,:,k)*penta(k,7) - y(:,:,k+2)*penta(k,5)*penta(k,7) - y(:,:,k+1)*penta(k,10)
        end do 
