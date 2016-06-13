    subroutine SolveZTriREAL(n,TriData,y,n1,n2)
        integer, intent(in) :: n, n1,n2
        real(rkind), intent(in), dimension(n,3) :: TriData
        real(rkind), dimension(n1,n2,n), intent(inout) :: y  
        integer :: k
  
        y(:,:,1) = y(:,:,1)*TriData(1,2)
        do k = 2,n
            y(:,:,k) = y(:,:,k)*TriData(k,2) - y(:,:,k-1)*TriData(k,1)
        end do
        
        do k = n-1,1,-1
            y(:,:,k) = y(:,:,k) - TriData(k,3)*y(:,:,k+1)
        end do 
        
    end subroutine

    subroutine SolveZTriCMPLX(n,TriData,y,n1,n2)
        integer, intent(in) :: n, n1,n2
        real(rkind), intent(in), dimension(n,3) :: TriData
        complex(rkind), dimension(n1,n2,n), intent(inout) :: y  
        integer :: k
  
        y(:,:,1) = y(:,:,1)*TriData(1,2)
        do k = 2,n
            y(:,:,k) = y(:,:,k)*TriData(k,2) - y(:,:,k-1)*TriData(k,1)
        end do
        
        do k = n-1,1,-1
            y(:,:,k) = y(:,:,k) - TriData(k,3)*y(:,:,k+1)
        end do 
        
    end subroutine
