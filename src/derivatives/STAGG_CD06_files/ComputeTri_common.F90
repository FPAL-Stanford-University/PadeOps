        cp(1) = dup(1)/dg(1)
        do i = 2,n-1
            cp(i) = dup(i)/(dg(i) - ddn(i)*cp(i-1))
        end do
            
        den(1) = one/dg(1)
        den(2:n) = one/(dg(2:n) - ddn(2:n)*cp(1:n-1))
