! Routines specific to 10th order Compact Finite Differencing scheme
! This file is included in the Derivatives module

#ifndef TEST_LU10
function InitLU10() result(ierr)

    integer :: ierr

    ! Based on  Neossi Nguetchue, Abelman (Appl. Math. & Comp. 2008)

    ! Allocate 1st derivative LU matrices.
    if(allocated( LU10X1 )) deallocate( LU10X1 ); allocate( LU10X1(nx,9) )
    if(allocated( LU10Y1 )) deallocate( LU10Y1 ); allocate( LU10Y1(ny,9) )
    if(allocated( LU10Z1 )) deallocate( LU10Z1 ); allocate( LU10Z1(nz,9) )

    ! Compute 1st derivative LU matrices
    if (nx .GE. 8) then
        call ComputeLU10(LU10X1,nx,beta10d1,alpha10d1,one,alpha10d1,beta10d1)
    else if (nx == 1) then
        LU10X1 = zero
    else
        ierr = 2
        return
    end if

    if (ny .GE. 8) then
        call ComputeLU10(LU10Y1,ny,beta10d1,alpha10d1,one,alpha10d1,beta10d1)
    else if (ny == 1) then
        LU10Y1 = zero
    else
        ierr = 2
        return
    end if

    if (nz .GE. 8) then
        call ComputeLU10(LU10Z1,nz,beta10d1,alpha10d1,one,alpha10d1,beta10d1)
    else if (nz == 1) then
        LU10Z1 = zero
    else
        ierr = 2
        return
    end if

    ! Allocate 2nd derivative LU matrices.
    if(allocated( LU10X2 )) deallocate( LU10X2 ); allocate( LU10X2(nx,9) )
    if(allocated( LU10Y2 )) deallocate( LU10Y2 ); allocate( LU10Y2(ny,9) )
    if(allocated( LU10Z2 )) deallocate( LU10Z2 ); allocate( LU10Z2(nz,9) )

    ! Compute 2nd derivative LU matrices
    if (nx .GE. 8) then
        call ComputeLU10(LU10X2,nx,beta10d2,alpha10d2,one,alpha10d2,beta10d2)
    else if (nx == 1) then
        LU10X2 = zero
    end if

    if (ny .GE. 8) then
        call ComputeLU10(LU10Y2,ny,beta10d2,alpha10d2,one,alpha10d2,beta10d2)
    else if (ny == 1) then
        LU10Y2 = zero
    end if

    if (nz .GE. 8) then
        call ComputeLU10(LU10Z2,nz,beta10d2,alpha10d2,one,alpha10d2,beta10d2)
    else if (nz == 1) then
        LU10Z2 = zero
    end if

    ! If everything passes
    ierr = 0

end function
#endif

subroutine ComputeLU10(LU,n,e,a,d,c,f)

#ifdef TEST_LU10
    use kind_parameters, only: rkind
#endif

    integer, intent(in) :: n
    real(rkind), intent(in) :: d,a,c,e,f
    real(rkind), dimension(n,9), intent(out) :: LU
    integer :: i

    LU = 0.0_rkind

    associate( b=>LU(:,1), eg=>LU(:,2), k=>LU(:,3),&
               l=>LU(:,4),  g=>LU(:,5), h=>LU(:,6),&
               ff=>LU(:,7),  v=>LU(:,8), w=>LU(:,9))
        
        ! Step 1       
        g(1) = d
        b(2) = a/g(1)
        h(1) = c
        k(1) = f/g(1)
        w(1) = a
        v(1) = e
        l(1) = c/g(1)
        g(2) = d - b(2)*h(1)
        k(2) = -k(1)*h(1)/g(2)
        w(2) = e - b(2)*w(1)
        v(2) = -b(2)*v(1)
        l(2) = (f - l(1)*h(1)) / g(2)
        h(2) = c - b(2)*f

        ! Step 2
        do i = 3,n-3
            b(i) = ( a - ( e/g(i-2) )*h(i-2) ) / g(i-1)
            h(i) = c - b(i)*f
            g(i) = d - ( e/g(i-2) )*f - b(i)*h(i-1)
        end do

        ! Step 3
        b(n-2) = ( a - ( e/g(n-4) )*h(n-4) ) / g(n-3)
        g(n-2) = d - ( e/g(n-4) )*f - b(n-2)*h(n-3)

        ! Step 4
        do i = 3,n-4
            k(i) = -( k(i-2)*f + k(i-1)*h(i-1) )/g(i)
            v(i) = -( e/g(i-2) )*v(i-2) - b(i)*v(i-1)
        end do

        ! Step 5
        k(n-3) = ( e - k(n-5)*f - k(n-4)*h(n-4) ) / g(n-3)
        k(n-2) = ( a - k(n-4)*f - k(n-3)*h(n-3) ) / g(n-2)
        v(n-3) = f - ( e/g(n-5) )*v(n-5) - b(n-3)*v(n-4)
        v(n-2) = c - ( e/g(n-4) )*v(n-4) - b(n-2)*v(n-3)
        g(n-1) = d - SUM( k(1:n-2)*v(1:n-2) )

        ! Step 6
        do i = 3,n-3
            w(i) = -( e/g(i-2) )*w(i-2) - b(i)*w(i-1)
            l(i) = -( l(i-2)*f + l(i-1)*h(i-1) ) / g(i)
        end do

        ! Step 7
        w(n-2) = f - ( e/g(n-4) )*w(n-4) - b(n-2)*w(n-3)
        w(n-1) = c - SUM( k(1:n-2)*w(1:n-2) )
        l(n-2) = ( e - l(n-4)*f - l(n-3)*h(n-3) ) / g(n-2)
        l(n-1) = ( a - SUM( l(1:n-2)*v(1:n-2) ) ) / g(n-1)
        g(n)   = d - SUM( l(1:n-1)*w(1:n-1) )

        ! Set eg(i) = e/g(i-2)
        eg(3:n-2) = e/g(1:n-4)

        ! Set ff = f
        ff(1:n-4) = f

    end associate

end subroutine


#ifndef TEST_LU10
function InitCD10( direction ) result(ierr)

    character(len=1), intent(in) :: direction
    integer :: ierr

    ierr = 0
end function
#endif
