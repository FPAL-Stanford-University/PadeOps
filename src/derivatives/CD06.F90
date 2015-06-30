! Routines specific to 6th order Compact Finite Differencing scheme
! This file is included in the Derivatives module

#ifndef TEST_LU06
function InitLU06() result(ierr)

    integer :: ierr

    ! Based on  Neossi Nguetchue, Abelman (Appl. Math. & Comp. 2008)

    ! Allocate 1st derivative LU matrices.
    if(allocated( LU06X1 )) deallocate( LU06X1 ); allocate( LU06X1(nx,5) )
    if(allocated( LU06Y1 )) deallocate( LU06Y1 ); allocate( LU06Y1(ny,5) )
    if(allocated( LU06Z1 )) deallocate( LU06Z1 ); allocate( LU06Z1(nz,5) )

    ! Compute 1st derivative LU matrices
    if (nx .GE. 6) then
        call ComputeLU06(LU06X1,nx,alpha06d1,one,alpha06d1)
    else if (nx == 1) then
        LU06X1 = zero
    else
        ierr = 2
        return
    end if

    if (ny .GE. 6) then
        call ComputeLU06(LU06Y1,ny,alpha06d1,one,alpha06d1)
    else if (ny == 1) then
        LU06Y1 = zero
    else
        ierr = 2
        return
    end if

    if (nz .GE. 6) then
        call ComputeLU06(LU06Z1,nz,alpha06d1,one,alpha06d1)
    else if (nz == 1) then
        LU06Z1 = zero
    else
        ierr = 2
        return
    end if

    ! Allocate 2nd derivative LU matrices.
    if(allocated( LU06X2 )) deallocate( LU06X2 ); allocate( LU06X2(nx,5) )
    if(allocated( LU06Y2 )) deallocate( LU06Y2 ); allocate( LU06Y2(ny,5) )
    if(allocated( LU06Z2 )) deallocate( LU06Z2 ); allocate( LU06Z2(nz,5) )

    ! Compute 2nd derivative LU matrices
    if (nx .GE. 6) then
        call ComputeLU06(LU06X2,nx,alpha06d2,one,alpha06d2)
    else if (nx == 1) then
        LU06X2 = zero
    end if

    if (ny .GE. 6) then
        call ComputeLU06(LU06Y2,ny,alpha06d2,one,alpha06d2)
    else if (ny == 1) then
        LU06Y2 = zero
    end if

    if (nz .GE. 6) then
        call ComputeLU06(LU06Z2,nz,alpha06d2,one,alpha06d2)
    else if (nz == 1) then
        LU06Z2 = zero
    end if

    ! If everything passes
    ierr = 0

end function
#endif

subroutine ComputeLU06(LU,n,b,d,a)

#ifdef TEST_LU06
    use kind_parameters, only: rkind
#endif

    integer, intent(in) :: n
    real(rkind), intent(in) :: d,a,b
    real(rkind), dimension(n,9), intent(out) :: LU
    integer :: i

    LU = 0.0_rkind

    associate ( bc=>LU(:,1), h=>LU(:,2), c=>LU(:,3), aa=>LU(:,4), v=>LU(:,5) )

        ! Step 0
        c(1) = d
        v(1) = b
        h(1) = a/c(1)

        ! Step 1
        do i = 2,n-1
            c(i) = d - ( b/c(i-1) )*a
        end do
        do i = 2,n-2
            v(i) = -( b/c(i-1) )*v(i-1)
            h(i) = -( a/c(i) )*h(i-1)
        end do
        v(n-1) = a - ( b/c(n-2) )*v(n-2)
        h(n-1) = ( b - h(n-2)*a ) / c(n-1)
        c(n)   = d - SUM( h(1:n-1)*v(1:n-1) )

        bc(2:n-1) = b / c(1:n-2)
        aa(1:n-2) = a

        ! Set c = 1/c
        c = 1._rkind/c

    end associate

end subroutine

subroutine SolveLU06(LU,k,n,nvec)

#ifdef TEST_LU06
    use kind_parameters, only: rkind
#endif

    integer, intent(in) :: n, nvec
    real(rkind), dimension(n,5), intent(in)  :: LU
    real(rkind), dimension(n,nvec), intent(inout) :: k  ! Take in RHS and put solution into it
    integer :: i

    associate ( bc=>LU(:,1), h=>LU(:,2), onebyc=>LU(:,3), a=>LU(:,4), v=>LU(:,5) )
        
        ! Step 2
        do i = 2,n-1
            k(i,:) = k(i,:) - bc(i)*k(i-1,:)
        end do
        do i=1,nvec
            k(n,i) = k(n,i) - SUM( h(1:n-1)*k(1:n-1,i) )
        end do

        ! Step 3
        k(n,:) = k(n,:) * onebyc(n)
        k(n-1,:) = ( k(n-1,:) - v(n-1)*k(n,:) ) * onebyc(n-1)
        do i = n-2,1,-1
            k(i,:) = ( k(i,:) - a(i)*k(i+1,:) - v(i)*k(n,:) ) * onebyc(i)
        end do

    end associate


end subroutine

#ifndef TEST_LU06
function InitCD06( direction ) result(ierr)

    character(len=1), intent(in) :: direction
    integer :: ierr

    ierr = 0
end function
#endif
