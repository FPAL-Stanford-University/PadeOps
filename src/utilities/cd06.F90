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

subroutine SolveLU06(LU,k,n)

#ifdef TEST_LU06
    use kind_parameters, only: rkind
#endif

    integer, intent(in) :: n
    real(rkind), dimension(n,5), intent(in)  :: LU
    real(rkind), dimension(n), intent(inout) :: k  ! Take in RHS and put solution into it
    integer :: i

    associate ( bc=>LU(:,1), h=>LU(:,2), onebyc=>LU(:,3), a=>LU(:,4), v=>LU(:,5) )
        
        ! Step 2
        do i = 2,n-1
            k(i) = k(i) - bc(i)*k(i-1)
        end do
        k(n) = k(n) - SUM( h(1:n-1)*k(1:n-1) )

        ! Step 3
        k(n) = k(n) * onebyc(n)
        k(n-1) = ( k(n-1) - v(n-1)*k(n) ) * onebyc(n-1)
        do i = n-2,1,-1
            k(i) = ( k(i) - a(i)*k(i+1) - v(i)*k(n) ) * onebyc(i)
        end do

    end associate


end subroutine

pure function CD06D1RHS(f,n,periodic) result (RHS)

#ifdef TEST_LU06
    use kind_parameters, only: rkind
#include "PadeCoeffs.F90"
#endif

    integer, intent(in) :: n
    real(rkind), dimension(n), intent(in) :: f
    logical, intent(in) :: periodic
    real(rkind), dimension(n) :: RHS
    integer :: i

    select case (periodic)
    case (.TRUE.)
        RHS(1) = a06d1 * ( f(2)   - f(n)   ) &
               + b06d1 * ( f(3)   - f(n-1) ) 
        RHS(2) = a06d1 * ( f(3)   - f(1)   ) &
               + b06d1 * ( f(4)   - f(n)   ) 
        do i = 3,n-2
            RHS(i) = a06d1 * ( f(i+1) - f(i-1) ) &
                   + b06d1 * ( f(i+2) - f(i-2) ) 
        end do
        RHS(n-1) = a06d1 * ( f(n)   - f(n-2) ) &
                 + b06d1 * ( f(1)   - f(n-3) ) 
        RHS(n)   = a06d1 * ( f(1)   - f(n-1) ) &
                 + b06d1 * ( f(2)   - f(n-2) )
    case (.FALSE.)
        RHS = zero
    end select

end function

pure function CD06D2RHS(f,n,periodic) result (RHS)

#ifdef TEST_LU06
    use kind_parameters, only: rkind
#include "PadeCoeffs.F90"
#endif

    integer, intent(in) :: n
    real(rkind), dimension(n), intent(in) :: f
    logical, intent(in) :: periodic
    real(rkind), dimension(n) :: RHS
    integer :: i

    select case (periodic)
    case (.TRUE.)
        RHS(1) = a06d2 * ( f(2)   - two*f(1) + f(n)   ) &
               + b06d2 * ( f(3)   - two*f(1) + f(n-1) ) 
        RHS(2) = a06d2 * ( f(3)   - two*f(2) + f(1)   ) &
               + b06d2 * ( f(4)   - two*f(2) + f(n)   ) 
        do i = 3,n-2
            RHS(i) = a06d2 * ( f(i+1) - two*f(i) + f(i-1) ) &
                   + b06d2 * ( f(i+2) - two*f(i) + f(i-2) ) 
        end do
        RHS(n-1) = a06d2 * ( f(n)   - two*f(n-1) + f(n-2) ) &
                 + b06d2 * ( f(1)   - two*f(n-1) + f(n-3) ) 
        RHS(n)   = a06d2 * ( f(1)   - two*f(n)   + f(n-1) ) &
                 + b06d2 * ( f(2)   - two*f(n)   + f(n-2) )
    case (.FALSE.)
        RHS = zero
    end select

end function

pure function d1xCD06(f) result(df)
    real(rkind), dimension(nx), intent(in) :: f
    real(rkind), dimension(nx) :: df
    real(rkind), dimension(nx) :: tmp

    tmp = CD06D1RHS(f,nx,periodicx)  
    if (periodicx) then
        call SolveLU06(LU06X1,tmp,nx)
        df = tmp*onebydx
    else
        df = zero
    end if 

end function

pure function d1yCD06(f) result(df)
    real(rkind), dimension(ny), intent(in) :: f
    real(rkind), dimension(ny) :: df
    real(rkind), dimension(ny) :: tmp

    tmp = CD06D1RHS(f,ny,periodicy)  
    if (periodicy) then
        call SolveLU06(LU06Y1,tmp,ny)
        df = tmp*onebydy
    else
        df = zero
    end if 

end function


pure function d1zCD06(f) result(df)
    real(rkind), dimension(nz), intent(in) :: f
    real(rkind), dimension(nz) :: df
    real(rkind), dimension(nz) :: tmp

    tmp = CD06D1RHS(f,nz,periodicz)  
    if (periodicz) then
        call SolveLU06(LU06Z1,tmp,nz)
        df = tmp*onebydz
    else
        df = zero
    end if 

end function


pure function d2xCD06(f) result(df)
    real(rkind), dimension(nx), intent(in) :: f
    real(rkind), dimension(nx) :: df
    real(rkind), dimension(nx) :: tmp

    tmp = CD06D2RHS(f,nx,periodicx)  
    if (periodicx) then
        call SolveLU06(LU06X2,tmp,nx)
        df = tmp*onebydx*onebydx
    else
        df = zero
    end if 

end function

pure function d2yCD06(f) result(df)
    real(rkind), dimension(ny), intent(in) :: f
    real(rkind), dimension(ny) :: df
    real(rkind), dimension(ny) :: tmp

    tmp = CD06D2RHS(f,ny,periodicy)  
    if (periodicy) then
        call SolveLU06(LU06Y2,tmp,ny)
        df = tmp*onebydy*onebydy
    else
        df = zero
    end if 

end function


pure function d2zCD06(f) result(df)
    real(rkind), dimension(nz), intent(in) :: f
    real(rkind), dimension(nz) :: df
    real(rkind), dimension(nz) :: tmp

    tmp = CD06D2RHS(f,nz,periodicz)  
    if (periodicz) then
        call SolveLU06(LU06Z2,tmp,nz)
        df = tmp*onebydz
    else
        df = zero
    end if 

end function

#ifndef TEST_LU06
function InitCD06( direction ) result(ierr)

    character(len=1), intent(in) :: direction
    integer :: ierr

    ierr = 0
end function
#endif
