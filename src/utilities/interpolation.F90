module interpolation
    use kind_parameters, only: rkind
    use constants, only: one
    implicit none
    contains

    pure subroutine binarysearch(length, array, value, val,delta)
        ! Given an array and a value, returns the index of the element that
        ! is closest to, but less than, the given value.
        ! Uses a binary search algorithm.
        ! "delta" is the tolerance used to determine if two values are equal
        ! if ( abs(x1 - x2) <= delta) then
        !    assume x1 = x2
        ! endif

        implicit none
        integer, intent(in) :: length
        real(rkind), dimension(length), intent(in) :: array
        real(rkind), intent(in) :: value
        real(rkind), intent(in), optional :: delta

        integer, intent(out) :: val

        integer :: left, middle, right
        real(rkind) :: d

        if (present(delta) .eqv. .true.) then
            d = delta
        else
            d = 1d-9
        endif
        
        left = 1
        right = length
        do
            if (left > right) then
                exit
            endif
            middle = nint((left+right) / 2.0)
            if ( abs(array(middle) - value) <= d) then
                val = middle
                return
            else if (array(middle) > value) then
                right = middle - 1
            else
                left = middle + 1
            end if
        end do
        val = right

    end subroutine binarysearch

    pure function interpolate2d(x_len, x_array, y_len, y_array, f, x, y) result(val)
        ! This function uses bilinear interpolation to estimate the value
        ! of a function f at point (x,y)
        ! f is assumed to be sampled on a regular grid, with the grid x values specified
        ! by x_array and the grid y values specified by y_array
        ! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation
        implicit none
        integer, intent(in) :: x_len, y_len           
        real(rkind), dimension(x_len), intent(in) :: x_array
        real(rkind), dimension(y_len), intent(in) :: y_array
        real(rkind), dimension(x_len, y_len), intent(in) :: f
        real(rkind), intent(in) :: x,y

        real(rkind) :: denom, x1, x2, y1, y2
        integer :: i,j
        real(rkind) :: val

        call binarysearch(x_len, x_array, x, i)
        call binarysearch(y_len, y_array, y, j)

        x1 = x_array(i)
        x2 = x_array(i+1)

        y1 = y_array(j)
        y2 = y_array(j+1)
        
        denom = (x2 - x1)*(y2 - y1)

        val = (f(i,j)*(x2-x)*(y2-y) + f(i+1,j)*(x-x1)*(y2-y) + &
            f(i,j+1)*(x2-x)*(y-y1) + f(i+1, j+1)*(x-x1)*(y-y1))/denom

    end function interpolate2d

    pure subroutine interpolate3d(x_len, x_array, y_len, y_array, z_len, z_array, f, x, y, z, val)
        ! This function uses trilinear interpolation to estimate the value
        ! of a function f at point (x,y,z)
        ! f is assumed to be sampled on a regular grid, with the grid x values specified
        ! by x_array and the grid y values specified by y_array
        ! Reference: http://en.wikipedia.org/wiki/Trilinear_interpolation
        implicit none
        integer, intent(in) :: x_len, y_len, z_len           
        real(rkind), dimension(x_len), intent(in) :: x_array
        real(rkind), dimension(y_len), intent(in) :: y_array
        real(rkind), dimension(z_len), intent(in) :: z_array
        real(rkind), dimension(x_len, y_len, z_len), intent(in) :: f
        real(rkind), intent(in) :: x,y,z

        real(rkind) :: c00, c10, c01, c11, c0, c1
        real(rkind) :: xd, yd, zd

        real(rkind) :: x1, x2, y1, y2, z1, z2
        integer :: i,j,k
        real(rkind), intent(out) :: val

        call binarysearch(x_len, x_array, x, i)
        call binarysearch(y_len, y_array, y, j)
        call binarysearch(z_len, z_array, z, k)

        if (i == x_len) then
            i = x_len - 1
        end if 
        if (j == y_len) then
            j = y_len - 1
        end if 
        if (k == z_len) then
            k = z_len - 1
        end if 

        x1 = x_array(i)
        x2 = x_array(i+1)

        y1 = y_array(j)
        y2 = y_array(j+1)

        z1 = z_array(k)
        z2 = z_array(k+1)

        xd = (x - x1)/(x2 - x1)
        yd = (y - y1)/(y2 - y1)
        zd = (z - z1)/(z2 - z1)

        c00 = f(i,j,k)     * (one - xd) + f(i+1,j,k)     * xd
        c10 = f(i,j+1,k)   * (one - xd) + f(i+1,j+1,k)   * xd
        c01 = f(i,j,k+1)   * (one - xd) + f(i+1,j,k+1)   * xd
        c11 = f(i,j+1,k+1) * (one - xd) + f(i+1,j+1,k+1) * xd 

        c0  = c00 * (one - yd) + c10 * yd
        c1  = c01 * (one - yd) + c11 * yd

        val = c0  * (one - zd) + c1  * zd

    end subroutine interpolate3d

    pure subroutine bilinearInterp(xF,yF, fF, xC, yC, fC)
        real(rkind), dimension(:,:,:), intent(in) :: xC, yC, fC
        real(rkind), dimension(:,:,:), intent(in) :: xF, yF
        real(rkind), dimension(:,:,:), intent(out) :: fF

        real(rkind), dimension(size(xC,1)) :: xline
        real(rkind), dimension(size(xC,2)) :: yline
        integer :: nxC, nyC, nzC, nxF, nyF, nzF
        integer :: i, j, k


        yline = yC(1,:,1)
        xline = xC(:,1,1)
        fF = 0._rkind

        nxC = size(xC,1)
        nyC = size(xC,2)
        nzC = size(xC,3)

        nxF = size(xF,1)
        nyF = size(xF,2)
        nzF = size(xF,3)

        do k = 1,nzF
            do j = 1,nyF
                do i = 1,nxF
                    fF(i,j,k) = interpolate2d(nxC, xline, nyC, yline, fC(:,:,1), xF(i,j,k), yF(i,j,k))    
                end do
            end do
        end do


    end subroutine

    pure subroutine trilinearInterp(xF, yF, zF, fF, xC, yC, zC, fC)
        real(rkind), dimension(:,:,:), intent(in) :: xC, yC, zC, fC
        real(rkind), dimension(:,:,:), intent(in) :: xF, yF, zF
        real(rkind), dimension(:,:,:), intent(out) :: fF

        real(rkind), dimension(size(xC,1)) :: xline
        real(rkind), dimension(size(xC,2)) :: yline
        real(rkind), dimension(size(xC,3)) :: zline
        integer :: nxC, nyC, nzC, nxF, nyF, nzF
        integer :: i, j, k

        
        yline = yC(1,:,1)
        xline = xC(:,1,1)
        zline = zC(1,1,:)
        fF = 0._rkind

        nxC = size(xC,1)
        nyC = size(xC,2)
        nzC = size(xC,3)

        nxF = size(xF,1)
        nyF = size(xF,2)
        nzF = size(xF,3)

        do k = 1,nzF
            do j = 1,nyF
                do i = 1,nxF
                    call interpolate3d(nxC, xline, nyC, yline, nzC, zline, fC, xF(i,j,k), yF(i,j,k), zF(i,j,k),fF(i,j,k))    
                end do
            end do
        end do

    end subroutine

   subroutine spline (x, y, b, c, d, n)
        !======================================================================
        !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
        !  for cubic spline interpolation
        !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
        !  for  x(i) <= x <= x(i+1)
        !  Alex G: January 2010
        !----------------------------------------------------------------------
        !  input..
        !  x = the arrays of data abscissas (in strictly increasing order)
        !  y = the arrays of data ordinates
        !  n = size of the arrays xi() and yi() (n>=2)
        !  output..
        !  b, c, d  = arrays of spline coefficients
        !  comments ...
        !  spline.f90 program is based on fortran version of program spline.f
        !  the accompanying function fspline can be used for interpolation
        !======================================================================
        implicit none
        integer n
        double precision x(n), y(n), b(n), c(n), d(n)
        integer i, j, gap
        double precision h
        
        gap = n-1
        ! check input
        if ( n < 2 ) return
        if ( n < 3 ) then
          b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
          c(1) = 0.
          d(1) = 0.
          b(2) = b(1)
          c(2) = 0.
          d(2) = 0.
          return
        end if
        !
        ! step 1: preparation
        !
        d(1) = x(2) - x(1)
        c(2) = (y(2) - y(1))/d(1)
        do i = 2, gap
          d(i) = x(i+1) - x(i)
          b(i) = 2.0*(d(i-1) + d(i))
          c(i+1) = (y(i+1) - y(i))/d(i)
          c(i) = c(i+1) - c(i)
        end do
        !
        ! step 2: end conditions 
        !
        b(1) = -d(1)
        b(n) = -d(n-1)
        c(1) = 0.0
        c(n) = 0.0
        if(n /= 3) then
          c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
          c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
          c(1) = c(1)*d(1)**2/(x(4)-x(1))
          c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
        end if
        !
        ! step 3: forward elimination 
        !
        do i = 2, n
          h = d(i-1)/b(i-1)
          b(i) = b(i) - h*d(i-1)
          c(i) = c(i) - h*c(i-1)
        end do
        !
        ! step 4: back substitution
        !
        c(n) = c(n)/b(n)
        do j = 1, gap
          i = n-j
          c(i) = (c(i) - d(i)*c(i+1))/b(i)
        end do
        !
        ! step 5: compute spline coefficients
        !
        b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
        do i = 1, gap
          b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
          d(i) = (c(i+1) - c(i))/d(i)
          c(i) = 3.*c(i)
        end do
        c(n) = 3.0*c(n)
        d(n) = d(n-1)
    end subroutine spline
        
    function ispline(u, x, y, b, c, d, n)
        !======================================================================
        ! function ispline evaluates the cubic spline interpolation at point z
        ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
        ! where  x(i) <= u <= x(i+1)
        !----------------------------------------------------------------------
        ! input..
        ! u       = the abscissa at which the spline is to be evaluated
        ! x, y    = the arrays of given data points
        ! b, c, d = arrays of spline coefficients computed by spline
        ! n       = the number of data points
        ! output:
        ! ispline = interpolated value at point u
        !=======================================================================
        implicit none
        double precision ispline
        integer n
        double precision  u, x(n), y(n), b(n), c(n), d(n)
        integer i, j, k
        double precision dx
        
        ! if u is ouside the x() interval take a boundary value (left or right)
        if(u <= x(1)) then
          ispline = y(1)
          return
        end if
        if(u >= x(n)) then
          ispline = y(n)
          return
        end if
        
        !*
        !  binary search for for i, such that x(i) <= u <= x(i+1)
        !*
        i = 1
        j = n+1
        do while (j > i+1)
          k = (i+j)/2
          if(u < x(k)) then
            j=k
            else
            i=k
           end if
        end do
        !*
        !  evaluate spline interpolation
        !*
        dx = u - x(i)
        ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
    end function ispline


    
end module interpolation
