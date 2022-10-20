module hygfx_mod

contains
    subroutine psi ( x, ps )
    
    !*****************************************************************************80
    !
    !! PSI computes the PSI function.
    !
    !  Licensing:
    !
    !    The original FORTRAN77 version of this routine is copyrighted by 
    !    Shanjie Zhang and Jianming Jin.  However, they give permission to 
    !    incorporate this routine into a user program that the copyright 
    !    is acknowledged.
    !
    !  Modified:
    !
    !    08 September 2007
    !
    !  Author:
    !
    !    Original FORTRAN77 by Shanjie Zhang, Jianming Jin.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Shanjie Zhang, Jianming Jin,
    !    Computation of Special Functions,
    !    Wiley, 1996,
    !    ISBN: 0-471-11963-6,
    !    LC: QA351.C45
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument.
    !
    !    Output, real ( kind = 8 ) PS, the value of the PSI function.
    !
      implicit none
    
      real ( kind = 8 ), parameter :: a1 = -0.83333333333333333D-01
      real ( kind = 8 ), parameter :: a2 =  0.83333333333333333D-02
      real ( kind = 8 ), parameter :: a3 = -0.39682539682539683D-02
      real ( kind = 8 ), parameter :: a4 =  0.41666666666666667D-02
      real ( kind = 8 ), parameter :: a5 = -0.75757575757575758D-02
      real ( kind = 8 ), parameter :: a6 =  0.21092796092796093D-01
      real ( kind = 8 ), parameter :: a7 = -0.83333333333333333D-01
      real ( kind = 8 ), parameter :: a8 =  0.4432598039215686D+00
      real ( kind = 8 ), parameter :: el = 0.5772156649015329D+00
      integer ( kind = 4 ) k
      integer ( kind = 4 ) n
      real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
      real ( kind = 8 ) ps
      real ( kind = 8 ) s
      real ( kind = 8 ) x
      real ( kind = 8 ) x2
      real ( kind = 8 ) xa
    
      xa = abs ( x )
      s = 0.0D+00
    
      if ( x == aint ( x ) .and. x <= 0.0D+00 ) then
    
        ps = 1.0D+300
        return
    
      else if ( xa == aint ( xa ) ) then
    
        n = int ( xa )
        do k = 1, n - 1
          s = s + 1.0D+00 / real ( k, kind = 8 )
        end do
    
        ps = - el + s
    
      else if ( xa + 0.5D+00 == aint ( xa + 0.5D+00 ) ) then
    
        n = int ( xa - 0.5D+00 )
    
        do k = 1, n
          s = s + 1.0D+00 / real ( 2 * k - 1, kind = 8 )
        end do
    
        ps = - el + 2.0D+00 * s - 1.386294361119891D+00
    
      else
    
        if ( xa < 10.0D+00 ) then
    
          n = 10 - int ( xa )
          do k = 0, n - 1
            s = s + 1.0D+00 / ( xa + real ( k, kind = 8 ) )
          end do
    
          xa = xa + real ( n, kind = 8 )
    
        end if
    
        x2 = 1.0D+00 / ( xa * xa )
    
        ps = log ( xa ) - 0.5D+00 / xa + x2 * ((((((( &
                 a8   &
          * x2 + a7 ) &
          * x2 + a6 ) &
          * x2 + a5 ) &
          * x2 + a4 ) &
          * x2 + a3 ) &
          * x2 + a2 ) &
          * x2 + a1 )
    
        ps = ps - s
    
      end if
    
      if ( x < 0.0D+00 ) then
        ps = ps - pi * cos ( pi * x ) / sin ( pi * x ) - 1.0D+00 / x
      end if
    
      return
    end subroutine
    
    subroutine gamma ( x, ga )
    
    !*****************************************************************************80
    !
    !! GAMMA evaluates the Gamma function.
    !
    !  Licensing:
    !
    !    The original FORTRAN77 version of this routine is copyrighted by 
    !    Shanjie Zhang and Jianming Jin.  However, they give permission to 
    !    incorporate this routine into a user program that the copyright 
    !    is acknowledged.
    !
    !  Modified:
    !
    !    08 September 2007
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Shanjie Zhang, Jianming Jin,
    !    Computation of Special Functions,
    !    Wiley, 1996,
    !    ISBN: 0-471-11963-6,
    !    LC: QA351.C45
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument.
    !    X must not be 0, or any negative integer.
    !
    !    Output, real ( kind = 8 ) GA, the value of the Gamma function.
    !
      implicit none
    
      real ( kind = 8 ), dimension ( 26 ) :: g = (/ &
        1.0D+00, &
        0.5772156649015329D+00, &
       -0.6558780715202538D+00, &
       -0.420026350340952D-01, &
        0.1665386113822915D+00, &
       -0.421977345555443D-01, &
       -0.96219715278770D-02, &
        0.72189432466630D-02, &
       -0.11651675918591D-02, &
       -0.2152416741149D-03, &
        0.1280502823882D-03, & 
       -0.201348547807D-04, &
       -0.12504934821D-05, &
        0.11330272320D-05, &
       -0.2056338417D-06, & 
        0.61160950D-08, &
        0.50020075D-08, &
       -0.11812746D-08, &
        0.1043427D-09, & 
        0.77823D-11, &
       -0.36968D-11, &
        0.51D-12, &
       -0.206D-13, &
       -0.54D-14, &
        0.14D-14, &
        0.1D-15 /)
      real ( kind = 8 ) ga
      real ( kind = 8 ) gr
      integer ( kind = 4 ) k
      integer ( kind = 4 ) m
      integer ( kind = 4 ) m1
      real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
      real ( kind = 8 ) r
      real ( kind = 8 ) x
      real ( kind = 8 ) z
    
      if ( x == aint ( x ) ) then
    
        if ( 0.0D+00 < x ) then
          ga = 1.0D+00
          m1 = int ( x ) - 1
          do k = 2, m1
            ga = ga * k
          end do
        else
          ga = 1.0D+300
        end if
    
      else
    
        if ( 1.0D+00 < abs ( x ) ) then
          z = abs ( x )
          m = int ( z )
          r = 1.0D+00
          do k = 1, m
            r = r * ( z - real ( k, kind = 8 ) )
          end do
          z = z - real ( m, kind = 8 )
        else
          z = x
        end if
    
        gr = g(26)
        do k = 25, 1, -1
          gr = gr * z + g(k)
        end do
    
        ga = 1.0D+00 / ( gr * z )
    
        if ( 1.0D+00 < abs ( x ) ) then
          ga = ga * r
          if ( x < 0.0D+00 ) then
            ga = - pi / ( x* ga * sin ( pi * x ) )
          end if
        end if
    
      end if
    
      return
    end

    subroutine hygfx ( a, b, c, x, hf )
    
    !*****************************************************************************80
    !
    !! HYGFX evaluates the hypergeometric function F(A,B,C,X).
    !
    !  Licensing:
    !
    !    The original FORTRAN77 version of this routine is copyrighted by 
    !    Shanjie Zhang and Jianming Jin.  However, they give permission to 
    !    incorporate this routine into a user program that the copyright 
    !    is acknowledged.
    !
    !  Modified:
    !
    !    08 September 2007
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Shanjie Zhang, Jianming Jin,
    !    Computation of Special Functions,
    !    Wiley, 1996,
    !    ISBN: 0-471-11963-6,
    !    LC: QA351.C45
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, C, X, the arguments of the function.
    !    C must not be equal to a nonpositive integer.
    !    X < 1.
    !
    !    Output, real HF, the value of the function.
    !
      implicit none
    
      real ( kind = 8 ) a
      real ( kind = 8 ) a0
      real ( kind = 8 ) aa
      real ( kind = 8 ) b
      real ( kind = 8 ) bb
      real ( kind = 8 ) c
      real ( kind = 8 ) c0
      real ( kind = 8 ) c1
      real ( kind = 8 ), parameter :: el = 0.5772156649015329D+00
      real ( kind = 8 ) eps
      real ( kind = 8 ) f0
      real ( kind = 8 ) f1
      real ( kind = 8 ) g0
      real ( kind = 8 ) g1
      real ( kind = 8 ) g2
      real ( kind = 8 ) g3
      real ( kind = 8 ) ga
      real ( kind = 8 ) gabc
      real ( kind = 8 ) gam
      real ( kind = 8 ) gb
      real ( kind = 8 ) gbm
      real ( kind = 8 ) gc
      real ( kind = 8 ) gca
      real ( kind = 8 ) gcab
      real ( kind = 8 ) gcb
      real ( kind = 8 ) gm
      real ( kind = 8 ) hf
      real ( kind = 8 ) hw
      integer ( kind = 4 ) j
      integer ( kind = 4 ) k
      logical l0
      logical l1
      logical l2
      logical l3
      logical l4
      logical l5
      integer ( kind = 4 ) m
      integer ( kind = 4 ) nm
      real ( kind = 8 ) pa
      real ( kind = 8 ) pb
      real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
      real ( kind = 8 ) r
      real ( kind = 8 ) r0
      real ( kind = 8 ) r1
      real ( kind = 8 ) rm
      real ( kind = 8 ) rp
      real ( kind = 8 ) sm
      real ( kind = 8 ) sp
      real ( kind = 8 ) sp0
      real ( kind = 8 ) x
      real ( kind = 8 ) x1
    
      l0 = ( c == aint ( c ) ) .and. ( c < 0.0D+00 )
      l1 = ( 1.0D+00 - x < 1.0D-15 ) .and. ( c - a - b <= 0.0D+00 )
      l2 = ( a == aint ( a ) ) .and. ( a < 0.0D+00 )
      l3 = ( b == aint ( b ) ) .and. ( b < 0.0D+00 )
      l4 = ( c - a == aint ( c - a ) ) .and. ( c - a <= 0.0D+00 )
      l5 = ( c - b == aint ( c - b ) ) .and. ( c - b <= 0.0D+00 )
    
      if ( l0 .or. l1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HYGFX - Fatal error!'
        write ( *, '(a)' ) '  The hypergeometric series is divergent.'
        return
      end if
    
      if ( 0.95D+00 < x ) then 
        eps = 1.0D-08
      else
        eps = 1.0D-15
      end if
    
      if ( x == 0.0D+00 .or. a == 0.0D+00 .or. b == 0.0D+00 ) then
    
        hf = 1.0D+00
        return
    
      else if ( 1.0D+00 - x == eps .and. 0.0D+00 < c - a - b ) then
    
        call gamma ( c, gc )
        call gamma ( c - a - b, gcab )
        call gamma ( c - a, gca )
        call gamma ( c - b, gcb )
        hf = gc * gcab /( gca *gcb )
        return
    
      else if ( 1.0D+00 + x <= eps .and. abs ( c - a + b - 1.0D+00 ) <= eps ) then
    
        g0 = sqrt ( pi ) * 2.0D+00**( - a )
        call gamma ( c, g1 )
        call gamma ( 1.0D+00 + a / 2.0D+00 - b, g2 )
        call gamma ( 0.5D+00 + 0.5D+00 * a, g3 )
        hf = g0 * g1 / ( g2 * g3 )
        return
    
      else if ( l2 .or. l3 ) then
    
        if ( l2 ) then
          nm = int ( abs ( a ) )
        end if
    
        if ( l3 ) then
          nm = int ( abs ( b ) )
        end if
    
        hf = 1.0D+00
        r = 1.0D+00
    
        do k = 1, nm
          r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( c + k - 1.0D+00 ) ) * x
          hf = hf + r
        end do
    
        return
    
      else if ( l4 .or. l5 ) then
    
        if ( l4 ) then
          nm = int ( abs ( c - a ) )
        end if
    
        if ( l5 ) then
          nm = int ( abs ( c - b ) )
        end if
    
        hf = 1.0D+00
        r  = 1.0D+00
        do k = 1, nm
          r = r * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
            / ( k * ( c + k - 1.0D+00 ) ) * x
          hf = hf + r
        end do
        hf = ( 1.0D+00 - x )**( c - a - b ) * hf
        return
    
      end if
    
      aa = a
      bb = b
      x1 = x
    !
    !  WARNING: ALTERATION OF INPUT ARGUMENTS A AND B, WHICH MIGHT BE CONSTANTS.
    !
      if ( x < 0.0D+00 ) then
        x = x / ( x - 1.0D+00 )
        if ( a < c .and. b < a .and. 0.0D+00 < b ) then
          a = bb
          b = aa
        end if
        b = c - b
      end if
    
      if ( 0.75D+00 <= x ) then
    
        gm = 0.0D+00
    
        if ( abs ( c - a - b - aint ( c - a - b ) ) < 1.0D-15 ) then
    
          m = int ( c - a - b )
          call gamma ( a, ga )
          call gamma ( b, gb )
          call gamma ( c, gc )
          call gamma ( a + m, gam )
          call gamma ( b + m, gbm )
          call psi ( a, pa )
          call psi ( b, pb )
    
          if ( m /= 0 ) then
            gm = 1.0D+00
          end if
    
          do j = 1, abs ( m ) - 1
            gm = gm * j
          end do
    
          rm = 1.0D+00
          do j = 1, abs ( m )
            rm = rm * j
          end do
    
          f0 = 1.0D+00
          r0 = 1.0D+00
          r1 = 1.0D+00
          sp0 = 0.0D+00
          sp = 0.0D+00
    
          if ( 0 <= m ) then
    
            c0 = gm * gc / ( gam * gbm )
            c1 = - gc * ( x - 1.0D+00 )**m / ( ga * gb * rm )
    
            do k = 1, m - 1
              r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
                / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
            end do
    
            do k = 1, m
              sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) &
                + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / real ( k, kind = 8 )
            end do
    
            f1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
            hw = f1
    
            do k = 1, 250
    
              sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) ) &
                + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )
    
              sm = 0.0D+00
              do j = 1, m
                sm = sm + ( 1.0D+00 - a ) &
                  / ( ( j + k ) * ( a + j + k - 1.0D+00 ) ) &
                  + 1.0D+00 / ( b + j + k - 1.0D+00 )
              end do
    
              rp = pa + pb + 2.0D+00 * el + sp + sm + log ( 1.0D+00 - x )
    
              r1 = r1 * ( a + m + k - 1.0D+00 ) * ( b + m + k - 1.0D+00 ) &
                / ( k * ( m + k ) ) * ( 1.0D+00 - x )
    
              f1 = f1 + r1 * rp
    
              if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
                exit
              end if
    
              hw = f1
    
            end do
    
            hf = f0 * c0 + f1 * c1
    
          else if ( m < 0 ) then
    
            m = - m
            c0 = gm * gc / ( ga * gb * ( 1.0D+00 - x )**m )
            c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )
    
            do k = 1, m - 1
              r0 = r0 * ( a - m + k - 1.0D+00 ) * ( b - m + k - 1.0D+00 ) &
                / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
            end do
    
            do k = 1, m
              sp0 = sp0 + 1.0D+00 / real ( k, kind = 8 )
            end do
    
            f1 = pa + pb - sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
    
            do k = 1, 250
    
              sp = sp + ( 1.0D+00 - a ) &
                / ( k * ( a + k - 1.0D+00 ) ) &
                + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )
    
              sm = 0.0D+00
              do j = 1, m
                sm = sm + 1.0D+00 / real ( j + k, kind = 8 )
              end do
    
              rp = pa + pb + 2.0D+00 * el + sp - sm + log ( 1.0D+00 - x )
    
              r1 = r1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
                / ( k * ( m + k ) ) * ( 1.0D+00 - x )
    
              f1 = f1 + r1 * rp
    
              if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
                exit
              end if
    
              hw = f1
    
            end do
    
            hf = f0 * c0 + f1 * c1
    
          end if
    
        else
    
          call gamma ( a, ga )
          call gamma ( b, gb )
          call gamma ( c, gc )
          call gamma ( c - a, gca )
          call gamma ( c - b, gcb )
          call gamma ( c - a - b, gcab )
          call gamma ( a + b - c, gabc )
          c0 = gc * gcab / ( gca * gcb )
          c1 = gc * gabc / ( ga * gb ) * ( 1.0D+00 - x )**( c - a - b )
          hf = 0.0D+00
          r0 = c0
          r1 = c1
    
          hw = 0.d0 
          do k = 1, 250
    
            r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
              / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - x )
    
            r1 = r1 * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
              / ( k * ( c - a - b + k ) ) * ( 1.0D+00 - x )
    
            hf = hf + r0 + r1
            
            if ( abs ( hf - hw ) < abs ( hf ) * eps ) then
              exit
            end if
    
            hw = hf
    
          end do
    
          hf = hf + c0 + c1
    
        end if
    
      else
    
        a0 = 1.0D+00
    
        if ( a < c .and. c < 2.0D+00 * a .and. b < c .and. c < 2.0D+00 * b ) then
    
          a0 = ( 1.0D+00 - x )**( c - a - b )
          a = c - a
          b = c - b
    
        end if
    
        hf = 1.0D+00
        r = 1.0D+00
    
        hw = 0.d0 
        do k = 1, 250
    
          r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( c + k - 1.0D+00 ) ) * x
    
          hf = hf + r
    
          if ( abs ( hf - hw ) <= abs ( hf ) * eps ) then
            exit
          end if
    
          hw = hf
    
        end do
    
        hf = a0 * hf
    
      end if
    
      if ( x1 < 0.0D+00 ) then
        x = x1
        c0 = 1.0D+00 / ( 1.0D+00 - x )**aa
        hf = c0 * hf
      end if
    
      a = aa
      b = bb
    
      if ( 120 < k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HYGFX - Warning!'
        write ( *, '(a)' ) '  A large number of iterations were needed.'
        write ( *, '(a)' ) '  The accuracy of the results should be checked.'
      end if
    
      return
    end subroutine

end module 
