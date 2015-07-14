      SUBROUTINE SolvePenta(N,E,A,D,C,F,B,X)

c   RESULTS:  matrix has 5 bands, EADCF, with D being the main diagonal,
c   E and A are the lower diagonals, and C and F are the upper diagonals.

c     E is defined for rows i = 3:N, but is defined as E(1) to E(N-2)
c     A is defined for rows i = 2:N, but is defined as A(1) to A(N-1)
c     D is defined for rows i = 1:N
c     C is defined for rows i = 1:N-1, but the last element isn't used
c     F is defined for rows i = 1:N-2, but the last 2 elements aren't used

c   B is the right-hand side
c   X is the solution vector

      implicit none
      integer i,n
      double precision E(N),A(N),D(N),C(N),F(N),B(N),X(N),XMULT
      DO 2 I = 2,N-1
        XMULT = A(I-1)/D(I-1)
        D(I) = D(I) - XMULT*C(I-1)
        C(I) = C(I) - XMULT*F(I-1)
        B(I) = B(I) - XMULT*B(I-1)
        XMULT = E(I-1)/D(I-1)
        A(I) = A(I) - XMULT*C(I-1)
        D(I+1) = D(I+1) - XMULT*F(I-1)
        B(I+1) = B(I+1) - XMULT*B(I-1)
   2  CONTINUE
      XMULT = A(N-1)/D(N-1)
      D(N) = D(N) - XMULT*C(N-1)
      X(N) = (B(N) - XMULT*B(N-1))/D(N)
      X(N-1) = (B(N-1) - C(N-1)*X(N))/D(N-1)
      DO 3 I = N-2,1,-1
        X(I) = (B(I) - F(I)*X(I+2) - C(I)*X(I+1))/D(I)
   3  CONTINUE
      RETURN
      END
