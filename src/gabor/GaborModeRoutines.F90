module GaborModeRoutines
  use kind_parameters, only: rkind
  use constants,       only: pi
  use fortran_assert,  only: assert
  implicit none
  
  logical :: doWarning = .true.
  real(rkind), parameter :: small = 1.d-10, big = 1.d12

  contains
  include "GaborModeRoutines_files/timeSteppingStuff.F90"
  include "GaborModeRoutines_files/interpolation.F90"
    
    function getModelSpectrum(k,KE,L,nk) result(E) 
      use exits, only: warning
      real(rkind), dimension(:), intent(in) :: k
      integer, intent(in) :: nk
      real(rkind), dimension(nk) :: kL, fL, feta
      real(rkind), intent(in) :: KE, L
      real(rkind) :: C
      real(rkind), dimension(nk) :: E

      ! Non-dimensional wavenumber
      kL = k*L

      ! Model coefficient chosen such that the integrated spectrum equals total
      ! kinetic energy
      C = 0.4843d0

      ! Energetic scalse
      fL = kL**4.d0/(1.d0 + kL*kL)**(17.d0/6.d0)

      ! TODO: Dissipative scales
      feta = 1.d0
      if (doWarning) then
        call warning("WARNING: Finite Re model spectrum is not implemented")
        doWarning = .false.
      end if

      ! Model spectrum
      E = C*2.d0*KE*L*fL*feta
    end function

    pure subroutine normalizeVec(x,y,z)
      real(rkind), dimension(:,:), intent(inout) :: x, y, z
      real(rkind), dimension(size(x,1),size(x,2)) :: mag
    
      mag = sqrt(x*x + y*y + z*z)
      x = x/mag
      y = y/mag
      z = z/mag
    end subroutine

    pure function isOrthogonal(a1,a2,a3,b1,b2,b3) result(TF)
      real(rkind), dimension(:), intent(in) :: a1, a2, a3, b1, b2, b3
      logical :: TF
      real(rkind) :: maxDiv, small

      small = 1.d-14
    
      maxDiv = maxval(a1*b1 + a2*b2 + a3*b3)
      TF = .false.
      if (maxDiv < small) TF = .true.
    end function
      
    pure function getNyquist(L,n) result(kNyq)
      real(rkind), intent(in) :: L
      integer, intent(in) :: n
      real(rkind) :: kNyq

      kNyq = real(n,rkind)*0.5d0*(2.d0*pi/L)
    end function
   
    subroutine computeKminKmax(Lx,Ly,Lz,nxLS,nyLS,nzLS,nxSS,nySS,nzSS,kmin,kmax)
      real(rkind), intent(in) :: Lx, Ly, Lz 
      integer, intent(in) :: nxLS, nyLS, nzLS, nxSS, nySS, nzSS
      real(rkind), intent(out) :: kmin, kmax
      real(rkind) :: kxNyqLES, kyNyqLES, kzNyqLES
      real(rkind) :: kxNyqF, kyNyqF, kzNyqF

      kxNyqLES = getNyquist(Lx,nxLS)
      kyNyqLES = getNyquist(Ly,nyLS)
      kzNyqLES = getNyquist(Lz,nzLS)

      kxNyqF = getNyquist(Lx,nxSS)
      kyNyqF = getNyquist(Ly,nySS)
      kzNyqF = getNyquist(Lz,nzSS)

      kmin = minval([kxNyqLES, kyNyqLES, kzNyqLES])
      kmax = minval([kxNyqF,   kyNyqF,   kzNyqF])
    end subroutine

    SUBROUTINE PFQ(X_ACT,Z_ACT)
      use hygfx_mod, only: HYGFX
      REAL(kind=8) :: A, B, C, X, HF
      REAL(kind=8) :: A_ACT, B_ACT, C_ACT
      REAL(rkind), intent(in) :: X_ACT
      real(rkind), intent(out) :: Z_ACT
      A_ACT = 0.333333d0
      B_ACT = 2.833333d0
      C_ACT = 1.333333d0
      A = A_ACT
      B = B_ACT
      C = C_ACT
      X = real(X_ACT,kind=8)
      
      CALL HYGFX(A,B,C,X,HF)
      Z_ACT = real(HF,rkind)
      RETURN
    END SUBROUTINE PFQ 

end module
