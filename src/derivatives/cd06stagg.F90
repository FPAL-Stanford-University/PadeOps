! Routines specific to 6th order (stagerred) Compact Finite Differencing scheme

module cd06staggstuff

    use kind_parameters, only: rkind
    use constants,       only : zero,one,two, three, four, ten
    use exits,           only : GracefulExit
    implicit none

    private
    public :: cd06stagg
  

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic 1st derivative evaluation !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
    ! Set the scheme for the edge nodes (Ref. for notation: Lele - JCP paper)
    ! 1st derivative 
    real(rkind), parameter                   :: alpha =  3._rkind
    real(rkind), parameter                   :: p     = 17._rkind / 6._rkind
    real(rkind), parameter                   :: q     =  3._rkind / 2._rkind
    real(rkind), parameter                   :: r     =  3._rkind / 2._rkind
    real(rkind), parameter                   :: s     = -1._rkind / 6._rkind 


    !!  Calculate the corressponding weights 
    ! Step 1: Assign the interior scheme
    real(rkind), parameter                   :: alpha06d1   =  1.0_rkind / 3.0_rkind
    real(rkind), parameter                   :: a06d1       = (14.0_rkind / 9.0_rkind) / 2.0_rkind
    real(rkind), parameter                   :: b06d1       = ( 1.0_rkind / 9.0_rkind) / 4.0_rkind
    real(rkind), parameter                   :: qhat        = a06d1
    real(rkind), parameter                   :: rhat        = b06d1
    real(rkind), parameter                   :: alpha_hat   = alpha06d1

    ! Step 2: Assign the scheme at node 2 to be Standard Pade (4th Order)
    real(rkind), parameter                   :: q_p          = 3._rkind/4._rkind
    real(rkind), parameter                   :: alpha_p     = 1._rkind/4._rkind
    
    ! Step 3: Get the scheme at the node 3
    real(rkind), parameter                   :: alpha_pp = ((40*alpha_hat - 1)*q  + 7*(4*alpha_hat &
                                                         -  1)*s)/(16*(alpha_hat + 2)*q + 8*(1      &
                                                         -  4*alpha_hat)*s)
    real(rkind), parameter                   :: q_pp     = (1._rkind/3._rkind)*(alpha_pp + 2)
    real(rkind), parameter                   :: r_pp     = (1._rkind/12._rkind)*(4*alpha_pp - 1)
    real(rkind), parameter                   :: s_pp     =  0._rkind

    ! Step 4: Get the weights
    real(rkind), parameter                   :: w1 = (2*alpha_hat + 1)/(2*(q + s))
    real(rkind), parameter                   :: w2 = ((8*alpha_hat + 7)*q - 6*(2*alpha_hat + 1)*r &
                                                   + (8*alpha_hat + 7)*s)/(9*(q + s))
    real(rkind), parameter                   :: w3 = (4*(alpha_hat + 2)*q + 2*(1 - 4*alpha_hat)*s) &
                                                   / (9*(q + s))

    ! Step 5: Weights for stagerred sided scheme
    real(rkind), parameter                   :: w0s = 223.d0/186.d0, w1s = 61.d0/62.d0    


    type cd06stagg

        private

        logical     :: AmIPeriodic = .false. 
        integer     :: n
        integer     :: nE
        real(rkind) :: dx
        real(rkind) :: onebydx
        real(rkind) :: onebydx2

        logical     :: isTopEven                          ! Boundary condition type. 
        logical     :: isBotEven                          ! Boundary condition type. 

        logical     :: isTopSided = .FALSE.
        logical     :: isBotSided = .FALSE.

        real(rkind), allocatable, dimension(:,:) :: TriD1_E2C
        real(rkind), allocatable, dimension(:,:) :: TriD1_C2E
        real(rkind), allocatable, dimension(:,:) :: TriD1_C2C
        real(rkind), allocatable, dimension(:,:) :: TriD1_E2E
        real(rkind), allocatable, dimension(:,:) :: TriD2_C2C
        real(rkind), allocatable, dimension(:,:) :: TriD2_E2E
        real(rkind), allocatable, dimension(:,:) :: TriInterp_E2C
        real(rkind), allocatable, dimension(:,:) :: TriInterp_C2E
        real(rkind), allocatable, dimension(:,:) :: LU_D1
        real(rkind), allocatable, dimension(:,:) :: LU_D2
        real(rkind), allocatable, dimension(:,:) :: LU_interp

        contains

        procedure, private :: init_nonperiodic
        procedure, private :: init_periodic
        generic            :: init => init_periodic, init_nonperiodic
        
        procedure, private :: ComputeD1RHS_E2C_REAL 
        procedure, private :: ComputeD1RHS_C2E_REAL
        procedure, private :: ComputeD1RHS_E2E_REAL
        procedure, private :: ComputeD1RHS_C2C_REAL
        procedure, private :: ComputeD2RHS_E2E_REAL
        procedure, private :: ComputeD2RHS_C2C_REAL
        procedure, private :: ComputeD1RHS_E2C_CMPLX 
        procedure, private :: ComputeD1RHS_C2E_CMPLX
        procedure, private :: ComputeD1RHS_E2E_CMPLX
        procedure, private :: ComputeD1RHS_C2C_CMPLX
        procedure, private :: ComputeD2RHS_E2E_CMPLX
        procedure, private :: ComputeD2RHS_C2C_CMPLX

        procedure, private :: ComputeInterpRHS_C2E_REAL
        procedure, private :: ComputeInterpRHS_C2E_CMPLX
        procedure, private :: ComputeInterpRHS_E2C_REAL
        procedure, private :: ComputeInterpRHS_E2C_CMPLX
        
        procedure, private :: ComputeTriD1_E2C
        procedure, private :: ComputeTriD1_C2E
        procedure, private :: ComputeTriD1_E2E
        procedure, private :: ComputeTriD1_C2C
        procedure, private :: ComputeTriD2_E2E
        procedure, private :: ComputeTriD2_C2C

        procedure, private :: ComputeZD1RHS_E2C_REAL_periodic
        procedure, private :: ComputeZD1RHS_E2C_CMPLX_periodic
        procedure, private :: ComputeZD1RHS_C2E_REAL_periodic
        procedure, private :: ComputeZD1RHS_C2E_CMPLX_periodic
        procedure, private :: ComputeZInterpRHS_E2C_REAL_periodic
        procedure, private :: ComputeZInterpRHS_E2C_CMPLX_periodic
        procedure, private :: ComputeZInterpRHS_C2E_REAL_periodic
        procedure, private :: ComputeZInterpRHS_C2E_CMPLX_periodic
        procedure, private :: SolveZLU_REAL
        procedure, private :: SolveZLU_CMPLX

        procedure, private :: ComputeTriInterp_E2C
        procedure, private :: ComputeTriInterp_C2E
        
        procedure, private :: ddz_E2C_REAL 
        procedure, private :: ddz_C2E_REAL
        procedure, private :: ddz_C2C_REAL
        procedure, private :: ddz_E2E_REAL
        procedure, private :: d2dz2_C2C_REAL
        procedure, private :: d2dz2_E2E_REAL
        procedure, private :: ddz_E2C_CMPLX 
        procedure, private :: ddz_C2E_CMPLX
        procedure, private :: ddz_C2C_CMPLX
        procedure, private :: ddz_E2E_CMPLX
        procedure, private :: d2dz2_C2C_CMPLX
        procedure, private :: d2dz2_E2E_CMPLX
        procedure          :: destroy
        
        generic :: ddz_E2C => ddz_E2C_REAL, ddz_E2C_CMPLX
        generic :: ddz_C2E => ddz_C2E_REAL, ddz_C2E_CMPLX
        generic :: ddz_E2E => ddz_E2E_REAL, ddz_E2E_CMPLX
        generic :: ddz_C2C => ddz_C2C_REAL, ddz_C2C_CMPLX
        generic :: d2dz2_E2E => d2dz2_E2E_REAL, d2dz2_E2E_CMPLX
        generic :: d2dz2_C2C => d2dz2_C2C_REAL, d2dz2_C2C_CMPLX

        procedure, private :: InterpZ_E2C_REAL 
        procedure, private :: InterpZ_C2E_REAL
        procedure, private :: InterpZ_E2C_CMPLX 
        procedure, private :: InterpZ_C2E_CMPLX

        generic :: InterpZ_E2C => InterpZ_E2C_REAL, InterpZ_E2C_CMPLX
        generic :: InterpZ_C2E => InterpZ_C2E_REAL, InterpZ_C2E_CMPLX

    end type

contains

    subroutine init_periodic(this, nx, dx)
        class( cd06stagg ), intent(inout) :: this
        integer, intent(in) :: nx
        real(rkind), intent(in) :: dx
        real(rkind), parameter :: alpha06stagg_d1     = 9._rkind/62._rkind 
        real(rkind), parameter :: alpha06stagg_d2     = 2._rkind/11._rkind
        real(rkind), parameter :: alpha06stagg_interp = 3._rkind/10._rkind

        this%n = nx; this%nE = nx + 1
        this%dx = dx; this%onebydx = one/dx; this%onebydx2 = this%onebydx/dx
        this%AmIPeriodic = .true. 
        
        if (nx .LE. 4) then
            call GracefulExit("CD06_stagg requires at least 4 points",21)
        end if  
         
        if(allocated( this%LU_D1 )) deallocate( this%LU_D1 ); allocate( this%LU_D1(this%n,5) )
        call ComputeLU(this%LU_D1,nx,alpha06stagg_d1,one,alpha06stagg_d1)

        if(allocated( this%LU_D2 )) deallocate( this%LU_D2 ); allocate( this%LU_D2(this%n,5) )
        call ComputeLU(this%LU_D2,nx,alpha06stagg_d2,one,alpha06stagg_d2)

        if(allocated( this%LU_interp )) deallocate( this%LU_interp ); allocate( this%LU_interp(this%n,5) )
        call ComputeLU(this%LU_interp,nx,alpha06stagg_interp,one,alpha06stagg_interp)

    end subroutine 

    subroutine init_nonperiodic(this, nx, dx, isTopEven, isBotEven, isTopSided, isBotSided) 
        class( cd06stagg ), intent(inout) :: this
        integer, intent(in) :: nx
        real(rkind), intent(in) :: dx
        logical, intent(in) :: isTopEven, isBotEven
        logical, intent(in), optional :: isTopSided, isBotSided
    
        this%n = nx; this%nE = nx + 1
        this%dx = dx; this%onebydx = one/dx; this%onebydx2 = this%onebydx/dx
    
        this%isTopEven = isTopEven; this%isBotEven = isBotEven
 
        if (present(isTopSided)) then
            this%isTopSided = isTopSided
        end if 

        if (present(isBotSided)) then
            this%isBotSided = isBotSided
        end if 

        if (nx .LE. 4) then
            call GracefulExit("CD06_stagg requires at least 4 points",21)
        end if  
        allocate( this%TriD1_E2C(nx  ,3) ); allocate( this%TriD1_C2E(nx+1,3) )
        allocate( this%TriD1_E2E(nx+1,3) ); allocate( this%TriD1_C2C(nx  ,3) )
        allocate( this%TriD2_E2E(nx+1,3) ); allocate( this%TriD2_C2C(nx  ,3) )
  
        allocate( this%TriInterp_E2C(nx  ,3) ); allocate( this%TriInterp_C2E(nx+1,3) )
        
        call this%ComputeTriD1_E2C(); call this%ComputeTriD1_C2E()
        call this%ComputeTriD1_E2E(); call this%ComputeTriD1_C2C()
        call this%ComputeTriD2_E2E(); call this%ComputeTriD2_C2C()
    
        call this%ComputeTriInterp_E2C(); call this%ComputeTriInterp_C2E()
        this%AmIPeriodic = .false. 
    end subroutine

    subroutine destroy(this)
        class( cd06stagg ), intent(inout) :: this
        if (this%AmIperiodic) then
            deallocate(this%LU_D1, this%LU_D2, this%LU_interp)
        else
            deallocate( this%TriD1_E2C, this%TriD1_C2E, this%TriD1_E2E, this%TriD1_C2C )
            deallocate( this%TriInterp_E2C, this%TriInterp_C2E)
            deallocate( this%TriD2_E2E, this%TriD2_C2C)
        end if 
    end subroutine
    
#include "STAGG_CD06_files/ComputeTri_allRoutines.F90"
#include "STAGG_CD06_files/TridiagSolver_allRoutines.F90"

    subroutine SolveZLU_REAL(this,y,n1,n2,LU)
        
        class (cd06stagg), intent(in) :: this
        integer, intent(in) :: n1,n2
        real(rkind), dimension(this%n,5), intent(in) :: LU
        real(rkind), dimension(n1,n2,this%n), intent(inout) :: y  ! Take in RHS and put solution into it
        integer ::  k
        real(rkind), dimension(n1,n2) :: sum1 

        ! Step 2
        sum1 = LU(1,2)*y(:,:,1)
        do k = 2,this%n-1
            y(:,:,k) = y(:,:,k) - LU(k,1)*y(:,:,k-1)
            sum1 = sum1 + LU(k,2)*y(:,:,k)
        end do
        y(:,:,this%n) = y(:,:,this%n) - sum1
    
        ! Step 3
        y(:,:,this%n)   = y(:,:,this%n) * LU(this%n,3)

        y(:,:,this%n-1) =  y(:,:,this%n-1) * LU(this%n-1,3) - y(:,:,this%n) * LU(this%n-1,5) 
        do k = this%n-2,1,-1
            y(:,:,k) =  y(:,:,k) * LU(k,3)- y(:,:,k+1) * LU(k,4)- y(:,:,this%n) * LU(k,5)
        end do
    
    end subroutine
      
    subroutine SolveZLU_CMPLX(this,y,n1,n2,LU)
        
        class (cd06stagg), intent(in) :: this
        integer, intent(in) :: n1,n2
        real(rkind), dimension(this%n,5), intent(in) :: LU
        complex(rkind), dimension(n1,n2,this%n), intent(inout) :: y  ! Take in RHS and put solution into it
        integer ::  k
        complex(rkind), dimension(n1,n2) :: sum1 

        ! Step 2
        sum1 = LU(1,2)*y(:,:,1)
        do k = 2,this%n-1
            y(:,:,k) = y(:,:,k) - LU(k,1)*y(:,:,k-1)
            sum1 = sum1 + LU(k,2)*y(:,:,k)
        end do
        y(:,:,this%n) = y(:,:,this%n) - sum1  
        ! Step 3
        y(:,:,this%n)   = y(:,:,this%n) * LU(this%n,3)

        y(:,:,this%n-1) =  y(:,:,this%n-1) * LU(this%n-1,3) - y(:,:,this%n) * LU(this%n-1,5) 
        do k = this%n-2,1,-1
            y(:,:,k) =  y(:,:,k) * LU(k,3)- y(:,:,k+1) * LU(k,4)- y(:,:,this%n) * LU(k,5)
        end do
    
    end subroutine
    
    pure subroutine ComputeZD1RHS_E2C_REAL_periodic(this,f, RHS, n1, n2) 
         
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%nE), intent(in)  :: f
        real(rkind), dimension(n1,n2,this%n) , intent(out) :: RHS
        real(rkind), parameter :: a = (63._rkind/62._rkind) , b = (17._rkind/62._rkind)/three
        real(rkind) :: a06, b06
        
        a06 = a * this%onebydx
        b06 = b * this%onebydx
        
        RHS(:,:,1         ) = a06 * ( f(:,:,2)          - f(:,:,1       ) ) &
                            + b06 * ( f(:,:,3)          - f(:,:,this%n  ) ) 
        
        RHS(:,:,2:this%n-1) = a06 * ( f(:,:,3:this%n  ) - f(:,:,2:this%n-1) ) &
                            + b06 * ( f(:,:,4:this%n+1) - f(:,:,1:this%n-2) ) 
        
        RHS(:,:,this%n    ) = a06 * ( f(:,:,this%n+1)   - f(:,:,this%n  ) ) &
                            + b06 * ( f(:,:,2)          - f(:,:,this%n-1) )

   end subroutine  

   pure subroutine ComputeZD1RHS_E2C_CMPLX_periodic(this,f, RHS, n1, n2) 
         
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%nE), intent(in) :: f
        complex(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
        real(rkind), parameter :: a = (63._rkind/62._rkind)/one , b = (17._rkind/62._rkind)/three
        real(rkind) :: a06, b06

        a06 = a * this%onebydx
        b06 = b * this%onebydx
        
        RHS(:,:,1         ) = a06 * ( f(:,:,2)          - f(:,:,1       ) ) &
                            + b06 * ( f(:,:,3)          - f(:,:,this%n  ) ) 
        
        RHS(:,:,2:this%n-1) = a06 * ( f(:,:,3:this%n  ) - f(:,:,2:this%n-1) ) &
                            + b06 * ( f(:,:,4:this%n+1) - f(:,:,1:this%n-2) ) 
        
        RHS(:,:,this%n    ) = a06 * ( f(:,:,this%n+1)   - f(:,:,this%n  ) ) &
                            + b06 * ( f(:,:,2)          - f(:,:,this%n-1) )

   end subroutine  



    pure subroutine ComputeZD1RHS_C2E_REAL_periodic(this,f, RHS, n1, n2) 
         
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in)  :: f
        real(rkind), dimension(n1,n2,this%nE) , intent(out) :: RHS
        real(rkind), parameter :: a = (63._rkind/62._rkind) , b = (17._rkind/62._rkind)/three
        real(rkind) :: a06, b06

        a06 = a * this%onebydx
        b06 = b * this%onebydx
        
        RHS(:,:,1         ) = a06 * ( f(:,:,1)          - f(:,:,this%n       ) ) &
                            + b06 * ( f(:,:,2)          - f(:,:,this%n-1 ) ) 
        
        RHS(:,:,2         ) = a06 * ( f(:,:,2 ) - f(:,:,1) ) &
                            + b06 * ( f(:,:,3) - f(:,:,this%n) ) 
        
        RHS(:,:,3:this%n-1) = a06 * ( f(:,:,3:this%n-1) - f(:,:,2:this%n-2) ) &
                            + b06 * ( f(:,:,4:this%n) - f(:,:,1:this%n-3) ) 
        
        RHS(:,:,this%n) = a06 * ( f(:,:,this%n) - f(:,:,this%n-1) ) &
                            + b06 * ( f(:,:,1) - f(:,:,this%n-2) ) 
        
        RHS(:,:,this%n+1) = a06 * ( f(:,:,1) - f(:,:,this%n) ) &
                            + b06 * ( f(:,:,2) - f(:,:,this%n-1) ) 
   end subroutine  

    pure subroutine ComputeZD1RHS_C2E_CMPLX_periodic(this,f, RHS, n1, n2) 
         
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in)  :: f
        complex(rkind), dimension(n1,n2,this%nE) , intent(out) :: RHS
        real(rkind), parameter :: a = (63._rkind/62._rkind) , b = (17._rkind/62._rkind)/three
        complex(rkind) :: a06, b06

        a06 = a * this%onebydx
        b06 = b * this%onebydx
        
        RHS(:,:,1         ) = a06 * ( f(:,:,1)          - f(:,:,this%n       ) ) &
                            + b06 * ( f(:,:,2)          - f(:,:,this%n-1 ) ) 
        
        RHS(:,:,2         ) = a06 * ( f(:,:,2 ) - f(:,:,1) ) &
                            + b06 * ( f(:,:,3) - f(:,:,this%n) ) 
        
        RHS(:,:,3:this%n-1) = a06 * ( f(:,:,3:this%n-1) - f(:,:,2:this%n-2) ) &
                            + b06 * ( f(:,:,4:this%n) - f(:,:,1:this%n-3) ) 
        
        RHS(:,:,this%n) = a06 * ( f(:,:,this%n) - f(:,:,this%n-1) ) &
                            + b06 * ( f(:,:,1) - f(:,:,this%n-2) ) 
        
        RHS(:,:,this%n+1) = a06 * ( f(:,:,1) - f(:,:,this%n) ) &
                            + b06 * ( f(:,:,2) - f(:,:,this%n-1) ) 
   end subroutine  

    pure subroutine ComputeZInterpRHS_E2C_REAL_periodic(this,f, RHS, n1, n2) 
         
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%nE), intent(in)  :: f
        real(rkind), dimension(n1,n2,this%n) , intent(out) :: RHS
        real(rkind), parameter :: a = (3._rkind/2._rkind) , b = (1._rkind/10._rkind)
        real(rkind) :: a06, b06
        
        a06 = a * (1._rkind/2._rkind)
        b06 = b * (1._rkind/2._rkind)

        RHS(:,:,1         ) = a06 * ( f(:,:,2)          + f(:,:,1       ) ) &
                            + b06 * ( f(:,:,3)          + f(:,:,this%n  ) ) 
        
        RHS(:,:,2:this%n-1) = a06 * ( f(:,:,3:this%n  ) + f(:,:,2:this%n-1) ) &
                            + b06 * ( f(:,:,4:this%n+1) + f(:,:,1:this%n-2) ) 
        
        RHS(:,:,this%n    ) = a06 * ( f(:,:,this%n+1)   + f(:,:,this%n  ) ) &
                            + b06 * ( f(:,:,2)          + f(:,:,this%n-1) )

   end subroutine  

    pure subroutine ComputeZInterpRHS_E2C_CMPLX_periodic(this,f, RHS, n1, n2) 
         
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%nE), intent(in)  :: f
        complex(rkind), dimension(n1,n2,this%n) , intent(out) :: RHS
        real(rkind), parameter :: a = (3._rkind/2._rkind) , b = (1._rkind/10._rkind)
        real(rkind) :: a06, b06
        
        a06 = a * (1._rkind/2._rkind)
        b06 = b * (1._rkind/2._rkind)

        RHS(:,:,1         ) = a06 * ( f(:,:,2)          + f(:,:,1       ) ) &
                            + b06 * ( f(:,:,3)          + f(:,:,this%n  ) ) 
        
        RHS(:,:,2:this%n-1) = a06 * ( f(:,:,3:this%n  ) + f(:,:,2:this%n-1) ) &
                            + b06 * ( f(:,:,4:this%n+1) + f(:,:,1:this%n-2) ) 
        
        RHS(:,:,this%n    ) = a06 * ( f(:,:,this%n+1)   + f(:,:,this%n  ) ) &
                            + b06 * ( f(:,:,2)          + f(:,:,this%n-1) )

   end subroutine  


    pure subroutine ComputeZInterpRHS_C2E_REAL_periodic(this,f,RHS,n1,n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in)  :: f
        real(rkind), dimension(n1,n2,this%nE) , intent(out) :: RHS
        real(rkind), parameter :: a = (3._rkind/2._rkind) , b = (1._rkind/10._rkind)
        real(rkind) :: a06, b06
        
        a06 = a * (1._rkind/2._rkind)
        b06 = b * (1._rkind/2._rkind)


        RHS(:,:,1         ) = a06 * ( f(:,:,1)          + f(:,:,this%n       ) ) &
                            + b06 * ( f(:,:,2)          + f(:,:,this%n-1 ) ) 
        
        RHS(:,:,2         ) = a06 * ( f(:,:,2 ) + f(:,:,1) ) &
                            + b06 * ( f(:,:,3) + f(:,:,this%n) ) 
        
        RHS(:,:,3:this%n-1) = a06 * ( f(:,:,3:this%n-1) + f(:,:,2:this%n-2) ) &
                            + b06 * ( f(:,:,4:this%n) + f(:,:,1:this%n-3) ) 
        
        RHS(:,:,this%n) = a06 * ( f(:,:,this%n) + f(:,:,this%n-1) ) &
                            + b06 * ( f(:,:,1) + f(:,:,this%n-2) ) 
        
        RHS(:,:,this%n+1) = a06 * ( f(:,:,1) + f(:,:,this%n) ) &
                            + b06 * ( f(:,:,2) + f(:,:,this%n-1) )
 
    end subroutine
   
    pure subroutine ComputeZInterpRHS_C2E_CMPLX_periodic(this,f,RHS,n1,n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in)  :: f
        complex(rkind), dimension(n1,n2,this%nE) , intent(out) :: RHS
        real(rkind), parameter :: a = (3._rkind/2._rkind) , b = (1._rkind/10._rkind)
        real(rkind) :: a06, b06
        
        a06 = a * (1._rkind/2._rkind)
        b06 = b * (1._rkind/2._rkind)


        RHS(:,:,1         ) = a06 * ( f(:,:,1)          + f(:,:,this%n       ) ) &
                            + b06 * ( f(:,:,2)          + f(:,:,this%n-1 ) ) 
        
        RHS(:,:,2         ) = a06 * ( f(:,:,2 ) + f(:,:,1) ) &
                            + b06 * ( f(:,:,3) + f(:,:,this%n) ) 
        
        RHS(:,:,3:this%n-1) = a06 * ( f(:,:,3:this%n-1) + f(:,:,2:this%n-2) ) &
                            + b06 * ( f(:,:,4:this%n) + f(:,:,1:this%n-3) ) 
        
        RHS(:,:,this%n) = a06 * ( f(:,:,this%n) + f(:,:,this%n-1) ) &
                            + b06 * ( f(:,:,1) + f(:,:,this%n-2) ) 
        
        RHS(:,:,this%n+1) = a06 * ( f(:,:,1) + f(:,:,this%n) ) &
                            + b06 * ( f(:,:,2) + f(:,:,this%n-1) )
 
    end subroutine
 
   subroutine ComputeLU(LU,n,b,d,a)
    
        integer, intent(in) :: n
        real(rkind), intent(in) :: d,a,b
        real(rkind), dimension(n,5), intent(out) :: LU
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
   
            ! Overwrite aa by aa*c
            aa = aa*c
            
            ! Overwrite v by v*c
            v = v*c 
        end associate
    
    end subroutine
    
    pure subroutine ComputeD1RHS_E2C_REAL(this, fE, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%nE), intent(in) :: fE
        real(rkind), dimension(n1,n2,this%n), intent(out) :: rhs
    
#include "STAGG_CD06_files/D1RHS_E2C_common.F90"
    end subroutine
   
    pure subroutine ComputeD1RHS_E2C_CMPLX(this, fE, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%nE), intent(in) :: fE
        complex(rkind), dimension(n1,n2,this%n), intent(out) :: rhs
    
#include "STAGG_CD06_files/D1RHS_E2C_common.F90"
    end subroutine
   

    pure subroutine ComputeD1RHS_C2E_REAL(this, fC, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in)   :: fC
        real(rkind), dimension(n1,n2,this%nE), intent(out) :: rhs
    
#include "STAGG_CD06_files/D1RHS_C2E_common.F90"
    end subroutine

    pure subroutine ComputeD1RHS_C2E_CMPLX(this, fC, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in)   :: fC
        complex(rkind), dimension(n1,n2,this%nE), intent(out) :: rhs
    
#include "STAGG_CD06_files/D1RHS_C2E_common.F90"
    end subroutine

    pure subroutine ComputeD1RHS_C2C_REAL(this, fC, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in)   :: fC
        real(rkind), dimension(n1,n2,this%n), intent(out) :: rhs
    
#include "STAGG_CD06_files/D1RHS_C2C_common.F90"
    end subroutine

    pure subroutine ComputeD1RHS_C2C_CMPLX(this, fC, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in)   :: fC
        complex(rkind), dimension(n1,n2,this%n), intent(out) :: rhs
    
#include "STAGG_CD06_files/D1RHS_C2C_common.F90"
    end subroutine

    pure subroutine ComputeD1RHS_E2E_REAL(this, fE, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%nE), intent(in)   :: fE
        real(rkind), dimension(n1,n2,this%nE), intent(out) :: rhs
    
#include "STAGG_CD06_files/D1RHS_E2E_common.F90"
    end subroutine

    pure subroutine ComputeD1RHS_E2E_CMPLX(this, fE, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%nE), intent(in)   :: fE
        complex(rkind), dimension(n1,n2,this%nE), intent(out) :: rhs
    
#include "STAGG_CD06_files/D1RHS_E2E_common.F90"
    end subroutine

    pure subroutine ComputeInterpRHS_E2C_REAL(this, fE, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%nE), intent(in) :: fE
        real(rkind), dimension(n1,n2,this%n), intent(out) :: rhs
    
#include "STAGG_CD06_files/InterpRHS_E2C_common.F90"
    end subroutine
   
    pure subroutine ComputeInterpRHS_E2C_CMPLX(this, fE, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%nE), intent(in) :: fE
        complex(rkind), dimension(n1,n2,this%n), intent(out) :: rhs
    
#include "STAGG_CD06_files/InterpRHS_E2C_common.F90"
    end subroutine


    pure subroutine ComputeInterpRHS_C2E_REAL(this, fC, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in) :: fC
        real(rkind), dimension(n1,n2,this%nE), intent(out) :: rhs
    
#include "STAGG_CD06_files/InterpRHS_C2E_common.F90"
    end subroutine
   
    pure subroutine ComputeInterpRHS_C2E_CMPLX(this, fC, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in) :: fC
        complex(rkind), dimension(n1,n2,this%nE), intent(out) :: rhs
    
#include "STAGG_CD06_files/InterpRHS_C2E_common.F90"
    end subroutine

    pure subroutine ComputeD2RHS_C2C_REAL(this, fC, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in)   :: fC
        real(rkind), dimension(n1,n2,this%n), intent(out) :: rhs
    
#include "STAGG_CD06_files/D2RHS_C2C_common.F90"
    end subroutine

    pure subroutine ComputeD2RHS_C2C_CMPLX(this, fC, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in)   :: fC
        complex(rkind), dimension(n1,n2,this%n), intent(out) :: rhs
    
#include "STAGG_CD06_files/D2RHS_C2C_common.F90"
    end subroutine

    pure subroutine ComputeD2RHS_E2E_REAL(this, fE, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%nE), intent(in)   :: fE
        real(rkind), dimension(n1,n2,this%nE), intent(out) :: rhs
    
#include "STAGG_CD06_files/D2RHS_E2E_common.F90"
    end subroutine

    pure subroutine ComputeD2RHS_E2E_CMPLX(this, fE, rhs, n1, n2) 
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%nE), intent(in)   :: fE
        complex(rkind), dimension(n1,n2,this%nE), intent(out) :: rhs
    
#include "STAGG_CD06_files/D2RHS_E2E_common.F90"
    end subroutine

    subroutine ddz_E2C_REAL(this, fE, dfC, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%nE), intent(in)  :: fE
        real(rkind), dimension(n1,n2,this%n) , intent(out) :: dfC

        if (this%AmIPeriodic) then
            call this%ComputeZD1RHS_E2C_REAL_periodic(fE(:,:,1:this%n), dfC, n1, n2)
            call this%SolveZLU_REAL(dfC,n1,n2,this%LU_D1)
        else
            call this%ComputeD1RHS_E2C_REAL(fE, dfC, n1, n2)
            call SolveZTriREAL(this%n,this%TriD1_E2C, dfC,n1,n2)
        end if
    end subroutine
   
    subroutine ddz_E2C_CMPLX(this, fE, dfC, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%nE), intent(in)  :: fE
        complex(rkind), dimension(n1,n2,this%n) , intent(out) :: dfC

        if (this%AmIPeriodic) then
            call this%ComputeZD1RHS_E2C_CMPLX_periodic(fE(:,:,1:this%n), dfC, n1, n2)
            call this%SolveZLU_CMPLX(dfC,n1,n2,this%LU_D1)
        else
            call this%ComputeD1RHS_E2C_CMPLX(fE, dfC, n1, n2)
            call SolveZTriCMPLX(this%n,this%TriD1_E2C, dfC,n1,n2)
        end if 
    end subroutine

    subroutine ddz_C2E_REAL(this, fC, dfE, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in)  :: fC
        real(rkind), dimension(n1,n2,this%nE) , intent(out) :: dfE

        if (this%AmIPeriodic) then
            call this%ComputeZD1RHS_C2E_REAL_periodic(fC(:,:,1:this%n), dfE, n1, n2)
            call this%SolveZLU_REAL(dfE,n1,n2,this%LU_D1)
            dfE(:,:,this%nE) = dfE(:,:,1) !needed so that SolveZLU can work for C2E and E2C
        else
            call this%ComputeD1RHS_C2E_REAL(fC, dfE, n1, n2)
            call SolveZTriREAL(this%nE,this%TriD1_C2E, dfE,n1,n2)
        end if 
    end subroutine
    
    subroutine ddz_C2E_CMPLX(this, fC, dfE, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in)  :: fC
        complex(rkind), dimension(n1,n2,this%nE) , intent(out) :: dfE

        if (this%AmIPeriodic) then
            call this%ComputeZD1RHS_C2E_CMPLX_periodic(fC(:,:,1:this%n), dfE, n1, n2)
            call this%SolveZLU_CMPLX(dfE,n1,n2,this%LU_D1)
            dfE(:,:,this%nE) = dfE(:,:,1) !needed so that SolveZLU can work for C2E and E2C
        else
            call this%ComputeD1RHS_C2E_CMPLX(fC, dfE, n1, n2)
            call SolveZTriCMPLX(this%nE,this%TriD1_C2E, dfE,n1,n2)
        end if 

    end subroutine

    subroutine ddz_C2C_REAL(this, fC, dfC, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in)  :: fC
        real(rkind), dimension(n1,n2,this%n) , intent(out) :: dfC

        call this%ComputeD1RHS_C2C_REAL(fC, dfC, n1, n2)
        call SolveZTriREAL(this%n,this%TriD1_C2C, dfC,n1,n2)

    end subroutine
    
    subroutine ddz_C2C_CMPLX(this, fC, dfC, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in)  :: fC
        complex(rkind), dimension(n1,n2,this%n) , intent(out) :: dfC

        call this%ComputeD1RHS_C2C_CMPLX(fC, dfC, n1, n2)
        call SolveZTriCMPLX(this%n,this%TriD1_C2C, dfC,n1,n2)

    end subroutine

    subroutine ddz_E2E_REAL(this, fE, dfE, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%nE), intent(in)  :: fE
        real(rkind), dimension(n1,n2,this%nE) , intent(out) :: dfE

        call this%ComputeD1RHS_E2E_REAL(fE, dfE, n1, n2)
        call SolveZTriREAL(this%nE,this%TriD1_E2E, dfE,n1,n2)

    end subroutine
    
    subroutine ddz_E2E_CMPLX(this, fE, dfE, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%nE), intent(in)  :: fE
        complex(rkind), dimension(n1,n2,this%nE) , intent(out) :: dfE

        call this%ComputeD1RHS_E2E_CMPLX(fE, dfE, n1, n2)
        call SolveZTriCMPLX(this%nE,this%TriD1_E2E, dfE,n1,n2)

    end subroutine


    subroutine InterpZ_E2C_REAL(this, fE, IfC, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%nE), intent(in)  :: fE
        real(rkind), dimension(n1,n2,this%n) , intent(out) :: IfC

        if (this%AmIPeriodic) then
            call this%ComputeZInterpRHS_E2C_REAL_periodic(fE(:,:,1:this%n), IfC, n1, n2)
            call this%SolveZLU_REAL(IfC,n1,n2,this%LU_Interp)
        else
            call this%ComputeInterpRHS_E2C_REAL(fE, IfC, n1, n2)
            call SolveZTriREAL(this%n,this%TriInterp_E2C, IfC,n1,n2)
        end if

    end subroutine
    
    subroutine InterpZ_E2C_CMPLX(this, fE, IfC, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%nE), intent(in)  :: fE
        complex(rkind), dimension(n1,n2,this%n) , intent(out) :: IfC

        if (this%AmIPeriodic) then
            call this%ComputeZInterpRHS_E2C_CMPLX_periodic(fE(:,:,1:this%n), IfC, n1, n2)
            call this%SolveZLU_CMPLX(IfC,n1,n2,this%LU_Interp)
        else
            call this%ComputeInterpRHS_E2C_CMPLX(fE, IfC, n1, n2)
            call SolveZTriCMPLX(this%n,this%TriInterp_E2C, IfC,n1,n2)
        end if

    end subroutine

    subroutine InterpZ_C2E_REAL(this, fC, IfE, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in)  :: fC
        real(rkind), dimension(n1,n2,this%nE) , intent(out) :: IfE

        if (this%AmIPeriodic) then
            call this%ComputeZInterpRHS_C2E_REAL_periodic(fC, IfE, n1, n2)
            call this%SolveZLU_REAL(IfE,n1,n2,this%LU_Interp)
            IfE(:,:,this%nE) = IfE(:,:,1) !needed so that SolveZLU can work for C2E and E2C
        else
            call this%ComputeInterpRHS_C2E_REAL(fC, IfE, n1, n2)
            call SolveZTriREAL(this%nE,this%TriInterp_C2E, IfE,n1,n2)
        end if

    end subroutine
    
    subroutine InterpZ_C2E_CMPLX(this, fC, IfE, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in)  :: fC
        complex(rkind), dimension(n1,n2,this%nE) , intent(out) :: IfE

        
        if (this%AmIPeriodic) then
            call this%ComputeZInterpRHS_C2E_CMPLX_periodic(fC(:,:,1:this%n), IfE, n1, n2)
            call this%SolveZLU_CMPLX(IfE,n1,n2,this%LU_Interp)
            IfE(:,:,this%nE) = IfE(:,:,1) !needed so that SolveZLU can work for C2E and E2C
        else
            call this%ComputeInterpRHS_C2E_CMPLX(fC, IfE, n1, n2)
            call SolveZTriCMPLX(this%nE,this%TriInterp_C2E, IfE,n1,n2)
        end if

    end subroutine

    subroutine d2dz2_C2C_REAL(this, fC, dfC, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in)  :: fC
        real(rkind), dimension(n1,n2,this%n) , intent(out) :: dfC

        call this%ComputeD2RHS_C2C_REAL(fC, dfC, n1, n2)
        call SolveZTriREAL(this%n,this%TriD2_C2C, dfC,n1,n2)
    
    end subroutine

    subroutine d2dz2_C2C_CMPLX(this, fC, dfC, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in)  :: fC
        complex(rkind), dimension(n1,n2,this%n) , intent(out) :: dfC

        call this%ComputeD2RHS_C2C_CMPLX(fC, dfC, n1, n2)
        call SolveZTriCMPLX(this%n,this%TriD2_C2C, dfC,n1,n2)
    
    end subroutine

    subroutine d2dz2_E2E_REAL(this, fE, dfE, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%nE), intent(in)  :: fE
        real(rkind), dimension(n1,n2,this%nE) , intent(out) :: dfE

        call this%ComputeD2RHS_E2E_REAL(fE, dfE, n1, n2)
        call SolveZTriREAL(this%nE,this%TriD2_E2E, dfE,n1,n2)
    
    end subroutine

    subroutine d2dz2_E2E_CMPLX(this, fE, dfE, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%nE), intent(in)  :: fE
        complex(rkind), dimension(n1,n2,this%nE) , intent(out) :: dfE

        call this%ComputeD2RHS_E2E_CMPLX(fE, dfE, n1, n2)
        call SolveZTriCMPLX(this%nE,this%TriD2_E2E, dfE,n1,n2)
    
    end subroutine


end module
