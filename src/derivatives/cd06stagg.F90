! Routines specific to 6th order (stagerred) Compact Finite Differencing scheme

module cd06staggstuff

    use kind_parameters, only: rkind
    use constants,       only : zero,one,two, three, four, ten
    use exits,           only : GracefulExit
    implicit none

    private
    public :: cd06stagg
    
    type cd06stagg

        private

        integer     :: n
        integer     :: nE
        real(rkind) :: dx
        real(rkind) :: onebydx
        real(rkind) :: onebydx2

        logical     :: isTopEven                          ! Boundary condition type. 
        logical     :: isBotEven                          ! Boundary condition type. 

        real(rkind), allocatable, dimension(:,:) :: TriD1_E2C
        real(rkind), allocatable, dimension(:,:) :: TriD1_C2E
        real(rkind), allocatable, dimension(:,:) :: TriD1_C2C
        real(rkind), allocatable, dimension(:,:) :: TriD1_E2E
        real(rkind), allocatable, dimension(:,:) :: TriInterp_E2C
        real(rkind), allocatable, dimension(:,:) :: TriInterp_C2E

        contains

        procedure :: init
        procedure :: destroy
        procedure, private :: ComputeD1RHS_E2C_REAL 
        procedure, private :: ComputeD1RHS_C2E_REAL
        procedure, private :: ComputeD1RHS_E2E_REAL
        procedure, private :: ComputeD1RHS_C2C_REAL
        procedure, private :: ComputeD1RHS_E2C_CMPLX 
        procedure, private :: ComputeD1RHS_C2E_CMPLX
        procedure, private :: ComputeD1RHS_E2E_CMPLX
        procedure, private :: ComputeD1RHS_C2C_CMPLX

        procedure, private :: ComputeInterpRHS_C2E_REAL
        procedure, private :: ComputeInterpRHS_C2E_CMPLX
        procedure, private :: ComputeInterpRHS_E2C_REAL
        procedure, private :: ComputeInterpRHS_E2C_CMPLX
        
        procedure, private :: ComputeTriD1_E2C
        procedure, private :: ComputeTriD1_C2E
        procedure, private :: ComputeTriD1_E2E
        procedure, private :: ComputeTriD1_C2C

        procedure, private :: ComputeTriInterp_E2C
        procedure, private :: ComputeTriInterp_C2E
        
        procedure, private :: ddz_E2C_REAL 
        procedure, private :: ddz_C2E_REAL
        procedure, private :: ddz_C2C_REAL
        procedure, private :: ddz_E2E_REAL
        procedure, private :: ddz_E2C_CMPLX 
        procedure, private :: ddz_C2E_CMPLX
        procedure, private :: ddz_C2C_CMPLX
        procedure, private :: ddz_E2E_CMPLX
        
        generic :: ddz_E2C => ddz_E2C_REAL, ddz_E2C_CMPLX
        generic :: ddz_C2E => ddz_C2E_REAL, ddz_C2E_CMPLX
        generic :: ddz_E2E => ddz_E2E_REAL, ddz_E2E_CMPLX
        generic :: ddz_C2C => ddz_C2C_REAL, ddz_C2C_CMPLX

        procedure, private :: InterpZ_E2C_REAL 
        procedure, private :: InterpZ_C2E_REAL
        procedure, private :: InterpZ_E2C_CMPLX 
        procedure, private :: InterpZ_C2E_CMPLX

        generic :: InterpZ_E2C => InterpZ_E2C_REAL, InterpZ_E2C_CMPLX
        generic :: InterpZ_C2E => InterpZ_C2E_REAL, InterpZ_C2E_CMPLX

    end type

contains
    subroutine init(this, nx, dx, isTopEven, isBotEven) 
        class( cd06stagg ), intent(inout) :: this
        integer, intent(in) :: nx
        real(rkind), intent(in) :: dx
        logical, intent(in) :: isTopEven, isBotEven
    
        this%n = nx; this%nE = nx + 1
        this%dx = dx; this%onebydx = one/dx; this%onebydx2 = this%onebydx/dx
    
        this%isTopEven = isTopEven; this%isBotEven = isBotEven
  
        if (nx .LE. 4) then
            call GracefulExit("CD06_stagg requires at least 4 points",21)
        end if  
        allocate( this%TriD1_E2C(nx  ,3) ); allocate( this%TriD1_C2E(nx+1,3) )
        allocate( this%TriD1_E2E(nx+1,3) ); allocate( this%TriD1_C2C(nx  ,3) )
  
        allocate( this%TriInterp_E2C(nx  ,3) ); allocate( this%TriInterp_C2E(nx+1,3) )
        
        call this%ComputeTriD1_E2C(); call this%ComputeTriD1_C2E()
        call this%ComputeTriD1_E2E(); call this%ComputeTriD1_C2C()
    
        call this%ComputeTriInterp_E2C(); call this%ComputeTriInterp_C2E()
    end subroutine

    subroutine destroy(this)
        class( cd06stagg ), intent(inout) :: this
        deallocate( this%TriD1_E2C, this%TriD1_C2E, this%TriD1_E2E, this%TriD1_C2C )
        deallocate( this%TriInterp_E2C, this%TriInterp_C2E)
    end subroutine
    
#include "STAGG_CD06_files/ComputeTri_allRoutines.F90"
#include "STAGG_CD06_files/TridiagSolver_allRoutines.F90"

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

   
    subroutine ddz_E2C_REAL(this, fE, dfC, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%nE), intent(in)  :: fE
        real(rkind), dimension(n1,n2,this%n) , intent(out) :: dfC

        call this%ComputeD1RHS_E2C_REAL(fE, dfC, n1, n2)
        call SolveZTriREAL(this%n,this%TriD1_E2C, dfC,n1,n2)

    end subroutine
    
    subroutine ddz_E2C_CMPLX(this, fE, dfC, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%nE), intent(in)  :: fE
        complex(rkind), dimension(n1,n2,this%n) , intent(out) :: dfC

        call this%ComputeD1RHS_E2C_CMPLX(fE, dfC, n1, n2)
        call SolveZTriCMPLX(this%n,this%TriD1_E2C, dfC,n1,n2)

    end subroutine

    subroutine ddz_C2E_REAL(this, fC, dfE, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in)  :: fC
        real(rkind), dimension(n1,n2,this%nE) , intent(out) :: dfE

        call this%ComputeD1RHS_C2E_REAL(fC, dfE, n1, n2)
        call SolveZTriREAL(this%nE,this%TriD1_C2E, dfE,n1,n2)

    end subroutine
    
    subroutine ddz_C2E_CMPLX(this, fC, dfE, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in)  :: fC
        complex(rkind), dimension(n1,n2,this%nE) , intent(out) :: dfE

        call this%ComputeD1RHS_C2E_CMPLX(fC, dfE, n1, n2)
        call SolveZTriCMPLX(this%nE,this%TriD1_C2E, dfE,n1,n2)

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

        call this%ComputeInterpRHS_E2C_REAL(fE, IfC, n1, n2)
        call SolveZTriREAL(this%n,this%TriInterp_E2C, IfC,n1,n2)

    end subroutine
    
    subroutine InterpZ_E2C_CMPLX(this, fE, IfC, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%nE), intent(in)  :: fE
        complex(rkind), dimension(n1,n2,this%n) , intent(out) :: IfC

        call this%ComputeInterpRHS_E2C_CMPLX(fE, IfC, n1, n2)
        call SolveZTriCMPLX(this%n,this%TriInterp_E2C, IfC,n1,n2)

    end subroutine


    subroutine InterpZ_C2E_REAL(this, fC, IfE, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in)  :: fC
        real(rkind), dimension(n1,n2,this%nE) , intent(out) :: IfE

        call this%ComputeInterpRHS_C2E_REAL(fC, IfE, n1, n2)
        call SolveZTriREAL(this%nE,this%TriInterp_C2E, IfE,n1,n2)

    end subroutine
    
    subroutine InterpZ_C2E_CMPLX(this, fC, IfE, n1, n2)
        class( cd06stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in)  :: fC
        complex(rkind), dimension(n1,n2,this%nE) , intent(out) :: IfE

        call this%ComputeInterpRHS_C2E_CMPLX(fC, IfE, n1, n2)
        call SolveZTriCMPLX(this%nE,this%TriInterp_C2E, IfE,n1,n2)

    end subroutine

end module
