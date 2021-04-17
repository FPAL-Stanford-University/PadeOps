 ! Routines specific to 2nd-order 15-diagonal Compact Filter
 ! Periodic LU based on  Neossi Nguetchue, Abelman (Appl. Math. & Comp. 2008)

 module cfo2D15stuff
   
   use kind_parameters, only: rkind
   use constants,       only: zero,one,two,four
   use exits,           only: GracefulExit

   implicit none

   private
   !public :: cfo2d15, alp0_d15, alp1_d15, alp2_d15, alp3_d15, alp4_d15, alp5_d15, alp6_d15, alp7_d15, a0_d15, a1_d15, a2_d15, a3_d15, a4_d15, a5_d15, a6_d15, a7_d15
   public :: cfo2d15, alp0_d15,alp1_d15,alp2_d15,alp3_d15,alp4_d15,alp5_d15,alp6_d15,alp7_d15,a0_d15,a1_d15,a2_d15,a3_d15,a4_d15,a5_d15,a6_d15,a7_d15,b1_alp0_d15,b1_a0_d15,b2_alp0_d15,b2_alp1_d15,b2_alp2_d15,b2_alp3_d15,b2_alp4_d15,b2_alp5_d15,b2_alp6_d15,b2_alp7_d15,b2_a0_d15,b2_a1_d15,b2_a2_d15,b2_a3_d15,b2_a4_d15,b2_a5_d15,b2_a6_d15,b2_a7_d15,b2_bet1_d15,b2_b1_d15,b3_alp0_d15,b3_alp1_d15,b3_alp2_d15,b3_alp3_d15,b3_alp4_d15,b3_alp5_d15,b3_alp6_d15,b3_alp7_d15,b3_a0_d15,b3_a1_d15,b3_a2_d15,b3_a3_d15,b3_a4_d15,b3_a5_d15,b3_a6_d15,b3_a7_d15,b3_bet1_d15,b3_bet2_d15,b3_b1_d15,b3_b2_d15,b4_alp0_d15,b4_alp1_d15,b4_alp2_d15,b4_alp3_d15,b4_alp4_d15,b4_alp5_d15,b4_alp6_d15,b4_alp7_d15,b4_a0_d15,b4_a1_d15,b4_a2_d15,b4_a3_d15,b4_a4_d15,b4_a5_d15,b4_a6_d15,b4_a7_d15,b4_bet1_d15,b4_bet2_d15,b4_bet3_d15,b4_b1_d15,b4_b2_d15,b4_b3_d15,b5_alp0_d15,b5_alp1_d15,b5_alp2_d15,b5_alp3_d15,b5_alp4_d15,b5_alp5_d15,b5_alp6_d15,b5_alp7_d15,b5_a0_d15,b5_a1_d15,b5_a2_d15,b5_a3_d15,b5_a4_d15,b5_a5_d15,b5_a6_d15,b5_a7_d15,b5_bet1_d15,b5_bet2_d15,b5_bet3_d15,b5_bet4_d15,b5_b1_d15,b5_b2_d15,b5_b3_d15,b5_b4_d15,b6_alp0_d15,b6_alp1_d15,b6_alp2_d15,b6_alp3_d15,b6_alp4_d15,b6_alp5_d15,b6_alp6_d15,b6_alp7_d15,b6_a0_d15,b6_a1_d15,b6_a2_d15,b6_a3_d15,b6_a4_d15,b6_a5_d15,b6_a6_d15,b6_a7_d15,b6_bet1_d15,b6_bet2_d15,b6_bet3_d15,b6_bet4_d15,b6_bet5_d15,b6_b1_d15,b6_b2_d15,b6_b3_d15,b6_b4_d15,b6_b5_d15,b7_alp0_d15,b7_alp1_d15,b7_alp2_d15,b7_alp3_d15,b7_alp4_d15,b7_alp5_d15,b7_alp6_d15,b7_alp7_d15,b7_a0_d15,b7_a1_d15,b7_a2_d15,b7_a3_d15,b7_a4_d15,b7_a5_d15,b7_a6_d15,b7_a7_d15,b7_bet1_d15,b7_bet2_d15,b7_bet3_d15,b7_bet4_d15,b7_bet5_d15,b7_bet6_d15,b7_b1_d15,b7_b2_d15,b7_b3_d15,b7_b4_d15,b7_b5_d15,b7_b6_d15  

   ! second-order (framework for additional RHS points for higher (14th) order included) 15-diagonal filter with adjustable coefficients -- (d^2 sf(0)) / (d omega^2) = zero -- set even derivatives of transfer funciton to zero at wavenumber zero to solve for LHS coefficients if second order

   !###############################################################

   !Interior points
   real(rkind) :: alp0_d15  = one
   real(rkind) :: alp1_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: alp2_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: alp3_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: alp4_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: alp5_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: alp6_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: alp7_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: a0_d15    = zero
   real(rkind) :: a1_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: a2_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: a3_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: a4_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: a5_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: a6_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: a7_d15    = zero ! Already divided by factor of 2 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !! NOTE : The following variables are used for non-periodic filter !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

   ! Point 1 is just identity
   real(rkind) :: b1_alp0_d15  = one
   real(rkind) :: b1_a0_d15    = one

   ! Point 2
   !Interior coefficients
   real(rkind) :: b2_alp0_d15  = one
   real(rkind) :: b2_alp1_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b2_alp2_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b2_alp3_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b2_alp4_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b2_alp5_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b2_alp6_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b2_alp7_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b2_a0_d15    = zero
   real(rkind) :: b2_a1_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b2_a2_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b2_a3_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b2_a4_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b2_a5_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b2_a6_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b2_a7_d15    = zero ! Already divided by factor of 2 
   !Exterior coefficients
   real(rkind) :: b2_bet1_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b2_b1_d15    = zero ! Already divided by factor of 2 

   ! Point 3
   !Interior coefficients
   real(rkind) :: b3_alp0_d15  = one
   real(rkind) :: b3_alp1_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b3_alp2_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b3_alp3_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b3_alp4_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b3_alp5_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b3_alp6_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b3_alp7_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b3_a0_d15    = zero
   real(rkind) :: b3_a1_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b3_a2_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b3_a3_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b3_a4_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b3_a5_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b3_a6_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b3_a7_d15    = zero ! Already divided by factor of 2 
   !Exterior coefficients
   real(rkind) :: b3_bet1_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b3_bet2_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b3_b1_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b3_b2_d15    = zero ! Already divided by factor of 2 

   ! Point 4
   !Interior coefficients
   real(rkind) :: b4_alp0_d15  = one
   real(rkind) :: b4_alp1_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b4_alp2_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b4_alp3_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b4_alp4_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b4_alp5_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b4_alp6_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b4_alp7_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b4_a0_d15    = zero
   real(rkind) :: b4_a1_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b4_a2_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b4_a3_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b4_a4_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b4_a5_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b4_a6_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b4_a7_d15    = zero ! Already divided by factor of 2 
   !Exterior coefficients
   real(rkind) :: b4_bet1_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b4_bet2_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b4_bet3_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b4_b1_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b4_b2_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b4_b3_d15    = zero ! Already divided by factor of 2 

   ! Point 5
   !Interior coefficients
   real(rkind) :: b5_alp0_d15  = one
   real(rkind) :: b5_alp1_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b5_alp2_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b5_alp3_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b5_alp4_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b5_alp5_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b5_alp6_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b5_alp7_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b5_a0_d15    = zero
   real(rkind) :: b5_a1_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b5_a2_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b5_a3_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b5_a4_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b5_a5_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b5_a6_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b5_a7_d15    = zero ! Already divided by factor of 2 
   !Exterior coefficients
   real(rkind) :: b5_bet1_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b5_bet2_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b5_bet3_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b5_bet4_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b5_b1_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b5_b2_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b5_b3_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b5_b4_d15    = zero ! Already divided by factor of 2 

   ! Point 6
   !Interior coefficients
   real(rkind) :: b6_alp0_d15  = one
   real(rkind) :: b6_alp1_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b6_alp2_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b6_alp3_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b6_alp4_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b6_alp5_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b6_alp6_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b6_alp7_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b6_a0_d15    = zero
   real(rkind) :: b6_a1_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b6_a2_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b6_a3_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b6_a4_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b6_a5_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b6_a6_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b6_a7_d15    = zero ! Already divided by factor of 2 
   !Exterior coefficients
   real(rkind) :: b6_bet1_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b6_bet2_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b6_bet3_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b6_bet4_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b6_bet5_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b6_b1_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b6_b2_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b6_b3_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b6_b4_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b6_b5_d15    = zero ! Already divided by factor of 2 

   ! Point 7
   !Interior coefficients
   real(rkind) :: b7_alp0_d15  = one
   real(rkind) :: b7_alp1_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b7_alp2_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b7_alp3_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b7_alp4_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b7_alp5_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b7_alp6_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b7_alp7_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b7_a0_d15    = zero
   real(rkind) :: b7_a1_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b7_a2_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b7_a3_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b7_a4_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b7_a5_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b7_a6_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b7_a7_d15    = zero ! Already divided by factor of 2 
   !Exterior coefficients
   real(rkind) :: b7_bet1_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b7_bet2_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b7_bet3_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b7_bet4_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b7_bet5_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b7_bet6_d15  = zero ! Already divided by factor of 2 
   real(rkind) :: b7_b1_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b7_b2_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b7_b3_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b7_b4_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b7_b5_d15    = zero ! Already divided by factor of 2 
   real(rkind) :: b7_b6_d15    = zero ! Already divided by factor of 2 

   !##################################################################


   type cfo2d15

      private

      integer     :: n

      logical     :: periodic=.TRUE.

      integer, allocatable, dimension(:) :: LU_ipiv
      integer, allocatable, dimension(:) :: d15_ipiv_nn, d15_ipiv_ns, d15_ipiv_na, d15_ipiv_sn, d15_ipiv_ss, d15_ipiv_sa, d15_ipiv_an, d15_ipiv_as, d15_ipiv_aa

      real(rkind), allocatable, dimension(:,:) :: LU
      real(rkind), allocatable, dimension(:,:) :: d15_nn
      real(rkind), allocatable, dimension(:,:) :: d15_ns
      real(rkind), allocatable, dimension(:,:) :: d15_na
      real(rkind), allocatable, dimension(:,:) :: d15_sn
      real(rkind), allocatable, dimension(:,:) :: d15_ss
      real(rkind), allocatable, dimension(:,:) :: d15_sa
      real(rkind), allocatable, dimension(:,:) :: d15_an
      real(rkind), allocatable, dimension(:,:) :: d15_as
      real(rkind), allocatable, dimension(:,:) :: d15_aa

    contains

      procedure :: init
      procedure :: destroy

      procedure, private :: ComputeXRHS
      procedure, private :: ComputeYRHS
      procedure, private :: ComputeZRHS_REAL
      procedure, private :: ComputeZRHS_CMPLX

      procedure, private :: SolveXLU
      procedure, private :: SolveYLU
      procedure, private :: SolveZLU

      procedure, private :: ComputeLU
      procedure, private :: Computed15

      procedure, private :: SolveXd15
      procedure, private :: SolveYd15
      procedure, private :: SolveZd15_REAL
      procedure, private :: SolveZd15_CMPLX

      procedure :: filter1
      procedure :: filter2
      procedure, private :: filter3_REAL
      procedure, private :: filter3_CMPLX
      generic :: filter3 => filter3_REAL, filter3_CMPLX 

   end type cfo2d15



 contains

   function init(this, n_, periodic_, alp1_d15_) result(ierr)
     class( cfo2d15 ), intent(inout) :: this
     integer, intent(in) :: n_
     logical, intent(in) :: periodic_
     real(rkind), intent(in) :: alp1_d15_
     integer :: ierr

     this%n = n_

     this%periodic = periodic_

     alp1_d15 = alp1_d15_

     !second order filter -- coeficients set from variable alp1_d15 -- LHS set by solving for higher even derivatives of transfer function equal zero for wavenumber zero
     
     !########################################################################################
     
     !Interior points
     alp0_d15  = one
     alp1_d15  = alp1_d15                           ! Already divided by factor of 2 
     alp2_d15  =  (7. /33.  )*( 1. - 2.*alp1_d15)   ! Already divided by factor of 2 
     alp3_d15  = -(7. /66.  )*( 1. - 2.*alp1_d15)   ! Already divided by factor of 2 
     alp4_d15  =  (14./363. )*( 1. - 2.*alp1_d15)   ! Already divided by factor of 2 
     alp5_d15  = -(7. /726. )*( 1. - 2.*alp1_d15)   ! Already divided by factor of 2 
     alp6_d15  =  (7. /4719.)*( 1. - 2.*alp1_d15)   ! Already divided by factor of 2 
     alp7_d15  = -(1. /9438.)*( 1. - 2.*alp1_d15)   ! Already divided by factor of 2 
     a0_d15    =  (1. /11.  )*( 7. + 8.*alp1_d15)  
     a1_d15    =  (1. /22.  )*( 7. + 8.*alp1_d15)   ! Already divided by factor of 2 
     a2_d15    = zero                               ! Already divided by factor of 2 
     a3_d15    = zero                               ! Already divided by factor of 2 
     a4_d15    = zero                               ! Already divided by factor of 2 
     a5_d15    = zero                               ! Already divided by factor of 2 
     a6_d15    = zero                               ! Already divided by factor of 2 
     a7_d15    = zero                               ! Already divided by factor of 2 

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     !! NOTE : The following variables are used for non-periodic filter !!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

     ! Point 1 is just identity
     b1_alp0_d15  = one
     b1_a0_d15    = one

     ! Point 2
     !Interior coefficients
     b2_alp0_d15  = one
     b2_alp1_d15  = alp1_d15                        ! Already divided by factor of 2 
     b2_alp2_d15  = zero                            ! Already divided by factor of 2 
     b2_alp3_d15  = zero                            ! Already divided by factor of 2 
     b2_alp4_d15  = zero                            ! Already divided by factor of 2 
     b2_alp5_d15  = zero                            ! Already divided by factor of 2 
     b2_alp6_d15  = zero                            ! Already divided by factor of 2 
     b2_alp7_d15  = zero                            ! Already divided by factor of 2 
     b2_a0_d15    =  (1./2. )*( 1. + 2.*alp1_d15)     
     b2_a1_d15    =  (1./4. )*( 1. + 2.*alp1_d15)   ! Already divided by factor of 2 
     b2_a2_d15    = zero                            ! Already divided by factor of 2 
     b2_a3_d15    = zero                            ! Already divided by factor of 2 
     b2_a4_d15    = zero                            ! Already divided by factor of 2 
     b2_a5_d15    = zero                            ! Already divided by factor of 2 
     b2_a6_d15    = zero                            ! Already divided by factor of 2 
     b2_a7_d15    = zero                            ! Already divided by factor of 2 
     !Exterior coefficients                        
     b2_bet1_d15  = alp1_d15                        ! Already divided by factor of 2 
     b2_b1_d15    =  (1./4. )*( 1. + 2.*alp1_d15)   ! Already divided by factor of 2 
                                                   
     ! Point 3                                     
     !Interior coefficients                        
     b3_alp0_d15  = one                            
     b3_alp1_d15  = alp1_d15                        ! Already divided by factor of 2 
     b3_alp2_d15  =  (1./14.)*( 1. - 2.*alp1_d15)   ! Already divided by factor of 2 
     b3_alp3_d15  = zero                            ! Already divided by factor of 2 
     b3_alp4_d15  = zero                            ! Already divided by factor of 2 
     b3_alp5_d15  = zero                            ! Already divided by factor of 2 
     b3_alp6_d15  = zero                            ! Already divided by factor of 2 
     b3_alp7_d15  = zero                            ! Already divided by factor of 2 
     b3_a0_d15    =  (2./7. )*( 2. + 3.*alp1_d15)    
     b3_a1_d15    =  (1./7. )*( 2. + 3.*alp1_d15)   ! Already divided by factor of 2 
     b3_a2_d15    = zero                            ! Already divided by factor of 2 
     b3_a3_d15    = zero                            ! Already divided by factor of 2 
     b3_a4_d15    = zero                            ! Already divided by factor of 2 
     b3_a5_d15    = zero                            ! Already divided by factor of 2 
     b3_a6_d15    = zero                            ! Already divided by factor of 2 
     b3_a7_d15    = zero                            ! Already divided by factor of 2 
     !Exterior coefficients                        
     b3_bet1_d15  = alp1_d15                        ! Already divided by factor of 2 
     b3_bet2_d15  =  (1./14.)*( 1. - 2.*alp1_d15)   ! Already divided by factor of 2 
     b3_b1_d15    =  (1./7. )*( 2. + 3.*alp1_d15)   ! Already divided by factor of 2 
     b3_b2_d15    = zero                            ! Already divided by factor of 2 
                                                   
     ! Point 4                                     
     !Interior coefficients                        
     b4_alp0_d15  = one                            
     b4_alp1_d15  = alp1_d15                        ! Already divided by factor of 2 
     b4_alp2_d15  =  (3./25.)*( 1. - 2.*alp1_d15)   ! Already divided by factor of 2 
     b4_alp3_d15  = -(1./50.)*( 1. - 2.*alp1_d15)   ! Already divided by factor of 2 
     b4_alp4_d15  = zero                            ! Already divided by factor of 2 
     b4_alp5_d15  = zero                            ! Already divided by factor of 2 
     b4_alp6_d15  = zero                            ! Already divided by factor of 2 
     b4_alp7_d15  = zero                            ! Already divided by factor of 2 
     b4_a0_d15    =  (1./5. )*( 3. + 4.*alp1_d15)  
     b4_a1_d15    =  (1./10.)*( 3. + 4.*alp1_d15)   ! Already divided by factor of 2 
     b4_a2_d15    = zero                            ! Already divided by factor of 2 
     b4_a3_d15    = zero                            ! Already divided by factor of 2 
     b4_a4_d15    = zero                            ! Already divided by factor of 2 
     b4_a5_d15    = zero                            ! Already divided by factor of 2 
     b4_a6_d15    = zero                            ! Already divided by factor of 2 
     b4_a7_d15    = zero                            ! Already divided by factor of 2 
     !Exterior coefficients                        
     b4_bet1_d15  = alp1_d15                        ! Already divided by factor of 2 
     b4_bet2_d15  =  (3./25.)*( 1. - 2.*alp1_d15)   ! Already divided by factor of 2 
     b4_bet3_d15  = -(1./50.)*( 1. - 2.*alp1_d15)   ! Already divided by factor of 2 
     b4_b1_d15    =  (1./10.)*( 3. + 4.*alp1_d15)   ! Already divided by factor of 2 
     b4_b2_d15    = zero                            ! Already divided by factor of 2 
     b4_b3_d15    = zero                            ! Already divided by factor of 2 
                                                   
     ! Point 5                                     
     !Interior coefficients                        
     b5_alp0_d15  = one                            
     b5_alp1_d15  = alp1_d15                        ! Already divided by factor of 2 
     b5_alp2_d15  =  (2./13. )*( 1. - 2.*alp1_d15)  ! Already divided by factor of 2 
     b5_alp3_d15  = -(4./91. )*( 1. - 2.*alp1_d15)  ! Already divided by factor of 2 
     b5_alp4_d15  =  (1./182.)*( 1. - 2.*alp1_d15)  ! Already divided by factor of 2 
     b5_alp5_d15  = zero                            ! Already divided by factor of 2 
     b5_alp6_d15  = zero                            ! Already divided by factor of 2 
     b5_alp7_d15  = zero                            ! Already divided by factor of 2 
     b5_a0_d15    =  (2./13. )*( 4. + 5.*alp1_d15) 
     b5_a1_d15    =  (1./13. )*( 4. + 5.*alp1_d15)  ! Already divided by factor of 2 
     b5_a2_d15    = zero                            ! Already divided by factor of 2 
     b5_a3_d15    = zero                            ! Already divided by factor of 2 
     b5_a4_d15    = zero                            ! Already divided by factor of 2 
     b5_a5_d15    = zero                            ! Already divided by factor of 2 
     b5_a6_d15    = zero                            ! Already divided by factor of 2 
     b5_a7_d15    = zero                            ! Already divided by factor of 2 
     !Exterior coefficients                        
     b5_bet1_d15  = alp1_d15                        ! Already divided by factor of 2 
     b5_bet2_d15  =  (2./13. )*( 1. - 2.*alp1_d15)  ! Already divided by factor of 2 
     b5_bet3_d15  = -(4./91. )*( 1. - 2.*alp1_d15)  ! Already divided by factor of 2 
     b5_bet4_d15  =  (1./182.)*( 1. - 2.*alp1_d15)  ! Already divided by factor of 2 
     b5_b1_d15    =  (1./13. )*( 4. + 5.*alp1_d15)  ! Already divided by factor of 2 
     b5_b2_d15    = zero                            ! Already divided by factor of 2 
     b5_b3_d15    = zero                            ! Already divided by factor of 2 
     b5_b4_d15    = zero                            ! Already divided by factor of 2 

     ! Point 6
     !Interior coefficients
     b6_alp0_d15  = one
     b6_alp1_d15  = alp1_d15                        ! Already divided by factor of 2 
     b6_alp2_d15  =  (5. /28. )*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b6_alp3_d15  = -(15./224.)*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b6_alp4_d15  =  (5. /336.)*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b6_alp5_d15  = -(1. /672.)*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b6_alp6_d15  = zero                            ! Already divided by factor of 2 
     b6_alp7_d15  = zero                            ! Already divided by factor of 2 
     b6_a0_d15    =  (1. /8.  )*( 5. + 6.*alp1_d15)
     b6_a1_d15    =  (1. /16. )*( 5. + 6.*alp1_d15) ! Already divided by factor of 2 
     b6_a2_d15    = zero                            ! Already divided by factor of 2 
     b6_a3_d15    = zero                            ! Already divided by factor of 2 
     b6_a4_d15    = zero                            ! Already divided by factor of 2 
     b6_a5_d15    = zero                            ! Already divided by factor of 2 
     b6_a6_d15    = zero                            ! Already divided by factor of 2 
     b6_a7_d15    = zero                            ! Already divided by factor of 2 
     !Exterior coefficients
     b6_bet1_d15  = alp1_d15                        ! Already divided by factor of 2 
     b6_bet2_d15  =  (5. /28. )*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b6_bet3_d15  = -(15./224.)*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b6_bet4_d15  =  (5. /336.)*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b6_bet5_d15  = -(1. /672.)*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b6_b1_d15    =  (1. /16. )*( 5. + 6.*alp1_d15) ! Already divided by factor of 2 
     b6_b2_d15    = zero                            ! Already divided by factor of 2 
     b6_b3_d15    = zero                            ! Already divided by factor of 2 
     b6_b4_d15    = zero                            ! Already divided by factor of 2 
     b6_b5_d15    = zero                            ! Already divided by factor of 2 

     ! Point 7
     !Interior coefficients
     b7_alp0_d15  = one
     b7_alp1_d15  = alp1_d15                        ! Already divided by factor of 2 
     b7_alp2_d15  =  (15./76. )*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b7_alp3_d15  = -(5. /57. )*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b7_alp4_d15  =  (1. /38. )*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b7_alp5_d15  = -(1. /209.)*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b7_alp6_d15  =  (1./2508.)*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b7_alp7_d15  = zero                            ! Already divided by factor of 2 
     b7_a0_d15    =  (2. /19. )*( 6. + 7.*alp1_d15)
     b7_a1_d15    =  (1. /19. )*( 6. + 7.*alp1_d15) ! Already divided by factor of 2 
     b7_a2_d15    = zero                            ! Already divided by factor of 2 
     b7_a3_d15    = zero                            ! Already divided by factor of 2 
     b7_a4_d15    = zero                            ! Already divided by factor of 2 
     b7_a5_d15    = zero                            ! Already divided by factor of 2 
     b7_a6_d15    = zero                            ! Already divided by factor of 2 
     b7_a7_d15    = zero                            ! Already divided by factor of 2 
     !Exterior coefficients
     b7_bet1_d15  = alp1_d15                        ! Already divided by factor of 2 
     b7_bet2_d15  =  (15./76. )*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b7_bet3_d15  = -(5. /57. )*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b7_bet4_d15  =  (1. /38. )*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b7_bet5_d15  = -(1. /209.)*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b7_bet6_d15  =  (1./2508.)*( 1. - 2.*alp1_d15) ! Already divided by factor of 2 
     b7_b1_d15    =  (1. /19. )*( 6. + 7.*alp1_d15) ! Already divided by factor of 2 
     b7_b2_d15    = zero                            ! Already divided by factor of 2 
     b7_b3_d15    = zero                            ! Already divided by factor of 2 
     b7_b4_d15    = zero                            ! Already divided by factor of 2 
     b7_b5_d15    = zero                            ! Already divided by factor of 2 
     b7_b6_d15    = zero                            ! Already divided by factor of 2 
 
     !########################################################################################

     if (periodic_) then 
        ! Allocate LU matrix
        if(allocated( this%LU )) deallocate( this%LU ); allocate( this%LU(n_,n_) )
        ! Allocate LU ipiv
        if(allocated( this%LU_ipiv )) deallocate( this%LU_ipiv ); allocate( this%LU_ipiv(n_) )

        ! Compute LU matrix
        if (n_ .GE. 15) then
           call this%ComputeLU(this%LU,this%LU_ipiv)

        else if (n_ == 1) then
           this%LU = one
        else
           ierr = 7
           return
        end if

     else 

        ! Allocate d15 matrix
        if(allocated( this%d15_nn )) deallocate( this%d15_nn ); allocate( this%d15_nn(n_,n_) )
        if(allocated( this%d15_ns )) deallocate( this%d15_ns ); allocate( this%d15_ns(n_,n_) )
        if(allocated( this%d15_na )) deallocate( this%d15_na ); allocate( this%d15_na(n_,n_) )
        if(allocated( this%d15_sn )) deallocate( this%d15_sn ); allocate( this%d15_sn(n_,n_) )
        if(allocated( this%d15_ss )) deallocate( this%d15_ss ); allocate( this%d15_ss(n_,n_) )
        if(allocated( this%d15_sa )) deallocate( this%d15_sa ); allocate( this%d15_sa(n_,n_) )
        if(allocated( this%d15_an )) deallocate( this%d15_an ); allocate( this%d15_an(n_,n_) )
        if(allocated( this%d15_as )) deallocate( this%d15_as ); allocate( this%d15_as(n_,n_) )
        if(allocated( this%d15_aa )) deallocate( this%d15_aa ); allocate( this%d15_aa(n_,n_) )
        ! Allocate d15 ipiv
        if(allocated( this%d15_ipiv_nn )) deallocate( this%d15_ipiv_nn ); allocate( this%d15_ipiv_nn(n_) )
        if(allocated( this%d15_ipiv_ns )) deallocate( this%d15_ipiv_ns ); allocate( this%d15_ipiv_ns(n_) )
        if(allocated( this%d15_ipiv_na )) deallocate( this%d15_ipiv_na ); allocate( this%d15_ipiv_na(n_) )
        if(allocated( this%d15_ipiv_sn )) deallocate( this%d15_ipiv_sn ); allocate( this%d15_ipiv_sn(n_) )
        if(allocated( this%d15_ipiv_ss )) deallocate( this%d15_ipiv_ss ); allocate( this%d15_ipiv_ss(n_) )
        if(allocated( this%d15_ipiv_sa )) deallocate( this%d15_ipiv_sa ); allocate( this%d15_ipiv_sa(n_) )
        if(allocated( this%d15_ipiv_an )) deallocate( this%d15_ipiv_an ); allocate( this%d15_ipiv_an(n_) )
        if(allocated( this%d15_ipiv_as )) deallocate( this%d15_ipiv_as ); allocate( this%d15_ipiv_as(n_) )
        if(allocated( this%d15_ipiv_aa )) deallocate( this%d15_ipiv_aa ); allocate( this%d15_ipiv_aa(n_) )

        if (n_ .GE. 15) then
           call this%ComputeD15(this%d15_nn, this%d15_ipiv_nn, 0, 0)
           call this%ComputeD15(this%d15_ns, this%d15_ipiv_ns, 0, 1)
           call this%ComputeD15(this%d15_na, this%d15_ipiv_na, 0,-1)
           call this%ComputeD15(this%d15_sn, this%d15_ipiv_sn, 1, 0)
           call this%ComputeD15(this%d15_ss, this%d15_ipiv_ss, 1, 1)
           call this%ComputeD15(this%d15_sa, this%d15_ipiv_sa, 1,-1)
           call this%ComputeD15(this%d15_an, this%d15_ipiv_an,-1, 0)
           call this%ComputeD15(this%d15_as, this%d15_ipiv_as,-1, 1)
           call this%ComputeD15(this%d15_aa, this%d15_ipiv_aa,-1,-1)
        else if (n_ .EQ. 1) then
           this%d15_nn = one
           this%d15_ns = one
           this%d15_na = one
           this%d15_sn = one
           this%d15_ss = one
           this%d15_sa = one
           this%d15_an = one
           this%d15_as = one
           this%d15_aa = one
        else
           ierr = 7
           return 
        end if

     end if

     ! If everything passes
     ierr = 0

   end function init

   subroutine destroy(this)

     class( cfo2d15 ), intent(inout) :: this

     ! Dellocate LU matrix.
     if(allocated( this%LU )) deallocate( this%LU )
     ! Dellocate LUipiv matrix.
     if(allocated( this%LU_ipiv )) deallocate( this%LU_ipiv )

     ! Dellocate d15 matrix.
     if(allocated( this%d15_nn )) deallocate( this%d15_nn )
     if(allocated( this%d15_ns )) deallocate( this%d15_ns )
     if(allocated( this%d15_na )) deallocate( this%d15_na )
     if(allocated( this%d15_sn )) deallocate( this%d15_sn )
     if(allocated( this%d15_ss )) deallocate( this%d15_ss )
     if(allocated( this%d15_sa )) deallocate( this%d15_sa )
     if(allocated( this%d15_an )) deallocate( this%d15_an )
     if(allocated( this%d15_as )) deallocate( this%d15_as )
     if(allocated( this%d15_aa )) deallocate( this%d15_aa )
     ! Dellocate d15 ipiv matrix.
     if(allocated( this%d15_ipiv_nn )) deallocate( this%d15_ipiv_nn )
     if(allocated( this%d15_ipiv_ns )) deallocate( this%d15_ipiv_ns )
     if(allocated( this%d15_ipiv_na )) deallocate( this%d15_ipiv_na )
     if(allocated( this%d15_ipiv_sn )) deallocate( this%d15_ipiv_sn )
     if(allocated( this%d15_ipiv_ss )) deallocate( this%d15_ipiv_ss )
     if(allocated( this%d15_ipiv_sa )) deallocate( this%d15_ipiv_sa )
     if(allocated( this%d15_ipiv_an )) deallocate( this%d15_ipiv_an )
     if(allocated( this%d15_ipiv_as )) deallocate( this%d15_ipiv_as )
     if(allocated( this%d15_ipiv_aa )) deallocate( this%d15_ipiv_aa )

   end subroutine destroy


   subroutine ComputeLU(this,LU,ipiv)
     class( cfo2D15 ), intent(in) :: this
     real(rkind), dimension(this%n,this%n), intent(out) :: LU
     integer, dimension(this%n), intent(out) :: ipiv
     integer :: i,info
     
     LU = zero

     do i=1,this%n
        ! i-7
        if (i-7 .ge. one   ) then
           LU(i,i-7) = alp7_d15
        else
           LU(i,i-7+this%n) = alp7_d15
        endif
        ! i-6
        if (i-6 .ge. one   ) then
           LU(i,i-6) = alp6_d15
        else
           LU(i,i-6+this%n) = alp6_d15
        endif
        ! i-5
        if (i-5 .ge. one   ) then
           LU(i,i-5) = alp5_d15
        else
           LU(i,i-5+this%n) = alp5_d15
        endif
        ! i-4
        if (i-4 .ge. one   ) then
           LU(i,i-4) = alp4_d15
        else
           LU(i,i-4+this%n) = alp4_d15
        endif
        ! i-3
        if (i-3 .ge. one   ) then
           LU(i,i-3) = alp3_d15
        else
           LU(i,i-3+this%n) = alp3_d15
        endif
        ! i-2
        if (i-2 .ge. one   ) then
           LU(i,i-2) = alp2_d15
        else
           LU(i,i-2+this%n) = alp2_d15
        endif
        ! i-1
        if (i-1 .ge. one   ) then
           LU(i,i-1) = alp1_d15
        else
           LU(i,i-1+this%n) = alp1_d15
        endif
        ! i
        LU(i,i  ) = alp0_d15
        ! i+1
        if (i+1 .le. this%n) then
           LU(i,i+1) = alp1_d15
        else
           LU(i,i+1-this%n) = alp1_d15
        endif
        ! i+2
        if (i+2 .le. this%n) then
           LU(i,i+2) = alp2_d15
        else
           LU(i,i+2-this%n) = alp2_d15
        endif
        ! i+3
        if (i+3 .le. this%n) then
           LU(i,i+3) = alp3_d15
        else
           LU(i,i+3-this%n) = alp3_d15
        endif
        ! i+4
        if (i+4 .le. this%n) then
           LU(i,i+4) = alp4_d15
        else
           LU(i,i+4-this%n) = alp4_d15
        endif
        ! i+5
        if (i+5 .le. this%n) then
           LU(i,i+5) = alp5_d15
        else
           LU(i,i+5-this%n) = alp5_d15
        endif
        ! i+6
        if (i+6 .le. this%n) then
           LU(i,i+6) = alp6_d15
        else
           LU(i,i+6-this%n) = alp6_d15
        endif
        ! i+7
        if (i+7 .le. this%n) then
           LU(i,i+7) = alp7_d15
        else
           LU(i,i+7-this%n) = alp7_d15
        endif
     enddo


     call dgetrf(this%n,this%n,LU,this%n,ipiv,info)
     
 end subroutine

 subroutine ComputeD15(this,d15,ipiv,bc1,bcn)
   class (cfo2D15), intent(in) :: this
   integer, intent(in) :: bc1,bcn
   real(rkind), dimension(this%n,this%n), intent(out) :: d15
   integer, dimension(this%n), intent(out) :: ipiv
   real(rkind), dimension(this%n,this%n) :: d15o
   integer             :: i,info,j
   
   d15 = zero
   d15o = zero
   
   ! interior points
   do i = 8,this%n-7
      d15o(i,i-7) = alp7_d15
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15
      d15o(i,i+2) = alp2_d15
      d15o(i,i+3) = alp3_d15
      d15o(i,i+4) = alp4_d15
      d15o(i,i+5) = alp5_d15
      d15o(i,i+6) = alp6_d15
      d15o(i,i+7) = alp7_d15
   enddo
   
   ! BC at 1
   select case(bc1)
   case(0)
      ! b7
      i = 7
      d15o(i,i-6) = b7_bet6_d15
      d15o(i,i-5) = b7_bet5_d15
      d15o(i,i-4) = b7_bet4_d15
      d15o(i,i-3) = b7_bet3_d15
      d15o(i,i-2) = b7_bet2_d15
      d15o(i,i-1) = b7_bet1_d15
      d15o(i,i  ) = b7_alp0_d15
      d15o(i,i+1) = b7_alp1_d15
      d15o(i,i+2) = b7_alp2_d15
      d15o(i,i+3) = b7_alp3_d15
      d15o(i,i+4) = b7_alp4_d15
      d15o(i,i+5) = b7_alp5_d15
      d15o(i,i+6) = b7_alp6_d15
      d15o(i,i+7) = b7_alp7_d15
      ! b6
      i = 6
      d15o(i,i-5) = b6_bet5_d15
      d15o(i,i-4) = b6_bet4_d15
      d15o(i,i-3) = b6_bet3_d15
      d15o(i,i-2) = b6_bet2_d15
      d15o(i,i-1) = b6_bet1_d15
      d15o(i,i  ) = b6_alp0_d15
      d15o(i,i+1) = b6_alp1_d15
      d15o(i,i+2) = b6_alp2_d15
      d15o(i,i+3) = b6_alp3_d15
      d15o(i,i+4) = b6_alp4_d15
      d15o(i,i+5) = b6_alp5_d15
      d15o(i,i+6) = b6_alp6_d15
      d15o(i,i+7) = b6_alp7_d15
      ! b5
      i = 5
      d15o(i,i-4) = b5_bet4_d15
      d15o(i,i-3) = b5_bet3_d15
      d15o(i,i-2) = b5_bet2_d15
      d15o(i,i-1) = b5_bet1_d15
      d15o(i,i  ) = b5_alp0_d15
      d15o(i,i+1) = b5_alp1_d15
      d15o(i,i+2) = b5_alp2_d15
      d15o(i,i+3) = b5_alp3_d15
      d15o(i,i+4) = b5_alp4_d15
      d15o(i,i+5) = b5_alp5_d15
      d15o(i,i+6) = b5_alp6_d15
      d15o(i,i+7) = b5_alp7_d15
      ! b4
      i = 4
      d15o(i,i-3) = b4_bet3_d15
      d15o(i,i-2) = b4_bet2_d15
      d15o(i,i-1) = b4_bet1_d15
      d15o(i,i  ) = b4_alp0_d15
      d15o(i,i+1) = b4_alp1_d15
      d15o(i,i+2) = b4_alp2_d15
      d15o(i,i+3) = b4_alp3_d15
      d15o(i,i+4) = b4_alp4_d15
      d15o(i,i+5) = b4_alp5_d15
      d15o(i,i+6) = b4_alp6_d15
      d15o(i,i+7) = b4_alp7_d15
      ! b3
      i = 3
      d15o(i,i-2) = b3_bet2_d15
      d15o(i,i-1) = b3_bet1_d15
      d15o(i,i  ) = b3_alp0_d15
      d15o(i,i+1) = b3_alp1_d15
      d15o(i,i+2) = b3_alp2_d15
      d15o(i,i+3) = b3_alp3_d15
      d15o(i,i+4) = b3_alp4_d15
      d15o(i,i+5) = b3_alp5_d15
      d15o(i,i+6) = b3_alp6_d15
      d15o(i,i+7) = b3_alp7_d15
      ! b2
      i = 2
      d15o(i,i-1) = b2_bet1_d15
      d15o(i,i  ) = b2_alp0_d15
      d15o(i,i+1) = b2_alp1_d15
      d15o(i,i+2) = b2_alp2_d15
      d15o(i,i+3) = b2_alp3_d15
      d15o(i,i+4) = b2_alp4_d15
      d15o(i,i+5) = b2_alp5_d15
      d15o(i,i+6) = b2_alp6_d15
      d15o(i,i+7) = b2_alp7_d15
      ! b1
      i = 1
      d15o(i,i  ) = b1_alp0_d15
   case(1)
      !b7
      i = 7
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15 + alp7_d15
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15
      d15o(i,i+2) = alp2_d15
      d15o(i,i+3) = alp3_d15
      d15o(i,i+4) = alp4_d15
      d15o(i,i+5) = alp5_d15
      d15o(i,i+6) = alp6_d15
      d15o(i,i+7) = alp7_d15
      !b6
      i = 6
      d15o(i,i-5) = alp5_d15
      d15o(i,i-4) = alp4_d15 + alp6_d15
      d15o(i,i-3) = alp3_d15 + alp7_d15
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15
      d15o(i,i+2) = alp2_d15
      d15o(i,i+3) = alp3_d15
      d15o(i,i+4) = alp4_d15
      d15o(i,i+5) = alp5_d15
      d15o(i,i+6) = alp6_d15
      d15o(i,i+7) = alp7_d15
      !b5
      i = 5
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15 + alp5_d15
      d15o(i,i-2) = alp2_d15 + alp6_d15
      d15o(i,i-1) = alp1_d15 + alp7_d15
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15
      d15o(i,i+2) = alp2_d15
      d15o(i,i+3) = alp3_d15
      d15o(i,i+4) = alp4_d15
      d15o(i,i+5) = alp5_d15
      d15o(i,i+6) = alp6_d15
      d15o(i,i+7) = alp7_d15
      !b4
      i = 4
      d15o(i,i-3) = alp3_d15 
      d15o(i,i-2) = alp2_d15 + alp4_d15
      d15o(i,i-1) = alp1_d15 + alp5_d15
      d15o(i,i  ) = alp0_d15 + alp6_d15
      d15o(i,i+1) = alp1_d15 + alp7_d15
      d15o(i,i+2) = alp2_d15
      d15o(i,i+3) = alp3_d15
      d15o(i,i+4) = alp4_d15
      d15o(i,i+5) = alp5_d15
      d15o(i,i+6) = alp6_d15
      d15o(i,i+7) = alp7_d15
      !b3
      i = 3
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15 + alp3_d15
      d15o(i,i  ) = alp0_d15 + alp4_d15
      d15o(i,i+1) = alp1_d15 + alp5_d15
      d15o(i,i+2) = alp2_d15 + alp6_d15
      d15o(i,i+3) = alp3_d15 + alp7_d15
      d15o(i,i+4) = alp4_d15
      d15o(i,i+5) = alp5_d15
      d15o(i,i+6) = alp6_d15
      d15o(i,i+7) = alp7_d15
      !b2
      i = 2
      d15o(i,i-1) = alp1_d15
      d15o(i,i  ) = alp0_d15 + alp2_d15
      d15o(i,i+1) = alp1_d15 + alp3_d15
      d15o(i,i+2) = alp2_d15 + alp4_d15
      d15o(i,i+3) = alp3_d15 + alp5_d15
      d15o(i,i+4) = alp4_d15 + alp6_d15
      d15o(i,i+5) = alp5_d15 + alp7_d15
      d15o(i,i+6) = alp6_d15
      d15o(i,i+7) = alp7_d15
      !b1
      i = 1
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15 + alp1_d15
      d15o(i,i+2) = alp2_d15 + alp2_d15
      d15o(i,i+3) = alp3_d15 + alp3_d15
      d15o(i,i+4) = alp4_d15 + alp4_d15
      d15o(i,i+5) = alp5_d15 + alp5_d15
      d15o(i,i+6) = alp6_d15 + alp6_d15
      d15o(i,i+7) = alp7_d15 + alp7_d15
   case(-1)
      !b7
      i = 7
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15 - alp7_d15
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15
      d15o(i,i+2) = alp2_d15
      d15o(i,i+3) = alp3_d15
      d15o(i,i+4) = alp4_d15
      d15o(i,i+5) = alp5_d15
      d15o(i,i+6) = alp6_d15
      d15o(i,i+7) = alp7_d15
      !b6
      i = 6
      d15o(i,i-5) = alp5_d15
      d15o(i,i-4) = alp4_d15 - alp6_d15
      d15o(i,i-3) = alp3_d15 - alp7_d15
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15
      d15o(i,i+2) = alp2_d15
      d15o(i,i+3) = alp3_d15
      d15o(i,i+4) = alp4_d15
      d15o(i,i+5) = alp5_d15
      d15o(i,i+6) = alp6_d15
      d15o(i,i+7) = alp7_d15
      !b5
      i = 5
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15 - alp5_d15
      d15o(i,i-2) = alp2_d15 - alp6_d15
      d15o(i,i-1) = alp1_d15 - alp7_d15
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15
      d15o(i,i+2) = alp2_d15
      d15o(i,i+3) = alp3_d15
      d15o(i,i+4) = alp4_d15
      d15o(i,i+5) = alp5_d15
      d15o(i,i+6) = alp6_d15
      d15o(i,i+7) = alp7_d15
      !b4
      i = 4
      d15o(i,i-3) = alp3_d15 
      d15o(i,i-2) = alp2_d15 - alp4_d15
      d15o(i,i-1) = alp1_d15 - alp5_d15
      d15o(i,i  ) = alp0_d15 - alp6_d15
      d15o(i,i+1) = alp1_d15 - alp7_d15
      d15o(i,i+2) = alp2_d15
      d15o(i,i+3) = alp3_d15
      d15o(i,i+4) = alp4_d15
      d15o(i,i+5) = alp5_d15
      d15o(i,i+6) = alp6_d15
      d15o(i,i+7) = alp7_d15
      !b3
      i = 3
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15 - alp3_d15
      d15o(i,i  ) = alp0_d15 - alp4_d15
      d15o(i,i+1) = alp1_d15 - alp5_d15
      d15o(i,i+2) = alp2_d15 - alp6_d15
      d15o(i,i+3) = alp3_d15 - alp7_d15
      d15o(i,i+4) = alp4_d15
      d15o(i,i+5) = alp5_d15
      d15o(i,i+6) = alp6_d15
      d15o(i,i+7) = alp7_d15
      !b2
      i = 2
      d15o(i,i-1) = alp1_d15
      d15o(i,i  ) = alp0_d15 - alp2_d15
      d15o(i,i+1) = alp1_d15 - alp3_d15
      d15o(i,i+2) = alp2_d15 - alp4_d15
      d15o(i,i+3) = alp3_d15 - alp5_d15
      d15o(i,i+4) = alp4_d15 - alp6_d15
      d15o(i,i+5) = alp5_d15 - alp7_d15
      d15o(i,i+6) = alp6_d15
      d15o(i,i+7) = alp7_d15
      !b1
      i = 1
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15 - alp1_d15
      d15o(i,i+2) = alp2_d15 - alp2_d15
      d15o(i,i+3) = alp3_d15 - alp3_d15
      d15o(i,i+4) = alp4_d15 - alp4_d15
      d15o(i,i+5) = alp5_d15 - alp5_d15
      d15o(i,i+6) = alp6_d15 - alp6_d15
      d15o(i,i+7) = alp7_d15 - alp7_d15
   end select

   ! BC at n    
   select case(bcn)
   case(0)
      !bn-6
      i = this%n-6
      d15o(i,i-7) = b7_alp7_d15
      d15o(i,i-6) = b7_alp6_d15
      d15o(i,i-5) = b7_alp5_d15
      d15o(i,i-4) = b7_alp4_d15
      d15o(i,i-3) = b7_alp3_d15
      d15o(i,i-2) = b7_alp2_d15
      d15o(i,i-1) = b7_alp1_d15
      d15o(i,i  ) = b7_alp0_d15
      d15o(i,i+1) = b7_bet1_d15
      d15o(i,i+2) = b7_bet2_d15
      d15o(i,i+3) = b7_bet3_d15
      d15o(i,i+4) = b7_bet4_d15
      d15o(i,i+5) = b7_bet5_d15
      d15o(i,i+6) = b7_bet6_d15
      !bn-5
      i = this%n-5
      d15o(i,i-7) = b6_alp7_d15
      d15o(i,i-6) = b6_alp6_d15
      d15o(i,i-5) = b6_alp5_d15
      d15o(i,i-4) = b6_alp4_d15
      d15o(i,i-3) = b6_alp3_d15
      d15o(i,i-2) = b6_alp2_d15
      d15o(i,i-1) = b6_alp1_d15
      d15o(i,i  ) = b6_alp0_d15
      d15o(i,i+1) = b6_bet1_d15
      d15o(i,i+2) = b6_bet2_d15
      d15o(i,i+3) = b6_bet3_d15
      d15o(i,i+4) = b6_bet4_d15
      d15o(i,i+5) = b6_bet5_d15
      !bn-4
      i = this%n-4
      d15o(i,i-7) = b5_alp7_d15
      d15o(i,i-6) = b5_alp6_d15
      d15o(i,i-5) = b5_alp5_d15
      d15o(i,i-4) = b5_alp4_d15
      d15o(i,i-3) = b5_alp3_d15
      d15o(i,i-2) = b5_alp2_d15
      d15o(i,i-1) = b5_alp1_d15
      d15o(i,i  ) = b5_alp0_d15
      d15o(i,i+1) = b5_bet1_d15
      d15o(i,i+2) = b5_bet2_d15
      d15o(i,i+3) = b5_bet3_d15
      d15o(i,i+4) = b5_bet4_d15
      !bn-3
      i = this%n-3
      d15o(i,i-7) = b4_alp7_d15
      d15o(i,i-6) = b4_alp6_d15
      d15o(i,i-5) = b4_alp5_d15
      d15o(i,i-4) = b4_alp4_d15
      d15o(i,i-3) = b4_alp3_d15
      d15o(i,i-2) = b4_alp2_d15
      d15o(i,i-1) = b4_alp1_d15
      d15o(i,i  ) = b4_alp0_d15
      d15o(i,i+1) = b4_bet1_d15
      d15o(i,i+2) = b4_bet2_d15
      d15o(i,i+3) = b4_bet3_d15
      !bn-2
      i = this%n-2
      d15o(i,i-7) = b3_alp7_d15
      d15o(i,i-6) = b3_alp6_d15
      d15o(i,i-5) = b3_alp5_d15
      d15o(i,i-4) = b3_alp4_d15
      d15o(i,i-3) = b3_alp3_d15
      d15o(i,i-2) = b3_alp2_d15
      d15o(i,i-1) = b3_alp1_d15
      d15o(i,i  ) = b3_alp0_d15
      d15o(i,i+1) = b3_bet1_d15
      d15o(i,i+2) = b3_bet2_d15
      !bn-1
      i = this%n-1
      d15o(i,i-7) = b2_alp7_d15
      d15o(i,i-6) = b2_alp6_d15
      d15o(i,i-5) = b2_alp5_d15
      d15o(i,i-4) = b2_alp4_d15
      d15o(i,i-3) = b2_alp3_d15
      d15o(i,i-2) = b2_alp2_d15
      d15o(i,i-1) = b2_alp1_d15
      d15o(i,i  ) = b2_alp0_d15
      d15o(i,i+1) = b2_bet1_d15
      !bn
      i = this%n
      d15o(i,i  ) = b1_alp0_d15
   case(1)
      !bn-6
      i = this%n-6
      d15o(i,i-7) = alp7_d15
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15
      d15o(i,i+2) = alp2_d15
      d15o(i,i+3) = alp3_d15
      d15o(i,i+4) = alp4_d15
      d15o(i,i+5) = alp5_d15 + alp7_d15
      d15o(i,i+6) = alp6_d15
      !bn-5
      i = this%n-5
      d15o(i,i-7) = alp7_d15
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15
      d15o(i,i+2) = alp2_d15
      d15o(i,i+3) = alp3_d15 + alp7_d15
      d15o(i,i+4) = alp4_d15 + alp6_d15
      d15o(i,i+5) = alp5_d15
      !bn-4
      i = this%n-4
      d15o(i,i-7) = alp7_d15
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15 + alp7_d15
      d15o(i,i+2) = alp2_d15 + alp6_d15
      d15o(i,i+3) = alp3_d15 + alp5_d15
      d15o(i,i+4) = alp4_d15
      !bn-3
      i = this%n-3
      d15o(i,i-7) = alp7_d15
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15 + alp7_d15
      d15o(i,i  ) = alp0_d15 + alp6_d15
      d15o(i,i+1) = alp1_d15 + alp5_d15
      d15o(i,i+2) = alp2_d15 + alp4_d15
      d15o(i,i+3) = alp3_d15
      !bn-2
      i = this%n-2
      d15o(i,i-7) = alp7_d15
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15 + alp7_d15
      d15o(i,i-2) = alp2_d15 + alp6_d15
      d15o(i,i-1) = alp1_d15 + alp5_d15
      d15o(i,i  ) = alp0_d15 + alp4_d15
      d15o(i,i+1) = alp1_d15 + alp3_d15
      d15o(i,i+2) = alp2_d15
      !bn-1
      i = this%n-1
      d15o(i,i-7) = alp7_d15
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15 + alp7_d15
      d15o(i,i-4) = alp4_d15 + alp6_d15
      d15o(i,i-3) = alp3_d15 + alp5_d15
      d15o(i,i-2) = alp2_d15 + alp4_d15
      d15o(i,i-1) = alp1_d15 + alp3_d15
      d15o(i,i  ) = alp0_d15 + alp2_d15
      d15o(i,i+1) = alp1_d15
      !bn
      i = this%n
      d15o(i,i-7) = alp7_d15 + alp7_d15
      d15o(i,i-6) = alp6_d15 + alp6_d15
      d15o(i,i-5) = alp5_d15 + alp5_d15
      d15o(i,i-4) = alp4_d15 + alp4_d15
      d15o(i,i-3) = alp3_d15 + alp3_d15
      d15o(i,i-2) = alp2_d15 + alp2_d15
      d15o(i,i-1) = alp1_d15 + alp1_d15
      d15o(i,i  ) = alp0_d15
   case(-1)
      !bn-6
      i = this%n-6
      d15o(i,i-7) = alp7_d15
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15
      d15o(i,i+2) = alp2_d15
      d15o(i,i+3) = alp3_d15
      d15o(i,i+4) = alp4_d15
      d15o(i,i+5) = alp5_d15 - alp7_d15
      d15o(i,i+6) = alp6_d15
      !bn-5
      i = this%n-5
      d15o(i,i-7) = alp7_d15
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15
      d15o(i,i+2) = alp2_d15
      d15o(i,i+3) = alp3_d15 - alp7_d15
      d15o(i,i+4) = alp4_d15 - alp6_d15
      d15o(i,i+5) = alp5_d15
      !bn-4
      i = this%n-4
      d15o(i,i-7) = alp7_d15
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15
      d15o(i,i  ) = alp0_d15
      d15o(i,i+1) = alp1_d15 - alp7_d15
      d15o(i,i+2) = alp2_d15 - alp6_d15
      d15o(i,i+3) = alp3_d15 - alp5_d15
      d15o(i,i+4) = alp4_d15
      !bn-3
      i = this%n-3
      d15o(i,i-7) = alp7_d15
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15
      d15o(i,i-2) = alp2_d15
      d15o(i,i-1) = alp1_d15 - alp7_d15
      d15o(i,i  ) = alp0_d15 - alp6_d15
      d15o(i,i+1) = alp1_d15 - alp5_d15
      d15o(i,i+2) = alp2_d15 - alp4_d15
      d15o(i,i+3) = alp3_d15
      !bn-2
      i = this%n-2
      d15o(i,i-7) = alp7_d15
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15
      d15o(i,i-4) = alp4_d15
      d15o(i,i-3) = alp3_d15 - alp7_d15
      d15o(i,i-2) = alp2_d15 - alp6_d15
      d15o(i,i-1) = alp1_d15 - alp5_d15
      d15o(i,i  ) = alp0_d15 - alp4_d15
      d15o(i,i+1) = alp1_d15 - alp3_d15
      d15o(i,i+2) = alp2_d15
      !bn-1
      i = this%n-1
      d15o(i,i-7) = alp7_d15
      d15o(i,i-6) = alp6_d15
      d15o(i,i-5) = alp5_d15 - alp7_d15
      d15o(i,i-4) = alp4_d15 - alp6_d15
      d15o(i,i-3) = alp3_d15 - alp5_d15
      d15o(i,i-2) = alp2_d15 - alp4_d15
      d15o(i,i-1) = alp1_d15 - alp3_d15
      d15o(i,i  ) = alp0_d15 - alp2_d15
      d15o(i,i+1) = alp1_d15
      !bn
      i = this%n
      d15o(i,i-7) = alp7_d15 - alp7_d15
      d15o(i,i-6) = alp6_d15 - alp6_d15
      d15o(i,i-5) = alp5_d15 - alp5_d15
      d15o(i,i-4) = alp4_d15 - alp4_d15
      d15o(i,i-3) = alp3_d15 - alp3_d15
      d15o(i,i-2) = alp2_d15 - alp2_d15
      d15o(i,i-1) = alp1_d15 - alp1_d15
      d15o(i,i  ) = alp0_d15
   end select

   !Rearrange d15 for dgbtrf banded form
   do j = 1,this%n
      do i = max(1,j-7),min(this%n,j+7)
         d15(15+i-j,j) = d15o(i,j)
      enddo
   enddo
   !end

   call dgbtrf(this%n,this%n,7,7,d15,this%n,ipiv,info)
   
 end subroutine ComputeD15

 subroutine SolveXLU(this,y,n2,n3)

   class( cfo2D15 ), intent(in) :: this
   integer, intent(in) :: n2,n3
   real(rkind), dimension(this%n,n2,n3), intent(inout) :: y  ! Take in RHS and put solution into it
   integer :: j,k,info

   do k=1,n3
      do j=1,n2
         call dgetrs('N',this%n,1,this%LU,this%n,this%LU_ipiv,y(:,j,k),this%n,info)
      end do
   end do

 end subroutine SolveXLU

 subroutine SolveYLU(this,y,n1,n3)

   class( cfo2D15 ), intent(in) :: this
   integer, intent(in) :: n1,n3
   real(rkind), dimension(n1,this%n,n3), intent(inout) :: y  ! Take in RHS and put solution into it
   integer :: i,k,info

   do k=1,n3
      do i=1,n1
         call dgetrs('N',this%n,1,this%LU,this%n,this%LU_ipiv,y(i,:,k),this%n,info)
      end do
   end do

 end subroutine SolveYLU

 subroutine SolveZLU(this,y,n1,n2)

   class( cfo2D15 ), intent(in) :: this
   integer, intent(in) :: n1,n2
   real(rkind), dimension(n1,n2,this%n), intent(inout) :: y  ! Take in RHS and put solution into it
   integer :: j,i,info

   do j=1,n2
      do i=1,n1
         call dgetrs('N',this%n,1,this%LU,this%n,this%LU_ipiv,y(i,j,:),this%n,info)
      end do
   end do

 end subroutine SolveZLU

 subroutine SolveXD15(this,d15,ipiv,y,n2,n3)

   class( cfo2D15 ), intent(in) :: this
   real(rkind), dimension(this%n,this%n), intent(in) :: d15
   integer, dimension(this%n), intent(in) :: ipiv
   integer, intent(in) :: n2,n3
   real(rkind), dimension(this%n,n2,n3), intent(inout) :: y
   integer :: j, k, info

   do k=1,n3
      do j=1,n2
         call dgbtrs('N',this%n,7,7,1,d15,this%n,ipiv,y(:,j,k),this%n,info)
      enddo
   enddo

 end subroutine SolveXD15

 subroutine SolveYD15(this,d15,ipiv,y,n1,n3)

   class( cfo2D15 ), intent(in) :: this
   real(rkind), dimension(this%n,this%n), intent(in) :: d15
   integer, dimension(this%n), intent(in) :: ipiv
   integer, intent(in) :: n1,n3
   real(rkind), dimension(n1,this%n,n3), intent(inout) :: y
   integer :: i, k, info

   do k=1,n3
      do i=1,n1
         call dgbtrs('N',this%n,7,7,1,d15,this%n,ipiv,y(i,:,k),this%n,info)
      enddo
   enddo

 end subroutine SolveYD15

 subroutine SolveZD15_REAL(this,d15,ipiv,y,n1,n2)

   class( cfo2D15 ), intent(in) :: this
   real(rkind), dimension(this%n,this%n), intent(in) :: d15
   integer, dimension(this%n), intent(in) :: ipiv
   integer, intent(in) :: n1,n2
   real(rkind), dimension(n1,n2,this%n), intent(inout) :: y
   integer :: i,j,info

#include "CFo2D15_files/SolveZD15_common.F90"        

 end subroutine SolveZD15_REAL

 subroutine SolveZD15_CMPLX(this,d15,ipiv,y,n1,n2)

   class( cfo2D15 ), intent(in) :: this
   real(rkind), dimension(this%n,this%n), intent(in) :: d15
   integer, dimension(this%n), intent(in) :: ipiv
   integer, intent(in) :: n1,n2
   complex(rkind), dimension(n1,n2,this%n), intent(inout) :: y
   integer :: i,j,info

#include "CFo2D15_files/SolveZD15_common.F90"        

 end subroutine SolveZD15_CMPLX


 pure subroutine ComputeXRHS(this, f, RHS, n2, n3, bc1, bcn)
 
   class( cfo2D15 ), intent(in) :: this
   integer, intent(in) :: n2, n3
   real(rkind), dimension(this%n,n2,n3), intent(in) :: f
   real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
   integer, intent(in) :: bc1, bcn
   integer :: j,k

   select case (this%periodic)
   case (.TRUE.)
      do k=1,n3
         do j=1,n2
            RHS(         1,j,k) = a0_d15 * ( f(         1,j,k) )       &
                 + a1_d15 * ( f(         2 ,j,k) + f(    this%n,j,k) ) &
                 + a2_d15 * ( f(         3 ,j,k) + f(  this%n-1,j,k) ) &
                 + a3_d15 * ( f(         4 ,j,k) + f(  this%n-2,j,k) ) &
                 + a4_d15 * ( f(         5 ,j,k) + f(  this%n-3,j,k) ) &
                 + a5_d15 * ( f(         6 ,j,k) + f(  this%n-4,j,k) ) &
                 + a6_d15 * ( f(         7 ,j,k) + f(  this%n-5,j,k) ) &
                 + a7_d15 * ( f(         8 ,j,k) + f(  this%n-6,j,k) )
            RHS(         2,j,k) = a0_d15 * ( f(         2,j,k) )       &
                 + a1_d15 * ( f(         3 ,j,k) + f(         1,j,k) ) &
                 + a2_d15 * ( f(         4 ,j,k) + f(    this%n,j,k) ) &
                 + a3_d15 * ( f(         5 ,j,k) + f(  this%n-1,j,k) ) &
                 + a4_d15 * ( f(         6 ,j,k) + f(  this%n-2,j,k) ) &
                 + a5_d15 * ( f(         7 ,j,k) + f(  this%n-3,j,k) ) &
                 + a6_d15 * ( f(         8 ,j,k) + f(  this%n-4,j,k) ) &
                 + a7_d15 * ( f(         9 ,j,k) + f(  this%n-5,j,k) )
            RHS(         3,j,k) = a0_d15 * ( f(         3,j,k) )       &
                 + a1_d15 * ( f(         4 ,j,k) + f(         2,j,k) ) &
                 + a2_d15 * ( f(         5 ,j,k) + f(         1,j,k) ) &
                 + a3_d15 * ( f(         6 ,j,k) + f(    this%n,j,k) ) &
                 + a4_d15 * ( f(         7 ,j,k) + f(  this%n-1,j,k) ) &
                 + a5_d15 * ( f(         8 ,j,k) + f(  this%n-2,j,k) ) &
                 + a6_d15 * ( f(         9 ,j,k) + f(  this%n-3,j,k) ) &
                 + a7_d15 * ( f(         10,j,k) + f(  this%n-4,j,k) )
            RHS(         4,j,k) = a0_d15 * ( f(         4,j,k) )       &
                 + a1_d15 * ( f(         5 ,j,k) + f(         3,j,k) ) &
                 + a2_d15 * ( f(         6 ,j,k) + f(         2,j,k) ) &
                 + a3_d15 * ( f(         7 ,j,k) + f(         1,j,k) ) &
                 + a4_d15 * ( f(         8 ,j,k) + f(    this%n,j,k) ) &
                 + a5_d15 * ( f(         9 ,j,k) + f(  this%n-1,j,k) ) &
                 + a6_d15 * ( f(         10,j,k) + f(  this%n-2,j,k) ) &
                 + a7_d15 * ( f(         11,j,k) + f(  this%n-3,j,k) )
            RHS(         5,j,k) = a0_d15 * ( f(         5,j,k) )       &
                 + a1_d15 * ( f(         6 ,j,k) + f(         4,j,k) ) &
                 + a2_d15 * ( f(         7 ,j,k) + f(         3,j,k) ) &
                 + a3_d15 * ( f(         8 ,j,k) + f(         2,j,k) ) &
                 + a4_d15 * ( f(         9 ,j,k) + f(         1,j,k) ) &
                 + a5_d15 * ( f(         10,j,k) + f(    this%n,j,k) ) &
                 + a6_d15 * ( f(         11,j,k) + f(  this%n-1,j,k) ) &
                 + a7_d15 * ( f(         12,j,k) + f(  this%n-2,j,k) )
            RHS(         6,j,k) = a0_d15 * ( f(         6,j,k) )       &
                 + a1_d15 * ( f(         7 ,j,k) + f(         5,j,k) ) &
                 + a2_d15 * ( f(         8 ,j,k) + f(         4,j,k) ) &
                 + a3_d15 * ( f(         9 ,j,k) + f(         3,j,k) ) &
                 + a4_d15 * ( f(         10,j,k) + f(         2,j,k) ) &
                 + a5_d15 * ( f(         11,j,k) + f(         1,j,k) ) &
                 + a6_d15 * ( f(         12,j,k) + f(    this%n,j,k) ) &
                 + a7_d15 * ( f(         13,j,k) + f(  this%n-1,j,k) )
            RHS(         7,j,k) = a0_d15 * ( f(         7,j,k) )       &
                 + a1_d15 * ( f(         8 ,j,k) + f(         6,j,k) ) &
                 + a2_d15 * ( f(         9 ,j,k) + f(         5,j,k) ) &
                 + a3_d15 * ( f(         10,j,k) + f(         4,j,k) ) &
                 + a4_d15 * ( f(         11,j,k) + f(         3,j,k) ) &
                 + a5_d15 * ( f(         12,j,k) + f(         2,j,k) ) &
                 + a6_d15 * ( f(         13,j,k) + f(         1,j,k) ) &
                 + a7_d15 * ( f(         14,j,k) + f(    this%n,j,k) )
            RHS(8:this%n-7,j,k) = a0_d15 * ( f(8:this%n-7,j,k) )        &
                 + a1_d15 * ( f(9 :this%n-6,j,k) + f(7:this%n-8 ,j,k) ) &
                 + a2_d15 * ( f(10:this%n-5,j,k) + f(6:this%n-9 ,j,k) ) &
                 + a3_d15 * ( f(11:this%n-4,j,k) + f(5:this%n-10,j,k) ) &
                 + a4_d15 * ( f(12:this%n-3,j,k) + f(4:this%n-11,j,k) ) &
                 + a5_d15 * ( f(13:this%n-2,j,k) + f(3:this%n-12,j,k) ) &
                 + a6_d15 * ( f(14:this%n-1,j,k) + f(2:this%n-13,j,k) ) &
                 + a7_d15 * ( f(15:this%n  ,j,k) + f(1:this%n-14,j,k) )
            RHS(  this%n-6,j,k) = a0_d15 * ( f(  this%n-6,j,k) )        &
                 + a1_d15 * ( f(   this%n-5,j,k) + f(  this%n-7 ,j,k) ) &
                 + a2_d15 * ( f(   this%n-4,j,k) + f(  this%n-8 ,j,k) ) &
                 + a3_d15 * ( f(   this%n-3,j,k) + f(  this%n-9 ,j,k) ) &
                 + a4_d15 * ( f(   this%n-2,j,k) + f(  this%n-10,j,k) ) &
                 + a5_d15 * ( f(   this%n-1,j,k) + f(  this%n-11,j,k) ) &
                 + a6_d15 * ( f(   this%n  ,j,k) + f(  this%n-12,j,k) ) &
                 + a7_d15 * ( f(          1,j,k) + f(  this%n-13,j,k) )
            RHS(  this%n-5,j,k) = a0_d15 * ( f(  this%n-5,j,k) )        &
                 + a1_d15 * ( f(   this%n-4,j,k) + f(  this%n-6 ,j,k) ) &
                 + a2_d15 * ( f(   this%n-3,j,k) + f(  this%n-7 ,j,k) ) &
                 + a3_d15 * ( f(   this%n-2,j,k) + f(  this%n-8 ,j,k) ) &
                 + a4_d15 * ( f(   this%n-1,j,k) + f(  this%n-9 ,j,k) ) &
                 + a5_d15 * ( f(   this%n  ,j,k) + f(  this%n-10,j,k) ) &
                 + a6_d15 * ( f(          1,j,k) + f(  this%n-11,j,k) ) &
                 + a7_d15 * ( f(          2,j,k) + f(  this%n-12,j,k) )
            RHS(  this%n-4,j,k) = a0_d15 * ( f(  this%n-4,j,k) )        &
                 + a1_d15 * ( f(   this%n-3,j,k) + f(  this%n-5 ,j,k) ) &
                 + a2_d15 * ( f(   this%n-2,j,k) + f(  this%n-6 ,j,k) ) &
                 + a3_d15 * ( f(   this%n-1,j,k) + f(  this%n-7 ,j,k) ) &
                 + a4_d15 * ( f(   this%n  ,j,k) + f(  this%n-8 ,j,k) ) &
                 + a5_d15 * ( f(          1,j,k) + f(  this%n-9 ,j,k) ) &
                 + a6_d15 * ( f(          2,j,k) + f(  this%n-10,j,k) ) &
                 + a7_d15 * ( f(          3,j,k) + f(  this%n-11,j,k) )
            RHS(  this%n-3,j,k) = a0_d15 * ( f(  this%n-3,j,k) )        &
                 + a1_d15 * ( f(   this%n-2,j,k) + f(  this%n-4 ,j,k) ) &
                 + a2_d15 * ( f(   this%n-1,j,k) + f(  this%n-5 ,j,k) ) &
                 + a3_d15 * ( f(   this%n  ,j,k) + f(  this%n-6 ,j,k) ) &
                 + a4_d15 * ( f(          1,j,k) + f(  this%n-7 ,j,k) ) &
                 + a5_d15 * ( f(          2,j,k) + f(  this%n-8 ,j,k) ) &
                 + a6_d15 * ( f(          3,j,k) + f(  this%n-9 ,j,k) ) &
                 + a7_d15 * ( f(          4,j,k) + f(  this%n-10,j,k) )
            RHS(  this%n-2,j,k) = a0_d15 * ( f(  this%n-2,j,k) )        &
                 + a1_d15 * ( f(   this%n-1,j,k) + f(  this%n-3 ,j,k) ) &
                 + a2_d15 * ( f(   this%n  ,j,k) + f(  this%n-4 ,j,k) ) &
                 + a3_d15 * ( f(          1,j,k) + f(  this%n-5 ,j,k) ) &
                 + a4_d15 * ( f(          2,j,k) + f(  this%n-6 ,j,k) ) &
                 + a5_d15 * ( f(          3,j,k) + f(  this%n-7 ,j,k) ) &
                 + a6_d15 * ( f(          4,j,k) + f(  this%n-8 ,j,k) ) &
                 + a7_d15 * ( f(          5,j,k) + f(  this%n-9 ,j,k) )
            RHS(  this%n-1,j,k) = a0_d15 * ( f(  this%n-1,j,k) )        &
                 + a1_d15 * ( f(   this%n  ,j,k) + f(  this%n-2 ,j,k) ) &
                 + a2_d15 * ( f(          1,j,k) + f(  this%n-3 ,j,k) ) &
                 + a3_d15 * ( f(          2,j,k) + f(  this%n-4 ,j,k) ) &
                 + a4_d15 * ( f(          3,j,k) + f(  this%n-5 ,j,k) ) &
                 + a5_d15 * ( f(          4,j,k) + f(  this%n-6 ,j,k) ) &
                 + a6_d15 * ( f(          5,j,k) + f(  this%n-7 ,j,k) ) &
                 + a7_d15 * ( f(          6,j,k) + f(  this%n-8 ,j,k) )
            RHS(  this%n  ,j,k) = a0_d15 * ( f(  this%n  ,j,k) )        &
                 + a1_d15 * ( f(          1,j,k) + f(  this%n-1 ,j,k) ) &
                 + a2_d15 * ( f(          2,j,k) + f(  this%n-2 ,j,k) ) &
                 + a3_d15 * ( f(          3,j,k) + f(  this%n-3 ,j,k) ) &
                 + a4_d15 * ( f(          4,j,k) + f(  this%n-4 ,j,k) ) &
                 + a5_d15 * ( f(          5,j,k) + f(  this%n-5 ,j,k) ) &
                 + a6_d15 * ( f(          6,j,k) + f(  this%n-6 ,j,k) ) &
                 + a7_d15 * ( f(          7,j,k) + f(  this%n-7 ,j,k) )
         end do
      end do
   
    case (.FALSE.)
      
       ! interior points
       do k=1,n3
          do j=1,n2
             RHS(8:this%n-7,j,k) = a0_d15 * ( f(8:this%n-7,j,k) )        &
                  + a1_d15 * ( f(9 :this%n-6,j,k) + f(7:this%n-8 ,j,k) ) &
                  + a2_d15 * ( f(10:this%n-5,j,k) + f(6:this%n-9 ,j,k) ) &
                  + a3_d15 * ( f(11:this%n-4,j,k) + f(5:this%n-10,j,k) ) &
                  + a4_d15 * ( f(12:this%n-3,j,k) + f(4:this%n-11,j,k) ) &
                  + a5_d15 * ( f(13:this%n-2,j,k) + f(3:this%n-12,j,k) ) &
                  + a6_d15 * ( f(14:this%n-1,j,k) + f(2:this%n-13,j,k) ) &
                  + a7_d15 * ( f(15:this%n  ,j,k) + f(1:this%n-14,j,k) )
          enddo
       enddo
       
       ! bc 1
       select case(bc1)
       case(0)
          
          do k=1,n3
             do j=1,n2
                RHS(         1,j,k) = b1_a0_d15 * ( f(         1,j,k) )                
                RHS(         2,j,k) = b2_a0_d15 * ( f(         2,j,k) )                &
                     + b2_a1_d15 * f(         3 ,j,k) + b2_b1_d15 * f(         1,j,k)  &
                     + b2_a2_d15 * f(         4 ,j,k)                                  &
                     + b2_a3_d15 * f(         5 ,j,k)                                  &
                     + b2_a4_d15 * f(         6 ,j,k)                                  &
                     + b2_a5_d15 * f(         7 ,j,k)                                  &
                     + b2_a6_d15 * f(         8 ,j,k)                                  &
                     + b2_a7_d15 * f(         9 ,j,k)                     
                RHS(         3,j,k) = b3_a0_d15 * ( f(         3,j,k) )                &
                     + b3_a1_d15 * f(         4 ,j,k) + b3_b1_d15 * f(         2,j,k)  &
                     + b3_a2_d15 * f(         5 ,j,k) + b3_b2_d15 * f(         1,j,k)  &
                     + b3_a3_d15 * f(         6 ,j,k)                                  &
                     + b3_a4_d15 * f(         7 ,j,k)                                  &
                     + b3_a5_d15 * f(         8 ,j,k)                                  &
                     + b3_a6_d15 * f(         9 ,j,k)                                  &
                     + b3_a7_d15 * f(         10,j,k)                     
                RHS(         4,j,k) = b4_a0_d15 * ( f(         4,j,k) )                &
                     + b4_a1_d15 * f(         5 ,j,k) + b4_b1_d15 * f(         3,j,k)  &
                     + b4_a2_d15 * f(         6 ,j,k) + b4_b2_d15 * f(         2,j,k)  &
                     + b4_a3_d15 * f(         7 ,j,k) + b4_b3_d15 * f(         1,j,k)  &
                     + b4_a4_d15 * f(         8 ,j,k)                                  &
                     + b4_a5_d15 * f(         9 ,j,k)                                  &
                     + b4_a6_d15 * f(         10,j,k)                                  &
                     + b4_a7_d15 * f(         11,j,k)                     
                RHS(         5,j,k) = b5_a0_d15 * ( f(         5,j,k) )                &
                     + b5_a1_d15 * f(         6 ,j,k) + b5_b1_d15 * f(         4,j,k)  &
                     + b5_a2_d15 * f(         7 ,j,k) + b5_b2_d15 * f(         3,j,k)  &
                     + b5_a3_d15 * f(         8 ,j,k) + b5_b3_d15 * f(         2,j,k)  &
                     + b5_a4_d15 * f(         9 ,j,k) + b5_b4_d15 * f(         1,j,k)  &
                     + b5_a5_d15 * f(         10,j,k)                                  &
                     + b5_a6_d15 * f(         11,j,k)                                  &
                     + b5_a7_d15 * f(         12,j,k)                     
                RHS(         6,j,k) = b6_a0_d15 * ( f(         6,j,k) )                &
                     + b6_a1_d15 * f(         7 ,j,k) + b6_b1_d15 * f(         5,j,k)  &
                     + b6_a2_d15 * f(         8 ,j,k) + b6_b2_d15 * f(         4,j,k)  &
                     + b6_a3_d15 * f(         9 ,j,k) + b6_b3_d15 * f(         3,j,k)  &
                     + b6_a4_d15 * f(         10,j,k) + b6_b4_d15 * f(         2,j,k)  &
                     + b6_a5_d15 * f(         11,j,k) + b6_b5_d15 * f(         1,j,k)  &
                     + b6_a6_d15 * f(         12,j,k)                                  &
                     + b6_a7_d15 * f(         13,j,k)                     
                RHS(         7,j,k) = b7_a0_d15 * ( f(         7,j,k) )                &
                     + b7_a1_d15 * f(         8 ,j,k) + b7_b1_d15 * f(         6,j,k)  &
                     + b7_a2_d15 * f(         9 ,j,k) + b7_b2_d15 * f(         5,j,k)  &
                     + b7_a3_d15 * f(         10,j,k) + b7_b3_d15 * f(         4,j,k)  &
                     + b7_a4_d15 * f(         11,j,k) + b7_b4_d15 * f(         3,j,k)  &
                     + b7_a5_d15 * f(         12,j,k) + b7_b5_d15 * f(         2,j,k)  &
                     + b7_a6_d15 * f(         13,j,k) + b7_b6_d15 * f(         1,j,k)  &
                     + b7_a7_d15 * f(         14,j,k)                     
             enddo
          enddo

       case(1)
          do k=1,n3
             do j=1,n2
                RHS(         1,j,k) = a0_d15 * ( f(         1,j,k) )       &
                     + a1_d15 * ( f(         2 ,j,k) + f(         2,j,k) ) &
                     + a2_d15 * ( f(         3 ,j,k) + f(         3,j,k) ) &
                     + a3_d15 * ( f(         4 ,j,k) + f(         4,j,k) ) &
                     + a4_d15 * ( f(         5 ,j,k) + f(         5,j,k) ) &
                     + a5_d15 * ( f(         6 ,j,k) + f(         6,j,k) ) &
                     + a6_d15 * ( f(         7 ,j,k) + f(         7,j,k) ) &
                     + a7_d15 * ( f(         8 ,j,k) + f(         8,j,k) )
                RHS(         2,j,k) = a0_d15 * ( f(         2,j,k) )       &
                     + a1_d15 * ( f(         3 ,j,k) + f(         1,j,k) ) &
                     + a2_d15 * ( f(         4 ,j,k) + f(         2,j,k) ) &
                     + a3_d15 * ( f(         5 ,j,k) + f(         3,j,k) ) &
                     + a4_d15 * ( f(         6 ,j,k) + f(         4,j,k) ) &
                     + a5_d15 * ( f(         7 ,j,k) + f(         5,j,k) ) &
                     + a6_d15 * ( f(         8 ,j,k) + f(         6,j,k) ) &
                     + a7_d15 * ( f(         9 ,j,k) + f(         7,j,k) )
                RHS(         3,j,k) = a0_d15 * ( f(         3,j,k) )       &
                     + a1_d15 * ( f(         4 ,j,k) + f(         2,j,k) ) &
                     + a2_d15 * ( f(         5 ,j,k) + f(         1,j,k) ) &
                     + a3_d15 * ( f(         6 ,j,k) + f(         2,j,k) ) &
                     + a4_d15 * ( f(         7 ,j,k) + f(         3,j,k) ) &
                     + a5_d15 * ( f(         8 ,j,k) + f(         4,j,k) ) &
                     + a6_d15 * ( f(         9 ,j,k) + f(         5,j,k) ) &
                     + a7_d15 * ( f(         10,j,k) + f(         6,j,k) )
                RHS(         4,j,k) = a0_d15 * ( f(         4,j,k) )       &
                     + a1_d15 * ( f(         5 ,j,k) + f(         3,j,k) ) &
                     + a2_d15 * ( f(         6 ,j,k) + f(         2,j,k) ) &
                     + a3_d15 * ( f(         7 ,j,k) + f(         1,j,k) ) &
                     + a4_d15 * ( f(         8 ,j,k) + f(         2,j,k) ) &
                     + a5_d15 * ( f(         9 ,j,k) + f(         3,j,k) ) &
                     + a6_d15 * ( f(         10,j,k) + f(         4,j,k) ) &
                     + a7_d15 * ( f(         11,j,k) + f(         5,j,k) )
                RHS(         5,j,k) = a0_d15 * ( f(         5,j,k) )       &
                     + a1_d15 * ( f(         6 ,j,k) + f(         4,j,k) ) &
                     + a2_d15 * ( f(         7 ,j,k) + f(         3,j,k) ) &
                     + a3_d15 * ( f(         8 ,j,k) + f(         2,j,k) ) &
                     + a4_d15 * ( f(         9 ,j,k) + f(         1,j,k) ) &
                     + a5_d15 * ( f(         10,j,k) + f(         2,j,k) ) &
                     + a6_d15 * ( f(         11,j,k) + f(         3,j,k) ) &
                     + a7_d15 * ( f(         12,j,k) + f(         4,j,k) )
                RHS(         6,j,k) = a0_d15 * ( f(         6,j,k) )       &
                     + a1_d15 * ( f(         7 ,j,k) + f(         5,j,k) ) &
                     + a2_d15 * ( f(         8 ,j,k) + f(         4,j,k) ) &
                     + a3_d15 * ( f(         9 ,j,k) + f(         3,j,k) ) &
                     + a4_d15 * ( f(         10,j,k) + f(         2,j,k) ) &
                     + a5_d15 * ( f(         11,j,k) + f(         1,j,k) ) &
                     + a6_d15 * ( f(         12,j,k) + f(         2,j,k) ) &
                     + a7_d15 * ( f(         13,j,k) + f(         3,j,k) )
                RHS(         7,j,k) = a0_d15 * ( f(         7,j,k) )       &
                     + a1_d15 * ( f(         8 ,j,k) + f(         6,j,k) ) &
                     + a2_d15 * ( f(         9 ,j,k) + f(         5,j,k) ) &
                     + a3_d15 * ( f(         10,j,k) + f(         4,j,k) ) &
                     + a4_d15 * ( f(         11,j,k) + f(         3,j,k) ) &
                     + a5_d15 * ( f(         12,j,k) + f(         2,j,k) ) &
                     + a6_d15 * ( f(         13,j,k) + f(         1,j,k) ) &
                     + a7_d15 * ( f(         14,j,k) + f(         2,j,k) )
             enddo
          enddo

       case(-1)
          do k=1,n3
             do j=1,n2
                RHS(         1,j,k) = a0_d15 * ( f(         1,j,k) )       &
                     + a1_d15 * ( f(         2 ,j,k) - f(         2,j,k) ) &
                     + a2_d15 * ( f(         3 ,j,k) - f(         3,j,k) ) &
                     + a3_d15 * ( f(         4 ,j,k) - f(         4,j,k) ) &
                     + a4_d15 * ( f(         5 ,j,k) - f(         5,j,k) ) &
                     + a5_d15 * ( f(         6 ,j,k) - f(         6,j,k) ) &
                     + a6_d15 * ( f(         7 ,j,k) - f(         7,j,k) ) &
                     + a7_d15 * ( f(         8 ,j,k) - f(         8,j,k) )
                RHS(         2,j,k) = a0_d15 * ( f(         2,j,k) )       &
                     + a1_d15 * ( f(         3 ,j,k) + f(         1,j,k) ) &
                     + a2_d15 * ( f(         4 ,j,k) - f(         2,j,k) ) &
                     + a3_d15 * ( f(         5 ,j,k) - f(         3,j,k) ) &
                     + a4_d15 * ( f(         6 ,j,k) - f(         4,j,k) ) &
                     + a5_d15 * ( f(         7 ,j,k) - f(         5,j,k) ) &
                     + a6_d15 * ( f(         8 ,j,k) - f(         6,j,k) ) &
                     + a7_d15 * ( f(         9 ,j,k) - f(         7,j,k) )
                RHS(         3,j,k) = a0_d15 * ( f(         3,j,k) )       &
                     + a1_d15 * ( f(         4 ,j,k) + f(         2,j,k) ) &
                     + a2_d15 * ( f(         5 ,j,k) + f(         1,j,k) ) &
                     + a3_d15 * ( f(         6 ,j,k) - f(         2,j,k) ) &
                     + a4_d15 * ( f(         7 ,j,k) - f(         3,j,k) ) &
                     + a5_d15 * ( f(         8 ,j,k) - f(         4,j,k) ) &
                     + a6_d15 * ( f(         9 ,j,k) - f(         5,j,k) ) &
                     + a7_d15 * ( f(         10,j,k) - f(         6,j,k) )
                RHS(         4,j,k) = a0_d15 * ( f(         4,j,k) )       &
                     + a1_d15 * ( f(         5 ,j,k) + f(         3,j,k) ) &
                     + a2_d15 * ( f(         6 ,j,k) + f(         2,j,k) ) &
                     + a3_d15 * ( f(         7 ,j,k) + f(         1,j,k) ) &
                     + a4_d15 * ( f(         8 ,j,k) - f(         2,j,k) ) &
                     + a5_d15 * ( f(         9 ,j,k) - f(         3,j,k) ) &
                     + a6_d15 * ( f(         10,j,k) - f(         4,j,k) ) &
                     + a7_d15 * ( f(         11,j,k) - f(         5,j,k) )
                RHS(         5,j,k) = a0_d15 * ( f(         5,j,k) )       &
                     + a1_d15 * ( f(         6 ,j,k) + f(         4,j,k) ) &
                     + a2_d15 * ( f(         7 ,j,k) + f(         3,j,k) ) &
                     + a3_d15 * ( f(         8 ,j,k) + f(         2,j,k) ) &
                     + a4_d15 * ( f(         9 ,j,k) + f(         1,j,k) ) &
                     + a5_d15 * ( f(         10,j,k) - f(         2,j,k) ) &
                     + a6_d15 * ( f(         11,j,k) - f(         3,j,k) ) &
                     + a7_d15 * ( f(         12,j,k) - f(         4,j,k) )
                RHS(         6,j,k) = a0_d15 * ( f(         6,j,k) )       &
                     + a1_d15 * ( f(         7 ,j,k) + f(         5,j,k) ) &
                     + a2_d15 * ( f(         8 ,j,k) + f(         4,j,k) ) &
                     + a3_d15 * ( f(         9 ,j,k) + f(         3,j,k) ) &
                     + a4_d15 * ( f(         10,j,k) + f(         2,j,k) ) &
                     + a5_d15 * ( f(         11,j,k) + f(         1,j,k) ) &
                     + a6_d15 * ( f(         12,j,k) - f(         2,j,k) ) &
                     + a7_d15 * ( f(         13,j,k) - f(         3,j,k) )
                RHS(         7,j,k) = a0_d15 * ( f(         7,j,k) )       &
                     + a1_d15 * ( f(         8 ,j,k) + f(         6,j,k) ) &
                     + a2_d15 * ( f(         9 ,j,k) + f(         5,j,k) ) &
                     + a3_d15 * ( f(         10,j,k) + f(         4,j,k) ) &
                     + a4_d15 * ( f(         11,j,k) + f(         3,j,k) ) &
                     + a5_d15 * ( f(         12,j,k) + f(         2,j,k) ) &
                     + a6_d15 * ( f(         13,j,k) + f(         1,j,k) ) &
                     + a7_d15 * ( f(         14,j,k) - f(         2,j,k) )
             enddo
          enddo

       end select

       ! bcn
       select case(bcn)
       case(0)
          do k=1,n3
             do j=1,n2
                RHS(  this%n-6,j,k) = a0_d15 * ( f(  this%n-6,j,k) )                    &
                     + b7_b1_d15 * f(   this%n-5,j,k) + b7_a1_d15 * f(  this%n-7 ,j,k)  &
                     + b7_b2_d15 * f(   this%n-4,j,k) + b7_a2_d15 * f(  this%n-8 ,j,k)  &
                     + b7_b3_d15 * f(   this%n-3,j,k) + b7_a3_d15 * f(  this%n-9 ,j,k)  &
                     + b7_b4_d15 * f(   this%n-2,j,k) + b7_a4_d15 * f(  this%n-10,j,k)  &
                     + b7_b5_d15 * f(   this%n-1,j,k) + b7_a5_d15 * f(  this%n-11,j,k)  &
                     + b7_b6_d15 * f(   this%n  ,j,k) + b7_a6_d15 * f(  this%n-12,j,k)  &
                                                      + b7_a7_d15 * f(  this%n-13,j,k) 
                RHS(  this%n-5,j,k) = a0_d15 * ( f(  this%n-5,j,k) )                    &
                     + b6_b1_d15 * f(   this%n-4,j,k) + b6_a1_d15 * f(  this%n-6 ,j,k)  &
                     + b6_b2_d15 * f(   this%n-3,j,k) + b6_a2_d15 * f(  this%n-7 ,j,k)  &
                     + b6_b3_d15 * f(   this%n-2,j,k) + b6_a3_d15 * f(  this%n-8 ,j,k)  &
                     + b6_b4_d15 * f(   this%n-1,j,k) + b6_a4_d15 * f(  this%n-9 ,j,k)  &
                     + b6_b5_d15 * f(   this%n  ,j,k) + b6_a5_d15 * f(  this%n-10,j,k)  &
                                                      + b6_a6_d15 * f(  this%n-11,j,k)  &
                                                      + b6_a7_d15 * f(  this%n-12,j,k) 
                RHS(  this%n-4,j,k) = a0_d15 * ( f(  this%n-4,j,k) )                    &
                     + b5_b1_d15 * f(   this%n-3,j,k) + b5_a1_d15 * f(  this%n-5 ,j,k)  &
                     + b5_b2_d15 * f(   this%n-2,j,k) + b5_a2_d15 * f(  this%n-6 ,j,k)  &
                     + b5_b3_d15 * f(   this%n-1,j,k) + b5_a3_d15 * f(  this%n-7 ,j,k)  &
                     + b5_b4_d15 * f(   this%n  ,j,k) + b5_a4_d15 * f(  this%n-8 ,j,k)  &
                                                      + b5_a5_d15 * f(  this%n-9 ,j,k)  &
                                                      + b5_a6_d15 * f(  this%n-10,j,k)  &
                                                      + b5_a7_d15 * f(  this%n-11,j,k) 
                RHS(  this%n-3,j,k) = a0_d15 * ( f(  this%n-3,j,k) )                    &
                     + b4_b1_d15 * f(   this%n-2,j,k) + b4_a1_d15 * f(  this%n-4 ,j,k)  &
                     + b4_b2_d15 * f(   this%n-1,j,k) + b4_a2_d15 * f(  this%n-5 ,j,k)  &
                     + b4_b3_d15 * f(   this%n  ,j,k) + b4_a3_d15 * f(  this%n-6 ,j,k)  &
                                                      + b4_a4_d15 * f(  this%n-7 ,j,k)  &
                                                      + b4_a5_d15 * f(  this%n-8 ,j,k)  &
                                                      + b4_a6_d15 * f(  this%n-9 ,j,k)  &
                                                      + b4_a7_d15 * f(  this%n-10,j,k) 
                RHS(  this%n-2,j,k) = a0_d15 * ( f(  this%n-2,j,k) )                    &
                     + b3_b1_d15 * f(   this%n-1,j,k) + b3_a1_d15 * f(  this%n-3 ,j,k)  &
                     + b3_b2_d15 * f(   this%n  ,j,k) + b3_a2_d15 * f(  this%n-4 ,j,k)  &
                                                      + b3_a3_d15 * f(  this%n-5 ,j,k)  &
                                                      + b3_a4_d15 * f(  this%n-6 ,j,k)  &
                                                      + b3_a5_d15 * f(  this%n-7 ,j,k)  &
                                                      + b3_a6_d15 * f(  this%n-8 ,j,k)  &
                                                      + b3_a7_d15 * f(  this%n-9 ,j,k) 
                RHS(  this%n-1,j,k) = a0_d15 * ( f(  this%n-1,j,k) )                    &
                     + b2_b1_d15 * f(   this%n  ,j,k) + b2_a1_d15 * f(  this%n-2 ,j,k)  &
                                                      + b2_a2_d15 * f(  this%n-3 ,j,k)  &
                                                      + b2_a3_d15 * f(  this%n-4 ,j,k)  &
                                                      + b2_a4_d15 * f(  this%n-5 ,j,k)  &
                                                      + b2_a5_d15 * f(  this%n-6 ,j,k)  &
                                                      + b2_a6_d15 * f(  this%n-7 ,j,k)  &
                                                      + b2_a7_d15 * f(  this%n-8 ,j,k) 
                RHS(  this%n  ,j,k) = a0_d15 * ( f(  this%n  ,j,k) )                    
             end do
          end do

       case(1)
          do k=1,n3
             do j=1,n2
                RHS(  this%n-6,j,k) = a0_d15 * ( f(  this%n-6,j,k) )        &
                     + a1_d15 * ( f(   this%n-5,j,k) + f(  this%n-7 ,j,k) ) &
                     + a2_d15 * ( f(   this%n-4,j,k) + f(  this%n-8 ,j,k) ) &
                     + a3_d15 * ( f(   this%n-3,j,k) + f(  this%n-9 ,j,k) ) &
                     + a4_d15 * ( f(   this%n-2,j,k) + f(  this%n-10,j,k) ) &
                     + a5_d15 * ( f(   this%n-1,j,k) + f(  this%n-11,j,k) ) &
                     + a6_d15 * ( f(   this%n  ,j,k) + f(  this%n-12,j,k) ) &
                     + a7_d15 * ( f(   this%n-1,j,k) + f(  this%n-13,j,k) )
                RHS(  this%n-5,j,k) = a0_d15 * ( f(  this%n-5,j,k) )        &
                     + a1_d15 * ( f(   this%n-4,j,k) + f(  this%n-6 ,j,k) ) &
                     + a2_d15 * ( f(   this%n-3,j,k) + f(  this%n-7 ,j,k) ) &
                     + a3_d15 * ( f(   this%n-2,j,k) + f(  this%n-8 ,j,k) ) &
                     + a4_d15 * ( f(   this%n-1,j,k) + f(  this%n-9 ,j,k) ) &
                     + a5_d15 * ( f(   this%n  ,j,k) + f(  this%n-10,j,k) ) &
                     + a6_d15 * ( f(   this%n-1,j,k) + f(  this%n-11,j,k) ) &
                     + a7_d15 * ( f(   this%n-2,j,k) + f(  this%n-12,j,k) )
                RHS(  this%n-4,j,k) = a0_d15 * ( f(  this%n-4,j,k) )        &
                     + a1_d15 * ( f(   this%n-3,j,k) + f(  this%n-5 ,j,k) ) &
                     + a2_d15 * ( f(   this%n-2,j,k) + f(  this%n-6 ,j,k) ) &
                     + a3_d15 * ( f(   this%n-1,j,k) + f(  this%n-7 ,j,k) ) &
                     + a4_d15 * ( f(   this%n  ,j,k) + f(  this%n-8 ,j,k) ) &
                     + a5_d15 * ( f(   this%n-1,j,k) + f(  this%n-9 ,j,k) ) &
                     + a6_d15 * ( f(   this%n-2,j,k) + f(  this%n-10,j,k) ) &
                     + a7_d15 * ( f(   this%n-3,j,k) + f(  this%n-11,j,k) )
                RHS(  this%n-3,j,k) = a0_d15 * ( f(  this%n-3,j,k) )        &
                     + a1_d15 * ( f(   this%n-2,j,k) + f(  this%n-4 ,j,k) ) &
                     + a2_d15 * ( f(   this%n-1,j,k) + f(  this%n-5 ,j,k) ) &
                     + a3_d15 * ( f(   this%n  ,j,k) + f(  this%n-6 ,j,k) ) &
                     + a4_d15 * ( f(   this%n-1,j,k) + f(  this%n-7 ,j,k) ) &
                     + a5_d15 * ( f(   this%n-2,j,k) + f(  this%n-8 ,j,k) ) &
                     + a6_d15 * ( f(   this%n-3,j,k) + f(  this%n-9 ,j,k) ) &
                     + a7_d15 * ( f(   this%n-4,j,k) + f(  this%n-10,j,k) )
                RHS(  this%n-2,j,k) = a0_d15 * ( f(  this%n-2,j,k) )        &
                     + a1_d15 * ( f(   this%n-1,j,k) + f(  this%n-3 ,j,k) ) &
                     + a2_d15 * ( f(   this%n  ,j,k) + f(  this%n-4 ,j,k) ) &
                     + a3_d15 * ( f(   this%n-1,j,k) + f(  this%n-5 ,j,k) ) &
                     + a4_d15 * ( f(   this%n-2,j,k) + f(  this%n-6 ,j,k) ) &
                     + a5_d15 * ( f(   this%n-3,j,k) + f(  this%n-7 ,j,k) ) &
                     + a6_d15 * ( f(   this%n-4,j,k) + f(  this%n-8 ,j,k) ) &
                     + a7_d15 * ( f(   this%n-5,j,k) + f(  this%n-9 ,j,k) )
                RHS(  this%n-1,j,k) = a0_d15 * ( f(  this%n-1,j,k) )        &
                     + a1_d15 * ( f(   this%n  ,j,k) + f(  this%n-2 ,j,k) ) &
                     + a2_d15 * ( f(   this%n-1,j,k) + f(  this%n-3 ,j,k) ) &
                     + a3_d15 * ( f(   this%n-2,j,k) + f(  this%n-4 ,j,k) ) &
                     + a4_d15 * ( f(   this%n-3,j,k) + f(  this%n-5 ,j,k) ) &
                     + a5_d15 * ( f(   this%n-4,j,k) + f(  this%n-6 ,j,k) ) &
                     + a6_d15 * ( f(   this%n-5,j,k) + f(  this%n-7 ,j,k) ) &
                     + a7_d15 * ( f(   this%n-6,j,k) + f(  this%n-8 ,j,k) )
                RHS(  this%n  ,j,k) = a0_d15 * ( f(  this%n  ,j,k) )        &
                     + a1_d15 * ( f(   this%n-1,j,k) + f(  this%n-1 ,j,k) ) &
                     + a2_d15 * ( f(   this%n-2,j,k) + f(  this%n-2 ,j,k) ) &
                     + a3_d15 * ( f(   this%n-3,j,k) + f(  this%n-3 ,j,k) ) &
                     + a4_d15 * ( f(   this%n-4,j,k) + f(  this%n-4 ,j,k) ) &
                     + a5_d15 * ( f(   this%n-5,j,k) + f(  this%n-5 ,j,k) ) &
                     + a6_d15 * ( f(   this%n-6,j,k) + f(  this%n-6 ,j,k) ) &
                     + a7_d15 * ( f(   this%n-7,j,k) + f(  this%n-7 ,j,k) )
             end do
          end do

       case(-1)
          do k=1,n3
             do j=1,n2
                RHS(  this%n-6,j,k) = a0_d15 * ( f(  this%n-6,j,k) )        &
                     + a1_d15 * ( f(   this%n-5,j,k) + f(  this%n-7 ,j,k) ) &
                     + a2_d15 * ( f(   this%n-4,j,k) + f(  this%n-8 ,j,k) ) &
                     + a3_d15 * ( f(   this%n-3,j,k) + f(  this%n-9 ,j,k) ) &
                     + a4_d15 * ( f(   this%n-2,j,k) + f(  this%n-10,j,k) ) &
                     + a5_d15 * ( f(   this%n-1,j,k) + f(  this%n-11,j,k) ) &
                     + a6_d15 * ( f(   this%n  ,j,k) + f(  this%n-12,j,k) ) &
                     + a7_d15 * (-f(   this%n-1,j,k) + f(  this%n-13,j,k) )
                RHS(  this%n-5,j,k) = a0_d15 * ( f(  this%n-5,j,k) )        &
                     + a1_d15 * ( f(   this%n-4,j,k) + f(  this%n-6 ,j,k) ) &
                     + a2_d15 * ( f(   this%n-3,j,k) + f(  this%n-7 ,j,k) ) &
                     + a3_d15 * ( f(   this%n-2,j,k) + f(  this%n-8 ,j,k) ) &
                     + a4_d15 * ( f(   this%n-1,j,k) + f(  this%n-9 ,j,k) ) &
                     + a5_d15 * ( f(   this%n  ,j,k) + f(  this%n-10,j,k) ) &
                     + a6_d15 * (-f(   this%n-1,j,k) + f(  this%n-11,j,k) ) &
                     + a7_d15 * (-f(   this%n-2,j,k) + f(  this%n-12,j,k) )
                RHS(  this%n-4,j,k) = a0_d15 * ( f(  this%n-4,j,k) )        &
                     + a1_d15 * ( f(   this%n-3,j,k) + f(  this%n-5 ,j,k) ) &
                     + a2_d15 * ( f(   this%n-2,j,k) + f(  this%n-6 ,j,k) ) &
                     + a3_d15 * ( f(   this%n-1,j,k) + f(  this%n-7 ,j,k) ) &
                     + a4_d15 * ( f(   this%n  ,j,k) + f(  this%n-8 ,j,k) ) &
                     + a5_d15 * (-f(   this%n-1,j,k) + f(  this%n-9 ,j,k) ) &
                     + a6_d15 * (-f(   this%n-2,j,k) + f(  this%n-10,j,k) ) &
                     + a7_d15 * (-f(   this%n-3,j,k) + f(  this%n-11,j,k) )
                RHS(  this%n-3,j,k) = a0_d15 * ( f(  this%n-3,j,k) )        &
                     + a1_d15 * ( f(   this%n-2,j,k) + f(  this%n-4 ,j,k) ) &
                     + a2_d15 * ( f(   this%n-1,j,k) + f(  this%n-5 ,j,k) ) &
                     + a3_d15 * ( f(   this%n  ,j,k) + f(  this%n-6 ,j,k) ) &
                     + a4_d15 * (-f(   this%n-1,j,k) + f(  this%n-7 ,j,k) ) &
                     + a5_d15 * (-f(   this%n-2,j,k) + f(  this%n-8 ,j,k) ) &
                     + a6_d15 * (-f(   this%n-3,j,k) + f(  this%n-9 ,j,k) ) &
                     + a7_d15 * (-f(   this%n-4,j,k) + f(  this%n-10,j,k) )
                RHS(  this%n-2,j,k) = a0_d15 * ( f(  this%n-2,j,k) )        &
                     + a1_d15 * ( f(   this%n-1,j,k) + f(  this%n-3 ,j,k) ) &
                     + a2_d15 * ( f(   this%n  ,j,k) + f(  this%n-4 ,j,k) ) &
                     + a3_d15 * (-f(   this%n-1,j,k) + f(  this%n-5 ,j,k) ) &
                     + a4_d15 * (-f(   this%n-2,j,k) + f(  this%n-6 ,j,k) ) &
                     + a5_d15 * (-f(   this%n-3,j,k) + f(  this%n-7 ,j,k) ) &
                     + a6_d15 * (-f(   this%n-4,j,k) + f(  this%n-8 ,j,k) ) &
                     + a7_d15 * (-f(   this%n-5,j,k) + f(  this%n-9 ,j,k) )
                RHS(  this%n-1,j,k) = a0_d15 * ( f(  this%n-1,j,k) )        &
                     + a1_d15 * ( f(   this%n  ,j,k) + f(  this%n-2 ,j,k) ) &
                     + a2_d15 * (-f(   this%n-1,j,k) + f(  this%n-3 ,j,k) ) &
                     + a3_d15 * (-f(   this%n-2,j,k) + f(  this%n-4 ,j,k) ) &
                     + a4_d15 * (-f(   this%n-3,j,k) + f(  this%n-5 ,j,k) ) &
                     + a5_d15 * (-f(   this%n-4,j,k) + f(  this%n-6 ,j,k) ) &
                     + a6_d15 * (-f(   this%n-5,j,k) + f(  this%n-7 ,j,k) ) &
                     + a7_d15 * (-f(   this%n-6,j,k) + f(  this%n-8 ,j,k) )
                RHS(  this%n  ,j,k) = a0_d15 * ( f(  this%n  ,j,k) )        &
                     + a1_d15 * (-f(   this%n-1,j,k) + f(  this%n-1 ,j,k) ) &
                     + a2_d15 * (-f(   this%n-2,j,k) + f(  this%n-2 ,j,k) ) &
                     + a3_d15 * (-f(   this%n-3,j,k) + f(  this%n-3 ,j,k) ) &
                     + a4_d15 * (-f(   this%n-4,j,k) + f(  this%n-4 ,j,k) ) &
                     + a5_d15 * (-f(   this%n-5,j,k) + f(  this%n-5 ,j,k) ) &
                     + a6_d15 * (-f(   this%n-6,j,k) + f(  this%n-6 ,j,k) ) &
                     + a7_d15 * (-f(   this%n-7,j,k) + f(  this%n-7 ,j,k) )
             end do
          end do

       end select
    
    end select

 end subroutine ComputeXRHS

 pure subroutine ComputeYRHS(this, f, RHS, n1, n3, bc1, bcn) 

   class( cfo2D15 ), intent(in) :: this
   integer, intent(in) :: n1, n3
   real(rkind), dimension(n1,this%n,n3), intent(in) :: f
   real(rkind), dimension(n1,this%n,n3), intent(out) :: RHS
   integer, intent(in) :: bc1, bcn
   integer :: k

   select case (this%periodic)
   case (.TRUE.)
      do k=1,n3
         RHS(:,         1,k) = a0_d15 * ( f(:,         1,k) )       &
              + a1_d15 * ( f(:,         2 ,k) + f(:,    this%n,k) ) &
              + a2_d15 * ( f(:,         3 ,k) + f(:,  this%n-1,k) ) &
              + a3_d15 * ( f(:,         4 ,k) + f(:,  this%n-2,k) ) &
              + a4_d15 * ( f(:,         5 ,k) + f(:,  this%n-3,k) ) &
              + a5_d15 * ( f(:,         6 ,k) + f(:,  this%n-4,k) ) &
              + a6_d15 * ( f(:,         7 ,k) + f(:,  this%n-5,k) ) &
              + a7_d15 * ( f(:,         8 ,k) + f(:,  this%n-6,k) )
         RHS(:,         2,k) = a0_d15 * ( f(:,         2,k) )       &
              + a1_d15 * ( f(:,         3 ,k) + f(:,         1,k) ) &
              + a2_d15 * ( f(:,         4 ,k) + f(:,    this%n,k) ) &
              + a3_d15 * ( f(:,         5 ,k) + f(:,  this%n-1,k) ) &
              + a4_d15 * ( f(:,         6 ,k) + f(:,  this%n-2,k) ) &
              + a5_d15 * ( f(:,         7 ,k) + f(:,  this%n-3,k) ) &
              + a6_d15 * ( f(:,         8 ,k) + f(:,  this%n-4,k) ) &
              + a7_d15 * ( f(:,         9 ,k) + f(:,  this%n-5,k) )
         RHS(:,         3,k) = a0_d15 * ( f(:,         3,k) )       &
              + a1_d15 * ( f(:,         4 ,k) + f(:,         2,k) ) &
              + a2_d15 * ( f(:,         5 ,k) + f(:,         1,k) ) &
              + a3_d15 * ( f(:,         6 ,k) + f(:,    this%n,k) ) &
              + a4_d15 * ( f(:,         7 ,k) + f(:,  this%n-1,k) ) &
              + a5_d15 * ( f(:,         8 ,k) + f(:,  this%n-2,k) ) &
              + a6_d15 * ( f(:,         9 ,k) + f(:,  this%n-3,k) ) &
              + a7_d15 * ( f(:,         10,k) + f(:,  this%n-4,k) )
         RHS(:,         4,k) = a0_d15 * ( f(:,         4,k) )       &
              + a1_d15 * ( f(:,         5 ,k) + f(:,         3,k) ) &
              + a2_d15 * ( f(:,         6 ,k) + f(:,         2,k) ) &
              + a3_d15 * ( f(:,         7 ,k) + f(:,         1,k) ) &
              + a4_d15 * ( f(:,         8 ,k) + f(:,    this%n,k) ) &
              + a5_d15 * ( f(:,         9 ,k) + f(:,  this%n-1,k) ) &
              + a6_d15 * ( f(:,         10,k) + f(:,  this%n-2,k) ) &
              + a7_d15 * ( f(:,         11,k) + f(:,  this%n-3,k) )
         RHS(:,         5,k) = a0_d15 * ( f(:,         5,k) )       &
              + a1_d15 * ( f(:,         6 ,k) + f(:,         4,k) ) &
              + a2_d15 * ( f(:,         7 ,k) + f(:,         3,k) ) &
              + a3_d15 * ( f(:,         8 ,k) + f(:,         2,k) ) &
              + a4_d15 * ( f(:,         9 ,k) + f(:,         1,k) ) &
              + a5_d15 * ( f(:,         10,k) + f(:,    this%n,k) ) &
              + a6_d15 * ( f(:,         11,k) + f(:,  this%n-1,k) ) &
              + a7_d15 * ( f(:,         12,k) + f(:,  this%n-2,k) )
         RHS(:,         6,k) = a0_d15 * ( f(:,         6,k) )       &
              + a1_d15 * ( f(:,         7 ,k) + f(:,         5,k) ) &
              + a2_d15 * ( f(:,         8 ,k) + f(:,         4,k) ) &
              + a3_d15 * ( f(:,         9 ,k) + f(:,         3,k) ) &
              + a4_d15 * ( f(:,         10,k) + f(:,         2,k) ) &
              + a5_d15 * ( f(:,         11,k) + f(:,         1,k) ) &
              + a6_d15 * ( f(:,         12,k) + f(:,    this%n,k) ) &
              + a7_d15 * ( f(:,         13,k) + f(:,  this%n-1,k) )
         RHS(:,         7,k) = a0_d15 * ( f(:,         7,k) )       &
              + a1_d15 * ( f(:,         8 ,k) + f(:,         6,k) ) &
              + a2_d15 * ( f(:,         9 ,k) + f(:,         5,k) ) &
              + a3_d15 * ( f(:,         10,k) + f(:,         4,k) ) &
              + a4_d15 * ( f(:,         11,k) + f(:,         3,k) ) &
              + a5_d15 * ( f(:,         12,k) + f(:,         2,k) ) &
              + a6_d15 * ( f(:,         13,k) + f(:,         1,k) ) &
              + a7_d15 * ( f(:,         14,k) + f(:,    this%n,k) )
         RHS(:,8:this%n-7,k) = a0_d15 * ( f(:,8:this%n-7,k) )        &
              + a1_d15 * ( f(:,9 :this%n-6,k) + f(:,7:this%n-8 ,k) ) &
              + a2_d15 * ( f(:,10:this%n-5,k) + f(:,6:this%n-9 ,k) ) &
              + a3_d15 * ( f(:,11:this%n-4,k) + f(:,5:this%n-10,k) ) &
              + a4_d15 * ( f(:,12:this%n-3,k) + f(:,4:this%n-11,k) ) &
              + a5_d15 * ( f(:,13:this%n-2,k) + f(:,3:this%n-12,k) ) &
              + a6_d15 * ( f(:,14:this%n-1,k) + f(:,2:this%n-13,k) ) &
              + a7_d15 * ( f(:,15:this%n  ,k) + f(:,1:this%n-14,k) )
         RHS(:,  this%n-6,k) = a0_d15 * ( f(:,  this%n-6,k) )       &
              + a1_d15 * ( f(:,  this%n-5,k) + f(:,  this%n-7 ,k) ) &
              + a2_d15 * ( f(:,  this%n-4,k) + f(:,  this%n-8 ,k) ) &
              + a3_d15 * ( f(:,  this%n-3,k) + f(:,  this%n-9 ,k) ) &
              + a4_d15 * ( f(:,  this%n-2,k) + f(:,  this%n-10,k) ) &
              + a5_d15 * ( f(:,  this%n-1,k) + f(:,  this%n-11,k) ) &
              + a6_d15 * ( f(:,  this%n  ,k) + f(:,  this%n-12,k) ) &
              + a7_d15 * ( f(:,         1,k) + f(:,  this%n-13,k) )
         RHS(:,  this%n-5,k) = a0_d15 * ( f(:,  this%n-5,k) )       &
              + a1_d15 * ( f(:,  this%n-4,k) + f(:,  this%n-6 ,k) ) &
              + a2_d15 * ( f(:,  this%n-3,k) + f(:,  this%n-7 ,k) ) &
              + a3_d15 * ( f(:,  this%n-2,k) + f(:,  this%n-8 ,k) ) &
              + a4_d15 * ( f(:,  this%n-1,k) + f(:,  this%n-9 ,k) ) &
              + a5_d15 * ( f(:,    this%n,k) + f(:,  this%n-10,k) ) &
              + a6_d15 * ( f(:,         1,k) + f(:,  this%n-11,k) ) &
              + a7_d15 * ( f(:,         2,k) + f(:,  this%n-12,k) )
         RHS(:,  this%n-4,k) = a0_d15 * ( f(:,  this%n-4,k) )       &
              + a1_d15 * ( f(:,  this%n-3,k) + f(:,  this%n-5 ,k) ) &
              + a2_d15 * ( f(:,  this%n-2,k) + f(:,  this%n-6 ,k) ) &
              + a3_d15 * ( f(:,  this%n-1,k) + f(:,  this%n-7 ,k) ) &
              + a4_d15 * ( f(:,    this%n,k) + f(:,  this%n-8 ,k) ) &
              + a5_d15 * ( f(:,         1,k) + f(:,  this%n-9 ,k) ) &
              + a6_d15 * ( f(:,         2,k) + f(:,  this%n-10,k) ) &
              + a7_d15 * ( f(:,         3,k) + f(:,  this%n-11,k) )
         RHS(:,  this%n-3,k) = a0_d15 * ( f(:,  this%n-3,k) )       &
              + a1_d15 * ( f(:,  this%n-2,k) + f(:,  this%n-4 ,k) ) &
              + a2_d15 * ( f(:,  this%n-1,k) + f(:,  this%n-5 ,k) ) &
              + a3_d15 * ( f(:,    this%n,k) + f(:,  this%n-6 ,k) ) &
              + a4_d15 * ( f(:,         1,k) + f(:,  this%n-7 ,k) ) &
              + a5_d15 * ( f(:,         2,k) + f(:,  this%n-8 ,k) ) &
              + a6_d15 * ( f(:,         3,k) + f(:,  this%n-9 ,k) ) &
              + a7_d15 * ( f(:,         4,k) + f(:,  this%n-10,k) )
         RHS(:,  this%n-2,k) = a0_d15 * ( f(:,  this%n-2,k) )       &
              + a1_d15 * ( f(:,  this%n-1,k) + f(:,  this%n-3 ,k) ) &
              + a2_d15 * ( f(:,    this%n,k) + f(:,  this%n-4 ,k) ) &
              + a3_d15 * ( f(:,         1,k) + f(:,  this%n-5 ,k) ) &
              + a4_d15 * ( f(:,         2,k) + f(:,  this%n-6 ,k) ) &
              + a5_d15 * ( f(:,         3,k) + f(:,  this%n-7 ,k) ) &
              + a6_d15 * ( f(:,         4,k) + f(:,  this%n-8 ,k) ) &
              + a7_d15 * ( f(:,         5,k) + f(:,  this%n-9 ,k) )
         RHS(:,  this%n-1,k) = a0_d15 * ( f(:,  this%n-1,k) )       &
              + a1_d15 * ( f(:,    this%n,k) + f(:,  this%n-2 ,k) ) &
              + a2_d15 * ( f(:,         1,k) + f(:,  this%n-3 ,k) ) &
              + a3_d15 * ( f(:,         2,k) + f(:,  this%n-4 ,k) ) &
              + a4_d15 * ( f(:,         3,k) + f(:,  this%n-5 ,k) ) &
              + a5_d15 * ( f(:,         4,k) + f(:,  this%n-6 ,k) ) &
              + a6_d15 * ( f(:,         5,k) + f(:,  this%n-7 ,k) ) &
              + a7_d15 * ( f(:,         6,k) + f(:,  this%n-8 ,k) )
         RHS(:,  this%n,k) = a0_d15 * ( f(:,  this%n,k) )           &
              + a1_d15 * ( f(:,         1,k) + f(:,  this%n-1 ,k) ) &
              + a2_d15 * ( f(:,         2,k) + f(:,  this%n-2 ,k) ) &
              + a3_d15 * ( f(:,         3,k) + f(:,  this%n-3 ,k) ) &
              + a4_d15 * ( f(:,         4,k) + f(:,  this%n-4 ,k) ) &
              + a5_d15 * ( f(:,         5,k) + f(:,  this%n-5 ,k) ) &
              + a6_d15 * ( f(:,         6,k) + f(:,  this%n-6 ,k) ) &
              + a7_d15 * ( f(:,         7,k) + f(:,  this%n-7 ,k) )
      end do
   case (.FALSE.)

      ! interior points
      do k=1,n3
         RHS(:,8:this%n-7,k) = a0_d15 * ( f(:,8:this%n-7,k) )        &
              + a1_d15 * ( f(:,9 :this%n-6,k) + f(:,7:this%n-8 ,k) ) &
              + a2_d15 * ( f(:,10:this%n-5,k) + f(:,6:this%n-9 ,k) ) &
              + a3_d15 * ( f(:,11:this%n-4,k) + f(:,5:this%n-10,k) ) &
              + a4_d15 * ( f(:,12:this%n-3,k) + f(:,4:this%n-11,k) ) &
              + a5_d15 * ( f(:,13:this%n-2,k) + f(:,3:this%n-12,k) ) &
              + a6_d15 * ( f(:,14:this%n-1,k) + f(:,2:this%n-13,k) ) &
              + a7_d15 * ( f(:,15:this%n  ,k) + f(:,1:this%n-14,k) )
      enddo

      ! bc 1
      select case(bc1)
      case(0)

         do k=1,n3
            RHS(:,         1,k) = b1_a0_d15 * ( f(:,         1,k) )                
            RHS(:,         2,k) = b2_a0_d15 * ( f(:,         2,k) )                &
                 + b2_a1_d15 * f(:,         3 ,k) + b2_b1_d15 * f(:,         1,k)  &
                 + b2_a2_d15 * f(:,         4 ,k)                                  &
                 + b2_a3_d15 * f(:,         5 ,k)                                  &
                 + b2_a4_d15 * f(:,         6 ,k)                                  &
                 + b2_a5_d15 * f(:,         7 ,k)                                  &
                 + b2_a6_d15 * f(:,         8 ,k)                                  &
                 + b2_a7_d15 * f(:,         9 ,k)                     
            RHS(:,         3,k) = b3_a0_d15 * ( f(:,         3,k) )                &
                 + b3_a1_d15 * f(:,         4 ,k) + b3_b1_d15 * f(:,         2,k)  &
                 + b3_a2_d15 * f(:,         5 ,k) + b3_b2_d15 * f(:,         1,k)  &
                 + b3_a3_d15 * f(:,         6 ,k)                                  &
                 + b3_a4_d15 * f(:,         7 ,k)                                  &
                 + b3_a5_d15 * f(:,         8 ,k)                                  &
                 + b3_a6_d15 * f(:,         9 ,k)                                  &
                 + b3_a7_d15 * f(:,         10,k)                     
            RHS(:,         4,k) = b4_a0_d15 * ( f(:,         4,k) )                &
                 + b4_a1_d15 * f(:,         5 ,k) + b4_b1_d15 * f(:,         3,k)  &
                 + b4_a2_d15 * f(:,         6 ,k) + b4_b2_d15 * f(:,         2,k)  &
                 + b4_a3_d15 * f(:,         7 ,k) + b4_b3_d15 * f(:,         1,k)  &
                 + b4_a4_d15 * f(:,         8 ,k)                                  &
                 + b4_a5_d15 * f(:,         9 ,k)                                  &
                 + b4_a6_d15 * f(:,         10,k)                                  &
                 + b4_a7_d15 * f(:,         11,k)                     
            RHS(:,         5,k) = b5_a0_d15 * ( f(:,         5,k) )                &
                 + b5_a1_d15 * f(:,         6 ,k) + b5_b1_d15 * f(:,         4,k)  &
                 + b5_a2_d15 * f(:,         7 ,k) + b5_b2_d15 * f(:,         3,k)  &
                 + b5_a3_d15 * f(:,         8 ,k) + b5_b3_d15 * f(:,         2,k)  &
                 + b5_a4_d15 * f(:,         9 ,k) + b5_b4_d15 * f(:,         1,k)  &
                 + b5_a5_d15 * f(:,         10,k)                                  &
                 + b5_a6_d15 * f(:,         11,k)                                  &
                 + b5_a7_d15 * f(:,         12,k)                     
            RHS(:,         6,k) = b6_a0_d15 * ( f(:,         6,k) )                &
                 + b6_a1_d15 * f(:,         7 ,k) + b6_b1_d15 * f(:,         5,k)  &
                 + b6_a2_d15 * f(:,         8 ,k) + b6_b2_d15 * f(:,         4,k)  &
                 + b6_a3_d15 * f(:,         9 ,k) + b6_b3_d15 * f(:,         3,k)  &
                 + b6_a4_d15 * f(:,         10,k) + b6_b4_d15 * f(:,         2,k)  &
                 + b6_a5_d15 * f(:,         11,k) + b6_b5_d15 * f(:,         1,k)  &
                 + b6_a6_d15 * f(:,         12,k)                                  &
                 + b6_a7_d15 * f(:,         13,k)                     
            RHS(:,         7,k) = b7_a0_d15 * ( f(:,         7,k) )                &
                 + b7_a1_d15 * f(:,         8 ,k) + b7_b1_d15 * f(:,         6,k)  &
                 + b7_a2_d15 * f(:,         9 ,k) + b7_b2_d15 * f(:,         5,k)  &
                 + b7_a3_d15 * f(:,         10,k) + b7_b3_d15 * f(:,         4,k)  &
                 + b7_a4_d15 * f(:,         11,k) + b7_b4_d15 * f(:,         3,k)  &
                 + b7_a5_d15 * f(:,         12,k) + b7_b5_d15 * f(:,         2,k)  &
                 + b7_a6_d15 * f(:,         13,k) + b7_b6_d15 * f(:,         1,k)  &
                 + b7_a7_d15 * f(:,         14,k)                     
         enddo

      case(1)
         do k=1,n3
            RHS(:,         1,k) = a0_d15 * ( f(:,         1,k) )       &
                 + a1_d15 * ( f(:,         2 ,k) + f(:,         2,k) ) &
                 + a2_d15 * ( f(:,         3 ,k) + f(:,         3,k) ) &
                 + a3_d15 * ( f(:,         4 ,k) + f(:,         4,k) ) &
                 + a4_d15 * ( f(:,         5 ,k) + f(:,         5,k) ) &
                 + a5_d15 * ( f(:,         6 ,k) + f(:,         6,k) ) &
                 + a6_d15 * ( f(:,         7 ,k) + f(:,         7,k) ) &
                 + a7_d15 * ( f(:,         8 ,k) + f(:,         8,k) )
            RHS(:,         2,k) = a0_d15 * ( f(:,         2,k) )       &
                 + a1_d15 * ( f(:,         3 ,k) + f(:,         1,k) ) &
                 + a2_d15 * ( f(:,         4 ,k) + f(:,         2,k) ) &
                 + a3_d15 * ( f(:,         5 ,k) + f(:,         3,k) ) &
                 + a4_d15 * ( f(:,         6 ,k) + f(:,         4,k) ) &
                 + a5_d15 * ( f(:,         7 ,k) + f(:,         5,k) ) &
                 + a6_d15 * ( f(:,         8 ,k) + f(:,         6,k) ) &
                 + a7_d15 * ( f(:,         9 ,k) + f(:,         7,k) )
            RHS(:,         3,k) = a0_d15 * ( f(:,         3,k) )       &
                 + a1_d15 * ( f(:,         4 ,k) + f(:,         2,k) ) &
                 + a2_d15 * ( f(:,         5 ,k) + f(:,         1,k) ) &
                 + a3_d15 * ( f(:,         6 ,k) + f(:,         2,k) ) &
                 + a4_d15 * ( f(:,         7 ,k) + f(:,         3,k) ) &
                 + a5_d15 * ( f(:,         8 ,k) + f(:,         4,k) ) &
                 + a6_d15 * ( f(:,         9 ,k) + f(:,         5,k) ) &
                 + a7_d15 * ( f(:,         10,k) + f(:,         6,k) )
            RHS(:,         4,k) = a0_d15 * ( f(:,         4,k) )       &
                 + a1_d15 * ( f(:,         5 ,k) + f(:,         3,k) ) &
                 + a2_d15 * ( f(:,         6 ,k) + f(:,         2,k) ) &
                 + a3_d15 * ( f(:,         7 ,k) + f(:,         1,k) ) &
                 + a4_d15 * ( f(:,         8 ,k) + f(:,         2,k) ) &
                 + a5_d15 * ( f(:,         9 ,k) + f(:,         3,k) ) &
                 + a6_d15 * ( f(:,         10,k) + f(:,         4,k) ) &
                 + a7_d15 * ( f(:,         11,k) + f(:,         5,k) )
            RHS(:,         5,k) = a0_d15 * ( f(:,         5,k) )       &
                 + a1_d15 * ( f(:,         6 ,k) + f(:,         4,k) ) &
                 + a2_d15 * ( f(:,         7 ,k) + f(:,         3,k) ) &
                 + a3_d15 * ( f(:,         8 ,k) + f(:,         2,k) ) &
                 + a4_d15 * ( f(:,         9 ,k) + f(:,         1,k) ) &
                 + a5_d15 * ( f(:,         10,k) + f(:,         2,k) ) &
                 + a6_d15 * ( f(:,         11,k) + f(:,         3,k) ) &
                 + a7_d15 * ( f(:,         12,k) + f(:,         4,k) )
            RHS(:,         6,k) = a0_d15 * ( f(:,         6,k) )       &
                 + a1_d15 * ( f(:,         7 ,k) + f(:,         5,k) ) &
                 + a2_d15 * ( f(:,         8 ,k) + f(:,         4,k) ) &
                 + a3_d15 * ( f(:,         9 ,k) + f(:,         3,k) ) &
                 + a4_d15 * ( f(:,         10,k) + f(:,         2,k) ) &
                 + a5_d15 * ( f(:,         11,k) + f(:,         1,k) ) &
                 + a6_d15 * ( f(:,         12,k) + f(:,         2,k) ) &
                 + a7_d15 * ( f(:,         13,k) + f(:,         3,k) )
            RHS(:,         7,k) = a0_d15 * ( f(:,         7,k) )       &
                 + a1_d15 * ( f(:,         8 ,k) + f(:,         6,k) ) &
                 + a2_d15 * ( f(:,         9 ,k) + f(:,         5,k) ) &
                 + a3_d15 * ( f(:,         10,k) + f(:,         4,k) ) &
                 + a4_d15 * ( f(:,         11,k) + f(:,         3,k) ) &
                 + a5_d15 * ( f(:,         12,k) + f(:,         2,k) ) &
                 + a6_d15 * ( f(:,         13,k) + f(:,         1,k) ) &
                 + a7_d15 * ( f(:,         14,k) + f(:,         2,k) )
         enddo

      case(-1)
         do k=1,n3
            RHS(:,         1,k) = a0_d15 * ( f(:,         1,k) )       &
                 + a1_d15 * ( f(:,         2 ,k) - f(:,         2,k) ) &
                 + a2_d15 * ( f(:,         3 ,k) - f(:,         3,k) ) &
                 + a3_d15 * ( f(:,         4 ,k) - f(:,         4,k) ) &
                 + a4_d15 * ( f(:,         5 ,k) - f(:,         5,k) ) &
                 + a5_d15 * ( f(:,         6 ,k) - f(:,         6,k) ) &
                 + a6_d15 * ( f(:,         7 ,k) - f(:,         7,k) ) &
                 + a7_d15 * ( f(:,         8 ,k) - f(:,         8,k) )
            RHS(:,         2,k) = a0_d15 * ( f(:,         2,k) )       &
                 + a1_d15 * ( f(:,         3 ,k) + f(:,         1,k) ) &
                 + a2_d15 * ( f(:,         4 ,k) - f(:,         2,k) ) &
                 + a3_d15 * ( f(:,         5 ,k) - f(:,         3,k) ) &
                 + a4_d15 * ( f(:,         6 ,k) - f(:,         4,k) ) &
                 + a5_d15 * ( f(:,         7 ,k) - f(:,         5,k) ) &
                 + a6_d15 * ( f(:,         8 ,k) - f(:,         6,k) ) &
                 + a7_d15 * ( f(:,         9 ,k) - f(:,         7,k) )
            RHS(:,         3,k) = a0_d15 * ( f(:,         3,k) )       &
                 + a1_d15 * ( f(:,         4 ,k) + f(:,         2,k) ) &
                 + a2_d15 * ( f(:,         5 ,k) + f(:,         1,k) ) &
                 + a3_d15 * ( f(:,         6 ,k) - f(:,         2,k) ) &
                 + a4_d15 * ( f(:,         7 ,k) - f(:,         3,k) ) &
                 + a5_d15 * ( f(:,         8 ,k) - f(:,         4,k) ) &
                 + a6_d15 * ( f(:,         9 ,k) - f(:,         5,k) ) &
                 + a7_d15 * ( f(:,         10,k) - f(:,         6,k) )
            RHS(:,         4,k) = a0_d15 * ( f(:,         4,k) )       &
                 + a1_d15 * ( f(:,         5 ,k) + f(:,         3,k) ) &
                 + a2_d15 * ( f(:,         6 ,k) + f(:,         2,k) ) &
                 + a3_d15 * ( f(:,         7 ,k) + f(:,         1,k) ) &
                 + a4_d15 * ( f(:,         8 ,k) - f(:,         2,k) ) &
                 + a5_d15 * ( f(:,         9 ,k) - f(:,         3,k) ) &
                 + a6_d15 * ( f(:,         10,k) - f(:,         4,k) ) &
                 + a7_d15 * ( f(:,         11,k) - f(:,         5,k) )
            RHS(:,         5,k) = a0_d15 * ( f(:,         5,k) )       &
                 + a1_d15 * ( f(:,         6 ,k) + f(:,         4,k) ) &
                 + a2_d15 * ( f(:,         7 ,k) + f(:,         3,k) ) &
                 + a3_d15 * ( f(:,         8 ,k) + f(:,         2,k) ) &
                 + a4_d15 * ( f(:,         9 ,k) + f(:,         1,k) ) &
                 + a5_d15 * ( f(:,         10,k) - f(:,         2,k) ) &
                 + a6_d15 * ( f(:,         11,k) - f(:,         3,k) ) &
                 + a7_d15 * ( f(:,         12,k) - f(:,         4,k) )
            RHS(:,         6,k) = a0_d15 * ( f(:,         6,k) )       &
                 + a1_d15 * ( f(:,         7 ,k) + f(:,         5,k) ) &
                 + a2_d15 * ( f(:,         8 ,k) + f(:,         4,k) ) &
                 + a3_d15 * ( f(:,         9 ,k) + f(:,         3,k) ) &
                 + a4_d15 * ( f(:,         10,k) + f(:,         2,k) ) &
                 + a5_d15 * ( f(:,         11,k) + f(:,         1,k) ) &
                 + a6_d15 * ( f(:,         12,k) - f(:,         2,k) ) &
                 + a7_d15 * ( f(:,         13,k) - f(:,         3,k) )
            RHS(:,         7,k) = a0_d15 * ( f(:,         7,k) )       &
                 + a1_d15 * ( f(:,         8 ,k) + f(:,         6,k) ) &
                 + a2_d15 * ( f(:,         9 ,k) + f(:,         5,k) ) &
                 + a3_d15 * ( f(:,         10,k) + f(:,         4,k) ) &
                 + a4_d15 * ( f(:,         11,k) + f(:,         3,k) ) &
                 + a5_d15 * ( f(:,         12,k) + f(:,         2,k) ) &
                 + a6_d15 * ( f(:,         13,k) + f(:,         1,k) ) &
                 + a7_d15 * ( f(:,         14,k) - f(:,         2,k) )
         enddo

      end select

      ! bcn
      select case(bcn)
      case(0)
         do k=1,n3
            RHS(:,  this%n-6,k) = a0_d15 * ( f(:,  this%n-6,k) )                    &
                 + b7_b1_d15 * f(:,   this%n-5,k) + b7_a1_d15 * f(:,  this%n-7 ,k)  &
                 + b7_b2_d15 * f(:,   this%n-4,k) + b7_a2_d15 * f(:,  this%n-8 ,k)  &
                 + b7_b3_d15 * f(:,   this%n-3,k) + b7_a3_d15 * f(:,  this%n-9 ,k)  &
                 + b7_b4_d15 * f(:,   this%n-2,k) + b7_a4_d15 * f(:,  this%n-10,k)  &
                 + b7_b5_d15 * f(:,   this%n-1,k) + b7_a5_d15 * f(:,  this%n-11,k)  &
                 + b7_b6_d15 * f(:,   this%n  ,k) + b7_a6_d15 * f(:,  this%n-12,k)  &
                                                  + b7_a7_d15 * f(:,  this%n-13,k) 
            RHS(:,  this%n-5,k) = a0_d15 * ( f(:,  this%n-5,k) )                    &
                 + b6_b1_d15 * f(:,   this%n-4,k) + b6_a1_d15 * f(:,  this%n-6 ,k)  &
                 + b6_b2_d15 * f(:,   this%n-3,k) + b6_a2_d15 * f(:,  this%n-7 ,k)  &
                 + b6_b3_d15 * f(:,   this%n-2,k) + b6_a3_d15 * f(:,  this%n-8 ,k)  &
                 + b6_b4_d15 * f(:,   this%n-1,k) + b6_a4_d15 * f(:,  this%n-9 ,k)  &
                 + b6_b5_d15 * f(:,   this%n  ,k) + b6_a5_d15 * f(:,  this%n-10,k)  &
                                                  + b6_a6_d15 * f(:,  this%n-11,k)  &
                                                  + b6_a7_d15 * f(:,  this%n-12,k) 
            RHS(:,  this%n-4,k) = a0_d15 * ( f(:,  this%n-4,k) )                    &
                 + b5_b1_d15 * f(:,   this%n-3,k) + b5_a1_d15 * f(:,  this%n-5 ,k)  &
                 + b5_b2_d15 * f(:,   this%n-2,k) + b5_a2_d15 * f(:,  this%n-6 ,k)  &
                 + b5_b3_d15 * f(:,   this%n-1,k) + b5_a3_d15 * f(:,  this%n-7 ,k)  &
                 + b5_b4_d15 * f(:,   this%n  ,k) + b5_a4_d15 * f(:,  this%n-8 ,k)  &
                                                  + b5_a5_d15 * f(:,  this%n-9 ,k)  &
                                                  + b5_a6_d15 * f(:,  this%n-10,k)  &
                                                  + b5_a7_d15 * f(:,  this%n-11,k) 
            RHS(:,  this%n-3,k) = a0_d15 * ( f(:,  this%n-3,k) )                    &
                 + b4_b1_d15 * f(:,   this%n-2,k) + b4_a1_d15 * f(:,  this%n-4 ,k)  &
                 + b4_b2_d15 * f(:,   this%n-1,k) + b4_a2_d15 * f(:,  this%n-5 ,k)  &
                 + b4_b3_d15 * f(:,   this%n  ,k) + b4_a3_d15 * f(:,  this%n-6 ,k)  &
                                                  + b4_a4_d15 * f(:,  this%n-7 ,k)  &
                                                  + b4_a5_d15 * f(:,  this%n-8 ,k)  &
                                                  + b4_a6_d15 * f(:,  this%n-9 ,k)  &
                                                  + b4_a7_d15 * f(:,  this%n-10,k) 
            RHS(:,  this%n-2,k) = a0_d15 * ( f(:,  this%n-2,k) )                    &
                 + b3_b1_d15 * f(:,   this%n-1,k) + b3_a1_d15 * f(:,  this%n-3 ,k)  &
                 + b3_b2_d15 * f(:,   this%n  ,k) + b3_a2_d15 * f(:,  this%n-4 ,k)  &
                                                  + b3_a3_d15 * f(:,  this%n-5 ,k)  &
                                                  + b3_a4_d15 * f(:,  this%n-6 ,k)  &
                                                  + b3_a5_d15 * f(:,  this%n-7 ,k)  &
                                                  + b3_a6_d15 * f(:,  this%n-8 ,k)  &
                                                  + b3_a7_d15 * f(:,  this%n-9 ,k) 
            RHS(:,  this%n-1,k) = a0_d15 * ( f(:,  this%n-1,k) )                    &
                 + b2_b1_d15 * f(:,   this%n  ,k) + b2_a1_d15 * f(:,  this%n-2 ,k)  &
                                                  + b2_a2_d15 * f(:,  this%n-3 ,k)  &
                                                  + b2_a3_d15 * f(:,  this%n-4 ,k)  &
                                                  + b2_a4_d15 * f(:,  this%n-5 ,k)  &
                                                  + b2_a5_d15 * f(:,  this%n-6 ,k)  &
                                                  + b2_a6_d15 * f(:,  this%n-7 ,k)  &
                                                  + b2_a7_d15 * f(:,  this%n-8 ,k) 
            RHS(:,  this%n  ,k) = a0_d15 * ( f(:,  this%n  ,k) )                    
         end do

         case(1)
            do k=1,n3
               RHS(:,  this%n-6,k) = a0_d15 * ( f(:,  this%n-6,k) )        &
                    + a1_d15 * ( f(:,   this%n-5,k) + f(:,  this%n-7 ,k) ) &
                    + a2_d15 * ( f(:,   this%n-4,k) + f(:,  this%n-8 ,k) ) &
                    + a3_d15 * ( f(:,   this%n-3,k) + f(:,  this%n-9 ,k) ) &
                    + a4_d15 * ( f(:,   this%n-2,k) + f(:,  this%n-10,k) ) &
                    + a5_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-11,k) ) &
                    + a6_d15 * ( f(:,   this%n  ,k) + f(:,  this%n-12,k) ) &
                    + a7_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-13,k) )
               RHS(:,  this%n-5,k) = a0_d15 * ( f(:,  this%n-5,k) )        &
                    + a1_d15 * ( f(:,   this%n-4,k) + f(:,  this%n-6 ,k) ) &
                    + a2_d15 * ( f(:,   this%n-3,k) + f(:,  this%n-7 ,k) ) &
                    + a3_d15 * ( f(:,   this%n-2,k) + f(:,  this%n-8 ,k) ) &
                    + a4_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-9 ,k) ) &
                    + a5_d15 * ( f(:,   this%n  ,k) + f(:,  this%n-10,k) ) &
                    + a6_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-11,k) ) &
                    + a7_d15 * ( f(:,   this%n-2,k) + f(:,  this%n-12,k) )
               RHS(:,  this%n-4,k) = a0_d15 * ( f(:,  this%n-4,k) )        &
                    + a1_d15 * ( f(:,   this%n-3,k) + f(:,  this%n-5 ,k) ) &
                    + a2_d15 * ( f(:,   this%n-2,k) + f(:,  this%n-6 ,k) ) &
                    + a3_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-7 ,k) ) &
                    + a4_d15 * ( f(:,   this%n  ,k) + f(:,  this%n-8 ,k) ) &
                    + a5_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-9 ,k) ) &
                    + a6_d15 * ( f(:,   this%n-2,k) + f(:,  this%n-10,k) ) &
                    + a7_d15 * ( f(:,   this%n-3,k) + f(:,  this%n-11,k) )
               RHS(:,  this%n-3,k) = a0_d15 * ( f(:,  this%n-3,k) )        &
                    + a1_d15 * ( f(:,   this%n-2,k) + f(:,  this%n-4 ,k) ) &
                    + a2_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-5 ,k) ) &
                    + a3_d15 * ( f(:,   this%n  ,k) + f(:,  this%n-6 ,k) ) &
                    + a4_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-7 ,k) ) &
                    + a5_d15 * ( f(:,   this%n-2,k) + f(:,  this%n-8 ,k) ) &
                    + a6_d15 * ( f(:,   this%n-3,k) + f(:,  this%n-9 ,k) ) &
                    + a7_d15 * ( f(:,   this%n-4,k) + f(:,  this%n-10,k) )
               RHS(:,  this%n-2,k) = a0_d15 * ( f(:,  this%n-2,k) )        &
                    + a1_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-3 ,k) ) &
                    + a2_d15 * ( f(:,   this%n  ,k) + f(:,  this%n-4 ,k) ) &
                    + a3_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-5 ,k) ) &
                    + a4_d15 * ( f(:,   this%n-2,k) + f(:,  this%n-6 ,k) ) &
                    + a5_d15 * ( f(:,   this%n-3,k) + f(:,  this%n-7 ,k) ) &
                    + a6_d15 * ( f(:,   this%n-4,k) + f(:,  this%n-8 ,k) ) &
                    + a7_d15 * ( f(:,   this%n-5,k) + f(:,  this%n-9 ,k) )
               RHS(:,  this%n-1,k) = a0_d15 * ( f(:,  this%n-1,k) )        &
                    + a1_d15 * ( f(:,   this%n  ,k) + f(:,  this%n-2 ,k) ) &
                    + a2_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-3 ,k) ) &
                    + a3_d15 * ( f(:,   this%n-2,k) + f(:,  this%n-4 ,k) ) &
                    + a4_d15 * ( f(:,   this%n-3,k) + f(:,  this%n-5 ,k) ) &
                    + a5_d15 * ( f(:,   this%n-4,k) + f(:,  this%n-6 ,k) ) &
                    + a6_d15 * ( f(:,   this%n-5,k) + f(:,  this%n-7 ,k) ) &
                    + a7_d15 * ( f(:,   this%n-6,k) + f(:,  this%n-8 ,k) )
               RHS(:,  this%n  ,k) = a0_d15 * ( f(:,  this%n  ,k) )        &
                    + a1_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-1 ,k) ) &
                    + a2_d15 * ( f(:,   this%n-2,k) + f(:,  this%n-2 ,k) ) &
                    + a3_d15 * ( f(:,   this%n-3,k) + f(:,  this%n-3 ,k) ) &
                    + a4_d15 * ( f(:,   this%n-4,k) + f(:,  this%n-4 ,k) ) &
                    + a5_d15 * ( f(:,   this%n-5,k) + f(:,  this%n-5 ,k) ) &
                    + a6_d15 * ( f(:,   this%n-6,k) + f(:,  this%n-6 ,k) ) &
                    + a7_d15 * ( f(:,   this%n-7,k) + f(:,  this%n-7 ,k) )
            end do

         case(-1)
            do k=1,n3
               RHS(:,  this%n-6,k) = a0_d15 * ( f(:,  this%n-6,k) )        &
                    + a1_d15 * ( f(:,   this%n-5,k) + f(:,  this%n-7 ,k) ) &
                    + a2_d15 * ( f(:,   this%n-4,k) + f(:,  this%n-8 ,k) ) &
                    + a3_d15 * ( f(:,   this%n-3,k) + f(:,  this%n-9 ,k) ) &
                    + a4_d15 * ( f(:,   this%n-2,k) + f(:,  this%n-10,k) ) &
                    + a5_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-11,k) ) &
                    + a6_d15 * ( f(:,   this%n  ,k) + f(:,  this%n-12,k) ) &
                    + a7_d15 * (-f(:,   this%n-1,k) + f(:,  this%n-13,k) )
               RHS(:,  this%n-5,k) = a0_d15 * ( f(:,  this%n-5,k) )        &
                    + a1_d15 * ( f(:,   this%n-4,k) + f(:,  this%n-6 ,k) ) &
                    + a2_d15 * ( f(:,   this%n-3,k) + f(:,  this%n-7 ,k) ) &
                    + a3_d15 * ( f(:,   this%n-2,k) + f(:,  this%n-8 ,k) ) &
                    + a4_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-9 ,k) ) &
                    + a5_d15 * ( f(:,   this%n  ,k) + f(:,  this%n-10,k) ) &
                    + a6_d15 * (-f(:,   this%n-1,k) + f(:,  this%n-11,k) ) &
                    + a7_d15 * (-f(:,   this%n-2,k) + f(:,  this%n-12,k) )
               RHS(:,  this%n-4,k) = a0_d15 * ( f(:,  this%n-4,k) )        &
                    + a1_d15 * ( f(:,   this%n-3,k) + f(:,  this%n-5 ,k) ) &
                    + a2_d15 * ( f(:,   this%n-2,k) + f(:,  this%n-6 ,k) ) &
                    + a3_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-7 ,k) ) &
                    + a4_d15 * ( f(:,   this%n  ,k) + f(:,  this%n-8 ,k) ) &
                    + a5_d15 * (-f(:,   this%n-1,k) + f(:,  this%n-9 ,k) ) &
                    + a6_d15 * (-f(:,   this%n-2,k) + f(:,  this%n-10,k) ) &
                    + a7_d15 * (-f(:,   this%n-3,k) + f(:,  this%n-11,k) )
               RHS(:,  this%n-3,k) = a0_d15 * ( f(:,  this%n-3,k) )        &
                    + a1_d15 * ( f(:,   this%n-2,k) + f(:,  this%n-4 ,k) ) &
                    + a2_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-5 ,k) ) &
                    + a3_d15 * ( f(:,   this%n  ,k) + f(:,  this%n-6 ,k) ) &
                    + a4_d15 * (-f(:,   this%n-1,k) + f(:,  this%n-7 ,k) ) &
                    + a5_d15 * (-f(:,   this%n-2,k) + f(:,  this%n-8 ,k) ) &
                    + a6_d15 * (-f(:,   this%n-3,k) + f(:,  this%n-9 ,k) ) &
                    + a7_d15 * (-f(:,   this%n-4,k) + f(:,  this%n-10,k) )
               RHS(:,  this%n-2,k) = a0_d15 * ( f(:,  this%n-2,k) )        &
                    + a1_d15 * ( f(:,   this%n-1,k) + f(:,  this%n-3 ,k) ) &
                    + a2_d15 * ( f(:,   this%n  ,k) + f(:,  this%n-4 ,k) ) &
                    + a3_d15 * (-f(:,   this%n-1,k) + f(:,  this%n-5 ,k) ) &
                    + a4_d15 * (-f(:,   this%n-2,k) + f(:,  this%n-6 ,k) ) &
                    + a5_d15 * (-f(:,   this%n-3,k) + f(:,  this%n-7 ,k) ) &
                    + a6_d15 * (-f(:,   this%n-4,k) + f(:,  this%n-8 ,k) ) &
                    + a7_d15 * (-f(:,   this%n-5,k) + f(:,  this%n-9 ,k) )
               RHS(:,  this%n-1,k) = a0_d15 * ( f(:,  this%n-1,k) )        &
                    + a1_d15 * ( f(:,   this%n  ,k) + f(:,  this%n-2 ,k) ) &
                    + a2_d15 * (-f(:,   this%n-1,k) + f(:,  this%n-3 ,k) ) &
                    + a3_d15 * (-f(:,   this%n-2,k) + f(:,  this%n-4 ,k) ) &
                    + a4_d15 * (-f(:,   this%n-3,k) + f(:,  this%n-5 ,k) ) &
                    + a5_d15 * (-f(:,   this%n-4,k) + f(:,  this%n-6 ,k) ) &
                    + a6_d15 * (-f(:,   this%n-5,k) + f(:,  this%n-7 ,k) ) &
                    + a7_d15 * (-f(:,   this%n-6,k) + f(:,  this%n-8 ,k) )
               RHS(:,  this%n  ,k) = a0_d15 * ( f(:,  this%n  ,k) )        &
                    + a1_d15 * (-f(:,   this%n-1,k) + f(:,  this%n-1 ,k) ) &
                    + a2_d15 * (-f(:,   this%n-2,k) + f(:,  this%n-2 ,k) ) &
                    + a3_d15 * (-f(:,   this%n-3,k) + f(:,  this%n-3 ,k) ) &
                    + a4_d15 * (-f(:,   this%n-4,k) + f(:,  this%n-4 ,k) ) &
                    + a5_d15 * (-f(:,   this%n-5,k) + f(:,  this%n-5 ,k) ) &
                    + a6_d15 * (-f(:,   this%n-6,k) + f(:,  this%n-6 ,k) ) &
                    + a7_d15 * (-f(:,   this%n-7,k) + f(:,  this%n-7 ,k) )
            end do

         end select
      end select

 end subroutine ComputeYRHS

 pure subroutine ComputeZRHS_REAL(this, f, RHS, n1, n2, bc1, bcn)

   class( cfo2D15 ), intent(in) :: this
   integer, intent(in) :: n1, n2
   real(rkind), dimension(n1,n2,this%n), intent(in) :: f
   real(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
   integer, intent(in) :: bc1, bcn

#include "CFo2D15_files/ComputeZRHS_common.F90"    

 end subroutine ComputeZRHS_REAL

 pure subroutine ComputeZRHS_CMPLX(this, f, RHS, n1, n2, bc1, bcn)

   class( cfo2D15 ), intent(in) :: this
   integer, intent(in) :: n1, n2
   complex(rkind), dimension(n1,n2,this%n), intent(in) :: f
   complex(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
   integer, intent(in) :: bc1, bcn

#include "CFo2D15_files/ComputeZRHS_common.F90"    

 end subroutine ComputeZRHS_CMPLX

 subroutine filter1(this, f, df, na, nb, bc1_, bcn_)
   class( cfo2D15 ), intent(in) :: this
   integer, intent(in) :: na, nb
   real(rkind), dimension(this%n,na,nb), intent(in) :: f
   real(rkind), dimension(this%n,na,nb), intent(out) :: df
   integer, optional, intent(in) :: bc1_, bcn_
   integer :: bc1, bcn

   if(this%n == 1) then
      df = f
      return
   end if

   if (present(bc1_)) then
      bc1 = bc1_
      if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
         call GracefulExit("Incorrect boundary specification for bc1 (should be 0, 1 or -1)", 324)
      end if
   else
      bc1 = 0
   end if

   if (present(bcn_)) then
      bcn = bcn_
      if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
         call GracefulExit("Incorrect boundary specification for bcn (should be 0, 1 or -1)", 324)
      end if
   else
      bcn = 0
   end if

   call this%ComputeXRHS(f, df, na, nb, bc1, bcn)

   select case (this%periodic)
   case(.TRUE.)
      call this%SolveXLU(df, na, nb)
   case(.FALSE.)
      select case(bc1)
      case(0) ! Normal non-periodic left boundary
         select case(bcn)
         case(0)  ! Normal non-periodic right boundary
            call this%SolveXD15(this%d15_nn, this%d15_ipiv_nn, df, na, nb)
         case(1) ! Symmetric right boundary
            call this%SolveXD15(this%d15_ns, this%d15_ipiv_ns, df, na, nb)
         case(-1) ! Antisymmetric right boundary
            call this%SolveXD15(this%d15_na, this%d15_ipiv_na, df, na, nb)
         end select
      case(1) ! Symmetric left boundary
         select case(bcn)
         case(0)  ! Normal non-periodic right boundary
            call this%SolveXD15(this%d15_sn, this%d15_ipiv_sn, df, na, nb)
         case(1)  ! Symmetric right boundary
            call this%SolveXD15(this%d15_ss, this%d15_ipiv_ss, df, na, nb)
         case(-1) ! Antisymmetric right boundary
            call this%SolveXD15(this%d15_sa, this%d15_ipiv_sa, df, na, nb)
         end select
      case(-1) ! Antisymmetric left boundary
         select case(bcn)
         case(0)  ! Normal non-periodic right boundary
            call this%SolveXD15(this%d15_an, this%d15_ipiv_an, df, na, nb)
         case(1) ! Symmetric right boundary
            call this%SolveXD15(this%d15_as, this%d15_ipiv_as, df, na, nb)
         case(-1) ! Antisymmetric right boundary
            call this%SolveXD15(this%d15_aa, this%d15_ipiv_aa, df, na, nb)
         end select
      end select
   end select

 end subroutine filter1

 subroutine filter2(this, f, df, na, nb, bc1_, bcn_)
   class( cfo2D15 ), intent(in) :: this
   integer, intent(in) :: na, nb
   real(rkind), dimension(na,this%n,nb), intent(in) :: f
   real(rkind), dimension(na,this%n,nb), intent(out) :: df
   integer, optional, intent(in) :: bc1_, bcn_
   integer :: bc1, bcn

   if(this%n == 1) then
      df = f
      return
   end if

   if (present(bc1_)) then
      bc1 = bc1_
      if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
         call GracefulExit("Incorrect boundary specification for bc1 (should be 0, 1 or -1)", 324)
      end if
   else
      bc1 = 0
   end if

   if (present(bcn_)) then
      bcn = bcn_
      if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
         call GracefulExit("Incorrect boundary specification for bcn (should be 0, 1 or -1)", 324)
      end if
   else
      bcn = 0
   end if

   call this%ComputeYRHS(f, df, na, nb, bc1, bcn)

   select case (this%periodic)
   case(.TRUE.)
      call this%SolveYLU(df, na, nb)
   case(.FALSE.)
      select case(bc1)
      case(0) ! Normal non-periodic left boundary
         select case(bcn)
         case(0)  ! Normal non-periodic right boundary
            call this%SolveYD15(this%d15_nn, this%d15_ipiv_nn, df, na, nb)
         case(1) ! Symmetric right boundary
            call this%SolveYD15(this%d15_ns, this%d15_ipiv_ns, df, na, nb)
         case(-1) ! Antisymmetric right boundary
            call this%SolveYD15(this%d15_na, this%d15_ipiv_na, df, na, nb)
         end select
      case(1) ! Symmetric left boundary
         select case(bcn)
         case(0)  ! Normal non-periodic right boundary
            call this%SolveYD15(this%d15_sn, this%d15_ipiv_sn, df, na, nb)
         case(1)  ! Symmetric right boundary
            call this%SolveYD15(this%d15_ss, this%d15_ipiv_ss, df, na, nb)
         case(-1) ! Antisymmetric right boundary
            call this%SolveYD15(this%d15_sa, this%d15_ipiv_sa, df, na, nb)
         end select
      case(-1) ! Antisymmetric left boundary
         select case(bcn)
         case(0)  ! Normal non-periodic right boundary
            call this%SolveYD15(this%d15_an, this%d15_ipiv_an, df, na, nb)
         case(1) ! Symmetric right boundary
            call this%SolveYD15(this%d15_as, this%d15_ipiv_as, df, na, nb)
         case(-1) ! Antisymmetric right boundary
            call this%SolveYD15(this%d15_aa, this%d15_ipiv_aa, df, na, nb)
         end select
      end select
   end select

 end subroutine filter2

 subroutine filter3_REAL(this, f, df, na, nb, bc1_, bcn_)
   class( cfo2D15 ), intent(in) :: this
   integer, intent(in) :: na, nb
   real(rkind), dimension(na,nb,this%n), intent(in) :: f
   real(rkind), dimension(na,nb,this%n), intent(out) :: df
   integer, optional, intent(in) :: bc1_, bcn_
   integer :: bc1, bcn

   if(this%n == 1) then
      df = f
      return
   end if

   if (present(bc1_)) then
      bc1 = bc1_
      if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
         call GracefulExit("Incorrect boundary specification for bc1 (should be 0, 1 or -1)", 324)
      end if
   else
      bc1 = 0
   end if

   if (present(bcn_)) then
      bcn = bcn_
      if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
         call GracefulExit("Incorrect boundary specification for bcn (should be 0, 1 or -1)", 324)
      end if
   else
      bcn = 0
   end if

   call this%ComputeZRHS_REAL(f, df, na, nb, bc1, bcn)

   select case (this%periodic)
   case(.TRUE.)
      call this%SolveZLU(df, na, nb)
   case(.FALSE.)
      select case(bc1)
      case(0) ! Normal non-periodic left boundary
         select case(bcn)
         case(0)  ! Normal non-periodic right boundary
            call this%SolveZD15_REAL(this%d15_nn, this%d15_ipiv_nn, df, na, nb)
         case(1) ! Symmetric right boundary
            call this%SolveZD15_REAL(this%d15_ns, this%d15_ipiv_ns, df, na, nb)
         case(-1) ! Antisymmetric right boundary
            call this%SolveZD15_REAL(this%d15_na, this%d15_ipiv_na, df, na, nb)
         end select
      case(1) ! Symmetric left boundary
         select case(bcn)
         case(0)  ! Normal non-periodic right boundary
            call this%SolveZD15_REAL(this%d15_sn, this%d15_ipiv_sn, df, na, nb)
         case(1)  ! Symmetric right boundary
            call this%SolveZD15_REAL(this%d15_ss, this%d15_ipiv_ss, df, na, nb)
         case(-1) ! Antisymmetric right boundary
            call this%SolveZD15_REAL(this%d15_sa, this%d15_ipiv_sa, df, na, nb)
         end select
      case(-1) ! Antisymmetric left boundary
         select case(bcn)
         case(0)  ! Normal non-periodic right boundary
            call this%SolveZD15_REAL(this%d15_an, this%d15_ipiv_an, df, na, nb)
         case(1) ! Symmetric right boundary
            call this%SolveZD15_REAL(this%d15_as, this%d15_ipiv_as, df, na, nb)
         case(-1) ! Antisymmetric right boundary
            call this%SolveZD15_REAL(this%d15_aa, this%d15_ipiv_aa, df, na, nb)
         end select
      end select
   end select

 end subroutine filter3_REAL

 subroutine filter3_CMPLX(this, f, df, na, nb, bc1_, bcn_)
   class( cfo2D15 ), intent(in) :: this
   integer, intent(in) :: na, nb
   complex(rkind), dimension(na,nb,this%n), intent(in) :: f
   complex(rkind), dimension(na,nb,this%n), intent(out) :: df
   integer, optional, intent(in) :: bc1_, bcn_
   integer :: bc1, bcn

   if(this%n == 1) then
      df = f
      return
   end if
   if (present(bcn_)) then
      bc1 = bc1_; bcn = bcn_
   else 
      bc1 = 0; bcn = 0
   end if

   call this%ComputeZRHS_CMPLX(f, df, na, nb, bc1, bcn)
   call this%SolveZD15_CMPLX(this%d15_nn, this%d15_ipiv_nn, df, na, nb)

 end subroutine filter3_CMPLX

 end module
