! Routines specific to 2nd order Explicit Midpoint Interpolation Scheme

module ei06stuff

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two, half
    use exits,           only: GracefulExit
    
    implicit none

    private
    public :: ei06, aI02 !TODO Rename this last var
    ! 2nd order first derivative explicit centeral difference coefficients
    real(rkind), parameter :: aI02     =  75.0d0/64.0d0
    real(rkind), parameter :: bI02     = -25.0d0/128.0d0
    real(rkind), parameter :: cI02     =   3.0d0/128.0d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic 1st derivative evaluation !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !TODO: implment boundary scheme 
    ! Interior scheme
    real(rkind), parameter                   :: q_hat           = aI02

    ! Scheme at node 1 (boundary): explicit one-sided 2nd order
    real(rkind), parameter                   :: p               = -3.0_rkind / 2.0_rkind
    real(rkind), parameter                   :: q               = 4.0_rkind / 2.0_rkind
    real(rkind), parameter                   :: r               = -1.0_rkind / 2.0_rkind
    real(rkind), parameter                   :: s               = 0.0_rkind / 2.0_rkind

    ! Scheme at node 2: explicit central 2nd order
    real(rkind), parameter                   :: q_p             = 1.0_rkind / 2.0_rkind
   
    ! Scheme at node 3: explicit central 2nd order
    real(rkind), parameter                   :: q_pp            = 1.0_rkind / 2.0_rkind

    ! Scheme at node 4: same as interior
    real(rkind), parameter                   :: q_ppp           = aI02

    ! Weights are now trivial, but fluxes do not telescope
    real(rkind), parameter                   :: w1              = 1.0_rkind
    real(rkind), parameter                   :: w2              = 1.0_rkind
    real(rkind), parameter                   :: w3              = 1.0_rkind
    real(rkind), parameter                   :: w4              = 1.0_rkind
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type ei06
        
        private
        
        integer     :: n
        real(rkind) :: dx
        real(rkind) :: onebydx
        real(rkind) :: onebydx2

        logical     :: periodic=.TRUE.
        integer     :: bc1=0                               ! Boundary condition type. 0=Dirichlet, 1=Neumann
        integer     :: bcn=0                               ! Boundary condition type. 0=Dirichlet, 1=Neumann 

        contains

        procedure :: init
        procedure :: destroy
        procedure :: GetSize
        procedure, private :: ComputeXIRHS !TODO: change these names to be I123?
        procedure, private :: ComputeYIRHS
        procedure, private :: ComputeZIRHS
        
        procedure :: iN2F1 
        procedure :: iN2F2
        procedure :: iN2F3

        procedure :: iF2N1 
        procedure :: iF2N2 
        procedure :: iF2N3 
        
    end type



contains
    
    pure function GetSize(this) result(val)
        class(ei06), intent(in) :: this
        integer  :: val 
        val = this%n
    end function

    function init(this, n_, dx_, periodic_, bc1_, bcn_) result(ierr)
   
        class( ei06 ), intent(inout) :: this
        integer, intent(in) :: n_
        real(rkind), intent(in) :: dx_
        logical, intent(in) :: periodic_
        integer, intent(in) :: bc1_, bcn_
        integer :: ierr
        
        this%n = n_
        this%dx = dx_
        this%onebydx = one/dx_
        this%onebydx2 = this%onebydx/dx_

        this%periodic = periodic_

        this%bc1 = bc1_
        this%bcn = bcn_

        ! If everything passes
        ierr = 0
    
    end function
    
    subroutine destroy(this)

        class( ei06 ), intent(inout) :: this

        return

    end subroutine


    pure subroutine ComputeXIRHS(this, f, RHS, dir, n2, n3, bc1, bcn) !TODO add an extra argument to account for forward or backward
    
        class( ei06 ), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10, b10, c10, a101, a102, a104, b104
        integer :: j,k
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4
        real(rkind) :: a_np_3
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1

        select case (this%periodic)
        case (.TRUE.)
            a10 = aI02 * half
            b10 = bI02 * half
            c10 = cI02 * half
            RHS = 0.0d0
            select case (dir)
                case ("N2F")
                    !interior    
                    RHS(3:this%n-3,:,:) = RHS(3:this%n-3,:,:) + a10 * (f(4:this%n-2,:,:) + f(3:this%n-3,:,:)) &
                                                              + b10 * (f(5:this%n-1,:,:) + f(2:this%n-4,:,:)) &
                                                              + c10 * (f(6:this%n,  :,:) + f(1:this%n-5,:,:))

                    !left boundary
                    RHS(1,:,:) = RHS(1,:,:) + a10 * (f(2,:,:) + f(1,       :,:)) &
                                            + b10 * (f(3,:,:) + f(this%n,  :,:)) &
                                            + c10 * (f(4,:,:) + f(this%n-1,:,:))
                    RHS(2,:,:) = RHS(2,:,:) + a10 * (f(3,:,:) + f(2,       :,:)) &
                                            + b10 * (f(4,:,:) + f(1,       :,:)) &
                                            + c10 * (f(5,:,:) + f(this%n,  :,:))
                    !right boundary (n-2:n)
                    RHS(this%n-2,:,:) = RHS(this%n-2,:,:) + a10 * (f(this%n-1,:,:) + f(this%n-2,  :,:)) &
                                                          + b10 * (f(this%n,  :,:) + f(this%n-3,  :,:)) &
                                                          + c10 * (f(1,       :,:) + f(this%n-4,  :,:))
                    RHS(this%n-1,:,:) = RHS(this%n-1,:,:) + a10 * (f(this%n,  :,:) + f(this%n-1,  :,:)) &
                                                          + b10 * (f(1,       :,:) + f(this%n-2,  :,:)) &
                                                          + c10 * (f(2,       :,:) + f(this%n-3,  :,:))
                    RHS(this%n,  :,:) = RHS(this%n,  :,:) + a10 * (f(1,       :,:) + f(this%n,    :,:)) &
                                                          + b10 * (f(2,       :,:) + f(this%n-1,  :,:)) &
                                                          + c10 * (f(3,       :,:) + f(this%n-2,  :,:))

                case ("F2N")
                    !interior    
                    RHS(4:this%n-2,:,:) = RHS(4:this%n-2,:,:) + a10 * (f(4:this%n-2,:,:) + f(3:this%n-3,:,:)) &
                                                              + b10 * (f(5:this%n-1,:,:) + f(2:this%n-4,:,:)) &
                                                              + c10 * (f(6:this%n,  :,:) + f(1:this%n-5,:,:))

                    !left boundary
                    RHS(1,:,:) = RHS(1,:,:) + a10 * (f(1,:,:) + f(this%n,  :,:)) &
                                            + b10 * (f(2,:,:) + f(this%n-1,:,:)) &
                                            + c10 * (f(3,:,:) + f(this%n-2,:,:))
                    RHS(2,:,:) = RHS(2,:,:) + a10 * (f(2,:,:) + f(1,       :,:)) &
                                            + b10 * (f(3,:,:) + f(this%n,  :,:)) &
                                            + c10 * (f(4,:,:) + f(this%n-1,:,:))
                    RHS(3,:,:) = RHS(3,:,:) + a10 * (f(3,:,:) + f(2,       :,:)) &
                                            + b10 * (f(4,:,:) + f(1,       :,:)) &
                                            + c10 * (f(5,:,:) + f(this%n,  :,:))
                    !right boundary (n-1:n)
                    RHS(this%n-1,:,:) = RHS(this%n-1,:,:) + a10 * (f(this%n-1,:,:) + f(this%n-2,  :,:)) &
                                                          + b10 * (f(this%n,  :,:) + f(this%n-3,  :,:)) &
                                                          + c10 * (f(1,       :,:) + f(this%n-4,  :,:))
                    RHS(this%n,  :,:) = RHS(this%n,  :,:) + a10 * (f(this%n,  :,:) + f(this%n-1,  :,:)) &
                                                          + b10 * (f(1,       :,:) + f(this%n-2,  :,:)) &
                                                          + c10 * (f(2,       :,:) + f(this%n-3,  :,:))
            end select

        case (.FALSE.)  
            a10 = aI02 * half
            b10 = bI02 * half
            c10 = cI02 * half
            a101 =  one
            a102 =  half
            a104 =  9.0d0/16.0d0
            b104 = -1.0d0/16.0d0
            RHS = 0.0d0
            select case (dir)
                case ("N2F")!TODO: implement better non-periodic BC: currently 2466...6642
                    !interior    
                    RHS(3:this%n-3,:,:) = RHS(3:this%n-3,:,:) + a10 * (f(4:this%n-2,:,:) + f(3:this%n-3,:,:)) &
                                                              + b10 * (f(5:this%n-1,:,:) + f(2:this%n-4,:,:)) &
                                                              + c10 * (f(6:this%n,  :,:) + f(1:this%n-5,:,:))
                    select case(bc1)
                        !left boundary (1:2)
                        case(1) !symm
                            RHS(1,:,:) = RHS(1,:,:) + a10 * (f(2,:,:) + f(1,:,:)) &
                                                    + b10 * (f(3,:,:) + f(2,:,:)) &
                                                    + c10 * (f(4,:,:) + f(3,:,:))
                            RHS(2,:,:) = RHS(2,:,:) + a10 * (f(3,:,:) + f(2,:,:)) &
                                                    + b10 * (f(4,:,:) + f(1,:,:)) &
                                                    + c10 * (f(5,:,:) + f(2,:,:))
                        case(-1) !anti-symm
                            RHS(1,:,:) = RHS(1,:,:) + a10 * (f(2,:,:) + f(1,:,:)) &
                                                    + b10 * (f(3,:,:) - f(2,:,:)) &
                                                    + c10 * (f(4,:,:) - f(3,:,:))
                            RHS(2,:,:) = RHS(2,:,:) + a10 * (f(3,:,:) + f(2,:,:)) &
                                                    + b10 * (f(4,:,:) + f(1,:,:)) &
                                                    + c10 * (f(5,:,:) - f(2,:,:))
                        case(0)
                            RHS(1,:,:) = RHS(1,:,:) + a102 * (f(2,:,:) + f(1,:,:))   !2nd order
                            RHS(2,:,:) = RHS(2,:,:) + a104 * (f(3,:,:) + f(2,:,:)) & !4th order
                                                    + b104 * (f(4,:,:) + f(1,:,:)) 
                    end select

                    select case(bcn)
                        !right boundary (n-2:n-1)
                        case(1)
                            RHS(this%n-2,:,:) = RHS(this%n-2,:,:) + a10 * (f(this%n-1,:,:) + f(this%n-2,:,:)) &
                                                                  + b10 * (f(this%n  ,:,:) + f(this%n-3,:,:)) &
                                                                  + c10 * (f(this%n-1,:,:) + f(this%n-4,:,:))
                            RHS(this%n-1,:,:) = RHS(this%n-1,:,:) + a10 * (f(this%n  ,:,:) + f(this%n-1,:,:)) &
                                                                  + b10 * (f(this%n-1,:,:) + f(this%n-2,:,:)) &
                                                                  + c10 * (f(this%n-2,:,:) + f(this%n-3,:,:))
                        case(-1)
                            RHS(this%n-2,:,:) = RHS(this%n-2,:,:) + a10 * ( f(this%n-1,:,:) + f(this%n-2,:,:)) &
                                                                  + b10 * ( f(this%n  ,:,:) + f(this%n-3,:,:)) &
                                                                  + c10 * (-f(this%n-1,:,:) + f(this%n-4,:,:))
                            RHS(this%n-1,:,:) = RHS(this%n-1,:,:) + a10 * ( f(this%n  ,:,:) + f(this%n-1,:,:)) &
                                                                  + b10 * (-f(this%n-1,:,:) + f(this%n-2,:,:)) &
                                                                  + c10 * (-f(this%n-2,:,:) + f(this%n-3,:,:))
                        case(0)
                            RHS(this%n-2,:,:) = RHS(this%n-2,:,:) + a104 * (f(this%n-1,:,:) + f(this%n-2,  :,:)) & !4th order
                                                                  + b104 * (f(this%n,  :,:) + f(this%n-3,  :,:))   
                            RHS(this%n-1,:,:) = RHS(this%n-1,:,:) + a102 * (f(this%n,  :,:) + f(this%n-1,  :,:))   !2nd order
                    end select

                case ("F2N")!TODO: implement better non-periodic BC: currently 12466...66421
                    !interior    
                    RHS(4:this%n-3,:,:) = RHS(4:this%n-3,:,:) + a10 * (f(4:this%n-3,:,:) + f(3:this%n-4,:,:)) &
                                                              + b10 * (f(5:this%n-2,:,:) + f(2:this%n-5,:,:)) &
                                                              + c10 * (f(6:this%n-1,:,:) + f(1:this%n-6,:,:))

                    select case(bc1)
                        ! left boundary (1:3)
                        case(1)
                            RHS(3,:,:) = RHS(3,:,:) + a10 * (f(3,:,:) + f(2,:,:)) &
                                                    + b10 * (f(4,:,:) + f(1,:,:)) &
                                                    + c10 * (f(5,:,:) + f(1,:,:))
                            RHS(2,:,:) = RHS(2,:,:) + a10 * (f(2,:,:) + f(1,:,:)) &
                                                    + b10 * (f(3,:,:) + f(1,:,:)) &
                                                    + c10 * (f(4,:,:) + f(2,:,:))
                            RHS(1,:,:) = RHS(1,:,:) + a10 * (f(1,:,:) + f(1,:,:)) &
                                                    + b10 * (f(2,:,:) + f(2,:,:)) &
                                                    + c10 * (f(3,:,:) + f(3,:,:))
                        case(-1)
                            RHS(3,:,:) = RHS(3,:,:) + a10 * (f(3,:,:) + f(2,:,:)) &
                                                    + b10 * (f(4,:,:) + f(1,:,:)) &
                                                    + c10 * (f(5,:,:) - f(1,:,:))
                            RHS(2,:,:) = RHS(2,:,:) + a10 * (f(2,:,:) + f(1,:,:)) &
                                                    + b10 * (f(3,:,:) - f(1,:,:)) &
                                                    + c10 * (f(4,:,:) - f(2,:,:))
                            RHS(1,:,:) = RHS(1,:,:) + a10 * (f(1,:,:) - f(1,:,:)) &
                                                    + b10 * (f(2,:,:) - f(2,:,:)) &
                                                    + c10 * (f(3,:,:) - f(3,:,:))
                        case(0)
                            RHS(1,:,:) = RHS(1,:,:) + a101 * (f(1,:,:)           )
                            RHS(2,:,:) = RHS(2,:,:) + a102 * (f(2,:,:) + f(1,:,:)) 
                            RHS(3,:,:) = RHS(3,:,:) + a104 * (f(3,:,:) + f(2,:,:)) &
                                                    + b104 * (f(4,:,:) + f(1,:,:)) 
                    end select

                    select case(bcn)
                        !right boundary (n-2:n)
                        case(1)
                            RHS(this%n-2,:,:) = RHS(this%n-2,:,:) + a10 * (f(this%n-2,:,:) + f(this%n-3,:,:)) &
                                                                  + b10 * (f(this%n-1,:,:) + f(this%n-4,:,:)) &
                                                                  + c10 * (f(this%n-1,:,:) + f(this%n-5,:,:))
                            RHS(this%n-1,:,:) = RHS(this%n-1,:,:) + a10 * (f(this%n-1,:,:) + f(this%n-2,:,:)) &
                                                                  + b10 * (f(this%n-1,:,:) + f(this%n-3,:,:)) &
                                                                  + c10 * (f(this%n-2,:,:) + f(this%n-4,:,:))
                            RHS(this%n,:,:) = RHS(this%n,:,:)     + a10 * (f(this%n-1,:,:) + f(this%n-1,:,:)) &
                                                                  + b10 * (f(this%n-2,:,:) + f(this%n-2,:,:)) &
                                                                  + c10 * (f(this%n-3,:,:) + f(this%n-3,:,:))
                        case(-1)
                            RHS(this%n-2,:,:) = RHS(this%n-2,:,:) + a10 * ( f(this%n-2,:,:) + f(this%n-3,:,:)) &
                                                                  + b10 * ( f(this%n-1,:,:) + f(this%n-4,:,:)) &
                                                                  + c10 * (-f(this%n-1,:,:) + f(this%n-5,:,:))
                            RHS(this%n-1,:,:) = RHS(this%n-1,:,:) + a10 * ( f(this%n-1,:,:) + f(this%n-2,:,:)) &
                                                                  + b10 * (-f(this%n-1,:,:) + f(this%n-3,:,:)) &
                                                                  + c10 * (-f(this%n-2,:,:) + f(this%n-4,:,:))
                            RHS(this%n,:,:) = RHS(this%n,:,:)     + a10 * (-f(this%n-1,:,:) + f(this%n-1,:,:)) &
                                                                  + b10 * (-f(this%n-2,:,:) + f(this%n-2,:,:)) &
                                                                  + c10 * (-f(this%n-3,:,:) + f(this%n-3,:,:))
                        case(0)
                            RHS(this%n,:,:)   = RHS(this%n  ,:,:) + a101 * (f(this%n-1,:,:)                  )
                            RHS(this%n-1,:,:) = RHS(this%n-1,:,:) + a102 * (f(this%n-1,:,:) + f(this%n-2,:,:)) 
                            RHS(this%n-2,:,:) = RHS(this%n-2,:,:) + a104 * (f(this%n-2,:,:) + f(this%n-3,:,:)) &
                                                                  + b104 * (f(this%n-1,:,:) + f(this%n-4,:,:)) 
                    end select

            end select

        end select
    
    end subroutine
    
    pure subroutine ComputeYIRHS(this, f, RHS, dir, n1, n3, bc1, bcn) 
    
        class( ei06 ), intent(in) :: this
        integer, intent(in) :: n1, n3
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10, b10, c10, a101, a102, a104, b104
        integer :: i,k
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4
        real(rkind) :: a_np_3
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1


        select case (this%periodic)
        case (.TRUE.)
            a10 = aI02 * half
            b10 = bI02 * half
            c10 = cI02 * half
            RHS = 0.0d0
            select case (dir)
                case ("N2F")
                    !interior    
                    RHS(:,3:this%n-3,:) = RHS(:,3:this%n-3,:) + a10 * (f(:,4:this%n-2,:) + f(:,3:this%n-3,:)) &
                                                              + b10 * (f(:,5:this%n-1,:) + f(:,2:this%n-4,:)) &
                                                              + c10 * (f(:,6:this%n,  :) + f(:,1:this%n-5,:))

                    !left boundary
                    RHS(:,1,:) = RHS(:,1,:) + a10 * (f(:,2,:) + f(:,1,       :)) &
                                            + b10 * (f(:,3,:) + f(:,this%n,  :)) &
                                            + c10 * (f(:,4,:) + f(:,this%n-1,:))
                    RHS(:,2,:) = RHS(:,2,:) + a10 * (f(:,3,:) + f(:,2,       :)) &
                                            + b10 * (f(:,4,:) + f(:,1,       :)) &
                                            + c10 * (f(:,5,:) + f(:,this%n,  :))
                    !right boundary (n-2:n)
                    RHS(:,this%n-2,:) = RHS(:,this%n-2,:) + a10 * (f(:,this%n-1,:) + f(:,this%n-2,  :)) &
                                                          + b10 * (f(:,this%n,  :) + f(:,this%n-3,  :)) &
                                                          + c10 * (f(:,1,       :) + f(:,this%n-4,  :))
                    RHS(:,this%n-1,:) = RHS(:,this%n-1,:) + a10 * (f(:,this%n,  :) + f(:,this%n-1,  :)) &
                                                          + b10 * (f(:,1,       :) + f(:,this%n-2,  :)) &
                                                          + c10 * (f(:,2,       :) + f(:,this%n-3,  :))
                    RHS(:,this%n,  :) = RHS(:,this%n,  :) + a10 * (f(:,1,       :) + f(:,this%n,    :)) &
                                                          + b10 * (f(:,2,       :) + f(:,this%n-1,  :)) &
                                                          + c10 * (f(:,3,       :) + f(:,this%n-2,  :))

                case ("F2N")
                    !interior    
                    RHS(:,4:this%n-2,:) = RHS(:,4:this%n-2,:) + a10 * (f(:,4:this%n-2,:) + f(:,3:this%n-3,:)) &
                                                              + b10 * (f(:,5:this%n-1,:) + f(:,2:this%n-4,:)) &
                                                              + c10 * (f(:,6:this%n,  :) + f(:,1:this%n-5,:))

                    !left boundary
                    RHS(:,1,:) = RHS(:,1,:) + a10 * (f(:,1,:) + f(:,this%n,  :)) &
                                            + b10 * (f(:,2,:) + f(:,this%n-1,:)) &
                                            + c10 * (f(:,3,:) + f(:,this%n-2,:))
                    RHS(:,2,:) = RHS(:,2,:) + a10 * (f(:,2,:) + f(:,1,       :)) &
                                            + b10 * (f(:,3,:) + f(:,this%n,  :)) &
                                            + c10 * (f(:,4,:) + f(:,this%n-1,:))
                    RHS(:,3,:) = RHS(:,3,:) + a10 * (f(:,3,:) + f(:,2,       :)) &
                                            + b10 * (f(:,4,:) + f(:,1,       :)) &
                                            + c10 * (f(:,5,:) + f(:,this%n,  :))
                    !right boundary (n-1:n)
                    RHS(:,this%n-1,:) = RHS(:,this%n-1,:) + a10 * (f(:,this%n-1,:) + f(:,this%n-2,  :)) &
                                                          + b10 * (f(:,this%n,  :) + f(:,this%n-3,  :)) &
                                                          + c10 * (f(:,1,       :) + f(:,this%n-4,  :))
                    RHS(:,this%n,  :) = RHS(:,this%n,  :) + a10 * (f(:,this%n,  :) + f(:,this%n-1,  :)) &
                                                          + b10 * (f(:,1,       :) + f(:,this%n-2,  :)) &
                                                          + c10 * (f(:,2,       :) + f(:,this%n-3,  :))
            end select

        case (.FALSE.)  !TODO: implement non-periodic BC
            a10 = aI02 * half
            b10 = bI02 * half
            c10 = cI02 * half
            a101 =  one
            a102 =  half
            a104 =  9.0d0/16.0d0
            b104 = -1.0d0/16.0d0
            RHS = 0.0d0
            select case (dir)
                case ("N2F")!TODO: implement better non-periodic BC: currently 2466...6642
                    !interior    
                    RHS(:,3:this%n-3,:) = RHS(:,3:this%n-3,:) + a10 * (f(:,4:this%n-2,:) + f(:,3:this%n-3,:)) &
                                                              + b10 * (f(:,5:this%n-1,:) + f(:,2:this%n-4,:)) &
                                                              + c10 * (f(:,6:this%n  ,:) + f(:,1:this%n-5,:))
                    select case(bc1)
                        !left boundary (1:2)
                        case(1) !symm
                            RHS(:,1,:) = RHS(:,1,:) + a10 * (f(:,2,:) + f(:,1,:)) &
                                                    + b10 * (f(:,3,:) + f(:,2,:)) &
                                                    + c10 * (f(:,4,:) + f(:,3,:))
                            RHS(:,2,:) = RHS(:,2,:) + a10 * (f(:,3,:) + f(:,2,:)) &
                                                    + b10 * (f(:,4,:) + f(:,1,:)) &
                                                    + c10 * (f(:,5,:) + f(:,2,:))
                        case(-1) !anti-symm
                            RHS(:,1,:) = RHS(:,1,:) + a10 * (f(:,2,:) + f(:,1,:)) &
                                                    + b10 * (f(:,3,:) - f(:,2,:)) &
                                                    + c10 * (f(:,4,:) - f(:,3,:))
                            RHS(:,2,:) = RHS(:,2,:) + a10 * (f(:,3,:) + f(:,2,:)) &
                                                    + b10 * (f(:,4,:) + f(:,1,:)) &
                                                    + c10 * (f(:,5,:) - f(:,2,:))
                        case(0)
                            RHS(:,1,:) = RHS(:,1,:) + a102 * (f(:,2,:) + f(:,1,:))   !2nd order
                            RHS(:,2,:) = RHS(:,2,:) + a104 * (f(:,3,:) + f(:,2,:)) & !4th order
                                                    + b104 * (f(:,4,:) + f(:,1,:)) 
                    end select

                    select case(bcn)
                        !right boundary (n-2:n-1)
                        case(1)
                            RHS(:,this%n-2,:) = RHS(:,this%n-2,:) + a10 * (f(:,this%n-1,:) + f(:,this%n-2,:)) &
                                                                  + b10 * (f(:,this%n  ,:) + f(:,this%n-3,:)) &
                                                                  + c10 * (f(:,this%n-1,:) + f(:,this%n-4,:))
                            RHS(:,this%n-1,:) = RHS(:,this%n-1,:) + a10 * (f(:,this%n  ,:) + f(:,this%n-1,:)) &
                                                                  + b10 * (f(:,this%n-1,:) + f(:,this%n-2,:)) &
                                                                  + c10 * (f(:,this%n-2,:) + f(:,this%n-3,:))
                        case(-1)
                            RHS(:,this%n-2,:) = RHS(:,this%n-2,:) + a10 * ( f(:,this%n-1,:) + f(:,this%n-2,:)) &
                                                                  + b10 * ( f(:,this%n  ,:) + f(:,this%n-3,:)) &
                                                                  + c10 * (-f(:,this%n-1,:) + f(:,this%n-4,:))
                            RHS(:,this%n-1,:) = RHS(:,this%n-1,:) + a10 * ( f(:,this%n  ,:) + f(:,this%n-1,:)) &
                                                                  + b10 * (-f(:,this%n-1,:) + f(:,this%n-2,:)) &
                                                                  + c10 * (-f(:,this%n-2,:) + f(:,this%n-3,:))
                        case(0)
                            RHS(:,this%n-2,:) = RHS(:,this%n-2,:) + a104 * (f(:,this%n-1,:) + f(:,this%n-2,  :)) & !4th order
                                                                  + b104 * (f(:,this%n,  :) + f(:,this%n-3,  :))   
                            RHS(:,this%n-1,:) = RHS(:,this%n-1,:) + a102 * (f(:,this%n,  :) + f(:,this%n-1,  :))   !2nd order
                    end select

                case ("F2N")!TODO: implement better non-periodic BC: currently 12466...66421
                    !interior    
                    RHS(:,4:this%n-3,:) = RHS(:,4:this%n-3,:) + a10 * (f(:,4:this%n-3,:) + f(:,3:this%n-4,:)) &
                                                              + b10 * (f(:,5:this%n-2,:) + f(:,2:this%n-5,:)) &
                                                              + c10 * (f(:,6:this%n-1,:) + f(:,1:this%n-6,:))

                    select case(bc1)
                        ! left boundary (1:3)
                        case(1)
                            RHS(:,3,:) = RHS(:,3,:) + a10 * (f(:,3,:) + f(:,2,:)) &
                                                    + b10 * (f(:,4,:) + f(:,1,:)) &
                                                    + c10 * (f(:,5,:) + f(:,1,:))
                            RHS(:,2,:) = RHS(:,2,:) + a10 * (f(:,2,:) + f(:,1,:)) &
                                                    + b10 * (f(:,3,:) + f(:,1,:)) &
                                                    + c10 * (f(:,4,:) + f(:,2,:))
                            RHS(:,1,:) = RHS(:,1,:) + a10 * (f(:,1,:) + f(:,1,:)) &
                                                    + b10 * (f(:,2,:) + f(:,2,:)) &
                                                    + c10 * (f(:,3,:) + f(:,3,:))
                        case(-1)
                            RHS(:,3,:) = RHS(:,3,:) + a10 * (f(:,3,:) + f(:,2,:)) &
                                                    + b10 * (f(:,4,:) + f(:,1,:)) &
                                                    + c10 * (f(:,5,:) - f(:,1,:))
                            RHS(:,2,:) = RHS(:,2,:) + a10 * (f(:,2,:) + f(:,1,:)) &
                                                    + b10 * (f(:,3,:) - f(:,1,:)) &
                                                    + c10 * (f(:,4,:) - f(:,2,:))
                            RHS(:,1,:) = RHS(:,1,:) + a10 * (f(:,1,:) - f(:,1,:)) &
                                                    + b10 * (f(:,2,:) - f(:,2,:)) &
                                                    + c10 * (f(:,3,:) - f(:,3,:))
                        case(0)
                            RHS(:,1,:) = RHS(:,1,:) + a101 * (f(:,1,:)           )
                            RHS(:,2,:) = RHS(:,2,:) + a102 * (f(:,2,:) + f(:,1,:)) 
                            RHS(:,3,:) = RHS(:,3,:) + a104 * (f(:,3,:) + f(:,2,:)) &
                                                    + b104 * (f(:,4,:) + f(:,1,:)) 
                    end select

                    select case(bcn)
                        !right boundary (n-2:n)
                        case(1)
                            RHS(:,this%n-2,:) = RHS(:,this%n-2,:) + a10 * (f(:,this%n-2,:) + f(:,this%n-3,:)) &
                                                                  + b10 * (f(:,this%n-1,:) + f(:,this%n-4,:)) &
                                                                  + c10 * (f(:,this%n-1,:) + f(:,this%n-5,:))
                            RHS(:,this%n-1,:) = RHS(:,this%n-1,:) + a10 * (f(:,this%n-1,:) + f(:,this%n-2,:)) &
                                                                  + b10 * (f(:,this%n-1,:) + f(:,this%n-3,:)) &
                                                                  + c10 * (f(:,this%n-2,:) + f(:,this%n-4,:))
                            RHS(:,this%n,:) = RHS(:,this%n,:)     + a10 * (f(:,this%n-1,:) + f(:,this%n-1,:)) &
                                                                  + b10 * (f(:,this%n-2,:) + f(:,this%n-2,:)) &
                                                                  + c10 * (f(:,this%n-3,:) + f(:,this%n-3,:))
                        case(-1)
                            RHS(:,this%n-2,:) = RHS(:,this%n-2,:) + a10 * ( f(:,this%n-2,:) + f(:,this%n-3,:)) &
                                                                  + b10 * ( f(:,this%n-1,:) + f(:,this%n-4,:)) &
                                                                  + c10 * (-f(:,this%n-1,:) + f(:,this%n-5,:))
                            RHS(:,this%n-1,:) = RHS(:,this%n-1,:) + a10 * ( f(:,this%n-1,:) + f(:,this%n-2,:)) &
                                                                  + b10 * (-f(:,this%n-1,:) + f(:,this%n-3,:)) &
                                                                  + c10 * (-f(:,this%n-2,:) + f(:,this%n-4,:))
                            RHS(:,this%n,:) = RHS(:,this%n,:)     + a10 * (-f(:,this%n-1,:) + f(:,this%n-1,:)) &
                                                                  + b10 * (-f(:,this%n-2,:) + f(:,this%n-2,:)) &
                                                                  + c10 * (-f(:,this%n-3,:) + f(:,this%n-3,:))
                        case(0)
                            RHS(:,this%n,:)   = RHS(:,this%n  ,:) + a101 * (f(:,this%n-1,:)                  )
                            RHS(:,this%n-1,:) = RHS(:,this%n-1,:) + a102 * (f(:,this%n-1,:) + f(:,this%n-2,:)) 
                            RHS(:,this%n-2,:) = RHS(:,this%n-2,:) + a104 * (f(:,this%n-2,:) + f(:,this%n-3,:)) &
                                                                  + b104 * (f(:,this%n-1,:) + f(:,this%n-4,:)) 
                    end select

            end select
        end select
    
    end subroutine

    !pure subroutine ComputeZIRHS(this, f, RHS, n1, n2, bc1, bcn)
    subroutine ComputeZIRHS(this, f, RHS, dir, n1, n2, bc1, bcn)
    
        class( ei06 ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in) :: f
        real(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10, b10, c10, a101, a102, a104, b104
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4
        real(rkind) :: a_np_3
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1
        integer :: j,i

        select case (this%periodic)
        case (.TRUE.)
            a10 = aI02 * half
            b10 = bI02 * half
            c10 = cI02 * half
            RHS = 0.0d0
            select case (dir)
                case ("N2F")
                    !interior    
                    RHS(:,:,3:this%n-3) = RHS(:,:,3:this%n-3) + a10 * (f(:,:,4:this%n-2) + f(:,:,3:this%n-3)) &
                                                              + b10 * (f(:,:,5:this%n-1) + f(:,:,2:this%n-4)) &
                                                              + c10 * (f(:,:,6:this%n  ) + f(:,:,1:this%n-5))

                    !left boundary
                    RHS(:,:,1) = RHS(:,:,1) + a10 * (f(:,:,2) + f(:,:,1       )) &
                                            + b10 * (f(:,:,3) + f(:,:,this%n  )) &
                                            + c10 * (f(:,:,4) + f(:,:,this%n-1))
                    RHS(:,:,2) = RHS(:,:,2) + a10 * (f(:,:,3) + f(:,:,2       )) &
                                            + b10 * (f(:,:,4) + f(:,:,1       )) &
                                            + c10 * (f(:,:,5) + f(:,:,this%n ))
                    !right boundary (n-2:n)
                    RHS(:,:,this%n-2) = RHS(:,:,this%n-2) + a10 * (f(:,:,this%n-1) + f(:,:,this%n-2)) &
                                                          + b10 * (f(:,:,this%n  ) + f(:,:,this%n-3)) &
                                                          + c10 * (f(:,:,1       ) + f(:,:,this%n-4))
                    RHS(:,:,this%n-1) = RHS(:,:,this%n-1) + a10 * (f(:,:,this%n  ) + f(:,:,this%n-1)) &
                                                          + b10 * (f(:,:,1       ) + f(:,:,this%n-2)) &
                                                          + c10 * (f(:,:,2       ) + f(:,:,this%n-3))
                    RHS(:,:,this%n  ) = RHS(:,:,this%n  ) + a10 * (f(:,:,1       ) + f(:,:,this%n  )) &
                                                          + b10 * (f(:,:,2       ) + f(:,:,this%n-1)) &
                                                          + c10 * (f(:,:,3       ) + f(:,:,this%n-2))

                case ("F2N")
                    !interior    
                    RHS(:,:,4:this%n-2) = RHS(:,:,4:this%n-2) + a10 * (f(:,:,4:this%n-2) + f(:,:,3:this%n-3)) &
                                                              + b10 * (f(:,:,5:this%n-1) + f(:,:,2:this%n-4)) &
                                                              + c10 * (f(:,:,6:this%n  ) + f(:,:,1:this%n-5))

                    !left boundary
                    RHS(:,:,1) = RHS(:,:,1) + a10 * (f(:,:,1) + f(:,:,this%n  )) &
                                            + b10 * (f(:,:,2) + f(:,:,this%n-1)) &
                                            + c10 * (f(:,:,3) + f(:,:,this%n-2))
                    RHS(:,:,2) = RHS(:,:,2) + a10 * (f(:,:,2) + f(:,:,1       )) &
                                            + b10 * (f(:,:,3) + f(:,:,this%n  )) &
                                            + c10 * (f(:,:,4) + f(:,:,this%n-1))
                    RHS(:,:,3) = RHS(:,:,3) + a10 * (f(:,:,3) + f(:,:,2       )) &
                                            + b10 * (f(:,:,4) + f(:,:,1       )) &
                                            + c10 * (f(:,:,5) + f(:,:,this%n  ))
                    !right boundary (n-1:n)
                    RHS(:,:,this%n-1) = RHS(:,:,this%n-1) + a10 * (f(:,:,this%n-1) + f(:,:,this%n-2)) &
                                                          + b10 * (f(:,:,this%n  ) + f(:,:,this%n-3)) &
                                                          + c10 * (f(:,:,1       ) + f(:,:,this%n-4))
                    RHS(:,:,this%n  ) = RHS(:,:,this%n  ) + a10 * (f(:,:,this%n  ) + f(:,:,this%n-1)) &
                                                          + b10 * (f(:,:,1       ) + f(:,:,this%n-2)) &
                                                          + c10 * (f(:,:,2       ) + f(:,:,this%n-3))
            end select

        case (.FALSE.)
            a10 = aI02 * half
            b10 = bI02 * half
            c10 = cI02 * half
            a101 =  one
            a102 =  half
            a104 =  9.0d0/16.0d0
            b104 = -1.0d0/16.0d0
            RHS = 0.0d0
            select case (dir)
                case ("N2F")!TODO: implement better non-periodic BC: currently 2466...6642
                    !interior    
                    RHS(:,:,3:this%n-3) = RHS(:,:,3:this%n-3) + a10 * (f(:,:,4:this%n-2) + f(:,:,3:this%n-3)) &
                                                              + b10 * (f(:,:,5:this%n-1) + f(:,:,2:this%n-4)) &
                                                              + c10 * (f(:,:,6:this%n  ) + f(:,:,1:this%n-5))
                    select case(bc1)
                        !left boundary (1:2)
                        case(1) !symm
                            RHS(:,:,1) = RHS(:,:,1) + a10 * (f(:,:,2) + f(:,:,1)) &
                                                    + b10 * (f(:,:,3) + f(:,:,2)) &
                                                    + c10 * (f(:,:,4) + f(:,:,3))
                            RHS(:,:,2) = RHS(:,:,2) + a10 * (f(:,:,3) + f(:,:,2)) &
                                                    + b10 * (f(:,:,4) + f(:,:,1)) &
                                                    + c10 * (f(:,:,5) + f(:,:,2))
                        case(-1) !anti-symm
                            RHS(:,:,1) = RHS(:,:,1) + a10 * (f(:,:,2) + f(:,:,1)) &
                                                    + b10 * (f(:,:,3) - f(:,:,2)) &
                                                    + c10 * (f(:,:,4) - f(:,:,3))
                            RHS(:,:,2) = RHS(:,:,2) + a10 * (f(:,:,3) + f(:,:,2)) &
                                                    + b10 * (f(:,:,4) + f(:,:,1)) &
                                                    + c10 * (f(:,:,5) - f(:,:,2))
                        case(0)
                            RHS(:,:,1) = RHS(:,:,1) + a102 * (f(:,:,2) + f(:,:,1))   !2nd order
                            RHS(:,:,2) = RHS(:,:,2) + a104 * (f(:,:,3) + f(:,:,2)) & !4th order
                                                    + b104 * (f(:,:,4) + f(:,:,1)) 
                    end select

                    select case(bcn)
                        !right boundary (n-2:n-1)
                        case(1)
                            RHS(:,:,this%n-2) = RHS(:,:,this%n-2) + a10 * (f(:,:,this%n-1) + f(:,:,this%n-2)) &
                                                                  + b10 * (f(:,:,this%n  ) + f(:,:,this%n-3)) &
                                                                  + c10 * (f(:,:,this%n-1) + f(:,:,this%n-4))
                            RHS(:,:,this%n-1) = RHS(:,:,this%n-1) + a10 * (f(:,:,this%n  ) + f(:,:,this%n-1)) &
                                                                  + b10 * (f(:,:,this%n-1) + f(:,:,this%n-2)) &
                                                                  + c10 * (f(:,:,this%n-2) + f(:,:,this%n-3))
                        case(-1)
                            RHS(:,:,this%n-2) = RHS(:,:,this%n-2) + a10 * ( f(:,:,this%n-1) + f(:,:,this%n-2)) &
                                                                  + b10 * ( f(:,:,this%n  ) + f(:,:,this%n-3)) &
                                                                  + c10 * (-f(:,:,this%n-1) + f(:,:,this%n-4))
                            RHS(:,:,this%n-1) = RHS(:,:,this%n-1) + a10 * ( f(:,:,this%n  ) + f(:,:,this%n-1)) &
                                                                  + b10 * (-f(:,:,this%n-1) + f(:,:,this%n-2)) &
                                                                  + c10 * (-f(:,:,this%n-2) + f(:,:,this%n-3))
                        case(0)
                            RHS(:,:,this%n-2) = RHS(:,:,this%n-2) + a104 * (f(:,:,this%n-1) + f(:,:,this%n-2)) & !4th order
                                                                  + b104 * (f(:,:,this%n) + f(:,:,this%n-3))   
                            RHS(:,:,this%n-1) = RHS(:,:,this%n-1) + a102 * (f(:,:,this%n) + f(:,:,this%n-1))   !2nd order
                    end select

                case ("F2N")!TODO: implement better non-periodic BC: currently 12466...66421
                    !interior    
                    RHS(:,:,4:this%n-3) = RHS(:,:,4:this%n-3) + a10 * (f(:,:,4:this%n-3) + f(:,:,3:this%n-4)) &
                                                              + b10 * (f(:,:,5:this%n-2) + f(:,:,2:this%n-5)) &
                                                              + c10 * (f(:,:,6:this%n-1) + f(:,:,1:this%n-6))

                    select case(bc1)
                        ! left boundary (1:3)
                        case(1)
                            RHS(:,:,3) = RHS(:,:,3) + a10 * (f(:,:,3) + f(:,:,2)) &
                                                    + b10 * (f(:,:,4) + f(:,:,1)) &
                                                    + c10 * (f(:,:,5) + f(:,:,1))
                            RHS(:,:,2) = RHS(:,:,2) + a10 * (f(:,:,2) + f(:,:,1)) &
                                                    + b10 * (f(:,:,3) + f(:,:,1)) &
                                                    + c10 * (f(:,:,4) + f(:,:,2))
                            RHS(:,:,1) = RHS(:,:,1) + a10 * (f(:,:,1) + f(:,:,1)) &
                                                    + b10 * (f(:,:,2) + f(:,:,2)) &
                                                    + c10 * (f(:,:,3) + f(:,:,3))
                        case(-1)
                            RHS(:,:,3) = RHS(:,:,3) + a10 * (f(:,:,3) + f(:,:,2)) &
                                                    + b10 * (f(:,:,4) + f(:,:,1)) &
                                                    + c10 * (f(:,:,5) - f(:,:,1))
                            RHS(:,:,2) = RHS(:,:,2) + a10 * (f(:,:,2) + f(:,:,1)) &
                                                    + b10 * (f(:,:,3) - f(:,:,1)) &
                                                    + c10 * (f(:,:,4) - f(:,:,2))
                            RHS(:,:,1) = RHS(:,:,1) + a10 * (f(:,:,1) - f(:,:,1)) &
                                                    + b10 * (f(:,:,2) - f(:,:,2)) &
                                                    + c10 * (f(:,:,3) - f(:,:,3))
                        case(0)
                            RHS(:,:,1) = RHS(:,:,1) + a101 * (f(:,:,1)           )
                            RHS(:,:,2) = RHS(:,:,2) + a102 * (f(:,:,2) + f(:,:,1)) 
                            RHS(:,:,3) = RHS(:,:,3) + a104 * (f(:,:,3) + f(:,:,2)) &
                                                    + b104 * (f(:,:,4) + f(:,:,1)) 
                    end select

                    select case(bcn)
                        !right boundary (n-2:n)
                        case(1)
                            RHS(:,:,this%n-2) = RHS(:,:,this%n-2) + a10 * (f(:,:,this%n-2) + f(:,:,this%n-3)) &
                                                                  + b10 * (f(:,:,this%n-1) + f(:,:,this%n-4)) &
                                                                  + c10 * (f(:,:,this%n-1) + f(:,:,this%n-5))
                            RHS(:,:,this%n-1) = RHS(:,:,this%n-1) + a10 * (f(:,:,this%n-1) + f(:,:,this%n-2)) &
                                                                  + b10 * (f(:,:,this%n-1) + f(:,:,this%n-3)) &
                                                                  + c10 * (f(:,:,this%n-2) + f(:,:,this%n-4))
                            RHS(:,:,this%n) = RHS(:,:,this%n)     + a10 * (f(:,:,this%n-1) + f(:,:,this%n-1)) &
                                                                  + b10 * (f(:,:,this%n-2) + f(:,:,this%n-2)) &
                                                                  + c10 * (f(:,:,this%n-3) + f(:,:,this%n-3))
                        case(-1)
                            RHS(:,:,this%n-2) = RHS(:,:,this%n-2) + a10 * ( f(:,:,this%n-2) + f(:,:,this%n-3)) &
                                                                  + b10 * ( f(:,:,this%n-1) + f(:,:,this%n-4)) &
                                                                  + c10 * (-f(:,:,this%n-1) + f(:,:,this%n-5))
                            RHS(:,:,this%n-1) = RHS(:,:,this%n-1) + a10 * ( f(:,:,this%n-1) + f(:,:,this%n-2)) &
                                                                  + b10 * (-f(:,:,this%n-1) + f(:,:,this%n-3)) &
                                                                  + c10 * (-f(:,:,this%n-2) + f(:,:,this%n-4))
                            RHS(:,:,this%n) = RHS(:,:,this%n)     + a10 * (-f(:,:,this%n-1) + f(:,:,this%n-1)) &
                                                                  + b10 * (-f(:,:,this%n-2) + f(:,:,this%n-2)) &
                                                                  + c10 * (-f(:,:,this%n-3) + f(:,:,this%n-3))
                        case(0)
                            RHS(:,:,this%n)   = RHS(:,:,this%n  ) + a101 * (f(:,:,this%n-1)                  )
                            RHS(:,:,this%n-1) = RHS(:,:,this%n-1) + a102 * (f(:,:,this%n-1) + f(:,:,this%n-2)) 
                            RHS(:,:,this%n-2) = RHS(:,:,this%n-2) + a104 * (f(:,:,this%n-2) + f(:,:,this%n-3)) &
                                                                  + b104 * (f(:,:,this%n-1) + f(:,:,this%n-4)) 
                    end select

            end select
        end select
    
    end subroutine
    
    
    subroutine iN2F1(this, fN, fF, na, nb, bc1_, bcn_)
        class( ei06 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(this%n,na,nb), intent(in)  :: fN !fNode
        real(rkind), dimension(this%n,na,nb), intent(out) :: fF !fFace

        if(this%n == 1) then
            fF = zero
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

        
        call this%ComputeXIRHS(fN, fF, "N2F", na, nb, bc1, bcn)

    end subroutine

    subroutine iN2F2(this, fN, fF, na, nb, bc1_, bcn_)
        class( ei06 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(na,this%n,nb), intent(in)  :: fN  !fNode
        real(rkind), dimension(na,this%n,nb), intent(out) :: fF  !fFace

        if(this%n == 1) then
            fF = zero
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

        
        call this%ComputeYIRHS(fN, fF, "N2F", na, nb, bc1, bcn)

    end subroutine

    subroutine iN2F3(this, fN, fF, na, nb, bc1_, bcn_)
        class( ei06 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(na,nb,this%n), intent(in)  :: fN !fNode
        real(rkind), dimension(na,nb,this%n), intent(out) :: fF !fFace

        if(this%n == 1) then
            fF = zero
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


        call this%ComputeZIRHS(fN, fF, "N2F", na, nb, bc1, bcn)

    end subroutine

    subroutine iF2N1(this, fF, fN, na, nb, bc1_, bcn_)
        class( ei06 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(this%n,na,nb), intent(in) :: fF !fFace
        real(rkind), dimension(this%n,na,nb), intent(out)  :: fN !fNode

        if(this%n == 1) then
            fN = zero
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

        
        call this%ComputeXIRHS(fF, fN, "F2N", na, nb, bc1, bcn)

    end subroutine

    subroutine iF2N2(this, fF, fN, na, nb, bc1_, bcn_)
        class( ei06 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(this%n,na,nb), intent(in) :: fF !fFace
        real(rkind), dimension(this%n,na,nb), intent(out)  :: fN !fNode

        if(this%n == 1) then
            fN = zero
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

        
        call this%ComputeYIRHS(fF, fN, "F2N", na, nb, bc1, bcn)

    end subroutine

    subroutine iF2N3(this, fF, fN, na, nb, bc1_, bcn_)
        class( ei06 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(this%n,na,nb), intent(in) :: fF !fFace
        real(rkind), dimension(this%n,na,nb), intent(out)  :: fN !fNode

        if(this%n == 1) then
            fN = zero
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

        
        call this%ComputeZIRHS(fF, fN, "F2N", na, nb, bc1, bcn)

    end subroutine

end module
