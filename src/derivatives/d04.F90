! Routines specific to 4th order Explicit Finite Difference Scheme

module d04stuff

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two
    use exits,           only: GracefulExit
    
    implicit none

    private
    public :: d04, aD04d1, bD04d1, aD04d2, bD04d2
    ! 4th order first derivative explicit centeral difference coefficients
    real(rkind), parameter :: aD04d1     = 2.0_rkind / 3.0_rkind
    real(rkind), parameter :: bD04d1     = -1.0_rkind / 12.0_rkind

    ! 4th order second derivative explicit centeral difference coefficients
    real(rkind), parameter :: aD04d2     = 4.0_rkind / 3.0_rkind
    real(rkind), parameter :: bD04d2     = -1.0_rkind / 12.0_rkind

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic 1st derivative evaluation !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
    ! Interior scheme
    real(rkind), parameter                   :: q_hat           = aD04d1
    real(rkind), parameter                   :: r_hat           = bD04d1

    ! Scheme at node 1 (boundary): explicit one-sided 2nd order
    real(rkind), parameter                   :: p               = -3.0_rkind / 2.0_rkind
    real(rkind), parameter                   :: q               = 4.0_rkind / 2.0_rkind
    real(rkind), parameter                   :: r               = -1.0_rkind / 2.0_rkind
    real(rkind), parameter                   :: s               = 0.0_rkind / 2.0_rkind

    ! Scheme at node 2: explicit central 2nd order
    real(rkind), parameter                   :: q_p             = 1.0_rkind / 2.0_rkind
   
    ! Scheme at node 3: explicit central 4th order
    real(rkind), parameter                   :: q_pp            = 2.0_rkind / 3.0_rkind
    real(rkind), parameter                   :: r_pp            = -1.0_rkind / 12.0_rkind

    ! Scheme at node 4: same as interior
    real(rkind), parameter                   :: q_ppp           = aD04d1
    real(rkind), parameter                   :: r_ppp           = bD04d1

    ! Weights are now trivial, but fluxes do not telescope
    real(rkind), parameter                   :: w1              = 1.0_rkind
    real(rkind), parameter                   :: w2              = 1.0_rkind
    real(rkind), parameter                   :: w3              = 1.0_rkind
    real(rkind), parameter                   :: w4              = 1.0_rkind
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic 2nd derivative evaluation !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   
    ! 1st point: explicit one-sided 2nd order  
    real(rkind), parameter :: b1_aD04d2     = 2.0_rkind
    real(rkind), parameter :: b1_bD04d2     = -5.0_rkind
    real(rkind), parameter :: b1_cD04d2     = 4.0_rkind
    real(rkind), parameter :: b1_dD04d2     = -1.0_rkind
    real(rkind), parameter :: b1_eD04d2     = 0.0_rkind

    ! 2nd point: explicit central 2nd order
    real(rkind), parameter :: b2_aD04d2 = 1.0_rkind

    ! 3rd point: explicit central 4th order
    real(rkind), parameter :: b3_aD04d2  = 4.0_rkind / 3.0_rkind
    real(rkind), parameter :: b3_bD04d2  = -1.0_rkind / 12.0_rkind

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type d04
        
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
        procedure, private :: ComputeXD1RHS
        procedure, private :: ComputeYD1RHS
        procedure, private :: ComputeZD1RHS
        
        procedure, private :: ComputeXD2RHS
        procedure, private :: ComputeYD2RHS
        procedure, private :: ComputeZD2RHS

        procedure :: dd1
        procedure :: dd2
        procedure :: dd3
        
        procedure :: d2d1
        procedure :: d2d2
        procedure :: d2d3

    end type



contains
    
    pure function GetSize(this) result(val)
        class(d04), intent(in) :: this
        integer  :: val 
        val = this%n
    end function

    function init(this, n_, dx_, periodic_, bc1_, bcn_) result(ierr)
   
        class( d04 ), intent(inout) :: this
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

        class( d04 ), intent(inout) :: this

        return

    end subroutine


    pure subroutine ComputeXD1RHS(this, f, RHS, n2, n3, bc1, bcn)
    
        class( d04 ), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10,b10
        integer :: j,k
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4, b_np_4
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1

        select case (this%periodic)
        case (.TRUE.)
            a10 = aD04d1 * this%onebydx
            b10 = bD04d1 * this%onebydx
            do k=1,n3
                do j=1,n2
                    RHS(1,j,k) = a10 * ( f(2,j,k)   - f(this%n  ,j,k) ) &
                           + b10 * ( f(3,j,k)   - f(this%n-1,j,k) )
                    RHS(2,j,k) = a10 * ( f(3,j,k)   - f(1       ,j,k) ) &
                           + b10 * ( f(4,j,k)   - f(this%n  ,j,k) )
                    RHS(3,j,k) = a10 * ( f(4,j,k)   - f(2       ,j,k) ) &
                           + b10 * ( f(5,j,k)   - f(1       ,j,k) )
                    RHS(4:this%n-3,j,k) = a10 * ( f(5:this%n-2,j,k) - f(3:this%n-4,j,k) ) &
                                    + b10 * ( f(6:this%n-1,j,k) - f(2:this%n-5,j,k) )
                    RHS(this%n-2,j,k) = a10 * ( f(this%n-1,j,k) - f(this%n-3,j,k) ) &
                                  + b10 * ( f(this%n  ,j,k) - f(this%n-4,j,k) )
                    RHS(this%n-1,j,k) = a10 * ( f(this%n  ,j,k) - f(this%n-2,j,k) ) &
                                  + b10 * ( f(1       ,j,k) - f(this%n-3,j,k) )
                    RHS(this%n  ,j,k) = a10 * ( f(1       ,j,k) - f(this%n-1,j,k) ) &
                                  + b10 * ( f(2       ,j,k) - f(this%n-2,j,k) )
                end do
            end do
        case (.FALSE.)
            a10    = q_hat * this%onebydx  
            b10    = r_hat * this%onebydx 

            a_np_4 = w4*q_ppp * this%onebydx  
            b_np_4 = w4*r_ppp * this%onebydx 
            
            a_np_3 = w3*q_pp * this%onebydx  
            b_np_3 = w3*r_pp * this%onebydx 

            a_np_2 = w2*q_p * this%onebydx
            
            a_np_1 = w1*( p * this%onebydx)
            b_np_1 = w1*( q * this%onebydx)
            c_np_1 = w1*( r * this%onebydx)
            d_np_1 = w1*( s * this%onebydx)

            do k = 1,n3
                do j = 1,n2
                    select case(bc1)
                    case(0)
                        RHS(1         ,j,k) =   a_np_1* f(1         ,j,k) +  b_np_1*f(2         ,j,k)   &
                                            +   c_np_1* f(3         ,j,k) +  d_np_1*f(4         ,j,k) 
                   
                        RHS(2         ,j,k) =   a_np_2*(f(3         ,j,k) -         f(1         ,j,k))
                        
                        RHS(3         ,j,k) =   a_np_3*(f(4         ,j,k) -         f(2         ,j,k)) &
                                            +   b_np_3*(f(5         ,j,k) -         f(1         ,j,k)) 
                        
                        RHS(4         ,j,k) =   a_np_4*(f(5         ,j,k) -         f(3         ,j,k)) &
                                            +   b_np_4*(f(6         ,j,k) -         f(2         ,j,k))
                    case(1)
                        RHS(1,j,k) =   zero
                   
                        RHS(2,j,k) =   a10   *(f(3,j,k) - f(1,j,k)) &
                                   +   b10   *(f(4,j,k) - f(2,j,k))
                        
                        RHS(3,j,k) =   a10   *(f(4,j,k) - f(2,j,k)) &
                                   +   b10   *(f(5,j,k) - f(1,j,k))
                    
                        RHS(4,j,k) =   a10   *(f(5,j,k) - f(3,j,k)) &
                                   +   b10   *(f(6,j,k) - f(2,j,k))
                    case(-1)
                        RHS(1,j,k) =   a10   *(f(2,j,k) + f(2,j,k)) &
                                   +   b10   *(f(3,j,k) + f(3,j,k))
                        
                        RHS(2,j,k) =   a10   *(f(3,j,k) - f(1,j,k)) &
                                   +   b10   *(f(4,j,k) + f(2,j,k))
                        
                        RHS(3,j,k) =   a10   *(f(4,j,k) - f(2,j,k)) &
                                   +   b10   *(f(5,j,k) - f(1,j,k))
                    
                        RHS(4,j,k) =   a10   *(f(5,j,k) - f(3,j,k)) &
                                   +   b10   *(f(6,j,k) - f(2,j,k))
                    end select
                    
                    RHS(5:this%n-4,j,k) =   a10   *(f(6:this%n-3,j,k) -         f(4:this%n-5,j,k)) &
                                        +   b10   *(f(7:this%n-2,j,k) -         f(3:this%n-6,j,k))
                    
                    select case(bcn)
                    case(0)
                        RHS(this%n-3  ,j,k) =   a_np_4*(f(this%n-2  ,j,k) -         f(this%n-4  ,j,k)) &
                                            +   b_np_4*(f(this%n-1  ,j,k) -         f(this%n-5  ,j,k))
        
                        RHS(this%n-2  ,j,k) =   a_np_3*(f(this%n-1  ,j,k) -         f(this%n-3  ,j,k)) &
                                            +   b_np_3*(f(this%n    ,j,k) -         f(this%n-4  ,j,k)) 
                        
                        RHS(this%n-1  ,j,k) =   a_np_2*(f(this%n    ,j,k) -         f(this%n-2  ,j,k))

                        RHS(this%n    ,j,k) =  -a_np_1* f(this%n    ,j,k) -  b_np_1*f(this%n-1  ,j,k)   &
                                            -   c_np_1* f(this%n-2  ,j,k) -  d_np_1*f(this%n-3  ,j,k)
                    case(1)
                        RHS(this%n-3,j,k) =   a10   *( f(this%n-2,j,k) - f(this%n-4,j,k)) &
                                          +   b10   *( f(this%n-1,j,k) - f(this%n-5,j,k))

                        RHS(this%n-2,j,k) =   a10   *( f(this%n-1,j,k) - f(this%n-3,j,k)) &
                                          +   b10   *( f(this%n  ,j,k) - f(this%n-4,j,k))

                        RHS(this%n-1,j,k) =   a10   *( f(this%n  ,j,k) - f(this%n-2,j,k)) &
                                          +   b10   *( f(this%n-1,j,k) - f(this%n-3,j,k))

                        RHS(this%n  ,j,k) =   zero

                    case(-1)
                        RHS(this%n-3,j,k) =   a10   *( f(this%n-2,j,k) - f(this%n-4,j,k)) &
                                          +   b10   *( f(this%n-1,j,k) - f(this%n-5,j,k))

                        RHS(this%n-2,j,k) =   a10   *( f(this%n-1,j,k) - f(this%n-3,j,k)) &
                                          +   b10   *( f(this%n  ,j,k) - f(this%n-4,j,k))

                        RHS(this%n-1,j,k) =   a10   *( f(this%n  ,j,k) - f(this%n-2,j,k)) &
                                          +   b10   *(-f(this%n-1,j,k) - f(this%n-3,j,k))

                        RHS(this%n  ,j,k) =   a10   *(-f(this%n-1,j,k) - f(this%n-1,j,k)) &
                                          +   b10   *(-f(this%n-2,j,k) - f(this%n-2,j,k))

                    end select
               end do 
            end do 
        end select
    
    end subroutine
    
    pure subroutine ComputeYD1RHS(this, f, RHS, n1, n3, bc1, bcn) 
    
        class( d04 ), intent(in) :: this
        integer, intent(in) :: n1, n3
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10,b10
        integer :: k
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4, b_np_4
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1


        select case (this%periodic)
        case (.TRUE.)
            a10 = aD04d1 * this%onebydx
            b10 = bD04d1 * this%onebydx
            do k=1,n3
                RHS(:,1,k) = a10 * ( f(:,2,k)   - f(:,this%n  ,k) ) &
                       + b10 * ( f(:,3,k)   - f(:,this%n-1,k) )
                RHS(:,2,k) = a10 * ( f(:,3,k)   - f(:,1       ,k) ) &
                       + b10 * ( f(:,4,k)   - f(:,this%n  ,k) )
                RHS(:,3,k) = a10 * ( f(:,4,k)   - f(:,2       ,k) ) &
                       + b10 * ( f(:,5,k)   - f(:,1       ,k) )
                RHS(:,4:this%n-3,k) = a10 * ( f(:,5:this%n-2,k) - f(:,3:this%n-4,k) ) &
                                + b10 * ( f(:,6:this%n-1,k) - f(:,2:this%n-5,k) )
                RHS(:,this%n-2,k) = a10 * ( f(:,this%n-1,k) - f(:,this%n-3,k) ) &
                              + b10 * ( f(:,this%n  ,k) - f(:,this%n-4,k) )
                RHS(:,this%n-1,k) = a10 * ( f(:,this%n  ,k) - f(:,this%n-2,k) ) &
                              + b10 * ( f(:,1       ,k) - f(:,this%n-3,k) )
                RHS(:,this%n  ,k) = a10 * ( f(:,1       ,k) - f(:,this%n-1,k) ) &
                              + b10 * ( f(:,2       ,k) - f(:,this%n-2,k) )
            end do
        case (.FALSE.)
            a10    = q_hat * this%onebydx  
            b10    = r_hat * this%onebydx 

            a_np_4 = w4*q_ppp * this%onebydx  
            b_np_4 = w4*r_ppp * this%onebydx 
            
            a_np_3 = w3*q_pp * this%onebydx  
            b_np_3 = w3*r_pp * this%onebydx 

            a_np_2 = w2*q_p * this%onebydx
            
            a_np_1 = w1*( p * this%onebydx)
            b_np_1 = w1*( q * this%onebydx)
            c_np_1 = w1*( r * this%onebydx)
            d_np_1 = w1*( s * this%onebydx)

            do k = 1,n3

                select case(bc1)
                case(0)    
                    RHS(:,1         ,k) =   a_np_1* f(:,1         ,k) +  b_np_1*f(:,2         ,k)   &
                                        +   c_np_1* f(:,3         ,k) +  d_np_1*f(:,4         ,k) 
                    
                    RHS(:,2         ,k) =   a_np_2*(f(:,3         ,k) -         f(:,1         ,k))
                    
                    RHS(:,3         ,k) =   a_np_3*(f(:,4         ,k) -         f(:,2         ,k)) &
                                        +   b_np_3*(f(:,5         ,k) -         f(:,1         ,k)) 
                    
                    RHS(:,4         ,k) =   a_np_4*(f(:,5         ,k) -         f(:,3         ,k)) &
                                        +   b_np_4*(f(:,6         ,k) -         f(:,2         ,k))
                case(1)
                    RHS(:,1,k) =   zero
                   
                    RHS(:,2,k) =   a10   *(f(:,3,k) - f(:,1,k)) &
                               +   b10   *(f(:,4,k) - f(:,2,k))
                    
                    RHS(:,3,k) =   a10   *(f(:,4,k) - f(:,2,k)) &
                               +   b10   *(f(:,5,k) - f(:,1,k))
                    
                    RHS(:,4,k) =   a10   *(f(:,5,k) - f(:,3,k)) &
                               +   b10   *(f(:,6,k) - f(:,2,k))
                case(-1)
                    RHS(:,1,k) =   a10   *(f(:,2,k) + f(:,2,k)) &
                               +   b10   *(f(:,3,k) + f(:,3,k))
                    
                    RHS(:,2,k) =   a10   *(f(:,3,k) - f(:,1,k)) &
                               +   b10   *(f(:,4,k) + f(:,2,k))
                    
                    RHS(:,3,k) =   a10   *(f(:,4,k) - f(:,2,k)) &
                               +   b10   *(f(:,5,k) - f(:,1,k))
                    
                    RHS(:,4,k) =   a10   *(f(:,5,k) - f(:,3,k)) &
                               +   b10   *(f(:,6,k) - f(:,2,k))
                end select
                
                RHS(:,5:this%n-4,k) =   a10   *(f(:,6:this%n-3,k) -         f(:,4:this%n-5,k)) &
                                    +   b10   *(f(:,7:this%n-2,k) -         f(:,3:this%n-6,k))
                
                select case(bcn)
                case(0)    
                    RHS(:,this%n-3  ,k) =   a_np_4*(f(:,this%n-2  ,k) -         f(:,this%n-4  ,k)) &
                                        +   b_np_4*(f(:,this%n-1  ,k) -         f(:,this%n-5  ,k))
        
                    RHS(:,this%n-2  ,k) =   a_np_3*(f(:,this%n-1  ,k) -         f(:,this%n-3  ,k)) &
                                        +   b_np_3*(f(:,this%n    ,k) -         f(:,this%n-4  ,k)) 
                    
                    RHS(:,this%n-1  ,k) =   a_np_2*(f(:,this%n    ,k) -         f(:,this%n-2  ,k))

                    RHS(:,this%n    ,k) =  -a_np_1* f(:,this%n    ,k) -  b_np_1*f(:,this%n-1  ,k)   &
                                        -   c_np_1* f(:,this%n-2  ,k) -  d_np_1*f(:,this%n-3  ,k)
                case(1)
                    RHS(:,this%n-3,k) =   a10   *( f(:,this%n-2,k) - f(:,this%n-4,k)) &
                                      +   b10   *( f(:,this%n-1,k) - f(:,this%n-5,k))

                    RHS(:,this%n-2,k) =   a10   *( f(:,this%n-1,k) - f(:,this%n-3,k)) &
                                      +   b10   *( f(:,this%n  ,k) - f(:,this%n-4,k))

                    RHS(:,this%n-1,k) =   a10   *( f(:,this%n  ,k) - f(:,this%n-2,k)) &
                                      +   b10   *( f(:,this%n-1,k) - f(:,this%n-3,k))

                    RHS(:,this%n  ,k) =   zero
                case(-1)
                    RHS(:,this%n-3,k) =   a10   *( f(:,this%n-2,k) - f(:,this%n-4,k)) &
                                      +   b10   *( f(:,this%n-1,k) - f(:,this%n-5,k))

                    RHS(:,this%n-2,k) =   a10   *( f(:,this%n-1,k) - f(:,this%n-3,k)) &
                                      +   b10   *( f(:,this%n  ,k) - f(:,this%n-4,k))

                    RHS(:,this%n-1,k) =   a10   *( f(:,this%n  ,k) - f(:,this%n-2,k)) &
                                      +   b10   *(-f(:,this%n-1,k) - f(:,this%n-3,k))

                    RHS(:,this%n  ,k) =   a10   *(-f(:,this%n-1,k) - f(:,this%n-1,k)) &
                                      +   b10   *(-f(:,this%n-2,k) - f(:,this%n-2,k))
                end select
            end do 
        end select
    
    end subroutine

    pure subroutine ComputeZD1RHS(this, f, RHS, n1, n2, bc1, bcn)
    
        class( d04 ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in) :: f
        real(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10,b10
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4, b_np_4
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1
        integer :: j


        select case (this%periodic)
        case (.TRUE.)
            a10 = aD04d1 * this%onebydx
            b10 = bD04d1 * this%onebydx
            do j=1,n2
                RHS(:,j,1) = a10 * ( f(:,j,2)   - f(:,j,this%n  ) ) &
                                   + b10 * ( f(:,j,3)   - f(:,j,this%n-1) )
                RHS(:,j,2) = a10 * ( f(:,j,3)   - f(:,j,1       ) ) &
                                   + b10 * ( f(:,j,4)   - f(:,j,this%n  ) )
                RHS(:,j,3) = a10 * ( f(:,j,4)   - f(:,j,2       ) ) &
                                   + b10 * ( f(:,j,5)   - f(:,j,1       ) )
                RHS(:,j,4:this%n-3) = a10 * ( f(:,j,5:this%n-2) - f(:,j,3:this%n-4) ) &
                                            + b10 * ( f(:,j,6:this%n-1) - f(:,j,2:this%n-5) )
                RHS(:,j,this%n-2) = a10 * ( f(:,j,this%n-1) - f(:,j,this%n-3) ) &
                                          + b10 * ( f(:,j,this%n  ) - f(:,j,this%n-4) )
                RHS(:,j,this%n-1) = a10 * ( f(:,j,this%n  ) - f(:,j,this%n-2) ) &
                                          + b10 * ( f(:,j,1       ) - f(:,j,this%n-3) )
                RHS(:,j,this%n  ) = a10 * ( f(:,j,1       ) - f(:,j,this%n-1) ) &
                                                  + b10 * ( f(:,j,2       ) - f(:,j,this%n-2) )
            end do
        case (.FALSE.)
            a10    = q_hat * this%onebydx  
            b10    = r_hat * this%onebydx 

            a_np_4 = w4*q_ppp * this%onebydx  
            b_np_4 = w4*r_ppp * this%onebydx 
            
            a_np_3 = w3*q_pp * this%onebydx  
            b_np_3 = w3*r_pp * this%onebydx 

            a_np_2 = w2*q_p * this%onebydx
            
            a_np_1 = w1*( p * this%onebydx)
            b_np_1 = w1*( q * this%onebydx)
            c_np_1 = w1*( r * this%onebydx)
            d_np_1 = w1*( s * this%onebydx)

                    
            select case(bc1)
            case(0)    
                RHS(:,:,1         ) =   a_np_1* f(:,:,1         ) +  b_np_1*f(:,:,2         )   &
                                    +   c_np_1* f(:,:,3         ) +  d_np_1*f(:,:,4         ) 
                
                RHS(:,:,2         ) =   a_np_2*(f(:,:,3         ) -         f(:,:,1         ))
                
                RHS(:,:,3         ) =   a_np_3*(f(:,:,4         ) -         f(:,:,2         )) &
                                    +   b_np_3*(f(:,:,5         ) -         f(:,:,1         )) 
                
                RHS(:,:,4         ) =   a_np_4*(f(:,:,5         ) -         f(:,:,3         )) &
                                    +   b_np_4*(f(:,:,6         ) -         f(:,:,2         ))
            case(1)
                RHS(:,:,1) =   zero
                
                RHS(:,:,2) =   a10   *(f(:,:,3) - f(:,:,1)) &
                           +   b10   *(f(:,:,4) - f(:,:,2))
                
                RHS(:,:,3) =   a10   *(f(:,:,4) - f(:,:,2)) &
                           +   b10   *(f(:,:,5) - f(:,:,1))
                
                RHS(:,:,4) =   a10   *(f(:,:,5) - f(:,:,3)) &
                           +   b10   *(f(:,:,6) - f(:,:,2))
            case(-1)
                RHS(:,:,1) =   a10   *(f(:,:,2) + f(:,:,2)) &
                           +   b10   *(f(:,:,3) + f(:,:,3))
                
                RHS(:,:,2) =   a10   *(f(:,:,3) - f(:,:,1)) &
                           +   b10   *(f(:,:,4) + f(:,:,2))
                
                RHS(:,:,3) =   a10   *(f(:,:,4) - f(:,:,2)) &
                           +   b10   *(f(:,:,5) - f(:,:,1))
                
                RHS(:,:,4) =   a10   *(f(:,:,5) - f(:,:,3)) &
                           +   b10   *(f(:,:,6) - f(:,:,2))
            end select
            
            RHS(:,:,5:this%n-4) =   a10   *(f(:,:,6:this%n-3) -         f(:,:,4:this%n-5)) &
                                +   b10   *(f(:,:,7:this%n-2) -         f(:,:,3:this%n-6))
            
            select case(bcn)
            case(0)    
                RHS(:,:,this%n-3  ) =   a_np_4*(f(:,:,this%n-2  ) -         f(:,:,this%n-4  )) &
                                    +   b_np_4*(f(:,:,this%n-1  ) -         f(:,:,this%n-5  ))
        
                RHS(:,:,this%n-2  ) =   a_np_3*(f(:,:,this%n-1  ) -         f(:,:,this%n-3  )) &
                                    +   b_np_3*(f(:,:,this%n    ) -         f(:,:,this%n-4  )) 
                
                RHS(:,:,this%n-1  ) =   a_np_2*(f(:,:,this%n    ) -         f(:,:,this%n-2  ))

                RHS(:,:,this%n    ) =  -a_np_1* f(:,:,this%n    ) -  b_np_1*f(:,:,this%n-1  )   &
                                    -   c_np_1* f(:,:,this%n-2  ) -  d_np_1*f(:,:,this%n-3  )
            case(1)
                RHS(:,:,this%n-3) =   a10   *( f(:,:,this%n-2) - f(:,:,this%n-4)) &
                                  +   b10   *( f(:,:,this%n-1) - f(:,:,this%n-5))

                RHS(:,:,this%n-2) =   a10   *( f(:,:,this%n-1) - f(:,:,this%n-3)) &
                                  +   b10   *( f(:,:,this%n  ) - f(:,:,this%n-4))

                RHS(:,:,this%n-1) =   a10   *( f(:,:,this%n  ) - f(:,:,this%n-2)) &
                                  +   b10   *( f(:,:,this%n-1) - f(:,:,this%n-3))

                RHS(:,:,this%n  ) =   zero
            case(-1)
                RHS(:,:,this%n-3) =   a10   *( f(:,:,this%n-2) - f(:,:,this%n-4)) &
                                  +   b10   *( f(:,:,this%n-1) - f(:,:,this%n-5))

                RHS(:,:,this%n-2) =   a10   *( f(:,:,this%n-1) - f(:,:,this%n-3)) &
                                  +   b10   *( f(:,:,this%n  ) - f(:,:,this%n-4))

                RHS(:,:,this%n-1) =   a10   *( f(:,:,this%n  ) - f(:,:,this%n-2)) &
                                  +   b10   *(-f(:,:,this%n-1) - f(:,:,this%n-3))

                RHS(:,:,this%n  ) =   a10   *(-f(:,:,this%n-1) - f(:,:,this%n-1)) &
                                  +   b10   *(-f(:,:,this%n-2) - f(:,:,this%n-2))
            end select
            
        end select
    
    end subroutine
    
   pure subroutine ComputeXD2RHS(this, f, RHS, n2, n3, bc1, bcn)
        class( d04 ), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10, b10
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1, e_np_1
        integer :: j,k 
    
        select case (this%periodic)
        case (.TRUE.)
            a10 = aD04d2 * this%onebydx2
            b10 = bD04d2 * this%onebydx2

            do k = 1,n3
                do j = 1,n2
                        RHS(1         ,j,k) = a10 * ( f(2          ,j,k)   - two*f(1            ,j,k) + f(this%n    ,j,k)) &
                                            + b10 * ( f(3          ,j,k)   - two*f(1            ,j,k) + f(this%n-1  ,j,k))
                        
                        RHS(2         ,j,k) = a10 * ( f(3          ,j,k)   - two*f(2            ,j,k) + f(1         ,j,k)) &
                                            + b10 * ( f(4          ,j,k)   - two*f(2            ,j,k) + f(this%n    ,j,k))
                        
                        RHS(3         ,j,k) = a10 * ( f(4          ,j,k)   - two*f(3            ,j,k) + f(2         ,j,k)) &
                                            + b10 * ( f(5          ,j,k)   - two*f(3            ,j,k) + f(1         ,j,k))
                    
                    RHS(4:this%n-3,j,k) = a10 * ( f(5:this%n-2 ,j,k)   - two*f(4:this%n-3   ,j,k) + f(3:this%n-4,j,k)) &
                                        + b10 * ( f(6:this%n-1 ,j,k)   - two*f(4:this%n-3   ,j,k) + f(2:this%n-5,j,k))
                   
                        RHS(this%n-2  ,j,k) = a10 * ( f(this%n-1   ,j,k)   - two*f(this%n-2     ,j,k) + f(this%n-3  ,j,k)) &
                                            + b10 * ( f(this%n     ,j,k)   - two*f(this%n-2     ,j,k) + f(this%n-4  ,j,k))
                        
                        RHS(this%n-1  ,j,k) = a10 * ( f(this%n     ,j,k)   - two*f(this%n-1     ,j,k) + f(this%n-2  ,j,k)) &
                                            + b10 * ( f(1          ,j,k)   - two*f(this%n-1     ,j,k) + f(this%n-3  ,j,k))
                        
                        RHS(this%n    ,j,k) = a10 * ( f(1          ,j,k)   - two*f(this%n       ,j,k) + f(this%n-1  ,j,k)) &
                                            + b10 * ( f(2          ,j,k)   - two*f(this%n       ,j,k) + f(this%n-2  ,j,k))
                end do 
            end do 

        case (.FALSE.)
            a10 = aD04d2 * this%onebydx2
            b10 = bD04d2 * this%onebydx2
            
            a_np_3 = b3_aD04d2 * this%onebydx2 
            b_np_3 = b3_bD04d2 * this%onebydx2 
       
            a_np_2 = b2_aD04d2 * this%onebydx2

            a_np_1 = b1_aD04d2 * this%onebydx2
            b_np_1 = b1_bD04d2 * this%onebydx2
            c_np_1 = b1_cD04d2 * this%onebydx2
            d_np_1 = b1_dD04d2 * this%onebydx2
            e_np_1 = b1_eD04d2 * this%onebydx2


            do k = 1,n3
                do j = 1,n2
                    select case(bc1)
                    case(0)
                        RHS(1         ,j,k) = a_np_1*f(1,j,k) + b_np_1*f(2,j,k) + &
                                            & c_np_1*f(3,j,k) + d_np_1*f(4,j,k) + &
                                            & e_np_1*f(5,j,k)
                        
                        RHS(2         ,j,k) = a_np_2 * ( f(3          ,j,k)   - two*f(2            ,j,k) + f(1         ,j,k)) 
                        
                        RHS(3         ,j,k) = a_np_3 * ( f(4          ,j,k)   - two*f(3            ,j,k) + f(2         ,j,k)) &
                                            + b_np_3 * ( f(5          ,j,k)   - two*f(3            ,j,k) + f(1         ,j,k)) 
                    case(1)
                        RHS(1,j,k) = a10 * ( f(2,j,k)   - two*f(1,j,k) + f(2,j,k)) &
                                   + b10 * ( f(3,j,k)   - two*f(1,j,k) + f(3,j,k))
                        
                        RHS(2,j,k) = a10 * ( f(3,j,k)   - two*f(2,j,k) + f(1,j,k)) &
                                   + b10 * ( f(4,j,k)   - two*f(2,j,k) + f(2,j,k))
                        
                        RHS(3,j,k) = a10 * ( f(4,j,k)   - two*f(3,j,k) + f(2,j,k)) &
                                   + b10 * ( f(5,j,k)   - two*f(3,j,k) + f(1,j,k))
                    case(-1)
                        RHS(1,j,k) = a10 * ( f(2,j,k)   - two*f(1,j,k) - f(2,j,k)) &
                                   + b10 * ( f(3,j,k)   - two*f(1,j,k) - f(3,j,k))
                        
                        RHS(2,j,k) = a10 * ( f(3,j,k)   - two*f(2,j,k) + f(1,j,k)) &
                                   + b10 * ( f(4,j,k)   - two*f(2,j,k) - f(2,j,k))
                        
                        RHS(3,j,k) = a10 * ( f(4,j,k)   - two*f(3,j,k) + f(2,j,k)) &
                                   + b10 * ( f(5,j,k)   - two*f(3,j,k) + f(1,j,k))
                    end select
                    
                    RHS(4:this%n-3,j,k) = a10 * ( f(5:this%n-2 ,j,k)   - two*f(4:this%n-3   ,j,k) + f(3:this%n-4,j,k)) &
                                        + b10 * ( f(6:this%n-1 ,j,k)   - two*f(4:this%n-3   ,j,k) + f(2:this%n-5,j,k))
                    
                    select case(bcn)
                    case(0) 
                        RHS(this%n-2  ,j,k) = a_np_3 * ( f(this%n-1   ,j,k)   - two*f(this%n-2     ,j,k) + f(this%n-3  ,j,k)) &
                                            + b_np_3 * ( f(this%n     ,j,k)   - two*f(this%n-2     ,j,k) + f(this%n-4  ,j,k)) 
                        
                        RHS(this%n-1  ,j,k) = a_np_2 * ( f(this%n     ,j,k)   - two*f(this%n-1     ,j,k) + f(this%n-2  ,j,k)) 
                        
                        RHS(this%n    ,j,k) =  a_np_1*f(this%n  ,j,k) + b_np_1*f(this%n-1,j,k) + &
                                            &  c_np_1*f(this%n-2,j,k) + d_np_1*f(this%n-3,j,k) + &
                                            &  e_np_1*f(this%n-4,j,k)
                    case(1)
                        RHS(this%n-2,j,k) = a10 * ( f(this%n-1,j,k)   - two*f(this%n-2,j,k) + f(this%n-3,j,k)) &
                                          + b10 * ( f(this%n  ,j,k)   - two*f(this%n-2,j,k) + f(this%n-4,j,k))
                        
                        RHS(this%n-1,j,k) = a10 * ( f(this%n  ,j,k)   - two*f(this%n-1,j,k) + f(this%n-2,j,k)) &
                                          + b10 * ( f(this%n-1,j,k)   - two*f(this%n-1,j,k) + f(this%n-3,j,k))
                        
                        RHS(this%n  ,j,k) = a10 * ( f(this%n-1,j,k)   - two*f(this%n  ,j,k) + f(this%n-1,j,k)) &
                                          + b10 * ( f(this%n-2,j,k)   - two*f(this%n  ,j,k) + f(this%n-2,j,k))
                    case(-1)
                        RHS(this%n-2,j,k) = a10 * ( f(this%n-1,j,k)   - two*f(this%n-2,j,k) + f(this%n-3,j,k)) &
                                          + b10 * ( f(this%n  ,j,k)   - two*f(this%n-2,j,k) + f(this%n-4,j,k))
                        
                        RHS(this%n-1,j,k) = a10 * ( f(this%n  ,j,k)   - two*f(this%n-1,j,k) + f(this%n-2,j,k)) &
                                          + b10 * (-f(this%n-1,j,k)   - two*f(this%n-1,j,k) + f(this%n-3,j,k))
                        
                        RHS(this%n  ,j,k) = a10 * (-f(this%n-1,j,k)   - two*f(this%n  ,j,k) + f(this%n-1,j,k)) &
                                          + b10 * (-f(this%n-2,j,k)   - two*f(this%n  ,j,k) + f(this%n-2,j,k))

                    end select

                end do 
            end do 
        end select
   
    end subroutine  


   pure subroutine ComputeYD2RHS(this, f, RHS, n2, n3, bc1, bcn)
        class( d04 ), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(n2,this%n,n3), intent(in) :: f
        real(rkind), dimension(n2,this%n,n3), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10, b10
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1, e_np_1
    
        integer :: k
    
        select case (this%periodic)
        case (.TRUE.)
            a10 = aD04d2 * this%onebydx2
            b10 = bD04d2 * this%onebydx2

            do k = 1,n3
                RHS(:,1         ,k) = a10 * ( f(:,2          ,k)   - two*f(:,1            ,k) + f(:,this%n    ,k)) &
                                    + b10 * ( f(:,3          ,k)   - two*f(:,1            ,k) + f(:,this%n-1  ,k))
                
                RHS(:,2         ,k) = a10 * ( f(:,3          ,k)   - two*f(:,2            ,k) + f(:,1         ,k)) &
                                    + b10 * ( f(:,4          ,k)   - two*f(:,2            ,k) + f(:,this%n    ,k))
                
                RHS(:,3         ,k) = a10 * ( f(:,4          ,k)   - two*f(:,3            ,k) + f(:,2         ,k)) &
                                    + b10 * ( f(:,5          ,k)   - two*f(:,3            ,k) + f(:,1         ,k))
                
                RHS(:,4:this%n-3,k) = a10 * ( f(:,5:this%n-2 ,k)   - two*f(:,4:this%n-3   ,k) + f(:,3:this%n-4,k)) &
                                    + b10 * ( f(:,6:this%n-1 ,k)   - two*f(:,4:this%n-3   ,k) + f(:,2:this%n-5,k))
                
                RHS(:,this%n-2  ,k) = a10 * ( f(:,this%n-1   ,k)   - two*f(:,this%n-2     ,k) + f(:,this%n-3  ,k)) &
                                    + b10 * ( f(:,this%n     ,k)   - two*f(:,this%n-2     ,k) + f(:,this%n-4  ,k))
                
                RHS(:,this%n-1  ,k) = a10 * ( f(:,this%n     ,k)   - two*f(:,this%n-1     ,k) + f(:,this%n-2  ,k)) &
                                    + b10 * ( f(:,1          ,k)   - two*f(:,this%n-1     ,k) + f(:,this%n-3  ,k))
                
                RHS(:,this%n    ,k) = a10 * ( f(:,1          ,k)   - two*f(:,this%n       ,k) + f(:,this%n-1  ,k)) &
                                    + b10 * ( f(:,2          ,k)   - two*f(:,this%n       ,k) + f(:,this%n-2  ,k))
            end do 

        case (.FALSE.)
            a10 = aD04d2 * this%onebydx2
            b10 = bD04d2 * this%onebydx2
            
            a_np_3 = b3_aD04d2 * this%onebydx2 
            b_np_3 = b3_bD04d2 * this%onebydx2 
       
            a_np_2 = b2_aD04d2 * this%onebydx2

            a_np_1 = b1_aD04d2 * this%onebydx2
            b_np_1 = b1_bD04d2 * this%onebydx2
            c_np_1 = b1_cD04d2 * this%onebydx2
            d_np_1 = b1_dD04d2 * this%onebydx2
            e_np_1 = b1_eD04d2 * this%onebydx2


            do k = 1,n3
                select case(bc1)
                case(0)
                    RHS(:,1         ,k) = a_np_1*f(:,1,k) + b_np_1*f(:,2,k) + &
                                        & c_np_1*f(:,3,k) + d_np_1*f(:,4,k) + &
                                        & e_np_1*f(:,5,k)
                    
                    RHS(:,2         ,k) = a_np_2 * ( f(:,3          ,k)   - two*f(:,2            ,k) + f(:,1         ,k)) 
                    
                    RHS(:,3         ,k) = a_np_3 * ( f(:,4          ,k)   - two*f(:,3            ,k) + f(:,2         ,k)) &
                                        + b_np_3 * ( f(:,5          ,k)   - two*f(:,3            ,k) + f(:,1         ,k)) 
                case(1)
                    RHS(:,1,k) = a10 * ( f(:,2,k)   - two*f(:,1,k) + f(:,2,k)) &
                               + b10 * ( f(:,3,k)   - two*f(:,1,k) + f(:,3,k))
                    
                    RHS(:,2,k) = a10 * ( f(:,3,k)   - two*f(:,2,k) + f(:,1,k)) &
                               + b10 * ( f(:,4,k)   - two*f(:,2,k) + f(:,2,k))
                    
                    RHS(:,3,k) = a10 * ( f(:,4,k)   - two*f(:,3,k) + f(:,2,k)) &
                               + b10 * ( f(:,5,k)   - two*f(:,3,k) + f(:,1,k))
                case(-1)
                    RHS(:,1,k) = a10 * ( f(:,2,k)   - two*f(:,1,k) - f(:,2,k)) &
                               + b10 * ( f(:,3,k)   - two*f(:,1,k) - f(:,3,k))
                    
                    RHS(:,2,k) = a10 * ( f(:,3,k)   - two*f(:,2,k) + f(:,1,k)) &
                               + b10 * ( f(:,4,k)   - two*f(:,2,k) - f(:,2,k))
                    
                    RHS(:,3,k) = a10 * ( f(:,4,k)   - two*f(:,3,k) + f(:,2,k)) &
                               + b10 * ( f(:,5,k)   - two*f(:,3,k) + f(:,1,k))
                end select
                    
                    RHS(:,4:this%n-3,k) = a10 * ( f(:,5:this%n-2 ,k)   - two*f(:,4:this%n-3   ,k) + f(:,3:this%n-4,k)) &
                                        + b10 * ( f(:,6:this%n-1 ,k)   - two*f(:,4:this%n-3   ,k) + f(:,2:this%n-5,k))
                    
                select case(bcn)
                case(0) 
                    RHS(:,this%n-2  ,k) = a_np_3 * ( f(:,this%n-1   ,k)   - two*f(:,this%n-2     ,k) + f(:,this%n-3  ,k)) &
                                        + b_np_3 * ( f(:,this%n     ,k)   - two*f(:,this%n-2     ,k) + f(:,this%n-4  ,k)) 
                    
                    RHS(:,this%n-1  ,k) = a_np_2 * ( f(:,this%n     ,k)   - two*f(:,this%n-1     ,k) + f(:,this%n-2  ,k)) 
                    
                    RHS(:,this%n    ,k) =  a_np_1*f(:,this%n  ,k) + b_np_1*f(:,this%n-1,k) + &
                                        &  c_np_1*f(:,this%n-2,k) + d_np_1*f(:,this%n-3,k) + &
                                        &  e_np_1*f(:,this%n-4,k)
                case(1)
                    RHS(:,this%n-2,k) = a10 * ( f(:,this%n-1,k)   - two*f(:,this%n-2,k) + f(:,this%n-3,k)) &
                                      + b10 * ( f(:,this%n  ,k)   - two*f(:,this%n-2,k) + f(:,this%n-4,k))
                    
                    RHS(:,this%n-1,k) = a10 * ( f(:,this%n  ,k)   - two*f(:,this%n-1,k) + f(:,this%n-2,k)) &
                                      + b10 * ( f(:,this%n-1,k)   - two*f(:,this%n-1,k) + f(:,this%n-3,k))
                    
                    RHS(:,this%n  ,k) = a10 * ( f(:,this%n-1,k)   - two*f(:,this%n  ,k) + f(:,this%n-1,k)) &
                                      + b10 * ( f(:,this%n-2,k)   - two*f(:,this%n  ,k) + f(:,this%n-2,k))
                case(-1)
                    RHS(:,this%n-2,k) = a10 * ( f(:,this%n-1,k)   - two*f(:,this%n-2,k) + f(:,this%n-3,k)) &
                                      + b10 * ( f(:,this%n  ,k)   - two*f(:,this%n-2,k) + f(:,this%n-4,k))
                    
                    RHS(:,this%n-1,k) = a10 * ( f(:,this%n  ,k)   - two*f(:,this%n-1,k) + f(:,this%n-2,k)) &
                                      + b10 * (-f(:,this%n-1,k)   - two*f(:,this%n-1,k) + f(:,this%n-3,k))
                    
                    RHS(:,this%n  ,k) = a10 * (-f(:,this%n-1,k)   - two*f(:,this%n  ,k) + f(:,this%n-1,k)) &
                                      + b10 * (-f(:,this%n-2,k)   - two*f(:,this%n  ,k) + f(:,this%n-2,k))
                end select

            end do 
        
        end select
   
    end subroutine  

    pure subroutine ComputeZD2RHS(this, f, RHS, n2, n3, bc1, bcn)
        class( d04 ), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(n2,n3,this%n), intent(in) :: f
        real(rkind), dimension(n2,n3,this%n), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10, b10
        real(rkind) :: a_np_3, b_np_3   
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1, e_np_1
   
        select case (this%periodic)
        case (.TRUE.)
            a10 = aD04d2 * this%onebydx2
            b10 = bD04d2 * this%onebydx2

            RHS(:,:,1         ) = a10 * ( f(:,:,2          )   - two*f(:,:,1            ) + f(:,:,this%n    )) &
                                + b10 * ( f(:,:,3          )   - two*f(:,:,1            ) + f(:,:,this%n-1  ))
            
            RHS(:,:,2         ) = a10 * ( f(:,:,3          )   - two*f(:,:,2            ) + f(:,:,1         )) &
                                + b10 * ( f(:,:,4          )   - two*f(:,:,2            ) + f(:,:,this%n    ))
            
            RHS(:,:,3         ) = a10 * ( f(:,:,4          )   - two*f(:,:,3            ) + f(:,:,2         )) &
                                + b10 * ( f(:,:,5          )   - two*f(:,:,3            ) + f(:,:,1         ))
            
            RHS(:,:,4:this%n-3) = a10 * ( f(:,:,5:this%n-2 )   - two*f(:,:,4:this%n-3   ) + f(:,:,3:this%n-4)) &
                                + b10 * ( f(:,:,6:this%n-1 )   - two*f(:,:,4:this%n-3   ) + f(:,:,2:this%n-5))
            
            RHS(:,:,this%n-2  ) = a10 * ( f(:,:,this%n-1   )   - two*f(:,:,this%n-2     ) + f(:,:,this%n-3  )) &
                                + b10 * ( f(:,:,this%n     )   - two*f(:,:,this%n-2     ) + f(:,:,this%n-4  ))
            
            RHS(:,:,this%n-1  ) = a10 * ( f(:,:,this%n     )   - two*f(:,:,this%n-1     ) + f(:,:,this%n-2  )) &
                                + b10 * ( f(:,:,1          )   - two*f(:,:,this%n-1     ) + f(:,:,this%n-3  ))
            
            RHS(:,:,this%n    ) = a10 * ( f(:,:,1          )   - two*f(:,:,this%n       ) + f(:,:,this%n-1  )) &
                                + b10 * ( f(:,:,2          )   - two*f(:,:,this%n       ) + f(:,:,this%n-2  ))

        case (.FALSE.)
            a10 = aD04d2 * this%onebydx2
            b10 = bD04d2 * this%onebydx2
            
            a_np_3 = b3_aD04d2 * this%onebydx2 
            b_np_3 = b3_bD04d2 * this%onebydx2 
       
            a_np_2 = b2_aD04d2 * this%onebydx2

            a_np_1 = b1_aD04d2 * this%onebydx2
            b_np_1 = b1_bD04d2 * this%onebydx2
            c_np_1 = b1_cD04d2 * this%onebydx2
            d_np_1 = b1_dD04d2 * this%onebydx2
            e_np_1 = b1_eD04d2 * this%onebydx2
                    
            select case(bc1)
            case(0)
            RHS(:,:,1         ) = a_np_1*f(:,:,1) + b_np_1*f(:,:,2) + &
                                & c_np_1*f(:,:,3) + d_np_1*f(:,:,4) + &
                                & e_np_1*f(:,:,5)
            
            RHS(:,:,2         ) = a_np_2 * ( f(:,:,3          )   - two*f(:,:,2            ) + f(:,:,1         )) 
            
            RHS(:,:,3         ) = a_np_3 * ( f(:,:,4          )   - two*f(:,:,3            ) + f(:,:,2         )) &
                                + b_np_3 * ( f(:,:,5          )   - two*f(:,:,3            ) + f(:,:,1         )) 
            case(1)
                RHS(:,:,1) = a10 * ( f(:,:,2)   - two*f(:,:,1) + f(:,:,2)) &
                           + b10 * ( f(:,:,3)   - two*f(:,:,1) + f(:,:,3))
                
                RHS(:,:,2) = a10 * ( f(:,:,3)   - two*f(:,:,2) + f(:,:,1)) &
                           + b10 * ( f(:,:,4)   - two*f(:,:,2) + f(:,:,2))
                
                RHS(:,:,3) = a10 * ( f(:,:,4)   - two*f(:,:,3) + f(:,:,2)) &
                           + b10 * ( f(:,:,5)   - two*f(:,:,3) + f(:,:,1))
            case(-1)
                RHS(:,:,1) = a10 * ( f(:,:,2)   - two*f(:,:,1) - f(:,:,2)) &
                           + b10 * ( f(:,:,3)   - two*f(:,:,1) - f(:,:,3))
                
                RHS(:,:,2) = a10 * ( f(:,:,3)   - two*f(:,:,2) + f(:,:,1)) &
                           + b10 * ( f(:,:,4)   - two*f(:,:,2) - f(:,:,2))
                
                RHS(:,:,3) = a10 * ( f(:,:,4)   - two*f(:,:,3) + f(:,:,2)) &
                           + b10 * ( f(:,:,5)   - two*f(:,:,3) + f(:,:,1))
            end select
            
            RHS(:,:,4:this%n-3) = a10 * ( f(:,:,5:this%n-2 )   - two*f(:,:,4:this%n-3   ) + f(:,:,3:this%n-4)) &
                                + b10 * ( f(:,:,6:this%n-1 )   - two*f(:,:,4:this%n-3   ) + f(:,:,2:this%n-5))
            
            select case(bcn)
            case(0)
            RHS(:,:,this%n-2  ) = a_np_3 * ( f(:,:,this%n-1   )   - two*f(:,:,this%n-2     ) + f(:,:,this%n-3  )) &
                                + b_np_3 * ( f(:,:,this%n     )   - two*f(:,:,this%n-2     ) + f(:,:,this%n-4  )) 
            
            RHS(:,:,this%n-1  ) = a_np_2 * ( f(:,:,this%n     )   - two*f(:,:,this%n-1     ) + f(:,:,this%n-2  )) 
            
            RHS(:,:,this%n    ) =  a_np_1*f(:,:,this%n  ) + b_np_1*f(:,:,this%n-1) + &
                                &  c_np_1*f(:,:,this%n-2) + d_np_1*f(:,:,this%n-3) + &
                                &  e_np_1*f(:,:,this%n-4)
            case(1)
                RHS(:,:,this%n-2) = a10 * ( f(:,:,this%n-1)   - two*f(:,:,this%n-2) + f(:,:,this%n-3)) &
                                  + b10 * ( f(:,:,this%n  )   - two*f(:,:,this%n-2) + f(:,:,this%n-4))
                
                RHS(:,:,this%n-1) = a10 * ( f(:,:,this%n  )   - two*f(:,:,this%n-1) + f(:,:,this%n-2)) &
                                  + b10 * ( f(:,:,this%n-1)   - two*f(:,:,this%n-1) + f(:,:,this%n-3))
                
                RHS(:,:,this%n  ) = a10 * ( f(:,:,this%n-1)   - two*f(:,:,this%n  ) + f(:,:,this%n-1)) &
                                  + b10 * ( f(:,:,this%n-2)   - two*f(:,:,this%n  ) + f(:,:,this%n-2))
            case(-1)
                RHS(:,:,this%n-2) = a10 * ( f(:,:,this%n-1)   - two*f(:,:,this%n-2) + f(:,:,this%n-3)) &
                                  + b10 * ( f(:,:,this%n  )   - two*f(:,:,this%n-2) + f(:,:,this%n-4))
                
                RHS(:,:,this%n-1) = a10 * ( f(:,:,this%n  )   - two*f(:,:,this%n-1) + f(:,:,this%n-2)) &
                                  + b10 * (-f(:,:,this%n-1)   - two*f(:,:,this%n-1) + f(:,:,this%n-3))
                
                RHS(:,:,this%n  ) = a10 * (-f(:,:,this%n-1)   - two*f(:,:,this%n  ) + f(:,:,this%n-1)) &
                                  + b10 * (-f(:,:,this%n-2)   - two*f(:,:,this%n  ) + f(:,:,this%n-2))
            end select

        end select
   
    end subroutine  
    
    subroutine dd1(this, f, df, na, nb, bc1_, bcn_)
        class( d04 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(this%n,na,nb), intent(in)  :: f
        real(rkind), dimension(this%n,na,nb), intent(out) :: df

        if(this%n == 1) then
            df = zero
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

        call this%ComputeXD1RHS(f, df, na, nb, bc1, bcn)
        
    end subroutine

    subroutine dd2(this, f, df, na, nb, bc1_, bcn_)
        class( d04 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(na,this%n,nb), intent(in) :: f
        real(rkind), dimension(na,this%n,nb), intent(out) :: df

        if(this%n == 1) then
            df = zero
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

        call this%ComputeYD1RHS(f, df, na, nb, bc1, bcn)

    end subroutine

    subroutine dd3(this, f, df, na, nb, bc1_, bcn_)
        class( d04 ), intent(in) :: this
        integer, intent(in) :: na, nb
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        real(rkind), dimension(na,nb,this%n), intent(in) :: f
        real(rkind), dimension(na,nb,this%n), intent(out) :: df

        if(this%n == 1) then
            df = zero
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

        call this%ComputeZD1RHS(f, df, na, nb, bc1, bcn)

    end subroutine

    subroutine d2d1(this, f, df, na, nb, bc1_, bcn_)
        class( d04 ), intent(in) :: this
        integer, intent(in) :: na, nb
        real(rkind), dimension(this%n,na,nb), intent(in) :: f
        real(rkind), dimension(this%n,na,nb), intent(out) :: df
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn

        if(this%n == 1) then
            df = zero
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

        call this%ComputeXD2RHS(f, df, na, nb, bc1, bcn)

    end subroutine

    subroutine d2d2(this, f, df, na, nb, bc1_, bcn_)
        class( d04 ), intent(in) :: this
        integer, intent(in) :: na, nb
        real(rkind), dimension(na,this%n,nb), intent(in) :: f
        real(rkind), dimension(na,this%n,nb), intent(out) :: df
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn

        if(this%n == 1) then
            df = zero
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

        call this%ComputeYD2RHS(f, df, na, nb, bc1, bcn)

    end subroutine

    subroutine d2d3(this, f, df, na, nb, bc1_, bcn_)
        class( d04 ), intent(in) :: this
        integer, intent(in) :: na, nb
        real(rkind), dimension(na,nb,this%n), intent(in) :: f
        real(rkind), dimension(na,nb,this%n), intent(out) :: df
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn

        if(this%n == 1) then
            df = zero
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

        call this%ComputeZD2RHS(f, df, na, nb, bc1, bcn)

    end subroutine


end module
