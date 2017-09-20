! Routines specific to 10th order Compact Finite Differencing scheme
! Periodic LU based on  Neossi Nguetchue, Abelman (Appl. Math. & Comp. 2008)

module ed02stuff

    use kind_parameters, only: rkind
    use constants,       only: zero,half,one,two
    use exits,           only: GracefulExit
    
    implicit none

    private
    public :: ed02 

    type ed02
        
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
        class(ed02), intent(in) :: this
        integer  :: val 
        val = this%n
    end function

    function init(this, n_, dx_, periodic_, bc1_, bcn_) result(ierr)
   
        class(ed02), intent(inout) :: this
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
        class(ed02), intent(inout) :: this

        this%periodic = .FALSE.

    end subroutine

    pure subroutine ComputeXD1RHS(this, f, RHS, n2, n3, bc1, bcn)
    
        class(ed02), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a02
        integer :: j,k
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_1, b_np_1

        a02 = half * this%onebydx

        select case (this%periodic)
        case (.TRUE.)
            do k=1,n3
                do j=1,n2
                    RHS(1         ,j,k) = a02 * ( f(2       ,j,k) - f(this%n    ,j,k) )
                    RHS(2:this%n-1,j,k) = a02 * ( f(3:this%n,j,k) - f(1:this%n-2,j,k) )
                    RHS(this%n    ,j,k) = a02 * ( f(1       ,j,k) - f(this%n-1  ,j,k) )
                end do
            end do
        case (.FALSE.)
            a_np_1 = -one * this%onebydx
            b_np_1 =  one * this%onebydx

            do k = 1,n3
                do j = 1,n2
                    select case(bc1)
                    case(0)
                        RHS(1,j,k) = a_np_1* f(1,j,k) +  b_np_1*f(2,j,k)
                    case(1)
                        RHS(1,j,k) = zero
                    case(-1)
                        RHS(1,j,k) = a02*   (f(2,j,k) +         f(2,j,k))
                    end select
                    
                    RHS(2:this%n-1,j,k) = a02*(f(3:this%n  ,j,k) - f(1:this%n-2,j,k))
                    
                    select case(bcn)
                    case(0)
                        RHS(this%n,j,k) =  -a_np_1* f(this%n   ,j,k) -  b_np_1*f(this%n-1,j,k)
                    case(1)
                        RHS(this%n,j,k) =   zero
                    case(-1)
                        RHS(this%n,j,k) =   a02   *(-f(this%n-1,j,k) -         f(this%n-1,j,k))
                    end select
               end do 
            end do 
        end select
    
    end subroutine
    
    pure subroutine ComputeYD1RHS(this, f, RHS, n1, n3, bc1, bcn)
    
        class(ed02), intent(in) :: this
        integer, intent(in) :: n1, n3
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a02
        integer :: k
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_1, b_np_1

        a02 = half * this%onebydx

        select case (this%periodic)
        case (.TRUE.)
            do k=1,n3
                RHS(:,1         ,k) = a02 * ( f(:,2       ,k) - f(:,this%n    ,k) )
                RHS(:,2:this%n-1,k) = a02 * ( f(:,3:this%n,k) - f(:,1:this%n-2,k) )
                RHS(:,this%n    ,k) = a02 * ( f(:,1       ,k) - f(:,this%n-1  ,k) )
            end do
        case (.FALSE.)
            a_np_1 = -one * this%onebydx
            b_np_1 =  one * this%onebydx

            do k = 1,n3
                select case(bc1)
                case(0)
                    RHS(:,1,k) = a_np_1* f(:,1,k) +  b_np_1*f(:,2,k)
                case(1)
                    RHS(:,1,k) = zero
                case(-1)
                    RHS(:,1,k) = a02*   (f(:,2,k) +         f(:,2,k))
                end select
                
                RHS(:,2:this%n-1,k) = a02*(f(:,3:this%n  ,k) - f(:,1:this%n-2,k))
                
                select case(bcn)
                case(0)
                    RHS(:,this%n,k) =  -a_np_1*  f(:,this%n   ,k) - b_np_1*f(:,this%n-1,k)
                case(1)
                    RHS(:,this%n,k) =   zero
                case(-1)
                    RHS(:,this%n,k) =   a02   *(-f(:,this%n-1,k) -         f(:,this%n-1,k))
                end select
            end do 
        end select
    
    end subroutine
    
    pure subroutine ComputeZD1RHS(this, f, RHS, n1, n2, bc1, bcn)
    
        class(ed02), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in) :: f
        real(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a02
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_1, b_np_1

        a02 = half * this%onebydx

        select case (this%periodic)
        case (.TRUE.)
            RHS(:,:,1         ) = a02 * ( f(:,:,2       ) - f(:,:,this%n    ) )
            RHS(:,:,2:this%n-1) = a02 * ( f(:,:,3:this%n) - f(:,:,1:this%n-2) )
            RHS(:,:,this%n    ) = a02 * ( f(:,:,1       ) - f(:,:,this%n-1  ) )
        case (.FALSE.)
            a_np_1 = -one * this%onebydx
            b_np_1 =  one * this%onebydx

            select case(bc1)
            case(0)
                RHS(:,:,1) = a_np_1* f(:,:,1) +  b_np_1*f(:,:,2)
            case(1)
                RHS(:,:,1) = zero
            case(-1)
                RHS(:,:,1) = a02*   (f(:,:,2) +         f(:,:,2))
            end select
            
            RHS(:,:,2:this%n-1) = a02*(f(:,:,3:this%n  ) - f(:,:,1:this%n-2))
            
            select case(bcn)
            case(0)
                RHS(:,:,this%n) =  -a_np_1*  f(:,:,this%n  ) - b_np_1*f(:,:,this%n-1)
            case(1)
                RHS(:,:,this%n) =   zero
            case(-1)
                RHS(:,:,this%n) =   a02   *(-f(:,:,this%n-1) -        f(:,:,this%n-1))
            end select
        end select
    
    end subroutine
    
    
   pure subroutine ComputeXD2RHS(this, f, RHS, n2, n3, bc1, bcn)
        class(ed02), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a02
        real(rkind) :: a_np_1, b_np_1, c_np_1
        integer :: j,k 

        a02 = one * this%onebydx2
    
        select case (this%periodic)
        case (.TRUE.)

            do k = 1,n3
                do j = 1,n2
                    RHS(1         ,j,k) = a02 * ( f(2          ,j,k)   - two*f(1            ,j,k) + f(this%n    ,j,k))
                    
                    RHS(2:this%n-1,j,k) = a02 * ( f(3:this%n   ,j,k)   - two*f(2:this%n-1   ,j,k) + f(1:this%n-2,j,k))
                   
                    RHS(this%n    ,j,k) = a02 * ( f(1          ,j,k)   - two*f(this%n       ,j,k) + f(this%n-1  ,j,k))
                end do 
            end do 

        case (.FALSE.)
            a_np_1 =  one * this%onebydx2
            b_np_1 = -two * this%onebydx2
            c_np_1 =  one * this%onebydx2

            do k = 1,n3
                do j = 1,n2
                    select case(bc1)
                    case(0)
                        RHS(1,j,k) = a_np_1*f(1,j,k) + b_np_1*f(2,j,k) + c_np_1*f(3,j,k)
                    case(1)
                        RHS(1,j,k) = a02 * (f(2,j,k) -    two*f(1,j,k) +        f(2,j,k))
                    case(-1)
                        RHS(1,j,k) = a02 * (f(2,j,k) -    two*f(1,j,k) -        f(2,j,k))
                    end select
                    
                    RHS(2:this%n-1,j,k) = a02 * ( f(3:this%n   ,j,k)   - two*f(2:this%n-1   ,j,k) + f(1:this%n-2,j,k))
                    
                    select case(bcn)
                    case(0) 
                        RHS(this%n  ,j,k) = a_np_1* f(this%n  ,j,k) + b_np_1*f(this%n-1,j,k) + c_np_1*f(this%n-2,j,k)
                    case(1)
                        RHS(this%n  ,j,k) = a02 * ( f(this%n-1,j,k) -    two*f(this%n  ,j,k) +        f(this%n-1,j,k))
                    case(-1)
                        RHS(this%n  ,j,k) = a02 * (-f(this%n-1,j,k) -    two*f(this%n  ,j,k) +        f(this%n-1,j,k))
                    end select

                end do 
            end do 
        end select
   
    end subroutine  

   pure subroutine ComputeYD2RHS(this, f, RHS, n1, n3, bc1, bcn)
        class(ed02), intent(in) :: this
        integer, intent(in) :: n1, n3
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a02
        real(rkind) :: a_np_1, b_np_1, c_np_1
        integer :: k 

        a02 = one * this%onebydx2
    
        select case (this%periodic)
        case (.TRUE.)

            do k = 1,n3
                RHS(:,1         ,k) = a02 * ( f(:,2          ,k)   - two*f(:,1            ,k) + f(:,this%n    ,k))

                RHS(:,2:this%n-1,k) = a02 * ( f(:,3:this%n   ,k)   - two*f(:,2:this%n-1   ,k) + f(:,1:this%n-2,k))

                RHS(:,this%n    ,k) = a02 * ( f(:,1          ,k)   - two*f(:,this%n       ,k) + f(:,this%n-1  ,k))
            end do 

        case (.FALSE.)
            a_np_1 =  one * this%onebydx2
            b_np_1 = -two * this%onebydx2
            c_np_1 =  one * this%onebydx2

            do k = 1,n3
                select case(bc1)
                case(0)
                    RHS(:,1,k) = a_np_1*f(:,1,k) + b_np_1*f(:,2,k) + c_np_1*f(:,3,k)
                case(1)
                    RHS(:,1,k) = a02 * (f(:,2,k) -    two*f(:,1,k) +        f(:,2,k))
                case(-1)
                    RHS(:,1,k) = a02 * (f(:,2,k) -    two*f(:,1,k) -        f(:,2,k))
                end select
                
                RHS(:,2:this%n-1,k) = a02 * ( f(:,3:this%n   ,k)   - two*f(:,2:this%n-1   ,k) + f(:,1:this%n-2,k))
                
                select case(bcn)
                case(0) 
                    RHS(:,this%n  ,k) = a_np_1* f(:,this%n  ,k) + b_np_1*f(:,this%n-1,k) + c_np_1*f(:,this%n-2,k)
                case(1)
                    RHS(:,this%n  ,k) = a02 * ( f(:,this%n-1,k) -    two*f(:,this%n  ,k) +        f(:,this%n-1,k))
                case(-1)
                    RHS(:,this%n  ,k) = a02 * (-f(:,this%n-1,k) -    two*f(:,this%n  ,k) +        f(:,this%n-1,k))
                end select

            end do 
        end select
   
    end subroutine  

   pure subroutine ComputeZD2RHS(this, f, RHS, n1, n2, bc1, bcn)
        class(ed02), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in) :: f
        real(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a02
        real(rkind) :: a_np_1, b_np_1, c_np_1

        a02 = one * this%onebydx2
    
        select case (this%periodic)
        case (.TRUE.)
            RHS(:,:,1         ) = a02 * ( f(:,:,2          )   - two*f(:,:,1            ) + f(:,:,this%n    ))

            RHS(:,:,2:this%n-1) = a02 * ( f(:,:,3:this%n   )   - two*f(:,:,2:this%n-1   ) + f(:,:,1:this%n-2))

            RHS(:,:,this%n    ) = a02 * ( f(:,:,1          )   - two*f(:,:,this%n       ) + f(:,:,this%n-1  ))
        case (.FALSE.)
            a_np_1 =  one * this%onebydx2
            b_np_1 = -two * this%onebydx2
            c_np_1 =  one * this%onebydx2

            select case(bc1)
            case(0)
                RHS(:,:,1) = a_np_1*f(:,:,1) + b_np_1*f(:,:,2) + c_np_1*f(:,:,3)
            case(1)
                RHS(:,:,1) = a02 * (f(:,:,2) -    two*f(:,:,1) +        f(:,:,2))
            case(-1)
                RHS(:,:,1) = a02 * (f(:,:,2) -    two*f(:,:,1) -        f(:,:,2))
            end select
            
            RHS(:,:,2:this%n-1) = a02 * ( f(:,:,3:this%n   )   - two*f(:,:,2:this%n-1   ) + f(:,:,1:this%n-2))
            
            select case(bcn)
            case(0) 
                RHS(:,:,this%n  ) = a_np_1* f(:,:,this%n  ) + b_np_1*f(:,:,this%n-1) + c_np_1*f(:,:,this%n-2)
            case(1)
                RHS(:,:,this%n  ) = a02 * ( f(:,:,this%n-1) -    two*f(:,:,this%n  ) +        f(:,:,this%n-1))
            case(-1)
                RHS(:,:,this%n  ) = a02 * (-f(:,:,this%n-1) -    two*f(:,:,this%n  ) +        f(:,:,this%n-1))
            end select
        end select
   
    end subroutine  


    subroutine dd1(this, f, df, na, nb, bc1_, bcn_)
        class(ed02), intent(in) :: this
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
        class(ed02), intent(in) :: this
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
        class(ed02), intent(in) :: this
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
        class(ed02), intent(in) :: this
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
        class(ed02), intent(in) :: this
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
        class(ed02), intent(in) :: this
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
