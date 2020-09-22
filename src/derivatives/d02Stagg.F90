! Routines specific to 2nd order Explicit Staggered Finite Difference Scheme

module d02Staggstuff

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two
    use exits,           only: GracefulExit
    
    implicit none

    private
    public :: d02Stagg, aD02d1
    ! 2nd order first derivative explicit centeral difference coefficients
    real(rkind), parameter :: aD02d1     = 1.0_rkind

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic 1st derivative evaluation !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
    ! Interior scheme
    real(rkind), parameter                   :: q_hat           = aD02d1

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
    real(rkind), parameter                   :: q_ppp           = aD02d1

    ! Weights are now trivial, but fluxes do not telescope
    real(rkind), parameter                   :: w1              = 1.0_rkind
    real(rkind), parameter                   :: w2              = 1.0_rkind
    real(rkind), parameter                   :: w3              = 1.0_rkind
    real(rkind), parameter                   :: w4              = 1.0_rkind
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type d02Stagg
        
        private
        
        integer     :: n
        real(rkind) :: dx
        real(rkind) :: onebydx

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
        
        procedure :: dd1N2F
        procedure :: dd2N2F
        procedure :: dd3N2F
        procedure :: dd1F2N
        procedure :: dd2F2N
        procedure :: dd3F2N
        
    end type



contains
    
    pure function GetSize(this) result(val)
        class(d02Stagg), intent(in) :: this
        integer  :: val 
        val = this%n
    end function

    function init(this, n_, dx_, periodic_, bc1_, bcn_) result(ierr)
   
        class( d02Stagg ), intent(inout) :: this
        integer, intent(in) :: n_
        real(rkind), intent(in) :: dx_
        logical, intent(in) :: periodic_
        integer, intent(in) :: bc1_, bcn_
        integer :: ierr
        
        this%n = n_
        this%dx = dx_
        this%onebydx = one/dx_

        this%periodic = periodic_

        this%bc1 = bc1_
        this%bcn = bcn_

        ! If everything passes
        ierr = 0
    
    end function
    
    subroutine destroy(this)

        class( d02Stagg ), intent(inout) :: this

        return

    end subroutine


    pure subroutine ComputeXD1RHS(this, f, RHS, dir, n2, n3, bc1, bcn)   !TODO add an extra argument to account for forward or backward
    
        class( d02Stagg ), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10
        integer :: j,k

        select case (this%periodic)
        case (.TRUE.)
            a10 = aD02d1 * this%onebydx
            select case (dir)
                case ("N2F")
                    RHS(1:this%n-1,:,:) = a10 * (f(2:this%n,:,:) - f(1:this%n-1,:,:))
                    RHS(this%n,:,:)     = a10 * (f(1,       :,:) - f(this%n    ,:,:))
                case ("F2N")
                    RHS(2:this%n,:,:)   = a10 * (f(2:this%n,:,:) - f(1:this%n-1,:,:))
                    RHS(1       ,:,:)   = a10 * (f(1,       :,:) - f(this%n    ,:,:))
            end select
        case (.FALSE.)
            a10    = q_hat * this%onebydx  
            select case (dir)
                case ("N2F")
                    !2nd Order                        
                    RHS(1:this%n-1,:,:) = a10 * (f(2:this%n,:,:) - f(1:this%n-1,:,:))
                case ("F2N")
                    !2nd Order interior, 1st order at boundary                        
                    RHS(2:this%n-1,:,:) = a10 * (f(2:this%n-1,:,:) - f(1:this%n-2,:,:))
                    RHS(1         ,:,:) = RHS(2,:,:)
                    RHS(this%n    ,:,:) = RHS(this%n-1,:,:)
            end select

        end select
    
    end subroutine
    
    pure subroutine ComputeYD1RHS(this, f, RHS, dir, n1, n3, bc1, bcn)  !TODO add an extra argument to account for forward or backward
    
        class( d02Stagg ), intent(in) :: this
        integer, intent(in) :: n1, n3
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10
        integer :: k
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4
        real(rkind) :: a_np_3
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1


        select case (this%periodic)
        case (.TRUE.)
            a10 = aD02d1 * this%onebydx
            select case (dir)
                case ("N2F")
                    RHS(:,1:this%n-1,:) = a10 * (f(:,2:this%n,:) - f(:,1:this%n-1,:))
                    RHS(:,this%n,    :) = a10 * (f(:,1,       :) - f(:,this%n    ,:))
                case ("F2N")
                    RHS(:,2:this%n,:)   = a10 * (f(:,2:this%n,:) - f(:,1:this%n-1,:))
                    RHS(:,1       ,:)   = a10 * (f(:,1,       :) - f(:,this%n    ,:))
            end select
        case (.FALSE.)
            a10    = q_hat * this%onebydx  
            select case (dir)
                case ("N2F")
                    !2nd Order                        
                    RHS(:,1:this%n-1,:) = a10 * (f(:,2:this%n,:) - f(:,1:this%n-1,:))
                case ("F2N")
                    !2nd Order interior, 1st order at boundary                        
                    RHS(:,2:this%n-1,:) = a10 * (f(:,2:this%n-1,:) - f(:,1:this%n-2,:))
                    RHS(:,1         ,:) = RHS(:,2,:)
                    RHS(:,this%n    ,:) = RHS(:,this%n-1,:)
            end select
        end select
    
    end subroutine

    !pure subroutine ComputeZD1RHS(this, f, RHS, n1, n2, bc1, bcn)
    subroutine ComputeZD1RHS(this, f, RHS, dir, n1, n2, bc1, bcn)   !TODO add an extra argument to account for forward or backward
    
        class( d02Stagg ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in) :: f
        real(rkind), dimension(n1,n2,this%n), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4
        real(rkind) :: a_np_3
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1

        select case (this%periodic)
        case (.TRUE.)
            a10 = aD02d1 * this%onebydx
            select case (dir)
                case ("N2F")
                    RHS(:,:,1:this%n-1) = a10 * (f(:,:,2:this%n) - f(:,:,1:this%n-1))
                    RHS(:,:,this%n    ) = a10 * (f(:,:,1       ) - f(:,:,this%n    ))
                case ("F2N")
                    RHS(:,:,2:this%n)   = a10 * (f(:,:,2:this%n) - f(:,:,1:this%n-1))
                    RHS(:,:,1       )   = a10 * (f(:,:,1       ) - f(:,:,this%n    ))
            end select
        case (.FALSE.)
            a10    = q_hat * this%onebydx  
            select case (dir)
                case ("N2F")
                    !2nd Order                        
                    RHS(:,:,1:this%n-1) = a10 * (f(:,:,2:this%n) - f(:,:,1:this%n-1))
                case ("F2N")
                    !2nd Order interior, 1st order at boundary                        
                    RHS(:,:,2:this%n-1) = a10 * (f(:,:,2:this%n-1) - f(:,:,1:this%n-2))
                    RHS(:,:,1         ) = RHS(:,:,2)
                    RHS(:,:,this%n    ) = RHS(:,:,this%n-1)
            end select
            
        end select
    
    end subroutine
    
    subroutine dd1N2F(this, f, df, na, nb, bc1_, bcn_)
        class( d02Stagg ), intent(in) :: this
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
            if ( (bc1 .eq. 1) .OR. (bc1 .eq. -1) ) then
                call GracefulExit("Symmetric/Antisymmetric BC not yet implemented for staggered scheme", 324)
            endif
            if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bc1 (should be 0, 1 or -1)", 324)
            end if
        else
            bc1 = 0
        end if

        if (present(bcn_)) then
            bcn = bcn_
            if ( (bcn .eq. 1) .OR. (bcn .eq. -1) ) then
                call GracefulExit("Symmetric/Antisymmetric BC not yet implemented for staggered scheme", 324)
            endif
            if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bcn (should be 0, 1 or -1)", 324)
            end if
        else
            bcn = 0
        end if
        
        call this%ComputeXD1RHS(f, df, "N2F", na, nb, bc1, bcn)

    end subroutine

    subroutine dd2N2F(this, f, df, na, nb, bc1_, bcn_)
        class( d02Stagg ), intent(in) :: this
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
            if ( (bc1 .eq. 1) .OR. (bc1 .eq. -1) ) then
                call GracefulExit("Symmetric/Antisymmetric BC not yet implemented for staggered scheme", 324)
            endif
            if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bc1 (should be 0, 1 or -1)", 324)
            end if
        else
            bc1 = 0
        end if

        if (present(bcn_)) then
            bcn = bcn_
            if ( (bcn .eq. 1) .OR. (bcn .eq. -1) ) then
                call GracefulExit("Symmetric/Antisymmetric BC not yet implemented for staggered scheme", 324)
            endif
            if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bcn (should be 0, 1 or -1)", 324)
            end if
        else
            bcn = 0
        end if
        
        
        call this%ComputeYD1RHS(f, df, "N2F", na, nb, bc1, bcn)

    end subroutine

    subroutine dd3N2F(this, f, df, na, nb, bc1_, bcn_)
        class( d02Stagg ), intent(in) :: this
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
            if ( (bc1 .eq. 1) .OR. (bc1 .eq. -1) ) then
                call GracefulExit("Symmetric/Antisymmetric BC not yet implemented for staggered scheme", 324)
            endif
            if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bc1 (should be 0, 1 or -1)", 324)
            end if
        else
            bc1 = 0
        end if

        if (present(bcn_)) then
            bcn = bcn_
            if ( (bcn .eq. 1) .OR. (bcn .eq. -1) ) then
                call GracefulExit("Symmetric/Antisymmetric BC not yet implemented for staggered scheme", 324)
            endif
            if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bcn (should be 0, 1 or -1)", 324)
            end if
        else
            bcn = 0
        end if
        

        call this%ComputeZD1RHS(f, df, "N2F", na, nb, bc1, bcn)

    end subroutine

    subroutine dd1F2N(this, f, df, na, nb, bc1_, bcn_)
        class( d02Stagg ), intent(in) :: this
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
            if ( (bc1 .eq. 1) .OR. (bc1 .eq. -1) ) then
                call GracefulExit("Symmetric/Antisymmetric BC not yet implemented for staggered scheme", 324)
            endif
            if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bc1 (should be 0, 1 or -1)", 324)
            end if
        else
            bc1 = 0
        end if

        if (present(bcn_)) then
            bcn = bcn_
            if ( (bcn .eq. 1) .OR. (bcn .eq. -1) ) then
                call GracefulExit("Symmetric/Antisymmetric BC not yet implemented for staggered scheme", 324)
            endif
            if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bcn (should be 0, 1 or -1)", 324)
            end if
        else
            bcn = 0
        end if
        
        
        call this%ComputeXD1RHS(f, df, "F2N", na, nb, bc1, bcn)

    end subroutine

    subroutine dd2F2N(this, f, df, na, nb, bc1_, bcn_)
        class( d02Stagg ), intent(in) :: this
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
            if ( (bc1 .eq. 1) .OR. (bc1 .eq. -1) ) then
                call GracefulExit("Symmetric/Antisymmetric BC not yet implemented for staggered scheme", 324)
            endif
            if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bc1 (should be 0, 1 or -1)", 324)
            end if
        else
            bc1 = 0
        end if

        if (present(bcn_)) then
            bcn = bcn_
            if ( (bcn .eq. 1) .OR. (bcn .eq. -1) ) then
                call GracefulExit("Symmetric/Antisymmetric BC not yet implemented for staggered scheme", 324)
            endif
            if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bcn (should be 0, 1 or -1)", 324)
            end if
        else
            bcn = 0
        end if
        
        
        call this%ComputeYD1RHS(f, df, "F2N", na, nb, bc1, bcn)

    end subroutine

    subroutine dd3F2N(this, f, df, na, nb, bc1_, bcn_)
        class( d02Stagg ), intent(in) :: this
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
            if ( (bc1 .eq. 1) .OR. (bc1 .eq. -1) ) then
                call GracefulExit("Symmetric/Antisymmetric BC not yet implemented for staggered scheme", 324)
            endif
            if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bc1 (should be 0, 1 or -1)", 324)
            end if
        else
            bc1 = 0
        end if

        if (present(bcn_)) then
            bcn = bcn_
            if ( (bcn .eq. 1) .OR. (bcn .eq. -1) ) then
                call GracefulExit("Symmetric/Antisymmetric BC not yet implemented for staggered scheme", 324)
            endif
            if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                call GracefulExit("Incorrect boundary specification for bcn (should be 0, 1 or -1)", 324)
            end if
        else
            bcn = 0
        end if
        

        call this%ComputeZD1RHS(f, df, "F2N", na, nb, bc1, bcn)

    end subroutine


end module
