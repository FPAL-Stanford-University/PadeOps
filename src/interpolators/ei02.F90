! Routines specific to 2nd order Explicit Midpoint Interpolation Scheme

module ei02stuff

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two
    use exits,           only: GracefulExit
    
    implicit none

    private
    public :: ei02, aI02 !TODO Rename this last var
    ! 2nd order first derivative explicit centeral difference coefficients
    real(rkind), parameter :: aI02     = 1.0_rkind / 2.0_rkind

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

    type ei02
        
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
        class(ei02), intent(in) :: this
        integer  :: val 
        val = this%n
    end function

    function init(this, n_, dx_, periodic_, bc1_, bcn_) result(ierr)
   
        class( ei02 ), intent(inout) :: this
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

        class( ei02 ), intent(inout) :: this

        return

    end subroutine


    pure subroutine ComputeXIRHS(this, f, RHS, dir, n2, n3, bc1, bcn) !TODO add an extra argument to account for forward or backward
    
        class( ei02 ), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10
        integer :: j,k

        RHS = 0.0d0
        select case (this%periodic)
        case (.TRUE.)
            a10 = aI02
            do k=1,n3
                do j=1,n2
                    select case(dir)
                    case("N2F")
                        RHS(1:this%n-1,j,k) = a10 * ( f(1:this%n-1,j,k) + f(2:this%n,j,k) )
                        RHS(this%n,j,k)     = a10 * ( f(this%n,j,k)     + f(1,j,k) )
                    case("F2N")
                        RHS(1,j,k)          = a10 * ( f(this%n,j,k)     + f(1,j,k) )
                        RHS(2:this%n,j,k)   = a10 * ( f(1:this%n-1,j,k) + f(2:this%n,j,k) )
                    end select
                end do
            end do

        case (.FALSE.) 
            a10    = q_hat 
            select case(dir)
            case("N2F")
                !2nd order    
                RHS(1:this%n-1,:,:) = a10 * ( f(1:this%n-1,:,:) + f(2:this%n,:,:) )  
            case("F2N")
                !2nd order interior, 1st order at boundary
                RHS(2:this%n-1,:,:) = a10 * ( f(1:this%n-2,:,:) + f(2:this%n-1,:,:) )  

                RHS(1,:,:)       = f(1,:,:)
                RHS(this%n,:,:) = f(this%n,:,:)

            end select

        end select
    
    end subroutine
    
    pure subroutine ComputeYIRHS(this, f, RHS, dir, n1, n3, bc1, bcn) 
    
        class( ei02 ), intent(in) :: this
        integer, intent(in) :: n1, n3
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: RHS
        character(len=*)  , intent(in)             :: dir
        integer, intent(in) :: bc1, bcn
        real(rkind) :: a10
        integer :: i,k
        ! Non-periodic boundary a, b and c
        real(rkind) :: a_np_4
        real(rkind) :: a_np_3
        real(rkind) :: a_np_2
        real(rkind) :: a_np_1, b_np_1, c_np_1, d_np_1


        select case (this%periodic)
        case (.TRUE.)
            a10 = aI02
            do k=1,n3
                do i=1,n1
                    select case(dir)
                    case("N2F")
                        RHS(i,1:this%n-1,k) = a10 * ( f(i,1:this%n-1,k) + f(i,2:this%n,k) )
                        RHS(i,this%n,k)     = a10 * ( f(i,this%n,k)     + f(i,1,k) )
                    case("F2N")
                        RHS(i,1,k)          = a10 * ( f(i,this%n,k)     + f(i,1,k) )
                        RHS(i,2:this%n,k)   = a10 * ( f(i,1:this%n-1,k) + f(i,2:this%n,k) )
                    end select
                end do
            end do

        case (.FALSE.)  !TODO: implement non-periodic BC
            a10    = q_hat 
            select case(dir)
            case("N2F")
                !2nd order    
                RHS(:,1:this%n-1,:) = a10 * ( f(:,1:this%n-1,:) + f(:,2:this%n,:) )  
            case("F2N")
                !2nd order interior, 1st order at boundary
                RHS(:,2:this%n-1,:) = a10 * ( f(:,1:this%n-2,:) + f(:,2:this%n-1,:) )  

                RHS(:,1,:)       = f(:,1,:)
                RHS(:,this%n,:) = f(:,this%n,:)

            end select
        end select
    
    end subroutine

    !pure subroutine ComputeZIRHS(this, f, RHS, n1, n2, bc1, bcn)
    subroutine ComputeZIRHS(this, f, RHS, dir, n1, n2, bc1, bcn)
    
        class( ei02 ), intent(in) :: this
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
        integer :: j,i

        select case (this%periodic)
        case (.TRUE.)
            a10 = aI02
            do j=1,n2
                do i=1,n1
                    select case(dir)
                    case("N2F")
                        RHS(i,j,1:this%n-1) = a10 * ( f(i,j,1:this%n-1) + f(i,j,2:this%n) )
                        RHS(i,j,this%n)     = a10 * ( f(i,j,this%n)     + f(i,j,1) )
                    case("F2N")
                        RHS(i,j,1)          = a10 * ( f(i,j,this%n)     + f(i,j,1) )
                        RHS(i,j,2:this%n)   = a10 * ( f(i,j,1:this%n-1) + f(i,j,2:this%n) )
                    end select
                end do
            end do

        case (.FALSE.)
            a10    = q_hat 
            select case(dir)
            case("N2F")
                !2nd order    
                RHS(:,:,1:this%n-1) = a10 * ( f(:,:,1:this%n-1) + f(:,:,2:this%n) )  
            case("F2N")
                !2nd order interior, 1st order at boundary
                RHS(:,:,2:this%n-1) = a10 * ( f(:,:,1:this%n-2) + f(:,:,2:this%n-1) )  

                RHS(:,:,1)       = f(:,:,1)
                RHS(:,:,this%n) = f(:,:,this%n)

            end select
            
        end select
    
    end subroutine
    
    
    subroutine iN2F1(this, fN, fF, na, nb, bc1_, bcn_)
        class( ei02 ), intent(in) :: this
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

        call this%ComputeXIRHS(fN, fF, "N2F", na, nb, bc1, bcn)

    end subroutine

    subroutine iN2F2(this, fN, fF, na, nb, bc1_, bcn_)
        class( ei02 ), intent(in) :: this
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

        call this%ComputeYIRHS(fN, fF, "N2F", na, nb, bc1, bcn)

    end subroutine

    subroutine iN2F3(this, fN, fF, na, nb, bc1_, bcn_)
        class( ei02 ), intent(in) :: this
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

        call this%ComputeZIRHS(fN, fF, "N2F", na, nb, bc1, bcn)

    end subroutine

    subroutine iF2N1(this, fF, fN, na, nb, bc1_, bcn_)
        class( ei02 ), intent(in) :: this
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

        call this%ComputeXIRHS(fF, fN, "F2N", na, nb, bc1, bcn)

    end subroutine

    subroutine iF2N2(this, fF, fN, na, nb, bc1_, bcn_)
        class( ei02 ), intent(in) :: this
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

        call this%ComputeYIRHS(fF, fN, "F2N", na, nb, bc1, bcn)

    end subroutine

    subroutine iF2N3(this, fF, fN, na, nb, bc1_, bcn_)
        class( ei02 ), intent(in) :: this
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

        call this%ComputeZIRHS(fF, fN, "F2N", na, nb, bc1, bcn)

    end subroutine

end module
