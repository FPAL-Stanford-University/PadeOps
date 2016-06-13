! Routines specific to 8th order Compact Filter
! Periodic LU based on  Neossi Nguetchue, Abelman (Appl. Math. & Comp. 2008)

module gaussianstuff

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two
    use exits,           only: GracefulExit

    implicit none

    private
    public :: gaussian, agf, bgf, cgf, dgf, egf
    
    ! Gaussian filter of width 4 \Delta
    real(rkind), parameter :: agf    = real(3565, rkind)/real( 10368, rkind) 
    real(rkind), parameter :: bgf    = real(3091, rkind)/real( 12960, rkind) 
    real(rkind), parameter :: cgf    = real(1997, rkind)/real( 25920, rkind) 
    real(rkind), parameter :: dgf    = real( 149, rkind)/real( 12960, rkind) 
    real(rkind), parameter :: egf    = real( 107, rkind)/real(103680, rkind) 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic filter !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
    ! First point
    real(rkind), parameter :: b1_agf    = real(   5, rkind)/real(     6, rkind)
    real(rkind), parameter :: b1_bgf    = real(   1, rkind)/real(     6, rkind) 

    ! Second point
    real(rkind), parameter :: b2_agf    = real(   2, rkind)/real(     3, rkind)
    real(rkind), parameter :: b2_bgf    = real(   1, rkind)/real(     6, rkind) 
    
    ! Third point
    real(rkind), parameter :: b3_agf    = real(  31, rkind)/real(    64, rkind)
    real(rkind), parameter :: b3_bgf    = real(   7, rkind)/real(    32, rkind) 
    real(rkind), parameter :: b3_cgf    = real(   5, rkind)/real(   128, rkind) 
   
    ! Fourth point
    real(rkind), parameter :: b4_agf    = real(  17, rkind)/real(    48, rkind)
    real(rkind), parameter :: b4_bgf    = real(  15, rkind)/real(    64, rkind) 
    real(rkind), parameter :: b4_cgf    = real(   7, rkind)/real(    96, rkind) 
    real(rkind), parameter :: b4_dgf    = real(   1, rkind)/real(    64, rkind) 
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    type gaussian
        
        private
        
        integer     :: n

        logical     :: periodic=.TRUE.
        logical     :: initialized=.FALSE.

        contains

        procedure :: init
        procedure :: destroy

        procedure :: filter1
        procedure :: filter2
        procedure, private :: filter3_real
        procedure, private :: filter3_cmplx
        generic   :: filter3 => filter3_real, filter3_cmplx
    end type



contains

    function init(this, n_, periodic_) result(ierr)
   
        class( gaussian ), intent(inout) :: this
        integer, intent(in) :: n_
        logical, intent(in) :: periodic_
        integer :: ierr
        
        if (this%initialized) then
            call this%destroy()
        end if
        this%initialized = .TRUE.

        this%n = n_

        this%periodic = periodic_

        ! If everything passes
        ierr = 0
    
    end function
    
    subroutine destroy(this)

        class( gaussian ), intent(inout) :: this

        this%initialized = .FALSE.
        this%periodic = .TRUE.

    end subroutine

    subroutine filter1(this, f, fil, n2, n3, bc1_, bcn_)
    
        class( gaussian ), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: fil
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        integer :: j,k

        if(this%n == 1) then
            fil = f
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

        select case (this%periodic)
        case (.TRUE.)
            do k=1,n3
                do j=1,n2
                    fil(         1,j,k) = agf * ( f(         1,j,k) )                     &
                                        + bgf * ( f(         2,j,k) + f(    this%n,j,k) ) &
                                        + cgf * ( f(         3,j,k) + f(  this%n-1,j,k) ) &
                                        + dgf * ( f(         4,j,k) + f(  this%n-2,j,k) ) &
                                        + egf * ( f(         5,j,k) + f(  this%n-3,j,k) )
                    fil(         2,j,k) = agf * ( f(         2,j,k) )                     &
                                        + bgf * ( f(         3,j,k) + f(         1,j,k) ) &
                                        + cgf * ( f(         4,j,k) + f(    this%n,j,k) ) &
                                        + dgf * ( f(         5,j,k) + f(  this%n-1,j,k) ) &
                                        + egf * ( f(         6,j,k) + f(  this%n-2,j,k) )
                    fil(         3,j,k) = agf * ( f(         3,j,k) )                     &
                                        + bgf * ( f(         4,j,k) + f(         2,j,k) ) &
                                        + cgf * ( f(         5,j,k) + f(         1,j,k) ) &
                                        + dgf * ( f(         6,j,k) + f(    this%n,j,k) ) &
                                        + egf * ( f(         7,j,k) + f(  this%n-1,j,k) )
                    fil(         4,j,k) = agf * ( f(         4,j,k) )                     &
                                        + bgf * ( f(         5,j,k) + f(         3,j,k) ) &
                                        + cgf * ( f(         6,j,k) + f(         2,j,k) ) &
                                        + dgf * ( f(         7,j,k) + f(         1,j,k) ) &
                                        + egf * ( f(         8,j,k) + f(    this%n,j,k) )
                    fil(5:this%n-4,j,k) = agf * ( f(5:this%n-4,j,k) )                     &
                                        + bgf * ( f(6:this%n-3,j,k) + f(4:this%n-5,j,k) ) &
                                        + cgf * ( f(7:this%n-2,j,k) + f(3:this%n-6,j,k) ) &
                                        + dgf * ( f(8:this%n-1,j,k) + f(2:this%n-7,j,k) ) &
                                        + egf * ( f(9:this%n  ,j,k) + f(1:this%n-8,j,k) )
                    fil(  this%n-3,j,k) = agf * ( f(  this%n-3,j,k) )                     &
                                        + bgf * ( f(  this%n-2,j,k) + f(  this%n-4,j,k) ) &
                                        + cgf * ( f(  this%n-1,j,k) + f(  this%n-5,j,k) ) &
                                        + dgf * ( f(    this%n,j,k) + f(  this%n-6,j,k) ) &
                                        + egf * ( f(         1,j,k) + f(  this%n-7,j,k) )
                    fil(  this%n-2,j,k) = agf * ( f(  this%n-2,j,k) )                     &
                                        + bgf * ( f(  this%n-1,j,k) + f(  this%n-3,j,k) ) &
                                        + cgf * ( f(    this%n,j,k) + f(  this%n-4,j,k) ) &
                                        + dgf * ( f(         1,j,k) + f(  this%n-5,j,k) ) &
                                        + egf * ( f(         2,j,k) + f(  this%n-6,j,k) )
                    fil(  this%n-1,j,k) = agf * ( f(  this%n-1,j,k) )                     &
                                        + bgf * ( f(    this%n,j,k) + f(  this%n-2,j,k) ) &
                                        + cgf * ( f(         1,j,k) + f(  this%n-3,j,k) ) &
                                        + dgf * ( f(         2,j,k) + f(  this%n-4,j,k) ) &
                                        + egf * ( f(         3,j,k) + f(  this%n-5,j,k) )
                    fil(    this%n,j,k) = agf * ( f(    this%n,j,k) )                     &
                                        + bgf * ( f(         1,j,k) + f(  this%n-1,j,k) ) &
                                        + cgf * ( f(         2,j,k) + f(  this%n-2,j,k) ) &
                                        + dgf * ( f(         3,j,k) + f(  this%n-3,j,k) ) &
                                        + egf * ( f(         4,j,k) + f(  this%n-4,j,k) )
                end do
            end do

        case (.FALSE.)

            do k = 1,n3
                do j = 1,n2
                    select case(bc1)
                    case(0)
                        fil(         1,j,k) = b1_agf * ( f(         1,j,k) )                     &
                                            + b1_bgf * ( f(         2,j,k) ) 

                        fil(         2,j,k) = b2_agf * ( f(         2,j,k) )                     &
                                            + b2_bgf * ( f(         3,j,k) + f(         1,j,k) ) 
                        
                        fil(         3,j,k) = b3_agf * ( f(         3,j,k) )                     &
                                            + b3_bgf * ( f(         4,j,k) + f(         2,j,k) ) &
                                            + b3_cgf * ( f(         5,j,k) + f(         1,j,k) )

                        fil(         4,j,k) = b4_agf * ( f(         4,j,k) )                     &
                                            + b4_bgf * ( f(         5,j,k) + f(         3,j,k) ) &
                                            + b4_cgf * ( f(         6,j,k) + f(         2,j,k) ) &
                                            + b4_dgf * ( f(         7,j,k) + f(         1,j,k) ) 
                    case(1)
                        fil(1,j,k) =    agf * ( f(1,j,k) )            &
                                   +    bgf * ( f(2,j,k) + f(2,j,k) ) &
                                   +    cgf * ( f(3,j,k) + f(3,j,k) ) &
                                   +    dgf * ( f(4,j,k) + f(4,j,k) ) &
                                   +    egf * ( f(5,j,k) + f(5,j,k) )

                        fil(2,j,k) =    agf * ( f(2,j,k) )            &
                                   +    bgf * ( f(3,j,k) + f(1,j,k) ) &
                                   +    cgf * ( f(4,j,k) + f(2,j,k) ) &
                                   +    dgf * ( f(5,j,k) + f(3,j,k) ) &
                                   +    egf * ( f(6,j,k) + f(4,j,k) )

                        fil(3,j,k) =    agf * ( f(3,j,k) )            &
                                   +    bgf * ( f(4,j,k) + f(2,j,k) ) &
                                   +    cgf * ( f(5,j,k) + f(1,j,k) ) &
                                   +    dgf * ( f(6,j,k) + f(2,j,k) ) &
                                   +    egf * ( f(7,j,k) + f(3,j,k) )

                        fil(4,j,k) =    agf * ( f(4,j,k) )            &
                                   +    bgf * ( f(5,j,k) + f(3,j,k) ) &
                                   +    cgf * ( f(6,j,k) + f(2,j,k) ) &
                                   +    dgf * ( f(7,j,k) + f(1,j,k) ) &
                                   +    egf * ( f(8,j,k) + f(2,j,k) )
                    case(-1)
                        fil(1,j,k) =    agf * ( f(1,j,k) )            &
                                   +    bgf * ( f(2,j,k) - f(2,j,k) ) &
                                   +    cgf * ( f(3,j,k) - f(3,j,k) ) &
                                   +    dgf * ( f(4,j,k) - f(4,j,k) ) &
                                   +    egf * ( f(5,j,k) - f(5,j,k) )

                        fil(2,j,k) =    agf * ( f(2,j,k) )            &
                                   +    bgf * ( f(3,j,k) + f(1,j,k) ) &
                                   +    cgf * ( f(4,j,k) - f(2,j,k) ) &
                                   +    dgf * ( f(5,j,k) - f(3,j,k) ) &
                                   +    egf * ( f(6,j,k) - f(4,j,k) )

                        fil(3,j,k) =    agf * ( f(3,j,k) )            &
                                   +    bgf * ( f(4,j,k) + f(2,j,k) ) &
                                   +    cgf * ( f(5,j,k) + f(1,j,k) ) &
                                   +    dgf * ( f(6,j,k) - f(2,j,k) ) &
                                   +    egf * ( f(7,j,k) - f(3,j,k) )

                        fil(4,j,k) =    agf * ( f(4,j,k) )            &
                                   +    bgf * ( f(5,j,k) + f(3,j,k) ) &
                                   +    cgf * ( f(6,j,k) + f(2,j,k) ) &
                                   +    dgf * ( f(7,j,k) + f(1,j,k) ) &
                                   +    egf * ( f(8,j,k) - f(2,j,k) )
                    end select

                    fil(5:this%n-4,j,k) =    agf * ( f(5:this%n-4,j,k) )                     &
                                        +    bgf * ( f(6:this%n-3,j,k) + f(4:this%n-5,j,k) ) &
                                        +    cgf * ( f(7:this%n-2,j,k) + f(3:this%n-6,j,k) ) &
                                        +    dgf * ( f(8:this%n-1,j,k) + f(2:this%n-7,j,k) ) &
                                        +    egf * ( f(9:this%n  ,j,k) + f(1:this%n-8,j,k) )

                    select case(bcn)
                    case(0)
                        fil(  this%n-3,j,k) = b4_agf * ( f(  this%n-3,j,k) )                     &
                                            + b4_bgf * ( f(  this%n-2,j,k) + f(  this%n-4,j,k) ) &
                                            + b4_cgf * ( f(  this%n-1,j,k) + f(  this%n-5,j,k) ) &
                                            + b4_dgf * ( f(    this%n,j,k) + f(  this%n-6,j,k) ) 

                        fil(  this%n-2,j,k) = b3_agf * ( f(  this%n-2,j,k) )                     &
                                            + b3_bgf * ( f(  this%n-1,j,k) + f(  this%n-3,j,k) ) &
                                            + b3_cgf * ( f(    this%n,j,k) + f(  this%n-4,j,k) ) 

                        fil(  this%n-1,j,k) = b2_agf * ( f(  this%n-1,j,k) )                     &
                                            + b2_bgf * ( f(    this%n,j,k) + f(  this%n-2,j,k) ) 

                        fil(    this%n,j,k) = b1_agf * ( f(    this%n,j,k) )                     &
                                            + b1_bgf * ( f(  this%n-1,j,k) ) 
                    case(1)
                        fil(this%n-3,j,k) =    agf * ( f(this%n-3,j,k) )                   &
                                          +    bgf * ( f(this%n-2,j,k) + f(this%n-4,j,k) ) &
                                          +    cgf * ( f(this%n-1,j,k) + f(this%n-5,j,k) ) &
                                          +    dgf * ( f(this%n  ,j,k) + f(this%n-6,j,k) ) &
                                          +    egf * ( f(this%n-1,j,k) + f(this%n-7,j,k) )

                        fil(this%n-2,j,k) =    agf * ( f(this%n-2,j,k) )                   &
                                          +    bgf * ( f(this%n-1,j,k) + f(this%n-3,j,k) ) &
                                          +    cgf * ( f(this%n  ,j,k) + f(this%n-4,j,k) ) &
                                          +    dgf * ( f(this%n-1,j,k) + f(this%n-5,j,k) ) &
                                          +    egf * ( f(this%n-2,j,k) + f(this%n-6,j,k) )

                        fil(this%n-1,j,k) =    agf * ( f(this%n-1,j,k) )                   &
                                          +    bgf * ( f(this%n  ,j,k) + f(this%n-2,j,k) ) &
                                          +    cgf * ( f(this%n-1,j,k) + f(this%n-3,j,k) ) &
                                          +    dgf * ( f(this%n-2,j,k) + f(this%n-4,j,k) ) &
                                          +    egf * ( f(this%n-3,j,k) + f(this%n-5,j,k) )

                        fil(this%n  ,j,k) =    agf * ( f(this%n  ,j,k) )                   &
                                          +    bgf * ( f(this%n-1,j,k) + f(this%n-1,j,k) ) &
                                          +    cgf * ( f(this%n-2,j,k) + f(this%n-2,j,k) ) &
                                          +    dgf * ( f(this%n-3,j,k) + f(this%n-3,j,k) ) &
                                          +    egf * ( f(this%n-4,j,k) + f(this%n-4,j,k) )
                    case(-1)
                        fil(this%n-3,j,k) =    agf * ( f(this%n-3,j,k) )                   &
                                          +    bgf * ( f(this%n-2,j,k) + f(this%n-4,j,k) ) &
                                          +    cgf * ( f(this%n-1,j,k) + f(this%n-5,j,k) ) &
                                          +    dgf * ( f(this%n  ,j,k) + f(this%n-6,j,k) ) &
                                          +    egf * (-f(this%n-1,j,k) + f(this%n-7,j,k) )

                        fil(this%n-2,j,k) =    agf * ( f(this%n-2,j,k) )                   &
                                          +    bgf * ( f(this%n-1,j,k) + f(this%n-3,j,k) ) &
                                          +    cgf * ( f(this%n  ,j,k) + f(this%n-4,j,k) ) &
                                          +    dgf * (-f(this%n-1,j,k) + f(this%n-5,j,k) ) &
                                          +    egf * (-f(this%n-2,j,k) + f(this%n-6,j,k) )

                        fil(this%n-1,j,k) =    agf * ( f(this%n-1,j,k) )                   &
                                          +    bgf * ( f(this%n  ,j,k) + f(this%n-2,j,k) ) &
                                          +    cgf * (-f(this%n-1,j,k) + f(this%n-3,j,k) ) &
                                          +    dgf * (-f(this%n-2,j,k) + f(this%n-4,j,k) ) &
                                          +    egf * (-f(this%n-3,j,k) + f(this%n-5,j,k) )

                        fil(this%n  ,j,k) =    agf * ( f(this%n  ,j,k) )                   &
                                          +    bgf * (-f(this%n-1,j,k) + f(this%n-1,j,k) ) &
                                          +    cgf * (-f(this%n-2,j,k) + f(this%n-2,j,k) ) &
                                          +    dgf * (-f(this%n-3,j,k) + f(this%n-3,j,k) ) &
                                          +    egf * (-f(this%n-4,j,k) + f(this%n-4,j,k) )
                    end select
               end do 
            end do 
        end select
    
    end subroutine
    
    subroutine filter2(this, f, fil, n1, n3, bc1_, bcn_) 
    
        class( gaussian ), intent(in) :: this
        integer, intent(in) :: n1, n3
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: fil
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
        integer :: k

        if(this%n == 1) then
            fil = f
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

        select case (this%periodic)
        case (.TRUE.)
            do k=1,n3
                fil(:,         1,k) = agf * ( f(:,         1,k) )                     &
                                    + bgf * ( f(:,         2,k) + f(:,    this%n,k) ) &
                                    + cgf * ( f(:,         3,k) + f(:,  this%n-1,k) ) &
                                    + dgf * ( f(:,         4,k) + f(:,  this%n-2,k) ) &
                                    + egf * ( f(:,         5,k) + f(:,  this%n-3,k) )
                fil(:,         2,k) = agf * ( f(:,         2,k) )                     &
                                    + bgf * ( f(:,         3,k) + f(:,         1,k) ) &
                                    + cgf * ( f(:,         4,k) + f(:,    this%n,k) ) &
                                    + dgf * ( f(:,         5,k) + f(:,  this%n-1,k) ) &
                                    + egf * ( f(:,         6,k) + f(:,  this%n-2,k) )
                fil(:,         3,k) = agf * ( f(:,         3,k) )                     &
                                    + bgf * ( f(:,         4,k) + f(:,         2,k) ) &
                                    + cgf * ( f(:,         5,k) + f(:,         1,k) ) &
                                    + dgf * ( f(:,         6,k) + f(:,    this%n,k) ) &
                                    + egf * ( f(:,         7,k) + f(:,  this%n-1,k) )
                fil(:,         4,k) = agf * ( f(:,         4,k) )                     &
                                    + bgf * ( f(:,         5,k) + f(:,         3,k) ) &
                                    + cgf * ( f(:,         6,k) + f(:,         2,k) ) &
                                    + dgf * ( f(:,         7,k) + f(:,         1,k) ) &
                                    + egf * ( f(:,         8,k) + f(:,    this%n,k) )
                fil(:,5:this%n-4,k) = agf * ( f(:,5:this%n-4,k) )                     &
                                    + bgf * ( f(:,6:this%n-3,k) + f(:,4:this%n-5,k) ) &
                                    + cgf * ( f(:,7:this%n-2,k) + f(:,3:this%n-6,k) ) &
                                    + dgf * ( f(:,8:this%n-1,k) + f(:,2:this%n-7,k) ) &
                                    + egf * ( f(:,9:this%n  ,k) + f(:,1:this%n-8,k) )
                fil(:,  this%n-3,k) = agf * ( f(:,  this%n-3,k) )                     &
                                    + bgf * ( f(:,  this%n-2,k) + f(:,  this%n-4,k) ) &
                                    + cgf * ( f(:,  this%n-1,k) + f(:,  this%n-5,k) ) &
                                    + dgf * ( f(:,    this%n,k) + f(:,  this%n-6,k) ) &
                                    + egf * ( f(:,         1,k) + f(:,  this%n-7,k) )
                fil(:,  this%n-2,k) = agf * ( f(:,  this%n-2,k) )                     &
                                    + bgf * ( f(:,  this%n-1,k) + f(:,  this%n-3,k) ) &
                                    + cgf * ( f(:,    this%n,k) + f(:,  this%n-4,k) ) &
                                    + dgf * ( f(:,         1,k) + f(:,  this%n-5,k) ) &
                                    + egf * ( f(:,         2,k) + f(:,  this%n-6,k) )
                fil(:,  this%n-1,k) = agf * ( f(:,  this%n-1,k) )                     &
                                    + bgf * ( f(:,    this%n,k) + f(:,  this%n-2,k) ) &
                                    + cgf * ( f(:,         1,k) + f(:,  this%n-3,k) ) &
                                    + dgf * ( f(:,         2,k) + f(:,  this%n-4,k) ) &
                                    + egf * ( f(:,         3,k) + f(:,  this%n-5,k) )
                fil(:,    this%n,k) = agf * ( f(:,    this%n,k) )                     &
                                    + bgf * ( f(:,         1,k) + f(:,  this%n-1,k) ) &
                                    + cgf * ( f(:,         2,k) + f(:,  this%n-2,k) ) &
                                    + dgf * ( f(:,         3,k) + f(:,  this%n-3,k) ) &
                                    + egf * ( f(:,         4,k) + f(:,  this%n-4,k) )
            end do
        case (.FALSE.)

            do k = 1,n3
                select case(bc1)
                case(0)
                    fil(:,         1,k) = b1_agf * ( f(:,         1,k) )                     &
                                        + b1_bgf * ( f(:,         2,k) ) 

                    fil(:,         2,k) = b2_agf * ( f(:,         2,k) )                     &
                                        + b2_bgf * ( f(:,         3,k) + f(:,         1,k) ) 
                    
                    fil(:,         3,k) = b3_agf * ( f(:,         3,k) )                     &
                                        + b3_bgf * ( f(:,         4,k) + f(:,         2,k) ) &
                                        + b3_cgf * ( f(:,         5,k) + f(:,         1,k) )

                    fil(:,         4,k) = b4_agf * ( f(:,         4,k) )                     &
                                        + b4_bgf * ( f(:,         5,k) + f(:,         3,k) ) &
                                        + b4_cgf * ( f(:,         6,k) + f(:,         2,k) ) &
                                        + b4_dgf * ( f(:,         7,k) + f(:,         1,k) ) 
                case(1)
                    fil(:,1,k) =    agf * ( f(:,1,k) )            &
                               +    bgf * ( f(:,2,k) + f(:,2,k) ) &
                               +    cgf * ( f(:,3,k) + f(:,3,k) ) &
                               +    dgf * ( f(:,4,k) + f(:,4,k) ) &
                               +    egf * ( f(:,5,k) + f(:,5,k) )

                    fil(:,2,k) =    agf * ( f(:,2,k) )            &
                               +    bgf * ( f(:,3,k) + f(:,1,k) ) &
                               +    cgf * ( f(:,4,k) + f(:,2,k) ) &
                               +    dgf * ( f(:,5,k) + f(:,3,k) ) &
                               +    egf * ( f(:,6,k) + f(:,4,k) )

                    fil(:,3,k) =    agf * ( f(:,3,k) )            &
                               +    bgf * ( f(:,4,k) + f(:,2,k) ) &
                               +    cgf * ( f(:,5,k) + f(:,1,k) ) &
                               +    dgf * ( f(:,6,k) + f(:,2,k) ) &
                               +    egf * ( f(:,7,k) + f(:,3,k) )

                    fil(:,4,k) =    agf * ( f(:,4,k) )            &
                               +    bgf * ( f(:,5,k) + f(:,3,k) ) &
                               +    cgf * ( f(:,6,k) + f(:,2,k) ) &
                               +    dgf * ( f(:,7,k) + f(:,1,k) ) &
                               +    egf * ( f(:,8,k) + f(:,2,k) )
                case(-1)  
                    fil(:,1,k) =    agf * ( f(:,1,k) )            &
                               +    bgf * ( f(:,2,k) - f(:,2,k) ) &
                               +    cgf * ( f(:,3,k) - f(:,3,k) ) &
                               +    dgf * ( f(:,4,k) - f(:,4,k) ) &
                               +    egf * ( f(:,5,k) - f(:,5,k) )

                    fil(:,2,k) =    agf * ( f(:,2,k) )            &
                               +    bgf * ( f(:,3,k) + f(:,1,k) ) &
                               +    cgf * ( f(:,4,k) - f(:,2,k) ) &
                               +    dgf * ( f(:,5,k) - f(:,3,k) ) &
                               +    egf * ( f(:,6,k) - f(:,4,k) )

                    fil(:,3,k) =    agf * ( f(:,3,k) )            &
                               +    bgf * ( f(:,4,k) + f(:,2,k) ) &
                               +    cgf * ( f(:,5,k) + f(:,1,k) ) &
                               +    dgf * ( f(:,6,k) - f(:,2,k) ) &
                               +    egf * ( f(:,7,k) - f(:,3,k) )

                    fil(:,4,k) =    agf * ( f(:,4,k) )            &
                               +    bgf * ( f(:,5,k) + f(:,3,k) ) &
                               +    cgf * ( f(:,6,k) + f(:,2,k) ) &
                               +    dgf * ( f(:,7,k) + f(:,1,k) ) &
                               +    egf * ( f(:,8,k) - f(:,2,k) )
                end select

                fil(:,5:this%n-4,k) =    agf * ( f(:,5:this%n-4,k) )                     &
                                    +    bgf * ( f(:,6:this%n-3,k) + f(:,4:this%n-5,k) ) &
                                    +    cgf * ( f(:,7:this%n-2,k) + f(:,3:this%n-6,k) ) &
                                    +    dgf * ( f(:,8:this%n-1,k) + f(:,2:this%n-7,k) ) &
                                    +    egf * ( f(:,9:this%n  ,k) + f(:,1:this%n-8,k) )

                select case(bcn)
                case(0)
                    fil(:,  this%n-3,k) = b4_agf * ( f(:,  this%n-3,k) )                     &
                                        + b4_bgf * ( f(:,  this%n-2,k) + f(:,  this%n-4,k) ) &
                                        + b4_cgf * ( f(:,  this%n-1,k) + f(:,  this%n-5,k) ) &
                                        + b4_dgf * ( f(:,    this%n,k) + f(:,  this%n-6,k) ) 

                    fil(:,  this%n-2,k) = b3_agf * ( f(:,  this%n-2,k) )                     &
                                        + b3_bgf * ( f(:,  this%n-1,k) + f(:,  this%n-3,k) ) &
                                        + b3_cgf * ( f(:,    this%n,k) + f(:,  this%n-4,k) ) 

                    fil(:,  this%n-1,k) = b2_agf * ( f(:,  this%n-1,k) )                     &
                                        + b2_bgf * ( f(:,    this%n,k) + f(:,  this%n-2,k) ) 

                    fil(:,    this%n,k) = b1_agf * ( f(:,    this%n,k) )                     &
                                        + b1_bgf * ( f(:,  this%n-1,k) ) 
                case(1)
                    fil(:,this%n-3,k) =    agf * ( f(:,this%n-3,k) )                   &
                                      +    bgf * ( f(:,this%n-2,k) + f(:,this%n-4,k) ) &
                                      +    cgf * ( f(:,this%n-1,k) + f(:,this%n-5,k) ) &
                                      +    dgf * ( f(:,this%n  ,k) + f(:,this%n-6,k) ) &
                                      +    egf * ( f(:,this%n-1,k) + f(:,this%n-7,k) )

                    fil(:,this%n-2,k) =    agf * ( f(:,this%n-2,k) )                   &
                                      +    bgf * ( f(:,this%n-1,k) + f(:,this%n-3,k) ) &
                                      +    cgf * ( f(:,this%n  ,k) + f(:,this%n-4,k) ) &
                                      +    dgf * ( f(:,this%n-1,k) + f(:,this%n-5,k) ) &
                                      +    egf * ( f(:,this%n-2,k) + f(:,this%n-6,k) )

                    fil(:,this%n-1,k) =    agf * ( f(:,this%n-1,k) )                   &
                                      +    bgf * ( f(:,this%n  ,k) + f(:,this%n-2,k) ) &
                                      +    cgf * ( f(:,this%n-1,k) + f(:,this%n-3,k) ) &
                                      +    dgf * ( f(:,this%n-2,k) + f(:,this%n-4,k) ) &
                                      +    egf * ( f(:,this%n-3,k) + f(:,this%n-5,k) )

                    fil(:,this%n  ,k) =    agf * ( f(:,this%n  ,k) )                   &
                                      +    bgf * ( f(:,this%n-1,k) + f(:,this%n-1,k) ) &
                                      +    cgf * ( f(:,this%n-2,k) + f(:,this%n-2,k) ) &
                                      +    dgf * ( f(:,this%n-3,k) + f(:,this%n-3,k) ) &
                                      +    egf * ( f(:,this%n-4,k) + f(:,this%n-4,k) )
                case(-1)   
                    fil(:,this%n-3,k) =    agf * ( f(:,this%n-3,k) )                   &
                                      +    bgf * ( f(:,this%n-2,k) + f(:,this%n-4,k) ) &
                                      +    cgf * ( f(:,this%n-1,k) + f(:,this%n-5,k) ) &
                                      +    dgf * ( f(:,this%n  ,k) + f(:,this%n-6,k) ) &
                                      +    egf * (-f(:,this%n-1,k) + f(:,this%n-7,k) )

                    fil(:,this%n-2,k) =    agf * ( f(:,this%n-2,k) )                   &
                                      +    bgf * ( f(:,this%n-1,k) + f(:,this%n-3,k) ) &
                                      +    cgf * ( f(:,this%n  ,k) + f(:,this%n-4,k) ) &
                                      +    dgf * (-f(:,this%n-1,k) + f(:,this%n-5,k) ) &
                                      +    egf * (-f(:,this%n-2,k) + f(:,this%n-6,k) )

                    fil(:,this%n-1,k) =    agf * ( f(:,this%n-1,k) )                   &
                                      +    bgf * ( f(:,this%n  ,k) + f(:,this%n-2,k) ) &
                                      +    cgf * (-f(:,this%n-1,k) + f(:,this%n-3,k) ) &
                                      +    dgf * (-f(:,this%n-2,k) + f(:,this%n-4,k) ) &
                                      +    egf * (-f(:,this%n-3,k) + f(:,this%n-5,k) )

                    fil(:,this%n  ,k) =    agf * ( f(:,this%n  ,k) )                   &
                                      +    bgf * (-f(:,this%n-1,k) + f(:,this%n-1,k) ) &
                                      +    cgf * (-f(:,this%n-2,k) + f(:,this%n-2,k) ) &
                                      +    dgf * (-f(:,this%n-3,k) + f(:,this%n-3,k) ) &
                                      +    egf * (-f(:,this%n-4,k) + f(:,this%n-4,k) )
                end select
                
            end do 
        end select
    
    end subroutine

    subroutine filter3_real(this, f, fil, n1, n2, bc1_, bcn_)
        class( gaussian ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in) :: f
        real(rkind), dimension(n1,n2,this%n), intent(out) :: fil
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
#include "gaussian_files/filz_common.F90" 
    end subroutine
    
    subroutine filter3_cmplx(this, f, fil, n1, n2, bc1_, bcn_)
        class( gaussian ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in) :: f
        complex(rkind), dimension(n1,n2,this%n), intent(out) :: fil
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: bc1, bcn
#include "gaussian_files/filz_common.F90" 
    end subroutine
end module
