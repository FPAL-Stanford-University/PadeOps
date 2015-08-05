! Routines specific to the Least Squares fiter

module lstsqstuff

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two
    use exits,           only: GracefulExit

    implicit none

    private
    public :: lstsq, als, bls, cls, dls, els
    
    ! Gaussian filter of width 4 \Delta
    real(rkind), parameter :: als    = real( 0.5      , rkind)
    real(rkind), parameter :: bls    = real( 0.6744132, rkind)/two 
    real(rkind), parameter :: cls    = real( 0.0      , rkind)/two 
    real(rkind), parameter :: dls    = real(-0.1744132, rkind)/two 
    real(rkind), parameter :: els    = real( 0.0      , rkind)/two 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! NOTE : The following variables are used for non-periodic filter !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
    ! First point
    real(rkind), parameter :: b1_als    = real(   5, rkind)/real(     6, rkind)
    real(rkind), parameter :: b1_bls    = real(   1, rkind)/real(     6, rkind) 

    ! Second point
    real(rkind), parameter :: b2_als    = real(   2, rkind)/real(     3, rkind)
    real(rkind), parameter :: b2_bls    = real(   1, rkind)/real(     6, rkind) 
    
    ! Third point
    real(rkind), parameter :: b3_als    = real(  31, rkind)/real(    64, rkind)
    real(rkind), parameter :: b3_bls    = real(   7, rkind)/real(    32, rkind) 
    real(rkind), parameter :: b3_cls    = real(   5, rkind)/real(   128, rkind) 
   
    ! Fourth point
    real(rkind), parameter :: b4_als    = real(  17, rkind)/real(    48, rkind)
    real(rkind), parameter :: b4_bls    = real(  15, rkind)/real(    64, rkind) 
    real(rkind), parameter :: b4_cls    = real(   7, rkind)/real(    96, rkind) 
    real(rkind), parameter :: b4_dls    = real(   1, rkind)/real(    64, rkind) 
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    type lstsq
        
        private
        
        integer     :: n

        logical     :: periodic=.TRUE.
        logical     :: initialized=.FALSE.

        contains

        procedure :: init
        procedure :: destroy

        procedure :: filter1
        procedure :: filter2
        procedure :: filter3
        
    end type



contains

    function init(this, n_, periodic_) result(ierr)
   
        class( lstsq ), intent(inout) :: this
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

        class( lstsq ), intent(inout) :: this

        this%initialized = .FALSE.
        this%periodic = .TRUE.

    end subroutine

    pure subroutine filter1(this, f, fil, n2, n3)
    
        class( lstsq ), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: fil
        integer :: j,k

        if(this%n == 1) then
            fil = f
            return
        end if
    
        select case (this%periodic)
        case (.TRUE.)
            do k=1,n3
                do j=1,n2
                    fil(         1,j,k) = als * ( f(         1,j,k) )                     &
                                        + bls * ( f(         2,j,k) + f(    this%n,j,k) ) &
                                        + cls * ( f(         3,j,k) + f(  this%n-1,j,k) ) &
                                        + dls * ( f(         4,j,k) + f(  this%n-2,j,k) ) &
                                        + els * ( f(         5,j,k) + f(  this%n-3,j,k) )
                    fil(         2,j,k) = als * ( f(         2,j,k) )                     &
                                        + bls * ( f(         3,j,k) + f(         1,j,k) ) &
                                        + cls * ( f(         4,j,k) + f(    this%n,j,k) ) &
                                        + dls * ( f(         5,j,k) + f(  this%n-1,j,k) ) &
                                        + els * ( f(         6,j,k) + f(  this%n-2,j,k) )
                    fil(         3,j,k) = als * ( f(         3,j,k) )                     &
                                        + bls * ( f(         4,j,k) + f(         2,j,k) ) &
                                        + cls * ( f(         5,j,k) + f(         1,j,k) ) &
                                        + dls * ( f(         6,j,k) + f(    this%n,j,k) ) &
                                        + els * ( f(         7,j,k) + f(  this%n-1,j,k) )
                    fil(         4,j,k) = als * ( f(         4,j,k) )                     &
                                        + bls * ( f(         5,j,k) + f(         3,j,k) ) &
                                        + cls * ( f(         6,j,k) + f(         2,j,k) ) &
                                        + dls * ( f(         7,j,k) + f(         1,j,k) ) &
                                        + els * ( f(         8,j,k) + f(    this%n,j,k) )
                    fil(5:this%n-4,j,k) = als * ( f(5:this%n-4,j,k) )                     &
                                        + bls * ( f(6:this%n-3,j,k) + f(4:this%n-5,j,k) ) &
                                        + cls * ( f(7:this%n-2,j,k) + f(3:this%n-6,j,k) ) &
                                        + dls * ( f(8:this%n-1,j,k) + f(2:this%n-7,j,k) ) &
                                        + els * ( f(9:this%n  ,j,k) + f(1:this%n-8,j,k) )
                    fil(  this%n-3,j,k) = als * ( f(  this%n-3,j,k) )                     &
                                        + bls * ( f(  this%n-2,j,k) + f(  this%n-4,j,k) ) &
                                        + cls * ( f(  this%n-1,j,k) + f(  this%n-5,j,k) ) &
                                        + dls * ( f(    this%n,j,k) + f(  this%n-6,j,k) ) &
                                        + els * ( f(         1,j,k) + f(  this%n-7,j,k) )
                    fil(  this%n-2,j,k) = als * ( f(  this%n-2,j,k) )                     &
                                        + bls * ( f(  this%n-1,j,k) + f(  this%n-3,j,k) ) &
                                        + cls * ( f(    this%n,j,k) + f(  this%n-4,j,k) ) &
                                        + dls * ( f(         1,j,k) + f(  this%n-5,j,k) ) &
                                        + els * ( f(         2,j,k) + f(  this%n-6,j,k) )
                    fil(  this%n-1,j,k) = als * ( f(  this%n-1,j,k) )                     &
                                        + bls * ( f(    this%n,j,k) + f(  this%n-2,j,k) ) &
                                        + cls * ( f(         1,j,k) + f(  this%n-3,j,k) ) &
                                        + dls * ( f(         2,j,k) + f(  this%n-4,j,k) ) &
                                        + els * ( f(         3,j,k) + f(  this%n-5,j,k) )
                    fil(    this%n,j,k) = als * ( f(    this%n,j,k) )                     &
                                        + bls * ( f(         1,j,k) + f(  this%n-1,j,k) ) &
                                        + cls * ( f(         2,j,k) + f(  this%n-2,j,k) ) &
                                        + dls * ( f(         3,j,k) + f(  this%n-3,j,k) ) &
                                        + els * ( f(         4,j,k) + f(  this%n-4,j,k) )
                end do
            end do

        case (.FALSE.)

            do k = 1,n3
                do j = 1,n2
                    
                    fil(         1,j,k) = b1_als * ( f(         1,j,k) )                     &
                                        + b1_bls * ( f(         2,j,k) ) 

                    fil(         2,j,k) = b2_als * ( f(         2,j,k) )                     &
                                        + b2_bls * ( f(         3,j,k) + f(         1,j,k) ) 
                    
                    fil(         3,j,k) = b3_als * ( f(         3,j,k) )                     &
                                        + b3_bls * ( f(         4,j,k) + f(         2,j,k) ) &
                                        + b3_cls * ( f(         5,j,k) + f(         1,j,k) )

                    fil(         4,j,k) = b4_als * ( f(         4,j,k) )                     &
                                        + b4_bls * ( f(         5,j,k) + f(         3,j,k) ) &
                                        + b4_cls * ( f(         6,j,k) + f(         2,j,k) ) &
                                        + b4_dls * ( f(         7,j,k) + f(         1,j,k) ) 

                    fil(5:this%n-4,j,k) =    als * ( f(5:this%n-4,j,k) )                     &
                                        +    bls * ( f(6:this%n-3,j,k) + f(4:this%n-5,j,k) ) &
                                        +    cls * ( f(7:this%n-2,j,k) + f(3:this%n-6,j,k) ) &
                                        +    dls * ( f(8:this%n-1,j,k) + f(2:this%n-7,j,k) ) &
                                        +    els * ( f(9:this%n  ,j,k) + f(1:this%n-8,j,k) )

                    fil(  this%n-3,j,k) = b4_als * ( f(  this%n-3,j,k) )                     &
                                        + b4_bls * ( f(  this%n-2,j,k) + f(  this%n-4,j,k) ) &
                                        + b4_cls * ( f(  this%n-1,j,k) + f(  this%n-5,j,k) ) &
                                        + b4_dls * ( f(    this%n,j,k) + f(  this%n-6,j,k) ) 

                    fil(  this%n-2,j,k) = b3_als * ( f(  this%n-2,j,k) )                     &
                                        + b3_bls * ( f(  this%n-1,j,k) + f(  this%n-3,j,k) ) &
                                        + b3_cls * ( f(    this%n,j,k) + f(  this%n-4,j,k) ) 

                    fil(  this%n-1,j,k) = b2_als * ( f(  this%n-1,j,k) )                     &
                                        + b2_bls * ( f(    this%n,j,k) + f(  this%n-2,j,k) ) 

                    fil(    this%n,j,k) = b1_als * ( f(    this%n,j,k) )                     &
                                        + b1_bls * ( f(  this%n-1,j,k) ) 
                
               end do 
            end do 
        end select
    
    end subroutine
    
    pure subroutine filter2(this, f, fil, n1, n3) 
    
        class( lstsq ), intent(in) :: this
        integer, intent(in) :: n1, n3
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: fil
        integer :: k

        if(this%n == 1) then
            fil = f
            return
        end if
    
        select case (this%periodic)
        case (.TRUE.)
            do k=1,n3
                fil(:,         1,k) = als * ( f(:,         1,k) )                     &
                                    + bls * ( f(:,         2,k) + f(:,    this%n,k) ) &
                                    + cls * ( f(:,         3,k) + f(:,  this%n-1,k) ) &
                                    + dls * ( f(:,         4,k) + f(:,  this%n-2,k) ) &
                                    + els * ( f(:,         5,k) + f(:,  this%n-3,k) )
                fil(:,         2,k) = als * ( f(:,         2,k) )                     &
                                    + bls * ( f(:,         3,k) + f(:,         1,k) ) &
                                    + cls * ( f(:,         4,k) + f(:,    this%n,k) ) &
                                    + dls * ( f(:,         5,k) + f(:,  this%n-1,k) ) &
                                    + els * ( f(:,         6,k) + f(:,  this%n-2,k) )
                fil(:,         3,k) = als * ( f(:,         3,k) )                     &
                                    + bls * ( f(:,         4,k) + f(:,         2,k) ) &
                                    + cls * ( f(:,         5,k) + f(:,         1,k) ) &
                                    + dls * ( f(:,         6,k) + f(:,    this%n,k) ) &
                                    + els * ( f(:,         7,k) + f(:,  this%n-1,k) )
                fil(:,         4,k) = als * ( f(:,         4,k) )                     &
                                    + bls * ( f(:,         5,k) + f(:,         3,k) ) &
                                    + cls * ( f(:,         6,k) + f(:,         2,k) ) &
                                    + dls * ( f(:,         7,k) + f(:,         1,k) ) &
                                    + els * ( f(:,         8,k) + f(:,    this%n,k) )
                fil(:,5:this%n-4,k) = als * ( f(:,5:this%n-4,k) )                     &
                                    + bls * ( f(:,6:this%n-3,k) + f(:,4:this%n-5,k) ) &
                                    + cls * ( f(:,7:this%n-2,k) + f(:,3:this%n-6,k) ) &
                                    + dls * ( f(:,8:this%n-1,k) + f(:,2:this%n-7,k) ) &
                                    + els * ( f(:,9:this%n  ,k) + f(:,1:this%n-8,k) )
                fil(:,  this%n-3,k) = als * ( f(:,  this%n-3,k) )                     &
                                    + bls * ( f(:,  this%n-2,k) + f(:,  this%n-4,k) ) &
                                    + cls * ( f(:,  this%n-1,k) + f(:,  this%n-5,k) ) &
                                    + dls * ( f(:,    this%n,k) + f(:,  this%n-6,k) ) &
                                    + els * ( f(:,         1,k) + f(:,  this%n-7,k) )
                fil(:,  this%n-2,k) = als * ( f(:,  this%n-2,k) )                     &
                                    + bls * ( f(:,  this%n-1,k) + f(:,  this%n-3,k) ) &
                                    + cls * ( f(:,    this%n,k) + f(:,  this%n-4,k) ) &
                                    + dls * ( f(:,         1,k) + f(:,  this%n-5,k) ) &
                                    + els * ( f(:,         2,k) + f(:,  this%n-6,k) )
                fil(:,  this%n-1,k) = als * ( f(:,  this%n-1,k) )                     &
                                    + bls * ( f(:,    this%n,k) + f(:,  this%n-2,k) ) &
                                    + cls * ( f(:,         1,k) + f(:,  this%n-3,k) ) &
                                    + dls * ( f(:,         2,k) + f(:,  this%n-4,k) ) &
                                    + els * ( f(:,         3,k) + f(:,  this%n-5,k) )
                fil(:,    this%n,k) = als * ( f(:,    this%n,k) )                     &
                                    + bls * ( f(:,         1,k) + f(:,  this%n-1,k) ) &
                                    + cls * ( f(:,         2,k) + f(:,  this%n-2,k) ) &
                                    + dls * ( f(:,         3,k) + f(:,  this%n-3,k) ) &
                                    + els * ( f(:,         4,k) + f(:,  this%n-4,k) )
            end do
        case (.FALSE.)

            do k = 1,n3
                fil(:,         1,k) = b1_als * ( f(:,         1,k) )                     &
                                    + b1_bls * ( f(:,         2,k) ) 

                fil(:,         2,k) = b2_als * ( f(:,         2,k) )                     &
                                    + b2_bls * ( f(:,         3,k) + f(:,         1,k) ) 
                
                fil(:,         3,k) = b3_als * ( f(:,         3,k) )                     &
                                    + b3_bls * ( f(:,         4,k) + f(:,         2,k) ) &
                                    + b3_cls * ( f(:,         5,k) + f(:,         1,k) )

                fil(:,         4,k) = b4_als * ( f(:,         4,k) )                     &
                                    + b4_bls * ( f(:,         5,k) + f(:,         3,k) ) &
                                    + b4_cls * ( f(:,         6,k) + f(:,         2,k) ) &
                                    + b4_dls * ( f(:,         7,k) + f(:,         1,k) ) 

                fil(:,5:this%n-4,k) =    als * ( f(:,5:this%n-4,k) )                     &
                                    +    bls * ( f(:,6:this%n-3,k) + f(:,4:this%n-5,k) ) &
                                    +    cls * ( f(:,7:this%n-2,k) + f(:,3:this%n-6,k) ) &
                                    +    dls * ( f(:,8:this%n-1,k) + f(:,2:this%n-7,k) ) &
                                    +    els * ( f(:,9:this%n  ,k) + f(:,1:this%n-8,k) )

                fil(:,  this%n-3,k) = b4_als * ( f(:,  this%n-3,k) )                     &
                                    + b4_bls * ( f(:,  this%n-2,k) + f(:,  this%n-4,k) ) &
                                    + b4_cls * ( f(:,  this%n-1,k) + f(:,  this%n-5,k) ) &
                                    + b4_dls * ( f(:,    this%n,k) + f(:,  this%n-6,k) ) 

                fil(:,  this%n-2,k) = b3_als * ( f(:,  this%n-2,k) )                     &
                                    + b3_bls * ( f(:,  this%n-1,k) + f(:,  this%n-3,k) ) &
                                    + b3_cls * ( f(:,    this%n,k) + f(:,  this%n-4,k) ) 

                fil(:,  this%n-1,k) = b2_als * ( f(:,  this%n-1,k) )                     &
                                    + b2_bls * ( f(:,    this%n,k) + f(:,  this%n-2,k) ) 

                fil(:,    this%n,k) = b1_als * ( f(:,    this%n,k) )                     &
                                    + b1_bls * ( f(:,  this%n-1,k) ) 
                
            end do 
        end select
    
    end subroutine

    pure subroutine filter3(this, f, fil, n1, n2)
    
        class( lstsq ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in) :: f
        real(rkind), dimension(n1,n2,this%n), intent(out) :: fil
    
        if(this%n == 1) then
            fil = f
            return
        end if
    
        select case (this%periodic)
        case (.TRUE.)
                fil(:,:,         1) = als * ( f(:,:,         1) )                     &
                                    + bls * ( f(:,:,         2) + f(:,:,    this%n) ) &
                                    + cls * ( f(:,:,         3) + f(:,:,  this%n-1) ) &
                                    + dls * ( f(:,:,         4) + f(:,:,  this%n-2) ) &
                                    + els * ( f(:,:,         5) + f(:,:,  this%n-3) )
                fil(:,:,         2) = als * ( f(:,:,         2) )                     &
                                    + bls * ( f(:,:,         3) + f(:,:,         1) ) &
                                    + cls * ( f(:,:,         4) + f(:,:,    this%n) ) &
                                    + dls * ( f(:,:,         5) + f(:,:,  this%n-1) ) &
                                    + els * ( f(:,:,         6) + f(:,:,  this%n-2) )
                fil(:,:,         3) = als * ( f(:,:,         3) )                     &
                                    + bls * ( f(:,:,         4) + f(:,:,         2) ) &
                                    + cls * ( f(:,:,         5) + f(:,:,         1) ) &
                                    + dls * ( f(:,:,         6) + f(:,:,    this%n) ) &
                                    + els * ( f(:,:,         7) + f(:,:,  this%n-1) )
                fil(:,:,         4) = als * ( f(:,:,         4) )                     &
                                    + bls * ( f(:,:,         5) + f(:,:,         3) ) &
                                    + cls * ( f(:,:,         6) + f(:,:,         2) ) &
                                    + dls * ( f(:,:,         7) + f(:,:,         1) ) &
                                    + els * ( f(:,:,         8) + f(:,:,    this%n) )
                fil(:,:,5:this%n-4) = als * ( f(:,:,5:this%n-4) )                     &
                                    + bls * ( f(:,:,6:this%n-3) + f(:,:,4:this%n-5) ) &
                                    + cls * ( f(:,:,7:this%n-2) + f(:,:,3:this%n-6) ) &
                                    + dls * ( f(:,:,8:this%n-1) + f(:,:,2:this%n-7) ) &
                                    + els * ( f(:,:,9:this%n  ) + f(:,:,1:this%n-8) )
                fil(:,:,  this%n-3) = als * ( f(:,:,  this%n-3) )                     &
                                    + bls * ( f(:,:,  this%n-2) + f(:,:,  this%n-4) ) &
                                    + cls * ( f(:,:,  this%n-1) + f(:,:,  this%n-5) ) &
                                    + dls * ( f(:,:,    this%n) + f(:,:,  this%n-6) ) &
                                    + els * ( f(:,:,         1) + f(:,:,  this%n-7) )
                fil(:,:,  this%n-2) = als * ( f(:,:,  this%n-2) )                     &
                                    + bls * ( f(:,:,  this%n-1) + f(:,:,  this%n-3) ) &
                                    + cls * ( f(:,:,    this%n) + f(:,:,  this%n-4) ) &
                                    + dls * ( f(:,:,         1) + f(:,:,  this%n-5) ) &
                                    + els * ( f(:,:,         2) + f(:,:,  this%n-6) )
                fil(:,:,  this%n-1) = als * ( f(:,:,  this%n-1) )                     &
                                    + bls * ( f(:,:,    this%n) + f(:,:,  this%n-2) ) &
                                    + cls * ( f(:,:,         1) + f(:,:,  this%n-3) ) &
                                    + dls * ( f(:,:,         2) + f(:,:,  this%n-4) ) &
                                    + els * ( f(:,:,         3) + f(:,:,  this%n-5) )
                fil(:,:,    this%n) = als * ( f(:,:,    this%n) )                     &
                                    + bls * ( f(:,:,         1) + f(:,:,  this%n-1) ) &
                                    + cls * ( f(:,:,         2) + f(:,:,  this%n-2) ) &
                                    + dls * ( f(:,:,         3) + f(:,:,  this%n-3) ) &
                                    + els * ( f(:,:,         4) + f(:,:,  this%n-4) )
        case (.FALSE.)
                    
            fil(:,:,         1) = b1_als * ( f(:,:,         1) )                     &
                                + b1_bls * ( f(:,:,         2) ) 

            fil(:,:,         2) = b2_als * ( f(:,:,         2) )                     &
                                + b2_bls * ( f(:,:,         3) + f(:,:,         1) ) 
            
            fil(:,:,         3) = b3_als * ( f(:,:,         3) )                     &
                                + b3_bls * ( f(:,:,         4) + f(:,:,         2) ) &
                                + b3_cls * ( f(:,:,         5) + f(:,:,         1) )

            fil(:,:,         4) = b4_als * ( f(:,:,         4) )                     &
                                + b4_bls * ( f(:,:,         5) + f(:,:,         3) ) &
                                + b4_cls * ( f(:,:,         6) + f(:,:,         2) ) &
                                + b4_dls * ( f(:,:,         7) + f(:,:,         1) ) 

            fil(:,:,5:this%n-4) =    als * ( f(:,:,5:this%n-4) )                     &
                                +    bls * ( f(:,:,6:this%n-3) + f(:,:,4:this%n-5) ) &
                                +    cls * ( f(:,:,7:this%n-2) + f(:,:,3:this%n-6) ) &
                                +    dls * ( f(:,:,8:this%n-1) + f(:,:,2:this%n-7) ) &
                                +    els * ( f(:,:,9:this%n  ) + f(:,:,1:this%n-8) )

            fil(:,:,  this%n-3) = b4_als * ( f(:,:,  this%n-3) )                     &
                                + b4_bls * ( f(:,:,  this%n-2) + f(:,:,  this%n-4) ) &
                                + b4_cls * ( f(:,:,  this%n-1) + f(:,:,  this%n-5) ) &
                                + b4_dls * ( f(:,:,    this%n) + f(:,:,  this%n-6) ) 

            fil(:,:,  this%n-2) = b3_als * ( f(:,:,  this%n-2) )                     &
                                + b3_bls * ( f(:,:,  this%n-1) + f(:,:,  this%n-3) ) &
                                + b3_cls * ( f(:,:,    this%n) + f(:,:,  this%n-4) ) 

            fil(:,:,  this%n-1) = b2_als * ( f(:,:,  this%n-1) )                     &
                                + b2_bls * ( f(:,:,    this%n) + f(:,:,  this%n-2) ) 

            fil(:,:,    this%n) = b1_als * ( f(:,:,    this%n) )                     &
                                + b1_bls * ( f(:,:,  this%n-1) ) 

        end select
    
    end subroutine
    
end module
