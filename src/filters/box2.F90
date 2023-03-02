! Routines specific to the Box filter with size 2*delta where delta is the
! grid spacing

module box2stuff

    use kind_parameters, only: rkind
    use constants,       only: zero,half,fourth,one,two
    use exits,           only: GracefulExit

    implicit none

    private
    public :: box2
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  Assumes box filter with width 2*delta evaluated using Composite Trapezoidal rule  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    type box2
        
        private
        
        integer     :: n

        logical     :: periodic=.TRUE.
        logical     :: initialized=.FALSE.

        contains

        procedure :: init
        procedure :: destroy

        procedure :: filter1
        procedure :: filter2
        procedure :: filter3_real
        procedure :: filter3_cmplx
        generic :: filter3 => filter3_real, filter3_cmplx
        
    end type



contains

    function init(this, n_, periodic_) result(ierr)
   
        class( box2 ), intent(inout) :: this
        integer, intent(in) :: n_
        logical, intent(in) :: periodic_
        integer :: ierr
        
        if (this%initialized) then
            call this%destroy()
        end if
        this%initialized = .TRUE.

        this%n = n_

        this%periodic = periodic_

        if(this%n < 3) then
            ! At least 3 points are needed for box2
            ierr = 1
        else
            ! 3 or more points are present so box2 can be used
            ierr = 0
        endif
    
    end function
    
    subroutine destroy(this)

        class( box2 ), intent(inout) :: this

        this%initialized = .FALSE.
        this%periodic = .TRUE.

    end subroutine

    subroutine filter1(this, f, fil, n2, n3, bc1_, bcn_)
    
        class( box2 ), intent(in) :: this
        integer, intent(in) :: n2, n3
        real(rkind), dimension(this%n,n2,n3), intent(in) :: f
        real(rkind), dimension(this%n,n2,n3), intent(out) :: fil
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: j, k, bc1, bcn

        select case (this%periodic)
        case (.TRUE.)
          do k=1,n3
            do j=1,n2
              fil(1         ,j,k) = fourth*f(  this%n  ,j,k) + half*f(         1,j,k) + fourth*f(       2,j,k)
              fil(2:this%n-1,j,k) = fourth*f(1:this%n-2,j,k) + half*f(2:this%n-1,j,k) + fourth*f(3:this%n,j,k)
              fil(this%n    ,j,k) = fourth*f(  this%n-1,j,k) + half*f(    this%n,j,k) + fourth*f(       1,j,k)
            end do
          end do

        case (.FALSE.)

          !! Ensure BC flags are set correctly
          if (present(bc1_)) then
              bc1 = bc1_
              if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                  call GracefulExit("Incorrect X-boundary specification for bc1 (should be 0, 1 or -1)", 324)
              end if
          else
              bc1 = 0
          end if

          if (present(bcn_)) then
              bcn = bcn_
              if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                  call GracefulExit("Incorrect X-boundary specification for bcn (should be 0, 1 or -1)", 324)
              end if
          else
              bcn = 0
          end if

          select case( bc1 )
          case( 0 )
            do k = 1,n3
              do j = 1,n2
                fil(1,j,k) = f(1,j,k) !! extrapolation
              enddo
            enddo
          case( 1 )
            do k = 1,n3
              do j = 1,n2
                fil(1,j,k) = half * (f(1,j,k) + f(2,j,k)) !! even extension
              enddo
            enddo
          case( -1 )
            do k = 1,n3
              do j = 1,n2
                fil(1,j,k) = f(1,j,k) !! extrapolation
              enddo
            enddo
          end select

          do k = 1,n3
            do j = 1,n2
              fil(2:this%n-1,j,k) = fourth*f(1:this%n-2,j,k) + half*f(2:this%n-1,j,k) + fourth*f(3:this%n,j,k)
            enddo
          enddo
                  
          select case( bcn )
          case( 0 )
            do k = 1,n3
              do j = 1,n2
                fil(this%n,j,k) = f(this%n,j,k) !! extrapolation
              enddo
            enddo
          case( 1 )
            do k = 1,n3
              do j = 1,n2
                fil(this%n,j,k) = half * (f(this%n-1,j,k) + f(this%n,j,k)) !! even extension
              enddo
            enddo
          case( -1 )
            do k = 1,n3
              do j = 1,n2
                fil(this%n,j,k) = f(this%n,j,k) !! extrapolation
              enddo
            enddo
          end select

        end select
    
    end subroutine
    
    subroutine filter2(this, f, fil, n1, n3, bc1_,  bcn_)
    
        class( box2 ), intent(in) :: this
        integer, intent(in) :: n1, n3
        real(rkind), dimension(n1,this%n,n3), intent(in) :: f
        real(rkind), dimension(n1,this%n,n3), intent(out) :: fil
        integer, optional, intent(in) :: bc1_, bcn_
        integer :: k, bc1, bcn

        select case (this%periodic)
        case (.TRUE.)
          do k=1,n3
            fil(:,1         ,k) = fourth*f(:,  this%n  ,k) + half*f(:,         1,k) + fourth*f(:,       2,k)
            fil(:,2:this%n-1,k) = fourth*f(:,1:this%n-2,k) + half*f(:,2:this%n-1,k) + fourth*f(:,3:this%n,k)
            fil(:,this%n    ,k) = fourth*f(:,  this%n-1,k) + half*f(:,    this%n,k) + fourth*f(:,       1,k)
          end do
        case (.FALSE.)

          !! Ensure BC flags are set correctly
          if (present(bc1_)) then
              bc1 = bc1_
              if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                  call GracefulExit("Incorrect Y-boundary specification for bc1 (should be 0, 1 or -1)", 324)
              end if
          else
              bc1 = 0
          end if

          if (present(bcn_)) then
              bcn = bcn_
              if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                  call GracefulExit("Incorrect Y-boundary specification for bcn (should be 0, 1 or -1)", 324)
              end if
          else
              bcn = 0
          end if

          select case( bc1 )
          case( 0 )
            do k = 1,n3
              fil(:,1,k) = f(:,1,k) !! extrapolation
            enddo
          case( 1 )
            do k = 1,n3
              fil(:,1,k) = half * (f(:,1,k) + f(:,2,k)) !! even extension
            enddo
          case( -1 )
            do k = 1,n3
              fil(:,1,k) = f(:,1,k) !! extrapolation
            enddo
          end select

          do k = 1,n3
            fil(:,2:this%n-1,k) = fourth*f(:,1:this%n-2,k) + half*f(:,2:this%n-1,k) + fourth*f(:,3:this%n,k)
          enddo

          select case( bcn )
          case( 0 )
            do k = 1,n3
              fil(:,this%n,k) = f(:,this%n,k) !! extrapolation
            enddo
          case( 1 )
            do k = 1,n3
              fil(:,this%n,k) = half * (f(:,this%n-1,k) + f(:,this%n,k)) !! even extension
            enddo
          case( -1 )
            do k = 1,n3
              fil(:,this%n,k) = f(:,this%n,k) !! extrapolation
            enddo
          end select

        end select
    
    end subroutine

    subroutine filter3_real(this, f, fil, n1, n2, bc1_, bcn_)
    
        class( box2 ), intent(in) :: this
        integer, intent(in) :: n1, n2
        real(rkind), dimension(n1,n2,this%n), intent(in) :: f
        real(rkind), dimension(n1,n2,this%n), intent(out) :: fil
        integer, optional, intent(in) :: bc1_, bcn_
    
        integer :: bc1, bcn

        select case (this%periodic)
        case (.TRUE.)
            fil(:,:,1         ) = fourth*f(:,:,  this%n  ) + half*f(:,:,         1) + fourth*f(:,:,       2)
            fil(:,:,2:this%n-1) = fourth*f(:,:,1:this%n-2) + half*f(:,:,2:this%n-1) + fourth*f(:,:,3:this%n)
            fil(:,:,this%n    ) = fourth*f(:,:,  this%n-1) + half*f(:,:,    this%n) + fourth*f(:,:,       1)
        case (.FALSE.)

          !! Ensure BC flags are set correctly
          if (present(bc1_)) then
              bc1 = bc1_
              if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                  call GracefulExit("Incorrect Z-boundary specification for bc1 (should be 0, 1 or -1)", 324)
              end if
          else
              bc1 = 0
          end if

          if (present(bcn_)) then
              bcn = bcn_
              if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                  call GracefulExit("Incorrect Z-boundary specification for bcn (should be 0, 1 or -1)", 324)
              end if
          else
              bcn = 0
          end if

          select case( bcn )
          case( 0 )
              fil(:,:,1) = f(:,:,1) !! extrapolation
          case( 1 )
              fil(:,:,1) = half * (f(:,:,1) + f(:,:,2)) !! even extension
          case( -1 )
              fil(:,:,1) = f(:,:,1) !! extrapolation
          end select
                    
          fil(:,:,2:this%n-1) = fourth*f(:,:,1:this%n-2) + half*f(:,:,2:this%n-1) + fourth*f(:,:,3:this%n)

          select case( bcn )
          case( 0 )
              fil(:,:,this%n) = f(:,:,this%n) !! extrapolation
          case( 1 )
              fil(:,:,this%n) = half * (f(:,:,this%n-1) + f(:,:,this%n)) !! even extension
          case( -1 )
              fil(:,:,this%n) = f(:,:,this%n) !! extrapolation
          end select

        end select
    
    end subroutine
    
    subroutine filter3_cmplx(this, f, fil, n1, n2, bc1_, bcn_)
    
        class( box2 ), intent(in) :: this
        integer, intent(in) :: n1, n2
        complex(rkind), dimension(n1,n2,this%n), intent(in) :: f
        complex(rkind), dimension(n1,n2,this%n), intent(out) :: fil
        integer, optional, intent(in) :: bc1_, bcn_
    
        integer :: bc1, bcn

        select case (this%periodic)
        case (.TRUE.)
            fil(:,:,1         ) = fourth*f(:,:,  this%n  ) + half*f(:,:,         1) + fourth*f(:,:,       2)
            fil(:,:,2:this%n-1) = fourth*f(:,:,1:this%n-2) + half*f(:,:,2:this%n-1) + fourth*f(:,:,3:this%n)
            fil(:,:,this%n    ) = fourth*f(:,:,  this%n-1) + half*f(:,:,    this%n) + fourth*f(:,:,       1)
        case (.FALSE.)

          !! Ensure BC flags are set correctly
          if (present(bc1_)) then
              bc1 = bc1_
              if ( (bc1 /= 0) .AND. (bc1 /= 1) .AND. (bc1 /= -1) ) then
                  call GracefulExit("Incorrect Z-boundary specification for bc1 (should be 0, 1 or -1)", 324)
              end if
          else
              bc1 = 0
          end if

          if (present(bcn_)) then
              bcn = bcn_
              if ( (bcn /= 0) .AND. (bcn /= 1) .AND. (bcn /= -1) ) then
                  call GracefulExit("Incorrect Z-boundary specification for bcn (should be 0, 1 or -1)", 324)
              end if
          else
              bcn = 0
          end if

          select case( bcn )
          case( 0 )
              fil(:,:,1) = f(:,:,1) !! extrapolation
          case( 1 )
              fil(:,:,1) = half * (f(:,:,1) + f(:,:,2)) !! even extension
          case( -1 )
              fil(:,:,1) = f(:,:,1) !! extrapolation
          end select
                    
          fil(:,:,2:this%n-1) = fourth*f(:,:,1:this%n-2) + half*f(:,:,2:this%n-1) + fourth*f(:,:,3:this%n)

          select case( bcn )
          case( 0 )
              fil(:,:,this%n) = f(:,:,this%n) !! extrapolation
          case( 1 )
              fil(:,:,this%n) = half * (f(:,:,this%n-1) + f(:,:,this%n)) !! even extension
          case( -1 )
              fil(:,:,this%n) = f(:,:,this%n) !! extrapolation
          end select

        end select
    
    end subroutine

end module
