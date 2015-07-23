module fftStuff

    use, intrinsic :: ISO_C_BINDING
    use kind_parameters, only: rkind
    use constants, only: zero,one,half,two,pi,imi
    implicit none

    private
    public :: ffts

    ! include "fftw3.F90"
    include "fftw3.f03"

    type ffts
        
        private
        
        character(len=4)                         :: fft_lib="FFTW"    ! FFT library to use. "FFTW", "FFTP"
        
        integer(kind=8) :: plan_fwd
        integer(kind=8) :: plan_bwd
        logical         :: initialized=.FALSE.
        logical         :: isEven=.TRUE.
        
        integer :: nd=1
        integer, dimension(:), allocatable :: npts
        
        real(rkind), dimension(:), allocatable :: dspace

        real(rkind), allocatable, dimension(:)       :: k1d, k1d_sq        ! Wavenumbers for 1D transform
        real(rkind), allocatable, dimension(:,:,:)   :: k2d                ! Wavenumbers for 2D transform
        real(rkind), allocatable, dimension(:,:,:,:) :: k3d                ! Wavenumbers for 3D transform



        integer, allocatable, dimension(:) :: split

        contains

        procedure, private :: init1d
        procedure, private :: init2d
        procedure, private :: init3d
        generic            :: init => init1d,init2d,init3d ! Initialize type
        procedure          :: destroy 

        procedure, private :: fft1d_half
        procedure, private :: fft1d_full
        generic :: fft => fft1d_half, fft1d_full!,fft2d,fft3d         ! Compute forward transform
        procedure, private :: ifft1d
        generic :: ifft => ifft1d!,ifft2d,ifft3d     ! Compute inverse transform
        
        procedure, private :: getk1d
        generic            :: k => getk1d!, getk2d, getk3d    ! Get wave number given indices
        !
    
        procedure, private :: nullOddBall
        procedure :: fourder1
        procedure :: fourder2
    end type


contains

    function nullOddBall(this,f) result(fout) 
        class( ffts ), intent(in) :: this
        complex(rkind), dimension(this%npts(1)), intent(in) :: f
        complex(rkind), dimension(this%npts(1)) :: fout

        fout = f
        if (this%isEven) fout(this%split(1)) = zero

    end function 

    function fourder1(this, f) result(df)
        class( ffts ), intent(in) :: this
        real(rkind), dimension(this%npts(1)), intent(in) :: f
        real(rkind), dimension(this%npts(1)) :: df 
        
        df = this%ifft( imi * this%k1d * this%nullOddBall( this%fft(f, .TRUE.) ) )
        !df = this%ifft( imi * this%k1d * this%fft(f) )
    end function 
    
    function fourder2(this, f) result(df)
        class( ffts ), intent(in) :: this
        real(rkind), dimension(this%npts(1)), intent(in) :: f
        real(rkind), dimension(this%npts(1)) :: df 

        df = this%ifft( -this%k1d_sq * this%nullOddBall( this%fft(f) ) )
    end function 

    function init1d(this, nx, dx) result(ierr)
        
        class( ffts ), intent(inout) :: this
        integer, intent(in) :: nx
        real(rkind), intent(in) :: dx
        integer :: ierr

        complex(rkind), dimension(:), allocatable :: arr_in, arr_out
        real(rkind), dimension(:), allocatable :: arr_in_real, arr_out_real

        if ( this%initialized ) then
            ierr = 100
            return
        end if
        
        this%nd = 1

        if ( .not. allocated(this%npts) ) allocate( this%npts(this%nd) )
        this%npts = [nx]

        if ( .not. allocated(this%split) ) allocate( this%split(this%nd) )
        this%split = this%npts/2+1

        if ( .not. allocated(this%dspace) ) allocate( this%dspace(this%nd) )
        this%dspace = [dx]

        if ( .not. allocated(this%k1d) ) allocate( this%k1d( this%npts(1) ) )
        this%k1d = GetWaveNums(nx,dx)
        if ( .not. allocated(this%k1d_sq) ) allocate( this%k1d_sq( this%npts(1) ) )
        this%k1d_sq = this%k1d*this%k1d 

        allocate( arr_in_real(nx) )
        allocate( arr_in(nx) )
        allocate( arr_out_real(nx) )
        allocate( arr_out(nx) )
        call dfftw_plan_dft_r2c_1d( this%plan_fwd, nx, arr_in_real, arr_out, FFTW_FORWARD,  FFTW_EXHAUSTIVE )
        call dfftw_plan_dft_c2r_1d( this%plan_bwd, nx, arr_in, arr_out_real, FFTW_BACKWARD, FFTW_EXHAUSTIVE )
        deallocate( arr_in_real )
        deallocate( arr_in  )
        deallocate( arr_out_real )
        deallocate( arr_out )

        select case (mod(nx,2))
        case (1)
            this%isEven = .false.
        case (0)
            this%isEven = .true.
        end select 

        this%initialized = .TRUE.
        ierr = 0
    end function

    function fft1d_half(this, f, half) result (fhat)
        class( ffts ), intent(in) :: this
        real(rkind), dimension(this%npts(1)), intent(in) :: f
        logical, intent(in) :: half
        complex(rkind), dimension(this%split(1)) :: fhat

        call dfftw_execute_dft(this%plan_fwd, f, fhat)

    end function

    function fft1d_full(this, f) result (fhat)
        class( ffts ), intent(in) :: this
        real(rkind), dimension(this%npts(1)), intent(in) :: f
        complex(rkind), dimension(this%npts(1)) :: fhat

        integer :: i

        call dfftw_execute_dft(this%plan_fwd, f, fhat(1:this%split(1)))

        do i=0,this%split(1)-2
            fhat(this%npts(1)-i) = conjg( fhat(2+i) )
        end do

    end function

    function ifft1d(this, fhat) result (f)
        class( ffts ), intent(in) :: this
        complex(rkind), dimension(this%split(1)), intent(in) :: fhat
        real(rkind), dimension(this%npts(1)) :: f

        call dfftw_execute_dft(this%plan_bwd, fhat, f)

        f = f / real(this%npts(1),rkind)

    end function

    function getk1d(this) result(karray)
        class( ffts ), intent(in) :: this
        real(rkind), dimension(this%npts(1)) :: karray

        karray = this%k1d
    end function
    
    function init2d(this, nx, dx, ny, dy) result(ierr)
        
        class( ffts ), intent(inout) :: this
        integer, intent(in) :: nx, ny
        real(rkind), intent(in) :: dx, dy
        integer :: ierr

        real(rkind), dimension(nx) :: kx
        real(rkind), dimension(ny) :: ky
        complex(rkind), dimension(:,:), allocatable :: arr_in, arr_out
        real(rkind), dimension(:,:), allocatable :: arr_in_real

        integer :: i,j

        if ( this%initialized ) then
            ierr = 100
            return
        end if
        
        this%nd = 2
        
        if ( .not. allocated(this%npts) ) allocate( this%npts(this%nd) )
        this%npts = [nx, ny]

        if ( .not. allocated(this%dspace) ) allocate( this%dspace(this%nd) )
        this%dspace = [dx, dy]

        if ( .not. allocated(this%k2d) ) allocate( this%k2d( this%npts(1),this%npts(2),this%nd ) )
        kx = GetWaveNums(nx,dx)
        ky = GetWaveNums(ny,dy)

        do j=1,ny
            do i=1,nx
                this%k2d(i,j,1) = kx(i)
                this%k2d(i,j,2) = ky(j)
            end do
        end do

        allocate( arr_in_real(nx,ny) )
        allocate( arr_in (nx,ny) )
        allocate( arr_out(nx,ny) )
        call dfftw_plan_dft_2d( this%plan_fwd, nx, ny, arr_in_real, arr_out, FFTW_FORWARD,  FFTW_EXHAUSTIVE )
        call dfftw_plan_dft_2d( this%plan_bwd, nx, ny, arr_in, arr_out, FFTW_BACKWARD, FFTW_EXHAUSTIVE )
        deallocate( arr_in  )
        deallocate( arr_out )


        this%initialized = .TRUE.
        ierr = 0
    end function

    function init3d(this, nx, dx, ny, dy, nz, dz) result(ierr)
        
        class( ffts ), intent(inout) :: this
        integer, intent(in) :: nx, ny, nz
        real(rkind), intent(in) :: dx, dy, dz
        integer :: ierr

        real(rkind), dimension(nx) :: kx
        real(rkind), dimension(ny) :: ky
        real(rkind), dimension(nz) :: kz
        complex(rkind), dimension(:,:,:), allocatable :: arr_in, arr_out
        real(rkind), dimension(:,:,:), allocatable :: arr_in_real

        integer :: i,j,k

        if ( this%initialized ) then
            ierr = 100
            return
        end if
        
        this%nd = 3
        
        if ( .not. allocated(this%npts) ) allocate( this%npts(this%nd) )
        this%npts = [nx, ny, nz]

        if ( .not. allocated(this%dspace) ) allocate( this%dspace(this%nd) )
        this%dspace = [dx, dy, dz]

        if ( .not. allocated(this%k3d) ) allocate( this%k3d( this%npts(1),this%npts(2),this%npts(3),this%nd ) )
        kx = GetWaveNums(nx,dx)
        ky = GetWaveNums(ny,dy)
        kz = GetWaveNums(nz,dz)

        do k=1,nz
            do j=1,ny
                do i=1,nx
                    this%k3d(i,j,k,1) = kx(i)
                    this%k3d(i,j,k,2) = ky(j)
                    this%k3d(i,j,k,3) = kz(k)
                end do
            end do
        end do

        allocate( arr_in_real(nx,ny,nz) )
        allocate( arr_in (nx,ny,nz) )
        allocate( arr_out(nx,ny,nz) )
        call dfftw_plan_dft_3d( this%plan_fwd, nx, ny, nz, arr_in_real, arr_out, FFTW_FORWARD,  FFTW_EXHAUSTIVE )
        call dfftw_plan_dft_3d( this%plan_bwd, nx, ny, nz, arr_in, arr_out, FFTW_BACKWARD, FFTW_EXHAUSTIVE )
        deallocate( arr_in  )
        deallocate( arr_out )


        this%initialized = .TRUE.
        ierr = 0
    end function

    pure function GetWaveNums(nx,dx) result(k)

        integer, intent(in) :: nx
        real(rkind), intent(in) :: dx
        real(rkind), dimension(nx) :: k

        integer :: i,dummy

        dummy = nx - MOD(nx,2)

        do i = 1,nx
            k(i) = ( -pi + (i-1)*two*pi/real(dummy,rkind) ) / dx
        end do

        k = ifftshift(k)

    end function

    pure function ifftshift(k) result(kshift)

        real(rkind), dimension(:), intent(in) :: k
        real(rkind), dimension(SIZE(k)) :: kshift
        integer :: n

        n = SIZE(k)

        select case ( MOD(n,2) )
        case (0)
            kshift(1:n/2) = k(n/2+1:n)
            kshift(n/2+1:n) = k(1:n/2)
        case (1)
            kshift(1:(n+1)/2) = k((n+1)/2:n)
            kshift((n+1)/2+1:n) = k(1:(n-1)/2)
        end select

    end function


    subroutine destroy(this)
        
        class(ffts), intent(inout) :: this

        if ( allocated(this%k1d) ) deallocate( this%k1d )
        if ( allocated(this%k2d) ) deallocate( this%k2d )
        if ( allocated(this%k3d) ) deallocate( this%k3d )
        if ( allocated(this%npts) ) deallocate( this%npts )
        if ( allocated(this%split) ) deallocate( this%split )
        if ( allocated(this%dspace) ) deallocate( this%dspace )

        call dfftw_destroy_plan ( this%plan_fwd )
        call dfftw_destroy_plan ( this%plan_bwd )

        this%initialized = .FALSE.

    end subroutine

end module
