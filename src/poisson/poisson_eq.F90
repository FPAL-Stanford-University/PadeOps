module poisson_eq

    use kind_parameters, only : rkind
    use fft_2d_stuff, only: fft_2d
    use decomp_2d
    use constants, only: r_eps

    implicit none 
    private
    public :: poisson

    type poisson
        type (fft_2d) :: fft_xy
        
        character(len=4) :: scheme_xy = "four"
        character(len=4) :: scheme_z  = "cd08"
      
        integer :: z_bc = 0 !BC in z direction; 0-none, 1-dirichlet, 2-neumann


        real(rkind) :: dx, dy, dz       ! Grid resolution
        real(rkind) :: nx_g, ny_g, nz_g ! Global dimensions for problem 
        
        real(rkind), dimension(:,:), allocatable :: k1, k2, kabs_sq ! Exact or modified wavenumbers

        real(rkind), dimension(:,:,:,:), allocatable :: Penta_z, Tri_z

        complex(rkind), dimension(:,:,:), allocatable :: rhs_xy_hat

        logical :: initialized = .false. 

        logical :: have_zero_wavenumber 
        contains 
            procedure :: init
            procedure :: destroy
            procedure :: Solve

            procedure, private :: generatePenta
            procedure, private :: generateTri 
    end type



contains

    function init(this,nx_, ny_, nz_,dx_, dy_, dz_,scheme_xy_, scheme_z_, z_bc_) result(ierr)
        class(poisson), intent(inout) :: this
        real(rkind), intent(in) :: dx_, dy_, dz_
        integer, intent(in) :: nx_, ny_, nz_
        character(len=4), intent(in) :: scheme_xy_
        character(len=4), intent(in) :: scheme_z_
        integer, intent(in) :: z_bc_

        integer :: ierr

        if (z_bc_ .ne. 0) then
            call decomp_2d_abort(400,"This BC is currently not supported / incomplete")
        end if

        this%dx = dx_
        this%dy = dy_
        this%dz = dz_
        this%nx_g = nx_ 
        this%ny_g = ny_ 
        this%nz_g = nz_ 
        
        ierr = this%fft_xy%init(nx_,ny_,nz_,"y",dx_, dy_, .true.)
        
        if (ierr .ne. 0) then
            return
        end if 
       
        call this%fft_xy%alloc_output(this%rhs_xy_hat) 
        allocate(this%k1(this%fft_xy%spectral%ysz(1),this%fft_xy%spectral%ysz(2)))            
        allocate(this%k2(this%fft_xy%spectral%ysz(1),this%fft_xy%spectral%ysz(2)))            
        allocate(this%kabs_sq(this%fft_xy%spectral%ysz(1),this%fft_xy%spectral%ysz(2)))            

        select case (scheme_xy_)
        case("four")
            this%k1 = this%fft_xy%k1
            this%k2 = this%fft_xy%k2
            this%kabs_sq = this%fft_xy%kabs_sq
            if (this%kabs_sq(1,1) == 0._rkind) then
                this%have_zero_wavenumber = .true.
            else
               this%have_zero_wavenumber = .false.  
            end if 
        case("cd10")
            ! Get the modified wavenumbers
            ! <incomplete>
            ierr = 2
            return
        case("cd06")
            ! Get the modified wavenumbers
            ! <incomplete>
            ierr = 2
            return 
        case default
            call decomp_2d_abort(401,"Verify that the selected scheme in xy plane is either: four, cd10 or cd06")
        end select
       
        this%scheme_xy = scheme_xy_

        select case (scheme_z_) 
        case ("cd08")
            allocate (this%Penta_z(size(this%k1,1),size(this%k1,2),nz_,11))
            call this%generatePenta(this%kabs_sq, this%Penta_z)
        case ("cd06")   
            allocate (this%Tri_z(size(this%k1,1),size(this%k1,2),nz_,3))
            call this%generateTri(this%kabs_sq, this%Tri_z)
        case default
            call decomp_2d_abort(402,"Verify that the selected scheme in z direction is either: cd08 or cd06")
        end select

       

        this%initialized = .true. 
        ierr = 0

    end function 


    subroutine destroy(this)
        class(poisson), intent(inout) :: this 
        
        call this%fft_xy%destroy
        if (allocated(this%Penta_z)) deallocate(this%Penta_z)
        if (allocated(this%Tri_z)) deallocate(this%Tri_z)

        if (allocated(this%k1)) deallocate(this%k1)
        if (allocated(this%k2)) deallocate(this%k2)
        if (allocated(this%kabs_sq)) deallocate(this%kabs_Sq)

        if (allocated(this%rhs_xy_hat)) deallocate(this%rhs_xy_hat)

        this%initialized = .false. 
    end subroutine

    subroutine generatePenta(this,kabs,Penta)
        class(poisson), intent(inout) :: this
        real(rkind), intent(in), dimension(:,:) :: kabs
        real(rkind), intent(out), dimension(:,:,:,:) :: Penta
    
        Penta = 0._rkind

    end subroutine

    subroutine generateTri(this,kabs,Tri)
        class(poisson), intent(inout) :: this
        real(rkind), intent(in), dimension(:,:) :: kabs
        real(rkind), intent(out), dimension(:,:,:,:) :: Tri

        Tri = 0._rkind
    
    end subroutine

    subroutine solve(this,rhs,fsol) 

        class(poisson), intent(inout) :: this
        real(rkind), dimension(this%fft_xy%physical%ysz(1),this%fft_xy%physical%ysz(2),this%fft_xy%physical%ysz(3)), intent(in) :: rhs
        real(rkind), dimension(this%fft_xy%physical%ysz(1),this%fft_xy%physical%ysz(2),this%fft_xy%physical%ysz(3)), intent(out) :: fsol
        integer :: k

        if (.not.this%initialized) then
            fsol = 0._rkind
        end if

        ! First get the xy planar FFT 
        call this%fft_xy%fft2(rhs,this%rhs_xy_hat)
       
        ! Trial - just 2d poisson equation
        do k = 1,this%fft_xy%physical%ysz(3)
            this%rhs_xy_hat(:,:,k) = -this%rhs_xy_hat(:,:,k)/(this%kabs_sq + r_eps)
        end do

        ! Set the mean value 
        if (this%have_zero_wavenumber) then
            this%rhs_xy_hat(1,1,:) = 0._rkind
        end if 

        ! Now invert the xy FFT transform  
        call this%fft_xy%ifft2(this%rhs_xy_hat,fsol)

20 format(1x,16D10.3)
    end subroutine
end module 
