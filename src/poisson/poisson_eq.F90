module poisson_eq

    use kind_parameters, only : rkind
    use fft_2d_stuff, only: fft_2d
    use decomp_2d
    use constants, only: r_eps, one

    implicit none 
    private
    public :: poisson

    type poisson
        private 
        type (fft_2d) :: fft_xy
        
        character(len=4) :: scheme_xy = "four"
        character(len=4) :: scheme_z  = "cd08"
      
        integer :: z_bc = 0 !BC in z direction; 0-none, 1-dirichlet, 2-neumann


        real(rkind) :: dx, dy, dz       ! Grid resolution
        integer     :: nx_g, ny_g, nz_g ! Global dimensions for problem 
        
        real(rkind), dimension(:,:), allocatable :: k1, k2, kabs_sq ! Exact or modified wavenumbers

        real(rkind), dimension(:,:,:,:), allocatable :: Penta_z, Tri_z

        complex(rkind), dimension(:,:,:), allocatable :: rhs_xy_hat_in_y
        complex(rkind), dimension(:,:,:), allocatable :: rhs_xy_hat_in_z
        complex(rkind), dimension(:,:,:), allocatable :: rhs_4_solve

        logical :: initialized = .false. 

        logical :: have_zero_wavenumber 
        contains 
            procedure :: init
            procedure :: destroy
            procedure :: Solve

            procedure, private :: solve_cd08
            procedure, private :: get_rhs_cd08
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
        real(rkind), dimension(:,:,:,:), allocatable ::  tmp_in_y, tmp_in_z
        integer :: ierr, k

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
       
        call this%fft_xy%alloc_output(this%rhs_xy_hat_in_y) 
        if (allocated(this%rhs_xy_hat_in_z)) deallocate(this%rhs_xy_hat_in_z)
        if (allocated(this%rhs_4_solve    )) deallocate(this%rhs_4_solve)
        allocate(this%rhs_xy_hat_in_z(this%fft_xy%spectral%zsz(1),this%fft_xy%spectral%zsz(2),this%fft_xy%spectral%zsz(3)))
        allocate(this%rhs_4_solve    (this%fft_xy%spectral%zsz(1),this%fft_xy%spectral%zsz(2),this%fft_xy%spectral%zsz(3)))

        allocate(tmp_in_y(this%fft_xy%spectral%ysz(1), this%fft_xy%spectral%ysz(2), this%fft_xy%spectral%ysz(3),2))
        allocate(tmp_in_z(this%fft_xy%spectral%zsz(1), this%fft_xy%spectral%zsz(2), this%fft_xy%spectral%zsz(3),2))
         
        allocate(this%k1(this%fft_xy%spectral%zsz(1),this%fft_xy%spectral%zsz(2)))            
        allocate(this%k2(this%fft_xy%spectral%zsz(1),this%fft_xy%spectral%zsz(2)))            
        allocate(this%kabs_sq(this%fft_xy%spectral%zsz(1),this%fft_xy%spectral%zsz(2)))            

        select case (scheme_xy_)
        case("four")
            ! Generate dummy 
            do k = 1,this%fft_xy%spectral%ysz(3)
                tmp_in_y(:,:,k,1) = this%fft_xy%k1
                tmp_in_y(:,:,k,2) = this%fft_xy%k2
            end do 
            
            ! Transpose dummy from y -> z
            call transpose_y_to_z(tmp_in_y(:,:,:,1),tmp_in_z(:,:,:,1),this%fft_xy%spectral)
            call transpose_y_to_z(tmp_in_y(:,:,:,2),tmp_in_z(:,:,:,2),this%fft_xy%spectral)

            this%k1 = tmp_in_z(:,:,1,1)
            this%k2 = tmp_in_z(:,:,1,2)
            this%kabs_sq = this%k1*this%k1 + this%k2*this%k2
            
            if (abs(this%kabs_sq(1,1)) .LT. real(0.0D-20,rkind)) then
                this%have_zero_wavenumber = .true.
            else
                this%have_zero_wavenumber = .false.  
            end if 
        case("cd10")
            ! Get the modified wavenumbers
            ! <incomplete>
            ierr = 2
            return
        case("cd04")
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
            allocate (this%Penta_z(size(this%k1,1),size(this%k1,2),nz_,10))
            call this%generatePenta
        case ("cd04")   
            allocate (this%Tri_z(size(this%k1,1),size(this%k1,2),nz_,3))
            call this%generateTri(this%kabs_sq, this%Tri_z)
        case default
            call decomp_2d_abort(402,"Verify that the selected scheme in z direction is either: cd08 or cd04")
        end select

        
        deallocate(tmp_in_y,tmp_in_z)

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
        if (allocated(this%kabs_sq)) deallocate(this%kabs_sq)

        if (allocated(this%rhs_xy_hat_in_y)) deallocate(this%rhs_xy_hat_in_y)
        if (allocated(this%rhs_xy_hat_in_z)) deallocate(this%rhs_xy_hat_in_z)
        if (allocated(this%rhs_4_solve)) deallocate(this%rhs_4_solve)

        this%initialized = .false. 
    end subroutine

    subroutine generatePenta(this)
   
        class(poisson), intent(inout) :: this
        real(rkind), parameter :: alpha  = 344._rkind/1179._rkind 
        real(rkind), parameter :: beta   =  23._rkind/2358._rkind 
        real(rkind), parameter :: alpha2 =   1._rkind/  10._rkind  
        real(rkind), parameter :: aa     = 320._rkind/ 393._rkind
        real(rkind), parameter :: bb     = 310._rkind/ 393._rkind
        real(rkind), parameter :: aa2    =   6._rkind/   5._rkind
        
        integer                :: k, i, j 

        associate (bt   => this%penta_z(:,:,:,1), b   => this%penta_z(:,:,:,2), d => this%penta_z(:,:,:,3),  &
                   a    => this%penta_z(:,:,:,4), at  => this%penta_z(:,:,:,5),                              &
                   e    => this%penta_z(:,:,:,6), obc => this%penta_z(:,:,:,7),                              &
                   f    => this%penta_z(:,:,:,8), g   => this%penta_z(:,:,:,9),                              &
                   eobc => this%penta_z(:,:,:,10)                                                            )
   
            do k = 1,this%fft_xy%spectral%zsz(3)
               at(:,:,k) = beta *this%kabs_sq - bb/(4._rkind*this%dz*this%dz)
               bt(:,:,k) = beta *this%kabs_sq - bb/(4._rkind*this%dz*this%dz)
               a (:,:,k) = alpha*this%kabs_sq - aa/(         this%dz*this%dz)
               b (:,:,k) = alpha*this%kabs_sq - aa/(         this%dz*this%dz)
               d (:,:,k) = this%kabs_sq + bb/(2._rkind*this%dz*this%dz) + 2._rkind*aa/(this%dz*this%dz)
            end do 

            bt(:,:,          1) =  0._rkind 
            b (:,:,          1) =  0._rkind
            d (:,:,          1) =  3._rkind/(2._rkind*this%dz)
            a (:,:,          1) = -2._rkind/this%dz
            at(:,:,          1) =  1._rkind/(2._rkind*this%dz)

            bt(:,:,          2) = 0._rkind 
            b (:,:,          2) = alpha2*this%kabs_sq - aa2/(         this%dz*this%dz)
            d (:,:,          2) = this%kabs_sq + 2._rkind*aa2/(this%dz*this%dz)
            a (:,:,          2) = alpha2*this%kabs_sq - aa2/(         this%dz*this%dz)
            at(:,:,          2) = 0._rkind 

            bt(:,:,this%nz_g-1) = 0._rkind 
            b (:,:,this%nz_g-1) = alpha2*this%kabs_sq - aa2/(         this%dz*this%dz)
            d (:,:,this%nz_g-1) = this%kabs_sq + 2._rkind*aa2/(this%dz*this%dz)
            a (:,:,this%nz_g-1) = alpha2*this%kabs_sq - aa2/(         this%dz*this%dz)
            at(:,:,this%nz_g-1) = 0._rkind 

            bt(:,:,  this%nz_g) = -1._rkind/(2._rkind*this%dz) 
            b (:,:,  this%nz_g) =  2._rkind/this%dz 
            d (:,:,  this%nz_g) = -3._rkind/(2._rkind*this%dz)
            a (:,:,  this%nz_g) =  0._rkind 
            at(:,:,  this%nz_g) =  0._rkind 
            
            if (this%have_zero_wavenumber) then
                do j = 1,this%fft_xy%spectral%zsz(2)
                    do i = 1,this%fft_xy%spectral%zsz(1)
                        if (abs(this%kabs_sq(i,j)) .LE. real(0.D-20,rkind)) then
                            a (i,j,1) = 0._rkind
                            b (i,j,1) = 0._rkind
                            at(i,j,1) = 0._rkind
                            bt(i,j,1) = 0._rkind
                            d (i,j,1) = 1._rkind
                        end if 
                    end do 
                end do 
            end if 


            ! Step 1
            obc(:,:,1) = one/d(:,:,1)

            ! Step 2
            obc(:,:,2) = one/(d(:,:,2) - b(:,:,2)*a(:,:,1)*obc(:,:,1))

            ! Step 3
            e(:,:,1) = a(:,:,1)
            f(:,:,2) = b(:,:,2)*obc(:,:,1)
            
            do i = 3,this%nz_g
                g(:,:,i) = bt(:,:,i)*obc(:,:,i-2)
                e(:,:,i-1) = a(:,:,i-1) - f(:,:,i-1)*at(:,:,i-2)
                f(:,:,i) = (b(:,:,i) - g(:,:,i)*e(:,:,i-2))*obc(:,:,i-1)
                obc(:,:,i) = one/(d(:,:,i) - f(:,:,i)*e(:,:,i-1) - g(:,:,i)*at(:,:,i-2))
            end do 

            eobc = e*obc
       
        end associate  
    end subroutine

    subroutine generateTri(this,kabs,Tri)
        class(poisson), intent(inout) :: this
        real(rkind), intent(in), dimension(:,:) :: kabs
        real(rkind), intent(out), dimension(:,:,:,:) :: Tri

        Tri = 0._rkind
    
    end subroutine

    pure subroutine get_rhs_cd08(this,in_rhs,out_rhs)
        class(poisson), intent(in) :: this
        complex(rkind), dimension(this%fft_xy%spectral%zsz(1),this%fft_xy%spectral%zsz(2),this%fft_xy%spectral%zsz(3)), intent(in) :: in_rhs
        complex(rkind), dimension(this%fft_xy%spectral%zsz(1),this%fft_xy%spectral%zsz(2),this%fft_xy%spectral%zsz(3)), intent(out) :: out_rhs
        real(rkind), parameter :: alpha  = 344._rkind/1179._rkind 
        real(rkind), parameter :: beta   =  23._rkind/2358._rkind 
        real(rkind), parameter :: alpha2 =   1._rkind/  10._rkind  
        integer :: k
       
        out_rhs(:,:,1) = 0._rkind
        out_rhs(:,:,2) = -alpha2*in_rhs(:,:,1) - in_rhs(:,:,2) - alpha2*in_rhs(:,:,3) 
        do k = 3,this%nz_g-2
            out_rhs(:,:,k) = -beta*in_rhs(:,:,k-2) - alpha*in_rhs(:,:,k-1) - in_rhs(:,:,k) - alpha*in_rhs(:,:,k+1) - beta*in_rhs(:,:,k+2)
        end do 
        out_rhs(:,:,this%nz_g-1) = -alpha2*in_rhs(:,:,this%nz_g-2) - in_rhs(:,:,this%nz_g-1) - alpha2*in_rhs(:,:,this%nz_g) 
        out_rhs(:,:,this%nz_g  ) = 0._rkind
         
    end subroutine

    subroutine solve_cd08(this,y)
        class(poisson), intent(in) :: this
        complex(rkind), dimension(this%fft_xy%spectral%zsz(1),this%fft_xy%spectral%zsz(2),this%fft_xy%spectral%zsz(3)), intent(inout) :: y
        integer :: k 
        
        ! Step 1
        y(:,:,2) = y(:,:,2) - this%Penta_z(:,:,2,8)*y(:,:,1)
        do k = 3,this%nz_g
            y(:,:,k) = y(:,:,k) - this%Penta_z(:,:,k,9)*y(:,:,k-2) - this%Penta_z(:,:,k,8)*y(:,:,k-1)
        end do 

        ! Step 2
        y(:,:,this%nz_g) = y(:,:,this%nz_g)*this%Penta_z(:,:,this%nz_g,7)
        
        y(:,:,this%nz_g-1) = y(:,:,this%nz_g-1)*this%Penta_z(:,:,this%nz_g-1,7) - this%Penta_z(:,:,this%nz_g-1,10)*y(:,:,this%nz_g)
        do k = this%nz_g-2,1,-1
            y(:,:,k) = y(:,:,k)*this%Penta_z(:,:,k,7) - y(:,:,k+2)*this%Penta_z(:,:,k,5)*this%Penta_z(:,:,k,7) - y(:,:,k+1)*this%Penta_z(:,:,k,10)
        end do 

    end subroutine 

    subroutine solve(this,rhs,fsol) 
        class(poisson), intent(inout) :: this
        real(rkind), dimension(this%fft_xy%physical%ysz(1),this%fft_xy%physical%ysz(2),this%fft_xy%physical%ysz(3)), intent(in) :: rhs
        real(rkind), dimension(this%fft_xy%physical%ysz(1),this%fft_xy%physical%ysz(2),this%fft_xy%physical%ysz(3)), intent(out) :: fsol
        
        if (.not.this%initialized) then
            fsol = 0._rkind
        end if

        ! First get the xy planar FFT 
        call this%fft_xy%fft2(rhs,this%rhs_xy_hat_in_y)
      
        ! Now transform to z decomp
        call transpose_y_to_z(this%rhs_xy_hat_in_y,this%rhs_xy_hat_in_z,this%fft_xy%spectral)

        ! Compute RHS for penta solve (in Z decomp)
        call this%get_rhs_cd08(this%rhs_xy_hat_in_z,this%rhs_4_solve)

        ! Solve the Pentadiagonal system (in Z decomp)
        call this%solve_cd08(this%rhs_4_solve)

        ! Now transform back from z to y
        call transpose_z_to_y(this%rhs_4_solve,this%rhs_xy_hat_in_y,this%fft_xy%spectral)
        
        ! Now invert the xy FFT transform  
        call this%fft_xy%ifft2(this%rhs_xy_hat_in_y,fsol)

        ! Done 
    end subroutine


end module 
