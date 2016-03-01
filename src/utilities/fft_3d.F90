module fft_3d_stuff

    use kind_parameters, only : rkind
    use decomp_2d
    use constants, only: zero  
    implicit none 
    private
    public :: fft_3d

    include "fftw3.f"
    
    type fft_3d
        private
        type(decomp_info),public    :: physical
        type(decomp_info),public    :: spectral
        
        character(len=1)  :: base_pencil
        real (rkind) :: normfactor, normfactor2d

        complex(rkind), dimension(:,:,:), allocatable :: f_yhat_in_yD
        complex(rkind), dimension(:,:,:), allocatable :: f_xyhat_in_xD
        complex(rkind), dimension(:,:,:), allocatable :: f_xyzhat_in_zD

        complex(rkind), dimension(:,:,:), allocatable :: f_xhat_in_xD
        complex(rkind), dimension(:,:,:), allocatable :: f_xyhat_in_yD
        complex(rkind), dimension(:,:,:), allocatable :: f_zyhat_in_yD

        complex(rkind), dimension(:,:,:), allocatable :: f_zhat_in_zD
        complex(rkind), dimension(:,:,:), allocatable :: f_xyzhat_in_xD

        real(rkind), dimension(:,:,:),allocatable, public :: k1
        real(rkind), dimension(:,:,:),allocatable, public :: k2
        real(rkind), dimension(:,:,:),allocatable, public :: k3
        real(rkind), dimension(:,:,:),allocatable, public :: kabs_sq

        logical :: allocK = .true. 
        integer(kind=8) :: plan_c2c_fwd_x
        integer(kind=8) :: plan_c2c_bwd_x
        integer(kind=8) :: plan_c2c_fwd_z
        integer(kind=8) :: plan_c2c_bwd_z
        integer(kind=8) :: plan_c2c_fwd_y
        integer(kind=8) :: plan_c2c_bwd_y
        integer(kind=8) :: plan_r2c_x
        integer(kind=8) :: plan_c2r_x
        integer(kind=8) :: plan_c2r_y
        integer(kind=8) :: plan_r2c_y
        integer(kind=8) :: plan_c2r_z
        integer(kind=8) :: plan_r2c_z
        integer(kind=8) :: plan_c2c_bwd_y_oop

        logical :: fixOddball = .false. 
        logical :: initialized = .false.

        integer :: fft_plan = FFTW_MEASURE
        contains
            procedure :: init
            procedure :: fft3_y2y
            procedure :: ifft3_y2y
            procedure :: fft3_x2z
            procedure :: ifft3_z2x
            procedure :: fft3_z2x
            procedure :: ifft3_x2z
            procedure :: ifft2_y2x
            procedure :: fft2_x2y
            procedure :: destroy
            procedure :: alloc_input
            procedure, private :: alloc_3d_real
            procedure, private :: alloc_3d_cmplx
            procedure, private :: alloc_4d 
            generic   :: alloc_output => alloc_3d_real, alloc_3d_cmplx, alloc_4d
    end type 


contains

    function init(this,nx_global,ny_global,nz_global,base_pencil_, dx, dy,dz, exhaustive, fixOddball_, allocK) result(ierr)
        class(fft_3d), intent(inout) :: this
        integer, intent(in) :: nx_global, ny_global, nz_global
        character(len=1), intent(in) :: base_pencil_
        real(rkind), intent(in) :: dx, dy,dz
        logical, optional :: exhaustive
        logical, optional :: fixOddball_
        logical, optional :: allocK
        integer :: ierr, i, j, k
        integer, dimension(2) :: dims, dummy_periods, dummy_coords
        real(rkind), dimension(:,:), allocatable :: real_arr_2d
        complex(rkind), dimension(:,:), allocatable :: cmplx_arr_2d, cmplx_arr_2d_oop
        integer :: n_sizeact, n_sizeinput, n_sizeoutput, n_chunk, n_howmany, n_jump
         
        real(rkind), dimension(:,:,:), allocatable :: k1_in_x
        real(rkind), dimension(:,:,:), allocatable :: k3_in_z
        real(rkind), dimension(:,:,:), allocatable :: k2_in_y
        real(rkind), dimension(:,:,:), allocatable :: temp
        complex(rkind), dimension(:,:,:), allocatable :: dummy_in_z
        complex(rkind), dimension(:,:,:), allocatable :: dummy_in_x

        ! Get the decompotion performed before calling this function
        call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, & 
                dims, dummy_periods, dummy_coords, ierr) 


        if ((base_pencil_ == "y") .or. (base_pencil_ == "z")) then
            if ( (mod(ny_global,2) .ne. 0) .or. (mod(nx_global,2) .ne. 0) .or. (mod(nz_global,2) .ne. 0) ) then
                call decomp_2d_abort(301,"Only even number of data are supported &
              & in X, Y and Z direction. NY_GLOBAL, NX_GLOBAL and NZ_GLOBAL must be an even number.") 
            end if 
        else
            if ( (mod(ny_global,2) .ne. 0) .or. (mod(nx_global,2) .ne. 0) ) then
                call decomp_2d_abort(301,"Only even number of data are supported &
              & in X, Y and Z direction. NY_GLOBAL, NX_GLOBAL and NZ_GLOBAL must be an even number.") 
            end if 
        end if 
        
        if (present(exhaustive)) then
            if (exhaustive) then 
                this%fft_plan = FFTW_EXHAUSTIVE
                if (nrank == 0) then
                    print*, "WARNING: FFT plans are being created using FFTW_EXHAUSTIVE mode. Intialization could take a while!"
                end if
            end if  
        end if 
       
        if (present(fixOddball_)) then
            this%fixOddball = fixOddball_
        end if 

        if (present(allocK)) then
            this%allocK = allocK
        end if 
        
        ! Generate the physical space decomposition
        call decomp_info_init(nx_global, ny_global, nz_global, this%physical)
       
       
        select case (base_pencil_)
        case ("y") 
            ! Generate the spectral space decomposition
            if (allocated(this%f_yhat_in_yD)) deallocate(this%f_yhat_in_yD)
            if (allocated(this%f_xyhat_in_xD)) deallocate(this%f_xyhat_in_xD)
            if (allocated(this%f_xyzhat_in_zD)) deallocate(this%f_xyzhat_in_zD)
      
            call decomp_info_init(nx_global, ny_global/2 + 1, nz_global, this%spectral)
            allocate(this%f_yhat_in_yD (this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)),STAT=ierr) 
            allocate(this%f_xyhat_in_xD(this%spectral%xsz(1),this%spectral%xsz(2),this%spectral%xsz(3)),STAT=ierr) 
            allocate(this%f_xyzhat_in_zD(this%spectral%zsz(1),this%spectral%zsz(2),this%spectral%zsz(3)),STAT=ierr) 

            ! 1d transform along y (ydecomp assumed)  
            allocate(real_arr_2d(this%physical%ysz(1), this%physical%ysz(2)),STAT=ierr)
            allocate(cmplx_arr_2d(this%spectral%ysz(1),this%spectral%ysz(2)),STAT=ierr)

            if (ierr .ne. 0) then
                call decomp_2d_abort(304,"Could not allocate memory for fft_3")
            end if 

            ! Create plan for fwd transform in y
            n_sizeact = this%physical%ysz(2)
            n_sizeinput = this%physical%ysz(2)
            n_sizeoutput = this%spectral%ysz(2)
            n_howmany = this%physical%ysz(1)
            n_jump = this%physical%ysz(1)
            n_chunk = 1
            call dfftw_plan_many_dft_r2c(this%plan_r2c_y, 1, n_sizeact, &
                    n_howmany, real_arr_2d, n_sizeinput, n_jump, n_chunk, cmplx_arr_2d, n_sizeoutput, &
                    n_jump, n_chunk, this%fft_plan)

            ! Create plan for bwd transform in y
            n_sizeact = this%physical%ysz(2)
            n_sizeinput = this%spectral%ysz(2)
            n_sizeoutput = this%physical%ysz(2)
            n_howmany = this%physical%ysz(1)
            n_jump = this%physical%ysz(1)
            n_chunk = 1
            call dfftw_plan_many_dft_c2r(this%plan_c2r_y, 1, n_sizeact, &
                    n_howmany, cmplx_arr_2d, n_sizeinput, n_jump, n_chunk, real_arr_2d, n_sizeoutput, &
                    n_jump, n_chunk, this%fft_plan)
            
            
            deallocate(real_arr_2d, cmplx_arr_2d)


            ! Create plan for fwd transform in x (in place transform)
             call dfftw_plan_many_dft(this%plan_c2c_fwd_x, 1, this%spectral%xsz(1),&  
                    this%spectral%xsz(2)*this%spectral%xsz(3), this%f_xyhat_in_xD, this%spectral%xsz(1), 1, &
                    this%spectral%xsz(1), this%f_xyhat_in_xD, this%spectral%xsz(1), 1, this%spectral%xsz(1), &   
                    FFTW_FORWARD, FFTW_MEASURE)!this%fft_plan)   
                  
 
            ! Create plan for bwd transform in x (in place transform)        
             call dfftw_plan_many_dft(this%plan_c2c_bwd_x, 1, this%spectral%xsz(1),&  
                    this%spectral%xsz(2)*this%spectral%xsz(3), this%f_xyhat_in_xD, this%spectral%xsz(1), 1, &
                    this%spectral%xsz(1), this%f_xyhat_in_xD, this%spectral%xsz(1), 1, this%spectral%xsz(1), &   
                    FFTW_BACKWARD, FFTW_MEASURE)!, this%fft_plan)   
        
                
            ! Create plan for bwd transform in z (in place transform)
             call dfftw_plan_many_dft(this%plan_c2c_bwd_z, 1, this%spectral%zsz(3),&  
                    this%spectral%zsz(2)*this%spectral%zsz(1), this%f_xyzhat_in_zD, this%spectral%zsz(3), &
                    this%spectral%zsz(1)*this%spectral%zsz(2), 1, this%f_xyzhat_in_zD, this%spectral%zsz(3), &
                    this%spectral%zsz(1)*this%spectral%zsz(2), 1, FFTW_BACKWARD, this%fft_plan)!this%fft_plan)   
            
            ! Create plan for fwd transform in z (in place transform)
            call dfftw_plan_many_dft(this%plan_c2c_fwd_z, 1, this%spectral%zsz(3),&  
                this%spectral%zsz(2)*this%spectral%zsz(1), this%f_xyzhat_in_zD, this%spectral%zsz(3), &
                this%spectral%zsz(1)*this%spectral%zsz(2), 1, this%f_xyzhat_in_zD, this%spectral%zsz(3), &
                this%spectral%zsz(1)*this%spectral%zsz(2), 1, FFTW_FORWARD, this%fft_plan)!this%fft_plan)   
                
        case ("x")
            if (allocated(this%f_xhat_in_xD)) deallocate(this%f_xhat_in_xD)
            if (allocated(this%f_xyhat_in_yD)) deallocate(this%f_xyhat_in_yD)
            if (allocated(this%f_xyzhat_in_zD)) deallocate(this%f_xyzhat_in_zD)

            call decomp_info_init(nx_global/2 + 1, ny_global, nz_global, this%spectral)
            allocate(this%f_xhat_in_xD(this%spectral%xsz(1),this%spectral%xsz(2),this%spectral%xsz(3)))
            allocate(this%f_xyhat_in_yD(this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)))
            allocate(this%f_xyzhat_in_zD(this%spectral%zsz(1),this%spectral%zsz(2),this%spectral%zsz(3)),STAT=ierr) 

            ! Define real -> complex transform in x
            allocate (temp(this%physical%xsz(1),this%physical%xsz(2),this%physical%xsz(3)))
            
            call dfftw_plan_many_dft_r2c(this%plan_r2c_x, 1, this%physical%xsz(1), &
                    this%physical%xsz(2)*this%physical%xsz(3), temp, this%physical%xsz(1), 1, &
                    this%physical%xsz(1), this%f_xhat_in_xD, this%spectral%xsz(1), 1, this%spectral%xsz(1), &
                    this%fft_plan)
            
            ! Define complex -> real transform in x
            call dfftw_plan_many_dft_c2r(this%plan_c2r_x, 1, this%physical%xsz(1), &
                    this%physical%xsz(2)*this%physical%xsz(3), this%f_xhat_in_xD, this%spectral%xsz(1), 1, &
                    this%spectral%xsz(1), temp, this%physical%xsz(1), 1, this%physical%xsz(1), &
                    this%fft_plan)
            
            deallocate (temp)
        
            ! Define complex -> complex transforms in y
            allocate(cmplx_arr_2d(this%spectral%ysz(1),this%spectral%ysz(2)),STAT=ierr)
            allocate(cmplx_arr_2d_oop(this%spectral%ysz(1),this%spectral%ysz(2)),STAT=ierr)

            ! fwd transform in y (in place transform)
             call dfftw_plan_many_dft(this%plan_c2c_fwd_y, 1, this%spectral%ysz(2),&  
                    this%spectral%ysz(1), cmplx_arr_2d, this%spectral%ysz(2),this%spectral%ysz(1), &
                    1, cmplx_arr_2d, this%spectral%ysz(2), this%spectral%ysz(1),1, &   
                    FFTW_FORWARD, this%fft_plan)

            ! bwd transform in y (in place transform)
             call dfftw_plan_many_dft(this%plan_c2c_bwd_y, 1, this%spectral%ysz(2),&  
                    this%spectral%ysz(1), cmplx_arr_2d, this%spectral%ysz(2),this%spectral%ysz(1), &
                    1, cmplx_arr_2d, this%spectral%ysz(2), this%spectral%ysz(1),1, &   
                    FFTW_BACKWARD, this%fft_plan)
            
             ! bwd transform in y (out of place transform)
             call dfftw_plan_many_dft(this%plan_c2c_bwd_y_oop, 1, this%spectral%ysz(2),&  
                    this%spectral%ysz(1), cmplx_arr_2d, this%spectral%ysz(2),this%spectral%ysz(1), &
                    1, cmplx_arr_2d_oop, this%spectral%ysz(2), this%spectral%ysz(1),1, &   
                    FFTW_BACKWARD, this%fft_plan)

            deallocate (cmplx_arr_2d,cmplx_arr_2d_oop)

            allocate(dummy_in_z(this%spectral%zsz(1),this%spectral%zsz(2),this%spectral%zsz(3)))        
            ! Create plan for bwd transform in z (out-of-place transform)
            call dfftw_plan_many_dft(this%plan_c2c_bwd_z, 1, this%spectral%zsz(3),&  
                this%spectral%zsz(2)*this%spectral%zsz(1),dummy_in_z, this%spectral%zsz(3), &
                this%spectral%zsz(1)*this%spectral%zsz(2), 1, this%f_xyzhat_in_zD, this%spectral%zsz(3), &
                this%spectral%zsz(1)*this%spectral%zsz(2), 1, FFTW_BACKWARD, this%fft_plan)!this%fft_plan)   

            deallocate(dummy_in_z)
        
            ! Create plan for fwd transform in z (in place transform)
            call dfftw_plan_many_dft(this%plan_c2c_fwd_z, 1, this%spectral%zsz(3),&  
                this%spectral%zsz(2)*this%spectral%zsz(1), this%f_xyzhat_in_zD, this%spectral%zsz(3), &
                this%spectral%zsz(1)*this%spectral%zsz(2), 1, this%f_xyzhat_in_zD, this%spectral%zsz(3), &
                this%spectral%zsz(1)*this%spectral%zsz(2), 1, FFTW_FORWARD, this%fft_plan)!this%fft_plan)   


        case ("z")
            call decomp_info_init(nx_global, ny_global, nz_global/2 + 1, this%spectral)
            allocate(this%f_xyzhat_in_xD(this%spectral%xsz(1),this%spectral%xsz(2),this%spectral%xsz(3)))   
            allocate(this%f_zyhat_in_yD(this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)))
            call dfftw_plan_many_dft(this%plan_c2c_fwd_x, 1, this%spectral%xsz(1), &
                this%spectral%xsz(2)*this%spectral%xsz(3), this%f_xyzhat_in_xD, this%spectral%xsz(1), 1, &
                this%spectral%xsz(1), this%f_xyzhat_in_xD, this%spectral%xsz(1), 1, this%spectral%xsz(1), &
                FFTW_FORWARD, this%fft_plan)

            allocate (dummy_in_x(this%spectral%xsz(1),this%spectral%xsz(2),this%spectral%xsz(3)))
            call dfftw_plan_many_dft(this%plan_c2c_bwd_x, 1, this%spectral%xsz(1), &
                this%spectral%xsz(2)*this%spectral%xsz(3), this%f_xyzhat_in_xD, this%spectral%xsz(1), 1, &
                this%spectral%xsz(1), dummy_in_x, this%spectral%xsz(1), 1, this%spectral%xsz(1), &
                FFTW_BACKWARD, this%fft_plan)
            deallocate(dummy_in_x)

            ! Define complex -> complex transforms in y
            allocate(cmplx_arr_2d(this%spectral%ysz(1),this%spectral%ysz(2)),STAT=ierr)

            ! fwd transform in y (in place transform)
             call dfftw_plan_many_dft(this%plan_c2c_fwd_y, 1, this%spectral%ysz(2),&  
                    this%spectral%ysz(1), cmplx_arr_2d, this%spectral%ysz(2),this%spectral%ysz(1), &
                    1, cmplx_arr_2d, this%spectral%ysz(2), this%spectral%ysz(1),1, &   
                    FFTW_FORWARD, this%fft_plan)

            ! bwd transform in y (in place transform)
             call dfftw_plan_many_dft(this%plan_c2c_bwd_y, 1, this%spectral%ysz(2),&  
                    this%spectral%ysz(1), cmplx_arr_2d, this%spectral%ysz(2),this%spectral%ysz(1), &
                    1, cmplx_arr_2d, this%spectral%ysz(2), this%spectral%ysz(1),1, &   
                    FFTW_BACKWARD, this%fft_plan)

            deallocate (cmplx_arr_2d)

    
            allocate (this%f_zhat_in_zD(this%spectral%zsz(1),this%spectral%zsz(2),this%spectral%zsz(3)))
            allocate (temp(this%physical%zsz(1),this%physical%zsz(2),this%physical%zsz(3)))
            call dfftw_plan_many_dft_r2c(this%plan_r2c_z, 1, this%physical%zsz(3), &
                 this%physical%zsz(1)*this%physical%zsz(2), temp, this%physical%zsz(3), &
                 this%physical%zsz(1)*this%physical%zsz(2), 1, this%f_zhat_in_zD, this%spectral%zsz(3), &
                 this%spectral%zsz(1)*this%spectral%zsz(2), 1, this%fft_plan)

            call dfftw_plan_many_dft_c2r(this%plan_c2r_z, 1, this%physical%zsz(3), &
                 this%physical%zsz(1)*this%physical%zsz(2), this%f_zhat_in_zD, this%spectral%zsz(3), &
                 this%spectral%zsz(1)*this%spectral%zsz(2), 1, temp, this%physical%zsz(3), &
                 this%physical%zsz(1)*this%physical%zsz(2), 1, this%fft_plan)
            deallocate(temp)
 
        end select 

                
         this%base_pencil = base_pencil_
         this%normfactor   = 1._rkind/(real(nx_global*ny_global*nz_global,rkind))
         this%normfactor2d = 1._rkind/(real(nx_global*ny_global,rkind))

       
         if (this%allocK) then 
            ! Make wavenumbers
            allocate (k1_in_x(nx_global           ,this%spectral%xsz(2),this%spectral%xsz(3)))
            allocate (k2_in_y(this%spectral%ysz(1),ny_global           ,this%spectral%ysz(3)))
            allocate (k3_in_z(this%spectral%zsz(1),this%spectral%zsz(2),nz_global           ))

            ! Generate full k1
            do k = 1,this%spectral%xsz(3)
               do j = 1,this%spectral%xsz(2)
                   k1_in_x(:,j,k) = GetWaveNums(nx_global,dx)
               end do 
            end do
            
            ! Generate full k2
            do k = 1,this%spectral%ysz(3)
               do i = 1,this%spectral%ysz(1)
                   k2_in_y(i,:,k) = GetWaveNums(ny_global,dy)
               end do 
            end do
               
            ! Generate full k3
            do j = 1,this%spectral%zsz(2)
               do i = 1,this%spectral%zsz(1)
                   k3_in_z(i,j,:) = GetWaveNums(nz_global,dz)
               end do
            end do
               
            if (this%fixOddball) then
               k1_in_x(nx_global/2+1,:,:) = zero 
               k2_in_y(:,ny_global/2+1,:) = zero
               k3_in_z(:,:,nz_global/2+1) = zero
            end if 
              
            select case (base_pencil_)   
            case ("y") 
               ! All wavenumbers are stored in y-decomp
               allocate(this%k1(this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)))
               allocate(this%k2(this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)))
               allocate(this%k3(this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)))
               allocate(this%kabs_sq(this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)))
               
               this%k2 = k2_in_y(1:this%spectral%ysz(1),1:this%spectral%ysz(2),1:this%spectral%ysz(3))

               ! Get k1 from x -> y
               call transpose_x_to_y(k1_in_x,this%k1,this%spectral)


               ! Get k3 from z -> y
               call transpose_z_to_y(k3_in_z,this%k3,this%spectral)

            case ("x")
               ! All wavenumbers are stored in z-decomp
               allocate(this%k1(this%spectral%zsz(1),this%spectral%zsz(2),this%spectral%zsz(3)))
               allocate(this%k2(this%spectral%zsz(1),this%spectral%zsz(2),this%spectral%zsz(3)))
               allocate(this%k3(this%spectral%zsz(1),this%spectral%zsz(2),this%spectral%zsz(3)))
               allocate(this%kabs_sq(this%spectral%zsz(1),this%spectral%zsz(2),this%spectral%zsz(3)))

               allocate(temp(this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)))
               
               this%k3 = k3_in_z

               ! Get k1 from x -> y -> z
               call transpose_x_to_y(k1_in_x(1:this%spectral%xsz(1),1:this%spectral%xsz(2),1:this%spectral%xsz(3)),temp,this%spectral)
               call transpose_y_to_z(temp,this%k1,this%spectral)

               ! Get k2 from y -> z
               call transpose_y_to_z(k2_in_y,this%k2,this%spectral)

               deallocate (temp)

            case ("z")
               ! All wavenumbers are stored in z-decomp
               allocate(this%k1(this%spectral%xsz(1),this%spectral%xsz(2),this%spectral%xsz(3)))
               allocate(this%k2(this%spectral%xsz(1),this%spectral%xsz(2),this%spectral%xsz(3)))
               allocate(this%k3(this%spectral%xsz(1),this%spectral%xsz(2),this%spectral%xsz(3)))
               allocate(this%kabs_sq(this%spectral%xsz(1),this%spectral%xsz(2),this%spectral%xsz(3)))

               allocate(temp(this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)))
               
               this%k1 = k1_in_x

               ! Get k3 from z -> y -> x
               call transpose_z_to_y(k3_in_z,temp,this%spectral)
               call transpose_y_to_x(temp,this%k3,this%spectral)

               ! Get k2 from y -> x
               call transpose_y_to_x(k2_in_y,this%k2,this%spectral)

               deallocate (temp)
               
            end select 

            deallocate (k1_in_x, k2_in_y, k3_in_z)

            this%kabs_sq = this%k1*this%k1 + this%k2*this%k2 + this%k3*this%k3

         end if 
         this%initialized = .true.

         ierr = 0

    end function


    subroutine destroy(this)
        class(fft_3d) , intent(inout) :: this

        if (this%initialized) then
            if (allocated(this%f_yhat_in_yD)) deallocate(this%f_yhat_in_yD)
            if (allocated(this%f_xhat_in_xD)) deallocate(this%f_xhat_in_xD)
            if (allocated(this%f_zhat_in_zD)) deallocate(this%f_zhat_in_zD)
            if (allocated(this%f_xyhat_in_yD)) deallocate(this%f_xyhat_in_yD)
            if (allocated(this%f_zyhat_in_yD)) deallocate(this%f_zyhat_in_yD)
            if (allocated(this%f_xyzhat_in_zD)) deallocate(this%f_xyzhat_in_zD)
            if (allocated(this%f_xyzhat_in_xD)) deallocate(this%f_xyzhat_in_xD)

            if (allocated(this%k1       )) deallocate( this%k1     )       
            if (allocated(this%k2       )) deallocate( this%k2     ) 
            if (allocated(this%k3       )) deallocate( this%k3     ) 
            if (allocated(this%kabs_sq  )) deallocate( this%kabs_sq) 
       
            select case (this%base_pencil)
            case ("y") 
                call dfftw_destroy_plan (this%plan_r2c_y    ) 
                call dfftw_destroy_plan (this%plan_c2c_fwd_x)
                call dfftw_destroy_plan (this%plan_c2c_bwd_x)
                call dfftw_destroy_plan (this%plan_c2c_fwd_z)
                call dfftw_destroy_plan (this%plan_c2c_bwd_z)
                call dfftw_destroy_plan (this%plan_c2r_y    )
            case ("x")
                call dfftw_destroy_plan (this%plan_c2c_fwd_y) 
                call dfftw_destroy_plan (this%plan_c2c_bwd_y) 
                call dfftw_destroy_plan (this%plan_c2c_fwd_z) 
                call dfftw_destroy_plan (this%plan_c2c_bwd_z) 
                call dfftw_destroy_plan (this%plan_r2c_x) 
                call dfftw_destroy_plan (this%plan_c2r_x) 
            case ("z")
                call dfftw_destroy_plan (this%plan_c2c_fwd_y) 
                call dfftw_destroy_plan (this%plan_c2c_bwd_y) 
                call dfftw_destroy_plan (this%plan_c2c_fwd_x) 
                call dfftw_destroy_plan (this%plan_c2c_bwd_x) 
                call dfftw_destroy_plan (this%plan_r2c_z) 
                call dfftw_destroy_plan (this%plan_c2r_z) 
            end select
 
            call decomp_info_finalize(this%spectral)
            call decomp_info_finalize(this%physical)

            this%initialized = .false.
        end if 

    end subroutine


    subroutine fft3_y2y(this,input, output)
        class(fft_3d), intent(inout) :: this
        real   (rkind), intent(in ), dimension(this%physical%ysz(1),this%physical%ysz(2),this%physical%ysz(3)) :: input
        complex(rkind), intent(out), dimension(this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)) :: output
        integer :: k

            ! First take transform in y  
            do k = 1,this%physical%ysz(3)
                call dfftw_execute_dft_r2c(this%plan_r2c_y, input(:,:,k), output(:,:,k))  
            end do 

            ! Now transform f_hat from y -> x
            call transpose_y_to_x(output, this%f_xyhat_in_xD, this%spectral)
            
            ! Now take in place transform along x (complex transform)
            call dfftw_execute_dft(this%plan_c2c_fwd_x, this%f_xyhat_in_xD, this%f_xyhat_in_xD)  
      
            ! Now tranform back from x -> y
            call transpose_x_to_y(this%f_xyhat_in_xD,output,this%spectral)
   
            ! Now transpose from y -> z
            call transpose_y_to_z(output,this%f_xyzhat_in_zD,this%spectral)

            ! Now take in place transform along z (complex transform)
            call dfftw_execute(this%plan_c2c_fwd_z,this%f_xyzhat_in_zD,this%f_xyzhat_in_zD)
      
            ! Now transpose the result back to y
            call transpose_z_to_y(this%f_xyzhat_in_zD,output,this%spectral)

    end subroutine 



    subroutine ifft3_y2y(this,input,output)
        class(fft_3d), intent(inout) :: this
        complex(rkind), intent(in ), dimension(this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)) :: input
        real   (rkind), intent(out), dimension(this%physical%ysz(1),this%physical%ysz(2),this%physical%ysz(3)) :: output
        complex(rkind), dimension(this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)) :: tmp_arr
        integer :: k 
  
        ! First transpose input y -> z
        call transpose_y_to_z(input,this%f_xyzhat_in_zD,this%spectral)

        ! Then take in place transform in z
        call dfftw_execute(this%plan_c2c_bwd_z,this%f_xyzhat_in_zD,this%f_xyzhat_in_zD)

        ! Now double-transform z -> x
        call transpose_z_to_y(this%f_xyzhat_in_zD,tmp_arr,this%spectral)
        call transpose_y_to_x(tmp_arr, this%f_xyhat_in_xD, this%spectral)

        ! Then take in place transform along x (complex transform)
        call dfftw_execute_dft(this%plan_c2c_bwd_x, this%f_xyhat_in_xD, this%f_xyhat_in_xD)  

        ! Now transform this output from x -> y
        call transpose_x_to_y(this%f_xyhat_in_xD,tmp_arr, this%spectral)

        ! Now take c2r transform along y 
        do k = 1,this%physical%ysz(3)
            call dfftw_execute_dft_c2r(this%plan_c2r_y, tmp_arr(:,:,k), output(:,:,k))
        end do 

        ! Now normalize the output
        output = output*this%normfactor
        
        ! Done 
     end subroutine

    subroutine fft3_x2z(this,input,output)
        class(fft_3d), intent(inout) :: this
        real(rkind), dimension(this%physical%xsz(1),this%physical%xsz(2),this%physical%xsz(3)), intent(in) :: input
        complex(rkind), dimension(this%spectral%zsz(1),this%spectral%zsz(2),this%spectral%zsz(3)), intent(out) :: output

        integer :: k

        ! First get the x tranform (r2c)
        call dfftw_execute_dft_r2c(this%plan_r2c_x, input, this%f_xhat_in_xD)  

        ! Now transpose from x-> y
        call transpose_x_to_y(this%f_xhat_in_xD,this%f_xyhat_in_yD,this%spectral)
        
        ! The take fwd transform in y (c2c, inplace)
        do k = 1,this%spectral%ysz(3)
            call dfftw_execute_dft(this%plan_c2c_fwd_y, this%f_xyhat_in_yD(:,:,k), this%f_xyhat_in_yD(:,:,k))  
        end do 

        ! Now transpose from y-> z
        call transpose_y_to_z(this%f_xyhat_in_yD,output,this%spectral) 

        ! Now take fwd transform in z (c2c, inplace)
        call dfftw_execute_dft(this%plan_c2c_fwd_z, output, output)  

        ! Done 
    end subroutine

    subroutine fft2_x2y(this,input,output)
        class(fft_3d), intent(inout) :: this
        real(rkind), dimension(this%physical%xsz(1),this%physical%xsz(2),this%physical%xsz(3)), intent(in) :: input
        complex(rkind), dimension(size(this%f_xyhat_in_yD,1),size(this%f_xyhat_in_yD,2),size(this%f_xyhat_in_yD,3)), intent(out) :: output
       
        integer :: k

        ! First get the x tranform (r2c)
        call dfftw_execute_dft_r2c(this%plan_r2c_x, input, this%f_xhat_in_xD)  

        ! Now transpose from x-> y
        call transpose_x_to_y(this%f_xhat_in_xD,output,this%spectral)
       
        ! The take fwd transform in y (c2c, inplace)
        do k = 1,this%spectral%ysz(3)
            call dfftw_execute_dft(this%plan_c2c_fwd_y, output(:,:,k), output(:,:,k))  
        end do 

    end subroutine

    subroutine ifft2_y2x(this,input,output,setOddBall)
        class(fft_3d), intent(inout) :: this
        complex(rkind), dimension(size(this%f_xyhat_in_yD,1),size(this%f_xyhat_in_yD,2),size(this%f_xyhat_in_yD,3)), intent(in) :: input
        real(rkind), dimension(this%physical%xsz(1),this%physical%xsz(2),this%physical%xsz(3)), intent(out) :: output
        
        integer :: k 
        logical :: setOddBall


        ! Then transform in y (c2c, out of place)
        do k = 1,this%spectral%ysz(3)
            call dfftw_execute_dft(this%plan_c2c_bwd_y_oop, input(:,:,k), this%f_xyhat_in_yD(:,:,k))  
        end do 

        ! Then transpose from y-> x
        call transpose_y_to_x(this%f_xyhat_in_yD,this%f_xhat_in_xD,this%spectral)
        
        if (setOddBall) then
            this%f_xhat_in_xD(this%physical%xsz(1)/2+1,:,:) = zero
        end if 

        ! Then transform in x (c2r transform)
        call dfftw_execute_dft_c2r(this%plan_c2r_x, this%f_xhat_in_xD, output)
        
        ! Normalize the tranform 
        output = output*this%normfactor2d

    end subroutine        
   
 
    subroutine ifft3_z2x(this,input,output)
        class(fft_3d), intent(inout) :: this
        complex(rkind), dimension(this%spectral%zsz(1),this%spectral%zsz(2),this%spectral%zsz(3)), intent(in) :: input
        real(rkind), dimension(this%physical%xsz(1),this%physical%xsz(2),this%physical%xsz(3)), intent(out) :: output
        
        integer :: k 
    
        call dfftw_execute_dft(this%plan_c2c_bwd_z, input, this%f_xyzhat_in_zD)  

        ! Then transpose from z -> y
        call transpose_z_to_y(this%f_xyzhat_in_zD,this%f_xyhat_in_yD,this%spectral)

        ! Then transform in y (c2c, inplace)
        do k = 1,this%spectral%ysz(3)
            call dfftw_execute_dft(this%plan_c2c_bwd_y, this%f_xyhat_in_yD(:,:,k), this%f_xyhat_in_yD(:,:,k))  
        end do 

        ! Then transpose from y-> x
        call transpose_y_to_x(this%f_xyhat_in_yD,this%f_xhat_in_xD,this%spectral)

        ! Then transform in x (c2r transform)
        call dfftw_execute_dft_c2r(this%plan_c2r_x, this%f_xhat_in_xD, output)

        ! Normalize the tranform 
        output = output*this%normfactor

    end subroutine
    
    subroutine fft3_z2x(this,input,output)
        class(fft_3d), intent(inout) :: this
        real(rkind), dimension(this%physical%zsz(1),this%physical%zsz(2),this%physical%zsz(3)), intent(in) :: input
        complex(rkind), dimension(this%spectral%xsz(1),this%spectral%xsz(2),this%spectral%xsz(3)), intent(out) :: output
        integer :: k 

        ! Take the z derivative

        call dfftw_execute_dft_r2c(this%plan_r2c_z, input, this%f_zhat_in_zD)  
        
        call transpose_z_to_y(this%f_zhat_in_zD,this%f_zyhat_in_yD,this%spectral)
        
        ! The take fwd transform in y (c2c, inplace)
        do k = 1,this%spectral%ysz(3)
            call dfftw_execute_dft(this%plan_c2c_fwd_y, this%f_zyhat_in_yD(:,:,k), this%f_zyhat_in_yD(:,:,k))  
        end do 

        call transpose_y_to_x(this%f_zyhat_in_yD,output,this%spectral)
        
        call dfftw_execute_dft(this%plan_c2c_fwd_x, output, output)  
        
    end subroutine


    subroutine ifft3_x2z(this,input,output)
        class(fft_3d), intent(inout) :: this
        complex(rkind), dimension(this%spectral%xsz(1),this%spectral%xsz(2),this%spectral%xsz(3)), intent(out) :: input
        real(rkind), dimension(this%physical%zsz(1),this%physical%zsz(2),this%physical%zsz(3)), intent(out) :: output
        integer :: k 

        call dfftw_execute_dft(this%plan_c2c_bwd_x, input, this%f_xyzhat_in_xD)  

        call transpose_x_to_y(this%f_xyzhat_in_xD,this%f_zyhat_in_yD,this%spectral)
        
        ! The take fwd transform in y (c2c, inplace)
        do k = 1,this%spectral%ysz(3)
            call dfftw_execute_dft(this%plan_c2c_bwd_y, this%f_zyhat_in_yD(:,:,k), this%f_zyhat_in_yD(:,:,k))  
        end do 
        
        call transpose_y_to_z(this%f_zyhat_in_yD,this%f_zhat_in_zD,this%spectral)

        
        ! Then transform in x (c2r transform)
        call dfftw_execute_dft_c2r(this%plan_c2r_z, this%f_zhat_in_zD, output)

        ! Normalize the tranform 
        output = output*this%normfactor

    end subroutine 

    subroutine alloc_3d_Real(this,arr_out) 
        class(fft_3d), intent(in) :: this
        real(rkind), dimension(:,:,:),allocatable, intent(inout) :: arr_out

        if (this%initialized) then
            if (allocated(arr_out)) deallocate(arr_out)
            select case (this%base_pencil)
            case ("y")
                allocate(arr_out(this%spectral%ysz(1), this%spectral%ysz(2), this%spectral%ysz(3)))
            case ("x")
                allocate(arr_out(this%spectral%zsz(1), this%spectral%zsz(2), this%spectral%zsz(3)))
            case ("z")
                allocate(arr_out(this%spectral%xsz(1), this%spectral%xsz(2), this%spectral%xsz(3)))
            end select 
            arr_out = 0._rkind
        else
            call decomp_2d_abort(305,"The fft_3d type has not been initialized")
        end if 

    end subroutine
    
    subroutine alloc_3d_Cmplx(this,arr_out) 
        class(fft_3d), intent(in) :: this
        complex(rkind), dimension(:,:,:),allocatable, intent(inout) :: arr_out

        if (this%initialized) then
            if (allocated(arr_out)) deallocate(arr_out)
            select case (this%base_pencil)
            case ("y")
                allocate(arr_out(this%spectral%ysz(1), this%spectral%ysz(2), this%spectral%ysz(3)))
            case ("x")
                allocate(arr_out(this%spectral%zsz(1), this%spectral%zsz(2), this%spectral%zsz(3)))
            case ("z")
                allocate(arr_out(this%spectral%xsz(1), this%spectral%xsz(2), this%spectral%xsz(3)))
            end select 
            arr_out = 0._rkind
        else
            call decomp_2d_abort(305,"The fft_3d type has not been initialized")
        end if 

    end subroutine
   
    subroutine alloc_input(this,arr_out)
        class(fft_3d), intent(in) :: this
        real(rkind), dimension(:,:,:), allocatable, intent(inout) :: arr_out

        if (this%initialized) then
            if (allocated(arr_out)) deallocate(arr_out)
            select case (this%base_pencil)
            case ("y")
                allocate(arr_out(this%physical%ysz(1), this%physical%ysz(2), this%physical%ysz(3)))
            case ("x")
                allocate(arr_out(this%physical%xsz(1), this%physical%xsz(2), this%physical%xsz(3)))
            case ("z")
                allocate(arr_out(this%physical%zsz(1), this%physical%zsz(2), this%physical%zsz(3)))
            end select 
        else
            call decomp_2d_abort(305,"The fft_3d type has not been initialized")
        end if 

    end subroutine  
    
    subroutine alloc_4d(this,arr_out,dim4) 
        class(fft_3d), intent(in) :: this
        integer, intent(in)       :: dim4
        complex(rkind), dimension(:,:,:,:),allocatable, intent(inout) :: arr_out

        if (this%initialized) then
            if (allocated(arr_out)) deallocate(arr_out)
            select case (this%base_pencil)
            case ("y")
                allocate(arr_out(this%spectral%ysz(1), this%spectral%ysz(2), this%spectral%ysz(3),dim4))
            case ("x")
                allocate(arr_out(this%spectral%zsz(1), this%spectral%zsz(2), this%spectral%zsz(3),dim4))
            case ("z")
                allocate(arr_out(this%spectral%xsz(1), this%spectral%xsz(2), this%spectral%xsz(3),dim4))
            end select 
            arr_out = 0._rkind
        else
            call decomp_2d_abort(305,"The fft_3d type has not been initialized")
        end if 

    end subroutine
    
    
    pure function GetWaveNums(nx,dx) result(k)
        use constants, only: pi, two
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

end module 
