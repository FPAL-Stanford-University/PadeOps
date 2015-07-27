module fft_2d_stuff

    use kind_parameters, only : rkind
    use decomp_2d
   
    implicit none 
    private
    public :: fft_2d

    include "fftw3.f"
    
    type fft_2d
        private
        type(decomp_info)           :: physical
        type(decomp_info),public    :: spectral
        
        character(len=1)  :: base_pencil
        real (rkind) :: normfactor

        complex(rkind), dimension(:,:,:), allocatable :: f_yhat_in_yD
        complex(rkind), dimension(:,:,:), allocatable :: f_xyhat_in_xD

        integer(kind=8) :: plan_r2c_y
        integer(kind=8) :: plan_c2c_fwd_x
        integer(kind=8) :: plan_c2c_bwd_x
        integer(kind=8) :: plan_c2r_y

        logical :: initialized = .false.

        integer :: fft_plan = FFTW_MEASURE
        contains
            procedure :: init
            procedure :: fft2
            procedure :: ifft2
            procedure :: destroy 

    end type 


contains

    function init(this,nx_global,ny_global,nz_global,base_pencil_, exhaustive) result(ierr)
        class(fft_2d), intent(inout) :: this
        integer, intent(in) :: nx_global, ny_global, nz_global
        character(len=1), intent(in) :: base_pencil_
        logical, optional :: exhaustive
        integer :: ierr
        integer, dimension(2) :: dims, dummy_periods, dummy_coords
        real(rkind), dimension(:,:), allocatable :: real_arr_2d
        complex(rkind), dimension(:,:), allocatable :: cmplx_arr_2d
        integer :: n_sizeact, n_sizeinput, n_sizeoutput, n_chunk, n_howmany, n_jump

        ! Get the decompotion performed before calling this function
        call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, & 
                dims, dummy_periods, dummy_coords, ierr) 


        if (base_pencil_ .ne. "y") then
            call decomp_2d_abort(300,"Currently only Y - baseline pencils are &
            supported for 2D ffts. FFT is always taken in XY plane for each Z &
            location. Input array must be provided in Y decomposition setting.")
        end if

        if ( mod(ny_global,2) .ne. 0) then
            call decomp_2d_abort(301,"Only even number of data are supported &
            in Y direction. NY_GLOBAL must be an even number.") 
        end if 
        
        if (((dims(1) == 1) .or. (dims(2) == 1)).and. (nproc > 1)) then
            !call decomp_2d_abort(303,"Only 2d decompositions supported. If the &
            !auto-tuner gives an effectively 1d decomposition, overwrite it to &
            !the next best 2d decomposition")
        end if 

        if (present(exhaustive)) then
            if (exhaustive) this%fft_plan = FFTW_EXHAUSTIVE
        end if 
        
        ! Generate the physical space decomposition
        call decomp_info_init(nx_global, ny_global, nz_global, this%physical)
        
        ! Generate the spectral space decompisition
        if (allocated(this%f_yhat_in_yD)) deallocate(this%f_yhat_in_yD)
        if (allocated(this%f_xyhat_in_xD)) deallocate(this%f_xyhat_in_xD)
      
        call decomp_info_init(nx_global, ny_global/2 + 1, nz_global, this%spectral)
        allocate(this%f_yhat_in_yD (this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)),STAT=ierr) 
        allocate(this%f_xyhat_in_xD(this%spectral%xsz(1),this%spectral%xsz(2),this%spectral%xsz(3)),STAT=ierr) 

        ! 1d transform along y (ydecomp assumed)  
        allocate(real_arr_2d(this%physical%ysz(1), this%physical%ysz(2)),STAT=ierr)
        allocate(cmplx_arr_2d(this%spectral%ysz(1),this%spectral%ysz(2)),STAT=ierr)

        if (ierr .ne. 0) then
            call decomp_2d_abort(304,"Could not allocate memory for fft_2")
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
                FFTW_FORWARD, this%fft_plan)   
               
        ! Create plan for bwd transform in x (in place transform)        
         call dfftw_plan_many_dft(this%plan_c2c_bwd_x, 1, this%spectral%xsz(1),&  
                this%spectral%xsz(2)*this%spectral%xsz(3), this%f_xyhat_in_xD, this%spectral%xsz(1), 1, &
                this%spectral%xsz(1), this%f_xyhat_in_xD, this%spectral%xsz(1), 1, this%spectral%xsz(1), &   
                FFTW_BACKWARD, this%fft_plan)   
        
                
         this%base_pencil = base_pencil_
         this%normfactor = 1._rkind
         this%initialized = .true.

         ierr = 0

    end function


    subroutine destroy(this)
        class(fft_2d) , intent(inout) :: this

        if (this%initialized) then
            if (allocated(this%f_yhat_in_yD)) deallocate(this%f_yhat_in_yD)
            if (allocated(this%f_xyhat_in_xD)) deallocate(this%f_xyhat_in_xD)
        

        call decomp_info_finalize(this%spectral)
        call decomp_info_finalize(this%physical)


            this%initialized = .false.
        end if 

    end subroutine


    subroutine fft2(this,input, output)
        class(fft_2d), intent(inout) :: this
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
    
       ! Done      
    end subroutine 



    subroutine ifft2(this,input,output)
        class(fft_2d), intent(inout) :: this
        complex(rkind), intent(in ), dimension(this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)) :: input
        real   (rkind), intent(out), dimension(this%physical%ysz(1),this%physical%ysz(2),this%physical%ysz(3)) :: output
        complex(rkind), dimension(this%spectral%ysz(1),this%spectral%ysz(2),this%spectral%ysz(3)) :: tmp_arr
        integer :: k 
        
        ! First transform input from y -> x
        call transpose_y_to_x(input, this%f_xyhat_in_xD, this%spectral)

        ! Now take in place transform along x (complex transform)
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


end module 
