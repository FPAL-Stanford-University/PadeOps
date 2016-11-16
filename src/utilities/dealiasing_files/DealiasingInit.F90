subroutine init_dealiasing_stuff(this,nx,ny,nz)
    use constants, only: two, three, zero, one
    class(fft_3d), intent(inout), target :: this
    integer, intent(in) :: nx, ny, nz
    integer :: nsize, howmany, inp_n_embed, out_n_embed, inp_stride, out_stride, inp_dist, out_dist
    real(rkind), dimension(:,:,:), allocatable, target :: tmp
    complex(rkind), pointer, dimension(:,:,:) :: c_arr_in, c_arr_out
    real(rkind), pointer, dimension(:,:,:) :: r_arr_in, r_arr_out

    this%nxF = ceiling((three*real(nx))/two) 
    this%nyF = ceiling((three*real(ny))/two)
    call decomp_info_init(this%nxF/2+1,ny,nz,this%D1)
    allocate(this%d11x(this%D1%xsz(1),this%D1%xsz(2),this%D1%xsz(3)))
    allocate(this%d11y(this%D1%ysz(1),this%D1%ysz(2),this%D1%ysz(3)))
    this%d11x = zero 

    call decomp_info_init(this%nxF/2+1,this%nyF,nz, this%D2)
    allocate(this%d21y(this%D2%ysz(1),this%D2%ysz(2),this%D2%ysz(3)))
    allocate(this%d21x(this%D2%xsz(1),this%D2%xsz(2),this%D2%xsz(3)))
   
    call decomp_info_init(this%nxF,this%nyF, nz, this%PhysicalUP)

    call decomp_info_init(nx/2+1,this%nyF, nz, this%D3)
    allocate(this%d31y(this%D3%ysz(1),this%D3%ysz(2),this%D3%ysz(3)))
    allocate(this%d31x(this%D3%xsz(1),this%D3%xsz(2),this%D3%xsz(3)))

    !call decomp_info_init(nx/2+1,ny,nz, this%D4) 

    this%D4 => this%spectral
    allocate(this%d41x(this%D4%xsz(1),this%D4%xsz(2),this%D4%xsz(3)))
    allocate(this%d41y(this%D4%ysz(1),this%D4%ysz(2),this%D4%ysz(3)))


    ! Plan for bwd_d21x (c2r)
    allocate(tmp(this%PhysicalUp%xsz(1),this%PhysicalUp%xsz(2),this%PhysicalUp%xsz(3)))
    nsize = this%PhysicalUp%xsz(1) 
    howmany = this%PhysicalUP%xsz(2)*this%PhysicalUP%xsz(3)
    c_arr_in => this%d21x
    inp_n_embed = this%D2%xsz(1)
    inp_stride = 1
    inp_dist = this%D2%xsz(1)
    r_arr_out => tmp 
    out_n_embed = this%PhysicalUP%xsz(1)
    out_stride = 1
    out_dist = this%PhysicalUP%xsz(1)
    
    call dfftw_plan_many_dft_c2r(this%plan_bwd_d21x,1,nsize,howmany,c_arr_in, inp_n_embed, inp_stride, &
                                inp_dist, r_arr_out , out_n_embed, out_stride, out_dist, this%fft_plan) 
    deallocate(tmp)

    ! Plan for bwd_d21y (inplace)
    nsize = this%D2%ysz(2) 
    howmany = this%D2%ysz(1)
    c_arr_in => this%d21y
    inp_n_embed = this%D2%ysz(2)
    inp_stride = this%D2%ysz(1)
    inp_dist = 1
    c_arr_out => this%d21y 
    out_n_embed = this%D2%ysz(2)
    out_stride = this%D2%ysz(1)
    out_dist = 1
    
    call dfftw_plan_many_dft(this%plan_bwd_d21y,1,nsize,howmany,c_arr_in, inp_n_embed, inp_stride, &
                                inp_dist, c_arr_out , out_n_embed, out_stride, out_dist, FFTW_BACKWARD, &
                                this%fft_plan) 
     
    ! Plan for fwd_d11y (inplace)                            
    nsize = this%D1%ysz(2) 
    howmany = this%D1%ysz(1)
    c_arr_in => this%d11y
    inp_n_embed = this%D1%ysz(2)
    inp_stride = this%D1%ysz(1)
    inp_dist = 1
    c_arr_out => this%d11y 
    out_n_embed = this%D1%ysz(2)
    out_stride = this%D1%ysz(1)
    out_dist = 1
    
    call dfftw_plan_many_dft(this%plan_fwd_d11y,1,nsize,howmany,c_arr_in, inp_n_embed, inp_stride, &
                                inp_dist, c_arr_out , out_n_embed, out_stride, out_dist, FFTW_FORWARD, &
                                this%fft_plan) 
    

    ! Plan for fwd_d21x (r2c)
    allocate(tmp(this%PhysicalUP%xsz(1),this%PhysicalUP%xsz(2),this%PhysicalUP%xsz(3)))
    nsize = this%PhysicalUP%xsz(1) 
    howmany = this%PhysicalUP%xsz(2)*this%PhysicalUP%xsz(3)
    r_arr_in => tmp
    inp_n_embed = this%PhysicalUP%xsz(1)
    inp_stride = 1
    inp_dist = this%PhysicalUP%xsz(1)
    c_arr_out => this%d21x
    out_n_embed = this%D2%xsz(1)
    out_stride = 1
    out_dist = this%D2%xsz(1)
    
    call dfftw_plan_many_dft_r2c(this%plan_fwd_d21x,1,nsize,howmany,r_arr_in, inp_n_embed, inp_stride, &
                                inp_dist, c_arr_out , out_n_embed, out_stride, out_dist, this%fft_plan) 
    deallocate(tmp)


    ! Plan for fwd_d31y (inplace)
    nsize = this%D3%ysz(2) 
    howmany = this%D3%ysz(1)
    c_arr_in => this%d31y
    inp_n_embed = this%D3%ysz(2)
    inp_stride = this%D3%ysz(1)
    inp_dist = 1
    c_arr_out => this%d31y 
    out_n_embed = this%D3%ysz(2)
    out_stride = this%D3%ysz(1)
    out_dist = 1
    
    call dfftw_plan_many_dft(this%plan_fwd_d31y,1,nsize,howmany,c_arr_in, inp_n_embed, inp_stride, &
                                inp_dist, c_arr_out , out_n_embed, out_stride, out_dist, FFTW_FORWARD, &
                                this%fft_plan) 
    
     
    ! Plan for bwd_d41y (inplace)
    !nsize = this%D4%ysz(2) 
    !howmany = this%D4%ysz(1)
    !c_arr_in => this%d41y
    !inp_n_embed = this%D4%ysz(2)
    !inp_stride = this%D4%ysz(1)
    !inp_dist = 1
    !c_arr_out => this%d41y 
    !out_n_embed = this%D4%ysz(2)
    !out_stride = this%D4%ysz(1)
    !out_dist = 1
    !
    !call dfftw_plan_many_dft(this%plan_fwd_d31y,1,nsize,howmany,c_arr_in, inp_n_embed, inp_stride, &
    !                            inp_dist, c_arr_out , out_n_embed, out_stride, out_dist, FFTW_FORWARD, &
    !                            this%fft_plan) 
     
    
    
    
    this%dealiasNormFact = one/(real(this%nxF, rkind) * real(this%nyF,rkind))                    

end subroutine


subroutine alloc_upsampledArr(this,arr2alloc)
    class(fft_3d), intent(in) :: this
    real(rkind), dimension(:,:,:), allocatable, intent(out) :: arr2alloc

    if (allocated(arr2alloc)) deallocate(arr2alloc)
    allocate(arr2alloc(this%PhysicalUP%xsz(1),this%PhysicalUP%xsz(2), this%PhysicalUP%xsz(3)))

end subroutine
