    subroutine upsample(this,input,output)
        use constants, only: three, two
        class(fft_3d), intent(inout) :: this
        real(rkind), dimension(this%physical%xsz(1),this%physical%xsz(2),this%physical%xsz(3)), intent(in) :: input
        real(rkind), dimension(this%physicalUP%xsz(1),this%physicalUP%xsz(2),this%physicalUP%xsz(3)), intent(out) :: output
        integer :: nx, ny, k
        
        
        nx = this%Physical%xsz(1); ny = this%Physical%ysz(2)

        call dfftw_execute_dft_r2c(this%plan_r2c_x, input, this%f_xhat_in_xD)  
        this%d11x(1:nx/2+1,:,:) = (three/two)*this%f_xhat_in_xD 
        this%d11x(nx/2+1:this%nxF/2+1,:,:) = zero

        call transpose_x_to_y(this%d11x, this%d11y,this%D1)

        do k = 1,this%D1%ysz(3)
            call dfftw_execute_dft(this%plan_fwd_d11y, this%d11y(:,:,k), this%d11y(:,:,k))  
            this%d21y(:,1:ny/2+1,k) = this%d11y(:,1:ny/2+1,k)
            this%d21y(:,ny/2+2:ny+1,k) = zero
            this%d21y(:,ny+2:this%nyF,k) = this%d11y(:,ny/2+2:ny,k)
        end do
        this%d21y = (three/two)*this%d21y

        do k = 1,this%D2%ysz(3)
            call dfftw_execute_dft(this%plan_bwd_d21y, this%d21y(:,:,k), this%d21y(:,:,k))  
        end do 
        
        call transpose_y_to_x(this%d21y,this%d21x,this%D2)
       
        call dfftw_execute_dft_c2r(this%plan_bwd_d21x, this%d21x, output) 

        output = output*this%dealiasNormFact

    end subroutine


    subroutine downsample(this,input,output)
        use constants, only: three, two
        class(fft_3d), intent(inout) :: this
        real(rkind), dimension(this%physicalUP%xsz(1),this%physicalUP%xsz(2),this%physicalUP%xsz(3)), intent(in) :: input
        real(rkind), dimension(this%physical%xsz(1),this%physical%xsz(2),this%physical%xsz(3)), intent(out) :: output
        integer :: nx, ny, k

        nx = this%Physical%xsz(1); ny = this%Physical%ysz(2)

        call dfftw_execute_dft_r2c(this%plan_fwd_d21x, input, this%d21x)
        this%d31x(1:nx/2+1,:,:) = this%d21x(1:nx/2+1,:,:)
        this%d31x = (two/three)*this%d31x

        call transpose_x_to_y(this%d31x,this%d31y,this%D3)
      
        do k = 1,this%D3%ysz(3)
            call dfftw_execute_dft(this%plan_fwd_d31y, this%d31y(:,:,k), this%d31y(:,:,k))  
            this%d41y(:,1:ny/2+1,k) = this%d31y(:,1:ny/2+1,k)
            this%d41y(:,ny/2+2:ny,k) = this%d31y(:,ny+2:this%nyF,k)
        end do  
        this%d41y = (two/three)*this%d41y

        do k = 1,this%spectral%ysz(3)
            call dfftw_execute_dft(this%plan_c2c_bwd_y, this%d41y(:,:,k), this%d41y(:,:,k))  
        end do 
       
        call transpose_y_to_x(this%d41y,this%d41x,this%spectral)
        call dfftw_execute_dft_c2r(this%plan_c2r_x, this%d41x, output) 
        
        output = output*this%normfactor2d 


    end subroutine
