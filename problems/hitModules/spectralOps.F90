module spectralOps
    ! NOTE: all input and output fields for this module are in spectral space. Z
    ! decomposition of the input is assumed.  
    
    use kind_parameters, only: rkind, mpirkind
    use variables, only: fieldsSpec, maxDivergence, buff_cmplx, buff_real, FT
    use decomp_2d
    use fft_3d_stuff, only: fft_3d
    use constants,  only: imi, zero
    implicit none

contains
    subroutine getRHS
        use variables, only: nu, rhs, fieldsPhys,fieldsSpec, Buff_real, Buff_cmplx

        ! Step 1: Compute the 6 Advection terms and add to RHS
        ! a) uu
        Buff_real = fieldsPhys(:,:,:,1)*fieldsPhys(:,:,:,1)
        call FT%fft3_x2z(Buff_real,Buff_cmplx)
        rhs(:,:,:,1) = FT%k1*Buff_cmplx

        ! b) uv
        Buff_real = fieldsPhys(:,:,:,1)*fieldsPhys(:,:,:,2)
        call FT%fft3_x2z(Buff_real,Buff_cmplx)
        rhs(:,:,:,1) = rhs(:,:,:,1) + FT%k2*Buff_cmplx(:,:,:)
        rhs(:,:,:,2) = FT%k1*Buff_cmplx
        
        ! c) uw 
        Buff_real = fieldsPhys(:,:,:,1)*fieldsPhys(:,:,:,3)
        call FT%fft3_x2z(Buff_real,Buff_cmplx)
        rhs(:,:,:,1) = rhs(:,:,:,1) + FT%k3*Buff_cmplx(:,:,:)
        rhs(:,:,:,3) = FT%k1*Buff_cmplx

        ! d) vv
        Buff_real = fieldsPhys(:,:,:,2)*fieldsPhys(:,:,:,2)
        call FT%fft3_x2z(Buff_real,Buff_cmplx)
        rhs(:,:,:,2) = rhs(:,:,:,2) + FT%k2*Buff_cmplx(:,:,:)
      
        ! e) vw
        Buff_real = fieldsPhys(:,:,:,2)*fieldsPhys(:,:,:,3)
        call FT%fft3_x2z(Buff_real,Buff_cmplx)
        rhs(:,:,:,2) = rhs(:,:,:,2) + FT%k3*Buff_cmplx(:,:,:)
        rhs(:,:,:,3) = rhs(:,:,:,3) + FT%k2*Buff_cmplx(:,:,:)

        ! f) ww
        Buff_real = fieldsPhys(:,:,:,3)*fieldsPhys(:,:,:,3)
        call FT%fft3_x2z(Buff_real,Buff_cmplx)
        rhs(:,:,:,3) = rhs(:,:,:,3) + FT%k3*Buff_cmplx(:,:,:)

        ! Multiply by i
        rhs = -imi*rhs

        ! Step 2: Add the diffusion term 
        rhs(:,:,:,1) = rhs(:,:,:,1) - nu*FT%kabs_sq*fieldsSpec(:,:,:,1)
        rhs(:,:,:,2) = rhs(:,:,:,2) - nu*FT%kabs_sq*fieldsSpec(:,:,:,2)
        rhs(:,:,:,3) = rhs(:,:,:,3) - nu*FT%kabs_sq*fieldsSpec(:,:,:,3)

        ! Done 
    end subroutine 

    subroutine pressureProjection
        use variables, only: Buff_cmplx, FT, fieldsSpec, rhs, oneByksq
        
        Buff_cmplx = FT%k1*fieldsSpec(:,:,:,1)
        Buff_cmplx = Buff_cmplx + FT%k2*fieldsSpec(:,:,:,2)
        Buff_cmplx = Buff_cmplx + FT%k3*fieldsSpec(:,:,:,3)
       
        Buff_cmplx = oneBykSq*Buff_cmplx 

        fieldsSpec(:,:,:,1) = fieldsSpec(:,:,:,1) - FT%k1*Buff_cmplx
        fieldsSpec(:,:,:,2) = fieldsSpec(:,:,:,2) - FT%k2*Buff_cmplx
        fieldsSpec(:,:,:,3) = fieldsSpec(:,:,:,3) - FT%k3*Buff_cmplx
    end subroutine    

    subroutine divergence 
        use mpi 
        real(rkind) :: mymax
        integer :: ierr 

        ! Create divergence in spectral space
        buff_cmplx = FT%k1*fieldsSpec(:,:,:,1) 
        buff_cmplx = buff_cmplx + FT%k2*fieldsSpec(:,:,:,2) 
        buff_cmplx = buff_cmplx + FT%k3*fieldsSpec(:,:,:,3) 
        buff_cmplx = imi*buff_cmplx

        ! Compute inverse FT
        call FT%ifft3_z2x(buff_cmplx,buff_real)

        ! Compute the local max value of divergence        
        mymax = maxval(abs(buff_real)) 
        call MPI_Reduce(mymax,maxDivergence, 1, mpirkind, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    end subroutine 

end module 
