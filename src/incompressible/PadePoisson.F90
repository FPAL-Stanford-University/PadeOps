module PadePoissonMod
    use kind_parameters, only: rkind
    use spectralMod, only: spectral, GetWaveNums 
    use constants, only: pi, zero, one, three, five, two, imi
    use exits, only: message, GracefulExit
    use decomp_2d
    use cd06staggstuff, only: cd06stagg
    use mpi
    use reductions, only: p_maxval   
 
    implicit none 

    private

    public :: Padepoisson
    
    include "fftw3.f"

    integer :: fft_planning = FFTW_PATIENT
    complex(rkind), parameter :: czero = zero + imi*zero
    
    type :: Padepoisson
        integer :: nzG, nx_in, ny_in, nz_in
        integer ::      nxE_in, nyE_in, nzE_in
        type(cd06stagg), pointer :: derZ
        type(spectral), pointer :: sp, spE
        type(decomp_info), pointer :: sp_gp, sp_gpE

        real(rkind), dimension(:,:,:), pointer :: k1_2d, k2_2d
        real(rkind), dimension(:,:,:), allocatable :: kradsq_inv

        complex(rkind), dimension(:,:,:), allocatable :: f2d, w2, f2dy
        complex(rkind), dimension(:,:,:), allocatable :: f2dext, wext

        complex(rkind), dimension(:), allocatable :: k3modcm, k3modcp

        real(rkind) :: mfact 
        integer(kind=8) :: plan_c2c_fwd_z
        integer(kind=8) :: plan_c2c_bwd_z

        contains
            procedure :: init
            procedure :: PressureProjection
            procedure :: destroy
            procedure :: DivergenceCheck
    end type

contains


    subroutine init(this, dx, dy, dz, sp, spE, derZ)
        class(padepoisson), intent(inout) :: this
        real(rkind), intent(in) :: dx, dy, dz
        class(cd06stagg), intent(in), target :: derZ
        class(spectral), intent(in), target :: sp, spE
        integer :: nxExt, nyExt, nzExt
        real(rkind), dimension(:), allocatable :: k1, k2, k3, k3mod
        complex(rkind), dimension(:), allocatable :: tfm, tfp
        real(rkind) :: kradsq
        integer     :: myxst, myyst, ii, jj, kk


        call  message("=========================================")
        call  message(0,"Initializing PADEPOISSON derived type")
        this%derZ => derZ
        this%sp => sp
        this%spE => spE
        this%sp_gp => sp%spectdecomp
        this%sp_gpE => spE%spectdecomp

        this%nzG = sp%spectdecomp%zsz(3)
        nxExt = sp%spectdecomp%zsz(1)
        nyExt = sp%spectdecomp%zsz(2)
        nzExt = 2*sp%spectdecomp%zsz(3)

        this%k1_2d => sp%k1; this%k2_2d => sp%k2
        allocate( this%f2dext(nxExt,nyExt,nzExt), this%wext(nxExt, nyExt, nzExt) )
   

        ! Prep modified wavenumbers
        allocate(k3(nzExt), k3mod(nzExt))
        allocate(tfm(nzExt), tfp(nzExt))
        allocate(k1(sp%nx_g), k2(sp%ny_g))

        k1 = GetWaveNums(sp%nx_g, dx)
        k2 = GetWaveNums(sp%ny_g, dy)
        k3 = GetWaveNums(nzExt, dz)
        call getmodCD06stagg(k3, dz, k3mod)
        tfm = exp(imi*(-dz/two)*k3); tfp = exp(imi*( dz/two)*k3)

        allocate(this%k3modcm(nzExt), this%k3modcp(nzExt))
        this%k3modcm = k3mod*tfm; this%k3modcp = k3mod*tfp

        allocate(this%kradsq_inv(nxExt, nyExt, nzExt))
        
        myxst = this%sp_gp%zst(1); myyst = this%sp_gp%zst(2)
        if ((myxst .ne. this%sp_gpE%zst(1)).or.(myyst .ne. this%sp_gpE%zst(2))) then
            call GracefulExit("Failed at initializing Padepoisson. sp_gp and sp_gpE &
                        have different x and y starts in z-decomp",423)
        end if 
        
        do kk = 1,nzExt
            do jj = 1,nyExt
                do ii = 1, nxExt
                    kradsq = k1(ii-1+myxst)**2 + k2(jj-1+myyst)**2 + k3mod(kk)**2
                    if (kradsq .le. 1d-14) then
                        this%kradsq_inv(ii,jj,kk) = zero
                    else
                        this%kradsq_inv(ii,jj,kk) = one/kradsq
                    end if  
                end do 
            end do 
        end do 
     
        deallocate(k1, k2, k3, k3mod, tfm, tfp)
        
        ! Prep for fft fwd and fft bwd
        call  message(1,"Generating Plans for FFT in PadePoisson. This could take a while ...")
        call dfftw_plan_many_dft(this%plan_c2c_fwd_z, 1, nzExt,nxExt*nyExt, this%f2dext, nzExt, &
                     nxExt*nyExt, 1, this%f2dext, nzExt,nxExt*nyExt, 1, FFTW_FORWARD, fft_planning)   

        call dfftw_plan_many_dft(this%plan_c2c_bwd_z, 1, nzExt,nxExt*nyExt, this%f2dext, nzExt, &
                     nxExt*nyExt, 1, this%f2dext, nzExt,nxExt*nyExt, 1, FFTW_BACKWARD, fft_planning)   
        
        
        this%mfact = one/real(nzExt,rkind)
        ! Prep for input array info and buffers 
        this%nx_in = this%sp_gp%ysz(1);  this%ny_in = this%sp_gp%ysz(2);  this%nz_in = this%sp_gp%ysz(3);  
        this%nxE_in = this%sp_gpE%ysz(1);  this%nyE_in = this%sp_gpE%ysz(2);  this%nzE_in = this%sp_gpE%ysz(3);  
        
        allocate(this%f2d(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)))
        allocate(this%f2dy(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)))
        allocate(this%w2(this%sp_gpE%zsz(1),this%sp_gpE%zsz(2),this%sp_gpE%zsz(3)))

        call  message(0,"PADEPOISSON derived type initialized successfully.")
        call  message("=========================================")

    end subroutine

    subroutine PressureProjection(this, uhat, vhat, what)
        class(PadePoisson), intent(inout) :: this
        complex(rkind), dimension(this%nx_in, this%ny_in, this%nz_in), intent(inout) :: uhat, vhat
        complex(rkind), dimension(this%nxE_in, this%nyE_in, this%nzE_in), intent(inout) :: what
        integer :: ii, jj, kk
        integer :: nz 
        
        nz = this%nzG

        ! Step 1: compute dudx + dvdy
        do kk = 1,size(uhat,3)
            do jj = 1,size(uhat,2)
                do ii = 1,size(uhat,1)
                    this%f2dy(ii,jj,kk) = this%k1_2d(ii,jj,kk)*uhat(ii,jj,kk)
                    this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + this%k2_2d(ii,jj,kk)*vhat(ii,jj,kk)
                end do
            end do 
        end do  
        this%f2dy = imi*this%f2dy

        ! Step 2: Transpose from y -> z
        call transpose_y_to_z(this%f2dy,this%f2d,this%sp_gp)
        call transpose_y_to_z(what, this%w2, this%sp_gpE)

        ! Step 3: Create Extensions
        do kk = 1,nz
            this%f2dext(:,:,kk) = this%f2d(:,:,nz-kk+1)
        end do
        this%f2dext(:,:,nz+1:2*nz) = this%f2d
        do kk = 2,nz
            this%wext(:,:,kk-1) = -this%w2(:,:,nz-kk+2)
        end do 
        this%wext(:,:,nz:2*nz) = this%w2

        ! Step 4: Take Fourier Transform 
        call dfftw_execute_dft(this%plan_c2c_fwd_z, this%f2dext, this%f2dext)  
        call dfftw_execute_dft(this%plan_c2c_fwd_z, this%wext  , this%wext  )  
      
        ! Step 5: Solve the Poisson System and project out w velocity
        do kk = 1,size(this%f2dext,3) 
            do jj = 1,size(this%f2dext,2) 
                do ii = 1,size(this%f2dext,1) 
                    this%f2dext(ii,jj,kk) = this%f2dext(ii,jj,kk) + imi*this%k3modcm(kk)*this%wext(ii,jj,kk)
                    this%f2dext(ii,jj,kk) = -this%f2dext(ii,jj,kk)*this%kradsq_inv(ii,jj,kk)
                end do 
            end do 
        end do 

        ! Step 6: Project out the w velocity
        do kk = 1,size(this%f2dext,3) 
            do jj = 1,size(this%f2dext,2) 
                do ii = 1,size(this%f2dext,1) 
                    this%wext(ii,jj,kk) = this%wext(ii,jj,kk) - imi*this%k3modcp(kk)*this%f2dext(ii,jj,kk)
                end do 
            end do 
        end do

        ! Step 7: Take inverse Fourier Transform
        call dfftw_execute_dft(this%plan_c2c_bwd_z, this%f2dext, this%f2dext)  
        call dfftw_execute_dft(this%plan_c2c_bwd_z, this%wext  , this%wext  )  
        this%f2dext = this%f2dext*this%mfact
        this%wext = this%wext*this%mfact

        this%f2d = this%f2dext(:,:,nz+1:2*nz)
        this%w2  = this%wext(:,:,nz:2*nz)
        this%w2(:,:,1) = czero; this%w2(:,:,nz+1) = czero; 


        ! Step 8: Transpose back from z -> y
        call transpose_z_to_y(this%w2,what,this%sp_gpE)
        call transpose_z_to_y(this%f2d,this%f2dy,this%sp_gp)

        ! Step 9: Project out u and v
        do kk = 1,size(uhat,3)
            do jj = 1,size(uhat,2)
                do ii = 1,size(uhat,1)
                    uhat(ii,jj,kk) = uhat(ii,jj,kk) - imi*this%k1_2d(ii,jj,kk)*this%f2dy(ii,jj,kk)
                end do 
            end do 
        end do 

        do kk = 1,size(vhat,3)
            do jj = 1,size(vhat,2)
                do ii = 1,size(vhat,1)
                    vhat(ii,jj,kk) = vhat(ii,jj,kk) - imi*this%k2_2d(ii,jj,kk)*this%f2dy(ii,jj,kk)
                end do 
            end do 
        end do 
       
    end subroutine 

    subroutine destroy(this)
        class(Padepoisson), intent(inout) :: this

        call dfftw_destroy_plan (this%plan_c2c_fwd_z) 
        call dfftw_destroy_plan (this%plan_c2c_bwd_z) 
        deallocate(this%kradsq_inv)
        deallocate(this%f2dext,  this%wext)
        deallocate(this%k3modcm, this%k3modcp)
        deallocate(this%f2d, this%f2dy, this%w2)
        nullify(this%derZ, this%sp_gp, this%sp_gpE, this%sp, this%spE)

    end subroutine


    subroutine DivergenceCheck(this,uhat,vhat,what,divergence, fixDiv)
        class(Padepoisson), intent(inout) :: this
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(inout) :: uhat,vhat
        complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: what

        real(rkind), dimension(this%sp_gp%xsz(1),this%sp_gp%xsz(2),this%sp_gp%xsz(3)), intent(out) :: divergence

        real(rkind) :: maxDiv
        integer :: ii, jj, kk, ierr
        logical, optional, intent(in) :: fixDiv
        logical :: fixDivergence
        
        fixDivergence = .false. 

        if (present(fixDiv)) then
            fixDivergence = fixDiv
        end if 

        ! Step 1: Compute dwdz
        call transpose_y_to_z(what, this%w2, this%sp_gpE)
        call this%derZ%ddz_E2C(this%w2,this%f2d, size(this%w2,1),size(this%w2,2))
        call transpose_z_to_y(this%f2d,this%f2dy, this%sp_gp)
        
        ! Step 2: compute and add dudx + dvdy
        do kk = 1,size(uhat,3)
            do jj = 1,size(uhat,2)
                do ii = 1,size(uhat,1)
                    this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + imi*this%k1_2d(ii,jj,kk)*uhat(ii,jj,kk)
                    this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + imi*this%k2_2d(ii,jj,kk)*vhat(ii,jj,kk)
                end do
            end do 
        end do  

        ! Step 3: Take inverse Fourier Transform
        call this%sp%ifft(this%f2dy,divergence)


        if (fixDivergence) then
            maxDiv = p_maxval(abs(divergence))

            if (maxDiv > 1d-13) then
                call this%PressureProjection(uhat,vhat,what)
                ! Step 1: Compute dwdz
                call transpose_y_to_z(what, this%w2, this%sp_gpE)
                call this%derZ%ddz_E2C(this%w2,this%f2d, size(this%w2,1),size(this%w2,2))
                call transpose_z_to_y(this%f2d,this%f2dy, this%sp_gp)
                ! Step 2: compute and add dudx + dvdy
                do kk = 1,size(uhat,3)
                    do jj = 1,size(uhat,2)
                        do ii = 1,size(uhat,1)
                            this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + imi*this%k1_2d(ii,jj,kk)*uhat(ii,jj,kk)
                            this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + imi*this%k2_2d(ii,jj,kk)*vhat(ii,jj,kk)
                        end do
                    end do 
                end do  
                ! Step 3: Take inverse Fourier Transform
                call this%sp%ifft(this%f2dy,divergence)
                maxDiv = p_maxval(abs(divergence))
                
                if (maxDiv > 1d-13) then
                    call message(0,"Maximum Divergence is found to be:", maxDiv)
                    call mpi_barrier(mpi_comm_world, ierr)
                    call GracefulExit("Divergence is not zero.",342)
                end if 

            end if 
        end if 

    end subroutine

    pure subroutine getmodCD06stagg(k,dx,kp)
        real(rkind), dimension(:), intent(in)  :: k
        real(rkind), intent(in) :: dx
        real(rkind), dimension(:), intent(out) :: kp
        real(rkind), dimension(:), allocatable :: omega
        real(rkind), parameter :: alpha = 9._rkind/62._rkind
        real(rkind), parameter :: beta = 0._rkind
        real(rkind), parameter :: a = 63._rkind/62._rkind
        real(rkind), parameter :: b = 17._rkind/62._rkind
        real(rkind), parameter :: c = 0._rkind

        allocate(omega(size(k)))
        omega = k*dx
        kp = (two*a*sin(omega/two) + (two/three)*b*sin(three*omega/two) + &
             (two/five)*c*sin(five*omega/two))/(one + two*alpha*cos(omega) +&
              two*beta*cos(two*omega))
        kp = kp/dx
        deallocate(omega)

    end subroutine

end module 
