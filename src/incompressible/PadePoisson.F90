module PadePoissonMod
    use kind_parameters, only: rkind
    use spectralMod, only: spectral, GetWaveNums, useExhaustiveFFT 
    use constants, only: pi, zero, one, three, five, two, imi
    use exits, only: message, GracefulExit, nancheck
    use decomp_2d
    use cd06staggstuff, only: cd06stagg
    use mpi
    use reductions, only: p_maxval   
    use gridtools, only: linspace 
    use timer, only: tic, toc
    use PadeDerOps, only: Pade6stagg  
    implicit none 

    private

    public :: Padepoisson
    
    include "fftw3.f"

    integer :: fft_planning 
    complex(rkind), parameter :: czero = zero + imi*zero
    
    type :: Padepoisson
        integer :: nzG, nx_in, ny_in, nz_in
        integer ::      nxE_in, nyE_in, nzE_in
        class(cd06stagg), allocatable :: derZ
        class(Pade6stagg), pointer :: derivZ 
        class(spectral), pointer :: sp, spE
        class(decomp_info), pointer :: sp_gp, sp_gpE

        real(rkind), dimension(:,:,:), pointer :: k1_2d, k2_2d
        real(rkind), dimension(:,:,:), allocatable :: kradsq_inv

        complex(rkind), dimension(:,:,:), allocatable :: f2d, w2, f2dy
        complex(rkind), dimension(:,:,:), allocatable :: f2dext, wext

        complex(rkind), dimension(:), allocatable :: k3modcm, k3modcp

        real(rkind), dimension(:,:), allocatable :: denFact, k1inZ, k2inZ, lambda
        complex(rkind), dimension(:,:), allocatable :: chat, phat, dpdzhat
        real(rkind), dimension(:,:,:), allocatable :: sinh_top, cosh_top, sinh_bot, cosh_bot
        real(rkind), dimension(:), allocatable :: zCell, zEdge
        real(rkind) :: Lz

        complex(rkind), dimension(:,:,:), allocatable :: phat_z1, phat_z2
        complex(rkind), dimension(:,:,:), allocatable :: uhatInZ, vhatInZ, whatInZ
      
        complex(rkind), dimension(:,:,:), allocatable :: dwdzHATz_Periodic, div
        real(rkind) :: mfact 
        integer(kind=8) :: plan_c2c_fwd_z
        integer(kind=8) :: plan_c2c_bwd_z

        type(decomp_info), pointer :: gpC
        logical :: computeStokesPressure = .false.
        
        logical :: PeriodicInZ = .false. 

        contains
            procedure :: init
            procedure :: PoissonSolver_HomogeneousNeumannBCz
            procedure, private :: PeriodicProjection
            procedure, private :: InitPeriodicPoissonSolver
            procedure :: PressureProjection
            procedure, private :: ProjectStokesPressure
            procedure, private :: GetStokesPressure
            procedure, private :: Periodic_getPressure
            procedure, private :: Periodic_getPressureAndUpdateRHS
            procedure :: destroy
            procedure :: DivergenceCheck
            procedure :: getPressure
            procedure :: getPressureAndUpdateRHS 
    end type

contains

    subroutine InitPeriodicPoissonSolver(this,dx,dy,dz)
         class(padepoisson), intent(inout) :: this
         real(rkind), intent(in) :: dx, dy, dz
         real(rkind), dimension(:), allocatable :: k1, k2, k3_loc_mod, k3_loc
         
         integer :: ii, jj, kk, nz, nxT, nyT, myxst, myyst
         real(rkind) :: kradsq

         allocate(this%kradsq_inv(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)))
         allocate(k3_loc(this%sp_gp%zsz(3)))
         allocate(k3_loc_mod(this%sp_gp%zsz(3)))
         allocate(k1(this%sp%nx_g), k2(this%sp%ny_g))
        
         allocate(this%dwdzHATz_Periodic(this%sp_gpE%zsz(1),this%sp_gpE%zsz(2),this%sp_gpE%zsz(3)))
         nxT = this%sp_gp%zsz(1)
         nyT = this%sp_gp%zsz(2)
         nz  = this%sp_gp%zsz(3)
         
         k1 = GetWaveNums(this%sp%nx_g, dx)
         k2 = GetWaveNums(this%sp%ny_g, dy)
         k3_loc = GetWaveNums(nz, dz)
         call this%sp%GetModifiedWavenumber_xy_ip(k1,dx)
         call this%sp%GetModifiedWavenumber_xy_ip(k2,dy)
         call this%derivZ%getModifiedWavenumbers(k3_loc,k3_loc_mod)
         
         myxst = this%sp_gp%zst(1); myyst = this%sp_gp%zst(2)
         do kk = 1,nz
             do jj = 1,nyT
                 do ii = 1,nxT
                     kradsq = k1(ii-1+myxst)**2 + k2(jj-1+myyst)**2 + k3_loc_mod(kk)**2
                     if (kradsq .le. 1d-14) then
                         this%kradsq_inv(ii,jj,kk) = zero
                     else
                         this%kradsq_inv(ii,jj,kk) = one/kradsq
                     end if 
                 end do 
             end do 
         end do 
         call  message(1,"Generating Plans for FFT in PadePoisson. This could take a while ...")
         
     
         call dfftw_plan_many_dft(this%plan_c2c_fwd_z, 1, nz,nxT*nyT, this%f2d, nz, &
                      nxT*nyT, 1, this%f2d, nz,nxT*nyT, 1, FFTW_FORWARD , fft_planning)   

         call dfftw_plan_many_dft(this%plan_c2c_bwd_z, 1, nz,nxT*nyT, this%f2d, nz, &
                      nxT*nyT, 1, this%f2d, nz,nxT*nyT, 1, FFTW_BACKWARD, fft_planning)   
         
         this%mfact = one/real(nz,rkind)
                
         allocate(this%uhatinZ(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)))
         allocate(this%vhatinZ(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)))
         
         deallocate(k1,k2,k3_loc, k3_loc_mod)

    end subroutine 

    subroutine init(this, dx, dy, dz, sp, spE, computeStokesPressure, Lz, storePressure, gpC, derivZ, PeriodicInZ, useTrueWavenumbers)
        class(padepoisson), intent(inout) :: this
        real(rkind), intent(in) :: dx, dy, dz
        !class(cd06stagg), intent(in), target :: derZ
        class(spectral), intent(in), target :: sp, spE
        integer :: nxExt, nyExt, nzExt
        real(rkind), dimension(:), allocatable :: k1, k2, k3, k3mod
        complex(rkind), dimension(:), allocatable :: tfm, tfp
        real(rkind) :: kradsq
        integer     :: myxst, myyst, ii, jj, kk
        logical, intent(in) :: computeStokesPressure, storePressure
        real(rkind), intent(in), optional :: Lz 
        integer :: i, j, k
        real(rkind), dimension(:,:), allocatable :: temp 
        type(decomp_info), intent(in), target :: gpC
        type(Pade6stagg), intent(in), target :: derivZ
        logical, intent(in) :: PeriodicInZ
        logical, intent(in), optional :: useTrueWavenumbers

        call  message("=========================================")
        call  message(0,"Initializing PADEPOISSON derived type")
        allocate(this%derZ)
        call this%derZ%init( gpC%zsz(3),  dz, isTopEven = .false., isBotEven = .false., & 
                             isTopSided = .false., isBotSided = .false.) 
        this%sp => sp
        this%spE => spE
        this%sp_gp => sp%spectdecomp
        this%sp_gpE => spE%spectdecomp
        this%derivZ => derivZ
        this%PeriodicInZ = PeriodicInZ

        this%nzG = sp%spectdecomp%zsz(3)

        this%computeStokesPressure = computeStokesPressure
        this%k1_2d => sp%k1; this%k2_2d => sp%k2
   
        if (useExhaustiveFFT) then
            fft_planning = FFTW_EXHAUSTIVE
        else
            fft_planning = FFTW_MEASURE
        end if
            
        allocate(this%f2d(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)))
        allocate(this%f2dy(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)))
        allocate(this%w2(this%sp_gpE%zsz(1),this%sp_gpE%zsz(2),this%sp_gpE%zsz(3)))
        ! Prep for input array info and buffers 
        this%nx_in = this%sp_gp%ysz(1);  this%ny_in = this%sp_gp%ysz(2);  this%nz_in = this%sp_gp%ysz(3);  
        this%nxE_in = this%sp_gpE%ysz(1);  this%nyE_in = this%sp_gpE%ysz(2);  this%nzE_in = this%sp_gpE%ysz(3);  
        this%gpC => gpC

        if (this%PeriodicInZ) then
            call this%InitPeriodicPoissonSolver(dx,dy,dz)
        else
            nxExt = sp%spectdecomp%zsz(1)
            nyExt = sp%spectdecomp%zsz(2)
            nzExt = 2*sp%spectdecomp%zsz(3)
            allocate( this%f2dext(nxExt,nyExt,nzExt), this%wext(nxExt, nyExt, nzExt) )
            ! Prep modified wavenumbers
            allocate(k3(nzExt), k3mod(nzExt))
            allocate(tfm(nzExt), tfp(nzExt))
            allocate(k1(sp%nx_g), k2(sp%ny_g))

            k1 = GetWaveNums(sp%nx_g, dx)
            k2 = GetWaveNums(sp%ny_g, dy)
            k3 = GetWaveNums(nzExt, dz)
            
            call this%sp%GetModifiedWavenumber_xy_ip(k1,dx)
            call this%sp%GetModifiedWavenumber_xy_ip(k2,dy)
            if (present(useTrueWavenumbers)) then
                if(useTrueWavenumbers) then
                    k3mod = k3
                    !k3mod(nzExt/2+1) = 0.d0
                else
                    call this%derivZ%getModifiedWavenumbers(k3,k3mod)
                end if 
            else
                call this%derivZ%getModifiedWavenumbers(k3,k3mod)
            end if 

            tfm = exp(imi*(-dz/two)*k3); tfp = exp(imi*( dz/two)*k3)

            allocate(this%k3modcm(nzExt), this%k3modcp(nzExt))
            this%k3modcm = k3mod*tfm; this%k3modcp = k3mod*tfp

            allocate(this%kradsq_inv(nxExt, nyExt, nzExt))
            
            myxst = this%sp_gp%zst(1); myyst = this%sp_gp%zst(2)
            if ((myxst .ne. this%sp_gpE%zst(1)).or.(myyst .ne. this%sp_gpE%zst(2))) then
                call GracefulExit("Failed at initializing Padepoisson. sp_gp and sp_gpE &
                            & have different x and y starts in z-decomp",423)
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
     
            
            ! Prep for fft fwd and fft bwd
            call  message(1,"Generating Plans for FFT in PadePoisson. This could take a while ...")
            call dfftw_plan_many_dft(this%plan_c2c_fwd_z, 1, nzExt,nxExt*nyExt, this%f2dext, nzExt, &
                         nxExt*nyExt, 1, this%f2dext, nzExt,nxExt*nyExt, 1, FFTW_FORWARD, fft_planning)   

            call dfftw_plan_many_dft(this%plan_c2c_bwd_z, 1, nzExt,nxExt*nyExt, this%f2dext, nzExt, &
                         nxExt*nyExt, 1, this%f2dext, nzExt,nxExt*nyExt, 1, FFTW_BACKWARD, fft_planning)   
            
            
            this%mfact = one/real(nzExt,rkind)
            

            if (this%computeStokesPressure) then
                allocate(this%lambda(this%sp_gp%zsz(1),this%sp_gp%zsz(2)))
                allocate(this%k1inZ(this%sp_gp%zsz(1),this%sp_gp%zsz(2)))
                allocate(this%k2inZ(this%sp_gp%zsz(1),this%sp_gp%zsz(2)))
                allocate(this%chat(this%sp_gp%zsz(1),this%sp_gp%zsz(2)))
                allocate(this%phat(this%sp_gp%zsz(1),this%sp_gp%zsz(2)))
                allocate(this%dpdzhat(this%sp_gp%zsz(1),this%sp_gp%zsz(2)))
                allocate(this%denFact(this%sp_gp%zsz(1),this%sp_gp%zsz(2)))
                allocate(temp(this%sp_gp%zsz(1),this%sp_gp%zsz(2)))
                
                allocate(this%sinh_top(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)+1))
                allocate(this%cosh_top(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)))
                allocate(this%sinh_bot(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)+1))
                allocate(this%cosh_bot(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)))
                this%Lz = Lz
                allocate(this%zCell(this%nzG), this%zEdge(this%nzG+1))
                this%zEdge = linspace(0.d0, Lz, this%nzG+1)
                this%zCell = 0.5d0*(this%zEdge(1:this%nzG) + this%zEdge(2:this%nzG+1))

                do j = 1,this%sp_gp%zsz(2)!this%sp_gp%zst(2),this%sp_gp%zen(2)
                    do i = 1,this%sp_gp%zsz(1)!this%sp_gp%zst(1),this%sp_gp%zen(1)
                        this%k1inZ(i,j) = k1(this%sp_gp%zst(1) - 1 + i)
                        this%k2inZ(i,j) = k2(this%sp_gp%zst(2) - 1 + j)
                    end do 
                end do
                this%lambda = sqrt(this%k1inZ**2 + this%k2inZ**2)
                temp = this%lambda*Lz
                where (temp < 500.d0) 
                    this%denFact = 1.d0/(this%lambda*sinh(this%lambda*Lz) + 1.d-13)
                elsewhere
                    this%denFact = zero
                end where
                where(this%denFact<1d-16)
                    this%denFact = 0.d0
                end where
                if (nrank == 0) this%denFact(1,1) = 0.d0
                

                do k = 1,this%sp_gp%zsz(3)
                    temp = this%lambda*(Lz - this%zCell(k))
                    where (temp < 32.d0)
                        this%cosh_bot(:,:,k) =  cosh(temp)
                    elsewhere
                        this%cosh_bot(:,:,k) = 4.d13
                    end where
                    temp = this%lambda*this%zCell(k)
                    where (temp < 32.d0) 
                        this%cosh_top(:,:,k) =  cosh(temp)
                    elsewhere
                        this%cosh_top(:,:,k) = 1.d13
                    end where
                end do 

                do k = 1,this%sp_gp%zsz(3)+1
                    temp = this%lambda*(Lz - this%zEdge(k))
                    where(temp<32.d0)
                        this%sinh_bot(:,:,k) = -this%lambda*sinh(temp)
                    elsewhere
                        this%sinh_bot(:,:,k) = -4.d13
                    end where
                    temp = this%lambda*(this%zEdge(k))
                    where(temp<32.d0) 
                        this%sinh_top(:,:,k) =  this%lambda*sinh(temp)
                    elsewhere
                        this%sinh_top(:,:,k) = 4.d13
                    end where
                end do 
                allocate(this%uhatinZ(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)))
                allocate(this%vhatinZ(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)))
                deallocate(temp)
                call message(1,"STOKES PRESSURE calculation enabled with the CD06 Poisson solver")
            end if 
            if (storePressure) then
                allocate(this%phat_z1(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)))
                allocate(this%phat_z2(this%sp_gp%zsz(1),this%sp_gp%zsz(2),this%sp_gp%zsz(3)))
                this%phat_z1 = cmplx(0.0_rkind, 0.0_rkind, rkind)
                this%phat_z2 = cmplx(0.0_rkind, 0.0_rkind, rkind)
            end if

            deallocate(k1, k2, k3, k3mod, tfm, tfp)
        end if 
        call  message(0,"PADEPOISSON derived type initialized successfully.")
        call  message("=========================================")

    end subroutine

    subroutine ProjectStokesPressure(this, uhat, vhat, what)
        class(PadePoisson), intent(inout) :: this
        complex(rkind), dimension(this%nx_in, this%ny_in, this%nz_in), intent(in) :: uhat, vhat
        complex(rkind), dimension(this%nxE_in, this%nyE_in, this%nzE_in), intent(in) :: what
        integer :: i,j,k

        ! First transpose to z 
        call transpose_y_to_z(uhat, this%uhatInZ, this%sp_gp)
        call transpose_y_to_z(vhat, this%vhatInZ, this%sp_gp)
        call transpose_y_to_z(what, this%w2, this%sp_gpE)
       
        ! Bottom BC
        do j = 1,this%sp_gp%zsz(2)
            !$omp simd
            do i = 1,this%sp_gp%zsz(1)
               this%chat(i,j) = -this%w2(i,j,1)*this%denFact(i,j)
            end do 
        end do 

        do k = 1,this%nzG
            do j = 1,this%sp_gp%zsz(2)
               !$omp simd
               do i = 1,this%sp_gp%zsz(1)
                  this%phat(i,j) = imi*this%chat(i,j) * this%cosh_bot(i,j,k)
                  this%uhatInZ(i,j,k) = this%uhatInZ(i,j,k) - this%k1inZ(i,j)*this%phat(i,j) 
                  this%vhatInZ(i,j,k) = this%vhatInZ(i,j,k) - this%k2inZ(i,j)*this%phat(i,j)
               end do 
            end do 
        end do 

        this%w2(:,:,1) = czero
        do k = 2,this%nzG+1
            do j = 1,this%sp_gp%zsz(2)
               !$omp simd
               do i = 1,this%sp_gp%zsz(1)
                  this%dpdzhat(i,j) = this%chat(i,j)*this%sinh_bot(i,j,k)
                  this%w2(i,j,k) = this%w2(i,j,k) - this%dpdzhat(i,j)
               end do 
            end do 
        end do 


        ! Top BC
        do j = 1,this%sp_gp%zsz(2)
            !$omp simd
            do i = 1,this%sp_gp%zsz(1)
               this%chat(i,j) = this%w2(i,j,this%nzG+1)*this%denFact(i,j)
            end do 
        end do 
        
        do k = 1,this%nzG
            do j = 1,this%sp_gp%zsz(2)
               !$omp simd
               do i = 1,this%sp_gp%zsz(1)
                  this%phat(i,j) = imi*this%chat(i,j) * this%cosh_top(i,j,k)
                  this%uhatInZ(i,j,k) = this%uhatInZ(i,j,k) - this%k1inZ(i,j)*this%phat(i,j) 
                  this%vhatInZ(i,j,k) = this%vhatInZ(i,j,k) - this%k2inZ(i,j)*this%phat(i,j)
                  this%dpdzhat(i,j) = this%chat(i,j)*this%sinh_top(i,j,k)
                  this%w2(i,j,k) = this%w2(i,j,k) - this%dpdzhat(i,j)
               end do 
            end do 
        end do 
        this%w2(:,:,this%nzG+1) = czero

    end subroutine

    subroutine PeriodicProjection(this, uhat, vhat, what)
        class(PadePoisson), intent(inout) :: this
        complex(rkind), dimension(this%nx_in, this%ny_in, this%nz_in), intent(inout) :: uhat, vhat
        complex(rkind), dimension(this%nxE_in, this%nyE_in, this%nzE_in), intent(inout) :: what
        integer :: ii, jj, kk

        do kk = 1,size(uhat,3)
            do jj = 1,size(uhat,2)
                !$omp simd
                do ii = 1,size(uhat,1)
                    this%f2dy(ii,jj,kk) = this%k1_2d(ii,jj,kk)*uhat(ii,jj,kk)
                    this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + this%k2_2d(ii,jj,kk)*vhat(ii,jj,kk)
                    this%f2dy(ii,jj,kk) = dcmplx(-dimag(this%f2dy(ii,jj,kk)),dreal(this%f2dy(ii,jj,kk)))
                end do
            end do 
        end do  

        call transpose_y_to_z(this%f2dy,this%uhatInZ,this%sp_gp)
        call transpose_y_to_z(what, this%w2, this%sp_gpE)

        call this%derivZ%ddz_E2C(this%w2, this%f2d, 0, 0)  ! Periodic dwdz
        
        this%f2d = this%f2d + this%uhatInZ
        call dfftw_execute_dft(this%plan_c2c_fwd_z, this%f2d    , this%f2d    )  

        this%f2d = -this%kradsq_inv*this%f2d
        call dfftw_execute_dft(this%plan_c2c_bwd_z, this%f2d    , this%f2d    )  
        this%f2d = this%f2d*this%mfact 

        call this%derivZ%ddz_C2E(this%f2d, this%dwdzHATz_Periodic, 0, 0)
        this%w2 = this%w2 - this%dwdzHATz_Periodic
        call transpose_z_to_y(this%w2,what, this%sp_gpE)

        call transpose_z_to_y(this%f2d, this%f2dy, this%sp_gp)


        do kk = 1,size(uhat,3)
           do jj = 1,size(uhat,2)
              !$omp simd
              do ii = 1,size(vhat,1)
                  uhat(ii,jj,kk) = uhat(ii,jj,kk) - imi*this%k1_2d(ii,jj,kk)*this%f2dy(ii,jj,kk)
                  vhat(ii,jj,kk) = vhat(ii,jj,kk) - imi*this%k2_2d(ii,jj,kk)*this%f2dy(ii,jj,kk)
              end do 
           end do 
        end do 

    end subroutine

    subroutine PressureProjection(this, uhat, vhat, what)
        class(PadePoisson), intent(inout) :: this
        complex(rkind), dimension(this%nx_in, this%ny_in, this%nz_in), intent(inout) :: uhat, vhat
        complex(rkind), dimension(this%nxE_in, this%nyE_in, this%nzE_in), intent(inout) :: what
        integer :: ii, jj, kk
     
        if (this%PeriodicInZ) then
            call this%PeriodicProjection(uhat, vhat, what)
        else
            ! Step 0: Project out the stokes pressure which will fix the wall BCs
            ! for velocity
            if (this%computeStokesPressure) then
                call this%ProjectStokesPressure(uhat, vhat, what)
                do kk = 1,this%nzG
                   do jj = 1,this%sp_gp%zsz(2)
                      !$omp simd
                      do ii = 1,this%sp_gp%zsz(1)
                         this%f2d(ii,jj,kk) = this%k1inZ(ii,jj)*this%uhatInZ(ii,jj,kk)
                         this%f2d(ii,jj,kk) = this%f2d(ii,jj,kk) + this%k2inZ(ii,jj)*this%vhatInZ(ii,jj,kk)
                         this%f2d(ii,jj,kk) = dcmplx(-dimag(this%f2d(ii,jj,kk)),dreal(this%f2d(ii,jj,kk)))
                      end do 
                   end do
                end do 
                !this%f2d = imi*this%f2d
            else
                ! Step 1: compute dudx + dvdy
                do kk = 1,size(uhat,3)
                    do jj = 1,size(uhat,2)
                        !$omp simd
                        do ii = 1,size(uhat,1)
                            this%f2dy(ii,jj,kk) = this%k1_2d(ii,jj,kk)*uhat(ii,jj,kk)
                            this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + this%k2_2d(ii,jj,kk)*vhat(ii,jj,kk)
                            this%f2dy(ii,jj,kk) = dcmplx(-dimag(this%f2dy(ii,jj,kk)),dreal(this%f2dy(ii,jj,kk)))
                        end do
                    end do 
                end do  
                !this%f2dy = imi*this%f2dy

                ! Step 2: Transpose from y -> z
                call transpose_y_to_z(this%f2dy,this%f2d,this%sp_gp)
                call transpose_y_to_z(what, this%w2, this%sp_gpE)
            end if

            ! Step 3: Create Extensions
            do kk = 1,this%nzG
                do jj = 1,size(this%f2dext,2)
                   !$omp simd
                   do ii = 1,size(this%f2dext,1)
                      this%f2dext(ii,jj,kk) = this%f2d(ii,jj,this%nzG-kk+1)
                   end do 
                end do 
            end do
            do kk = 1,size(this%f2d,3)
                do jj = 1,size(this%f2d,2)
                   !$omp simd
                   do ii = 1,size(this%f2d,1)
                      !this%f2dext(ii,jj,this%nzG+1:2*this%nzG) = this%f2d
                      this%f2dext(ii,jj,this%nzG+kk) = this%f2d(ii,jj,kk)
                   end do
                end do
             end do

            do kk = 2,this%nzG
                do jj = 1,size(this%wext,2)
                   !$omp simd
                   do ii = 1,size(this%wext,1)
                      this%wext(ii,jj,kk-1) = -this%w2(ii,jj,this%nzG-kk+2)
                   end do 
                end do
            end do 
            do kk = 1,size(this%w2,3)
                do jj = 1,size(this%w2,2)
                   !$omp simd
                   do ii = 1,size(this%w2,1)
                      !this%wext(:,:,this%nzG:2*this%nzG) = this%w2
                      this%wext(ii,jj,this%nzG+kk-1) = this%w2(ii,jj,kk)
                   end do
                end do
             end do

            
            ! Step 4: Take Fourier Transform 
            call dfftw_execute_dft(this%plan_c2c_fwd_z, this%f2dext, this%f2dext)  
            call dfftw_execute_dft(this%plan_c2c_fwd_z, this%wext  , this%wext  )  
      

            ! Step 5: Solve the Poisson System 
            do kk = 1,size(this%f2dext,3) 
                do jj = 1,size(this%f2dext,2) 
                    !$omp simd
                    do ii = 1,size(this%f2dext,1) 
                        this%f2dext(ii,jj,kk) = this%f2dext(ii,jj,kk) + imi*this%k3modcm(kk)*this%wext(ii,jj,kk)
                        this%f2dext(ii,jj,kk) = -this%f2dext(ii,jj,kk)*this%kradsq_inv(ii,jj,kk)
                    end do 
                end do 
            end do 

            ! Step 6: Project out the w velocity
            do kk = 1,size(this%f2dext,3) 
                do jj = 1,size(this%f2dext,2) 
                    !$omp simd
                    do ii = 1,size(this%f2dext,1) 
                        this%wext(ii,jj,kk) = this%wext(ii,jj,kk) - imi*this%k3modcp(kk)*this%f2dext(ii,jj,kk)
                    end do 
                end do 
            end do

            ! Step 7: Take inverse Fourier Transform
            call dfftw_execute_dft(this%plan_c2c_bwd_z, this%f2dext, this%f2dext)  
            call dfftw_execute_dft(this%plan_c2c_bwd_z, this%wext  , this%wext  )  
            do kk = 1,size(this%f2dext,3)
                do jj = 1,size(this%f2dext,2)
                   !$omp simd
                   do ii = 1,size(this%f2dext,1)
                      this%f2dext(ii,jj,kk) = this%mfact*this%f2dext(ii,jj,kk)
                   end do
                end do
            end do

            do kk = 1,size(this%wext,3)
                do jj = 1,size(this%wext,2)
                   !$omp simd
                   do ii = 1,size(this%wext,1)
                      this%wext(ii,jj,kk) = this%wext(ii,jj,kk)*this%mfact
                   end do
                end do
            end do

            
            do kk = 1,size(this%f2d,3)
                do jj = 1,size(this%f2d,2)
                   !$omp simd
                   do ii = 1,size(this%f2d,1)
                      this%f2d(ii,jj,kk) = this%f2dext(ii,jj,this%nzG+kk)
                   end do
                end do
            end do
            
            do kk = 1,size(this%w2,3)
                do jj = 1,size(this%w2,2)
                   !$omp simd
                   do ii = 1,size(this%w2,1)
                      this%w2(ii,jj,kk)  = this%wext(ii,jj,this%nzG+kk-1)
                   end do
                end do
            end do

            this%w2(:,:,1) = czero; 
            this%w2(:,:,this%nzG+1) = czero; 


            ! Step 8: Transpose back from z -> y
            call transpose_z_to_y(this%w2,what,this%sp_gpE)
            
            if (this%computeStokesPressure) then
                !this%f2d = imi*this%f2d
                do kk = 1,this%nzG
                   do jj = 1,size(this%uhatinZ,2)
                      !$omp simd
                      do ii = 1,size(this%uhatinZ,1)
                         this%f2d(ii,jj,kk) = dcmplx(-dimag(this%f2d(ii,jj,kk)),dreal(this%f2d(ii,jj,kk)))
                         this%uhatinZ(ii,jj,kk) = this%uhatinZ(ii,jj,kk) - this%f2d(ii,jj,kk)*this%k1inZ(ii,jj)
                         this%vhatinZ(ii,jj,kk) = this%vhatinZ(ii,jj,kk) - this%f2d(ii,jj,kk)*this%k2inZ(ii,jj)
                      end do
                   end do
                end do 
                call transpose_z_to_y(this%uhatinZ,uhat,this%sp_gp) 
                call transpose_z_to_y(this%vhatinZ,vhat,this%sp_gp) 
            else
                call transpose_z_to_y(this%f2d,this%f2dy,this%sp_gp)
                ! Step 9: Project out u and v
                do kk = 1,size(uhat,3)
                    do jj = 1,size(uhat,2)
                        !$omp simd
                        do ii = 1,size(uhat,1)
                            uhat(ii,jj,kk) = uhat(ii,jj,kk) - imi*this%k1_2d(ii,jj,kk)*this%f2dy(ii,jj,kk)
                        end do 
                    end do 
                end do 

                do kk = 1,size(vhat,3)
                    do jj = 1,size(vhat,2)
                        !$omp simd
                        do ii = 1,size(vhat,1)
                            vhat(ii,jj,kk) = vhat(ii,jj,kk) - imi*this%k2_2d(ii,jj,kk)*this%f2dy(ii,jj,kk)
                        end do 
                    end do 
                end do 
            end if
        end if
    end subroutine 

    subroutine destroy(this)
        class(Padepoisson), intent(inout) :: this

        call dfftw_destroy_plan (this%plan_c2c_fwd_z) 
        call dfftw_destroy_plan (this%plan_c2c_bwd_z) 
        deallocate(this%kradsq_inv)
        deallocate(this%f2dext,  this%wext)
        deallocate(this%k3modcm, this%k3modcp)
        deallocate(this%f2d, this%f2dy, this%w2)
        if (allocated(this%phat_z1)) deallocate(this%phat_z1)
        if (allocated(this%phat_z2)) deallocate(this%phat_z2)
        nullify( this%sp_gp, this%sp_gpE, this%sp, this%spE)
        deallocate(this%derZ)
    end subroutine

    subroutine GetStokesPressure(this,uhat,vhat,what)
        !! NOTE :: uhat is really urhshat, vhat is really vrhshat, ...
        class(PadePoisson), intent(inout) :: this
        complex(rkind), dimension(this%nx_in, this%ny_in, this%nz_in), intent(in) :: uhat, vhat
        complex(rkind), dimension(this%nxE_in, this%nyE_in, this%nzE_in), intent(in) :: what
        integer :: i,j,k

        ! First transpose to z 
        call transpose_y_to_z(uhat, this%uhatInZ, this%sp_gp)
        call transpose_y_to_z(vhat, this%vhatInZ, this%sp_gp)
        call transpose_y_to_z(what, this%w2, this%sp_gpE)
   
        ! Bottom BC
        do j = 1,size(this%chat,2)
           !$omp simd
           do i = 1,size(this%chat,1)
               this%chat(i,j) = -this%w2(i,j,1)*this%denFact(i,j)
           end do
        end do

        do k = 1,this%nzG
           do j = 1,size(this%chat,2)
              !$omp simd
              do i = 1,size(this%chat,1)
                  this%phat_z1(i,j,k) = imi*this%chat(i,j) * this%cosh_bot(i,j,k)
                  this%uhatInZ(i,j,k) = this%uhatInZ(i,j,k) - this%k1inZ(i,j)*this%phat_z1(i,j,k)
                  this%vhatInZ(i,j,k) = this%vhatInZ(i,j,k) - this%k2inZ(i,j)*this%phat_z1(i,j,k)
               end do
            end do
        end do 

        this%w2(:,:,1) = czero
        do k = 2,this%nzG+1
            do j = 1,size(this%chat,2)
               !$omp simd
               do i = 1,size(this%chat,1)
                  this%dpdzhat(i,j) = this%chat(i,j)*this%sinh_bot(i,j,k)
                  this%w2(i,j,k) = this%w2(i,j,k) - this%dpdzhat(i,j)
               end do
            end do
        end do 


        ! Top BC
        do j = 1,size(this%chat,2)
            !$omp simd
            do i = 1,size(this%chat,1)
               this%chat(i,j) = this%w2(i,j,this%nzG+1)*this%denFact(i,j)
            end do
        end do

        do k = 1,this%nzG
            do j = 1,size(this%chat,2)
               !$omp simd
               do i = 1,size(this%chat,1)
                  this%phat_z2(i,j,k) = imi*this%chat(i,j) * this%cosh_top(i,j,k)
                  this%uhatInZ(i,j,k) = this%uhatInZ(i,j,k) - this%k1inZ(i,j)*this%phat_z2(i,j,k) 
                  this%vhatInZ(i,j,k) = this%vhatInZ(i,j,k) - this%k2inZ(i,j)*this%phat_z2(i,j,k)
               end do
            end do
        end do 

        do k = 1,this%nzG
            do j = 1,size(this%chat,2)
               !$omp simd
               do i = 1,size(this%chat,1)
                  this%dpdzhat(i,j) = this%chat(i,j)*this%sinh_top(i,j,k)
                  this%w2(i,j,k) = this%w2(i,j,k) - this%dpdzhat(i,j)
               end do
            end do
        end do 
        this%w2(:,:,this%nzG+1) = czero
        
    end subroutine

    subroutine Periodic_getPressure(this, uhat, vhat, what, pressure)
        class(PadePoisson), intent(inout) :: this
        complex(rkind), dimension(this%nx_in, this%ny_in, this%nz_in), intent(in) :: uhat, vhat
        complex(rkind), dimension(this%nxE_in, this%nyE_in, this%nzE_in), intent(in) :: what
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(inout) :: pressure

        integer :: ii, jj, kk

        do kk = 1,size(uhat,3)
            do jj = 1,size(uhat,2)
                !$omp simd
                do ii = 1,size(uhat,1)
                    this%f2dy(ii,jj,kk) = this%k1_2d(ii,jj,kk)*uhat(ii,jj,kk)
                    this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + this%k2_2d(ii,jj,kk)*vhat(ii,jj,kk)
                    this%f2dy(ii,jj,kk) = dcmplx(-dimag(this%f2dy(ii,jj,kk)),dreal(this%f2dy(ii,jj,kk)))
                end do
            end do 
        end do  

        call transpose_y_to_z(this%f2dy,this%uhatInZ,this%sp_gp)
        call transpose_y_to_z(what, this%w2, this%sp_gpE)
        
        call this%derivZ%ddz_E2C(this%w2, this%f2d, 0, 0)  ! Periodic dwdz
        
        this%f2d = this%f2d + this%uhatInZ
        call dfftw_execute_dft(this%plan_c2c_fwd_z, this%f2d    , this%f2d    )  

        this%f2d = -this%kradsq_inv*this%f2d
        call dfftw_execute_dft(this%plan_c2c_bwd_z, this%f2d    , this%f2d    )  
        this%f2d = this%f2d*this%mfact 

        call transpose_z_to_y(this%f2d, this%f2dy, this%sp_gp)
        call this%sp%ifft(this%f2dy, pressure)

    end subroutine 

    subroutine getPressure(this, uhat, vhat, what, pressure)
        !! NOTE :: uhat is really urhshat, vhat is really vrhshat, ...
        class(PadePoisson), intent(inout) :: this
        complex(rkind), dimension(this%nx_in, this%ny_in, this%nz_in), intent(in) :: uhat, vhat
        complex(rkind), dimension(this%nxE_in, this%nyE_in, this%nzE_in), intent(in) :: what
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(inout) :: pressure
        integer :: ii, jj, kk

        if (this%PeriodicInZ) then
            call this%Periodic_getPressure(uhat, vhat, what, pressure)
        else
            ! Step 0: compute the stokes pressure which will fix the wall BCs
            ! for velocity rhs
            if (this%computeStokesPressure) then
                call this%GetStokesPressure(uhat, vhat, what)
                do kk = 1,this%nzG
                   do jj = 1,size(this%f2d,2)
                      !$omp simd
                      do ii = 1,size(this%f2d,1)
                         this%f2d(ii,jj,kk) = this%k1inZ(ii,jj)*this%uhatInZ(ii,jj,kk)
                         this%f2d(ii,jj,kk) = this%f2d(ii,jj,kk) + this%k2inZ(ii,jj)*this%vhatInZ(ii,jj,kk)
                         this%f2d(ii,jj,kk) = dcmplx(-dimag(this%f2d(ii,jj,kk)),dreal(this%f2d(ii,jj,kk)))
                      end do
                   end do 
                end do 
            else
                ! Step 1: compute dudx + dvdy
                do kk = 1,size(uhat,3)
                    do jj = 1,size(uhat,2)
                        !$omp simd
                        do ii = 1,size(uhat,1)
                            this%f2dy(ii,jj,kk) = this%k1_2d(ii,jj,kk)*uhat(ii,jj,kk)
                            this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + this%k2_2d(ii,jj,kk)*vhat(ii,jj,kk)
                            this%f2dy(ii,jj,kk) = dcmplx(-dimag(this%f2dy(ii,jj,kk)),dreal(this%f2dy(ii,jj,kk)))
                        end do
                    end do 
                end do  

                ! Step 2: Transpose from y -> z
                call transpose_y_to_z(this%f2dy,this%f2d,this%sp_gp)
                call transpose_y_to_z(what, this%w2, this%sp_gpE)
            end if

            ! Step 3: Create Extensions
            !do kk = 1,this%nzG
            !    this%f2dext(:,:,kk) = this%f2d(:,:,this%nzG-kk+1)
            !end do
            !this%f2dext(:,:,this%nzG+1:2*this%nzG) = this%f2d
            !do kk = 2,this%nzG
            !    this%wext(:,:,kk-1) = -this%w2(:,:,this%nzG-kk+2)
            !end do 
            !this%wext(:,:,this%nzG:2*this%nzG) = this%w2

            ! Step 3: Create Extensions
            do kk = 1,this%nzG
                do jj = 1,size(this%f2dext,2)
                   !$omp simd
                   do ii = 1,size(this%f2dext,1)
                      this%f2dext(ii,jj,kk) = this%f2d(ii,jj,this%nzG-kk+1)
                   end do 
                end do 
            end do
            do kk = 1,size(this%f2d,3)
                do jj = 1,size(this%f2d,2)
                   !$omp simd
                   do ii = 1,size(this%f2d,1)
                      !this%f2dext(ii,jj,this%nzG+1:2*this%nzG) = this%f2d
                      this%f2dext(ii,jj,this%nzG+kk) = this%f2d(ii,jj,kk)
                   end do
                end do
             end do

            do kk = 2,this%nzG
                do jj = 1,size(this%wext,2)
                   !$omp simd
                   do ii = 1,size(this%wext,1)
                      this%wext(ii,jj,kk-1) = -this%w2(ii,jj,this%nzG-kk+2)
                   end do 
                end do
            end do 
            do kk = 1,size(this%w2,3)
                do jj = 1,size(this%w2,2)
                   !$omp simd
                   do ii = 1,size(this%w2,1)
                      !this%wext(:,:,this%nzG:2*this%nzG) = this%w2
                      this%wext(ii,jj,this%nzG+kk-1) = this%w2(ii,jj,kk)
                   end do
                end do
             end do

            ! Step 4: Take Fourier Transform 
            call dfftw_execute_dft(this%plan_c2c_fwd_z, this%f2dext, this%f2dext)  
            call dfftw_execute_dft(this%plan_c2c_fwd_z, this%wext  , this%wext  )  
      
            ! Step 5: Solve the Poisson System and project out w velocity
            do kk = 1,size(this%f2dext,3) 
                do jj = 1,size(this%f2dext,2) 
                    !$omp simd
                    do ii = 1,size(this%f2dext,1) 
                        this%f2dext(ii,jj,kk) = this%f2dext(ii,jj,kk) + imi*this%k3modcm(kk)*this%wext(ii,jj,kk)
                        this%f2dext(ii,jj,kk) = -this%f2dext(ii,jj,kk)*this%kradsq_inv(ii,jj,kk)
                    end do 
                end do 
            end do 

            ! Step 6: Take inverse Fourier Transform
            call dfftw_execute_dft(this%plan_c2c_bwd_z, this%f2dext, this%f2dext)  
            !this%f2dext = this%f2dext*this%mfact
            do kk = 1,size(this%f2dext,3)
                do jj = 1,size(this%f2dext,2)
                   !$omp simd
                   do ii = 1,size(this%f2dext,1)
                      this%f2dext(ii,jj,kk) = this%mfact*this%f2dext(ii,jj,kk)
                   end do
                end do
            end do
            
            
            !this%f2d = this%f2dext(:,:,this%nzG+1:2*this%nzG)
            do kk = 1,size(this%f2d,3)
                do jj = 1,size(this%f2d,2)
                   !$omp simd
                   do ii = 1,size(this%f2d,1)
                      this%f2d(ii,jj,kk) = this%f2dext(ii,jj,this%nzG+kk)
                   end do
                end do
            end do

            if (this%computeStokesPressure) then
               do kk = 1,size(this%f2d,3)
                   do jj = 1,size(this%f2d,2)
                      !$omp simd
                      do ii = 1,size(this%f2d,1)
                         this%f2d(ii,jj,kk) = this%f2d(ii,jj,kk) + this%phat_z1(ii,jj,kk) + this%phat_z2(ii,jj,kk)
                      end do
                   end do
               end do
            end if 

            ! Step 8: Transpose back from z -> y
            call transpose_z_to_y(this%f2d,this%f2dy,this%sp_gp)

            ! Step 9: Get pressure by ifft
            call this%sp%ifft(this%f2dy,pressure)

         end if 
    end subroutine 

    subroutine PoissonSolver_HomogeneousNeumannBCz(this, frhs, pressure) 
        class(PadePoisson), intent(inout) :: this
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(out) :: pressure
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in) :: frhs
        integer :: ii, jj, kk
      
        ! Step 1: Take fft_xy for rhs and transpose to z
        call this%sp%fft(frhs,this%f2dy)
        call transpose_y_to_z(this%f2dy,this%f2d,this%sp_gp)

         ! Step 2: Create Extensions
         do kk = 1,this%nzG
             do jj = 1,size(this%f2dext,2)
                !$omp simd
                do ii = 1,size(this%f2dext,1)
                   this%f2dext(ii,jj,kk) = this%f2d(ii,jj,this%nzG-kk+1)
                end do 
             end do 
         end do
         do kk = 1,size(this%f2d,3)
            do jj = 1,size(this%f2d,2)
               !$omp simd
               do ii = 1,size(this%f2d,1)
                  !this%f2dext(ii,jj,this%nzG+1:2*this%nzG) = this%f2d
                  this%f2dext(ii,jj,this%nzG+kk) = this%f2d(ii,jj,kk)
               end do
            end do
         end do

         ! Step 3: Take Fourier Transform 
         call dfftw_execute_dft(this%plan_c2c_fwd_z, this%f2dext, this%f2dext)  

         ! Step 4: Solve the Poisson equation 
         do kk = 1,size(this%f2dext,3) 
             do jj = 1,size(this%f2dext,2) 
                 !$omp simd
                 do ii = 1,size(this%f2dext,1) 
                     this%f2dext(ii,jj,kk) = -this%f2dext(ii,jj,kk)*this%kradsq_inv(ii,jj,kk)
                 end do 
             end do 
         end do 

         ! Step 5: Take inverse Fourier Transform
         call dfftw_execute_dft(this%plan_c2c_bwd_z, this%f2dext, this%f2dext)  
         !this%f2dext = this%f2dext*this%mfact
         do kk = 1,size(this%f2dext,3)
             do jj = 1,size(this%f2dext,2)
                !$omp simd
                do ii = 1,size(this%f2dext,1)
                   this%f2dext(ii,jj,kk) = this%mfact*this%f2dext(ii,jj,kk)
                end do
             end do
         end do

         ! Step 6: Extract the top half
         do kk = 1,size(this%f2d,3)
             do jj = 1,size(this%f2d,2)
                !$omp simd
                do ii = 1,size(this%f2d,1)
                   this%f2d(ii,jj,kk) = this%f2dext(ii,jj,this%nzG+kk)
                end do
             end do
         end do

         ! Step 7: transpose to y and take inverse fourier transform
         call transpose_z_to_y(this%f2d,this%f2dy,this%sp_gp)
         call this%sp%ifft(this%f2dy,pressure)
    end subroutine 
    
    subroutine Periodic_getPressureAndUpdateRHS(this, uhat, vhat, what, pressure)
        class(PadePoisson), intent(inout) :: this
        complex(rkind), dimension(this%nx_in, this%ny_in, this%nz_in), intent(inout) :: uhat, vhat
        complex(rkind), dimension(this%nxE_in, this%nyE_in, this%nzE_in), intent(inout) :: what
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(inout) :: pressure

        integer :: ii, jj, kk

        do kk = 1,size(uhat,3)
            do jj = 1,size(uhat,2)
                !$omp simd
                do ii = 1,size(uhat,1)
                    this%f2dy(ii,jj,kk) = this%k1_2d(ii,jj,kk)*uhat(ii,jj,kk)
                    this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + this%k2_2d(ii,jj,kk)*vhat(ii,jj,kk)
                    this%f2dy(ii,jj,kk) = dcmplx(-dimag(this%f2dy(ii,jj,kk)),dreal(this%f2dy(ii,jj,kk)))
                end do
            end do 
        end do  

        call transpose_y_to_z(this%f2dy,this%uhatInZ,this%sp_gp)
        call transpose_y_to_z(what, this%w2, this%sp_gpE)
        
        call this%derivZ%ddz_E2C(this%w2, this%f2d, 0, 0)  ! Periodic dwdz
        
        this%f2d = this%f2d + this%uhatInZ
        call dfftw_execute_dft(this%plan_c2c_fwd_z, this%f2d    , this%f2d    )  

        this%f2d = -this%kradsq_inv*this%f2d
        call dfftw_execute_dft(this%plan_c2c_bwd_z, this%f2d    , this%f2d    )  
        this%f2d = this%f2d*this%mfact 

        call this%derivZ%ddz_C2E(this%f2d, this%dwdzHATz_Periodic, 0, 0)
        this%w2 = this%w2 - this%dwdzHATz_Periodic
        call transpose_z_to_y(this%w2,what, this%sp_gpE)
        call transpose_z_to_y(this%f2d, this%f2dy, this%sp_gp)


        do kk = 1,size(uhat,3)
           do jj = 1,size(uhat,2)
              !$omp simd
              do ii = 1,size(uhat,1)
                  uhat(ii,jj,kk) = uhat(ii,jj,kk) - imi*this%k1_2d(ii,jj,kk)*this%f2dy(ii,jj,kk)
                  vhat(ii,jj,kk) = vhat(ii,jj,kk) - imi*this%k2_2d(ii,jj,kk)*this%f2dy(ii,jj,kk)
              end do 
           end do 
        end do 

        call this%sp%ifft(this%f2dy, pressure)

    end subroutine 


    subroutine getPressureAndUpdateRHS(this, uhat, vhat, what, pressure)
        !! NOTE :: uhat is really urhshat, vhat is really vrhshat, ...
        class(PadePoisson), intent(inout) :: this
        complex(rkind), dimension(this%nx_in, this%ny_in, this%nz_in), intent(inout) :: uhat, vhat
        complex(rkind), dimension(this%nxE_in, this%nyE_in, this%nzE_in), intent(inout) :: what
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(inout) :: pressure
        integer :: ii, jj, kk
        complex(rkind) :: ctemp
        
        if (this%PeriodicInZ) then
            call this%Periodic_getPressureAndUpdateRHS(uhat, vhat, what, pressure)
        else

            ! Step 0: Project out the stokes pressure which will fix the wall BCs
            ! for velocity
            if (this%computeStokesPressure) then
                call this%ProjectStokesPressure(uhat, vhat, what)
                do kk = 1,this%nzG
                   do jj = 1,this%sp_gp%zsz(2)
                      !$omp simd
                      do ii = 1,this%sp_gp%zsz(1)
                         this%f2d(ii,jj,kk) = this%k1inZ(ii,jj)*this%uhatInZ(ii,jj,kk)
                         this%f2d(ii,jj,kk) = this%f2d(ii,jj,kk) + this%k2inZ(ii,jj)*this%vhatInZ(ii,jj,kk)
                         this%f2d(ii,jj,kk) = dcmplx(-dimag(this%f2d(ii,jj,kk)),dreal(this%f2d(ii,jj,kk)))
                      end do 
                   end do
                end do 
                !this%f2d = imi*this%f2d
            else
                ! Step 1: compute dudx + dvdy
                do kk = 1,size(uhat,3)
                    do jj = 1,size(uhat,2)
                        !$omp simd
                        do ii = 1,size(uhat,1)
                            this%f2dy(ii,jj,kk) = this%k1_2d(ii,jj,kk)*uhat(ii,jj,kk)
                            this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + this%k2_2d(ii,jj,kk)*vhat(ii,jj,kk)
                            this%f2dy(ii,jj,kk) = dcmplx(-dimag(this%f2dy(ii,jj,kk)),dreal(this%f2dy(ii,jj,kk)))
                        end do
                    end do 
                end do  
                !this%f2dy = imi*this%f2dy

                ! Step 2: Transpose from y -> z
                call transpose_y_to_z(this%f2dy,this%f2d,this%sp_gp)
                call transpose_y_to_z(what, this%w2, this%sp_gpE)
            end if

            ! Step 3: Create Extensions
            do kk = 1,this%nzG
                do jj = 1,size(this%f2dext,2)
                   !$omp simd
                   do ii = 1,size(this%f2dext,1)
                      this%f2dext(ii,jj,kk) = this%f2d(ii,jj,this%nzG-kk+1)
                   end do 
                end do 
            end do
            do kk = 1,size(this%f2d,3)
                do jj = 1,size(this%f2d,2)
                   !$omp simd
                   do ii = 1,size(this%f2d,1)
                      !this%f2dext(ii,jj,this%nzG+1:2*this%nzG) = this%f2d
                      this%f2dext(ii,jj,this%nzG+kk) = this%f2d(ii,jj,kk)
                   end do
                end do
             end do

            do kk = 2,this%nzG
                do jj = 1,size(this%wext,2)
                   !$omp simd
                   do ii = 1,size(this%wext,1)
                      this%wext(ii,jj,kk-1) = -this%w2(ii,jj,this%nzG-kk+2)
                   end do 
                end do
            end do 
            do kk = 1,size(this%w2,3)
                do jj = 1,size(this%w2,2)
                   !$omp simd
                   do ii = 1,size(this%w2,1)
                      !this%wext(:,:,this%nzG:2*this%nzG) = this%w2
                      this%wext(ii,jj,this%nzG+kk-1) = this%w2(ii,jj,kk)
                   end do
                end do
             end do

            
            ! Step 4: Take Fourier Transform 
            call dfftw_execute_dft(this%plan_c2c_fwd_z, this%f2dext, this%f2dext)  
            call dfftw_execute_dft(this%plan_c2c_fwd_z, this%wext  , this%wext  )  
      

            ! Step 5: Solve the Poisson System and project out w velocity
            do kk = 1,size(this%f2dext,3) 
                do jj = 1,size(this%f2dext,2) 
                    !$omp simd
                    do ii = 1,size(this%f2dext,1) 
                        this%f2dext(ii,jj,kk) = this%f2dext(ii,jj,kk) + imi*this%k3modcm(kk)*this%wext(ii,jj,kk)
                        this%f2dext(ii,jj,kk) = -this%f2dext(ii,jj,kk)*this%kradsq_inv(ii,jj,kk)
                    end do 
                end do 
            end do 

            ! Step 6: add -dpdz to wrhs 
            do kk = 1,size(this%f2dext,3) 
                do jj = 1,size(this%f2dext,2) 
                    !$omp simd
                    do ii = 1,size(this%f2dext,1) 
                        this%wext(ii,jj,kk) = this%wext(ii,jj,kk) - imi*this%k3modcp(kk)*this%f2dext(ii,jj,kk)
                    end do 
                end do 
            end do
            
            ! Step 7: Take inverse Fourier Transform
            call dfftw_execute_dft(this%plan_c2c_bwd_z, this%f2dext, this%f2dext)  
            call dfftw_execute_dft(this%plan_c2c_bwd_z, this%wext  , this%wext  )  
            do kk = 1,size(this%f2dext,3)
                do jj = 1,size(this%f2dext,2)
                   !$omp simd
                   do ii = 1,size(this%f2dext,1)
                      this%f2dext(ii,jj,kk) = this%mfact*this%f2dext(ii,jj,kk)
                   end do
                end do
            end do

            do kk = 1,size(this%wext,3)
                do jj = 1,size(this%wext,2)
                   !$omp simd
                   do ii = 1,size(this%wext,1)
                      this%wext(ii,jj,kk) = this%wext(ii,jj,kk)*this%mfact
                   end do
                end do
            end do

            
            do kk = 1,size(this%f2d,3)
                do jj = 1,size(this%f2d,2)
                   !$omp simd
                   do ii = 1,size(this%f2d,1)
                      this%f2d(ii,jj,kk) = this%f2dext(ii,jj,this%nzG+kk)
                   end do
                end do
            end do
            
            do kk = 1,size(this%w2,3)
                do jj = 1,size(this%w2,2)
                   !$omp simd
                   do ii = 1,size(this%w2,1)
                      this%w2(ii,jj,kk)  = this%wext(ii,jj,this%nzG+kk-1)
                   end do
                end do
            end do

            this%w2(:,:,1) = czero; 
            this%w2(:,:,this%nzG+1) = czero; 


            ! Step 8: Transpose back from z -> y
            call transpose_z_to_y(this%w2,what,this%sp_gpE)
            
            ! Step 9: Update the RHS for u and v
            if (this%computeStokesPressure) then
                !this%f2d = imi*this%f2d
                do kk = 1,this%nzG
                   do jj = 1,size(this%uhatinZ,2)
                      !$omp simd
                      do ii = 1,size(this%uhatinZ,1)
                         !this%f2d(ii,jj,kk) = dcmplx(-dimag(this%f2d(ii,jj,kk)),dreal(this%f2d(ii,jj,kk)))
                         ctemp  = dcmplx(-dimag(this%f2d(ii,jj,kk)),dreal(this%f2d(ii,jj,kk)))
                         this%uhatinZ(ii,jj,kk) = this%uhatinZ(ii,jj,kk) - ctemp*this%k1inZ(ii,jj)
                         this%vhatinZ(ii,jj,kk) = this%vhatinZ(ii,jj,kk) - ctemp*this%k2inZ(ii,jj)
                      end do
                   end do
                end do 
                call transpose_z_to_y(this%uhatinZ,uhat,this%sp_gp) 
                call transpose_z_to_y(this%vhatinZ,vhat,this%sp_gp) 
            else
                call transpose_z_to_y(this%f2d,this%f2dy,this%sp_gp)
                do kk = 1,size(uhat,3)
                    do jj = 1,size(uhat,2)
                        !$omp simd
                        do ii = 1,size(uhat,1)
                            uhat(ii,jj,kk) = uhat(ii,jj,kk) - imi*this%k1_2d(ii,jj,kk)*this%f2dy(ii,jj,kk)
                        end do 
                    end do 
                end do 

                do kk = 1,size(vhat,3)
                    do jj = 1,size(vhat,2)
                        !$omp simd
                        do ii = 1,size(vhat,1)
                            vhat(ii,jj,kk) = vhat(ii,jj,kk) - imi*this%k2_2d(ii,jj,kk)*this%f2dy(ii,jj,kk)
                        end do 
                    end do 
                end do 
            end if 

            ! Step 10: Get pressure by ifft
            if (this%computeStokesPressure) then
               do kk = 1,size(this%f2d,3)
                   do jj = 1,size(this%f2d,2)
                      !$omp simd
                      do ii = 1,size(this%f2d,1)
                         this%f2d(ii,jj,kk) = this%f2d(ii,jj,kk) + this%phat_z1(ii,jj,kk) + this%phat_z2(ii,jj,kk)
                      end do
                   end do
                end do
                call transpose_z_to_y(this%f2d,this%f2dy,this%sp_gp)
            end if 
            !print *, "ph1 : ", nancheck(real(this%phat_z1)), nancheck(aimag(this%phat_z1))
            !print *, "ph2 : ", nancheck(real(this%phat_z2)), nancheck(aimag(this%phat_z2))
            !print *, "f2d : ", nancheck(real(this%f2d)), nancheck(aimag(this%f2d))
            !print *, "f2dy: ", nancheck(real(this%f2dy)), nancheck(aimag(this%f2dy))
            !call message(1, " ++++Calling ifft")
            !print *, 'size of pressure: ', size(pressure)
            !call message(1, " ++++Calling ifft")
            call this%sp%ifft(this%f2dy,pressure)

         end if 

    end subroutine 

    subroutine DivergenceCheck(this,uhat,vhat,what,divergence, fixDiv)
        class(Padepoisson), intent(inout) :: this
        complex(rkind), dimension(this%sp_gp%ysz(1),this%sp_gp%ysz(2),this%sp_gp%ysz(3)), intent(inout) :: uhat,vhat
        complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: what

        real(rkind), dimension(this%sp_gp%xsz(1),this%sp_gp%xsz(2),this%sp_gp%xsz(3)), intent(out) :: divergence
        !real(rkind), dimension(:,:,:), intent(out) :: divergence

        real(rkind) :: maxDiv, myMaxDiv
        integer :: ii, jj, kk, ierr
        logical, optional, intent(in) :: fixDiv
        logical :: fixDivergence
        
        fixDivergence = .false. 

        if (present(fixDiv)) then
            fixDivergence = fixDiv
        end if 

        ! Step 1: Compute dwdz
        call transpose_y_to_z(what, this%w2, this%sp_gpE)
        !call this%derZ%ddz_E2C(this%w2,this%f2d, size(this%w2,1),size(this%w2,2))
        call this%derivZ%ddz_E2C(this%w2,this%f2d,-1,-1)
        call transpose_z_to_y(this%f2d,this%f2dy, this%sp_gp)
        
        ! Step 2: compute and add dudx + dvdy
        do kk = 1,size(uhat,3)
            do jj = 1,size(uhat,2)
                !$omp simd
                do ii = 1,size(uhat,1)
                    this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + imi*this%k1_2d(ii,jj,kk)*uhat(ii,jj,kk) &
                                        + imi*this%k2_2d(ii,jj,kk)*vhat(ii,jj,kk)
                end do
            end do 
        end do  

        ! Step 3: Take inverse Fourier Transform
        call this%sp%ifft(this%f2dy,divergence)
        !divergence = abs(divergence)

        if (fixDivergence) then
            myMaxDiv = maxval(divergence)
            maxDiv = p_maxval(myMaxDiv)

            if (maxDiv > 1d-13) then
                call this%PressureProjection(uhat,vhat,what)
                ! Step 1: Compute dwdz
                call transpose_y_to_z(what, this%w2, this%sp_gpE)
                !call this%derZ%ddz_E2C(this%w2,this%f2d, size(this%w2,1),size(this%w2,2))
                call this%derivZ%ddz_E2C(this%w2,this%f2d,-1,-1)
                call transpose_z_to_y(this%f2d,this%f2dy, this%sp_gp)
                ! Step 2: compute and add dudx + dvdy
                do kk = 1,size(uhat,3)
                    do jj = 1,size(uhat,2)
                        !$omp simd
                        do ii = 1,size(uhat,1)
                            !this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + imi*this%k1_2d(ii,jj,kk)*uhat(ii,jj,kk)
                            !this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + imi*this%k2_2d(ii,jj,kk)*vhat(ii,jj,kk)
                            this%f2dy(ii,jj,kk) = this%f2dy(ii,jj,kk) + imi*this%k1_2d(ii,jj,kk)*uhat(ii,jj,kk) &
                                        + imi*this%k2_2d(ii,jj,kk)*vhat(ii,jj,kk)
                        end do
                    end do 
                end do  
                ! Step 3: Take inverse Fourier Transform
                call this%sp%ifft(this%f2dy,divergence)
                !divergence = abs(divergence)
                myMaxDiv = maxval(divergence)
                maxDiv = p_maxval(myMaxDiv)
                
                if (maxDiv > 1d-10) then
                    call this%PressureProjection(uhat,vhat,what)
                    call message(0,"Maximum Divergence is found to be:", maxDiv)
                    call mpi_barrier(mpi_comm_world, ierr)
                    !call GracefulExit("Divergence is not zero.",342)
                end if 

            end if 
        end if 

    end subroutine

    !pure subroutine getmodCD06stagg(k,dx,kp)
    !    real(rkind), dimension(:), intent(in)  :: k
    !    real(rkind), intent(in) :: dx
    !    real(rkind), dimension(:), intent(out) :: kp
    !    real(rkind), dimension(:), allocatable :: omega
    !    real(rkind), parameter :: alpha = 9._rkind/62._rkind
    !    real(rkind), parameter :: beta = 0._rkind
    !    real(rkind), parameter :: a = 63._rkind/62._rkind
    !    real(rkind), parameter :: b = 17._rkind/62._rkind
    !    real(rkind), parameter :: c = 0._rkind

    !    allocate(omega(size(k)))
    !    omega = k*dx
    !    kp = (two*a*sin(omega/two) + (two/three)*b*sin(three*omega/two) + &
    !         (two/five)*c*sin(five*omega/two))/(one + two*alpha*cos(omega) +&
    !          two*beta*cos(two*omega))
    !    kp = kp/dx
    !    deallocate(omega)

    !end subroutine

end module 
