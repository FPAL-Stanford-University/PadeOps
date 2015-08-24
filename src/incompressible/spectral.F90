module spectralMod
    use kind_parameters, only: rkind
    use decomp_2d, only: decomp_info, decomp_info_init, &
                    transpose_x_to_y, transpose_y_to_x, &
                    transpose_y_to_z, transpose_z_to_y 
    use decomp_2d_fft, only: decomp_2d_fft_init, decomp_2d_fft_finalize, decomp_2d_fft_get_size
    use exits, only: GracefulExit, message 
    use constants, only: one, zero 
    use fft_3d_stuff, only: fft_3d 
    implicit none
    private
    public :: spectral 

    logical :: useExhaustiveFFT = .false. 

    type :: spectral
        private
        real(rkind), dimension(:,:,:), allocatable, public :: k1, k2, k3, kabs_sq, k_der2, one_by_kabs_sq
        integer, dimension(3) :: fft_start, fft_end, fft_size
        real(rkind), dimension(:,:,:), allocatable, public :: Gdealias
        integer :: rPencil ! Pencil dimension for the real input
        logical :: is3dFFT = .true. ! use 3d FFTs
        logical :: isInitialized = .false.
        logical :: fixOddball = .true. 
        integer :: nx_g, ny_g, nz_g
        integer :: nx_r,ny_r,nz_r
        integer :: nx_c,ny_c,nz_c
        type(fft_3d), allocatable :: FT
        logical :: use2decompFFT = .true. 
        logical :: useConsrvD2 = .true. 
        real(rkind) :: normfact 

        contains
            procedure           :: init
            procedure           :: destroy
            procedure, private  :: alloc_r2c_out_Rank3
            procedure, private  :: alloc_r2c_out_Rank4
            generic             :: alloc_r2c_out => alloc_r2c_out_Rank4, alloc_r2c_out_Rank3
            procedure           :: fft
            procedure           :: ifft
            procedure, private  :: initializeAllPeriodic
            procedure, private  :: initializeTwoPeriodic
    
    end type

contains
    subroutine init(this,pencil, nx_g, ny_g, nz_g, dx, dy, dz, scheme, filt, dimTransform, fixOddball, use2decompFFT, useConsrvD2) 
        class(spectral),  intent(inout)         :: this
        character(len=1), intent(in)            :: pencil              ! PHYSICAL decomposition direction
        integer,          intent(in)            :: nx_g, ny_g, nz_g    ! Global data size
        real(rkind),      intent(in)            :: dx, dy, dz          ! PHYSICAL grid spacing 
        character(len=*), intent(in)            :: scheme              ! Scheme used for modified wavenumbers
        character(len=*), intent(in)            :: filt                ! Scheme used for modified wavenumbers
        integer,          intent(in), optional  :: dimTransform        ! 2 or 3 (number of periodic directions) - default is 3
        logical,                      optional  :: fixOddball          ! Fix the oddball wavenumber - default is TRUE
        logical,                      optional  :: use2decompFFT       ! Use the 3d fft procedures defined in 2decomp - default is TRUE
        logical,                      optional  :: useConsrvD2         ! Use the conservative form of 2nd derivative - default is TRUE
        


        if (this%isInitialized) then
            call GracefulExit("You are trying to reinitialize SPECTRAL derived type. This is not allowed", 111)
        end if
        
        this%nx_g = nx_g
        this%ny_g = ny_g
        this%nz_g = nz_g

        if (present(fixOddball)) this%fixOddball = fixOddball
        if (present(use2decompFFT)) this%use2decompFFT = use2decompFFT 
        if (present(useConsrvD2)) this%useConsrvD2 = useConsrvD2

        if (present(dimTransform)) then
            if (dimTransform == 2) then
                call this%initializeTwoPeriodic(pencil,nx_g,ny_g,nz_g,dx,dy,dz,scheme,filt)
                this%is3dFFT = .false.
            else if (dimTransform == 3) then
                call this%initializeAllPeriodic(pencil,nx_g,ny_g,nz_g,dx,dy,dz,scheme,filt)
            else
                call GracefulExit("Incorrect choice for DIMTRANSFORM while initializing SPECTRAL derived type. &
                                   & Available options include 2 and 3", 102)
            end if 
        else
            call this%initializeAllPeriodic(pencil,nx_g,ny_g,nz_g,dx,dy,dz,scheme,filt)
        end if 
      
        this%normfact = one/real(nx_g)/real(ny_g)/real(nz_g) 
        this%isInitialized = .true.  
 
    end subroutine 


    subroutine initializeAllPeriodic(this,pencil,nx_g,ny_g,nz_g,dx,dy,dz,scheme,filt)
        class(spectral),  intent(inout)         :: this
        character(len=1), intent(in)            :: pencil              ! PHYSICAL decomposition direction
        integer,          intent(in)            :: nx_g, ny_g, nz_g    ! Global data size
        real(rkind),      intent(in)            :: dx, dy, dz          ! PHYSICAL grid spacing 
        character(len=*), intent(in)            :: scheme              ! Scheme used for modified wavenumbers
        character(len=*), intent(in)            :: filt                ! Scheme used for modified wavenumbers

        real(rkind), dimension(nx_g) :: k1_1d 
        real(rkind), dimension(ny_g) :: k2_1d 
        real(rkind), dimension(nz_g) :: k3_1d 
        type(decomp_info), allocatable :: spectdecomp
        real(rkind), dimension(:,:,:), allocatable :: tmp1, tmp2
        integer :: i, j, k, rPencil, ierr  
       
        ! STEP 0: Figure out what the input decomposition is going to be  
        select case (pencil)
        case("x")
            rPencil = 1
        case("z")
            rPencil = 3
        case("y")
            call GracefulExit("Y - input is currently unsupported",100) 
        case default
            call GracefulExit("Incorrect pencil direction chosen during initializing SPECTRAL derived type",101)
        end select

        this%rPencil = rPencil

        call message("===========================================================")
        call message("Now Generating the SPECTRAL - Derived Type for the problem.")

        ! STEP 1: Start the 3d FFT engine
        if (this%use2decompFFT) then
            call message(2, "WARNING: Using the 2decomp 3dFFT. Performance will be terrible if &
                        & 2decomp was not compiled using FFTW")
            call message(2, "Heading to 2DECOMP for an appropriate FFT Library.")
            call decomp_2d_fft_init(rPencil, nx_g, ny_g, nz_g)   
        else
            call message ( " ***** Using the FFTW (version 3.x) engine &
                    & (found in: src/utilities/fft_3d.F90)  ***** ")
            allocate(this%FT)
            ierr = this%FT%init(nx_g,ny_g,nz_g,pencil, dx, &
                dy,dz, useExhaustiveFFT, .false., .false.)  
            if (ierr == 0) then
                call message (3, "Successfully initialized!")
            else
                call GracefulExit("Couldn't initialize 3d FFT inside SPECTRAL derived type",123)
            end if    
        end if 

        ! STEP 2: Allocate wavenumbers and temporary decomp for spectral transposes
        allocate(spectdecomp)
        select case (rPencil)
        case(1)
           call decomp_info_init(nx_g/2+1, ny_g, nz_g, spectdecomp) 
           this%fft_size(1) = spectdecomp%zsz(1)    
           this%fft_size(2) = spectdecomp%zsz(2)    
           this%fft_size(3) = spectdecomp%zsz(3)    
        case(3)
           call decomp_info_init(nx_g, ny_g, nz_g/2+1, spectdecomp) 
           this%fft_size(1) = spectdecomp%xsz(1)    
           this%fft_size(2) = spectdecomp%xsz(2)    
           this%fft_size(3) = spectdecomp%xsz(3)    
        end select 
        
        
        if (allocated(this%k1)) deallocate(this%k1)
        allocate (this%k1(this%fft_size(1),this%fft_size(2),this%fft_size(3)))     
        if (allocated(this%k2)) deallocate(this%k2)
        allocate (this%k2(this%fft_size(1),this%fft_size(2),this%fft_size(3)))     
        if (allocated(this%k3)) deallocate(this%k3)
        allocate (this%k3(this%fft_size(1),this%fft_size(2),this%fft_size(3)))     
        if (allocated(this%kabs_sq)) deallocate(this%kabs_sq)
        allocate (this%kabs_sq(this%fft_size(1),this%fft_size(2),this%fft_size(3)))     
        if (allocated(this%k_der2)) deallocate(this%k_der2)
        allocate (this%k_der2(this%fft_size(1),this%fft_size(2),this%fft_size(3)))     
        if (allocated(this%one_by_kabs_sq)) deallocate(this%one_by_kabs_sq)
        allocate(this%one_by_kabs_sq(this%fft_size(1),this%fft_size(2),this%fft_size(3)))
        if (allocated(this%Gdealias)) deallocate(this%Gdealias)
        allocate (this%Gdealias(this%fft_size(1),this%fft_size(2),this%fft_size(3)))     
        

        ! STEP 3: Generate 1d wavenumbers 
        k1_1d = GetWaveNums(nx_g,dx) 
        k2_1d = GetWaveNums(ny_g,dy) 
        k3_1d = GetWaveNums(nz_g,dz) 
       
        if (this%fixOddball) then
            k1_1d(nx_g/2+1) = zero
            k2_1d(ny_g/2+1) = zero
            k3_1d(nz_g/2+1) = zero
        end if 

       
        ! STEP 4: Create temporary array for k1 and transpose it to the appropriate dimension
        allocate(tmp1(spectdecomp%xsz(1),spectdecomp%xsz(2),spectdecomp%xsz(3)))
        do k = 1,size(tmp1,3)
            do j = 1,size(tmp1,2)
                tmp1(1:spectdecomp%xsz(1),j,k) = k1_1d(1:spectdecomp%xsz(1))
            end do 
        end do 
        
        select case (rPencil)
        case(1)
            allocate(tmp2(spectdecomp%ysz(1),spectdecomp%ysz(2),spectdecomp%ysz(3)))
            call transpose_x_to_y(tmp1,tmp2,spectdecomp)
            call transpose_y_to_z(tmp2,this%k1,spectdecomp)
            deallocate(tmp2)
        case(3)
            this%k1 = tmp1
        end select     
        deallocate(tmp1)

        ! STEP 5: Create temporary array for k2 and transpose it to the appropriate dimension
        allocate(tmp1(spectdecomp%ysz(1),spectdecomp%ysz(2),spectdecomp%ysz(3))) 
        do k = 1,size(tmp1,3)
            do i = 1,size(tmp1,1)
                tmp1(i,1:spectdecomp%ysz(2),k) = k2_1d(1:spectdecomp%ysz(2))
            end do 
        end do 

        select case (rPencil)
        case (1)
            call transpose_y_to_z(tmp1,this%k2,spectdecomp)
        case(3)
            call transpose_y_to_x(tmp1,this%k2,spectdecomp)
        end select 
        deallocate(tmp1)

        ! STEP 6: Create temporary array for k3 and transpose it to the appropriate dimension
        allocate(tmp1(spectdecomp%zsz(1),spectdecomp%zsz(2),spectdecomp%zsz(3)))
        do j = 1,size(tmp1,2)
            do i = 1,size(tmp1,1)
                tmp1(i,j,1:spectdecomp%zsz(3)) = k3_1d(1:spectdecomp%zsz(3))
            end do 
        end do 

        select case (rPencil)
        case (1)
            this%k3 = tmp1
        case (3)
            allocate(tmp2(spectdecomp%ysz(1),spectdecomp%ysz(2),spectdecomp%ysz(3)))
            call transpose_z_to_y(tmp1,tmp2,spectdecomp)
            call transpose_y_to_x(tmp2,this%k3,spectdecomp)
            deallocate(tmp2)
        end select
        deallocate(tmp1)
        deallocate (spectdecomp) 
        this%kabs_sq = this%k1**2 + this%k2**2 + this%k3**2 

        ! STEP 7: Create the dealiasing filter transfer function 
        select case (filt)
        case ("2/3rd")
            this%Gdealias = TwoThirdsRule(nx_g,ny_g,nz_g,this%kabs_sq)
        case ("cf90")
            allocate(tmp1(size(this%Gdealias,1),size(this%Gdealias,2),size(this%Gdealias,3)))
            tmp1 = GetCF90TransferFunction(this%k1,dx)
            this%Gdealias = tmp1
            tmp1 = GetCF90TransferFunction(this%k2,dy)
            this%Gdealias = this%Gdealias*tmp1
            tmp1 = GetCF90TransferFunction(this%k3,dz)
            this%Gdealias = this%Gdealias*tmp1
            deallocate(tmp1)
        case default
            call GracefulExit("The dealiasing filter specified is incorrect.",104)
        end select
       
        ! STEP 8: Correct the wavenumber to be the modified wavenumber based on the scheme
        allocate(tmp1(size(this%k1,1),size(this%k1,2),size(this%k1,3)))
        select case (trim(scheme))
        case ("cd10")
            tmp1 = GetCD10ModWaveNum(this%k1,dx)
            this%k1 = tmp1 
            tmp1 = GetCD10ModWaveNum(this%k2,dy)
            this%k2 = tmp1 
            tmp1 = GetCD10ModWaveNum(this%k3,dz)
            this%k3 = tmp1 
        case ("cd08")
            tmp1 = GetCD06ModWaveNum(this%k1,dx)
            this%k1 = tmp1 
            tmp1 = GetCD06ModWaveNum(this%k2,dy)
            this%k2 = tmp1 
            tmp1 = GetCD06ModWaveNum(this%k3,dz)
            this%k3 = tmp1 
        case ("four")
            ! Do nothing !
        end select 
        deallocate(tmp1)        

        ! STEP 9: Get the other wavenumber dependent attributes
        this%kabs_sq = this%k1**2 + this%k2**2 + this%k3**2 
        do k = 1,size(this%kabs_sq,3)
            do j = 1,size(this%kabs_sq,2)  
                do i = 1,size(this%kabs_sq,1)
                    if ((this%kabs_sq(i,j,k)) .lt. real(1d-16,rkind)) then
                        this%one_by_kabs_sq(i,j,k) = 1d-16
                    else
                        this%one_by_kabs_sq(i,j,k) = one/this%kabs_sq(i,j,k) 
                    end if 
                end do 
            end do 
        end do 
        
        if (this%useConsrvD2) then
            this%k_der2 = this%kabs_sq
        else
            call GracefulExit("CODE INCOMPLETE: Non-conservative second derivative is not available right now",123)
        end if  
        
        ! STEP 10: Determine the sizes of the fft and ifft input arrays
        this%nx_c = this%fft_size(1); this%ny_c = this%fft_size(2); this%nz_c = this%fft_size(3)
        allocate(spectdecomp)
        call decomp_info_init(nx_g, ny_g, nz_g, spectdecomp) 
        select case (rPencil)
        case (1)
            this%nx_r = spectdecomp%xsz(1)
            this%ny_r = spectdecomp%xsz(2)
            this%nz_r = spectdecomp%xsz(3)
        case (3)
            this%nx_r = spectdecomp%zsz(1)
            this%ny_r = spectdecomp%zsz(2)
            this%nz_r = spectdecomp%zsz(3)
        end select
        deallocate(spectdecomp)

        call message("SPECTRAL - Derived Type for the problem generated successfully.")
        call message("===============================================================")
        ! Finished !
    end subroutine
    
    subroutine initializeTwoPeriodic(this,pencil,nx_g,ny_g,nz_g,dx,dy,dz,scheme,filt)
        class(spectral),  intent(inout)         :: this
        character(len=1), intent(in)            :: pencil              ! PHYSICAL decomposition direction
        integer,          intent(in)            :: nx_g, ny_g, nz_g    ! Global data size
        real(rkind),      intent(in)            :: dx, dy, dz          ! PHYSICAL grid spacing 
        character(len=*), intent(in)            :: scheme              ! Scheme used for modified wavenumbers
        character(len=*), intent(in)            :: filt                ! Scheme used for modified wavenumbers
    
        call GracefulExit("CODE INCOMPLETE: Currently 2d FFT is not supported. Only fully periodic problems allowed", 103)
        ! Rubbish code that avoids compile time warnings. 
        print*, this%isInitialized
        print*, nx_g, ny_g, nz_g
        print*, dx, dy, dz
        print*, scheme
        print*, filt
        print*, pencil 
    end subroutine

    subroutine destroy(this)
        class(spectral), intent(inout) :: this
      
        if (.not. this%isInitialized) then
            call GracefulExit("You are trying to destroy a SPECTRAL derived type befire initializing it",110)
        end if 
        if (this%use2decompFFT) then
            call decomp_2d_fft_finalize 
            this%use2decompFFT = .false. 
        else
            call this%FT%destroy
            deallocate(this%FT)
        end if 
        deallocate(this%Gdealias)
        deallocate( this%k1, this%k2, this%k3, this%kabs_sq, this%k_der2, this%one_by_kabs_sq)
        this%isInitialized = .false. 
    end subroutine 

    subroutine fft(this,arr_in,arr_out)
        use decomp_2d_fft, only: decomp_2d_fft_3d
        class(spectral), intent(inout) :: this
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r) :: arr_in
        complex(rkind), dimension(this%nx_c,this%ny_c,this%nz_c) :: arr_out

        if (.not. this%use2decompFFT) then
            select case (this%rPencil) 
            case (1)
                call this%FT%fft3_x2z(arr_in,arr_out)
            case (3)
                call this%FT%fft3_z2x(arr_in,arr_out)
            end select
        else 
            call decomp_2d_fft_3d(arr_in,arr_out)
        end if 

    end subroutine

    subroutine ifft(this,arr_in,arr_out)
        use decomp_2d_fft, only: decomp_2d_fft_3d
        class(spectral), intent(inout) :: this
        complex(rkind), dimension(this%nx_c,this%ny_c,this%nz_c) :: arr_in
        real(rkind), dimension(this%nx_r,this%ny_r,this%nz_r) :: arr_out 

        if (.not. this%use2decompFFT) then
            select case (this%rPencil)
            case (1)
                call this%FT%ifft3_z2x(arr_in,arr_out)
            case (3)
                call this%FT%ifft3_x2z(arr_in,arr_out)
            end select
        else  
            call decomp_2d_fft_3d(arr_in,arr_out)
            arr_out = arr_out*this%normfact
        end if 
    end subroutine

    subroutine alloc_r2c_out_Rank3(this,arr_inout)
        class(spectral), intent(in) :: this
        complex(rkind), dimension(:,:,:), allocatable, intent(out) :: arr_inout
        
        if (.not. this%isInitialized) then
            call GracefulExit("You cannot use alloc_r2c_out before initializing the SPECTRAL derived type",112)
        end if 
        allocate(arr_inout(this%fft_size(1),this%fft_size(2),this%fft_size(3)))
    end subroutine
    
    subroutine alloc_r2c_out_Rank4(this,arr_inout,numvars)
        class(spectral), intent(in) :: this
        integer, intent(in) :: numvars
        complex(rkind), dimension(:,:,:,:), allocatable, intent(out) :: arr_inout
        
        if (.not. this%isInitialized) then
            call GracefulExit("You cannot use alloc_r2c_out before initializing the SPECTRAL derived type",112)
        end if 
        allocate(arr_inout(this%fft_size(1),this%fft_size(2),this%fft_size(3),numvars))
    end subroutine
    
    function TwoThirdsRule(Nx,Ny,Nz,kin) result(Tf)
        use constants, only: one, zero, three
        real(rkind), intent(in), dimension(:,:,:) :: kin
        integer, intent(in) :: Nx, Ny, Nz
        real(rkind), dimension(size(kin,1),size(kin,2),size(kin,3)) :: Tf

        real(rkind) :: kd1, kd2, kd3, kdealias, kdealias_sq
        
        kd1 = real(floor(real(Nx)/three))
        kd2 = real(floor(real(Ny)/three))
        kd3 = real(floor(real(Nz)/three))
        kdealias = min(kd1,kd2,kd3)
        where (kin .gt. kdealias**2) 
            Tf = zero
        elsewhere
            Tf = one
        end where

    end function


    pure elemental function GetCF90TransferFunction(kin,dx) result(T)
        use constants, only: two
        use cf90stuff, only: alpha90, beta90, a90, b90, c90, d90, e90
        real(rkind), intent(in) :: kin,dx
        real(rkind) :: k
        real(rkind) :: T
    
        k = kin*dx
        T = (a90 + two* b90*COS(k) + two*c90*COS(two*k) +two*d90*COS(3._rkind*k) + two*e90*COS(4._rkind*k) ) &
          / (1._rkind + two*alpha90*COS(k) + two*beta90*COS(two*k) )
    end function 

    pure elemental function GetCD10ModWaveNum(kin,dx) result(kp)
        use constants, only: two
        use cd10stuff, only: alpha10d1, beta10d1, a10d1, b10d1, c10d1
        real(rkind), intent(in) :: kin,dx
        real(rkind) :: k
        real(rkind) :: kp
    
        k = kin*dx
        kp = ( two*a10d1*sin(k) + two*b10d1*sin(two*k) + two*c10d1*sin(3._rkind*k) ) / (1._rkind + two*alpha10d1*cos(k) + two*beta10d1*cos(two*k))
        kp = kp/dx 

    end function
    
    pure elemental function GetCD06ModWaveNum(kin,dx) result(kp)
        use constants, only: two
        use cd06stuff, only: alpha06d1, a06d1, b06d1
        real(rkind), intent(in) :: kin,dx
        real(rkind) :: k
        real(rkind) :: kp
    
        k = kin*dx
        kp = ( two*a06d1*sin(k) + two*b06d1*sin(two*k) ) / (1._rkind + two*alpha06d1*cos(k))
        kp = kp/dx 
    
    end function
    
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
