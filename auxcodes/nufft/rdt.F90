module rdtMod
    use kind_parameters, only: rkind
    use exits, only: GracefulExit, warning 
    use constants, only: imi, four, six, two, zero, one, pi
    use nufftMod, only: nufft
    implicit none

    private
    public :: rdt

    integer, parameter :: ntestRNG = 2000
    real(rkind), parameter :: kminRNG = 1.d-1, kmaxRNG = 1.d3
    real(rkind), parameter :: C_vk   = 1.4528_rkind   ! KE spectrum (Kolmogorov constant)
    real(rkind), parameter :: C_vkt = 0.3962_rkind ! PE spectrum (Corrsin constant?)
    real(rkind) :: C_mann = 9.4_rkind ! Mann constant for shearing
    real(rkind), parameter :: c_nut = 0.2d0/5.d0
    type :: rdt
        ! All attributes are public (for performance)

        integer :: nwaves           ! Number of isotropic shells
        
        real(rkind), dimension(:,:), allocatable :: wavenums
        complex(rkind), dimension(:,:), allocatable :: fieldhat

        real(rkind), dimension(:), pointer :: k1, k2, k3
        complex(rkind), dimension(:), pointer :: uhat, vhat, what, thetahat 
        
        real(rkind), dimension(:), allocatable :: buf
        real(rkind), dimension(:), allocatable :: buf1, buf2
        complex(rkind), dimension(:), allocatable :: sweepFact

        complex(rkind), dimension(:), allocatable :: ctmp1
        complex(rkind), dimension(:,:,:), allocatable :: ctmp2

        real(rkind), dimension(:), allocatable :: omega
        logical :: isStratified = .false. 
        logical :: isInitialized = .false. 

        logical :: allowFastCalculations = .false. 
        type(nufft), allocatable :: FT  

        real(rkind) :: dx, dy, dz
        real(rkind) :: factx, facty, factz 
        integer :: Nx, Ny, Nz
        real(rkind) :: dispx, dispy, dispz


        !! RNG terms
        logical :: createRNGterms = .false. 
        real(rkind), dimension(:), allocatable :: ktest, Etest, nuttest,lifetimetest, dnu, bnu, cnu, dtime, btime, ctime 
        real(rkind), dimension(:), allocatable :: lifetimes, nut, localtimer

        contains
            procedure, private :: getIsotropicModes
            procedure, private :: init_basic
            procedure, private :: init_copy
            generic :: init => init_copy, init_basic
            procedure :: destroy
            procedure          :: doSlowSum
            procedure, private :: doFastSum
            procedure, private :: doFastSumPressure
            generic :: getFields => doSlowSum, doFastSum, doFastSumPressure
            procedure, private :: getduidxj_fast
            procedure, private :: getduidxj_slow
            generic :: getduidxj => getduidxj_fast, getduidxj_slow
            procedure :: extractRange
            procedure :: sweep
            procedure, private :: initializeRNG
            procedure :: setRNG
            procedure :: setOmega
    end type

contains

    subroutine init_copy(this, wavenums, fieldhat, Nx, Ny, Nz, xLims, yLims, zLims, eps_nufft, allocExtraBufs, createRNGterms)

        class(rdt), intent(inout), target :: this
        real(rkind), dimension(:,:), intent(in) :: wavenums
        complex(rkind), dimension(:,:), intent(in) :: fieldhat
        integer, intent(in) :: Nx, Ny, Nz 
        real(rkind), dimension(2), intent(in) :: xLims, yLims, zLims
        real(rkind), intent(in) :: eps_nufft
        logical, intent(in), optional :: allocExtraBufs, createRNGterms
        integer :: nfields

        if (this%isInitialized) then
            call warning("WARNING: you are reinitializing an already initialized &
            & RDT derived type. Will destroy the existing type and reinitialize")
            call this%destroy()
        end if  
     
        this%nwaves = size(wavenums,1)
        nfields = size(fieldhat,2)

        if (nfields == 4) this%isStratified = .true.
        
        allocate(this%wavenums(this%nwaves,3))
        allocate(this%fieldhat(this%nwaves,nfields))
        allocate(this%buf(this%nwaves))

        this%Nx = Nx 
        this%Ny = Ny
        this%Nz = Nz
        
        this%wavenums = wavenums
        this%fieldhat = fieldhat
        this%k1 => this%wavenums(:,1) 
        this%k2 => this%wavenums(:,2) 
        this%k3 => this%wavenums(:,3) 
        this%uhat => this%fieldhat(:,1) 
        this%vhat => this%fieldhat(:,2) 
        this%what => this%fieldhat(:,3) 

        if (present(createRNGterms)) this%createRNGterms = createRNGterms

        if (this%isStratified) then
            this%thetahat => this%fieldhat(:,4)
        end if 


        if (.not. allocated(this%FT)) allocate(this%FT)
        call this%FT%init(eps_nufft,Nx,Ny,Nz)

        if (Nx .le. 1) then
            this%dx = 1.d-13
        else
            this%dx = (xLims(2) - xLims(1))/real(Nx - 1,rkind)   
        end if

        if (Ny .le. 1) then
            this%dy = 1.d-13
        else
            this%dy = (yLims(2) - yLims(1))/real(Ny - 1,rkind) 
        end if 

        if (Nz .le. 1) then
            this%dz = 1.d-13
        else
            this%dz = (zLims(2) - zLims(1))/real(Nz - 1,rkind)  
        end if 

        this%factx = (xLims(1) + this%dx*(real(Nx)/2)) 
        this%facty = (yLims(1) + this%dy*(real(Ny)/2))
        this%factz = (zLims(1) + this%dz*(real(Nz)/2))

        this%dispx = zero
        this%dispy = zero
        this%dispz = zero
        
        allocate(this%sweepFact(this%nwaves))
        this%sweepFact = exp(-imi*(this%k1*this%dispx + this%k2*this%dispy + this%k3*this%dispz))

        allocate(this%ctmp1(this%nwaves))
        allocate(this%ctmp2(this%nx,this%ny,this%nz))

        this%allowfastcalculations = .true. 
        this%isInitialized = .true.  

        if (present(allocExtraBufs)) then
            if (allocExtraBufs) then
                allocate(this%buf1(this%nwaves))
                allocate(this%buf2(this%nwaves))
            end if 
        end if

        if (this%createRNGterms) then
            call this%initializeRNG()
        end if 
 
    end subroutine

    subroutine initializeRNG(this)
        use gridtools, only: linspace, mytrapz
        use interpolation, only: spline
        class(rdt), intent(inout), target :: this
        
        integer :: idx 

        if (allocated(this%ktest)) then
            deallocate(this%ktest, this%Etest, this%nuttest, & 
                     this%lifetimetest, this%dnu, this%bnu, &
                     this%cnu, this%dtime, this%btime, & 
                     this%ctime)
        end if 
        allocate(this%ktest(ntestRNG), this%Etest(ntestRNG), this%nuttest(ntestRNG), & 
                 this%lifetimetest(ntestRNG), this%dnu(ntestRNG), this%bnu(ntestRNG), &
                 this%cnu(ntestRNG), this%dtime(ntestRNG), this%btime(ntestRNG), & 
                 this%ctime(ntestRNG))
        this%ktest = linspace(kminRNG,kmaxRNG,ntestRNG)
        this%Etest = 1.4528d0*(this%ktest**4)/((1.d0 + this%ktest**2)**(17.d0/6.d0))
        this%nuttest = 0.d0
        this%Etest = (this%ktest**(-2.d0))*this%Etest
        do idx = 1,ntestRNG-1
            this%nuttest(idx) = sqrt(c_nut*mytrapz(this%ktest(idx:ntestRNG),this%Etest(idx:ntestRNG)));
        end do
        this%nuttest(ntestRNG) = this%nuttest(ntestRNG-1) 
        this%lifetimetest = 1.d0/((this%ktest**2)*this%nuttest)
        call spline (this%ktest, this%nuttest, this%bnu, this%cnu, this%dnu, ntestRNG)        
        call spline (this%ktest, this%lifetimetest, this%btime, this%ctime, this%dtime, ntestRNG)        
        

    end subroutine 

    subroutine setOmega(this, omega, nut)
        class(rdt), intent(inout), target :: this
        real(rkind), dimension(this%nwaves), intent(in) :: omega, nut
        
        if(allocated(this%omega)) deallocate(this%omega)
        if(allocated(this%nut)) deallocate(this%nut)
        allocate(this%omega(this%nwaves))
        allocate(this%nut(this%nwaves))
        this%omega = omega
        this%nut = nut
    end subroutine
    
    subroutine setRNG(this)
        use interpolation, only: ispline
        class(rdt), intent(inout), target :: this
        integer :: idx
        real(rkind) :: kabs

        if (allocated(this%lifetimes)) deallocate(this%lifetimes)
        allocate(this%lifetimes(this%nwaves))
        
        if (allocated(this%nut)) deallocate(this%nut)
        allocate(this%nut(this%nwaves))

        if (allocated(this%localtimer)) deallocate(this%localtimer)
        allocate(this%localtimer(this%nwaves))
        this%localtimer = 0.d0        


        do idx = 1,this%nwaves
            kabs = sqrt(this%k1(idx)**2 + this%k2(idx)**2 + this%k3(idx)**2)
            if (kabs < kminRNG) then
                this%lifetimes(idx) = this%lifetimetest(1)
                this%nut(idx) = this%nuttest(1) 
            elseif (kabs > kmaxRNG) then
                this%lifetimes(idx) = this%lifetimetest(ntestRNG)
                this%nut(idx) = this%nuttest(ntestRNG) 
            else
                this%lifetimes(idx) = ispline(kabs, this%ktest, this%lifetimetest, this%btime, this%ctime, this%dtime, ntestRNG)
                this%nut(idx)       = ispline(kabs, this%ktest, this%nuttest, this%bnu, this%cnu, this%dnu, ntestRNG)
            end if 
        end do 

    end subroutine


    subroutine init_basic(this, nk, ntheta, kmin, kmax, isStratified, seed1, seed2, seed3, Nx, Ny, Nz, xLims, yLims, zLims, eps_nufft, createRNGterms)
        class(rdt), intent(inout), target :: this
        integer, intent(in) :: nk, ntheta
        real(rkind), intent(in) :: kmin, kmax
        integer :: nwaves
        integer, intent(in), optional :: seed1, seed2, seed3
        logical, intent(in), optional :: isStratified, createRNGterms
        integer, intent(in), optional :: Nx, Ny, Nz 
        real(rkind), dimension(2), intent(in), optional :: xLims, yLims, zLims
        real(rkind), intent(in), optional :: eps_nufft

        if (this%isInitialized) then
            call warning("WARNING: you are reinitializing an already initialized &
            & RDT derived type. Will destroy the existing type and reinitialize")
            call this%destroy()
        end if  
     

        nwaves = nk*ntheta 
        this%nwaves = nwaves 
        allocate(this%wavenums(nwaves,3))
        allocate(this%buf(nwaves))
        if (present(isStratified)) then
            this%isStratified = isStratified
        end if 

        if (this%isStratified) then
            allocate(this%fieldhat(nwaves,4))
        else
            allocate(this%fieldhat(nwaves,3))
        end if 
      
        if (present(createRNGterms)) this%createRNGterms = createRNGterms
        
        this%k1 => this%wavenums(:,1) 
        this%k2 => this%wavenums(:,2) 
        this%k3 => this%wavenums(:,3) 
        this%uhat => this%fieldhat(:,1) 
        this%vhat => this%fieldhat(:,2) 
        this%what => this%fieldhat(:,3) 
        if (this%isStratified) then
            this%thetahat => this%fieldhat(:,4)
        end if 

        

        if (present(seed1) .and. present(seed2) .and. present(seed3)) then
            call this%getIsotropicModes(nk,ntheta,kmin,kmax,seed1,seed2,seed3)
        else
            call this%getIsotropicModes(nk,ntheta,kmin,kmax)
        end if 

        if ( (present(Nx)) .and. (present(Ny)) .and. (present(Nz)) .and. &
                (present(xLims)) .and. present(yLims) .and. present(zLims) .and. &
                    present(eps_nufft) ) then

            this%Nx = Nx 
            this%Ny = Ny
            this%Nz = Nz
            
            if (.not. allocated(this%FT)) allocate(this%FT)
            call this%FT%init(eps_nufft,Nx,Ny,Nz)

            if (Nx .le. 1) then
                this%dx = 1.d-13
            else
                this%dx = (xLims(2) - xLims(1))/real(Nx - 1,rkind)   
            end if

            if (Ny .le. 1) then
                this%dy = 1.d-13
            else
                this%dy = (yLims(2) - yLims(1))/real(Ny - 1,rkind) 
            end if 

            if (Nz .le. 1) then
                this%dz = 1.d-13
            else
                this%dz = (zLims(2) - zLims(1))/real(Nz - 1,rkind)  
            end if 

            this%factx = (xLims(1) + this%dx*(real(Nx)/2)) 
            this%facty = (yLims(1) + this%dy*(real(Ny)/2))
            this%factz = (zLims(1) + this%dz*(real(Nz)/2))

            this%allowfastcalculations = .true. 
        end if 
        
        this%dispx = zero
        this%dispy = zero
        this%dispz = zero
        allocate(this%sweepFact(this%nwaves))
        this%sweepFact = exp(-imi*(this%k1*this%dispx + this%k2*this%dispy + this%k3*this%dispz))
        
        allocate(this%ctmp1(this%nwaves))
        allocate(this%ctmp2(this%nx,this%ny,this%nz))
        
        if (this%createRNGterms) then
            call this%initializeRNG()
        end if 
        
        this%isInitialized = .true.  
    end subroutine 

    subroutine getIsotropicModes(this,nk,ntheta,kmin,kmax,seed1,seed2,seed3)
        use gridtools, only: logspace
        use random, only: uniform_random, gaussian_random  
        class(rdt), intent(inout), target :: this
        logical, parameter :: DNSmode = .true.
        real(rkind), parameter :: eta_by_L = 0.001d0, cnu = 1.3d0
        integer, intent(in) :: nk, ntheta
        real(rkind), intent(in) :: kmin, kmax
        integer, intent(in), optional :: seed1, seed2, seed3

        real(rkind), dimension(:,:,:), allocatable, target :: randomNums
        real(rkind), dimension(:), allocatable, target :: k_abs, Ek, Sk
        real(rkind), dimension(:), allocatable :: k1, k2, k3, amag, bmag
        real(rkind), dimension(:,:), allocatable:: a, b, ua, ub
        integer :: idx, str_idx, end_idx

        real(rkind) :: dk, uabs, theta_iso 
        real(rkind), pointer :: mag, zeta(:), theta(:), alpha(:) 
        complex(rkind), pointer, dimension(:) :: uhat, vhat, what, thetahat

        ! Step 1: Generate shell radii
        allocate(k_abs(nk))
        k_abs = logspace(real(log10(kmin),rkind),real(log10(kmax),rkind),nk)
        
        ! Step 2: Allocate all the other local scope variables
        allocate(Ek(nk))
        if (this%isStratified) then
            allocate(randomNums(ntheta, nk, 6))
            allocate(Sk(nk)) 
        else
            allocate(randomNums(ntheta, nk, 4)) 
        end if 

        ! Step 3: Create the energy spectra for KE and PE
        Ek = C_vk*(k_abs**four)/((one + k_abs**two)**(17.0_rkind/six))
        if (DNSmode) then
            Ek = Ek*exp(-cnu*k_abs*eta_by_L)
        end if 
        if (this%isStratified) then
            Sk = C_vkt*(k_abs**two)/((one + k_abs**two)**(11.0_rkind/six))
            if (DNSmode) then
                Sk = Sk*exp(-cnu*k_abs*eta_by_L)
            end if 
        end if 

        ! Step 4: Generate the random numbers
        if (.not. present(seed3)) then 
            call uniform_random(randomNums(:,:,1),-one,one)  
            call uniform_random(randomNums(:,:,2:4),zero,two*pi)  
            if (this%isStratified) then
                call gaussian_random(randomNums(:,:,5:6),zero,one)
            end if 
        else
            call uniform_random(randomNums(:,:,1),-one,one,seed1)  
            call uniform_random(randomNums(:,:,2:4),zero,two*pi,seed2)  
            if (this%isStratified) then
                call gaussian_random(randomNums(:,:,5:6),zero,one,seed3)
            end if 
        end if 

        ! Step 5: Allocate local quantities for each shell
        allocate(k1(ntheta),k2(ntheta),k3(ntheta))
        allocate(a(ntheta,3),b(ntheta,3), amag(ntheta),bmag(ntheta))
        allocate(ua(ntheta,3),ub(ntheta,3))
        uhat => this%fieldhat(:,1)
        vhat => this%fieldhat(:,2)
        what => this%fieldhat(:,3)
        if (this%isStratified) thetahat => this%fieldhat(:,4)

        
        ! Step 6: Create the shells
        do idx = 1,nk
            mag => k_abs(idx)
            zeta => randomNums(:,idx,1)
            theta => randomNums(:,idx,2)
            k1 = mag*sqrt(1 - zeta**2)*cos(theta)
            k2 = mag*sqrt(1 - zeta**2)*sin(theta)
            k3 = mag*zeta

            if (idx .eq. 1) then
                dk = (k_abs(2) - k_abs(1))/two
            elseif (idx .eq. nk) then
                dk = (k_abs(nk) - k_abs(nk-1))/two
            else
                dk = (k_abs(idx+1) - k_abs(idx))/two + (k_abs(idx) - k_abs(idx-1))/two
            end if

            uabs = sqrt(two*Ek(idx)*dk/real(ntheta,rkind))
            a(:,1) =  zero
            a(:,2) = -k3
            a(:,3) =  k2
            b(:,1) =  k2**2 + k3**2
            b(:,2) = -k1*k2
            b(:,3) = -k1*k3
            amag = sqrt(a(:,1)**2 + a(:,2)**2 + a(:,3)**2)
            bmag = sqrt(b(:,1)**2 + b(:,2)**2 + b(:,3)**2)
            
            a(:,1) = a(:,1)/amag
            a(:,2) = a(:,2)/amag
            a(:,3) = a(:,3)/amag
            
            b(:,1) = b(:,1)/bmag
            b(:,2) = b(:,2)/bmag
            b(:,3) = b(:,3)/bmag

            alpha => randomNums(:,idx,3)
            ua(:,1) = uabs*(cos(alpha)*a(:,1) + sin(alpha)*b(:,1))
            ua(:,2) = uabs*(cos(alpha)*a(:,2) + sin(alpha)*b(:,2))
            ua(:,3) = uabs*(cos(alpha)*a(:,3) + sin(alpha)*b(:,3))

            alpha => randomNums(:,idx,4)
            ub(:,1) = uabs*(cos(alpha)*a(:,1) + sin(alpha)*b(:,1))
            ub(:,2) = uabs*(cos(alpha)*a(:,2) + sin(alpha)*b(:,2))
            ub(:,3) = uabs*(cos(alpha)*a(:,3) + sin(alpha)*b(:,3))

            str_idx = (idx - 1)*ntheta + 1
            end_idx = str_idx + ntheta - 1
            this%wavenums(str_idx:end_idx,1) = k1
            this%wavenums(str_idx:end_idx,2) = k2
            this%wavenums(str_idx:end_idx,3) = k3
            uhat(str_idx:end_idx) = cmplx(ua(:,1),ub(:,1),rkind) 
            vhat(str_idx:end_idx) = cmplx(ua(:,2),ub(:,2),rkind) 
            what(str_idx:end_idx) = cmplx(ua(:,3),ub(:,3),rkind) 

            if (this%isStratified) then
                theta_iso = dsqrt(two*Sk(idx)*dk/real(ntheta,rkind))
                thetahat(str_idx:end_idx) = cmplx(theta_iso*randomNums(:,idx,5),theta_iso*randomNums(:,idx,6),rkind) 
            end if 
        end do  


        ! Final Step: Deallocate all local scope variables
        deallocate(k_abs, Ek, ua, ub, amag, bmag, a, b)
        if (allocated(Sk)) deallocate(Sk)
        nullify(uhat, vhat, what, alpha, theta, zeta)
        if (associated(thetahat)) nullify(thetahat)
        deallocate(randomNums)
    end subroutine 

    subroutine destroy(this)
        class(rdt), intent(inout) :: this

        if (.not. this%isInitialized) then
            call GracefulExit("Calling destroy procedure for RDT derived type &
                & before initializing it", 1)
        end if 

        deallocate(this%fieldhat)
        deallocate(this%wavenums)
        deallocate(this%buf)
        deallocate(this%sweepFact)
        nullify(this%uhat, this%vhat,this%what,this%k1, this%k2, this%k3)
        this%isStratified = .false. 
        this%allowfastcalculations = .false. 
        this%isInitialized = .false.
        this%dispx = zero
        this%dispy = zero
        this%dispz = zero

        if (allocated(this%ctmp1)) deallocate(this%ctmp1)
        if (allocated(this%ctmp2)) deallocate(this%ctmp2)
        if (allocated(this%buf1)) deallocate(this%buf1)
        if (allocated(this%buf2)) deallocate(this%buf2)
        if (allocated(this%FT)) then
            call this%FT%destroy
            deallocate(this%FT)
        end if  
        
    end subroutine

    subroutine doSlowSum(this,x,y,z,u,v,w)
        class(rdt), intent(in) :: this
        real(rkind), intent(in) :: x, y, z
        real(rkind), intent(out) :: u, v, w
        integer :: idx
        real(rkind), dimension(this%nwaves) :: buf

        !this%uhat = this%uhat*this%sweepFact
        !this%vhat = this%vhat*this%sweepFact
        !this%what = this%what*this%sweepFact
        !if (this%isStratified) then
        !    this%thetahat = this%thetahat*this%sweepFact
        !end if 

        buf = this%k1*x + this%k2*y + this%k3*z
       
        u = zero
        do idx = 1,this%nwaves
            u = u + real(this%uhat(idx),rkind)*cos(buf(idx))
            u = u - dimag(this%uhat(idx))*sin(buf(idx))
        end do 

        v = zero
        do idx = 1,this%nwaves
            v = v + real(this%vhat(idx),rkind)*cos(buf(idx))
            v = v - dimag(this%vhat(idx))*sin(buf(idx))
        end do 

        w = zero
        do idx = 1,this%nwaves
            w = w + real(this%what(idx),rkind)*cos(buf(idx))
            w = w - dimag(this%what(idx))*sin(buf(idx))
        end do 
        
        !this%uhat = this%uhat/this%sweepFact
        !this%vhat = this%vhat/this%sweepFact
        !this%what = this%what/this%sweepFact
        !if (this%isStratified) then
        !    this%thetahat = this%thetahat/this%sweepFact
        !end if 

    end subroutine


    subroutine doFastSumPressure(this,phat,p)
        use constants, only: imi, zero 
        class(rdt), intent(inout), target :: this
        complex(rkind), dimension(this%nwaves), intent(inout) :: phat
        real(rkind), dimension(this%Nx,this%Ny,this%Nz), intent(out) :: p
        complex(rkind), dimension(:), pointer :: tempC1
        complex(rkind), dimension(:,:,:), pointer :: tempC2

        tempC1 => this%ctmp1
        tempC2 => this%ctmp2

        phat = phat*this%sweepFact
        tempC1 = exp(imi*(this%k1*this%factx + this%k2*this%facty + this%k3*this%factz))
        phat = phat*tempC1
        this%k1 = this%k1 * this%dx
        this%k2 = this%k2 * this%dy
        this%k3 = this%k3 * this%dz  

        if ( this%allowfastcalculations) then
            call this%FT%fft(this%nwaves,this%k1,this%k2,this%k3,phat,tempC2)
            p = real(tempC2,rkind)*this%nwaves
        else
            p = zero
        end if 

        this%k1 = this%k1/this%dx
        this%k2 = this%k2/this%dy
        this%k3 = this%k3/this%dz
        phat = phat/tempC1
        phat = phat/this%sweepFact


        nullify(tempC2, tempC1)
    end subroutine

    subroutine doFastSum(this,u,v,w,theta)
        use constants, only: imi 
        class(rdt), intent(inout), target :: this
        real(rkind), dimension(this%Nx,this%Ny,this%Nz), intent(out) :: u, v, w 
        real(rkind), dimension(this%Nx,this%Ny,this%Nz), intent(out), optional :: theta 
        complex(rkind), dimension(:), pointer :: tempC1
        complex(rkind), dimension(:,:,:), pointer :: tempC2

        tempC1 => this%ctmp1
        tempC2 => this%ctmp2

        
        this%uhat = this%uhat*this%sweepFact
        this%vhat = this%vhat*this%sweepFact
        this%what = this%what*this%sweepFact
        if (this%isStratified) then
            this%thetahat = this%thetahat*this%sweepFact
        end if 

        if ( this%allowfastcalculations) then
            tempC1 = exp(imi*(this%k1*this%factx + this%k2*this%facty + this%k3*this%factz))
            this%uhat = this%uhat * tempC1 
            this%vhat = this%vhat * tempC1
            this%what = this%what * tempC1
            this%k1 = this%k1 * this%dx
            this%k2 = this%k2 * this%dy
            this%k3 = this%k3 * this%dz  
            
            call this%FT%fft(this%nwaves,this%k1,this%k2,this%k3,this%uhat,tempC2)
            u = real(tempC2,rkind)*this%nwaves

            call this%FT%fft(this%nwaves,this%k1,this%k2,this%k3,this%vhat,tempC2)
            v = real(tempC2,rkind)*this%nwaves

            call this%FT%fft(this%nwaves,this%k1,this%k2,this%k3,this%what,tempC2)
            w = real(tempC2,rkind)*this%nwaves


            if ((present(theta)) .and. (this%isStratified)) then
                this%thetahat = this%thetahat * tempC1
                call this%FT%fft(this%nwaves,this%k1,this%k2,this%k3,this%thetahat,tempC2)
                theta = real(tempC2,rkind)*this%nwaves
                this%thetahat = this%thetahat / tempC1
            end if 

            this%k1 = this%k1/this%dx
            this%k2 = this%k2/this%dy
            this%k3 = this%k3/this%dz
            
            this%uhat = this%uhat / tempC1  
            this%vhat = this%vhat / tempC1
            this%what = this%what / tempC1

        else
            u = zero
            v = zero
            w = zero
        end if 
        
        this%uhat = this%uhat/this%sweepFact
        this%vhat = this%vhat/this%sweepFact
        this%what = this%what/this%sweepFact
        if (this%isStratified) then
            this%thetahat = this%thetahat/this%sweepFact
        end if 

        nullify(tempC2, tempC1)
    end subroutine 

    subroutine extractRange(this,wavenums_out, fieldhat_out, kleft, kright, success)
        class(rdt), intent(inout) :: this
        real(rkind), dimension(:,:), allocatable, intent(out) :: wavenums_out
        complex(rkind), dimension(:,:), allocatable, intent(out) :: fieldhat_out
        logical, intent(out), optional :: success
        real(rkind), intent(in) :: kleft, kright
        integer :: nwaves_out, idx, i
        integer, dimension(this%nwaves) :: marker

        marker = 0
        this%buf = sqrt(this%k1*this%k1 + this%k2*this%k2 + this%k3*this%k3)
        
        where (this%buf .ge. kleft)
            marker = 1
        end where

        where (this%buf .gt. kright) 
            marker = 0
        end where
        
        nwaves_out = sum(marker)

        if (nwaves_out < 1) then
            if (present(success)) success = .false. 
            return
        end if 
        allocate(wavenums_out(nwaves_out,3))
        allocate(fieldhat_out(nwaves_out,size(this%fieldhat,2)))

        idx = 0
        do i = 1,this%nwaves
            if (marker(i) == 1) then
                idx = idx + 1
                wavenums_out(idx,:) = this%wavenums(i,:)
                fieldhat_out(idx,:) = this%fieldhat(i,:)
            end if 
        end do 

        if (present(success)) success = .true. 

    end subroutine

    subroutine getduidxj_slow(this,duidxj,x,y,z)
        class(rdt), intent(in) :: this
        real(rkind), dimension(9), intent(out) :: duidxj
        real(rkind), intent(in) :: x, y, z

        integer :: idx
        real(rkind) :: tmp
        real(rkind), dimension(this%nwaves) :: buf
        
        real(rkind), dimension(:), pointer :: rtmp
        complex(rkind), dimension(:), pointer :: ctmp


        this%uhat = this%uhat*this%sweepFact
        this%vhat = this%vhat*this%sweepFact
        this%what = this%what*this%sweepFact
        
        buf = this%k1*x + this%k2*y + this%k3*z
        
        ! dudx
        rtmp => this%k1; ctmp => this%uhat; tmp = zero
        do idx = 1,this%nwaves
            tmp = tmp - rtmp(idx)* real(ctmp(idx),rkind)*sin(buf(idx))
            tmp = tmp - rtmp(idx)*dimag(ctmp(idx))*cos(buf(idx))
        end do 
        duidxj(1) = tmp 

        ! dudy
        rtmp => this%k2; ctmp => this%uhat; tmp = zero
        do idx = 1,this%nwaves
            tmp = tmp - rtmp(idx)* real(ctmp(idx),rkind)*sin(buf(idx))
            tmp = tmp - rtmp(idx)*dimag(ctmp(idx))*cos(buf(idx))
        end do 
        duidxj(2) = tmp 

        ! dudz
        rtmp => this%k3; ctmp => this%uhat; tmp = zero
        do idx = 1,this%nwaves
            tmp = tmp - rtmp(idx)* real(ctmp(idx),rkind)*sin(buf(idx))
            tmp = tmp - rtmp(idx)*dimag(ctmp(idx))*cos(buf(idx))
        end do 
        duidxj(3) = tmp 

        ! dvdx
        rtmp => this%k1; ctmp => this%vhat; tmp = zero
        do idx = 1,this%nwaves
            tmp = tmp - rtmp(idx)* real(ctmp(idx),rkind)*sin(buf(idx))
            tmp = tmp - rtmp(idx)*dimag(ctmp(idx))*cos(buf(idx))
        end do 
        duidxj(4) = tmp 

        ! dvdy
        rtmp => this%k2; ctmp => this%vhat; tmp = zero
        do idx = 1,this%nwaves
            tmp = tmp - rtmp(idx)* real(ctmp(idx),rkind)*sin(buf(idx))
            tmp = tmp - rtmp(idx)*dimag(ctmp(idx))*cos(buf(idx))
        end do 
        duidxj(5) = tmp 

        ! dvdz
        rtmp => this%k3; ctmp => this%vhat; tmp = zero
        do idx = 1,this%nwaves
            tmp = tmp - rtmp(idx)* real(ctmp(idx),rkind)*sin(buf(idx))
            tmp = tmp - rtmp(idx)*dimag(ctmp(idx))*cos(buf(idx))
        end do 
        duidxj(6) = tmp 

        ! dwdx
        rtmp => this%k1; ctmp => this%what; tmp = zero
        do idx = 1,this%nwaves
            tmp = tmp - rtmp(idx)* real(ctmp(idx),rkind)*sin(buf(idx))
            tmp = tmp - rtmp(idx)*dimag(ctmp(idx))*cos(buf(idx))
        end do 
        duidxj(7) = tmp 

        ! dwdy
        rtmp => this%k2; ctmp => this%what; tmp = zero
        do idx = 1,this%nwaves
            tmp = tmp - rtmp(idx)* real(ctmp(idx),rkind)*sin(buf(idx))
            tmp = tmp - rtmp(idx)*dimag(ctmp(idx))*cos(buf(idx))
        end do 
        duidxj(8) = tmp 

        ! dwdz
        rtmp => this%k3; ctmp => this%what; tmp = zero
        do idx = 1,this%nwaves
            tmp = tmp - rtmp(idx)* real(ctmp(idx),rkind)*sin(buf(idx))
            tmp = tmp - rtmp(idx)*dimag(ctmp(idx))*cos(buf(idx))
        end do 
        duidxj(9) = tmp 

        nullify (rtmp, ctmp)
        
        this%uhat = this%uhat/this%sweepFact
        this%vhat = this%vhat/this%sweepFact
        this%what = this%what/this%sweepFact
    end subroutine

    subroutine getduidxj_fast(this,duidxj)
        use constants, only: imi
        class(rdt), intent(inout) :: this
        real(rkind), dimension(this%Nx,this%Ny,this%Nz,9), target, intent(out) :: duidxj
        complex(rkind), dimension(this%nwaves) :: tempC1, tempC3
        complex(rkind), dimension(this%Nx,this%Ny,this%Nz) :: tempC2
        real(rkind), dimension(:,:,:), pointer :: dudx, dudy, dudz, dvdx, dvdy, dvdz,dwdx, dwdy, dwdz


        this%uhat = this%uhat*this%sweepFact
        this%vhat = this%vhat*this%sweepFact
        this%what = this%what*this%sweepFact

        dudx => duidxj(:,:,:,1); dudy => duidxj(:,:,:,2); dudz => duidxj(:,:,:,3)
        dvdx => duidxj(:,:,:,4); dvdy => duidxj(:,:,:,5); dvdz => duidxj(:,:,:,6)
        dwdx => duidxj(:,:,:,7); dwdy => duidxj(:,:,:,8); dwdz => duidxj(:,:,:,9)
        
        if ( this%allowfastcalculations) then
            tempC1 = exp(imi*(this%k1*this%factx + this%k2*this%facty + this%k3*this%factz))
           
            ! dudx  
            tempC3  = imi*this%k1*this%uhat
            tempC3  = tempC3*tempC1
            this%k1 = this%k1 * this%dx
            this%k2 = this%k2 * this%dy
            this%k3 = this%k3 * this%dz             
            call this%FT%fft(this%nwaves,this%k1,this%k2,this%k3,tempC3,tempC2)
            this%k1 = this%k1/this%dx
            this%k2 = this%k2/this%dy
            this%k3 = this%k3/this%dz
            dudx = real(tempC2,rkind)*this%nwaves
            
            ! dudy  
            tempC3  = imi*this%k2*this%uhat
            tempC3  = tempC3*tempC1
            this%k1 = this%k1 * this%dx
            this%k2 = this%k2 * this%dy
            this%k3 = this%k3 * this%dz             
            call this%FT%fft(this%nwaves,this%k1,this%k2,this%k3,tempC3,tempC2)
            this%k1 = this%k1/this%dx
            this%k2 = this%k2/this%dy
            this%k3 = this%k3/this%dz
            dudy = real(tempC2,rkind)*this%nwaves
           
            ! dudz  
            tempC3  = imi*this%k3*this%uhat
            tempC3  = tempC3*tempC1
            this%k1 = this%k1 * this%dx
            this%k2 = this%k2 * this%dy
            this%k3 = this%k3 * this%dz             
            call this%FT%fft(this%nwaves,this%k1,this%k2,this%k3,tempC3,tempC2)
            this%k1 = this%k1/this%dx
            this%k2 = this%k2/this%dy
            this%k3 = this%k3/this%dz
            dudz = real(tempC2,rkind)*this%nwaves

            ! dvdx  
            tempC3  = imi*this%k1*this%vhat
            tempC3  = tempC3*tempC1
            this%k1 = this%k1 * this%dx
            this%k2 = this%k2 * this%dy
            this%k3 = this%k3 * this%dz             
            call this%FT%fft(this%nwaves,this%k1,this%k2,this%k3,tempC3,tempC2)
            this%k1 = this%k1/this%dx
            this%k2 = this%k2/this%dy
            this%k3 = this%k3/this%dz
            dvdx = real(tempC2,rkind)*this%nwaves

            ! dvdy  
            tempC3  = imi*this%k2*this%vhat
            tempC3  = tempC3*tempC1
            this%k1 = this%k1 * this%dx
            this%k2 = this%k2 * this%dy
            this%k3 = this%k3 * this%dz             
            call this%FT%fft(this%nwaves,this%k1,this%k2,this%k3,tempC3,tempC2)
            this%k1 = this%k1/this%dx
            this%k2 = this%k2/this%dy
            this%k3 = this%k3/this%dz
            dvdy = real(tempC2,rkind)*this%nwaves

            ! dvdz  
            tempC3  = imi*this%k3*this%vhat
            tempC3  = tempC3*tempC1
            this%k1 = this%k1 * this%dx
            this%k2 = this%k2 * this%dy
            this%k3 = this%k3 * this%dz             
            call this%FT%fft(this%nwaves,this%k1,this%k2,this%k3,tempC3,tempC2)
            this%k1 = this%k1/this%dx
            this%k2 = this%k2/this%dy
            this%k3 = this%k3/this%dz
            dvdz = real(tempC2,rkind)*this%nwaves

            ! dwdx 
            tempC3  = imi*this%k1*this%what
            tempC3  = tempC3*tempC1
            this%k1 = this%k1 * this%dx
            this%k2 = this%k2 * this%dy
            this%k3 = this%k3 * this%dz             
            call this%FT%fft(this%nwaves,this%k1,this%k2,this%k3,tempC3,tempC2)
            this%k1 = this%k1/this%dx
            this%k2 = this%k2/this%dy
            this%k3 = this%k3/this%dz
            dwdx = real(tempC2,rkind)*this%nwaves

            ! dwdy 
            tempC3  = imi*this%k2*this%what
            tempC3  = tempC3*tempC1
            this%k1 = this%k1 * this%dx
            this%k2 = this%k2 * this%dy
            this%k3 = this%k3 * this%dz             
            call this%FT%fft(this%nwaves,this%k1,this%k2,this%k3,tempC3,tempC2)
            this%k1 = this%k1/this%dx
            this%k2 = this%k2/this%dy
            this%k3 = this%k3/this%dz
            dwdy = real(tempC2,rkind)*this%nwaves

            ! dwdz 
            tempC3  = imi*this%k3*this%what
            tempC3  = tempC3*tempC1
            this%k1 = this%k1 * this%dx
            this%k2 = this%k2 * this%dy
            this%k3 = this%k3 * this%dz             
            call this%FT%fft(this%nwaves,this%k1,this%k2,this%k3,tempC3,tempC2)
            this%k1 = this%k1/this%dx
            this%k2 = this%k2/this%dy
            this%k3 = this%k3/this%dz
            dwdz = real(tempC2,rkind)*this%nwaves
        else
            duidxj = zero
        end if 

        nullify(dudx, dudy, dudz, dvdx, dvdy, dvdz,dwdx, dwdy, dwdz)

        this%uhat = this%uhat/this%sweepFact
        this%vhat = this%vhat/this%sweepFact
        this%what = this%what/this%sweepFact
    end subroutine


    subroutine sweep(this,xmove, ymove, zmove)
        use constants, only: imi
        class(rdt), intent(inout) :: this
        real(rkind), intent(in) :: xmove, ymove, zmove
        !complex(rkind) :: tmp
        !integer :: idx2

        this%dispx = this%dispx + xmove
        this%dispy = this%dispy + ymove
        this%dispz = this%dispz + zmove
        
        this%sweepFact = exp(-imi*(this%k1*this%dispx + this%k2*this%dispy + this%k3*this%dispz))

    end subroutine

    
end module 
