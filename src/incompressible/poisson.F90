module poissonMod
    use kind_parameters, only: rkind
    use spectralMod, only: spectral, GetWaveNums 
    use constants, only: zero, one, two
    use exits, only: GracefulExit
    use decomp_2d

    implicit none 

    private
    public :: poisson

    type :: poisson
        private
        integer :: nx_in, ny_in, nz_in
        complex(rkind), dimension(:,:,:,:), allocatable :: tmp
        real(rkind), dimension(:,:,:), allocatable :: mfact1, mfact2, mfact3
        logical :: zperiodic = .true. 

        type(decomp_info), pointer, public :: sp_gp

        ! Stuff for thomas algorithm
        integer :: nz_inZ, nx_inZ, ny_inZ
        real(rkind), dimension(:), allocatable :: k1_inZ, k2_inZ
        real(rkind) :: dz, dzsq
        integer :: xst, xen, yst, yen, zst, zen
            
        contains
            procedure :: init
            procedure :: destroy
            procedure :: PressureProj
            procedure :: PoissonSolveZ    
            procedure :: allocArrZ    
    end type 

contains

    subroutine allocArrZ(this,arr)
        class(poisson), intent(in) :: this
        complex(rkind),dimension(:,:,:), allocatable, intent(out) :: arr
        type(decomp_info), pointer :: gp

        gp => this%sp_gp
        allocate(arr(gp%zsz(1),gp%zsz(2),gp%zsz(3)))
    
    end subroutine


    subroutine init(this,spect,PeriodicVertical,dx,dy,dz)
        class(poisson), intent(inout) :: this
        class(spectral), intent(in), target :: spect
        logical, intent(in), optional :: PeriodicVertical
        real(rkind), optional :: dx, dy, dz
       
        if (present(PeriodicVertical)) then
            this%zperiodic = PeriodicVertical
        end if  

        if (this%zperiodic) then
            call spect%alloc_r2c_out(this%tmp,2)
            this%nx_in = size(this%tmp,1)
            this%ny_in = size(this%tmp,2)
            this%nz_in = size(this%tmp,3)

            allocate(this%mfact1(size(spect%k1,1),size(spect%k1,2),size(spect%k1,3)))
            allocate(this%mfact2(size(spect%k1,1),size(spect%k1,2),size(spect%k1,3)))
            allocate(this%mfact3(size(spect%k1,1),size(spect%k1,2),size(spect%k1,3)))

            this%mfact1 = spect%k1*spect%one_by_kabs_sq
            this%mfact2 = spect%k2*spect%one_by_kabs_sq
            this%mfact3 = spect%k3*spect%one_by_kabs_sq
        else
            if (.not.present(dz)) call GracefulExit("Need to send in dz as input to poisson init",31)
            this%dz = dz 
            this%dzsq = dz**2
            this%nz_inZ = spect%nz_g
            this%nx_inZ = size(spect%k1,1)
            this%ny_inZ = size(spect%k1,2)
            allocate(this%k1_inZ(spect%nx_g))
            allocate(this%k2_inZ(spect%ny_g))
            this%k1_inZ = GetWaveNums(spect%nx_g,dx)
            this%k2_inZ = GetWaveNums(spect%ny_g,dy)
            this%sp_gp => spect%spectdecomp
            this%xst = this%sp_gp%zst(1)
            this%xen = this%sp_gp%zen(1)
            this%yst = this%sp_gp%zst(2)
            this%yen = this%sp_gp%zen(2)
            this%zst = this%sp_gp%zst(3)
            this%zen = this%sp_gp%zen(3)
        end if 

    end subroutine

    subroutine destroy(this)
        class(poisson), intent(inout) :: this
        
        if (allocated(this%tmp)) deallocate(this%tmp) 
        if (allocated(this%mfact1)) deallocate(this%mfact1) 
        if (allocated(this%mfact2)) deallocate(this%mfact2) 
        if (allocated(this%mfact3)) deallocate(this%mfact3) 
        if (allocated(this%k1_inZ)) deallocate(this%k1_inZ)
        if (allocated(this%k2_inZ)) deallocate(this%k2_inZ)
    end subroutine


    subroutine PressureProj(this,Sfields,spect)
        class(poisson), target, intent(inout) :: this
        class(spectral), target, intent(inout) :: spect
        complex(rkind), target, dimension(this%nx_in,this%ny_in,this%nz_in,3), intent(inout) :: Sfields
        complex(rkind), dimension(:,:,:), pointer :: tmp1, tmp2, uhat, vhat, what
        real(rkind), dimension(:,:,:), pointer :: k1, k2, k3

        tmp1 => this%tmp(:,:,:,1)
        tmp2 => this%tmp(:,:,:,2)
        uhat => Sfields(:,:,:,1)
        vhat => Sfields(:,:,:,2)
        what => Sfields(:,:,:,3)
        k1 => spect%k1
        k2 => spect%k2
        k3 => spect%k3

        if (this%zperiodic) then
            ! STEP 1: COMPUTE DIVERGENCE 
            tmp1 = k1*uhat 
            tmp1 = tmp1 + k2*vhat
            tmp1 = tmp1 + k3*what 

            ! STEP 2: PROJECT U VELOCITY
            tmp2 = this%mfact1*tmp1
            uhat = uhat - tmp2

            ! STEP 3: PROJECT V VELOCITY
            tmp2 = this%mfact2*tmp1
            vhat = vhat - tmp2

            ! STEP 4: PROJECT U VELOCITY
            tmp2 = this%mfact3*tmp1
            what = what - tmp2
        else
            print*, "Incomplete"    
        end if 

        ! DONE - the new velocity fields is divergence free
    end subroutine

    subroutine PoissonSolveZ(this,fhat,phat)
        ! Assuming that everything is in z-decomp
        class(poisson), intent(in) :: this 
        complex(rkind), dimension(this%xst:this%xen,this%yst:this%yen,this%zst:this%zen), intent(in) :: fhat
        complex(rkind), dimension(this%xst:this%xen,this%yst:this%yen,this%zst:this%zen), intent(out)::phat
        real(rkind), dimension(this%nz_inZ) :: a, b, c
        integer :: i, j
        real(rkind) :: k1, k2, aa       
        complex(rkind), dimension(this%nz_inZ) :: y, rhs
       
        do j = this%yst,this%yen
            do i = this%xst,this%xen
                k1 = this%k1_inZ(i)
                k2 = this%k2_inZ(j)
                
                rhs = this%dzsq*fhat(i,j,:)
                
                if ((i == 1).and. (j == 1)) then
                    call genTridiag2ndOrder(k1,k2,this%dz,this%nz_inZ,a,b,c)
                    a(1) = one; b(1) = zero; c(1) = zero
                    rhs(1) = zero
                    call solveTridiag_Poiss(a,b,c,rhs,y,this%nz_inZ)
                    phat(i,j,:) = y
                else
                    aa = (-(k1*this%dz)**2 - (k2*this%dz)**2 - two)
                    call solveTridiag_Poiss_InPlace_quick(aa,one,one,rhs,this%nz_inZ)
                    phat(i,j,:) = rhs
                end if 
            end do 
        end do 

    end subroutine    


    pure subroutine genTridiag2ndOrder(k1,k2,dz,n,a,b,c)
        real(rkind), intent(in) :: k1, k2, dz
        integer, intent(in) :: n
        real(rkind), dimension(n), intent(out) :: a, b, c

        a = (-(k1*dz)**2 - (k2*dz)**2 - 2)
        b = 1._rkind
        c = 1._rkind
        a(1) = a(1) + 1
        a(n) = a(n) + 1
        
    end subroutine


    pure subroutine solveTridiag_Poiss(a,b,c,f,y,n) 
        integer,intent(in) :: n
        real(rkind),dimension(n),intent(in) :: a,b,c
        complex(rkind), dimension(n), intent(in) :: f
        complex(rkind), dimension(n), intent(out) :: y
        real(rkind),dimension(n) :: v
        real(rkind) :: w
        integer :: i, j
        
        w = a(1)
        y(1) = f(1)/w
        
        do i = 2,n
            v(i-1) = c(i-1)/w
            w = a(i) - b(i)*v(i-1)
            y(i) = ( f(i) - b(i)*y(i-1))/w
        end do 

        do j = n-1,1,-1
            y(j) = y(j) - v(j)*y(j+1)
        end do 

    end subroutine
    
    subroutine solveTridiag_Poiss_InPlace_quick(a,b,c,f,n) 
        integer,intent(in) :: n
        real(rkind),intent(in) :: a,b,c
        complex(rkind), dimension(n), intent(inout) :: f
        !complex(rkind), dimension(n), intent(out) :: y
        real(rkind),dimension(n) :: v
        real(rkind) :: w
        integer :: i, j
        
        w = a + one
        f(1) = f(1)/w
        
        do i = 2,n-1
            v(i-1) = c/w
            w = a - b*v(i-1)
            f(i) = ( f(i) - b*f(i-1))/w
        end do
        v(n-1) = c/w
        w = (a + one) - b*v(n-1)
        f(n) = (f(n) - b*f(n-1))/w 

        do j = n-1,1,-1
            f(j) = f(j) - v(j)*f(j+1)
        end do 

    end subroutine
    
end module 
