module poissonMod
    use kind_parameters, only: rkind
    use spectralMod, only: spectral 

    implicit none 

    private
    public :: poisson

    type :: poisson
        private
        integer :: nx_in, ny_in, nz_in
        complex(rkind), dimension(:,:,:,:), allocatable :: tmp
        real(rkind), dimension(:,:,:), allocatable :: mfact1, mfact2, mfact3

        contains
            procedure :: init
            procedure :: destroy
            procedure :: PressureProj
            
    end type 

contains

    subroutine init(this,spect)
        class(poisson), intent(inout) :: this
        class(spectral), intent(in) :: spect

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

    end subroutine

    subroutine destroy(this)
        class(poisson), intent(inout) :: this
        
        if (allocated(this%tmp)) deallocate(this%tmp) 
        if (allocated(this%mfact1)) deallocate(this%mfact1) 
        if (allocated(this%mfact2)) deallocate(this%mfact2) 
        if (allocated(this%mfact3)) deallocate(this%mfact3) 

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

        ! DONE - the new velocity fields is divergence free
    end subroutine
end module 
