module standard_smagorinsky
    use kind_parameters, only: rkind
    use spectralMod, only: spectral
    use decomp_2d, only: decomp_info
    use constants, only: zero, three, imi, half, two 
    implicit none
    private

    public :: smag

    complex(rkind), parameter :: himi = half*imi 

    real(rkind), parameter :: c_s = 0.17_rkind 
    type :: smag
        private

        real(rkind), dimension(:,:,:), pointer :: k1, k2, k3 

        real(rkind), dimension(:,:,:,:), allocatable :: rbuff
        complex(rkind), dimension(:,:,:,:), allocatable :: cbuff
       
        real(rkind), dimension(:,:,:), allocatable :: S_ij 

        type(spectral), pointer :: spec

        real(rkind) :: deltaF

        logical :: isZperiodic

        integer :: nxhat, nyhat, nzhat
        integer :: nx, ny, nz

        contains
            procedure :: init
            procedure :: destroy
            procedure :: get_source 
    end type 

contains

    subroutine init(this, spec, decomp, dx, dy, dz, isZperiodic)
        class(smag), intent(inout) :: this
        class(spectral), intent(in), target :: spec
        class(decomp_info) :: decomp
        real(rkind), intent(in) :: dx, dy, dz
        logical, intent(in) :: isZperiodic

        this%isZperiodic = isZperiodic
        this%deltaF = (dx*dy*dz)**(1/three)
        
        this%k1 => spec%k1
        this%k2 => spec%k2
        this%k3 => spec%k3

        this%nx = decomp%xsz(1)
        this%ny = decomp%xsz(2)
        this%nz = decomp%xsz(3)
        this%nxhat = size(spec%kabs_sq,1)
        this%nyhat = size(spec%kabs_sq,2)
        this%nzhat = size(spec%kabs_sq,3)

        this%spec => spec

        allocate(this%rbuff(this%nx,this%ny,this%nz,2))
        allocate(this%cbuff(this%nxhat,this%nyhat,this%nzhat,2))

    end subroutine


    subroutine destroy(this)
        class(smag), intent(inout) :: this
        deallocate(this%rbuff)
        deallocate(this%cbuff)

        nullify(this%k1)
        nullify(this%k3)
        nullify(this%k2)
        nullify(this%spec)
    end subroutine 

    subroutine get_source(this, Sfield, duidxj, sourceHat, nut)
        class(smag), intent(inout), target :: this
        complex(rkind), dimension(this%nxhat, this%nyhat, this%nzhat,3), intent(in), target :: Sfield
        complex(rkind), dimension(this%nxhat, this%nyhat, this%nzhat,9), intent(in), target :: duidxj 
        complex(rkind), dimension(this%nxhat, this%nyhat, this%nzhat,3), intent(out) :: sourceHat
        complex(rkind), dimension(this%nx, this%ny, this%nz,3), intent(out) :: nut
        real(rkind), dimension(:,:,:), pointer :: rtmp1, rtmp2
        complex(rkind), dimension(:,:,:), pointer :: ctmp1, ctmp2 
        complex(rkind), dimension(:,:,:), pointer :: uhat, vhat, what
        

        uhat => Sfield(:,:,:,1)
        vhat => Sfield(:,:,:,2)
        what => Sfield(:,:,:,3)


        rtmp1 => this%rbuff(:,:,:,1)
        rtmp2 => this%rbuff(:,:,:,2)
        ctmp1 => this%cbuff(:,:,:,1)
        ctmp2 => this%cbuff(:,:,:,2)

        nut = zero  


        sourceHat = Sfield - Sfield 

    end subroutine 
end module 
