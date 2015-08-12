module poisson
    use kind_parameters, only: rkind
    use fft_3d_stuff, only: fft_3d 
    use decomp_2d, only: decomp_info
    use exits, only: GracefulExit, message
    use derivativesMod, only: derivatives 
    use filtersMod, only: filters

    implicit none
    private
    public :: poisson

    type poisson
        real(rkind), dimension(:,:,:), allocatable :: k1_inX, k2_inX, k3_inX, oneByksq_inX
        real(rkind), dimension(:,:,:), allocatable :: k1_inZ, k2_inZ, k3_inZ, oneByksq_inZ
        type(fft_3d), allocatable :: FTx, FTz

        logical :: initialized = .false. 
        contains
        procedure :: init
        procedure :: destroy
        
        procedure :: projectPressure_X
        procedure :: projectPressure_Z
        
    end type

subroutine init(this, gp, dx, dy, dz, der, fil, exhaustive)
    class(poisson), intent(inout) :: this
    class(derivatives), intent(in) :: der
    class(filters), intent(in) :: fil 
    class(decomp_info), intent(in) :: gp 
    integer, intent(in) :: nx, ny, nz
    character(len=1), intent(in) :: base_dec 
    real(rkind), intent(in) :: dx, dy, dz
    logical, intent(in) :: exhaustive


    integer :: ierr
   
    if (this%initialized) then
        call message("WARNING: POISSON class already initialized. Destroying &
        & previous definition and reinitializing")
        call this%destroy
    end if 

    allocate(FTx)
    ierr = FTx%init(gp%xsz(1),gp%ysz(2),gp%zsz(3),"x",dx, dy, dz, exhaustive)

    if (ierr .ne. 0) then
        call GracefulExit("FFT_3 class inside POISSON class could not be initialized",1434)
    end if 
    
    allocate(FTz)
    ierr = FTz%init(gp%xsz(1),gp%ysz(2),gp%zsz(3),"z",dx, dy, dz, exhaustive)

    if (ierr .ne. 0) then
        call GracefulExit("FFT_3 class inside POISSON class could not be initialized",1432)
    end if 

    allocate(this%k1_inX(size(FTx%k1,1),size(FTx%k1,2),size(FTx%k1,3)))
    allocate(this%k2_inX(size(FTx%k2,1),size(FTx%k2,2),size(FTx%k2,3)))
    allocate(this%k3_inX(size(FTx%k3,1),size(FTx%k3,2),size(FTx%k3,3)))
    allocate(this%k1_inZ(size(FTz%k1,1),size(FTz%k1,2),size(FTz%k1,3)))
    allocate(this%k2_inZ(size(FTz%k2,1),size(FTz%k2,2),size(FTz%k2,3)))
    allocate(this%k3_inZ(size(FTz%k3,1),size(FTz%k3,2),size(FTz%k3,3)))

    select case (der%getMethodx)
    case ("cd10")
        call genKCD10(FTx,this%k1_inX)
    case ("cd06")
        call genKCD06(FTx,this%k1_inX)
    case ("four")
        call genKFOUR(FTx,this%k1_inX)
    case ("cheb")
        call GracefulExit("Chebyshev is incomplete", 1021)
    end select 
    
    
end subroutine 


end module 
