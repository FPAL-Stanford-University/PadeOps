module variables
    use kind_parameters, only: rkind
    use fft_3d_Stuff, only: fft_3d 
    use decomp_2d,    only: decomp_info
    use constants,    only: pi, two

    integer                                             :: Nx
    integer                                             :: Ny
    integer                                             :: Nz
    real(rkind)                                         :: dx
    real(rkind)                                         :: dy
    real(rkind)                                         :: dz

    integer                                             :: NT
    logical                                             :: useCFL
    real(rkind)                                         :: CFL
    real(rkind)                                         :: dt 

    real(rkind)                                         :: REY
    real(rkind)                                         :: nu

    logical                                             :: RESTART  
    integer                                             :: Dealias
    logical                                             :: useHit3dinit
    character(len=128)                                  :: HIT3dinitDir
    character(len=128)                                  :: RESTARTDir
    character(len=128)                                  :: OutputDir

    logical                                             :: FFT_PLAN_EXHAUSTIVE
    real(rkind)   , dimension(:,:,:,:), allocatable     :: fieldsPhys
    complex(rkind), dimension(:,:,:,:), allocatable     :: fieldsSpec
    real(rkind)   , dimension(:,:,:)  , allocatable     :: DealiasMat
    real(rkind)   , dimension(:,:,:)  , allocatable     :: oneBykSq


    real(rkind)                                         :: kdealias_x, kdealias_y, kdealias_z 

    real(rkind)                                         :: maxDivergence 

    real(rkind)                                         :: time

    integer                                             :: TSTEP_RESTART
    integer                                             :: TSTEP_DUMP

    integer                                             :: runIDX

    complex(rkind)   , dimension(:,:,:,:), allocatable  :: rhs 
    complex(rkind)   , dimension(:,:,:,:), allocatable  :: rhsOld

    type(decomp_info), allocatable                      :: gp 
    type(fft_3d)     , allocatable                      :: ft
    
    
    real(rkind)     , dimension(:,:,:)  , allocatable   :: Buff_Real
    complex(rkind)  , dimension(:,:,:), allocatable     :: Buff_Cmplx

end module
