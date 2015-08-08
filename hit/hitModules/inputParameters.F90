module inputParameters
    use kind_parameters, only: rkind
    use constants,       only: pi, two
    implicit none
    
    ! Initial Conditions
    logical             :: initfromHIT3d = .true.       ! Get initial Conditions from HIT3d?
    character(len=64)   :: inputDir      = "./inputData"! Input Directory in which hit3d files are stored

    ! Grid Parameters
    integer             :: Nx            =  128     ! Unused if initFromHIT3d == .true.
    integer             :: Ny            =  128     ! Unused if initFromHIT3d == .true.
    integer             :: Nz            =  128     ! Unused if initFromHIT3d == .true.
    
    real(rkind)         :: Lx            = two*pi   ! Domain size in x 
    real(rkind)         :: Ly            = two*pi   ! Domain size in x 
    real(rkind)         :: Lz            = two*pi   ! Domain size in x 
   
    ! Simulation Time parameters 
    integer             :: tStop         =  1000    !Number of time steps
    logical             :: useCFL        =  .false. !If true, then variable time stepping is used 
    real(rkind)         :: CFL           =  0.1     ! CFL number used for variable time stepping (used if useCFL == .true.)
    real(rkind)         :: dt            =  0.001   ! Time step (used if useCFL == .false.)

    ! Subgrid Scale models
    logical             :: LESmodeON     = .false.
    character(len=64)   :: SubgridModel  = "STSMAG" ! Options: 1. STSMAG (standard smagorinsky) 2. DYSMAG (Dynamic Smagorinsky)

    ! Number of active fields
    integer             :: numScalars    = 0        ! Number of passive scalars

end module 

