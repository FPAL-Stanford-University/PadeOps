    
    ! 6th order first derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha06d1=  1.0_rkind / 3.0_rkind
    real(rkind), parameter :: a06d1    = (14.0_rkind / 9.0_rkind) / 2.0_rkind
    real(rkind), parameter :: b06d1    = ( 1.0_rkind / 9.0_rkind) / 4.0_rkind
    
    ! 6th order second derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha06d2=  2.0_rkind / 11.0_rkind
    real(rkind), parameter :: a06d2    = (12.0_rkind / 11.0_rkind) / 1.0_rkind
    real(rkind), parameter :: b06d2    = ( 3.0_rkind / 11.0_rkind) / 4.0_rkind
    
    ! 10th order first derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha10d1=  1.0_rkind /  2.0_rkind
    real(rkind), parameter :: beta10d1 =  1.0_rkind / 20.0_rkind
    real(rkind), parameter :: a10d1    =( 17.0_rkind / 12.0_rkind) / 2.0_rkind
    real(rkind), parameter :: b10d1    =(101.0_rkind /150.0_rkind) / 4.0_rkind
    real(rkind), parameter :: c10d1    =(  1.0_rkind /100.0_rkind) / 6.0_rkind
    
    ! 10th order second derivative coefficients (See Lele (1992) for explanation)
    real(rkind), parameter :: alpha10d2=334.0_rkind /899.0_rkind
    real(rkind), parameter :: beta10d2 = 43.0_rkind /1798.0_rkind
    real(rkind), parameter :: a10d2    =(1065.0_rkind /1798.0_rkind) / 1.0_rkind
    real(rkind), parameter :: b10d2    =(1038.0_rkind / 899.0_rkind) / 4.0_rkind
    real(rkind), parameter :: c10d2    =( 79.0_rkind /1798.0_rkind ) / 9.0_rkind
