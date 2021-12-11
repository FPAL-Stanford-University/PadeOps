program test_gaborMode_MMS
    use kind_parameters, only: rkind
    use GaborModeMod, only: gaborMode, MMS, nutMMS

    implicit none
    type(gaborMode) :: gm
    complex(kind=8), parameter :: uhat = (1.d0,1.d0), vhat = (0.d0,0.d0), what = (1.d0,1.d0)
    real(rkind), parameter :: kx = 2.d0, ky = 0.d0, kz = 0.d0
    real(rkind), parameter :: xMode = 1.d0, yMode = 1.d0, zMode = 1.d0
    real(rkind), parameter :: wSupport = 1.d0, kmin = 0.d0, epsKE = 1.d0
    real(rkind), parameter :: U = 1.d0, V = 1.d0, W = 1.d0
    real(rkind), parameter :: dudx = 1.1d0, dudy = 1.2d0, dudz = 1.3d0
    real(rkind), parameter :: dvdx = 2.1d0, dvdy = 2.2d0, dvdz = 2.3d0
    real(rkind), parameter :: dwdx = 3.1d0, dwdy = 3.2d0, dwdz = 3.3d0
    real(rkind) :: g1, g2, g3, f1, f2, f3
    real(rkind), parameter :: dt = 0.01d0
    real(rkind), parameter :: t0 = 0.d0, tf = 5.d0
    integer :: tsteps
     
    MMS = .true.
    nutMMS = 0.1d0
    call gm%init(uhat,vhat,what,xMode,yMode,zMode,kx,ky,kz,wSupport)
 
    tsteps = nint((tf - t0)/dt)
    do tid = 1,tsteps
        call gm%evolve(U, V, W, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, epsKE, dt, kmin)
    end do


end program test_gaborMode_MMS
