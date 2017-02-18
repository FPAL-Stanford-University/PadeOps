program test_Pade6ops
    use kind_parameters, only: rkind, clen
    use decomp_2d 
    use decomp_2d_io
    use spectralMod, only: spectral
    use mpi
    use constants
    use PadeDerOps, only: Pade6stagg
    use cd06staggstuff, only: cd06stagg
    use exits
    implicit none

    real(rkind), dimension(:,:,:), allocatable :: uC, zC, zE, uE, d2udz2C, dudzE, d2udz2E
    type(decomp_info) :: gpC
    type(decomp_info), pointer :: sp_gpC
    type(spectral), target :: spectC
    integer :: nx = 16, ny = 16, nz = 64
    integer :: ierr, k 
    real(rkind) :: dx, dy, dz
    type(Pade6stagg) :: Pade6opsZ
    type(cd06stagg) :: cd06ViscOp

    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(gpC)

    if (nproc>1) call gracefulExit("Meant to be a serial execution",13)
    
    dx = two*pi/real(nx,rkind)
    dy = two*pi/real(ny,rkind)
    dz = two/real(nz,rkind)

    call spectC%init("x", nx, ny, nz  , dx, dy, dz, "four", "2/3rd", 2, .false.)
    sp_gpC => spectC%spectdecomp

    allocate(uC(nx,ny,nz), zC(nx,ny,nz), d2udz2C(nx,ny,nz), dudzE(nx,ny,nz+1), zE(nx,ny,nz+1))
    allocate(uE(nx,ny,nz+1), d2udz2E(nx,ny,nz+1))
    do k = 1,nz
      zC(:,:,k) = k*dz
    end do 
    zE = 0.d0
    zE(:,:,2:nz+1) = zC
    zC = zC - dz/2.d0
   
    !uC = zC*(2.d0 - zC)
    uC = sin(2.d0*pi*zC)
    uE = sin(2.d0*pi*zE)
    call Pade6opsZ%init(gpC, sp_gpC, dz)
    call Pade6opsZ%ddz_C2E(uC,dudzE,-1,-1)
    
    !print*, maxval(abs(dudzE(1,1,:) - (2.d0 - 2.d0*zE(1,1,:))))
    
    call Pade6opsZ%ddz_E2C(dudzE,d2udz2C,0,0)
    !print*, maxval(abs(d2udz2C + 2.d0))

    call cd06ViscOp%init( nz, dz, .false., .false., .false., .false.) 
    call Pade6opsZ%d2dz2_C2C(uC, d2udz2C,-1,-1)
    call Pade6opsZ%d2dz2_E2E(uE, d2udz2E,-1,-1)
    print*, "**************************************************"
    print*, d2udz2C(1,1,:) - (-((2.d0*pi)**2)*sin(2.d0*pi*zC(1,1,:)))
    print*, d2udz2E(1,1,:) - (-((2.d0*pi)**2)*sin(2.d0*pi*zE(1,1,:)))


end program
