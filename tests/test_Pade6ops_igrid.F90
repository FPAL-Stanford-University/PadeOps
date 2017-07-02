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
    use PadePoissonMod, only: Padepoisson 
    implicit none

    integer, parameter :: nx = 16, ny = 16, nz = 64
    real(rkind), dimension(nx,ny,nz) :: x, y, z, u, v, wC, rbuffC1, rbuffC2, urhs
    real(rkind), dimension(nx,ny,nz+1) :: w, xE, zE, uE, rbuffE1, rbuffE2, wrhs

    complex(rkind), dimension(nx/2+1,ny,nz) :: uhat, wChat, vhat, cbuffC1, cbuffC2
    complex(rkind), dimension(nx/2+1,ny,nz+1) :: what, cbuffE1, cbuffE2
    type(padepoisson) :: padepoiss
    type(decomp_info) :: gpC, gpE
    type(decomp_info), pointer :: sp_gpC, sp_gpE
    type(spectral), target :: spectC, spectE
    integer :: ierr, k, i, j 
    real(rkind) :: dx, dy, dz
    type(Pade6stagg) :: Pade6opsZ
    real(rkind), parameter :: Re = 10.d0, dt = 1.d-4
    real(rkind) :: t 
    integer, parameter :: scheme = 1

    !! STEP 0: initialize mesh and fields
    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(gpC)
    call decomp_info_init(nx,ny,nz+1,gpE)

    if (nproc>1) call gracefulExit("Meant to be a serial execution",13)
    
    dx = two*pi/real(nx,rkind)
    dy = two*pi/real(ny,rkind)
    dz = two*pi/real(nz,rkind)

    call spectC%init("x", nx, ny, nz  , dx, dy, dz, "four", "2/3rd", 2, .false.)
    call spectE%init("x", nx, ny, nz+1, dx, dy, dz, "four", "2/3rd", 2, .false.)
    sp_gpC => spectC%spectdecomp
    sp_gpE => spectE%spectdecomp
    call Pade6opsZ%init(gpC, sp_gpC, gpE, sp_gpE, dz, scheme)
    
    do k=1,nz
        do j=1,ny
            do i=1,nz
                x(i,j,k) = real(  i , rkind ) * dx
                y(i,j,k) = real(  j , rkind ) * dy
                z(i,j,k) = real(  k , rkind ) * dz + dz/two
            end do
        end do
    end do

    x = x - dx; y = y - dy; z = z - dz 
    zE(:,:,1:nz) = z - dz/2.d0; zE(:,:,nz+1) = zE(:,:,nz) + dz;
    xE(:,:,1:nz) = x; xE(:,:,nz+1) = xE(:,:,nz)

    u  = sin(x)*cos(z)
    wC = -cos(x)*sin(z) 
    v  = 0.d0
    call Pade6opsZ%interpz_C2E(wC,w,-1,-1)
    call padepoiss%init(dx,dy,dz, spectC, spectE, .true., two*pi, .true., gpC, Pade6opsZ, .false.) 

    call spectC%fft(u,uhat)
    call spectC%fft(v,vhat)
    call spectE%fft(w,what)
    call padepoiss%PressureProjection(uhat,vhat,what)
    call spectC%ifft(uhat,u)
    call spectC%ifft(vhat,v)
    call spectE%ifft(what,w)
    call Pade6opsZ%interpz_C2E(u,uE, 1, 1)
    t = 0.d0

    print*, "Advection Term:"
    !! STEP 1: Compute the advection flux for u
    rbuffC1 = 0.5d0*u*u
    call spectC%fft(rbuffC1,cbuffC1)
    call spectC%mtimes_ik1_ip(cbuffC1)
    call spectC%ifft(cbuffC1, urhs)
    rbuffE1 = 0.5d0*uE*w
    call Pade6opsZ%ddz_E2C(rbuffE1,rbuffC1,-1,-1)
    urhs = urhs + rbuffC1
    call spectC%mtimes_ik1_oop(uhat,cbuffC1)
    call spectC%ifft(cbuffC1,rbuffC1)
    urhs = urhs + 0.5d0*(rbuffC1*u)
    call Pade6opsZ%ddz_C2E(u,rbuffE1, 1, 1)
    call Pade6opsZ%interpz_E2C(w*rbuffE1, rbuffC1, 1, 1)
    urhs = urhs + 0.5d0*rbuffC1

    print*, "urhs error:", maxval(abs(urhs - (0.5d0*sin(2.d0*x))))


    !! STEP 2: Compute the advection flux for w
    rbuffE1 = 0.5d0*uE*w
    call spectE%fft(rbuffE1,cbuffE1)
    call spectE%mtimes_ik1_ip(cbuffE1)
    call spectE%ifft(cbuffE1,wrhs)
    rbuffC1 = 0.5d0*wC*wC
    call Pade6opsZ%ddz_C2E(rbuffC1,rbuffE1, 1, 1)
    wrhs = wrhs + rbuffE1
    call spectE%mtimes_ik1_oop(what,cbuffE1)
    call spectE%ifft(cbuffE1,rbuffE1)
    wrhs = wrhs + 0.5d0*uE*rbuffE1
    call Pade6opsZ%ddz_C2E(wC,rbuffE1,-1, -1)
    wrhs = wrhs + 0.5d0*w*rbuffE1

    print*, "wrhs error:", maxval(abs(wrhs - (0.5d0*sin(2.d0*zE))))
    
    urhs = -urhs
    wrhs = -wrhs

    print*, "Diffusion Term:"
    !! STEP 3: Compute diffusion term for u
    call Pade6opsZ%d2dz2_C2C(u,rbuffC1,1,1)
    urhs = urhs + (1.d0/Re)*rbuffC1
    cbuffC1 = -spectC%kabs_sq*uhat
    call spectC%ifft(cbuffC1,rbuffC1)
    urhs = urhs + (1.d0/Re)*rbuffC1
    print*, "urhs error:", maxval(abs(urhs - (- sin(2.d0*x)/2.d0 - (cos(z)*sin(x))/5.d0)))


    !! STEP 4: Compute diffusion term for w
    call Pade6opsZ%d2dz2_E2E(w,rbuffE1,-1,-1)
    wrhs = wrhs + (1.d0/Re)*rbuffE1
    cbuffE1 = -spectE%kabs_sq*what
    call spectE%ifft(cbuffE1,rbuffE1)
    wrhs = wrhs + (1.d0/Re)*rbuffE1

    print*, "wrhs error:", maxval(abs(wrhs - ((cos(xE)*sin(zE))/5.d0 - sin(2.d0*zE)/2.d0)))


    print*, "Post-time stepping and projection:"
    ! STEP 5: Do time advance
    u = u + dt*urhs
    w = w + dt*wrhs
    t = t + dt
    call spectC%fft(u,uhat)
    call spectE%fft(w,what)

    ! STEP 6: Do the poisson projection
    call padepoiss%PressureProjection(uhat,vhat,what)
    call spectC%ifft(uhat,u)
    call spectE%ifft(what,w)

    print*, "Error in u:", maxval(abs(u - ( sin(x )*cos(z )*exp(-2.d0*t/Re))))
    print*, "Error in w:", maxval(abs(w - (-cos(xE)*sin(zE)*exp(-2.d0*t/Re))))


end program
