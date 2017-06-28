program test_pressureCalculation
    use kind_parameters, only: rkind
    use decomp_2d
    use staggOpsMod, only: staggOps
    use spectralMod, only: spectral
    use constants, only: pi
    use poissonMod, only: poisson
    use basic_io, only: write_3d_binary 
    use cd06staggstuff, only: cd06stagg
    use PadePoissonMod, only: padepoisson
    use decomp_2d_io, only: decomp_2d_write_one
    use reductions, only: p_maxval
    use exits, only: message
    use PadeDerOps, only: Pade6stagg

    implicit none

    type(spectral) :: spectC, spectE
    type(decomp_info) :: gpC, gpE
    integer :: nx = 64, ny = 64, nz = 64
    real(rkind), dimension(:,:,:), allocatable :: x, y, zC, zE, urhs, vrhs, wrhs, p, pexact
    complex(rkind), dimension(:,:,:), allocatable :: urhshat, vrhshat, wrhshat
    real(rkind), dimension(:,:,:), allocatable :: tmpCy, tmpCz, tmpEy, tmpEz
    integer :: i, j, k, ierr
    type(staggOps) :: ops
    type(poisson) :: poiss2ndOrder
    type(padepoisson) :: poiss6thOrder
    type(cd06stagg) :: derZ
    type(Pade6stagg) :: Pade6opsZ

    integer :: scheme = 1
    real(rkind) :: dx, dy, dz, maxE

    call MPI_Init(ierr)
    
    call decomp_2d_init(nx, ny, nz, 0, 0)

    call get_decomp_info(gpC)
    call decomp_info_init(nx, ny, nz+1, gpE)

    allocate(x(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(y(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(zC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(urhs(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(vrhs(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(wrhs(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
    allocate(zE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
    allocate(p(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))

    allocate(tmpCy(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3)))
    allocate(tmpCz(gpC%zsz(1),gpC%zsz(2),gpC%zsz(3)))
    allocate(tmpEy(gpE%ysz(1),gpE%ysz(2),gpE%ysz(3)))
    allocate(tmpEz(gpE%zsz(1),gpE%zsz(2),gpE%zsz(3)))
    
    dx = 2.d0*pi/real(nx,rkind) 
    dy = 2.d0*pi/real(ny,rkind) 
    dz = 2.d0*pi/real(nz,rkind) 

    vrhs = 0.d0; wrhs = 0.d0
    do k = gpC%xst(3),gpC%xen(3) 
        do j = gpC%xst(2),gpC%xen(2)
            do i = gpC%xst(1),gpC%xen(1)
                x (i-gpC%xst(1)+1,j-gpC%xst(2)+1,k) = (i-1)*dx
                y (i-gpC%xst(1)+1,j-gpC%xst(2)+1,k) = (j-1)*dy
                zC(i-gpC%xst(1)+1,j-gpC%xst(2)+1,k) = (k-1)*dz + dz/2.d0
            end do 
        end do 
    end do 

    call transpose_x_to_y(zC,tmpCy,gpC)
    call transpose_y_to_z(tmpCy,tmpCz,gpC)
    tmpEz = 0.d0
    tmpEz(:,:,2:nz) = 0.5d0*(tmpCz(:,:,1:nz-1) + tmpCz(:,:,2:nz))
    tmpEz(:,:,nz+1) = tmpEz(:,:,nz) + dz
    call transpose_z_to_y(tmpEz,tmpEy,gpE)
    call transpose_y_to_x(tmpEy,zE,gpE)
    !zC = 0.5d0*(zE(:,:,2:nz+1) + zE(:,:,1:nz))
    

    urhs = sin(2.d0*x)/2.d0
    wrhs = sin(2.d0*zE)/2.d0
    vrhs = 0.d0
    
    pexact = -0.25d0*(cos(2.d0*x) + cos(2.d0*zC))

    call derZ%init(nz,dz,.false.,.false.,.false.,.false.)
    call spectC%init("x", nx, ny, nz  , dx, dy, dz, "four","2/3rd", 2 , .false.)
    call spectE%init("x", nx, ny, nz+1, dx, dy, dz, "four","2/3rd", 2 , .false.)

    call Ops%init(gpC,gpE,0,dx,dy,dz,spectC%spectdecomp,spectE%spectdecomp, .false., .false.)
    call spectC%alloc_r2c_out(urhshat)
    call spectC%alloc_r2c_out(vrhshat)
    call spectE%alloc_r2c_out(wrhshat)

    call spectC%fft(urhs,urhshat)
    call spectC%fft(vrhs,vrhshat)
    call spectE%fft(wrhs,wrhshat)

    call Pade6opsZ%init(gpC, spectC%spectdecomp, gpE, spectE%spectdecomp, dz, scheme)
    call poiss2ndOrder%init(spectC,.false.,dx,dy,dz,Ops,spectE,.true.,gpC) 
    call poiss6thOrder%init(dx, dy, dz, spectC, spectE,.true.,2.d0*pi,.true.,gpC, Pade6opsZ,.false.)

    call poiss2ndOrder%getPressure(urhshat,vrhshat,wrhshat,p)
    maxE = p_maxval(abs(p - pexact))
    call message(1,"max error (2nd Order):", maxE)

    call poiss6thOrder%getPressure(urhshat,vrhshat,wrhshat,p)
    maxE = p_maxval(abs(p - pexact))
    call message(1,"max error (6th Order):", maxE) 

    !call write_3d_binary(p,"TaylorGreenPressure.dat")
    call decomp_2d_write_one(1,p,"TaylorGreenPressure.dat",gpC)
    call MPI_Finalize(ierr)

end program
