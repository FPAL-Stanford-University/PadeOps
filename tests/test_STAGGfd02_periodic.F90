program testSTAGGfd02
    use kind_parameters, only: rkind
    use staggOpsMod, only: staggOps
    use constants, only: pi, two, imi, one, zero
    use timer, only: tic, toc 
    use decomp_2d
    use mpi
        
    real(rkind), dimension(:,:,:), allocatable :: zE, zC
    complex(rkind), dimension(:,:,:), allocatable :: fE, fC, dfEt, dfE, dfCt, dfC, IfC, IfE
    complex(rkind), dimension(:,:,:), allocatable :: d2fEt, d2fCt, d2fC, d2fE
    complex(rkind) :: temp = 1.0d0
    integer :: nx = 1, ny = 1, nz = 128, stagg_scheme = 0 
       ! stagg_scheme = 0 for 2nd order FD
    real(rkind) :: omega = 1._rkind, dx, dy, dz
    logical :: isTopEven, isBotEven, isTopSided, isBotSided, isPeriodic
    type(staggOps), allocatable :: der
    integer :: i, j, k, ierr
    real(rkind) :: zzst = 0._rkind, zzend = 1._rkind 
    type(decomp_info) :: gpC, gpE, sp_gpC, sp_gpE
        
    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, 0, 0)
    
    call get_decomp_info(gpC)
    call decomp_info_init(nx, ny, nz+1, gpE)
    call decomp_info_init(nx, ny, nz+1, sp_gpE)
    call decomp_info_init(nx, ny, nz, sp_gpC)

    isTopEven = .true.
    isBotEven = .true. 
    isTopSided = .false.
    isBotSided = .false.
    isPeriodic = .true.

    allocate(zE(nx, ny, nz + 1), zC(nx, ny, nz), fE(nx, ny, nz+1), fC(nx, ny, nz))
    allocate(IfE(nx, ny, nz+1), IfC(nx, ny, nz))
    allocate(dfE(nx,ny,nz+1), dfC(nx,ny,nz),d2fE(nx,ny,nz+1), d2fC(nx,ny,nz))  
    allocate(dfEt(nx,ny,nz+1), dfCt(nx,ny,nz), d2fEt(nx,ny,nz+1), d2fCt(nx,ny,nz))

    dz = (zzend - zzst)/real(nz,rkind)
    dx = dz !value does not matter for this test
    dy = dz !value does not matter for this test
    do k = 1,nz+1
        do j = 1,ny
            do i = 1,nx
                zE(i,j,k) = zzst + real((k - 1),rkind)*dz
            end do 
        end do 
    end do 

    zC = 0.5_rkind*(zE(:,:,2:nz+1) + zE(:,:,1:nz))

    fE = -1*cos(two*pi*omega*zE)! + 1*imi*cos(two*pi*omega*zE)
    fC = -1*cos(two*pi*omega*zC)! + 1*imi*cos(two*pi*omega*zC)

    dfEt = 1*(omega*two*pi)*sin(two*pi*omega*zE)! + -1*imi*(omega*two*pi*sin(two*pi*omega*zE))
    dfCt = 1*(omega*two*pi)*sin(two*pi*omega*zC)! + -1*imi*(omega*two*pi*sin(two*pi*omega*zC))

    d2fEt = 1*((omega**2)*(two**2)*(pi**2))*cos(two*pi*omega*zE)!+1*imi*(-(omega**2)*(two**2)*(pi**2))*cos(two*pi*omega*zE)
    d2fCt = 1*((omega**2)*(two**2)*(pi**2))*cos(two*pi*omega*zC)!+1*imi*(-(omega**2)*(two**2)*(pi**2))*cos(two*pi*omega*zC)

   !Compute Derivatives 
    allocate(der)
    call der%init(gpC, gpE, stagg_scheme, dx, dy, dz, sp_gpC, sp_gpE, isTopSided, isBotSided, isPeriodic)
    call der%ddz_E2C(fE,dfC)
    call der%ddz_C2E(fC, dfE, isTopSided, isBotSided)
    call der%d2dz2_C2C(fC,d2fC,isTopSided,isBotSided)
    call der%d2dz2_E2E(fE,d2fE,isTopSided,isBotSided)

    print*, "-------------------"
    print*, "Max Error for 1st Deriv: E2C:", maxval(abs(dfC - dfCt))
    print*, "Max Error for 1st Deriv: C2E:", maxval(abs(dfE - dfEt))
    print*, "-------------------"
    print*, "-------------------"
    print*, "Max Error for 2nd Deriv: C2C:", maxval(abs(d2fC - d2fCt))
    print*, "Max Error for 2nd Deriv: E2E:", maxval(abs(d2fE - d2fEt))
    print*, "-------------------"

   ! Compute Interpolations
    call der%InterpZ_Edge2Cell(fE,IfC)
    call der%InterpZ_Cell2Edge(fC,IfE,temp,temp)
  
    print*, "-------------------"
    print*, "Max Error for Interp: E2C:", maxval(abs(IfC - fC))
    print*, "Max Error for Interp: C2E:", maxval(abs(IfE - fE))
    print*, "-------------------"
    print*, "==================================================="

    call der%destroy()
    deallocate(der)
    deallocate(zE, zC, fE, fC, dfEt, dfE, dfCt, dfC, IfC, IfE, d2fC, d2fE, d2fEt, d2fCt)

    call MPI_Finalize(ierr)

end program 
