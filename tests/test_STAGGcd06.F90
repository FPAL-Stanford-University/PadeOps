program testSTAGGcd06
    use kind_parameters, only: rkind
    use cd06staggstuff, only: cd06stagg
    use cd06stuff, only: cd06
    use constants, only: pi, two, imi
    use timer, only: tic, toc 
    use staggOpsMod, only: staggops
    use decomp_2d
    use mpi
        
    real(rkind), dimension(:,:,:), allocatable :: zE, zC
    real(rkind), dimension(:,:,:), allocatable :: fE, fC, dfEt, dfE, dfCt, dfC, IfC, IfE
    integer :: nx = 1, ny = 1, nz = 64
    real(rkind) :: omega = 1._rkind, dz
    logical :: isTopEven, isBotEven
    type(cd06stagg), allocatable :: der
    !type(cd06     ), allocatable :: der2
    integer :: i, j, k, ierr
    real(rkind) :: zzst = 0._rkind, zzend = 1._rkind 
    !type(staggops) :: ops
    type(decomp_info) :: gp, gpE, sp_gp, sp_gpE


    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, 0, 0)
    
    call get_decomp_info(gp)
    call decomp_info_init(nx, ny, nz+1, gpE)
    call decomp_info_init(nx, ny, nz+1, sp_gpE)
    call decomp_info_init(nx, ny, nz, sp_gp)

    isTopEven = .true.
    isBotEven = .true. 

    allocate(zE(nx, ny, nz + 1), zC(nx, ny, nz), fE(nx, ny, nz+1), fC(nx, ny, nz))
    allocate(IfE(nx, ny, nz+1), IfC(nx, ny, nz))
    allocate(dfE(nx,ny,nz+1), dfC(nx,ny,nz))  
    allocate(dfEt(nx,ny,nz+1), dfCt(nx,ny,nz))  

    dz = (zzend - zzst)/real(nz,rkind)
    do k = 1,nz+1
        do j = 1,ny
            do i = 1,nx
                zE(i,j,k) = zzst + real((k - 1),rkind)*dz
            end do 
        end do 
    end do 

    zC = 0.5_rkind*(zE(:,:,2:nz+1) + zE(:,:,1:nz))

    print*, zC
    fE = 1*cos(two*pi*omega*zE) !+ 1*imi*cos(two*pi*omega*zE)
    fC = 1*cos(two*pi*omega*zC) !+ 1*imi*cos(two*pi*omega*zC)

    dfEt = 1*(-omega*two*pi)*sin(two*pi*omega*zE) !+ 1*imi*(-omega*two*pi*sin(two*pi*omega*zE))
    dfCt = 1*(-omega*two*pi)*sin(two*pi*omega*zC) !+ 1*imi*(-omega*two*pi*sin(two*pi*omega*zC))

   !Compute Derivatives 
    allocate(der)
    !call der%init(nz, dz, isTopEven, isBotEven,.true.,.true.)
    call der%init(nz, dz)
    call der%ddz_E2C(fE,dfC,nx,ny)
    call der%ddz_C2E(fC,dfE,nx,ny)

!    print*,'computed derivative at cells'
!    print*, dfC(1,1,:)
!    print*, "-------------------"
!    print*, "computed derivative at edges"
!    print*, dfE(1,1,:)

    print*, "-------------------"
    print*, "Max Error for 1st Deriv: E2C:", maxval(abs(dfC - dfCt))
    print*, "Max Error for 1st Deriv: C2E:", maxval(abs(dfE - dfEt))
    print*, "-------------------"

   ! Compute Interpolations
    call der%InterpZ_E2C(fE,IfC,nx,ny)
    call der%InterpZ_C2E(fC,IfE,nx,ny)
  
!    print*,'interpolated value at cells'
!    print*, IfC(1,1,:)
!    print*, "-------------------"
!    print*,'interpolated value at edges'
!    print*, IfE(1,1,:)
    print*, "-------------------"
    print*, "Max Error for Interp: E2C:", maxval(abs(IfC - fC))
    print*, "Max Error for Interp: C2E:", maxval(abs(IfE - fE))
    print*, "-------------------"
  !print*, "==================================================="
    !print*, dfCt(1,1,:)
    !print*, "==================================================="
    !print*, dfCt(1,1,:) - dfC(1,1,:)
    !print*, "==================================================="
    !print*, "Second order" 
    !call ops%ddz_C2C(fC, dfC, .true., .true.)

    !print*, dfC(1,1,:)
    !print*, "==================================================="
    !print*, dfCt(1,1,:)
    !print*, "==================================================="
    !print*, dfCt(1,1,:) - dfC(1,1,:)
    !!print*, dfC
    !!print*, maxval(abs(real(dfCt,rkind) - real(dfC,rkind)))

    !print*, "============================="
    !call der2%dd3(real(fC,rkind), zC, nx, ny)
    !print*, zC(1,1,1:5)
    !print*, zC

    !print*, dfC - dfCt

    call der%destroy()
    deallocate(der)
    deallocate(zE, zC, fE, fC, dfEt, dfE, dfCt, dfC)

end program 
