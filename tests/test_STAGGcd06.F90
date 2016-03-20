program testSTAGGcd06
    use kind_parameters, only: rkind
    use cd06staggstuff, only: cd06stagg
    use constants, only: one, pi, two, imi
    use timer, only: tic, toc 

    real(rkind), dimension(:,:,:), allocatable :: zE, zC
    complex(rkind), dimension(:,:,:), allocatable :: fE, fC, dfE, dfC
    integer :: nx = 64, ny = 64, nz = 512
    real(rkind) :: omega = 2._rkind, dz
    logical :: isTopEven, isBotEven
    type(cd06stagg), allocatable :: der
    integer :: i, j, k

    isTopEven = .true.
    isBotEven = .true. 

    dz = one/real(nz,rkind)
    allocate(zE(nx, ny, nz + 1), zC(nx, ny, nz), fE(nx, ny, nz+1), fC(nx, ny, nz))
    allocate(dfE(nx,ny,nz+1), dfC(nx,ny,nz))  

    do k = 1,nz+1
        do j = 1,ny
            do i = 1,nx
                zE(i,j,k) = real((k - 1),rkind)*dz
            end do 
        end do 
    end do 

    zC = 0.5_rkind*(zE(:,:,2:nz+1) + zE(:,:,1:nz))

    fE = cos(two*pi*omega*zE) + imi*cos(two*pi*omega*zE)
    fC = cos(two*pi*omega*zC) + imi*cos(two*pi*omega*zC)

    allocate(der)
    call der%init(nz, dz, isTopEven, isBotEven)

    call tic()
    call der%ddz_E2E(fE,dfE,nx,ny)
    call toc()



    call der%destroy()
    deallocate(der)

end program 
