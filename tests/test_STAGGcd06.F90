program testSTAGGcd06
    use kind_parameters, only: rkind
    use cd06staggstuff, only: cd06stagg
    use cd06stuff, only: cd06
    use constants, only: one, pi, two, imi
    use timer, only: tic, toc 

    real(rkind), dimension(:,:,:), allocatable :: zE, zC
    complex(rkind), dimension(:,:,:), allocatable :: fE, fC, dfEt, dfE, dfCt, dfC
    integer :: nx = 1, ny = 1, nz = 36
    real(rkind) :: omega = 2._rkind, dz
    logical :: isTopEven, isBotEven
    type(cd06stagg), allocatable :: der
    type(cd06     ), allocatable :: der2
    integer :: i, j, k

    isTopEven = .true.
    isBotEven = .true. 

    dz = one/real(nz,rkind)
    allocate(zE(nx, ny, nz + 1), zC(nx, ny, nz), fE(nx, ny, nz+1), fC(nx, ny, nz))
    allocate(dfE(nx,ny,nz+1), dfC(nx,ny,nz))  
    allocate(dfEt(nx,ny,nz+1), dfCt(nx,ny,nz))  

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
    call der%init(nz, dz, isTopEven, isBotEven,.false.,.true.)
    allocate(der2)
    i =  der2%init(nz,dz, .false., 0, 0)



    call der%ddz_C2C(fC,dfC,nx,ny)

    dfEt = -omega*two*pi*sin(two*pi*omega*zE) + imi*(-omega*two*pi*sin(two*pi*omega*zE))
    dfCt = -omega*two*pi*sin(two*pi*omega*zC) + imi*(-omega*two*pi*sin(two*pi*omega*zC))

    !print*, dfC
    print*, maxval(abs(real(dfCt,rkind) - real(dfC,rkind)))

    print*, "============================="
    call der2%dd3(real(fC,rkind), zC, nx, ny)
    !print*, zC
    print*, maxval(abs(zC - real(dfCt,rkind)))

    call der%destroy()
    deallocate(der)
    deallocate(zE, zC, fE, fC, dfEt, dfE, dfCt, dfC)

end program 
