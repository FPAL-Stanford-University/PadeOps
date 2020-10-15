program test_interpolator 
    use kind_parameters
    use decomp_2d 
    use interpolatorMod, only: interpolator 
    use gridtools, only: linspace
    use reductions, only: p_maxval 
    use exits, only: message 
    
    implicit none 

    type(decomp_info) :: gpSource, gpDest 
    type(interpolator) :: interp 

    real(rkind), dimension(:,:,:), allocatable :: fSource, fDest, fDestExact 
    real(rkind), dimension(:), allocatable :: xSource, xDest, ySource, yDest, zSource, zDest 
    integer :: i, j, k, ii, jj, kk, ierr  

    integer, parameter :: nxSource = 97, nySource = 23, nzSource = 49
    integer, parameter :: nxDest = 68, nyDest = 54, nzDest = 32 
    real(rkind) :: err 


    call MPI_Init(ierr)

    ! Create the two grids 
    allocate(xSource(nxSource),ySource(nySource), zSource(nzSource))
    allocate(xDest(nxDest), yDest(nyDest), zDest(nzDest))

    xSource = linspace(-5.d0, 7.d0, nxSource)
    ySource = linspace(-25.d0, 8.d0, nySource)
    zSource = linspace(-5.d0, 3.d0, nzSource)

    xDest = linspace(-2.5d0, 4.d0, nxDest)
    yDest = linspace(-13.d0, 6.d0, nyDest)
    zDest = linspace(-1.0d0, 0.d0, nzDest)
   
   
    ! Create the two gp's 
    call decomp_2d_init(nxSource, nySource, nzSource, 0, 0)
    call get_decomp_info(gpSource)
    call decomp_info_init(nxDest,nyDest,nzDest,gpDest)

    ! Create the interpolator 
    call interp%init(gpSource, gpDest, xSource, ySource, zSource, xDest, yDest, zDest) 

    ! Create the testing data 
    allocate(fSource(gpSource%xsz(1),gpSource%xsz(2),gpSource%xsz(3)))
    allocate(fDest(gpDest%xsz(1),gpDest%xsz(2),gpDest%xsz(3)))
    allocate(fDestExact(gpDest%xsz(1),gpDest%xsz(2),gpDest%xsz(3)))

    kk = 1
    do k = gpSource%xst(3),gpSource%xen(3)
        jj = 1
        do j = gpSource%xst(2),gpSource%xen(2)
            ii = 1
            do i = gpSource%xst(1),gpSource%xen(1)
                fSource(ii,jj,kk) = xSource(i) + ySource(j) + zSource(k)  ! test with a linear function (returns exact interpolation)
                ii = ii + 1
            end do 
            jj = jj + 1
        end do 
        kk = kk + 1
    end do
    
    kk = 1
    do k = gpDest%xst(3),gpDest%xen(3)
        jj = 1
        do j = gpDest%xst(2),gpDest%xen(2)
            ii = 1
            do i = gpDest%xst(1),gpDest%xen(1)
                fDestExact(ii,jj,kk) = xDest(i) + yDest(j) + zDest(k)  
                ii = ii + 1
            end do 
            jj = jj + 1
        end do 
        kk = kk + 1
    end do
   

    ! Interpolate
    call interp%LinInterp3D(fSource,fDest)

    ! Check error
    err = p_maxval(abs(fDest - fDestExact))
    if (err > 1E-12) then 
        call message(0,"TEST FAILED.")        
    else 
        call message(1,"TEST PASSED.")        
    end if 


    ! Finalize 
    call interp%destroy()
    call MPI_Finalize(ierr)
end program 