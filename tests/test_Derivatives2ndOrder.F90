program test_deriv2ndOrder

    use mpi
    use kind_parameters, only : rkind
    use decomp_2d
    use constants, only: pi, two, one, imi, zero
    use reductions, only: p_maxval
    use timer, only: tic, toc
    use staggOpsMod, only: staggOps

    implicit none
   
    real(rkind), dimension(:,:,:), allocatable :: w, wC, x, y, z, zE, tmp1, tmp2
    real(rkind), dimension(:,:,:), allocatable :: tmp3, tmp4
    type(decomp_info) :: gpC, gpE
    type(staggOps), allocatable :: Ops
    integer :: nx = 8, ny = 8, nz = 16
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k, ii, jj, kk
    real(rkind) :: dx, dy, dz, zbot, fmean
    integer :: dimTransform = 2

    call MPI_Init(ierr)
    
    call decomp_2d_init(nx, ny, nz, prow, pcol)

    call get_decomp_info(gpC)
    
    call decomp_info_init(nx, ny, nz+1, gpE)
 
    dx= two*pi/nx; dy = two*pi/ny; dz = one/nz
    allocate(Ops)
    call Ops%Init(gpC,gpE,0, dx, dy, dz)

    allocate( x ( gpC%zsz(1), gpC%zsz(2), gpC%zsz(3) ) )
    allocate( y ( gpC%zsz(1), gpC%zsz(2), gpC%zsz(3) ) )
    allocate( z ( gpC%zsz(1), gpC%zsz(2), gpC%zsz(3) ) )
    allocate( zE( gpE%zsz(1), gpE%zsz(2), gpE%zsz(3) ) )
    allocate( w ( gpE%zsz(1), gpE%zsz(2), gpE%zsz(3) ) )
    allocate( wC ( gpC%zsz(1), gpC%zsz(2), gpC%zsz(3) ) )


    allocate( tmp1 ( gpC%zsz(1), gpC%zsz(2), gpC%zsz(3) ) )
    allocate( tmp2 ( gpC%zsz(1), gpC%zsz(2), gpC%zsz(3) ) )
    allocate( tmp3 ( gpE%zsz(1), gpE%zsz(2), gpE%zsz(3) ) )
    allocate( tmp4 ( gpE%zsz(1), gpE%zsz(2), gpE%zsz(3) ) )
   
      
    kk = 1
    do k = gpC%zst(3),gpC%zen(3)
        jj = 1
        do j = gpC%zst(2),gpC%zen(2)
            ii = 1
            do i = gpC%zst(1),gpC%zen(1)
                x(ii,jj,kk) = (i-1)*dx
                y(ii,jj,kk) = (j-1)*dy
                z(ii,jj,kk) = (k-1)*dz + dz/two
                zE(ii,jj,kk) = (k-1)*dz
                ii = ii + 1
            end do 
            jj = jj + 1
        end do 
        kk = kk + 1
    end do 

    zE(:,:,gpE%zen(3)) = zE(:,:,gpE%zen(3)-1) + dz 

    wC = cos(two*pi*z)
    w  = cos(two*pi*zE)

    !! First Derivative Check 
    ! Edge to center
    call Ops%ddz_E2C(w,tmp1)
    tmp2 = -two*pi*sin(two*pi*z)
    print*, "First Deivative E->C"
    print*, maxval(abs(tmp1 - tmp2))

    ! Center to Center 
    call Ops%ddz_C2C(wC,tmp1,.true.,.true.)
    tmp2 = -two*pi*sin(two*pi*z)
    print*, "First Deivative C->C"
    print*, maxval(abs(tmp1 - tmp2))


    ! Center to Edge
    call Ops%ddz_C2E(wC,tmp3,.true.,.true.)
    tmp4 = -two*pi*sin(two*pi*zE)
    print*, "First Deivative C->E"
    print*, maxval(abs(tmp4 - tmp3))


    !! Second Derivative Check 
    ! Edge to Edge
    call Ops%d2dz2_E2E(w,tmp3,.true.,.true.)
    tmp4 = -4._rkind*(pi**2)*cos(two*pi*zE)
    print*, "Second Derivative E->E"
    print*, maxval(abs(tmp4 - tmp3))
  
    call Ops%d2dz2_C2C(wC,tmp1,.true.,.true.)
    tmp2 = -4._rkind*(pi**2)*cos(two*pi*z)
    print*, "Second Derivative C->C"
    print*, maxval(abs(tmp2 - tmp1))
  
    call Ops%destroy
    deallocate(Ops)

    call MPI_Finalize(ierr)
    
end program 
