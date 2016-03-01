program test_interp

    use mpi
    use kind_parameters, only : rkind
    use decomp_2d
    use constants, only: pi, two, one, imi, zero
    use reductions, only: p_maxval
    use timer, only: tic, toc
    use staggOpsMod, only: staggOps

    implicit none
   
    real(rkind), dimension(:,:,:), allocatable :: u, v, w, wC, x, y, z, tmp1, tmp2
    real(rkind), dimension(:,:,:), allocatable :: tmp3, tmp4
    type(decomp_info) :: gpC, gpE
    type(staggOps), allocatable :: Ops
    integer :: nx = 8, ny = 8, nz = 32
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

    allocate( u ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( v ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( w ( gpE%xsz(1), gpE%xsz(2), gpE%xsz(3) ) )

    allocate( wC( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( x ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( y ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( z ( gpC%xsz(1), gpC%xsz(2), gpC%xsz(3) ) )
    allocate( tmp1 ( gpC%ysz(1), gpC%ysz(2), gpC%ysz(3) ) )
    allocate( tmp2 ( gpC%zsz(1), gpC%zsz(2), gpC%zsz(3) ) )
    allocate( tmp3 ( gpE%zsz(1), gpE%zsz(2), gpE%zsz(3) ) )
    allocate( tmp4 ( gpE%ysz(1), gpE%ysz(2), gpE%ysz(3) ) )
   
      
    kk = 1
    do k = gpC%xst(3),gpC%xen(3)
        jj = 1
        do j = gpC%xst(2),gpC%xen(2)
            ii = 1
            do i = gpC%xst(1),gpC%xen(1)
                x(ii,jj,kk) = (i-1)*dx
                y(ii,jj,kk) = (j-1)*dy
                z(ii,jj,kk) = (k-1)*dz + dz/two
                ii = ii + 1
            end do 
            jj = jj + 1
        end do 
        kk = kk + 1
    end do 

    wC = cosh(z)
   
    call transpose_x_to_y(wC,tmp1,gpC)
    call transpose_y_to_z(tmp1,tmp2,gpC) 

    call tic()
    call Ops%InterpZ_Cell2Edge(tmp2,tmp3,cosh(0._rkind),cosh(1._rkind))
    call toc()

    call tic()
    call Ops%InterpZ_Edge2Cell(tmp3,tmp2)
    call toc()

    call transpose_z_to_y(tmp3,tmp4,gpE)
    call transpose_y_to_x(tmp4,w,gpE)

    print*,nrank
    call sleep(nrank)
    print*, w(1,1,:)
    print*, "------------"
    print*, wC(1,1,:)
    print*, "------------"
    print*, "------------"

    
    call Ops%destroy
    deallocate(Ops)

    call MPI_Finalize(ierr)
    
end program 
