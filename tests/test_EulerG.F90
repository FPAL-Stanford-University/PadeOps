program test_EulerG
    use kind_parameters, only: rkind 
    use decomp_2d  
    use EulerG_mod, only: EulerG 
    use mpi 
    use gridTools, only: linspace 
    use reductions, only: p_maxval 
    use exits, only: message 

    integer :: nx = 24, ny = 16, nz = 8  ! Base
    integer, parameter :: levels = 4
    type(decomp_info) :: gp, gpFinal
    real(rkind), dimension(2), parameter :: DomX = [0.d0, 1.d0], DomY = [2.d0, 5.d0], DomZ = [-3.d0, 5.d0]
    type(EulerG) :: Egrid
    real(rkind), dimension(:), allocatable :: x, y, z
    real(rkind), dimension(:,:,:), pointer :: uBase, vBase, wBase, uF, vF, wF 
    real(rkind), dimension(:,:,:), allocatable :: uFinal, vFinal, wFinal 
    integer :: ierr, i, j, k, ii, jj, kk 
    real(rkind) :: errU, errV, errW 

    call MPI_Init(ierr)

    call decomp_2d_init(nx+1, ny+1, nz+1, 0, 0)
    call get_decomp_info(gp)
   
    call  Egrid%init( gp, DomX, DomY, DomZ, levels)
    

    allocate(x(nx+1), y(ny+1), z(nz+1))

    x = linspace(DomX(1), DomX(2), nx+1)
    y = linspace(DomY(1), DomY(2), ny+1)
    z = linspace(DomZ(1), DomZ(2), nz+1)

    call Egrid%allocateFieldForLevel(uFinal,levels)
    call Egrid%allocateFieldForLevel(vFinal,levels)
    call Egrid%allocateFieldForLevel(wFinal,levels)

    call Egrid%getPointersToFields(uBase, vBase, wBase, 1)
    call Egrid%getPointersToFields(uF, vF, wF, levels)

    kk = 1
    do k = gp%xst(3),gp%xen(3)
        jj = 1
        do j = gp%xst(2),gp%xen(2)
            ii = 1
            do i = gp%xst(1),gp%xen(1)
                uBase(ii,jj,kk) = x(i) + 3*y(j) + 7*z(k)  ! test with a linear function (returns exact interpolation)
                vBase(ii,jj,kk) = 2*x(i) + y(j) + z(k)    ! test with a linear function (returns exact interpolation)
                wBase(ii,jj,kk) = x(i) + 2*y(j) + 5*z(k)  ! test with a linear function (returns exact interpolation)
                ii = ii + 1
            end do 
            jj = jj + 1
        end do 
        kk = kk + 1
    end do


    deallocate(x, y, z)
    allocate(x(2**(levels-1)*nx + 1))
    allocate(y(2**(levels-1)*ny + 1))
    allocate(z(2**(levels-1)*nz + 1))
    
    x = linspace(DomX(1), DomX(2),2**(levels-1)*nx + 1)
    y = linspace(DomY(1), DomY(2),2**(levels-1)*ny + 1)
    z = linspace(DomZ(1), DomZ(2),2**(levels-1)*nz + 1)

    call decomp_info_init(2**(levels-1)*nx + 1,2**(levels-1)*ny + 1,2**(levels-1)*nz + 1,gpFinal)
    kk = 1
    do k = gpFinal%xst(3),gpFinal%xen(3)
        jj = 1
        do j = gpFinal%xst(2),gpFinal%xen(2)
            ii = 1
            do i = gpFinal%xst(1),gpFinal%xen(1)
                uFinal(ii,jj,kk) = x(i) + 3*y(j) + 7*z(k)  ! test with a linear function (returns exact interpolation)
                vFinal(ii,jj,kk) = 2*x(i) + y(j) + z(k)    ! test with a linear function (returns exact interpolation)
                wFinal(ii,jj,kk) = x(i) + 2*y(j) + 5*z(k)  ! test with a linear function (returns exact interpolation)
                ii = ii + 1
            end do 
            jj = jj + 1
        end do 
        kk = kk + 1
    end do

    call Egrid%agglomerate()


    errU = p_maxval(abs(uFinal - uF))
    errV = p_maxval(abs(vFinal - vF))
    errW = p_maxval(abs(wFinal - wF))
    
    if ((errU < 1E-10) .and. (errV < 1E-10) .and. (errW < 1E-10)) then 
        call message(0, "Test passed.")
    else
        print*, nrank, maxval(abs(uFinal - uF)), maxval(abs(vFinal - vF)), maxval(abs(wFinal - wF))
        call message(0, "Test Failed.")
    end if 
    
    call decomp_2d_finalize
    call MPI_Finalize(ierr)


end program 