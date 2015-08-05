program test_Cd10_Der2
    use kind_parameters, only: rkind
    use cd10stuff, only: cd10
    use constants, only: pi

    type(cd10) :: xcd10
    integer :: nx =32, ny = 1, nz = 1
    real(rkind), dimension(:,:,:), allocatable :: x, f, d2f, d2f_Exact
    real(rkind) :: dx
    integer :: i, ierr

    dx = 2.*pi/(nx)

    allocate(x(nx,ny,nz),f(nx,ny,nz),d2f(nx,ny,nz),d2f_Exact(nx,ny,nz))
    
    do i = 1,nx
        x(i,:,:) = (i - 1)*dx
    end do

    f = cos(x)
    d2f_exact = -cos(x)

    ierr =  xcd10%init ( nx, dx , .false., 0, 0) 

    call xcd10%d2d1(f,d2f,ny,nz)


    print*, maxval(abs(d2f - d2f_Exact))

    call xcd10%dd1(f,d2f,ny,nz)
    d2f_exact = -sin(x)
    print*, maxval(abs(d2f - d2f_Exact))

end program
