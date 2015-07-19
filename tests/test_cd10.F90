program test_cd10

    use kind_parameters, only: rkind
    use constants,       only: zero,one,two,pi
    use timer,           only: tic,toc
    use cd10stuff,       only: cd10
    implicit none

    integer:: nx = 8192, ny=128, nz=128
    logical, parameter :: periodic = .TRUE.

    type( cd10 ) :: mycd10
    real(rkind), dimension(:,:,:), allocatable :: x,f,df,df_exact
    real(rkind) :: dx!,t0,t1
    real(rkind), parameter :: omega = 1._rkind

    integer :: y1,yn,z1,zn
    integer :: i,ierr

    allocate( x(nx,ny,nz) )
    allocate( f(nx,ny,nz) )
    allocate( df(nx,ny,nz) )
    allocate( df_exact(nx,ny,nz) )

    dx = two*pi/real(nx,rkind)

    do i=1,nx
        x(i,:,:) = real(i-1,rkind)*two*pi/real(nx,rkind)
        f(i,:,:) = sin(omega * x(i,:,:))
        df_exact(i,:,:) = omega * cos( omega * x(i,:,:))
    end do

    ierr = mycd10%init( nx, dx, periodic, 0, 0)
    print*, ierr

    call tic() 
    y1 = 1;yn = ny; z1 = 1; zn = nz
    call mycd10 % cd10der1(f,df,ny,nz,y1,yn,z1,zn)
    call toc ("Time to get the x derivative:")

    print*, "Maximum error = ", MAXVAL( ABS(df - df_exact))

    deallocate( x )
    deallocate( f )
    deallocate( df )
    deallocate( df_exact )

end program
