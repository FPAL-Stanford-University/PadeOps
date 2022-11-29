program test_cd06

    use kind_parameters, only: rkind
    use constants,       only: two,pi
!    use timer,           only: tic,toc
    use cd06stuff,       only: cd06
    implicit none

    integer:: nx = 256, ny=256, nz=256
    logical, parameter :: periodic = .FALSE.

    type( cd06 ) :: xcd06, ycd06, zcd06
    real(rkind), dimension(:,:,:), allocatable :: x,f,df,df_exact
    real(rkind) :: dx, dy, dz
    real(rkind), parameter :: omega = 1._rkind

    integer :: i,ierr

    allocate( x(nx,ny,nz) )
    allocate( f(nx,ny,nz) )
    allocate( df(nx,ny,nz) )
    allocate( df_exact(nx,ny,nz) )

    dx = two*pi/real(nx,rkind)
    dy = dx
    dz = dx
    do i=1,nx
        x(i,:,:) = real(i-1,rkind)*two*pi/real(nx,rkind)
        f(i,:,:) = sin(omega * x(i,:,:))
        df_exact(i,:,:) = omega * cos( omega * x(i,:,:))
    end do

    ierr = xcd06%init( nx, dx, periodic, 0, 0)
    ierr = ycd06%init( ny, dy, periodic, 0, 0)
    ierr = zcd06%init( nz, dz, periodic, 0, 0)

!    call tic() 
    call xcd06 % dd1(f,df,ny,nz)
!    call toc ("Time to get the x derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - df_exact))

    df = 0._rkind
!    call tic() 
    call ycd06 % dd2(f,df,nx,nz)
!    call toc ("Time to get the y derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df ))
    
    df = 0._rkind
!    call tic() 
    call zcd06 % dd3(f,df,nx,ny)
!    call toc ("Time to get the z derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df ))
    
    
    deallocate( x )
    deallocate( f )
    deallocate( df )
    deallocate( df_exact )

end program
