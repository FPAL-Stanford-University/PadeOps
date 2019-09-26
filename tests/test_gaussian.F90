program test_gaussian

    use kind_parameters, only: rkind
    use constants,       only: two,pi
    use timer,           only: tic,toc
    use gaussianstuff,   only: gaussian
    implicit none

    integer :: nx = 32, ny=32, nz=32

    logical, parameter :: periodic = .FALSE.

    type( gaussian ) :: xgaussian, ygaussian, zgaussian
    real(rkind), dimension(:,:,:), allocatable :: x,y,z,f,df,df_x,df_y,df_z
    real(rkind) :: dx, dy, dz
    real(rkind) :: omega = 12._rkind
    real(rkind) :: TF_x, k_norm_x
    real(rkind) :: TF_y, k_norm_y
    real(rkind) :: TF_z, k_norm_z

    integer :: iounit = 17

    integer :: i,j,k,ierr

    allocate( x(nx,ny,nz) )
    allocate( y(nx,ny,nz) )
    allocate( z(nx,ny,nz) )
    allocate( f(nx,ny,nz) )
    allocate( df(nx,ny,nz) )
    allocate( df_x(nx,ny,nz) )
    allocate( df_y(nx,ny,nz) )
    allocate( df_z(nx,ny,nz) )

    if (periodic) then
        dx = two*pi/real(nx,rkind)
        dy = two*pi/real(ny,rkind)
        dz = two*pi/real(nz,rkind)
    else
        dx = two*pi/real(nx-1,rkind)
        dy = two*pi/real(ny-1,rkind)
        dz = two*pi/real(nz-1,rkind)
    end if

    k_norm_x = omega * dx
    TF_x = GetTransferFunction(k_norm_x)
    k_norm_y = omega * dy
    TF_y = GetTransferFunction(k_norm_y)
    k_norm_z = omega * dz
    TF_z = GetTransferFunction(k_norm_z)

    print*, "Wavenumber = ", omega
    print*, "Normalized wavenumber = ", k_norm_x
    print*, "Transfer Function value = ", TF_x

    do k=1,nz
        do j=1,ny
            do i=1,nx
                x(i,j,k) = real(i-1,rkind)*dx
                y(i,j,k) = real(j-1,rkind)*dy
                z(i,j,k) = real(k-1,rkind)*dz
                
                f(i,j,k) = cos(omega * x(i,j,k)) * cos(omega * y(i,j,k)) * cos(omega * z(i,j,k))
                
                df_x(i,j,k) = TF_x * f(i,j,k)
                df_y(i,j,k) = TF_y * f(i,j,k)
                df_z(i,j,k) = TF_z * f(i,j,k)
            end do
        end do
    end do
    print*, "Created initial data"

    ierr = xgaussian%init( nx, periodic, .false.)
    ierr = ygaussian%init( ny, periodic, .false.)
    ierr = zgaussian%init( nz, periodic, .false.)

    call tic() 
    call xgaussian % filter1(f,df,ny,nz)
    call toc ("Time to get the x filter:")
    print*, "Maximum error = ", MAXVAL( ABS(df - df_x))

    call tic() 
    call ygaussian % filter2(f,df,nx,nz)
    call toc ("Time to get the y filter:")
    print*, "Maximum error = ", MAXVAL( ABS(df - df_y))

    OPEN(UNIT=iounit, FILE="gaussian.txt", FORM='FORMATTED')
    do i=1,ny
        WRITE(iounit,'(3ES24.16)') y(1,i,1), f(1,i,1), df(1,i,1)
    end do
    CLOSE(iounit)

    call tic() 
    call zgaussian % filter3(f,df,nx,ny)
    call toc ("Time to get the z filter:")
    print*, "Maximum error = ", MAXVAL( ABS(df - df_z))

    deallocate( x )
    deallocate( y )
    deallocate( z )
    deallocate( f )
    deallocate( df )
    deallocate( df_x )
    deallocate( df_y )
    deallocate( df_z )

contains

    function GetTransferFunction(k) result(T)
        use gaussianstuff, only: agf, bgf, cgf, dgf, egf
        real(rkind), intent(in) :: k
        real(rkind) :: T

        T = (agf + two* bgf*COS(k) + two*cgf*COS(two*k) +two*dgf*COS(3._rkind*k) + two*egf*COS(4._rkind*k) )

    end function

end program
