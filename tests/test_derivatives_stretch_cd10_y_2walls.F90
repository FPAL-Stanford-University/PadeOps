module subroutines


contains

subroutine write_file_x(fname, nx, x, y, z, f, df, dfdx_exact)
    use kind_parameters, only: rkind

    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in) :: nx
    real(rkind), dimension(:,:,:), intent(in) :: x, y, z, f, df, dfdx_exact
    integer :: i, j, k

    j = 7; k = 9;
    open(10,file=fname,status='unknown')
    do i = 1, nx
      write(10,'(6(e19.12,1x))') x(i,j,k), y(i,j,k), z(i,j,k), f(i,j,k), df(i,j,k), dfdx_exact(i,j,k)
    enddo
    close(10)
end subroutine write_file_x

subroutine write_file_y(fname, ny, x, y, z, f, df, dfdx_exact)
    use kind_parameters, only: rkind

    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in) :: ny
    real(rkind), dimension(:,:,:), intent(in) :: x, y, z, f, df, dfdx_exact
    integer :: i, j, k

    i = 7; k = 9;
    open(10,file=fname,status='unknown')
    do j = 1, ny
      write(10,'(6(e19.12,1x))') x(i,j,k), y(i,j,k), z(i,j,k), f(i,j,k), df(i,j,k), dfdx_exact(i,j,k)
    enddo
    close(10)

end subroutine write_file_y

end module

program test_derivatives

    use kind_parameters, only: rkind
    use constants,       only: two,pi,half,zero,one
    use timer,           only: tic,toc
    use DerivativesMod,  only: derivatives
    use subroutines
    implicit none

    integer :: nx=16, ny=16, nz=16


    type( derivatives ) :: method1, method2, method3
    real(rkind), dimension(:,:,:), allocatable :: x,y,z,f,df,dfdx_exact,dfdy_exact,dfdz_exact,d2fdx2_exact,d2fdy2_exact,d2fdz2_exact
    real(rkind), dimension(:,:,:), allocatable :: xi, eta, zeta
    real(rkind) :: dx, dy, dz
    real(rkind), parameter :: omega = 3._rkind

    real(rkind) :: beta = 1.054d0, alpha = 0.50d0 , Lx = two*pi, Ly = two*pi, Lz = two*pi, yh = two*pi, ystart = -pi
    real(rkind) :: num, num1, den, BB, BB2
    real(rkind) :: ybyyfocm1, yuniform_loc, yfocus_loc
    real(rkind) :: linf, l2norm, err1, err2, err3, err4, errcen
    integer :: i,j,k,ierr
    logical :: xmetric=.false., ymetric=.false., zmetric=.false.
    integer :: xmetric_flag = 1, ymetric_flag = 2, zmetric_flag = 1 
    real(rkind), allocatable, dimension(:,:) :: metric_params

    namelist /INPUT/ nx, ny, nz, xmetric, ymetric, zmetric
    namelist /METRICS/ xmetric_flag, ymetric_flag, zmetric_flag, metric_params
    open(unit=123, file='input.dat', form='FORMATTED', iostat=ierr)
    read(unit=123, NML=INPUT)
    close(123)

    allocate(metric_params(3,5))    ! 3 :: (x,y,z); 5 :: max no of parameters
    metric_params = zero
    open(unit=15, file='input.dat', form='FORMATTED')
    read(unit=15, NML=METRICS)
    close(15)

    allocate( x(nx,ny,nz) )
    allocate( y(nx,ny,nz) )
    allocate( z(nx,ny,nz) )
    allocate( f(nx,ny,nz) )
    allocate( df(nx,ny,nz) )
    allocate( dfdx_exact(nx,ny,nz) )
    allocate( dfdy_exact(nx,ny,nz) )
    allocate( dfdz_exact(nx,ny,nz) )
    allocate( d2fdx2_exact(nx,ny,nz) )
    allocate( d2fdy2_exact(nx,ny,nz) )
    allocate( d2fdz2_exact(nx,ny,nz) )

    if(xmetric) then
        allocate( xi(nx,ny,nz) )
    else
        allocate(  xi( 0, 0, 0) )
    endif
    if(ymetric) then
        allocate( eta(nx,ny,nz) )
    else
        allocate( eta( 0, 0, 0) )
    endif
    if(zmetric) then
        allocate(zeta(nx,ny,nz) )
    else
        allocate(zeta( 0, 0, 0) )
    endif

    dx = Lx/real(max(nx-1,1),rkind)
    dy = Ly/real(max(ny-1,1),rkind)
    dz = Lz/real(max(nz-1,1),rkind)

    ! concentrate towards the center -- Pletcher, Tannehill, Anderson
    ! (Section 5.6, Transformation 3, pg. 332)
    
    alpha  = metric_params(2,1);  beta   = metric_params(2,2); ystart = metric_params(2,3); yh = metric_params(2,4)
    BB = (beta + 1) / (beta - 1)
    print *, '>>towards walls<<', alpha, beta, ystart, yh

    do k=1,nz
    do j=1,ny
    do i=1,nx
        x(i,j,k) = real(i-1,rkind)*dx
        y(i,j,k) = ystart + real(j-1,rkind)*dy
        z(i,j,k) = real(k-1,rkind)*dz

        ! stretched grid
        if(ymetric) then
            ! first set the uniform coordinate
            eta(i,j,k) = ystart + real(j-1,rkind)*dy

            ! then set the non-uniform coordinate
            yuniform_loc = (eta(i,j,k) - ystart) / yh
            BB2 = BB ** ( (yuniform_loc-alpha) / (1-alpha) )
            num = ((beta+2*alpha)*BB2 - beta + 2*alpha ) * yh
            y(i,j,k) = num/( (2*alpha+1)*(1+BB2) )   + ystart
        endif

        !! sin
        f(i,j,k) = 5.0_rkind + sin(omega * x(i,j,k)) + sin(omega * y(i,j,k)) + sin(omega * z(i,j,k))
        dfdx_exact(i,j,k) = omega * cos( omega * x(i,j,k))
        dfdy_exact(i,j,k) = omega * cos( omega * y(i,j,k)) 
        dfdz_exact(i,j,k) = omega * cos( omega * z(i,j,k)) 
        d2fdx2_exact(i,j,k) = -omega * omega*sin( omega * x(i,j,k))
        d2fdy2_exact(i,j,k) = -omega * omega*sin( omega * y(i,j,k)) 
        d2fdz2_exact(i,j,k) = -omega * omega*sin( omega * z(i,j,k)) 

        if(ymetric .and. i==3 .and. k==19) then
            !print '(10(e19.12,1x))', y(i,j,k), eta(i,j,k), f(i,j,k), dfdx_exact(i,j,k)
        endif
        !!! cos
        !f(i,j,k) = cos(omega * x(i,j,k)) + cos(omega * y(i,j,k)) + cos(omega * z(i,j,k))
        !dfdx_exact(i,j,k) = -omega * sin( omega * x(i,j,k))
        !dfdy_exact(i,j,k) = -omega * sin( omega * y(i,j,k)) 
        !dfdz_exact(i,j,k) = -omega * sin( omega * z(i,j,k)) 
        !d2fdx2_exact(i,j,k) = -omega * omega*cos( omega * x(i,j,k))
        !d2fdy2_exact(i,j,k) = -omega * omega*cos( omega * y(i,j,k)) 
        !d2fdz2_exact(i,j,k) = -omega * omega*cos( omega * z(i,j,k)) 
    end do
    end do 
    end do 

    !!print*, "Created initial data"
    !!do j=1,ny
    !!  write(*,'(i5,1x,2(e19.12,1x))') j, eta(1,j,1), y(1,j,1)
    !!enddo

    call method1%init(          nx,     ny,    nz, &
                                dx,     dy,    dz, &
                        .false., .false., .false., &
                           "cd10", "cd10", "cd10", &
                                x,      y,      z, &
     !!                   .false., .false., .false., &
     !!                   .false.,  .true., .false., &
                        xmetric, ymetric, zmetric, &
                             .false., 'input.dat', &
                               xi,    eta,   zeta  )

    !!call method2%init(          nx,     ny,    nz, &
    !!                            dx,     dy,    dz, &
    !!                    .false., .false., .false., &
    !!                       "cd10", "cd10", "cd10", &
    !!                            x,      y,      z, &
    !!                    .false.,  .true., .false., &
    !!                         .false., "input.dat", &
    !!                           xi,    eta,   zeta  )

    print*, "Initialized all methods"
  
    print*, "FIRST DERIVATIVE TESTS" 
   
    print*, "==========================================="
    print*, "Now trying METHOD 1: CD10"
    print*, "==========================================="
    !! -- 0, 0 -- !!
    call method1 % ddy(f,df,0,0)
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))
    linf = MAXVAL( ABS(df - dfdy_exact))
    l2norm = SQRT(SUM(ABS(df - dfdy_exact)**2)/(nx*ny*nz))
    errcen =  ABS(df(nx/2,ny/2,nz/2) - dfdy_exact(nx/2,ny/2,nz/2))
    open(unit=123, file='err_00.dat', action='write',iostat=ierr,position='append')
    !write(123,'(i5,1x,e19.12,1x)') ny, MAXVAL( ABS(df - dfdy_exact))
    write(123,'(i5,1x,e19.12,1x,e19.12,1x,e19.12,1x)') ny, linf, l2norm, errcen
    close(123)
    call write_file_y('dfdy_00.dat', ny, x, y, z, f, df, dfdy_exact)

    !! -- 0, 1 -- !!
    call method1 % ddy(f,df,0,1)
    print *, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))
    open(unit=123, file='err_01.dat', action='write',iostat=ierr,position='append')
    write(123,'(i5,1x,e19.12,1x)') ny, MAXVAL( ABS(df - dfdy_exact))
    close(123)
    call write_file_y('dfdy_01.dat', ny, x, y, z, f, df, dfdy_exact)

    !! -- 0, -1 -- !!
    call method1 % ddy(f,df,0,-1)
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))
    open(unit=123, file='err_0m1.dat', action='write',iostat=ierr,position='append')
    write(123,'(i5,1x,e19.12,1x)') ny, MAXVAL( ABS(df - dfdy_exact))
    close(123)
    call write_file_y('dfdy_0m1.dat', ny, x, y, z, f, df, dfdy_exact)

    !! -- 1, 0 -- !!
    call method1 % ddy(f,df,1,0)
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))
    open(unit=123, file='err_10.dat', action='write',iostat=ierr,position='append')
    write(123,'(i5,1x,e19.12,1x)') ny, MAXVAL( ABS(df - dfdy_exact))
    close(123)
    call write_file_y('dfdy_10.dat', ny, x, y, z, f, df, dfdy_exact)

    !! -- 1, 1 -- !!
    call method1 % ddy(f,df,1,1)
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))
    open(unit=123, file='err_11.dat', action='write',iostat=ierr,position='append')
    write(123,'(i5,1x,e19.12,1x)') ny, MAXVAL( ABS(df - dfdy_exact))
    close(123)
    call write_file_y('dfdy_11.dat', ny, x, y, z, f, df, dfdy_exact)

    !! -- 1, -1 -- !!
    call method1 % ddy(f,df,1,-1)
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))
    open(unit=123, file='err_1m1.dat', action='write',iostat=ierr,position='append')
    write(123,'(i5,1x,e19.12,1x)') ny, MAXVAL( ABS(df - dfdy_exact))
    close(123)
    call write_file_y('dfdy_1m1.dat', ny, x, y, z, f, df, dfdy_exact)

    !! -- -1, 0 -- !!
    call method1 % ddy(f,df,-1,0)
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))
    open(unit=123, file='err_m10.dat', action='write',iostat=ierr,position='append')
    write(123,'(i5,1x,e19.12,1x)') ny, MAXVAL( ABS(df - dfdy_exact))
    close(123)
    call write_file_y('dfdy_m10.dat', ny, x, y, z, f, df, dfdy_exact)

    !! -- -1, 1 -- !!
    call method1 % ddy(f,df,-1,1)
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))
    open(unit=123, file='err_m11.dat', action='write',iostat=ierr,position='append')
    write(123,'(i5,1x,e19.12,1x)') ny, MAXVAL( ABS(df - dfdy_exact))
    close(123)
    call write_file_y('dfdy_m11.dat', ny, x, y, z, f, df, dfdy_exact)

    !! -- -1, -1 -- !!
    call method1 % ddy(f,df,-1,-1)
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))
    open(unit=123, file='err_m1m1.dat', action='write',iostat=ierr,position='append')
    write(123,'(i5,1x,e19.12,1x)') ny, MAXVAL( ABS(df - dfdy_exact))
    close(123)
    call write_file_y('dfdy_m1m1.dat', ny, x, y, z, f, df, dfdy_exact)

    call tic() 
    call method1 % ddy(f,df,0,0)
    !call toc ("Time to get the y derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))

    !print*, "==========================================="
    !print*, "Now trying METHOD 2: CD06"
    !print*, "==========================================="
    !call tic() 
    !call method2 % ddx(f,df)
    !call toc ("Time to get the x derivative:")
    !print*, "Maximum error = ", MAXVAL( ABS(df - dfdx_exact))

    !call tic() 
    !call method2 % ddy(f,df)
    !call toc ("Time to get the y derivative:")
    !print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))

    !call tic() 
    !call method2 % ddz(f,df)
    !call toc ("Time to get the z derivative:")
    !print*, "Maximum error = ", MAXVAL( ABS(df - dfdz_exact))


    !print*, "==========================================="
    !print*, "Now trying METHOD 3: FOUR"
    !print*, "==========================================="
    !call tic() 
    !call method3 % ddx(f,df)
    !call toc ("Time to get the x derivative:")
    !print*, "Maximum error = ", MAXVAL( ABS(df - dfdx_exact))

    !call tic() 
    !call method3 % ddy(f,df)
    !call toc ("Time to get the y derivative:")
    !print*, "Maximum error = ", MAXVAL( ABS(df - dfdy_exact))

    !call tic() 
    !call method3 % ddz(f,df)
    !call toc ("Time to get the z derivative:")
    !print*, "Maximum error = ", MAXVAL( ABS(df - dfdz_exact))



    print*, "---------------------------------------------"


    print*, "SECOND DERIVATIVE TESTS" 
    
    print*, "==========================================="
    print*, "Now trying METHOD 1: CD10"
    print*, "==========================================="
    call tic() 
    call method1 % d2dx2(f,df,0,0)
    !call toc ("Time to get the x derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - d2fdx2_exact))
   

    call tic() 
    call method1 % d2dy2(f,df,0,0)
    !call toc ("Time to get the y derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - d2fdy2_exact))
    linf = MAXVAL( ABS(df - d2fdy2_exact))
    l2norm = SQRT(SUM(ABS(df - d2fdy2_exact)**2)/(nx*ny*nz))
    errcen =  ABS(df(nx/2,ny/2,nz/2) - d2fdy2_exact(nx/2,ny/2,nz/2))
    open(unit=123, file='err_00_2ndder.dat', action='write',iostat=ierr,position='append')
    write(123,'(i5,1x,e19.12,1x,e19.12,1x,e19.12,1x)') ny, linf, l2norm, errcen
    close(123)
    call write_file_y('d2fdy2_00.dat', ny, x, y, z, f, df, d2fdy2_exact)

    call tic() 
    call method1 % d2dz2(f,df,0,0)
    !call toc ("Time to get the z derivative:")
    print*, "Maximum error = ", MAXVAL( ABS(df - d2fdz2_exact))


    !!print*, "==========================================="
    !!print*, "Now trying METHOD 2: CD06"
    !!print*, "==========================================="
    !!call tic() 
    !!call method2 % d2dx2(f,df)
    !!call toc ("Time to get the x derivative:")
    !!print*, "Maximum error = ", MAXVAL( ABS(df - d2fdx2_exact))

    !!call tic() 
    !!call method2 % d2dy2(f,df)
    !!call toc ("Time to get the y derivative:")
    !!print*, "Maximum error = ", MAXVAL( ABS(df - d2fdy2_exact))

    !!call tic() 
    !!call method2 % d2dz2(f,df)
    !!call toc ("Time to get the z derivative:")
    !!print*, "Maximum error = ", MAXVAL( ABS(df - d2fdz2_exact))


    !print*, "==========================================="
    !print*, "Now trying METHOD 3: FOUR"
    !print*, "==========================================="
    !call tic() 
    !call method3 % d2dx2(f,df)
    !call toc ("Time to get the x derivative:")
    !print*, "Maximum error = ", MAXVAL( ABS(df - d2fdx2_exact))

    !call tic() 
    !call method3 % d2dy2(f,df)
    !call toc ("Time to get the y derivative:")
    !print*, "Maximum error = ", MAXVAL( ABS(df - d2fdy2_exact))

    !call tic() 
    !call method3 % d2dz2(f,df)
    !call toc ("Time to get the z derivative:")
    !print*, "Maximum error = ", MAXVAL( ABS(df - d2fdz2_exact))
    
    
    print*, "=========================================="
    print*, " Destroying everything"
    print*, "=========================================="

    deallocate( x )
    deallocate( y )
    deallocate( z )
    deallocate( f )
    deallocate( df )
    deallocate( dfdx_exact )
    deallocate( dfdy_exact )
    deallocate( dfdz_exact )
    deallocate( d2fdx2_exact )
    deallocate( d2fdy2_exact )
    deallocate( d2fdz2_exact )
    call method1%destroy
    call method2%destroy
    call method3%destroy
    

end program

