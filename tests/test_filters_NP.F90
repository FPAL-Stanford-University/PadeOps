module transfer_funcs

    use kind_parameters, only: rkind
    use constants,       only: two
    implicit none

contains
    
    function GetTransferFunctionCF90(k) result(T)
        use cf90stuff, only: alpha90, beta90, a90, b90, c90, d90, e90
        real(rkind), intent(in) :: k
        real(rkind) :: T

        T = (a90 + two* b90*COS(k) + two*c90*COS(two*k) +two*d90*COS(3._rkind*k) + two*e90*COS(4._rkind*k) ) &
          / (1._rkind + two*alpha90*COS(k) + two*beta90*COS(two*k) )

    end function
    
    function GetTransferFunctionGaussian(k) result(T)
        use gaussianstuff, only: agf, bgf, cgf, dgf, egf
        real(rkind), intent(in) :: k
        real(rkind) :: T

        T = (agf + two* bgf*COS(k) + two*cgf*COS(two*k) +two*dgf*COS(3._rkind*k) + two*egf*COS(4._rkind*k) )

    end function

    function GetTransferFunctionLstsq(k) result(T)
        use lstsqstuff, only: als, bls, cls, dls, els
        real(rkind), intent(in) :: k
        real(rkind) :: T

        T = (als + two* bls*COS(k) + two*cls*COS(two*k) +two*dls*COS(3._rkind*k) + two*els*COS(4._rkind*k) )

    end function

end module


program test_filters_NP

    use kind_parameters, only: rkind
    use constants,       only: two,pi,zero,half
    use exits,           only: message
    use timer,           only: tic,toc
    use FiltersMod,      only: filters
    use transfer_funcs,  only: GetTransferFunctionCF90,GetTransferFunctionGaussian,GetTransferFunctionLstsq
    implicit none

    integer :: nx = 256, ny=256, nz=256
    logical :: xcentered=.true.

    type( filters ) :: method1, method2, method3
    real(rkind), dimension(:,:,:), allocatable :: x,y,z,f,df
    real(rkind) :: dx, dy, dz
    real(rkind), parameter :: omega = 1._rkind

    real(rkind) :: xfst, x1err, x2err, y1err, y2err, z1err, z2err, Lx, x3err, y3err, z3err
    integer :: xa=0,xb=0,iprob=1
    character(len=15) :: fname
    character :: cent_str

    real(rkind) :: k_norm_x,k_norm_y,k_norm_z,TF1_x,TF1_y,TF1_z,TF2_x,TF2_y,TF2_z,TF3_x,TF3_y,TF3_z

    integer :: i,j,k

    open(10,file='input.in',status='old',action='read')
    read(10,*) nx, ny, nz
    read(10,*) xcentered
    read(10,*) iprob
    close(10)

    allocate( x(nx,ny,nz) )
    allocate( y(nx,ny,nz) )
    allocate( z(nx,ny,nz) )
    allocate( f(nx,ny,nz) )
    allocate( df(nx,ny,nz) )

        !if(iprob==1 .or. iprob==2) then
        !    Lx = pi/two
        !else
        !    Lx = pi
        !endif
        !if(xcentered) then
        !    dx = Lx/real(nx,rkind)
        !    xfst = half*dx
        !else
        !    dx = Lx/real(nx-1,rkind)
        !    xfst = zero
        !endif

    if(xcentered) then
        dx = two*pi/real(nx,rkind)
        xfst = half*dx
    else
        dx = two*pi/real(nx-1,rkind)
        xfst = zero
    endif
    dy = two*pi/real(ny-1,rkind)
    dz = two*pi/real(nz-1,rkind)

    if(iprob==0) then
      xa = 1
      xb = 1
    else
      xa = -1
      xb = -1
    endif

    do k=1,nz
        do j=1,ny
            do i=1,nx
                x(i,j,k) = real(i-1,rkind)*dx + xfst
                y(i,j,k) = real(j-1,rkind)*dy
                z(i,j,k) = real(k-1,rkind)*dz
                f(i,j,k) = (iprob*sin(omega * x(i,j,k)) +(1-iprob)*cos(omega*x(i,j,k))) * sin(omega * y(i,j,k)) * sin(omega * z(i,j,k))
            end do
        end do 
    end do 

    k_norm_x = omega * dx
    TF1_x = GetTransferFunctionCF90    (k_norm_x)
    TF2_x = GetTransferFunctionGaussian(k_norm_x)
    TF3_x = GetTransferFunctionLstsq   (k_norm_x)
    
    k_norm_y = omega * dy
    TF1_y = GetTransferFunctionCF90    (k_norm_y)
    TF2_y = GetTransferFunctionGaussian(k_norm_y)
    TF3_y = GetTransferFunctionLstsq   (k_norm_y)
    
    k_norm_z = omega * dz
    TF1_z = GetTransferFunctionCF90    (k_norm_z)
    TF2_z = GetTransferFunctionGaussian(k_norm_z)
    TF3_z = GetTransferFunctionLstsq   (k_norm_z)

    print*, "Created initial data"

    call method1%init(          nx,     ny,    nz, &
                            .FALSE., .FALSE., .FALSE., &
                            "cf90", "cf90", "cf90", &
                            xcentered, .FALSE., .FALSE. )

    call method2%init(          nx,     ny,    nz, &
                            .FALSE., .FALSE.,  .FALSE., &
                            "gaussian", "gaussian", "gaussian", &
                            xcentered, .FALSE., .FALSE. )
    
    call method3%init(          nx,     ny,    nz, &
                            .FALSE., .FALSE.,  .FALSE., &
                            "lstsq", "lstsq", "lstsq", &
                            .FALSE., .FALSE., .FALSE. )

    print*, "Initialized all methods"
   
   
    print*, "==========================================="
    print*, "Now trying METHOD 1: CF90"
    print*, "==========================================="
    call tic() 
    call method1 % filterx(f,df,xa,xb)
    call toc ("Time to get the x filter:")
    call message("Maximum error", MAXVAL(df - TF1_x*f))
    x1err = MAXVAL(df-TF1_x*f)

    call tic() 
    call method1 % filtery(f,df,-1,-1)
    call toc ("Time to get the y filter:")
    call message("Maximum error", MAXVAL(df - TF1_y*f))
    y1err = MAXVAL(df-TF1_y*f)

    call tic() 
    call method1 % filterz(f,df,-1,-1)
    call toc ("Time to get the z filter:")
    call message("Maximum error", MAXVAL(df - TF1_z*f))
    z1err = MAXVAL(df-TF1_z*f)

    print*, "==========================================="
    print*, "Now trying METHOD 2: GAUSSIAN"
    print*, "==========================================="
    call tic() 
    call method2 % filterx(f,df,xa,xb)
    call toc ("Time to get the x filter:")
    call message("Maximum error", MAXVAL(df - TF2_x*f))
    x2err = MAXVAL(df-TF2_x*f)

    call tic() 
    call method2 % filtery(f,df,-1,-1)
    call toc ("Time to get the y filter:")
    call message("Maximum error", MAXVAL(df - TF2_y*f))
    y2err = MAXVAL(df-TF2_y*f)

    call tic() 
    call method2 % filterz(f,df,-1,-1)
    call toc ("Time to get the z filter:")
    call message("Maximum error", MAXVAL(df - TF2_z*f))
    z2err = MAXVAL(df-TF2_z*f)


    print*, "==========================================="
    print*, "Now trying METHOD 3: LSTSQ"
    print*, "==========================================="
    call tic() 
    call method3 % filterx(f,df,xa,xb)
    call toc ("Time to get the x filter:")
    call message("Maximum error", MAXVAL(df - TF3_x*f))
    x3err = MAXVAL(df-TF3_x*f)

    call tic() 
    call method3 % filtery(f,df,-1,-1)
    call toc ("Time to get the y filter:")
    call message("Maximum error", MAXVAL(df - TF3_y*f))
    y3err = MAXVAL(df-TF3_y*f)

    call tic() 
    call method3 % filterz(f,df,-1,-1)
    call toc ("Time to get the z filter:")
    call message("Maximum error", MAXVAL(df - TF3_z*f))
    z3err = MAXVAL(df-TF3_z*f)

    cent_str = 'F'
    if(xcentered) cent_str='T'
    write(fname,'(a7,i2.2,2(a1),a4)') 'errors_', iprob, '_', cent_str, '.dat'
    open(10,file=fname,status='unknown',action='write',position='append')
    write(10,'(3(i4.4,1x),9(e12.5,1x))') nx, ny, nz, x1err, x2err, x3err, y1err, y2err, y3err, z1err, z2err, z3err
    close(10)
    print *, fname

    print*, "=========================================="
    print*, "Destroying everything"
    print*, "=========================================="

    call method1%destroy
    call method2%destroy
    call method3%destroy
    deallocate( x )
    deallocate( y )
    deallocate( z )
    deallocate( f )
    deallocate( df )
    

end program
