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


program test_filters_rangeomega

    use kind_parameters, only: rkind,clen
    use constants,       only: two,pi
    use exits,           only: message
    use timer,           only: tic,toc
    use FiltersMod,      only: filters
    use transfer_funcs,  only: GetTransferFunctionCF90,GetTransferFunctionGaussian,GetTransferFunctionLstsq
    implicit none

    integer :: nx = 48 , ny=48 , nz=48 , ii, jj, kk, iwv=1, ierr


    type( filters ) :: method1, method2, method3
    real(rkind), dimension(:,:,:), allocatable :: x,y,z,f,df
    real(rkind) :: dx, dy, dz
    real(rkind) :: omega = 1.0_rkind, om_min = 1.0_rkind, om_max = 24.0_rkind, om_step = 1.0_rkind
    character(len=clen) :: outputfile

    real(rkind) :: k_norm_x,k_norm_y,k_norm_z,TF1_x,TF1_y,TF1_z,TF2_x,TF2_y,TF2_z,TF3_x,TF3_y,TF3_z

    integer :: i,j,k, om_num, iounit2=234, iounit3=345

    allocate( x(nx,ny,nz) )
    allocate( y(nx,ny,nz) )
    allocate( z(nx,ny,nz) )
    allocate( f(nx,ny,nz) )
    allocate( df(nx,ny,nz) )

    dx = two*pi/real(nx,rkind)/two
    dy = two*pi/real(ny,rkind)
    dz = two*pi/real(nz,rkind)

    do k=1,nz
        do j=1,ny
            do i=1,nx
                x(i,j,k) = real(i-1,rkind)*dx
                y(i,j,k) = real(j-1,rkind)*dy
                z(i,j,k) = real(k-1,rkind)*dz
            end do
        end do 
    end do 

    call method1%init(          nx,     ny,    nz, &
                            .TRUE., .TRUE., .TRUE., &
                            "cf90", "cf90", "cf90" )

    call method2%init(          nx,     ny,    nz, &
                            .TRUE., .TRUE.,  .TRUE., &
                            "gaussian", "gaussian", "gaussian" )
    
    call method3%init(          nx,     ny,    nz, &
                            .TRUE., .TRUE.,  .TRUE., &
                            "lstsq", "lstsq", "lstsq" )

    print*, "Initialized all methods"

    write(outputfile,'(a,i4.4,a,i4.4,a)') 'xfil_cf90_nx',nx,'_err_per.dat'
    open(unit=iounit2, file=outputfile, action='write',iostat=ierr,position='append')

    write(outputfile,'(a,i4.4,a,i4.4,a)') 'xfil_gaus_nx',nx,'_err_per.dat'
    open(unit=iounit3, file=outputfile, action='write',iostat=ierr,position='append')

    om_num = nint((om_max-om_min)/(om_step + 1.0d-12)) + 1
    do iwv = 1, om_num
      ! set the wavenumber
      omega = om_min + (iwv-1) * om_step

      ! set the function with the single wavenumber
      do k=1,nz
        do j=1,ny
          do i=1,nx
            f(i,j,k) = cos(omega * x(i,j,k)) * cos(omega * y(i,j,k)) * cos(omega * z(i,j,k))
          end do
        end do 
      end do 

      ! get transfer functions
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

      print '(a,1x,i4.4,1x,e19.12)', "Created initial data", iwv, omega
   
      print*, "==========================================="
      print*, "Now trying METHOD 1: CF90"
      print*, "==========================================="
      call tic() 
      call method1 % filterx(f,df)
      call toc ("Time to get the x filter:")
      call message("Maximum error", MAXVAL(df - TF1_x*f))
      write(iounit2,'(i4.4,1x,2(e19.12,1x))') iwv, omega, MAXVAL(df - TF1_x*f)

      ! write x, f, df along an arbitrary line
      write(outputfile,'(a,i4.4,a,i4.4,a)') 'xfil_cf90_nx',nx,'_iomega',iwv,'_per.dat'
      jj = 4; kk = 7
      open(unit=123, file=outputfile, action='write',iostat=ierr,position='append')
      do ii = 1, nx
        write(123,'(4(e19.12,1x))') x(ii,jj,kk), f(ii,jj,kk), df(ii,jj,kk), omega
      enddo
      close(123)
      

      call tic() 
      call method1 % filtery(f,df)
      call toc ("Time to get the y filter:")
      call message("Maximum error", MAXVAL(df - TF1_y*f))

      call tic() 
      call method1 % filterz(f,df)
      call toc ("Time to get the z filter:")
      call message("Maximum error", MAXVAL(df - TF1_z*f))


      print*, "==========================================="
      print*, "Now trying METHOD 2: GAUSSIAN"
      print*, "==========================================="
      call tic() 
      call method2 % filterx(f,df)
      call toc ("Time to get the x filter:")
      call message("Maximum error", MAXVAL(df - TF2_x*f))
      write(iounit3,'(i4.4,1x,2(e19.12,1x))') iwv, omega, MAXVAL(df - TF2_x*f)

      ! write x, f, df along an arbitrary line
      write(outputfile,'(a,i4.4,a,i4.4,a)') 'xfil_gaus_nx',nx,'_iomega',iwv,'_per.dat'
      jj = 4; kk = 7
      open(unit=123, file=outputfile, action='write',iostat=ierr,position='append')
      do ii = 1, nx
        write(123,'(4(e19.12,1x))') x(ii,jj,kk), f(ii,jj,kk), df(ii,jj,kk), omega
      enddo
      close(123)
      
      call tic() 
      call method2 % filtery(f,df)
      call toc ("Time to get the y filter:")
      call message("Maximum error", MAXVAL(df - TF2_y*f))

      call tic() 
      call method2 % filterz(f,df)
      call toc ("Time to get the z filter:")
      call message("Maximum error", MAXVAL(df - TF2_z*f))


      print*, "==========================================="
      print*, "Now trying METHOD 3: LSTSQ"
      print*, "==========================================="
      call tic() 
      call method3 % filterx(f,df)
      call toc ("Time to get the x filter:")
      call message("Maximum error", MAXVAL(df - TF3_x*f))

      call tic() 
      call method3 % filtery(f,df)
      call toc ("Time to get the y filter:")
      call message("Maximum error", MAXVAL(df - TF3_y*f))

      call tic() 
      call method3 % filterz(f,df)
      call toc ("Time to get the z filter:")
      call message("Maximum error", MAXVAL(df - TF3_z*f))

    enddo

    close(iounit2)
    close(iounit3)

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
