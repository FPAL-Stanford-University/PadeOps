module subroutines

  use kind_parameters, only: rkind, clen
  use constants      , only: two, zero
  use exits,           only: gracefulExit, message
  implicit none

contains

subroutine write_file_z(fname, z, f, df, dfdz_exact, ii, jj)
    use kind_parameters, only: rkind, clen

    implicit none
    character(len=clen), intent(in) :: fname
    integer, intent(in) :: ii, jj
    real(rkind), dimension(:,:,:), intent(in) :: z, f, df, dfdz_exact
    integer :: k

    write(*,*) 'size(z,3) = ', size(z,3)
    open(10,file=fname,status='unknown')
    do k = 1, size(z,3)
      write(10,'(6(e22.15,1x))') z(ii,jj,k), f(ii,jj,k), df(ii,jj,k), dfdz_exact(ii,jj,k)
    enddo
    close(10)

end subroutine write_file_z

subroutine write_file_y(fname, y, f, df, dfdy_exact, ii, kk)
    use kind_parameters, only: rkind, clen

    implicit none
    character(len=clen), intent(in) :: fname
    integer, intent(in) :: ii, kk
    real(rkind), dimension(:,:,:), intent(in) :: y, f, df, dfdy_exact
    integer :: j

    write(*,*) 'size(y,2) = ', size(y,2)
    open(10,file=fname,status='unknown')
    do j = 1, size(y,2)
      write(10,'(6(e22.15,1x))') y(ii,j,kk), f(ii,j,kk), df(ii,j,kk), dfdy_exact(ii,j,kk)
    enddo
    close(10)

end subroutine write_file_y

subroutine write_file_x(fname, x, f, df, dfdx_exact, jj, kk)

    implicit none
    character(len=clen), intent(in) :: fname
    integer, intent(in) :: jj, kk
    real(rkind), dimension(:,:,:), intent(in) :: x, f, df, dfdx_exact
    integer :: i

    open(10,file=fname,status='unknown')
    do i = 1, size(x,1)
      write(10,'(6(e22.15,1x))') x(i,jj,kk), f(i,jj,kk), df(i,jj,kk), dfdx_exact(i,jj,kk)
    enddo
    close(10)
end subroutine write_file_x

function GetTransferFunction(k, method) result(T)
    use cf90stuff, only: alpha90, beta90, a90, b90, c90, d90, e90
    use gaussianstuff, only: agf, bgf, cgf, dgf, egf
    use lstsqstuff, only: als, bls, cls, dls, els
    real(rkind), intent(in) :: k
    character(len=*) , intent(in) :: method
    real(rkind) :: T

    select case (method)
    case("cf90")
      T = (a90 + two* b90*COS(k) + two*c90*COS(two*k) +two*d90*COS(3._rkind*k) + two*e90*COS(4._rkind*k) ) &
        / (1._rkind + two*alpha90*COS(k) + two*beta90*COS(two*k) )
    case("gaussian")
      T = (agf + two* bgf*COS(k) + two*cgf*COS(two*k) +two*dgf*COS(3._rkind*k) + two*egf*COS(4._rkind*k) )
    case("lstsq")
      T = (als + two* bls*COS(k) + two*cls*COS(two*k) +two*dls*COS(3._rkind*k) + two*els*COS(4._rkind*k) )
    case default
        call GracefulExit("Transfer function not implemented for your filter method", 52)
    end select 

end function

    
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

program test_multiblock_fil
    use mpi
    use kind_parameters, only : rkind, clen
    use decomp_2d
    use decomp_2d_io, only: decomp_2d_write_one
    use MultiBlockTopologyMod, only: multiblocktopol
    use constants, only: pi
    use FiltersMod,      only: filters
    use reductions, only: P_MAXVAL, p_sum
    use subroutines

    implicit none

    real(rkind), dimension(:,:,:,:), allocatable, target :: mesh, mesh_curv, xbuf, zbuf
    real(rkind), dimension(:,:,:),   allocatable :: f, df_x, df_y, df_z, ff_x, ff_y, ff_z
    real(rkind), dimension(:,:,:), pointer :: x, y, z, xi, eta, zeta, xtmp1, xtmp2, ztmp1, ztmp2
    real(rkind), dimension(:,:,:), pointer :: x_in_X, ff_x_in_X, z_in_Z, ff_z_in_Z
    type(decomp_info) :: gp
    type(multiblocktopol), allocatable :: mbtopology
    type(filters),     allocatable :: fil

    integer :: nx = 256, ny = 256, nz = 256
    real(rkind) :: Lx = 10.0d0, Ly = 2.0d0, Lz = 4.0d0
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k
    real(rkind) :: dx, dy, dz
    character(len=clen) :: inputfile, outputfile
    logical :: periodicx=.true., periodicy=.true., periodicz=.true.
    logical :: xmetric=.false., ymetric=.false., zmetric=.false.
    character(len=clen) :: filter_x="cf90", filter_y="cf90", filter_z = "cf90"
    real(rkind) :: linf, l2norm, errcen, omega = 12.0_rkind, TF1_x, TF1_y, TF1_z
    integer :: ii, kk, jj, x_bc(2)=(/0, 0/), y_bc(2)=(/0, 0/), z_bc(2)=(/0, 0/)
    integer :: rank_debug = 0, num_x_freq = 3, num_y_freq = 1, num_z_freq = 1
    real(rkind) :: om_x_arr(3) = (/1.0d0, 4.0d0, 8.0d0 /), amp_x_arr(3) = (/1.0d0, 0.2d0, 0.3d0/)
    real(rkind) :: om_y_arr(3) = (/1.0d0, 4.0d0, 8.0d0 /), amp_y_arr(3) = (/1.0d0, 0.2d0, 0.3d0/)
    real(rkind) :: om_z_arr(3) = (/1.0d0, 4.0d0, 8.0d0 /), amp_z_arr(3) = (/1.0d0, 0.2d0, 0.3d0/)
    real(rkind), allocatable, dimension(:) :: TF_x_arr, TF_y_arr, TF_z_arr
    integer :: iquery=-1,jquery=-1, kquery=-1, ifx, ify, ifz
    logical :: test_x=.true., test_y=.false., test_z=.false.

    namelist /INPUT/ nx, ny, nz, Lx, Ly, Lz, prow, pcol, rank_debug, periodicx, periodicy, periodicz, &
                     om_x_arr, amp_x_arr, om_y_arr, amp_y_arr, om_z_arr, amp_z_arr, iquery, jquery, kquery, &
                     x_bc, y_bc, z_bc, test_x, test_y, test_z, num_x_freq, num_y_freq, num_z_freq, &
                     filter_x, filter_y, filter_z

    call MPI_Init(ierr)

    write(inputfile,'(a)') 'input_test_multiblock_filter.dat'
    open(unit=123, file=inputfile, form='FORMATTED', iostat=ierr)
    read(unit=123, NML=INPUT)
    close(123)

    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)

    ! Initialize mesh arrays
    allocate(mesh(gp%ysz(1), gp%ysz(2), gp%ysz(3), 3) )
    allocate(xbuf(gp%xsz(1), gp%xsz(2), gp%xsz(3), 4) )
    allocate(zbuf(gp%zsz(1), gp%zsz(2), gp%zsz(3), 4) )
    allocate(mesh_curv(gp%ysz(1), gp%ysz(2), gp%ysz(3), 3) )
    allocate(   f(gp%ysz(1), gp%ysz(2), gp%ysz(3)))
    allocate(df_x(gp%ysz(1), gp%ysz(2), gp%ysz(3)))  ! numerical filtered field
    allocate(df_y(gp%ysz(1), gp%ysz(2), gp%ysz(3)))  ! numerical filtered field
    allocate(df_z(gp%ysz(1), gp%ysz(2), gp%ysz(3)))  ! numerical filtered field
    allocate(ff_x(gp%ysz(1), gp%ysz(2), gp%ysz(3)))  ! exact     filtered field
    allocate(ff_y(gp%ysz(1), gp%ysz(2), gp%ysz(3)))  ! exact     filtered field
    allocate(ff_z(gp%ysz(1), gp%ysz(2), gp%ysz(3)))  ! exact     filtered field

    x => mesh(:,:,:,1);  y => mesh(:,:,:,2);  z => mesh(:,:,:,3)
    xi => mesh_curv(:,:,:,1);  eta => mesh_curv(:,:,:,2);  zeta => mesh_curv(:,:,:,3)
    xtmp1 => xbuf(:,:,:,1);  xtmp2 => xbuf(:,:,:,2)
    ztmp1 => zbuf(:,:,:,1);  ztmp2 => zbuf(:,:,:,2)
    x_in_X => xbuf(:,:,:,3);  ff_x_in_X => xbuf(:,:,:,4)
    z_in_Z => zbuf(:,:,:,3);  ff_z_in_Z => zbuf(:,:,:,4)

    if(periodicx) then
        dx = Lx/real(max(nx,1),rkind)
    else
        dx = Lx/real(max(nx-1,1),rkind)
    endif

    if(periodicy) then
        dy = Ly/real(max(ny,1),rkind)
    else
        dy = Ly/real(max(ny-1,1),rkind)
    endif

    if(periodicz) then
        dz = Lz/real(max(nz,1),rkind)
    else
        dz = Lz/real(max(nz-1,1),rkind)
    endif


    ! Generate mesh
    do k = 1,gp%ysz(3)
        do j = 1,gp%ysz(2)
            do i = 1,gp%ysz(1)
                x(i,j,k) = real(gp%yst(1) - 1 + i - 1, rkind)*dx
                y(i,j,k) = real(gp%yst(2) - 1 + j - 1, rkind)*dy
                z(i,j,k) = real(gp%yst(3) - 1 + k - 1, rkind)*dz
            end do
        end do
    end do

    f = zero
    if(test_x) then
      do k = 1, gp%ysz(3)
       do j = 1, gp%ysz(2)
        do i = 1, gp%ysz(1)
         do ifx = 1, num_x_freq
           f(i,j,k) = f(i,j,k) + amp_x_arr(ifx) * cos(om_x_arr(ifx) * x(i,j,k)) * cos(omega * y(i,j,k)) * cos(omega * z(i,j,k))
         end do 
        end do
       end do
      end do
    endif
    if(test_y) then
      do k = 1, gp%ysz(3)
       do j = 1, gp%ysz(2)
        do i = 1, gp%ysz(1)
         do ify = 1, num_y_freq
           f(i,j,k) = f(i,j,k) + amp_y_arr(ify) * cos(om_y_arr(ify) * y(i,j,k)) * cos(omega * x(i,j,k)) * cos(omega * z(i,j,k))
         end do 
        end do
       end do
      end do
    endif
    if(test_z) then
      do k = 1, gp%ysz(3)
       do j = 1, gp%ysz(2)
        do i = 1, gp%ysz(1)
         do ifz = 1, num_z_freq
           f(i,j,k) = f(i,j,k) + amp_z_arr(ifz) * cos(om_z_arr(ifz) * z(i,j,k)) * cos(omega * x(i,j,k)) * cos(omega * y(i,j,k))
         end do 
        end do
       end do
      end do
    endif
    call decomp_2d_write_one(2, f,   "fval_original.out", gp)

    ! get transfer functions for discrete frequencies and use them to get the
    ! `exact' filtered solution. This is the exact filter only for periodic 
    ! boundaries and integer frequencies
    allocate(TF_x_arr(num_x_freq))
    allocate(TF_y_arr(num_y_freq))
    allocate(TF_z_arr(num_z_freq))

    if(test_x) then
      do ifx = 1, num_x_freq
        TF_x_arr(ifx) = GetTransferFunction (om_x_arr(ifx) * dx, filter_x)
      end do
    else
      TF_x_arr(1) = GetTransferFunction (omega * dx, filter_x)
    endif
    if(test_y) then
      do ify = 1, num_y_freq
        TF_y_arr(ify) = GetTransferFunction (om_y_arr(ify) * dy, filter_y)
      end do
    else
      TF_y_arr(1) = GetTransferFunction (omega * dy, filter_y)
    endif
    if(test_z) then
      do ifz = 1, num_z_freq
        TF_z_arr(ifz) = GetTransferFunction (om_z_arr(ifz) * dz, filter_z)
      end do
    else
      TF_z_arr(1) = GetTransferFunction (omega * dz, filter_z)
    endif

    ff_x = zero
    if(test_x) then
      do k = 1, gp%ysz(3)
       do j = 1, gp%ysz(2)
        do i = 1, gp%ysz(1)
         do ifx = 1, num_x_freq
           ff_x(i,j,k) = ff_x(i,j,k) + TF_x_arr(ifx) * amp_x_arr(ifx) * cos(om_x_arr(ifx) * x(i,j,k)) * cos(omega * y(i,j,k)) * cos(omega * z(i,j,k))
         end do 
        end do
       end do
      end do
    else
      ff_x = TF_x_arr(1) * f
    endif

    ff_y = zero
    if(test_y) then
      do k = 1, gp%ysz(3)
       do j = 1, gp%ysz(2)
        do i = 1, gp%ysz(1)
         do ify = 1, num_y_freq
           ff_y(i,j,k) = ff_y(i,j,k) + TF_y_arr(ify) * amp_y_arr(ify) * cos(om_y_arr(ify) * y(i,j,k)) * cos(omega * x(i,j,k)) * cos(omega * z(i,j,k))
         end do 
        end do
       end do
      end do
    else
      ff_y = TF_y_arr(1) * f
    endif

    ff_z = zero
    if(test_z) then
      do k = 1, gp%ysz(3)
       do j = 1, gp%ysz(2)
        do i = 1, gp%ysz(1)
         do ifz = 1, num_z_freq
           ff_z(i,j,k) = ff_z(i,j,k) + TF_z_arr(ifz) * amp_z_arr(ifz) * cos(om_z_arr(ifz) * z(i,j,k)) * cos(omega * y(i,j,k)) * cos(omega * x(i,j,k))
         end do 
        end do
       end do
      end do
    else
      ff_z = TF_z_arr(1) * f
    endif
    
    allocate(mbtopology)
    call mbtopology%init(gp, mesh, inputfile, xbuf, zbuf)

    ff_x = ff_x * mbtopology%mask
    ff_y = ff_y * mbtopology%mask
    ff_z = ff_z * mbtopology%mask


    allocate(fil)
    call fil%init(gp, periodicx, periodicy, periodicz, filter_x, filter_y, filter_z, mbtopology)

    !!!! Test the filtering along x direction !!!!!!
    if(test_x) then
        call transpose_y_to_x(x,x_in_X,gp)
        call transpose_y_to_x(ff_x,ff_x_in_X,gp)

        ! x-filtering 
        xtmp2 = zero
        call transpose_y_to_x(f, xtmp1, gp)
        call fil%filterx(xtmp1, xtmp2, x_bc(1), x_bc(2))
        call transpose_x_to_y(xtmp2, df_x,gp)
        df_x = df_x * mbtopology%mask

        !!! quantify the error in x-filtering !!!
        ! if not assigned through inputfile, set query location !
        if(iquery == -1) iquery = gp%xsz(1)/2
        if(jquery == -1) jquery = gp%xsz(2)/2
        if(kquery == -1) kquery = gp%xsz(3)/2

        ! write out the entire profile along one line !
        if(nrank==rank_debug) then
          if(periodicx) then
            write(outputfile,'(a,a,a,3(i4.4,a),a)') 'filx_per_',trim(filter_x),'_',nrank,'_',jquery,'_',kquery,'.dat'
          else
            write(outputfile,'(a,a,a,3(i4.4,a),a)') 'filx_nonper_',trim(filter_x),'_',nrank,'_',jquery,'_',kquery,'.dat'
          endif
          call write_file_x(outputfile, x_in_X, xtmp1, xtmp2, ff_x_in_X, jquery, kquery)
        endif

        ! write out error norms and error at one point !
        linf = P_MAXVAL( MAXVAL(ABS(df_x - ff_x)))
        l2norm = SQRT(P_SUM(ABS(df_x - ff_x)**2)/(nx*ny*nz))
        ii = gp%ysz(1)/2; jj = gp%ysz(2)/2; kk = gp%ysz(3)/2
        if(nrank==rank_debug) then
            errcen =  ABS(df_x(ii,jj,kk) - ff_x(ii,jj,kk))
            if(periodicx) then
              write(outputfile,'(a,a,a)') 'xfil_err_per_',trim(filter_x),'.dat'
            else
              write(outputfile,'(a,a,a)') 'xfil_err_nonper_',trim(filter_x),'.dat'
            endif
            open(unit=123, file=outputfile, action='write',iostat=ierr,position='append')
            write(123,'(i5,1x,e19.12,1x,e19.12,1x,e19.12,1x)') nx, linf, l2norm, errcen
            close(123)
        endif

        ! write the entire 3D field !
        call decomp_2d_write_one(2, df_x,"f_x_numerical.out", gp)
        call decomp_2d_write_one(2, ff_x,"f_x_transfunc.out", gp)
    endif
    
    !!!! Test the filtering along y direction !!!!!!
    if(test_y) then
        ! y-filtering 
        call fil%filtery(f, df_y, y_bc(1), y_bc(2))
        df_y = df_y * mbtopology%mask

        !!! quantify the error in y-filtering !!!
        ! if not assigned through inputfile, set query location !
        if(iquery == -1) iquery = gp%ysz(1)/2
        if(jquery == -1) jquery = gp%ysz(2)/2
        if(kquery == -1) kquery = gp%ysz(3)/2

        ! write out the entire profile along one line !
        if(nrank==rank_debug) then
          if(periodicy) then
            write(outputfile,'(a,a,a,3(i4.4,a),a)') 'fily_per_',trim(filter_y),'_',nrank,'_',iquery,'_',kquery,'.dat'
          else
            write(outputfile,'(a,a,a,3(i4.4,a),a)') 'fily_nonper_',trim(filter_y),'_',nrank,'_',iquery,'_',kquery,'.dat'
          endif
          call write_file_y(outputfile, y, f, df_y, ff_y, iquery, kquery)
        endif

        ! write out error norms and error at one point !
        linf = P_MAXVAL( MAXVAL(ABS(df_y - ff_y)))
        l2norm = SQRT(P_SUM(ABS(df_y - ff_y)**2)/(nx*ny*nz))
        ii = max(min(gp%ysz(1), iquery), 1)
        jj = max(min(gp%ysz(2), jquery), 1)
        kk = max(min(gp%ysz(3), kquery), 1)
        if(nrank==rank_debug) then
            errcen =  ABS(df_y(ii,jj,kk) - ff_y(ii,jj,kk))
            if(periodicy) then
              write(outputfile,'(a,a,a)') 'yfil_err_per_',trim(filter_y),'.dat'
            else
              write(outputfile,'(a,a,a)') 'yfil_err_nonper_',trim(filter_y),'.dat'
            endif
            open(unit=123, file=outputfile, action='write',iostat=ierr,position='append')
            write(123,'(i5,1x,e19.12,1x,e19.12,1x,e19.12,1x)') ny, linf, l2norm, errcen
            close(123)
        endif

        ! write the entire 3D field !
        call decomp_2d_write_one(2, df_y,"f_y_numerical.out", gp)
        call decomp_2d_write_one(2, ff_y,"f_y_transfunc.out", gp)
    endif

    !!!! Test the filtering along z direction !!!!!!
    if(test_z) then
        call transpose_y_to_z(z,z_in_Z,gp)
        call transpose_y_to_z(ff_z,ff_z_in_Z,gp)

        ! x-filtering 
        ztmp2 = zero
        call transpose_y_to_z(f, ztmp1, gp)
        call fil%filterz(ztmp1, ztmp2, z_bc(1), z_bc(2))
        call transpose_z_to_y(ztmp2, df_z,gp)
        df_z = df_z * mbtopology%mask

        !!! quantify the error in x-filtering !!!
        ! if not assigned through inputfile, set query location !
        if(iquery == -1) iquery = gp%zsz(1)/2
        if(jquery == -1) jquery = gp%zsz(2)/2
        if(kquery == -1) kquery = gp%zsz(3)/2

        ! write out the entire profile along one line !
        if(nrank==rank_debug) then
          if(periodicz) then
            write(outputfile,'(a,a,a,3(i4.4,a),a)') 'filz_per_',trim(filter_z),'_',nrank,'_',iquery,'_',jquery,'.dat'
          else
            write(outputfile,'(a,a,a,3(i4.4,a),a)') 'filz_nonper_',trim(filter_z),'_',nrank,'_',iquery,'_',jquery,'.dat'
          endif
          call write_file_z(outputfile, z_in_Z, ztmp1, ztmp2, ff_z_in_Z, iquery, jquery)
        endif

        ! write out error norms and error at one point !
        linf = P_MAXVAL( MAXVAL(ABS(df_z - ff_z)))
        l2norm = SQRT(P_SUM(ABS(df_z - ff_z)**2)/(nx*ny*nz))
        ii = gp%ysz(1)/2; jj = gp%ysz(2)/2; kk = gp%ysz(3)/2
        if(nrank==rank_debug) then
            errcen =  ABS(df_z(ii,jj,kk) - ff_z(ii,jj,kk))
            if(periodicz) then
              write(outputfile,'(a,a,a)') 'zfil_err_per_',trim(filter_z),'.dat'
            else
              write(outputfile,'(a,a,a)') 'zfil_err_nonper_',trim(filter_z),'.dat'
            endif
            open(unit=123, file=outputfile, action='write',iostat=ierr,position='append')
            write(123,'(i5,1x,e19.12,1x,e19.12,1x,e19.12,1x)') nz, linf, l2norm, errcen
            close(123)
        endif

        ! write the entire 3D field !
        call decomp_2d_write_one(2, df_z,"f_z_numerical.out", gp)
        call decomp_2d_write_one(2, ff_z,"f_z_transfunc.out", gp)
    endif
    
    call fil%destroy()
    deallocate(fil)

    call mbtopology%destroy()
    deallocate(mbtopology)

    nullify(x, y, z, xi, eta, zeta, xtmp1, xtmp2, ztmp1, ztmp2)
    nullify(z_in_Z, ff_z_in_Z)
    nullify(x_in_X, ff_x_in_X)
    deallocate(TF_x_arr)
    deallocate(TF_y_arr)
    deallocate(TF_z_arr)
    deallocate(f, df_x, df_y, df_z, ff_x, ff_y, ff_z)
    deallocate(zbuf, xbuf, mesh_curv, mesh)
    call decomp_2d_finalize
    call MPI_Finalize(ierr)

end program

