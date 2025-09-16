module subroutines


contains

!subroutine write_file_y(fname, ny, x, y, z, f, df, dfdx_exact)
!    use kind_parameters, only: rkind, clen
!
!    implicit none
!    character(len=clen), intent(in) :: fname
!    integer, intent(in) :: ny
!    real(rkind), dimension(:,:,:), intent(in) :: x, y, z, f, df, dfdx_exact
!    integer :: i, j, k
!
!    i = 7; k = 9;
!    open(10,file=fname,status='unknown')
!    do j = 1, ny
!      write(10,'(6(e19.12,1x))') x(i,j,k), y(i,j,k), z(i,j,k), f(i,j,k), df(i,j,k), dfdx_exact(i,j,k)
!    enddo
!    close(10)
!
!end subroutine write_file_y

subroutine write_file_x(fname, x, f, df, dfdx_exact, jj, kk)
    use kind_parameters, only: rkind, clen

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

end module

program test_multiblock_setup
    use mpi
    use kind_parameters, only : rkind, clen
    use decomp_2d
    use MultiBlockTopologyMod, only: multiblocktopol
    use constants, only: pi
    use DerivativesMod, only: derivatives
    use reductions, only: P_MAXVAL, p_sum
    use subroutines

    implicit none

    real(rkind), dimension(:,:,:,:), allocatable, target :: mesh, mesh_curv, xbuf, zbuf
    real(rkind), dimension(:,:,:),   allocatable :: f, dfdx, dfdy, dfdz, fxex, fyex, fzex
    real(rkind), dimension(:,:,:), pointer :: x, y, z, xi, eta, zeta, xtmp1, xtmp2, ztmp1, ztmp2
    real(rkind), dimension(:,:,:), pointer :: x_in_X, fxex_in_X
    type(decomp_info) :: gp
    type(multiblocktopol), allocatable :: mbtopology
    type(derivatives),     allocatable :: der

    integer :: nx = 256, ny = 256, nz = 256
    real(rkind) :: Lx = 10.0d0, Ly = 2.0d0, Lz = 4.0d0
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k
    real(rkind) :: dx, dy, dz
    character(len=clen) :: inputfile, outputfile
    logical :: periodicx=.true., periodicy=.true., periodicz=.true.
    logical :: xmetric=.false., ymetric=.false., zmetric=.false.
    character(len=clen) :: derivative_x="cd10", derivative_y="cd10", derivative_z = "cd10"
    real(rkind) :: linf, l2norm, errcen, omega
    integer :: ii, kk, jj, x_bc(2), y_bc(2), z_bc(2), rank_debug = 0

    namelist /INPUT/ nx, ny, nz, Lx, Ly, Lz, prow, pcol, rank_debug, periodicx, periodicy, periodicz, omega

    call MPI_Init(ierr)

    write(inputfile,'(a)') 'input_test_multiblock_setup.dat'
    open(unit=123, file=inputfile, form='FORMATTED', iostat=ierr)
    read(unit=123, NML=INPUT)
    close(123)

    call decomp_2d_init(nx, ny, nz, prow, pcol)
    call get_decomp_info(gp)

    ! Initialize mesh arrays
    allocate(mesh(gp%ysz(1), gp%ysz(2), gp%ysz(3), 3) )
    allocate(xbuf(gp%xsz(1), gp%xsz(2), gp%xsz(3), 4) )
    allocate(zbuf(gp%zsz(1), gp%zsz(2), gp%zsz(3), 2) )
    allocate(mesh_curv(gp%ysz(1), gp%ysz(2), gp%ysz(3), 3) )
    allocate(   f(gp%ysz(1), gp%ysz(2), gp%ysz(3)))
    allocate(dfdx(gp%ysz(1), gp%ysz(2), gp%ysz(3)))
    allocate(dfdy(gp%ysz(1), gp%ysz(2), gp%ysz(3)))
    allocate(dfdz(gp%ysz(1), gp%ysz(2), gp%ysz(3)))
    allocate(fxex(gp%ysz(1), gp%ysz(2), gp%ysz(3)))
    allocate(fyex(gp%ysz(1), gp%ysz(2), gp%ysz(3)))
    allocate(fzex(gp%ysz(1), gp%ysz(2), gp%ysz(3)))

    x => mesh(:,:,:,1);  y => mesh(:,:,:,2);  z => mesh(:,:,:,3)
    xi => mesh_curv(:,:,:,1);  eta => mesh_curv(:,:,:,2);  zeta => mesh_curv(:,:,:,3)
    xtmp1 => xbuf(:,:,:,1);  xtmp2 => xbuf(:,:,:,2)
    ztmp1 => zbuf(:,:,:,1);  ztmp2 => zbuf(:,:,:,2)
    x_in_X => xbuf(:,:,:,3);  fxex_in_X => xbuf(:,:,:,4)

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

    allocate(mbtopology)
    call mbtopology%init(gp, mesh, inputfile, xbuf, zbuf)

    allocate(der)
    call der%init(gp, dx, dy, dz, periodicx, periodicy, periodicz, &
           derivative_x, derivative_y, derivative_z, x, y, z, xmetric, ymetric, zmetric, &
           .false., inputfile, xi, eta, zeta, xbuf, zbuf, mbtopology)
    f = sin((4*pi/Lx)*x) * cos((2*pi/Ly)*y) * sin((4*pi/Lz)*z)
    fxex = (4*pi/Lx) * cos((4*pi/Lx)*x) * cos((2*pi/Ly)*y) * sin((4*pi/Lz)*z)
    fyex = -(2*pi/Ly) * sin((4*pi/Lx)*x) * sin((2*pi/Ly)*y) * sin((4*pi/Lz)*z)
    fzex = (4*pi/Lz) * sin((4*pi/Lx)*x) * cos((2*pi/Ly)*y) * cos((4*pi/Lz)*z)

    fxex = fxex * mbtopology%mask

    call transpose_y_to_x(x,x_in_X,gp)
    call transpose_y_to_x(fxex,fxex_in_X,gp)

    x_bc(1) = 0; x_bc(2) = 0  ! immaterial for periodic
    y_bc(1) = 0; y_bc(2) = 0  ! immaterial for periodic
    z_bc(1) = 0; z_bc(2) = 0  ! immaterial for periodic

    ! x-derivative
    call transpose_y_to_x(f,xtmp1,gp)
    call der%ddx(xtmp1,xtmp2,x_bc(1),x_bc(2))
    call transpose_x_to_y(xtmp2,dfdx,gp)
    ! quantify the error in x-derivative
    if(nrank==rank_debug) then
      if(periodicx) then
        write(outputfile,'(a)') 'dfdx_per.dat'
      else
        write(outputfile,'(a)') 'dfdx_nonper.dat'
      endif
      call write_file_x(outputfile, x_in_X, xtmp1, xtmp2, fxex_in_X, gp%xsz(2)/2, gp%xsz(3)/2)
    endif
    linf = P_MAXVAL( MAXVAL(ABS(dfdx - fxex)))
    l2norm = SQRT(P_SUM(ABS(dfdx - fxex)**2)/(nx*ny*nz))
    ii = gp%ysz(1)*3/4; jj = gp%ysz(2)/2; kk = gp%ysz(3)/2
    if(nrank==rank_debug) then
        errcen =  ABS(dfdx(ii,jj,kk) - fxex(ii,jj,kk))
        if(periodicx) then
          write(outputfile,'(a)') 'xder_err_per.dat'
        else
          write(outputfile,'(a)') 'xder_err_nonper.dat'
        endif
        open(unit=123, file=outputfile, action='write',iostat=ierr,position='append')
        write(123,'(i5,1x,e19.12,1x,e19.12,1x,e19.12,1x)') nx, linf, l2norm, errcen
        close(123)
    endif

    call der%destroy()
    deallocate(der)

    call mbtopology%destroy()
    deallocate(mbtopology)

    nullify(x, y, z, xi, eta, zeta, xtmp1, xtmp2, ztmp1, ztmp2)
    nullify(x_in_X, fxex_in_X)
    deallocate(f, dfdx, dfdy, dfdz, fxex, fyex, fzex)
    deallocate(zbuf, xbuf, mesh_curv, mesh)
    call decomp_2d_finalize
    call MPI_Finalize(ierr)

end program

