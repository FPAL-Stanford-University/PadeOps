! subroutine print_array(a, aname)
!     use kind_parameters, only: stdout
!     real(rkind), dimension(:,:,:), intent(in) :: a
!     character(len=*), intent(in) :: aname
!     integer :: i, j, k
! 
!     call sleep(gp%rank3d)
! 
!     write(stdout,'(A,I0,A)') 'Rank ', gp%rank3d, ', Array: '//trim(aname)
!     do k = 1,size(a,3)
!         write(stdout,'(A,I0)') 'k = ', k
!         do j = 1,size(a,2)
!             write(stdout,*) ( a(i,j,k), i=1,size(a,1) )
!         end do
!         write(stdout,'(A)') ' '
!     end do
!     write(stdout,'(A)') '-----------------------------------'
! 
! end subroutine

subroutine dowork(xinput, xoutput)
    use kind_parameters, only: rkind
    real(rkind), dimension(:,:,:), intent(in)  :: xinput
    real(rkind), dimension(:,:,:), intent(out) :: xoutput

    xoutput = exp(xinput)
end subroutine

program test_t3d_NB
    use mpi
    use kind_parameters, only: rkind
    use constants,       only: two, pi
    use exits,           only: GracefulExit
    use reductions,      only: P_MAXVAL
    use t3dMod,          only: t3d
    use cd10stuff,       only: cd10
    implicit none

    type(t3d) :: gp
    type(cd10) :: xcd10, ycd10, zcd10
    real(rkind), dimension(:,:,:), allocatable :: input, output, ninput, der, exder, eyder, ezder
    real(rkind), dimension(:), allocatable :: buffx, buffy, buffz, buff3d 
    integer :: nx = 1024, ny = 1024, nz = 1024
    integer :: px = 4, py = 4, pz = 4
    integer :: i, j, k, ierr
    logical :: fail
    real(rkind) :: dererr, omega = 1._rkind, dx, dy, dz
    real(rkind) :: start, tx, txd, t3dx, ty, tyd, t3dy, tz, tzd, t3dz
    logical :: optimize = .false.
    integer :: request
    integer, dimension(mpi_status_size) :: status

    interface 
        subroutine dowork(xinput, xoutput)
            import :: rkind
            real(rkind), dimension(:,:,:), intent(in)  :: xinput
            real(rkind), dimension(:,:,:), intent(out) :: xoutput
        end subroutine
    end interface

    call MPI_Init(ierr)

    if (.not. optimize) then
        gp = t3d(MPI_COMM_WORLD, nx, ny, nz, px, py, pz, [.TRUE., .TRUE., .TRUE.], .TRUE., fail)
        if (fail) call GracefulExit("t3d initialization failed",45)
    else
        gp = t3d(MPI_COMM_WORLD, nx, ny, nz, [.TRUE., .TRUE., .TRUE.])
    end if

    dx = two*pi/real(nx,rkind)
    dy = two*pi/real(ny,rkind)
    dz = two*pi/real(nz,rkind)
    ierr = xcd10%init( nx, dx, .TRUE., 0, 0 )
    ierr = ycd10%init( ny, dy, .TRUE., 0, 0 )
    ierr = zcd10%init( nz, dz, .TRUE., 0, 0 )

    allocate(  input(gp%sz3D(1), gp%sz3D(2), gp%sz3D(3)) )
    allocate(  exder(gp%sz3D(1), gp%sz3D(2), gp%sz3D(3)) )
    allocate(  eyder(gp%sz3D(1), gp%sz3D(2), gp%sz3D(3)) )
    allocate(  ezder(gp%sz3D(1), gp%sz3D(2), gp%sz3D(3)) )
    allocate( ninput(gp%sz3D(1), gp%sz3D(2), gp%sz3D(3)) )
    allocate( output(gp%szX (1), gp%szX (2), gp%szX (3)) )
    allocate(   der(gp%szX (1), gp%szX (2), gp%szX (3)) )

    allocate( buff3d(gp%sz3D(1) * gp%sz3D(2) * gp%sz3D(3)) )
    allocate( buffx (gp%szX (1) * gp%szX (2) * gp%szX (3)) )
    allocate( buffy (gp%szY (1) * gp%szY (2) * gp%szY (3)) )
    allocate( buffz (gp%szZ (1) * gp%szZ (2) * gp%szZ (3)) )

    ! call gp%print_summary()

    do k=1,gp%sz3d(3)
        do j=1,gp%sz3d(2)  
            do i=1,gp%sz3d(1)
                ! input(i,j,k) = (gp%st3d(1) - 1 + i - 1) + nx*(gp%st3d(2) - 1 + j - 1) + nx*ny*(gp%st3d(3) - 1 + k - 1)
                input(i,j,k) = sin( omega * (gp%st3d(1) - 1 + i - 1) * dx ) + &
                               sin( omega * (gp%st3d(2) - 1 + j - 1) * dy ) + &
                               sin( omega * (gp%st3d(3) - 1 + k - 1) * dz )
                exder(i,j,k) = omega * cos( omega * (gp%st3d(1) - 1 + i - 1) * dx )
                eyder(i,j,k) = omega * cos( omega * (gp%st3d(2) - 1 + j - 1) * dy )
                ezder(i,j,k) = omega * cos( omega * (gp%st3d(3) - 1 + k - 1) * dz )
            end do
        end do
    end do

    ! Warm-up dowork
    call dowork(input, eyder)
    call dowork(exder, ezder)
    
    start = gp%time(barrier=.false.)
    call dowork(input, eyder)
    tx = gp%time(start,reduce=.false.); if (gp%rank3d == 0) print*, "Time for dowork = ", tx

    start = gp%time(barrier=.true.)
    call gp%transpose_3D_to_x( input, output )
    tx = gp%time(start,reduce=.true.); if (gp%rank3d == 0) print*, "Time for blocking 3D to X transpose = ", tx

    ! X transpose and derivative
    start = gp%time(barrier=.true.)
    call gp%initiate_transpose_3D_to_x( input, buff3d, buffx, request)
    call dowork(exder, ezder)
    call gp%wait_transpose_3D_to_x( output, buffx, request, status)
    tx = gp%time(start,reduce=.true.); if (gp%rank3d == 0) print*, "Time for 3D to X transpose = ", tx

    start = gp%time(barrier=.true.)
    call xcd10%dd1(output, der, size(der,2), size(der,3))
    txd = gp%time(start,reduce=.true.)
    if (gp%rank3d == 0) print*, "Time for X derivative      = ", txd
    
    start = gp%time(barrier=.true.)
    call gp%transpose_x_to_3D(der,ninput)
    t3dx = gp%time(start,reduce=.true.)
    if (gp%rank3d == 0) print*, "Time for X to 3D transpose = ", t3dx

    dererr = P_MAXVAL(maxval(abs(ninput-exder)))
    if(gp%rank3d == 0) print*, " >>>> Error in x derivative = ", dererr

    deallocate(output, der);
    allocate( output(gp%szY (1), gp%szY (2), gp%szY (3)) )
    allocate(    der(gp%szY (1), gp%szY (2), gp%szY (3)) )

    ! ! Y transpose and derivative
    ! start = gp%time(barrier=.true.)
    ! call gp%transpose_3D_to_y(input,output)
    ! ty = gp%time(start,reduce=.false.)
    ! if (gp%rank3d == 0) print*, "Time for 3D to Y transpose = ", ty

    ! start = gp%time(barrier=.false.)
    ! call ycd10%dd2(output, der, size(der,1), size(der,3))
    ! tyd = gp%time(start,reduce=.false.)
    ! if (gp%rank3d == 0) print*, "Time for Y derivative      = ", tyd
    ! 
    ! start = gp%time(barrier=.false.)
    ! call gp%transpose_y_to_3D(der,ninput)
    ! t3dy = gp%time(start,reduce=.false.)
    ! if (gp%rank3d == 0) print*, "Time for X to 3D transpose = ", t3dy

    ! dererr = P_MAXVAL(maxval(abs(ninput-eyder)))
    ! if(gp%rank3d == 0) print*, ">>>> Error in y derivative  = ", dererr

    ! deallocate(output, der);
    ! allocate( output(gp%szZ (1), gp%szZ (2), gp%szZ (3)) )
    ! allocate(    der(gp%szZ (1), gp%szZ (2), gp%szZ (3)) )

    ! ! Z transpose and derivative
    ! start = gp%time(barrier=.true.)
    ! call gp%transpose_3D_to_z(input,output)
    ! tz = gp%time(start,reduce=.false.)
    ! if (gp%rank3d == 0) print*, "Time for 3D to Z transpose = ", tz

    ! start = gp%time(barrier=.false.)
    ! call zcd10%dd3(output, der, size(der,1), size(der,2))
    ! tzd = gp%time(start,reduce=.false.)
    ! if (gp%rank3d == 0) print*, "Time for Z derivative      = ", tzd
    ! 
    ! start = gp%time(barrier=.false.)
    ! call gp%transpose_z_to_3D(der,ninput)
    ! t3dz = gp%time(start,reduce=.false.)
    ! if (gp%rank3d == 0) print*, "Time for Z to 3D transpose = ", t3dz

    ! dererr = P_MAXVAL(maxval(abs(ninput-ezder)))
    ! if(gp%rank3d == 0) print*, ">>>> Error in z derivative  = ", dererr

    deallocate(  input )
    deallocate( ninput )
    deallocate( output )
    deallocate(   der  )
    deallocate(  exder )
    deallocate(  eyder )
    deallocate(  ezder )

    deallocate( buff3d )
    deallocate( buffx  )
    deallocate( buffy  )
    deallocate( buffz  )

    call MPI_Finalize(ierr)

end program
