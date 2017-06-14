program test_pcd06
    use mpi
    use kind_parameters, only: rkind
    use constants,       only: two, pi
    use exits,           only: GracefulExit
    use reductions,      only: P_MAXVAL
    use t3dMod,          only: t3d
    use pcd06stuff,      only: pcd06
    implicit none

    type(t3d) :: gp
    type(pcd06) :: xpcd06
    real(rkind), dimension(:,:,:), allocatable :: input
    integer :: nx = 16, ny = 16, nz = 16
    integer :: px = 4, py = 1, pz = 1
    integer :: i, j, k, ierr
    logical :: fail
    real(rkind) :: omega = 1._rkind, dx, dy, dz
    logical :: optimize = .false.

    call MPI_Init(ierr)

    if (.not. optimize) then
        gp = t3d(MPI_COMM_WORLD, nx, ny, nz, px, py, pz, [.TRUE., .TRUE., .TRUE.], .TRUE., fail)
        if (fail) call GracefulExit("t3d initialization failed",45)
    else
        gp = t3d(MPI_COMM_WORLD, nx, ny, nz, [.TRUE., .TRUE., .TRUE.], nghosts=[1,1,1])
    end if

    dx = two*pi/real(nx,rkind)
    dy = two*pi/real(ny,rkind)
    dz = two*pi/real(nz,rkind)

    ierr = xpcd06%init(gp, dx, .true., 0, 0)

    allocate(  input(gp%sz3D(1), gp%sz3D(2), gp%sz3D(3)) )

    do k=1,gp%sz3d(3)
        do j=1,gp%sz3d(2)  
            do i=1,gp%sz3d(1)
                input(i,j,k) = sin( omega * (gp%st3d(1) - 1 + i - 1) * dx ) + &
                               sin( omega * (gp%st3d(2) - 1 + j - 1) * dy ) + &
                               sin( omega * (gp%st3d(3) - 1 + k - 1) * dz )
            end do
        end do
    end do

    deallocate(  input )

    call MPI_Finalize(ierr)
contains

    subroutine print_array(a, aname)
        use kind_parameters, only: stdout
        real(rkind), dimension(:,:,:), intent(in) :: a
        character(len=*), intent(in) :: aname
        integer :: i, j, k

        call sleep(gp%rank3d)

        write(stdout,'(A,I0,A)') 'Rank ', gp%rank3d, ', Array: '//trim(aname)
        do k = 1,size(a,3)
            write(stdout,'(A,I0)') 'k = ', k
            do j = 1,size(a,2)
                write(stdout,*) ( a(i,j,k), i=1,size(a,1) )
            end do
            write(stdout,'(A)') ' '
        end do
        write(stdout,'(A)') '-----------------------------------'

    end subroutine

end program
