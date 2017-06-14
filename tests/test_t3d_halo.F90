program test_t3d_halo
    use mpi
    use kind_parameters, only: rkind
    use constants,       only: two, pi
    use exits,           only: GracefulExit
    use reductions,      only: P_MAXVAL
    use t3dMod,          only: t3d, square_factor, roundrobin_split
    use cd10stuff,       only: cd10
    implicit none

    type(t3d) :: gp
    real(rkind), dimension(:,:,:), allocatable :: input
    integer :: nx = 6, ny = 4, nz = 4
    integer :: px = 3, py = 1, pz = 1
    integer :: i, j, k, ierr
    real(rkind) :: dx, dy, dz
    logical :: fail
    logical :: optimize = .false.

    call MPI_Init(ierr)

    if (.not. optimize) then
        gp = t3d(MPI_COMM_WORLD, nx, ny, nz, px, py, pz, [.TRUE., .TRUE., .TRUE.], .TRUE., fail, nghosts=[1,1,1])
        if (fail) call GracefulExit("t3d initialization failed",45)
    else
        gp = t3d(MPI_COMM_WORLD, nx, ny, nz, [.TRUE., .TRUE., .TRUE.], nghosts=[1,1,1])
    end if

    allocate( input(gp%st3Dg(1):gp%en3Dg(1), gp%st3D(2):gp%en3D(2), gp%st3D(3):gp%en3D(3)) )
    input = 0._rkind

    dx = two*pi/real(nx,rkind)
    dy = two*pi/real(ny,rkind)
    dz = two*pi/real(nz,rkind)

    do k=gp%st3D(3),gp%en3D(3)
        do j=gp%st3D(2),gp%en3D(2)  
            do i=gp%st3D(1),gp%en3D(1)
                input(i,j,k) = (i-1) + (j-1)*nx + (k-1)*nx*ny
            end do
        end do
    end do

    call print_array(input, "input, before halo exchange")

    call mpi_barrier(MPI_COMM_WORLD, ierr)

    call gp%fill_halo_x( input )

    call print_array(input, "input, after halo exchange")

    deallocate( input )

    call MPI_Finalize(ierr)
contains

    subroutine print_array(a, aname)
        use kind_parameters, only: stdout
        real(rkind), dimension(:,:,:), intent(in) :: a
        character(len=*), intent(in) :: aname
        integer :: i, j, k

        call sleep(gp%rank3d)

        write(stdout,'(A,I0,A)') 'Rank ', gp%rank3d, ', Array: '//trim(aname)
        do k = lbound(a,3),ubound(a,3)
            write(stdout,'(A,I0)') 'k = ', k
            do j = lbound(a,2),ubound(a,2)
                write(stdout,*) ( a(i,j,k), i=lbound(a,1),ubound(a,1) )
            end do
            write(stdout,'(A)') ' '
        end do
        write(stdout,'(A)') '-----------------------------------'

    end subroutine

end program
