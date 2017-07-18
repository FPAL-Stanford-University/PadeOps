program test_t3d_halo
    use mpi
    use kind_parameters, only: rkind
    use constants,       only: eps, two, pi
    use exits,           only: GracefulExit
    use t3dMod,          only: t3d
    implicit none

    type(t3d) :: gp
    real(rkind), dimension(:,:,:), allocatable :: input
    integer :: nx = 16, ny = 16, nz = 16
    integer :: px = 2, py = 2, pz = 2
    integer, dimension(3) :: nghosts = [3,2,1]
    logical, dimension(3) :: periodic = [.false., .true., .true.]
    integer :: i, j, k, ierr
    integer :: ii, jj, kk
    real(rkind) :: dx, dy, dz
    logical :: fail
    logical :: optimize = .false.
    logical mycorrect, correct

    call MPI_Init(ierr)

    if (.not. optimize) then
        gp = t3d(MPI_COMM_WORLD, nx, ny, nz, px, py, pz, periodic, .TRUE., fail, nghosts=nghosts)
        if (fail) call GracefulExit("t3d initialization failed",45)
    else
        gp = t3d(MPI_COMM_WORLD, nx, ny, nz, periodic, nghosts=nghosts)
    end if

    dx = two*pi/real(nx,rkind)
    dy = two*pi/real(ny,rkind)
    dz = two*pi/real(nz,rkind)

    if (gp%rank3D == 0) then
        print *, "st3D  = ", gp%st3D
        print *, "en3D  = ", gp%en3D
        print *, "st3Dg = ", gp%st3Dg
        print *, "en3Dg = ", gp%en3Dg
    end if

    !!!! ============================== !!!!
    !!!! X-direction halo communication
    allocate( input(gp%st3Dg(1):gp%en3Dg(1), gp%st3D(2):gp%en3D(2), gp%st3D(3):gp%en3D(3)) )
    input = 0._rkind

    do k=gp%st3D(3),gp%en3D(3)
        do j=gp%st3D(2),gp%en3D(2)  
            do i=gp%st3D(1),gp%en3D(1)
                input(i,j,k) = (i-1) + (j-1)*nx + (k-1)*nx*ny
            end do
        end do
    end do

    call gp%fill_halo_x( input )

    mycorrect = .true.
    do kk=gp%st3D(3),gp%en3D(3)
        do jj=gp%st3D(2),gp%en3D(2)  
            do ii=gp%st3Dg(1),gp%en3Dg(1)
                i = mod(ii-1+nx,nx)+1
                j = mod(jj-1+ny,ny)+1
                k = mod(kk-1+nz,nz)+1
                if ( abs(input(ii,jj,kk) - ( (i-1) + (j-1)*nx + (k-1)*nx*ny ) ) > eps ) then
                    mycorrect = .false.
                end if
            end do
        end do
    end do
    call MPI_Reduce(mycorrect, correct, 1, MPI_LOGICAL, MPI_LAND, 0, MPI_COMM_WORLD, ierr)

    if (gp%rank3D == 0) then
        if (correct .eqv. .true.) then
            print '(A)', "Halo cells in X communicated correctly! :)"
        else
            print '(A)', "ERROR: Halo cells in X not communicated correctly!"
        end if
    end if
    deallocate( input )
    !!!! ============================== !!!!

    !!!! ============================== !!!!
    !!!! Y-direction halo communication
    allocate( input(gp%st3D(1):gp%en3D(1), gp%st3Dg(2):gp%en3Dg(2), gp%st3D(3):gp%en3D(3)) )
    input = 0._rkind

    do k=gp%st3D(3),gp%en3D(3)
        do j=gp%st3D(2),gp%en3D(2)  
            do i=gp%st3D(1),gp%en3D(1)
                input(i,j,k) = (i-1) + (j-1)*nx + (k-1)*nx*ny
            end do
        end do
    end do

    call gp%fill_halo_y( input )

    mycorrect = .true.
    do kk=gp%st3D(3),gp%en3D(3)
        do jj=gp%st3Dg(2),gp%en3Dg(2)  
            do ii=gp%st3D(1),gp%en3D(1)
                i = mod(ii-1+nx,nx)+1
                j = mod(jj-1+ny,ny)+1
                k = mod(kk-1+nz,nz)+1
                if ( abs(input(ii,jj,kk) - ( (i-1) + (j-1)*nx + (k-1)*nx*ny ) ) > eps ) then
                    mycorrect = .false.
                end if
            end do
        end do
    end do
    call MPI_Reduce(mycorrect, correct, 1, MPI_LOGICAL, MPI_LAND, 0, MPI_COMM_WORLD, ierr)

    if (gp%rank3D == 0) then
        if (correct .eqv. .true.) then
            print '(A)', "Halo cells in Y communicated correctly! :)"
        else
            print '(A)', "ERROR: Halo cells in Y not communicated correctly!"
        end if
    end if
    deallocate( input )
    !!!! ============================== !!!!

    !!!! ============================== !!!!
    !!!! Z-direction halo communication
    allocate( input(gp%st3D(1):gp%en3D(1), gp%st3D(2):gp%en3D(2), gp%st3Dg(3):gp%en3Dg(3)) )
    input = 0._rkind

    do k=gp%st3D(3),gp%en3D(3)
        do j=gp%st3D(2),gp%en3D(2)  
            do i=gp%st3D(1),gp%en3D(1)
                input(i,j,k) = (i-1) + (j-1)*nx + (k-1)*nx*ny
            end do
        end do
    end do

    call gp%fill_halo_z( input )

    mycorrect = .true.
    do kk=gp%st3Dg(3),gp%en3Dg(3)
        do jj=gp%st3D(2),gp%en3D(2)  
            do ii=gp%st3D(1),gp%en3D(1)
                i = mod(ii-1+nx,nx)+1
                j = mod(jj-1+ny,ny)+1
                k = mod(kk-1+nz,nz)+1
                if ( abs(input(ii,jj,kk) - ( (i-1) + (j-1)*nx + (k-1)*nx*ny ) ) > eps ) then
                    mycorrect = .false.
                end if
            end do
        end do
    end do
    call MPI_Reduce(mycorrect, correct, 1, MPI_LOGICAL, MPI_LAND, 0, MPI_COMM_WORLD, ierr)

    if (gp%rank3D == 0) then
        if (correct .eqv. .true.) then
            print '(A)', "Halo cells in Z communicated correctly! :)"
        else
            print '(A)', "ERROR: Halo cells in Z not communicated correctly!"
        end if
    end if
    deallocate( input )
    !!!! ============================== !!!!

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
