program test_t3d
    use mpi
    use kind_parameters, only: rkind
    use exits,           only: GracefulExit
    use t3dMod,          only: t3d, square_factor, roundrobin_split
    implicit none

    type(t3d) :: gp
    real(rkind), dimension(:,:,:), allocatable :: input, output, ninput
    integer :: nx = 5, ny = 5, nz = 4
    integer :: px = 4, py = 1, pz = 1
    integer :: ierr
    logical :: fail

    call MPI_Init(ierr)

    gp = t3d(MPI_COMM_WORLD, nx, ny, nz, px, py, pz, [.TRUE., .TRUE., .TRUE.], .TRUE., fail)
    if (fail) call GracefulExit("t3d initialization failed",45)

    allocate(  input(gp%sz3D(1), gp%sz3D(2), gp%sz3D(3)) )
    allocate( ninput(gp%sz3D(1), gp%sz3D(2), gp%sz3D(3)) )
    allocate( output(gp%szX (1), gp%szX (2), gp%szX (3)) )

    ! call gp%print_summary()

    call gp%transpose_3D_to_x(input,output)
    call gp%transpose_x_to_3D(output,ninput)
   
    print*, gp%rank3D, maxval(abs(ninput-input))

    deallocate(  input )
    deallocate( ninput )
    deallocate( output )

    call MPI_Finalize(ierr)

end program
