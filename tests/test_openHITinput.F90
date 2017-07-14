program test_openHITinput

    use decomp_2d
    use decomp_2d_io
    use kind_parameters, only: rkind
    use mpi 

    implicit none


    integer :: nx = 256, ny = 256, nz = 256
    type(decomp_info) :: gp
    integer :: ierr

    real(rkind), dimension(:,:,:), allocatable :: u 

    
    call mpi_init(ierr)
    call decomp_2d_init(nx,ny,nz,0,0)
    call get_decomp_info(gp)

    allocate(u(gp%xsz(1), gp%xsz(2), gp%xsz(3)))

    call decomp_2d_read_one(1,u,'/home/aditya90/Codes/PadeOps/data/HIT_testing/u_HIT_init_256.dat')

    print*, u(3,4,2)
    
    call mpi_finalize(ierr)  
    
end program 

