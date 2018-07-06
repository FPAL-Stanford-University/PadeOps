program test_openHITinput

    use decomp_2d
    use decomp_2d_io
    use kind_parameters, only: rkind
    use basic_io, only: read_2D_ascii 
    use mpi 

    implicit none


    integer :: nx = 32, ny = 32, nz = 32
    type(decomp_info) :: gp
    real(rkind), dimension(:,:), allocatable :: dat
    integer :: ierr

    real(rkind), dimension(:,:,:), allocatable :: u, v, w 

    
    call mpi_init(ierr)
    call decomp_2d_init(nx,ny,nz,1,1)
    call get_decomp_info(gp)

    allocate(dat(nx*ny*nz,3))
    allocate(u(gp%xsz(1), gp%xsz(2), gp%xsz(3)))
    allocate(v(gp%xsz(1), gp%xsz(2), gp%xsz(3)))
    allocate(w(gp%xsz(1), gp%xsz(2), gp%xsz(3)))

    call read_2D_ascii(dat,"/home1/04076/tg833754/HITinit32.dat")

    u = reshape(dat(:,1),[nx,ny,nz])
    v = reshape(dat(:,2),[nx,ny,nz])
    w = reshape(dat(:,3),[nx,ny,nz])


    call decomp_2d_write_one(1,u,"/home1/04076/tg833754/Run01_uVel_t000000.out",gp) 
    call decomp_2d_write_one(1,v,"/home1/04076/tg833754/Run01_vVel_t000000.out",gp) 
    call decomp_2d_write_one(1,w,"/home1/04076/tg833754/Run01_wVel_t000000.out",gp) 
    call mpi_finalize(ierr)  
    
end program 

