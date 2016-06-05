program test_KSprep
    use kind_parameters, only: rkind, clen
    use decomp_2d 
    use decomp_2d_io
    use spectralMod, only: spectral
    use mpi
    use constants
    use kspreprocessing, only: ksprep
    
    implicit none

    real(rkind), dimension(:,:,:), allocatable :: uC, vC, wE, uE, vE
    type(decomp_info) :: gpE, gpC
    real(rkind), dimension(:,:,:), allocatable :: buffyC, buffyE, buffzC, buffzE
    type(spectral) :: spectE
    integer :: nx = 512, ny = 256, nz = 128
    character(len=clen) :: outputdir = "/home/aditya90/Codes/PadeOps/data/LES2KS", fname
    integer :: ierr 
    real(rkind) :: dx, dy, dz
    type(ksprep) :: les2ks
    integer, dimension(4) :: planes2dumpC = [2,8,16,32]
    integer, dimension(4) :: planes2dumpF = [4,16,32,64]

    call MPI_Init(ierr)
    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(gpC)
    call decomp_info_init(nx,ny,nz+1, gpE)

    dx = four*pi/real(nx,rkind)
    dy = two*pi/real(ny,rkind)
    dz = one/real(nz,rkind)

    call spectE%init("x", nx, ny, nz+1, dx, dy, dz, "four", "2/3rd", 2, .false.)

    allocate(uC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(vC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(buffyC(gpC%ysz(1),gpC%ysz(2),gpC%ysz(3)))
    allocate(buffzC(gpC%zsz(1),gpC%zsz(2),gpC%zsz(3)))
    allocate(buffyE(gpE%ysz(1),gpE%ysz(2),gpE%ysz(3)))
    allocate(buffzE(gpE%zsz(1),gpE%zsz(2),gpE%zsz(3)))
    allocate(uE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
    allocate(vE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
    allocate(wE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
    fname = "/home/aditya90/Codes/PadeOps/data/LES2KS/RESTART_Run20_u.096000"
    call decomp_2d_read_one(1,uC,trim(fname), gpC)
    fname = "/home/aditya90/Codes/PadeOps/data/LES2KS/RESTART_Run20_v.096000"
    call decomp_2d_read_one(1,vC,trim(fname), gpC)
    fname = "/home/aditya90/Codes/PadeOps/data/LES2KS/RESTART_Run20_w.096000"
    call decomp_2d_read_one(1,wE,trim(fname), gpE)

    call transpose_x_to_y(uC,buffyC,gpC)
    call transpose_y_to_z(buffyC,buffzC,gpC)
    buffzE(:,:,2:nz) = half*(buffzC(:,:,1:nz-1) + buffzC(:,:,2:nz))
    buffzE(:,:,nz+1) = buffzE(:,:,nz)
    buffzE(:,:,1) = two*buffzC(:,:,1) - buffzE(:,:,2)
    call transpose_z_to_y(buffzE,buffyE,gpE)
    call transpose_y_to_x(buffyE,uE,gpE)

    call transpose_x_to_y(vC,buffyC,gpC)
    call transpose_y_to_z(buffyC,buffzC,gpC)
    buffzE(:,:,2:nz) = half*(buffzC(:,:,1:nz-1) + buffzC(:,:,2:nz))
    buffzE(:,:,nz+1) = buffzE(:,:,nz)
    buffzE(:,:,1) = two*buffzC(:,:,1) - buffzE(:,:,2)
    call transpose_z_to_y(buffzE,buffyE,gpE)
    call transpose_y_to_x(buffyE,vE,gpE)


    call les2ks%init(nx,ny,nz,spectE, gpE, outputdir,10, planes2dumpC, planes2dumpF)

    call les2ks%dumpKSprep(uE,vE,wE,0)
    call les2ks%destroy()



end program
