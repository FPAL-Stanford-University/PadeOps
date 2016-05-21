program upsapleFields
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa 
    use gridtools, only: alloc_buffs, destroy_buffs
    use decomp_2d
    use exits, only: GracefulExit, message
    use mpi 
    use timer, only: tic, toc
    use decomp_2d_io
    implicit none

    integer :: nx, ny, nz
    real(rkind), dimension(:,:,:) :: f, fup, fxup, fyup, fzup
    type(decomp_info) :: gpC, gpE, gpCL, gpEL
    namelist /INPUT/ nx, ny, nz, inputdir, outputdir, inputFile_TID, inputFile_RID, &
                     outputFile_TID, outputFile_RID 

    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file

    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(gpC)
    call decomp_info_init(nx,ny,nz+1,gpE)
    call decomp_info_init(2*nx,2*ny,2*nz,gpCL)
    call decomp_info_init(2*nx,2*ny,2*nz+1,gpEL)
    
    call allocate(f(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    call allocate(fL(gpCL%xsz(1),gpCL%xsz(2),gpCL%xsz(3)))





end program 
