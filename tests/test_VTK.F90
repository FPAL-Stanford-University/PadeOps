program test_VTK
    use mpi
    use kind_parameters, only : rkind
    use decomp_2d, only: decomp_info, get_decomp_info, decomp_2d_init, decomp_2d_finalize, &
                    transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y,&
                    update_halo, nrank, nproc
    use Lib_VTK_IO
    use IR_Precision
    use constants, only: pi, two

    implicit none

    real(rkind), dimension(:,:,:), allocatable :: x, y, z, f, g
    real(rkind), dimension(:,:,:), allocatable :: tmp1,tmp2,tmp3
    type(decomp_info) :: gp

    logical :: periodicx = .TRUE., periodicy = .TRUE., periodicz = .TRUE.

    integer :: nx = 128, ny = 128, nz = 128
    integer :: nxp, nyp, nzp
    integer :: nx1, ny1, nz1
    integer :: nx2, ny2, nz2
    integer :: nn
    integer :: prow = 0, pcol = 0
    integer :: ierr, i, j, k

    integer :: E_IO
    
    real(rkind) :: dx, dy, dz

    double precision :: t0, t1

    integer, dimension(MPI_STATUS_SIZE) :: mpistatus

    call MPI_Init(ierr)
    
    call decomp_2d_init(nx, ny, nz, prow, pcol, [periodicx, periodicy, periodicz])
    call get_decomp_info(gp)
   
    nxp = gp%ysz(1); nyp = gp%ysz(2); nzp = gp%ysz(3)

    ! Initialize input data (y decomposition) 
    allocate( x ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( y ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( z ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( f ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    allocate( g ( gp%ysz(1), gp%ysz(2), gp%ysz(3) ) )
    
    dx = two*pi/nx
    dy = two*pi/ny
    dz = two*pi/nz

    ! Generate input data
    do k = 1,nzp
        do j = 1,nyp
            do i = 1,nxp
                x(i,j,k) = real(gp%yst(1) - 1 + i - 1, rkind)*dx
                y(i,j,k) = real(gp%yst(2) - 1 + j - 1, rkind)*dy 
                z(i,j,k) = real(gp%yst(3) - 1 + k - 1, rkind)*dz
            end do 
        end do 
    end do 

    ! Create the fields
    f = sin(x)*sin(y)*cos(z)
    g = real(nrank,rkind)
    
    if (nrank == 0) print*, "Created coordinates and fields"

    nx1 = gp%yst(1)
    nx2 = gp%yen(1)+1
    ny1 = gp%yst(2)
    ny2 = gp%yen(2)      ! no +1 here since we're in the y decomposition
    nz1 = gp%yst(3)
    nz2 = gp%yen(3)+1

    ! No overlap point for boundary processors
    if ( (gp%yen(1) == nx) ) then
        nx2 = nx
    end if
    if ( (gp%yen(3) == nz) ) then
        nz2 = nz
    end if

    nn = (nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1)
    
    E_IO = VTK_INI_XML_WRITE(fformat='binary', filename='XML_MPI_f'//trim(strz(3,nrank))//'.vts', &
                       mesh_topology='StructuredGrid', nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2)
    
    ! Halo update for x, y and z
    call update_halo(x,tmp1,1,gp,.FALSE.)
    call update_halo(y,tmp2,1,gp,.FALSE.)
    call update_halo(z,tmp3,1,gp,.FALSE.)

    if (nrank == 0) print*, "Finished halo comm."

    E_IO = VTK_GEO_XML_WRITE(nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2,NN=nn,&
                             X=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1),          &
                             Y=tmp2(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1),          &
                             Z=tmp3(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1)           )

    deallocate(tmp1,tmp2,tmp3)

    E_IO = VTK_DAT_XML(var_location='node',var_block_action='open')
    
    call update_halo(f,tmp1,1,gp,.FALSE.)
    E_IO = VTK_VAR_XML(NC_NN=nn,varname='field',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
    deallocate(tmp1)

    call update_halo(g,tmp1,1,gp,.FALSE.)
    E_IO = VTK_VAR_XML(NC_NN=nn,varname='rank',var=tmp1(1:nx2-nx1+1,1:ny2-ny1+1,1:nz2-nz1+1))
    deallocate(tmp1)

    E_IO = VTK_DAT_XML(var_location='node',var_block_action='close')
    E_IO = VTK_GEO_XML_WRITE()
    
    E_IO = VTK_END_XML()

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (nrank == 0) print*, "Finished writing individual data files."

    if (nrank == 0) then
        ! First process saves also the composite .pvts file
        print*, "Now writing encapsulating pvts file"
        E_IO = PVTK_INI_XML(filename = 'XML_MPI.pvts', mesh_topology = 'PStructuredGrid',&
                            nx1=1, nx2=nx, ny1=1, ny2=ny, nz1=1, nz2=nz, tp='Float64')
        do i=0,nproc-1
            if (i .NE. 0) then
                call MPI_RECV(nx1,1,MPI_INTEGER,i,i        ,MPI_COMM_WORLD,mpistatus,ierr)
                call MPI_RECV(nx2,1,MPI_INTEGER,i,i+  nproc,MPI_COMM_WORLD,mpistatus,ierr)
                call MPI_RECV(ny1,1,MPI_INTEGER,i,i+2*nproc,MPI_COMM_WORLD,mpistatus,ierr)
                call MPI_RECV(ny2,1,MPI_INTEGER,i,i+3*nproc,MPI_COMM_WORLD,mpistatus,ierr)
                call MPI_RECV(nz1,1,MPI_INTEGER,i,i+4*nproc,MPI_COMM_WORLD,mpistatus,ierr)
                call MPI_RECV(nz2,1,MPI_INTEGER,i,i+5*nproc,MPI_COMM_WORLD,mpistatus,ierr)
            end if
            E_IO = PVTK_GEO_XML(nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2,&
                                source='XML_MPI_f'//trim(strz(3,i))//'.vts')
        enddo

        E_IO = PVTK_DAT_XML(var_location='node',var_block_action='open')
        E_IO = PVTK_VAR_XML(varname='field',tp='Float64')
         
        E_IO = PVTK_VAR_XML(varname='rank',tp='Float64')
        E_IO = PVTK_DAT_XML(var_location='node',var_block_action='close')
        
        E_IO = PVTK_END_XML()
    else
        call MPI_SEND(nx1,1,MPI_INTEGER,0,nrank        ,MPI_COMM_WORLD,ierr)
        call MPI_SEND(nx2,1,MPI_INTEGER,0,nrank+  nproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(ny1,1,MPI_INTEGER,0,nrank+2*nproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(ny2,1,MPI_INTEGER,0,nrank+3*nproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(nz1,1,MPI_INTEGER,0,nrank+4*nproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(nz2,1,MPI_INTEGER,0,nrank+5*nproc,MPI_COMM_WORLD,ierr)
    end if
  
    deallocate(x, y, z, f, g)
    call decomp_2d_finalize
    call MPI_Finalize(ierr)

end program 
