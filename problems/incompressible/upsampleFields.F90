module RoutinesUpsampling
    use kind_parameters, only: rkind
    
    implicit none
    contains
    subroutine upsampleX(arrIn, arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: i, j, k, nx

        nx = size(arrIn,1)
        do k = 1,size(arrIn,3)
            do j = 1,size(arrIn,2)
                do i = 1,nx
                    arrOut(2*i-1,j,k) = arrIn(i,j,k)
                end do 
            end do 
        end do 

        do k = 1,size(arrIn,3)
            do j = 1,size(arrIn,2)
                do i = 2,2*nx-2,2
                    arrOut(i,j,k) = 0.5d0*(arrOut(i - 1,j,k) + arrOut(i + 1,j,k))
                end do 
                arrOut(2*nx,j,k) = 0.5d0*(arrOut(2*nx-1,j,k) + arrOut(1,j,k))
            end do 
        end do 

    end subroutine

    subroutine upsampleY(arrIn, arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: j, k, ny

        ny = size(arrIn,2)
        do k = 1,size(arrIn,3)
            do j = 1,ny
                    arrOut(:,2*j-1,k) = arrIn(:,j,k)
            end do 
        end do 

        do k = 1,size(arrIn,3)
            do j = 2,2*ny-2,2
                    arrOut(:,j,k) = 0.5d0*(arrOut(:,j - 1,k) + arrOut(:,j + 1,k))
            end do 
            arrOut(:,2*ny,k) = 0.5d0*(arrOut(:,2*ny-1,k) + arrOut(:,1,k))
        end do 

    end subroutine

    subroutine upsampleZ_cells(arrIn,arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: k, nz, idx
        real(rkind), parameter :: ratEven = 0.25d0, ratOdd = 0.75d0

        nz = size(arrIn,3)
        do k = 1,nz-1
            arrOut(:,:,2*k+1) = arrIn(:,:,k) + (arrIn(:,:,k+1) - arrIn(:,:,k))*ratOdd 
        end do 
        
        idx = 1
        do k = 2,2*nz-1,2
            arrOut(:,:,k) = arrIn(:,:,idx) + (arrIn(:,:,idx + 1) - arrIn(:,:,idx))*ratEven
            idx = idx + 1
        end do
   
        arrOut(:,:,1) = 2.d0*arrOut(:,:,2) - arrOut(:,:,3)
        arrOut(:,:,2*nz) = arrOut(:,:,2*nz - 1) 

    end subroutine


    subroutine upsampleZ_edges(arrIn,arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: k, nz

        nz = size(arrIn,3) - 1

        do k = 1,nz
            arrOut(:,:,2*k-1) = arrIn(:,:,k)
        end do 
        arrOut(:,:,2*nz+1) = arrOut(:,:,nz+1)

        do k = 1,nz
            arrOut(:,:,2*k) = 0.5d0*(arrOut(:,:,2*k-1) + arrOut(:,:,2*k+1))
        end do 

        arrOut(:,:,2*nz+1) = 0.d0
        arrOut(:,:,1) = 0.d0
    end subroutine

end module



program upsampleFields
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth, pi, kappa 
    use gridtools, only: alloc_buffs, destroy_buffs
    use decomp_2d
    use exits, only: GracefulExit, message
    use mpi 
    use timer, only: tic, toc
    use decomp_2d_io
    use RoutinesUpsampling
   
    implicit none

    character(len=clen) :: inputfile 
    character(len=clen) :: outputdir, inputdir
    integer :: ioUnit, nx, ny, nz, ierr, inputFile_TID, inputFile_RID, outputFile_TID, outputFile_RID
    logical :: upsampleInZ = .false., isStratified = .false. 
    type(decomp_info) :: gpC, gpE, gpC_upX, gpC_upXY, gpC_upXYZ, gpE_upX, gpE_upXY, gpE_upXYZ 
    real(rkind), dimension(:,:,:), allocatable :: f
    real(rkind), dimension(:,:,:), allocatable :: fxup_inX, fxup_inY, fxyup_inY 
    real(rkind), dimension(:,:,:), allocatable :: fxyup_inZ, fxyup_inX, fxyzup_inZ, fxyzup_inY, fxyzup_inX
    character(len=clen) :: tempname, fname
    real(rkind) :: tsim 
    namelist /INPUT/ nx, ny, nz, inputdir, outputdir, inputFile_TID, inputFile_RID, &
    outputFile_TID, outputFile_RID, UpsampleInZ, isStratified 

    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=INPUT)
    close(ioUnit)


    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(gpC)
    call decomp_info_init(nx,ny,nz+1,gpE)
    
    call decomp_info_init(2*nx,ny,nz,gpC_upX)
    call decomp_info_init(2*nx,2*ny,nz,gpC_upXY)
    call decomp_info_init(2*nx,2*ny,2*nz,gpC_upXYZ)
    call decomp_info_init(2*nx,ny,nz+1,gpE_upX)
    call decomp_info_init(2*nx,2*ny,nz+1,gpE_upXY)
    call decomp_info_init(2*nx,2*ny,2*nz+1,gpE_upXYZ)
  

    !!!!!! READ HEADER !!!!!!!!!!!
    write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",inputFile_RID, "_info.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    open(unit=10,file=fname,access='sequential',form='formatted')
    read (10, *)  tsim
    close(10)
    call message(0, "Upsampling File data dumped at tSIM=", tsim)

    !!!!!!!!!!!!! CELL FIELDS !!!!!!!!!!!!!!!
    allocate(f(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(fxup_inX(gpC_upX%xsz(1),gpC_upX%xsz(2),gpC_upX%xsz(3))) 
    allocate(fxup_inY(gpC_upX%ysz(1),gpC_upX%ysz(2),gpC_upX%ysz(3))) 
    allocate(fxyup_inY(gpC_upXY%ysz(1),gpC_upXY%ysz(2),gpC_upXY%ysz(3))) 
    
    if (UpsampleInZ) then
        allocate(fxyup_inZ(gpC_upXY%zsz(1),gpC_upXY%zsz(2),gpC_upXY%zsz(3)))
        allocate(fxyzup_inX(gpC_upXYZ%xsz(1),gpC_upXYZ%xsz(2),gpC_upXYZ%xsz(3))) 
        allocate(fxyzup_inY(gpC_upXYZ%ysz(1),gpC_upXYZ%ysz(2),gpC_upXYZ%ysz(3))) 
        allocate(fxyzup_inZ(gpC_upXYZ%zsz(1),gpC_upXYZ%zsz(2),gpC_upXYZ%zsz(3)))
    else
        allocate(fxyup_inX(gpC_upXY%xsz(1),gpC_upXY%xsz(2),gpC_upXY%xsz(3))) 
    end if 

    ! Read and Write u - field
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_u.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,f,fname, gpC)
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_u.",outputFile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)

    call upsampleX(f,fxup_inX)
    call transpose_x_to_y(fxup_inX,fxup_inY,gpC_upX)
    call upsampleY(fxup_inY,fxyup_inY)

    if (UpsampleInZ) then
        call transpose_y_to_z(fxyup_inY,fxyup_inZ,gpC_upXY)
        call upsampleZ_cells(fxyup_inZ,fxyzup_inZ)
        call transpose_z_to_y(fxyzup_inZ,fxyzup_inY,gpC_upXYZ)
        call transpose_y_to_x(fxyzup_inY,fxyzup_inX,gpC_upXYZ)
        call decomp_2d_write_one(1,fxyzup_inX,fname, gpC_upXYZ)
    else
        call transpose_y_to_x(fxyup_inY,fxyup_inX,gpC_upXY)
        call decomp_2d_write_one(1,fxyup_inX,fname, gpC_upXY)
    end if 
    
    ! Read and Write v - field
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_v.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,f,fname, gpC)
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_v.",outputFile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)

    call upsampleX(f,fxup_inX)
    call transpose_x_to_y(fxup_inX,fxup_inY,gpC_upX)
    call upsampleY(fxup_inY,fxyup_inY)

    if (UpsampleInZ) then
        call transpose_y_to_z(fxyup_inY,fxyup_inZ,gpC_upXY)
        call upsampleZ_cells(fxyup_inZ,fxyzup_inZ)
        call transpose_z_to_y(fxyzup_inZ,fxyzup_inY,gpC_upXYZ)
        call transpose_y_to_x(fxyzup_inY,fxyzup_inX,gpC_upXYZ)
        call decomp_2d_write_one(1,fxyzup_inX,fname, gpC_upXYZ)
    else
        call transpose_y_to_x(fxyup_inY,fxyup_inX,gpC_upXY)
        call decomp_2d_write_one(1,fxyup_inX,fname, gpC_upXY)
    end if 
    
    if (isStratified) then
        ! Read and Write T - field
        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_T.",inputFile_TID
        fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
        call decomp_2d_read_one(1,f,fname, gpC)
        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_T.",outputFile_TID
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)

        call upsampleX(f,fxup_inX)
        call transpose_x_to_y(fxup_inX,fxup_inY,gpC_upX)
        call upsampleY(fxup_inY,fxyup_inY)

        if (UpsampleInZ) then
            call transpose_y_to_z(fxyup_inY,fxyup_inZ,gpC_upXY)
            call upsampleZ_cells(fxyup_inZ,fxyzup_inZ)
            call transpose_z_to_y(fxyzup_inZ,fxyzup_inY,gpC_upXYZ)
            call transpose_y_to_x(fxyzup_inY,fxyzup_inX,gpC_upXYZ)
            call decomp_2d_write_one(1,fxyzup_inX,fname, gpC_upXYZ)
        else
            call transpose_y_to_x(fxyup_inY,fxyup_inX,gpC_upXY)
            call decomp_2d_write_one(1,fxyup_inX,fname, gpC_upXY)
        end if 
    end if 
    
    ! < SCALARS GO HERE >
    deallocate(f, fxup_inX, fxup_inY, fxyup_inY)
    if (UpsampleInZ) then
        deallocate(fxyup_inZ, fxyzup_inX, fxyzup_inY, fxyzup_inZ)
    else
        deallocate(fxyup_inX)    
    end if 

    !!!!!!!!!!!!! EDGE FIELDS !!!!!!!!!!!!!!!
    allocate(f(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
    allocate(fxup_inX(gpE_upX%xsz(1),gpE_upX%xsz(2),gpE_upX%xsz(3))) 
    allocate(fxup_inY(gpE_upX%ysz(1),gpE_upX%ysz(2),gpE_upX%ysz(3))) 
    allocate(fxyup_inY(gpE_upXY%ysz(1),gpE_upXY%ysz(2),gpE_upXY%ysz(3))) 
    
    if (UpsampleInZ) then
        allocate(fxyup_inZ(gpE_upXY%zsz(1),gpE_upXY%zsz(2),gpE_upXY%zsz(3)))
        allocate(fxyzup_inX(gpE_upXYZ%xsz(1),gpE_upXYZ%xsz(2),gpE_upXYZ%xsz(3))) 
        allocate(fxyzup_inY(gpE_upXYZ%ysz(1),gpE_upXYZ%ysz(2),gpE_upXYZ%ysz(3))) 
        allocate(fxyzup_inZ(gpE_upXYZ%zsz(1),gpE_upXYZ%zsz(2),gpE_upXYZ%zsz(3)))
    else
        allocate(fxyup_inX(gpE_upXY%xsz(1),gpE_upXY%xsz(2),gpE_upXY%xsz(3))) 
    end if 

    ! Read and Write w - field
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_w.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,f,fname, gpE)
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_w.",outputFile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)

    call upsampleX(f,fxup_inX)
    call transpose_x_to_y(fxup_inX,fxup_inY,gpE_upX)
    call upsampleY(fxup_inY,fxyup_inY)

    if (UpsampleInZ) then
        call transpose_y_to_z(fxyup_inY,fxyup_inZ,gpE_upXY)
        call upsampleZ_edges(fxyup_inZ,fxyzup_inZ)
        call transpose_z_to_y(fxyzup_inZ,fxyzup_inY,gpE_upXYZ)
        call transpose_y_to_x(fxyzup_inY,fxyzup_inX,gpE_upXYZ)
        call decomp_2d_write_one(1,fxyzup_inX,fname, gpE_upXYZ)
    else
        call transpose_y_to_x(fxyup_inY,fxyup_inX,gpE_upXY)
        call decomp_2d_write_one(1,fxyup_inX,fname, gpE_upXY)
    end if 
    
    deallocate(f, fxup_inX, fxup_inY, fxyup_inY)
    if (UpsampleInZ) then
        deallocate(fxyup_inZ, fxyzup_inX, fxyzup_inY, fxyzup_inZ)
    else
        deallocate(fxyup_inX)    
    end if 

    !!!!!! WRITE HEADER !!!!!!!!!!!
    if (nrank == 0) then
        write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",outputFile_RID, "_info.",outputFile_TID
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
        OPEN(UNIT=10, FILE=trim(fname))
        write(10,"(100g15.5)") tsim
        close(10)
    end if 

    call MPI_Finalize(ierr)           !<-- Terminate MPI 

end program 
