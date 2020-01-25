module RoutinesDuplication
    use kind_parameters, only: rkind
    use exits, only: GracefulExit
    use constants, only : one
    
    implicit none
    contains
    subroutine duplicate_X(arrIn, arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: i, j, k, nx, nxf, i1, i2
        real(rkind) :: interpfac, dx, dxf
        real(rkind), allocatable, dimension(:) :: x, xfine

        nx = size(arrIn, 1); nxf = size(arrOut, 1)

        if(nx > nxf) then
            call GracefulExit("Check details. nx > nxf", 111)
        else
            do k = 1,size(arrIn,3)
                do j = 1,size(arrIn,2)
                    do i = 1,nx
                        arrOut(i,j,k) = arrIn(i,j,k)
                    end do 
                    do i = nx+1,size(arrOut,1)
                        arrOut(i,j,k) = arrIn(i-nx,j,k)
                    end do 
                end do 
            end do 
        endif
    end subroutine

    subroutine duplicate_Y(arrIn, arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: j, k, ny, nyf, j1, j2
        real(rkind) :: interpfac, dy, dyf
        real(rkind), allocatable, dimension(:) :: y, yfine

        ny = size(arrIn, 2); nyf = size(arrOut, 2)

        if(ny > nyf) then
            call GracefulExit("Check details. ny > nyf", 111)
        else
            do k = 1,size(arrIn,3)
                do j = 1,ny
                    arrOut(:,j,k) = arrIn(:,j,k)
                end do 
                do j = ny+1,size(arrOut,2)
                    arrOut(:,j,k) = arrIn(:,j-ny,k)
                end do 
            end do 
        endif
    end subroutine

    subroutine duplicate_Z(arrIn,arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: k, nz, idx, nzf, k1, k2
        real(rkind), parameter :: ratEven = 0.25d0, ratOdd = 0.75d0
        real(rkind), allocatable, dimension(:) :: z, zfine
        real(rkind) :: dz, dzf, interpfac

        nz = size(arrIn, 3); nzf = size(arrOut, 3)

        if(nz > nzf) then
            call GracefulExit("Check details. nz > nzf", 111)
        else
            do k = 1,nz
                arrOut(:,:,k) = arrIn(:,:,k)
            end do 
            do k = nz+1,size(arrOut,3)
                arrOut(:,:,k) = arrIn(:,:,k-nz)
            end do 
        endif
    end subroutine


    !subroutine duplicateZ_edges(arrIn,arrOut)
    !    real(rkind), dimension(:,:,:), intent(in) :: arrIn
    !    real(rkind), dimension(:,:,:), intent(out) :: arrOut
    !    integer :: k, nz, nzE, nzEf, k1, k2
    !    real(rkind), allocatable, dimension(:) :: zE, zEfine
    !    real(rkind) :: dzE, dzEf, interpfac

    !    nzE = size(arrIn, 3); nzEf = size(arrOut, 3)

    !    if(nzEf==nzE) then
    !        arrOut = arrIn
    !    elseif(nzEf==2*nzE+1) then
    !        nz = size(arrIn,3) - 1

    !        do k = 1,nz
    !            arrOut(:,:,2*k-1) = arrIn(:,:,k)
    !        end do 
    !        arrOut(:,:,2*nz+1) = arrOut(:,:,nz+1)

    !        do k = 1,nz
    !            arrOut(:,:,2*k) = 0.5d0*(arrOut(:,:,2*k-1) + arrOut(:,:,2*k+1))
    !        end do 

    !        !arrOut(:,:,2*nz+1) = 0.d0
    !        !arrOut(:,:,1) = 0.d0
    !    else
    !        allocate(zE(nzE), zEfine(nzEf))

    !        dzE  = one/real(nzE-1 , rkind)
    !        do k = 1, nzE
    !            zE(k) = real(k-1, rkind)*dzE
    !        enddo

    !        dzEf = one/real(nzEf-1, rkind)
    !        do k = 1, nzEf
    !            zEfine(k) = real(k-1, rkind)*dzEf
    !        enddo

    !        arrOut(:,:,1) = arrIn(:,:,1)
    !        do k = 2, nzEf-1
    !            k1 =  minloc(abs(zEfine(k)-zE),1)
    !            if(zEfine(k) .lt. zE(k1)) k1 = k1-1
    !            if(k1<1) then
    !                call GracefulExit("k1 must be between 1 and nzE. How did you get here!!??",999)
    !            endif
    !            if(k1==nzE) then
    !                k1 = nzE-1
    !            endif
    !            k2 = k1+1
    !            interpfac = (zEfine(k)-zE(k1))/(zE(k2)-zE(k1))
    !            arrOut(:,:,k) = arrIn(:,:,k1) + interpfac*(arrIn(:,:,k2) - arrIn(:,:,k1))
    !        enddo
    !        arrOut(:,:,nzEf) = arrIn(:,:,nzE)

    !        deallocate(zE, zEfine)
    !    endif

    !end subroutine

    !subroutine duplicateZ_periodic(fC, fE, fupC, fupE)
    !    real(rkind), dimension(:,:,:), intent(in) :: fC, fE
    !    real(rkind), dimension(:,:,:), intent(out) :: fupC, fupE
    !    integer :: k, nz, nzE, nzEf, k1, k2 
    !    real(rkind), allocatable, dimension(:) :: zE, zEfine
    !    real(rkind) :: dzE, dzEf, interpfac

    !    nzE = size(fE, 3); nzEf = size(fupE, 3)

    !    if(nzEf==nzE) then
    !        fupC = fC; fupE = fE
    !    elseif(nzEf==2*nzE+1) then
    !        nz = size(fC, 3)
    !        do k = 1,nz+1
    !            fupE(:,:,2*k-1) = fE(:,:,k)
    !        end do 

    !        do k = 1,nz
    !            fupE(:,:,2*k) = fC(:,:,k)
    !        end do 

    !        fupC = 0.5d0*(fupE(:,:,2:2*nz+1) + fupE(:,:,1:2*nz))
    !    else
    !        allocate(zE(nzE), zEfine(nzEf))

    !        dzE  = one/real(nzE-1 , rkind)
    !        do k = 1, nzE
    !            zE(k) = real(k-1, rkind)*dzE
    !        enddo

    !        dzEf = one/real(nzEf-1, rkind)
    !        do k = 1, nzEf
    !            zEfine(k) = real(k-1, rkind)*dzEf
    !        enddo

    !        fupE(:,:,1) = fE(:,:,1)
    !        do k = 2, nzEf-1
    !            k1 =  minloc(abs(zEfine(k)-zE),1)
    !            if(zEfine(k) .lt. zE(k1)) k1 = k1-1
    !            if(k1<1 .or. k1>nzE-1) then
    !                call GracefulExit("k1 must be between 1 and nzE-1. How did you get here!!??",999)
    !            endif
    !            k2 = k1+1
    !            interpfac = (zEfine(k)-zE(k1))/(zE(k2)-zE(k1))
    !            fupE(:,:,k) = fE(:,:,k1) + interpfac*(fE(:,:,k2) - fE(:,:,k1))
    !        enddo
    !        fupE(:,:,nzEf) = fE(:,:,nzE)
    !    endif

    !end subroutine 
end module



program duplicateFields
    use kind_parameters, only: rkind, clen
    use gridtools, only: alloc_buffs, destroy_buffs
    use decomp_2d
    use exits, only: GracefulExit, message
    use mpi 
    use timer, only: tic, toc
    use decomp_2d_io
    use RoutinesDuplication
   
    implicit none

    character(len=clen) :: inputfile 
    character(len=clen) :: outputdir, inputdir
    integer :: ioUnit, nx, ny, nz, ierr, inputFile_TID, inputFile_RID, outputFile_TID, outputFile_RID
    logical :: duplicateInX, duplicateInY, duplicateInZ, isStratified = .false. 
    type(decomp_info) :: gpC, gpE, gpC_dupX, gpC_dupXY, gpC_dupXYZ, gpE_dupX, gpE_dupXY, gpE_dupXYZ 
    real(rkind), dimension(:,:,:), allocatable :: f, fxydupE_inZ, fxyzdupE_inZ
    real(rkind), dimension(:,:,:), allocatable :: fxdup_inX, fxdup_inY, fxydup_inY 
    real(rkind), dimension(:,:,:), allocatable :: fxydup_inZ, fxydup_inX, fxyzdup_inZ, fxyzdup_inY, fxyzdup_inX
    character(len=clen) :: tempname, fname
    real(rkind) :: tsim 
    logical :: periodicInZ = .false. 
    integer :: nxf, nyf, nzf
    namelist /INPUT/ nx, ny, nz, inputdir, outputdir, inputFile_TID, inputFile_RID, &
    outputFile_TID, outputFile_RID, duplicateInX, duplicateInY, duplicateInZ, isStratified, nxf, nyf, nzf

    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=INPUT)
    close(ioUnit)

    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(gpC)
    call decomp_info_init(nx,ny,nz+1,gpE)
    if (duplicateInX) then   
        if(nxf < nx) then
          call GracefulExit("nxf cannot be less than nx", 111)
        endif
    endif
    if (duplicateInY) then   
        if(nyf < ny) then
          call GracefulExit("nyf cannot be less than ny", 111)
        endif
    endif
    if (duplicateInZ) then   
        if(nzf < nz) then
          call GracefulExit("nzf cannot be less than nz", 111)
        endif
    endif
    call decomp_info_init(nxf, ny , nz   , gpC_dupX  )
    call decomp_info_init(nxf, nyf, nz   , gpC_dupXY )
    call decomp_info_init(nxf, nyf, nzf  , gpC_dupXYZ)
    call decomp_info_init(nxf, ny , nz+1 , gpE_dupX  )
    call decomp_info_init(nxf, nyf, nz+1 , gpE_dupXY )
    call decomp_info_init(nxf, nyf, nzf+1, gpE_dupXYZ)
  

    !!!!!! READ HEADER !!!!!!!!!!!
    write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",inputFile_RID, "_info.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    open(unit=10,file=fname,access='sequential',form='formatted')
    read (10, *)  tsim
    close(10)
    call message(0, "Upsampling File data dumped at tSIM=", tsim)

    !!!!!!!!!!!!! CELL FIELDS !!!!!!!!!!!!!!!
    allocate(f(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(fxdup_inX(gpC_dupX%xsz(1),gpC_dupX%xsz(2),gpC_dupX%xsz(3)))
    allocate(fxdup_inY(gpC_dupX%ysz(1),gpC_dupX%ysz(2),gpC_dupX%ysz(3)))
    allocate(fxydup_inY(gpC_dupXY%ysz(1),gpC_dupXY%ysz(2),gpC_dupXY%ysz(3)))
    allocate(fxydup_inZ(gpC_dupXY%zsz(1),gpC_dupXY%zsz(2),gpC_dupXY%zsz(3)))
    allocate(fxyzdup_inX(gpC_dupXYZ%xsz(1),gpC_dupXYZ%xsz(2),gpC_dupXYZ%xsz(3))) 
    allocate(fxyzdup_inY(gpC_dupXYZ%ysz(1),gpC_dupXYZ%ysz(2),gpC_dupXYZ%ysz(3))) 
    allocate(fxyzdup_inZ(gpC_dupXYZ%zsz(1),gpC_dupXYZ%zsz(2),gpC_dupXYZ%zsz(3)))

    ! Read and Write u - field
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_u.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,f,fname, gpC)
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_u.",outputFile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)

    call duplicate_X(f,fxdup_inX)
    call transpose_x_to_y(fxdup_inX,fxdup_inY,gpC_dupX)
    call duplicate_Y(fxdup_inY,fxydup_inY)
    call transpose_y_to_z(fxydup_inY,fxydup_inZ,gpC_dupXY)
    call duplicate_Z(fxydup_inZ,fxyzdup_inZ)
    call transpose_z_to_y(fxyzdup_inZ,fxyzdup_inY,gpC_dupXYZ)
    call transpose_y_to_x(fxyzdup_inY,fxyzdup_inX,gpC_dupXYZ)
    call decomp_2d_write_one(1,fxyzdup_inX,fname, gpC_dupXYZ)
    
    ! Read and Write v - field
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_v.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,f,fname, gpC)
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_v.",outputFile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)

    call duplicate_X(f,fxdup_inX)
    call transpose_x_to_y(fxdup_inX,fxdup_inY,gpC_dupX)
    call duplicate_Y(fxdup_inY,fxydup_inY)
    call transpose_y_to_z(fxydup_inY,fxydup_inZ,gpC_dupXY)
    call duplicate_Z(fxydup_inZ,fxyzdup_inZ)
    call transpose_z_to_y(fxyzdup_inZ,fxyzdup_inY,gpC_dupXYZ)
    call transpose_y_to_x(fxyzdup_inY,fxyzdup_inX,gpC_dupXYZ)
    call decomp_2d_write_one(1,fxyzdup_inX,fname, gpC_dupXYZ)
    
    if (isStratified) then
        ! Read and Write T - field
        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_T.",inputFile_TID
        fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
        call decomp_2d_read_one(1,f,fname, gpC)
        write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_T.",outputFile_TID
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)

        call duplicate_X(f,fxdup_inX)
        call transpose_x_to_y(fxdup_inX,fxdup_inY,gpC_dupX)
        call duplicate_Y(fxdup_inY,fxydup_inY)
        call transpose_y_to_z(fxydup_inY,fxydup_inZ,gpC_dupXY)
        call duplicate_Z(fxydup_inZ,fxyzdup_inZ)
        call transpose_z_to_y(fxyzdup_inZ,fxyzdup_inY,gpC_dupXYZ)
        call transpose_y_to_x(fxyzdup_inY,fxyzdup_inX,gpC_dupXYZ)
        call decomp_2d_write_one(1,fxyzdup_inX,fname, gpC_dupXYZ)
    end if 
    
    ! < SCALARS GO HERE >
    deallocate(f, fxdup_inX, fxdup_inY, fxydup_inY)
    deallocate(fxydup_inZ, fxyzdup_inX, fxyzdup_inY, fxyzdup_inZ)

    !!!!!!!!!!!!! EDGE FIELDS !!!!!!!!!!!!!!!
    allocate(f(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
    allocate(fxdup_inX(gpE_dupX%xsz(1),gpE_dupX%xsz(2),gpE_dupX%xsz(3))) 
    allocate(fxdup_inY(gpE_dupX%ysz(1),gpE_dupX%ysz(2),gpE_dupX%ysz(3))) 
    allocate(fxydup_inY(gpE_dupXY%ysz(1),gpE_dupXY%ysz(2),gpE_dupXY%ysz(3))) 
    allocate(fxydup_inZ(gpE_dupXY%zsz(1),gpE_dupXY%zsz(2),gpE_dupXY%zsz(3)))
    allocate(fxyzdup_inX(gpE_dupXYZ%xsz(1),gpE_dupXYZ%xsz(2),gpE_dupXYZ%xsz(3))) 
    allocate(fxyzdup_inY(gpE_dupXYZ%ysz(1),gpE_dupXYZ%ysz(2),gpE_dupXYZ%ysz(3))) 
    allocate(fxyzdup_inZ(gpE_dupXYZ%zsz(1),gpE_dupXYZ%zsz(2),gpE_dupXYZ%zsz(3)))

    ! Read and Write w - field
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_w.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,f,fname, gpE)
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_w.",outputFile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)

    call duplicate_X(f,fxdup_inX)
    call transpose_x_to_y(fxdup_inX,fxdup_inY,gpE_dupX)
    call duplicate_Y(fxdup_inY,fxydup_inY)
    call transpose_y_to_z(fxydup_inY,fxydup_inZ,gpE_dupXY)
    call duplicate_Z(fxydup_inZ,fxyzdup_inZ)
    call transpose_z_to_y(fxyzdup_inZ,fxyzdup_inY,gpE_dupXYZ)
    call transpose_y_to_x(fxyzdup_inY,fxyzdup_inX,gpE_dupXYZ)
    call decomp_2d_write_one(1,fxyzdup_inX,fname, gpE_dupXYZ)
    
    deallocate(f, fxdup_inX, fxdup_inY, fxydup_inY)
    deallocate(fxydup_inZ, fxyzdup_inX, fxyzdup_inY, fxyzdup_inZ)

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
