module RoutinesUpsampling
    use kind_parameters, only: rkind
    use exits, only: GracefulExit
    use constants, only : one
    
    implicit none
    contains
    subroutine upsampleX(arrIn, arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: i, j, k, nx, nxf, i1, i2
        real(rkind) :: interpfac, dx, dxf
        real(rkind), allocatable, dimension(:) :: x, xfine

        nx = size(arrIn, 1); nxf = size(arrOut, 1)

        if(nx==nxf) then
            arrOut = arrIn
        elseif(nxf==2*nx) then
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
        else
            allocate(x(nx), xfine(nxf))

            dx  = one/real(nx , rkind)
            do i = 1, nx
                x(i) = real(i-1, rkind)*dx
            enddo

            dxf = one/real(nxf, rkind)
            do i = 1, nxf
                xfine(i) = real(i-1, rkind)*dxf
            enddo

            arrOut(1,:,:) = arrIn(1,:,:)
            do i = 2, nxf
                i1 =  minloc(abs(xfine(i)-x),1)
                if(xfine(i) .lt. x(i1)) i1 = i1-1
                if(i1<1) then
                     call GracefulExit("i1 must be between 1 and nx. How did you get here!!??",999)
                endif
                i2 = i1 + 1; if(i2>nx) i2 = 1
                interpfac = (xfine(i)-x(i1))/(x(i2)-x(i1))
                arrOut(i,:,:) = arrIn(i1,:,:) + interpfac*(arrIn(i2,:,:) - arrIn(i1,:,:))
            enddo

            deallocate(x, xfine)
        endif


    end subroutine

    subroutine upsampleY(arrIn, arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: j, k, ny, nyf, j1, j2
        real(rkind) :: interpfac, dy, dyf
        real(rkind), allocatable, dimension(:) :: y, yfine

        ny = size(arrIn, 2); nyf = size(arrOut, 2)

        if(ny==nyf) then
            arrOut = arrIn
        elseif(nyf==2*ny) then
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
        else
            allocate(y(ny), yfine(nyf))

            dy  = one/real(ny , rkind)
            do j = 1, ny
                y(j) = real(j-1, rkind)*dy
            enddo

            dyf = one/real(nyf, rkind)
            do j = 1, nyf
                yfine(j) = real(j-1, rkind)*dyf
            enddo

            do k = 1, size(arrOut,3)
                arrOut(:,1,k) = arrIn(:,1,k)
            enddo
            do j = 2, nyf
                j1 =  minloc(abs(yfine(j)-y),1)
                if(yfine(j) .lt. y(j1)) j1 = j1-1
                if(j1<1) then
                     call GracefulExit("j1 must be between 1 and ny. How did you get here!!??",999)
                endif
                j2 = j1 + 1; if(j2>ny) j2 = 1
                interpfac = (yfine(j)-y(j1))/(y(j2)-y(j1))
                do k = 1, size(arrOut,3)
                    arrOut(:,j,k) = arrIn(:,j1,k) + interpfac*(arrIn(:,j2,k) - arrIn(:,j1,k))
                enddo
            enddo

            deallocate(y, yfine)
        endif

    end subroutine

    subroutine upsampleZ_cells(arrIn,arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: k, nz, idx, nzf, k1, k2
        real(rkind), parameter :: ratEven = 0.25d0, ratOdd = 0.75d0
        real(rkind), allocatable, dimension(:) :: z, zfine
        real(rkind) :: dz, dzf, interpfac

        nz = size(arrIn, 3); nzf = size(arrOut, 3)

        if(nzf==nz) then
            arrOut = arrIn
        elseif(nzf==2*nz) then
            do k = 1,nz-1
                arrOut(:,:,2*k+1) = arrIn(:,:,k) + (arrIn(:,:,k+1) - arrIn(:,:,k))*ratOdd 
            end do 
            
            idx = 1
            do k = 2,2*nz-1,2
                arrOut(:,:,k) = arrIn(:,:,idx) + (arrIn(:,:,idx + 1) - arrIn(:,:,idx))*ratEven
                idx = idx + 1
            end do
   
            arrOut(:,:,1) = 2.d0*arrOut(:,:,2) - arrOut(:,:,3)
            arrOut(:,:,2*nz) = 2.d0*arrOut(:,:,2*nz - 1) - arrOut(:,:,2*nz - 2)
        else
            allocate(z(nz), zfine(nzf))

            dz  = one/real(nz , rkind)
            do k = 1, nz
                z(k) = (real(k, rkind)-0.5_rkind)*dz
            enddo

            dzf = one/real(nzf, rkind)
            do k = 1, nzf
                zfine(k) = (real(k, rkind)-0.5_rkind)*dzf
            enddo

            do k = 1, nzf
                k1 =  minloc(abs(zfine(k)-z),1)
                if(zfine(k) .lt. z(k1)) k1 = k1-1
                k2 = k1+1
                if(k1<1) then
                    k1 = 1; k2 = 2
                endif
                if(k1==nz) then
                    k1 = nz-1; k2 = nz
                endif
                interpfac = (zfine(k)-z(k1))/(z(k2)-z(k1))
                arrOut(:,:,k) = arrIn(:,:,k1) + interpfac*(arrIn(:,:,k2) - arrIn(:,:,k1))
            enddo

            deallocate(z, zfine)
        endif

    end subroutine


    subroutine upsampleZ_edges(arrIn,arrOut)
        real(rkind), dimension(:,:,:), intent(in) :: arrIn
        real(rkind), dimension(:,:,:), intent(out) :: arrOut
        integer :: k, nz, nzE, nzEf, k1, k2
        real(rkind), allocatable, dimension(:) :: zE, zEfine
        real(rkind) :: dzE, dzEf, interpfac

        nzE = size(arrIn, 3); nzEf = size(arrOut, 3)

        if(nzEf==nzE) then
            arrOut = arrIn
        elseif(nzEf==2*nzE+1) then
            nz = size(arrIn,3) - 1

            do k = 1,nz
                arrOut(:,:,2*k-1) = arrIn(:,:,k)
            end do 
            arrOut(:,:,2*nz+1) = arrOut(:,:,nz+1)

            do k = 1,nz
                arrOut(:,:,2*k) = 0.5d0*(arrOut(:,:,2*k-1) + arrOut(:,:,2*k+1))
            end do 

            !arrOut(:,:,2*nz+1) = 0.d0
            !arrOut(:,:,1) = 0.d0
        else
            allocate(zE(nzE), zEfine(nzEf))

            dzE  = one/real(nzE-1 , rkind)
            do k = 1, nzE
                zE(k) = real(k-1, rkind)*dzE
            enddo

            dzEf = one/real(nzEf-1, rkind)
            do k = 1, nzEf
                zEfine(k) = real(k-1, rkind)*dzEf
            enddo

            arrOut(:,:,1) = arrIn(:,:,1)
            do k = 2, nzEf-1
                k1 =  minloc(abs(zEfine(k)-zE),1)
                if(zEfine(k) .lt. zE(k1)) k1 = k1-1
                if(k1<1) then
                    call GracefulExit("k1 must be between 1 and nzE. How did you get here!!??",999)
                endif
                if(k1==nzE) then
                    k1 = nzE-1
                endif
                k2 = k1+1
                interpfac = (zEfine(k)-zE(k1))/(zE(k2)-zE(k1))
                arrOut(:,:,k) = arrIn(:,:,k1) + interpfac*(arrIn(:,:,k2) - arrIn(:,:,k1))
            enddo
            arrOut(:,:,nzEf) = arrIn(:,:,nzE)

            deallocate(zE, zEfine)
        endif

    end subroutine

    subroutine upsampleZ_periodic(fC, fE, fupC, fupE)
        real(rkind), dimension(:,:,:), intent(in) :: fC, fE
        real(rkind), dimension(:,:,:), intent(out) :: fupC, fupE
        integer :: k, nz, nzE, nzEf, k1, k2 
        real(rkind), allocatable, dimension(:) :: zE, zEfine
        real(rkind) :: dzE, dzEf, interpfac

        nzE = size(fE, 3); nzEf = size(fupE, 3)

        if(nzEf==nzE) then
            fupC = fC; fupE = fE
        elseif(nzEf==2*nzE+1) then
            nz = size(fC, 3)
            do k = 1,nz+1
                fupE(:,:,2*k-1) = fE(:,:,k)
            end do 

            do k = 1,nz
                fupE(:,:,2*k) = fC(:,:,k)
            end do 

            fupC = 0.5d0*(fupE(:,:,2:2*nz+1) + fupE(:,:,1:2*nz))
        else
            allocate(zE(nzE), zEfine(nzEf))

            dzE  = one/real(nzE-1 , rkind)
            do k = 1, nzE
                zE(k) = real(k-1, rkind)*dzE
            enddo

            dzEf = one/real(nzEf-1, rkind)
            do k = 1, nzEf
                zEfine(k) = real(k-1, rkind)*dzEf
            enddo

            fupE(:,:,1) = fE(:,:,1)
            do k = 2, nzEf-1
                k1 =  minloc(abs(zEfine(k)-zE),1)
                if(zEfine(k) .lt. zE(k1)) k1 = k1-1
                if(k1<1 .or. k1>nzE-1) then
                    call GracefulExit("k1 must be between 1 and nzE-1. How did you get here!!??",999)
                endif
                k2 = k1+1
                interpfac = (zEfine(k)-zE(k1))/(zE(k2)-zE(k1))
                fupE(:,:,k) = fE(:,:,k1) + interpfac*(fE(:,:,k2) - fE(:,:,k1))
            enddo
            fupE(:,:,nzEf) = fE(:,:,nzE)
        endif

    end subroutine 
end module



program upsampleFields
    use kind_parameters, only: rkind, clen
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
    integer :: ioUnit, nx, ny, nz, ierr, inputFile_TID, inputFile_RID, outputFile_TID, outputFile_RID, i1, j1, k1, i, j, k
    logical :: upsampleInZ = .false., isStratified = .false. 
    type(decomp_info) :: gpC, gpE, gpC_upX, gpC_upXY, gpC_upXYZ, gpE_upX, gpE_upXY, gpE_upXYZ 
    real(rkind), dimension(:,:,:), allocatable :: f, fxyupE_inZ, fxyzupE_inZ
    real(rkind), dimension(:,:,:), allocatable :: fxup_inX, fxup_inY, fxyup_inY 
    real(rkind), dimension(:,:,:), allocatable :: fxyup_inZ, fxyup_inX, fxyzup_inZ, fxyzup_inY, fxyzup_inX
    character(len=clen) :: tempname, fname
    real(rkind) :: tsim 
    logical :: periodicInZ = .false. 
    integer :: nxf, nyf, nzf
    namelist /INPUT/ nx, ny, nz, inputdir, outputdir, inputFile_TID, inputFile_RID, &
    outputFile_TID, outputFile_RID, UpsampleInZ, isStratified, PeriodicInZ, &
    nxf, nyf, nzf

    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=INPUT)
    close(ioUnit)

    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(gpC)
    call decomp_info_init(nx,ny,nz+1,gpE)
    
    call decomp_info_init(nxf, ny , nz   , gpC_upX  )
    call decomp_info_init(nxf, nyf, nz   , gpC_upXY )
    call decomp_info_init(nxf, nyf, nzf  , gpC_upXYZ)
    call decomp_info_init(nxf, ny , nz+1 , gpE_upX  )
    call decomp_info_init(nxf, nyf, nz+1 , gpE_upXY )
    call decomp_info_init(nxf, nyf, nzf+1, gpE_upXYZ)
  

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
        allocate(fxyupE_inZ(gpC_upXY%zsz(1),gpC_upXY%zsz(2),gpC_upXY%zsz(3)+1))
        allocate(fxyzup_inX(gpC_upXYZ%xsz(1),gpC_upXYZ%xsz(2),gpC_upXYZ%xsz(3))) 
        allocate(fxyzup_inY(gpC_upXYZ%ysz(1),gpC_upXYZ%ysz(2),gpC_upXYZ%ysz(3))) 
        allocate(fxyzup_inZ(gpC_upXYZ%zsz(1),gpC_upXYZ%zsz(2),gpC_upXYZ%zsz(3)))
        allocate(fxyzupE_inZ(gpE_upXYZ%zsz(1),gpE_upXYZ%zsz(2),gpE_upXYZ%zsz(3)))
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
!!!!!!!CHECKING!!!!!!!!
    j1 = 3; k1 = 4;
    open(101,file='check_upsample_x.dat',status='unknown',action='write')
    write(101,*) 'VARIABLES="x","val"'
    write(101,*) 'ZONE T="Sample X", F=POINT, I=', size(f,1)
    do i = 1, size(f,1)
    write(101,*) real(i-1,rkind)/real(size(f,1),rkind), f(i,j1,k1)
    enddo
    write(101,*) 'ZONE T="Upsample X", F=POINT, I=', size(fxup_inX,1)
    do i = 1, size(fxup_inX,1)
    write(101,*) real(i-1,rkind)/real(size(fxup_inX,1),rkind), fxup_inX(i,j1,k1)
    enddo
    close(101)
!!!!!!!CHECKING!!!!!!!!
    call transpose_x_to_y(fxup_inX,fxup_inY,gpC_upX)
    call upsampleY(fxup_inY,fxyup_inY)
!!!!!!!CHECKING!!!!!!!!
    i1 = 3; k1 = 4;
    open(101,file='check_upsample_y.dat',status='unknown',action='write')
    write(101,*) 'VARIABLES="y","val"'
    write(101,*) 'ZONE T="Sample Y", F=POINT, I=', size(fxup_inY,2)
    do j = 1, size(fxup_inY,2)
    write(101,*) real(j-1,rkind)/real(size(fxup_inY,2),rkind), fxup_inY(i1,j,k1)
    enddo
    write(101,*) 'ZONE T="Upsample Y", F=POINT, I=', size(fxyup_inY,2)
    do j = 1, size(fxyup_inY,2)
    write(101,*) real(j-1,rkind)/real(size(fxyup_inY,2),rkind), fxyup_inY(i1,j,k1)
    enddo
    close(101)
!!!!!!!CHECKING!!!!!!!!

    if (UpsampleInZ) then
        call transpose_y_to_z(fxyup_inY,fxyup_inZ,gpC_upXY)
        fxyupE_inZ(:,:,2:nz) = 0.5d0*(fxyup_inZ(:,:,1:nz-1) + fxyup_inZ(:,:,2:nz))
        fxyupE_inZ(:,:,1)    = 0.5d0*(fxyup_inZ(:,:,1) + fxyup_inZ(:,:,nz))
        fxyupE_inZ(:,:,nz+1) = fxyupE_inZ(:,:,1)

        if (PeriodicInZ) then
           call upsampleZ_periodic(fxyup_inZ, fxyupE_inZ, fxyzup_inZ, fxyzupE_inZ)
        else
           call upsampleZ_cells(fxyup_inZ,fxyzup_inZ)
        end if 
!!!!!!!CHECKING!!!!!!!!
    i1 = 3; j1 = 4;
    open(101,file='check_upsample_z.dat',status='unknown',action='write')
    write(101,*) 'VARIABLES="z","val"'
    write(101,*) 'ZONE T="Sample Z", F=POINT, I=', size(fxyup_inZ,3)
    do k = 1, size(fxyup_inZ,3)
    write(101,*) (real(k,rkind)-0.5d0)/real(size(fxyup_inZ,3),rkind), fxyup_inZ(i1,j1,k)
    enddo
    write(101,*) 'ZONE T="Upsample Z", F=POINT, I=', size(fxyzup_inZ,3)
    do k = 1, size(fxyzup_inZ,3)
    write(101,*) (real(k,rkind)-0.5d0)/real(size(fxyzup_inZ,3),rkind), fxyzup_inZ(i1,j1,k)
    enddo
    write(101,*) 'ZONE T="Sample Z E", F=POINT, I=', size(fxyupE_inZ,3)
    do k = 1, size(fxyupE_inZ,3)
    write(101,*) (real(k,rkind)-0.5d0)/real(size(fxyupE_inZ,3),rkind), fxyupE_inZ(i1,j1,k)
    enddo
    write(101,*) 'ZONE T="Upsample Z E", F=POINT, I=', size(fxyzupE_inZ,3)
    do k = 1, size(fxyzupE_inZ,3)
    write(101,*) (real(k,rkind)-0.5d0)/real(size(fxyzupE_inZ,3),rkind), fxyzupE_inZ(i1,j1,k)
    enddo
    close(101)
!!!!!!!CHECKING!!!!!!!!

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
        fxyupE_inZ(:,:,2:nz) = 0.5d0*(fxyup_inZ(:,:,1:nz-1) + fxyup_inZ(:,:,2:nz))
        fxyupE_inZ(:,:,1)    = 0.5d0*(fxyup_inZ(:,:,1) + fxyup_inZ(:,:,nz))
        fxyupE_inZ(:,:,nz+1) = fxyupE_inZ(:,:,1)

        if (PeriodicInZ) then
           call upsampleZ_periodic(fxyup_inZ, fxyupE_inZ, fxyzup_inZ, fxyzupE_inZ)
        else
           call upsampleZ_cells(fxyup_inZ,fxyzup_inZ)
        end if 
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
            fxyupE_inZ(:,:,2:nz) = 0.5d0*(fxyup_inZ(:,:,1:nz-1) + fxyup_inZ(:,:,2:nz))
            fxyupE_inZ(:,:,1)    = 0.5d0*(fxyup_inZ(:,:,1) + fxyup_inZ(:,:,nz))
            fxyupE_inZ(:,:,nz+1) = fxyupE_inZ(:,:,1)

            if (PeriodicInZ) then
               call upsampleZ_periodic(fxyup_inZ, fxyupE_inZ, fxyzup_inZ, fxyzupE_inZ)
            else
               call upsampleZ_cells(fxyup_inZ,fxyzup_inZ)
            end if 
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
