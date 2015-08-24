module hitCD_IO

    use kind_parameters, only: rkind, clen
    use decomp_2d,       only: decomp_info,nrank,nproc
    use exits,           only: GracefulExit, message
    implicit none

    integer, dimension(:,:), allocatable        :: xst,xen,xsz
    integer :: headerfid = 101
    integer :: NumDumps 

    integer :: RunIDX = 2
contains

    subroutine write_matlab_header(OutputDir)
        use mpi 
        character(len=clen) :: fname
        character(len=clen) :: tempname
        character(len=clen) :: command
        character(len=*), intent(in) :: OutputDir
        integer :: ierr, system 
        logical :: isThere
        integer :: idx 

        inquire(FILE=trim(OutputDir), exist=isThere)
        if (nrank == 0) then
            if (.not. isThere) then
                print*, "============================================="
                print*, "WARNING: Output directory not found. Creating a new one."
                print*, "============================================="
                command = "mkdir "//trim(OutputDir)
                ierr = system(trim(command))
            end if 
            write(tempname,"(A3,I2.2,A6,A4)") "Run", RunIDX, "HEADER",".txt"
            fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)

            open (headerfid, file=trim(fname), FORM='formatted', STATUS='replace',ACTION='write')
            write(headerfid,*)"========================================================================="
            write(headerfid,*)"---------------------  Header file for MATLAB ---------------------------"
            write(headerfid,"(A9,A10,A10,A10,A10,A10,A10)") "PROC", "xst", "xen", "yst", "yen","zst","zen"
            write(headerfid,*)"-------------------------------------------------------------------------"
            do idx = 0,nproc-1
                write(headerfid,"(I8,6I10)") idx, xst(idx,1), xen(idx,1), xst(idx,2), xen(idx,2), xst(idx,3), xen(idx,3)
            end do 
            write(headerfid,*)"-------------------------------------------------------------------------"
            write(headerfid,*)"Dumps made at:"
        end if
        numDumps = 0 
        call mpi_barrier(mpi_comm_world,ierr)
    end subroutine

    subroutine dumpData4Matlab(tid, outputdir,fieldsinY,decomp) 
        use GridMod,            only: alloc_buffs
        use decomp_2d,        only: transpose_y_to_x
        integer, intent(in) :: tid
        type(decomp_info), intent(in) :: decomp
        character(len=clen) :: fname
        character(len=clen) :: tempname
        real(rkind), dimension(:,:,:,:), intent(in) :: fieldsinY
        character(len=*), intent(in) :: OutputDir
        real(rkind), dimension(:,:,:,:), allocatable :: fieldsPhys
        integer :: fid = 1234


        call alloc_buffs(fieldsPhys,3,'x',decomp)
        ! transpose from y to x
        call transpose_y_to_x(fieldsinY(:,:,:,1),fieldsPhys(:,:,:,1),decomp)
        call transpose_y_to_x(fieldsinY(:,:,:,2),fieldsPhys(:,:,:,2),decomp)
        call transpose_y_to_x(fieldsinY(:,:,:,3),fieldsPhys(:,:,:,3),decomp)

        write(tempname,"(A3,I2.2,A2,I4.4,A2,I5.5,A5,A4)") "Run", RunIDX, "_p",nrank,"_t",tid,"_uVEL",".out"
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
        open(fid,file=trim(fname),form='unformatted',status='replace')
        write(fid) fieldsPhys(:,:,:,1)
        close(fid)

        write(tempname,"(A3,I2.2,A2,I4.4,A2,I5.5,A5,A4)") "Run", RunIDX, "_p",nrank,"_t",tid,"_vVEL",".out"
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
        open(fid,file=trim(fname),form='unformatted',status='replace')
        write(fid) fieldsPhys(:,:,:,2)
        close(fid)

        write(tempname,"(A3,I2.2,A2,I4.4,A2,I5.5,A5,A4)") "Run", RunIDX, "_p",nrank,"_t",tid,"_wVEL",".out"
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
        open(fid,file=trim(fname),form='unformatted',status='replace')
        write(fid) fieldsPhys(:,:,:,3)
        close(fid)
        
        if (nrank == 0) then
            write(headerfid,"(I8)") tid
        end if  
        NumDumps = NumDumps + 1

        deallocate(fieldsPhys)
    end subroutine 

    subroutine getHit3d_uvw(Nx,Ny,Nz,fieldsPhys,gp,dir)
        use kind_parameters, only: mpirkind
        use mpi
        class( decomp_info ), intent(in)                  :: gp
        integer, intent(in)                               :: Nx,Ny,Nz
        real(rkind), dimension(:,:,:,:), intent(inout)    :: fieldsPhys
        character(len=*), intent(in)                      :: dir
        character(len=8),parameter                        :: uFile = "U.000000"
        character(len=8),parameter                        :: vFile = "V.000000"
        character(len=8),parameter                        :: wFile = "W.000000"
        character(:), allocatable                         :: fname
        real(rkind), dimension(:,:,:), allocatable        :: full_Field
        integer :: fid = 1234, i, j, k
        integer :: tag, idx, status(MPI_STATUS_SIZE), ierr
        integer :: sizes(3), chunksize

        call message("Reading in Hit3D U, V, W data from directory: "//trim(dir))

        ! Create data sharing info
        if (nrank == 0) then
            allocate(xst(0:nproc-1,3),xen(0:nproc-1,3),xsz(0:nproc-1,3))
            allocate(full_Field(Nx,Ny,Nz))
        end if


        ! communicate local processor grid info (Assume x-decomposition)
        if (nrank == 0) then
            xst(0,:) = gp%xst
            xen(0,:) = gp%xen
            
            tag = 0
            do idx = 1,nproc-1
                call MPI_RECV(xst(idx,:), 3, MPI_INTEGER, idx, tag,&
                    MPI_COMM_WORLD, status, ierr)
            end do 
            tag = 1
            do idx = 1,nproc-1
                call MPI_RECV(xen(idx,:), 3, MPI_INTEGER, idx, tag,&
                    MPI_COMM_WORLD, status, ierr)
            end do
           tag = 2
            do idx = 1,nproc-1
                call MPI_RECV(xsz(idx,:), 3, MPI_INTEGER, idx, tag,&
                    MPI_COMM_WORLD, status, ierr)
            end do

        else
            tag = 0
            call MPI_SEND(gp%xst, 3, MPI_INTEGER, 0, tag, &
                 &      MPI_COMM_WORLD, ierr)
            tag = 1
            call MPI_SEND(gp%xen, 3, MPI_INTEGER, 0, tag, &
                 &      MPI_COMM_WORLD, ierr)
            tag = 2
            call MPI_SEND(gp%xsz, 3, MPI_INTEGER, 0, tag, &
                 &      MPI_COMM_WORLD, ierr)

        end if 
        
        ! Read the u data  
        if (nrank == 0) then
            ! U data 
            fname = dir(:len_trim(dir))//'/'//uFile
            open(fid,file=fname,form='unformatted', access='stream')  
            read(fid) sizes(1:3)
            if ((sizes(1) == Nx) .and. (sizes(2) == Ny) .and. (sizes(3) == Nz)) then
                read(fid) (((full_field(i,j,k),i=1,Nx),j=1,Ny),k=1,Nz)
            else
                close(fid)
                print*, "Problem size      :", Nx, Ny, Nz
                print*, "Hit 3d input size :", sizes 
                call GracefulExit("Problem size mismatch with HIT_3d input file:"//uFile,1001) 
            end if
            tag = 10
            fieldsPhys(:,:,:,1) = full_field(xst(0,1):xen(0,1), xst(0,2):xen(0,2), xst(0,3):xen(0,3))
            do idx = 1,nproc-1
                chunksize = xsz(idx,1)*xsz(idx,2)*xsz(idx,3)
                call MPI_SEND(full_field(xst(idx,1):xen(idx,1), xst(idx,2):xen(idx,2),xst(idx,3):xen(idx,3)), chunksize,&
                mpirkind, idx, tag, MPI_COMM_WORLD, ierr)
            end do

            ! V data 
            fname = dir(:len_trim(dir))//'/'//vFile
            open(fid,file=fname,form='unformatted', access='stream')  
            read(fid) sizes(1:3)
            if ((sizes(1) == Nx) .and. (sizes(2) == Ny) .and. (sizes(3) == Nz)) then
                read(fid) (((full_field(i,j,k),i=1,Nx),j=1,Ny),k=1,Nz)
            else
                close(fid)
                print*, "Problem size      :", Nx, Ny, Nz
                print*, "Hit 3d input size :", sizes 
                call GracefulExit("Problem size mismatch with HIT_3d input file:"//vFile,1002) 
            end if
            tag = 20
            fieldsPhys(:,:,:,2) = full_field(xst(0,1):xen(0,1), xst(0,2):xen(0,2), xst(0,3):xen(0,3))
            do idx = 1,nproc-1
                chunksize = xsz(idx,1)*xsz(idx,2)*xsz(idx,3)
                call MPI_SEND(full_field(xst(idx,1):xen(idx,1), xst(idx,2):xen(idx,2),xst(idx,3):xen(idx,3)), chunksize,&
                mpirkind, idx, tag, MPI_COMM_WORLD, ierr)
            end do

            ! W data 
            fname = dir(:len_trim(dir))//'/'//wFile
            open(fid,file=fname,form='unformatted', access='stream')  
            read(fid) sizes(1:3)
            if ((sizes(1) == Nx) .and. (sizes(2) == Ny) .and. (sizes(3) == Nz)) then
                read(fid) (((full_field(i,j,k),i=1,Nx),j=1,Ny),k=1,Nz)
            else
                close(fid)
                print*, "Problem size      :", Nx, Ny, Nz
                print*, "Hit 3d input size :", sizes 
                call GracefulExit("Problem size mismatch with HIT_3d input file:"//wFile,1003) 
            end if
            tag = 30
            fieldsPhys(:,:,:,3) = full_field(xst(0,1):xen(0,1), xst(0,2):xen(0,2), xst(0,3):xen(0,3))
            do idx = 1,nproc-1
                chunksize = xsz(idx,1)*xsz(idx,2)*xsz(idx,3)
                call MPI_SEND(full_field(xst(idx,1):xen(idx,1), xst(idx,2):xen(idx,2),xst(idx,3):xen(idx,3)), chunksize,&
                mpirkind, idx, tag, MPI_COMM_WORLD, ierr)
            end do

        else
            tag = 10
            chunksize = gp%xsz(1)*gp%xsz(2)*gp%xsz(3)
            call MPI_RECV(fieldsPhys(:,:,:,1), chunksize, mpirkind, 0, tag,&
                    MPI_COMM_WORLD, status, ierr)
            tag = 20
            call MPI_RECV(fieldsPhys(:,:,:,2), chunksize, mpirkind, 0, tag,&
                    MPI_COMM_WORLD, status, ierr)
            tag = 30
            call MPI_RECV(fieldsPhys(:,:,:,3), chunksize, mpirkind, 0, tag,&
                    MPI_COMM_WORLD, status, ierr)
            
        end if 

        if (nrank == 0) then
            deallocate(full_Field)
        end if
    end subroutine 

    subroutine closeMatlabFile
        if (nrank == 0) then
            write(headerfid,*) "--------------------------------------------------------------"
            write(headerfid,*) "------------------ END OF HEADER FILE ------------------------"
            close(headerfid)
        end if 
    end subroutine 
end module 
