module io
    use kind_parameters, only: rkind
    use decomp_2d
    use exits, only: GracefulExit 
    use variables
    
    implicit none
    

    integer :: fileinput  = 1
    integer :: fileoutput = 2

    integer, dimension(:,:), allocatable        :: xst,xen,xsz
    character(len=11) :: eof
    integer :: NumDumps 
    integer :: headerfid = 101
contains
    subroutine get_input(inputFile)
        character(len=*), intent(in) :: inputFile
        logical                      :: isThere
         

        inquire(file=trim(inputFile), exist=isThere)
        if (.not. isThere) then
            call GracefulExit("Input file:" // trim(inputFile)// "is missing",1001)
        end if 

        open(fileinput,file=trim(inputFile),form="formatted")
        
        include "hitSpectral_files/inputReadOrder.F90"
        
        close(fileinput)
      
        ! Fix all the string inputs  
        call fixIOstring(HIT3dinitDir)
        call fixIOstring(RESTARTDir)
        call fixIOstring(OutputDir)

        if (eof == "END OF FILE") then
            return         
        end if

2000    continue
        call GracefulExit("Error occured while reading the input file", 1002)
        

    end subroutine 

    subroutine write_matlab_header
        use mpi 
        character(len=256) :: fname
        character(len=32) :: tempname
        character(len=32) :: command
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
                ierr = system(command)
            end if 
            write(tempname,"(A3,I2.2,A6,A4)") "Run", RunIDX, "HEADER",".txt"
            fname = OutputDir(:len_trim(OutputDir))//trim(tempname)
            
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

    subroutine dumpData4Matlab(tid) 
        integer, intent(in) :: tid
        character(len=256) :: fname
        character(len=32) :: tempname
        integer :: fid = 1234

        write(tempname,"(A3,I2.2,A2,I4.4,A2,I5.5,A5,A4)") "Run", RunIDX, "_p",nrank,"_t",tid,"_uVEL",".out"
        fname = OutputDir(:len_trim(OutputDir))//trim(tempname)
        open(fid,file=trim(fname),form='unformatted',status='replace')
        write(fid) fieldsPhys(:,:,:,1)
        close(fid)

        write(tempname,"(A3,I2.2,A2,I4.4,A2,I5.5,A5,A4)") "Run", RunIDX, "_p",nrank,"_t",tid,"_vVEL",".out"
        fname = OutputDir(:len_trim(OutputDir))//trim(tempname)
        open(fid,file=trim(fname),form='unformatted',status='replace')
        write(fid) fieldsPhys(:,:,:,2)
        close(fid)

        write(tempname,"(A3,I2.2,A2,I4.4,A2,I5.5,A5,A4)") "Run", RunIDX, "_p",nrank,"_t",tid,"_wVEL",".out"
        fname = OutputDir(:len_trim(OutputDir))//trim(tempname)
        open(fid,file=trim(fname),form='unformatted',status='replace')
        write(fid) fieldsPhys(:,:,:,3)
        close(fid)
        
        if (nrank == 0) then
            write(headerfid,"(I8)") tid
        end if  
        NumDumps = NumDumps + 1
    end subroutine 

    subroutine closeMatlabFile
        if (nrank == 0) then
            write(headerfid,*) "--------------------------------------------------------------"
            write(headerfid,*) "------------------ END OF HEADER FILE ------------------------"
            close(headerfid)
        end if 
    end subroutine 
    
    subroutine createRestartFile 


    end subroutine

    subroutine fixIOstring(string)
        character(len=*), intent(inout) :: string
        integer :: trim_idx

        trim_idx = index(string,"|") 
        string = trim(string(1:trim_idx-1))
    end subroutine 

    subroutine getHit3d_uvw(dir)
        use mpi 
        use kind_parameters, only: mpirkind 
        character(len=*), intent(in)                :: dir
        character(len=8),parameter                  :: uFile = "U.000000"
        character(len=8),parameter                  :: vFile = "V.000000"
        character(len=8),parameter                  :: wFile = "W.000000"
        character(:), allocatable                   :: fname
        real(rkind), dimension(:,:,:), allocatable  :: full_Field
        integer :: fid = 1234, i, j, k
        integer :: tag, idx, status(MPI_STATUS_SIZE), ierr
        integer :: sizes(3), chunksize

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
            fname = dir(:len_trim(dir))//uFile
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
            fname = dir(:len_trim(dir))//vFile
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
            fname = dir(:len_trim(dir))//wFile
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
 
end module io
 
