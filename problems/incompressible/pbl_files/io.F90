module pbl_IO

    use kind_parameters, only: rkind, clen
    use decomp_2d,       only: decomp_info,nrank,nproc
    use exits,           only: GracefulExit, message
    implicit none

    integer, dimension(:,:), allocatable        :: xst,xen,xsz
    integer :: headerfid = 101
    integer :: NumDumps 

contains

    subroutine start_io(gp)
        use IncompressibleGridNP, only: igrid
        use mpi 
        
        class(igrid), target, intent(in) :: gp 
        character(len=clen) :: fname
        character(len=clen) :: tempname
        !character(len=clen) :: command
        character(len=clen) :: OutputDir
        !integer :: system 
        integer :: runIDX 
        logical :: isThere
        integer :: tag, idx, status(MPI_STATUS_SIZE), ierr

        ! Create data sharing info
        if (nrank == 0) then
            allocate(xst(0:nproc-1,3),xen(0:nproc-1,3),xsz(0:nproc-1,3))
        end if


        ! communicate local processor grid info (Assume x-decomposition)
        if (nrank == 0) then
            xst(0,:) = gp%gpC%xst
            xen(0,:) = gp%gpC%xen
            
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
            call MPI_SEND(gp%gpC%xst, 3, MPI_INTEGER, 0, tag, &
                 &      MPI_COMM_WORLD, ierr)
            tag = 1
            call MPI_SEND(gp%gpC%xen, 3, MPI_INTEGER, 0, tag, &
                 &      MPI_COMM_WORLD, ierr)
            tag = 2
            call MPI_SEND(gp%gpC%xsz, 3, MPI_INTEGER, 0, tag, &
                 &      MPI_COMM_WORLD, ierr)

        end if 

        OutputDir = gp%outputdir
        runIDX = gp%runID
        
        inquire(FILE=trim(OutputDir), exist=isThere)
        if (nrank == 0) then
            !if (.not. isThere) then
            !    print*, "============================================="
            !    print*, "WARNING: Output directory not found. Creating a new one."
            !    print*, "============================================="
            !    command = "mkdir "//trim(OutputDir)
            !    ierr = system(trim(command))
            !end if 
            write(tempname,"(A3,I2.2,A6,A4)") "Run", runIDX, "HEADER",".txt"
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

        ! Now perform the initializing data dump
        call dumpData4Matlab(gp)
    end subroutine

    subroutine dumpData4Matlab(gp)
        use IncompressibleGridNP, only: igrid
        use gridtools,          only: alloc_buffs
        use decomp_2d,        only: transpose_y_to_x
        
        class(igrid), target, intent(in) :: gp 
        integer :: tid, runIDX
        character(len=clen) :: fname
        character(len=clen) :: tempname
        character(len=clen) :: OutputDir
        real(rkind), dimension(:,:,:,:), pointer :: fieldsPhys
        integer :: fid = 1234

        OutputDir = gp%outputdir
        fieldsPhys => gp%PfieldsC
        runIDX = gp%runID
        tid = gp%step 

        write(tempname,"(A3,I2.2,A2,I4.4,A2,I6.6,A5,A4)") "Run", RunIDX, "_p",nrank,"_t",tid,"_uVEL",".out"
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
        open(fid,file=trim(fname),form='unformatted',status='replace')
        write(fid) fieldsPhys(:,:,:,1)
        close(fid)

        write(tempname,"(A3,I2.2,A2,I4.4,A2,I6.6,A5,A4)") "Run", RunIDX, "_p",nrank,"_t",tid,"_vVEL",".out"
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
        open(fid,file=trim(fname),form='unformatted',status='replace')
        write(fid) fieldsPhys(:,:,:,2)
        close(fid)

        write(tempname,"(A3,I2.2,A2,I4.4,A2,I6.6,A5,A4)") "Run", RunIDX, "_p",nrank,"_t",tid,"_wVEL",".out"
        fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
        open(fid,file=trim(fname),form='unformatted',status='replace')
        write(fid) fieldsPhys(:,:,:,3)
        close(fid)
        
        if (nrank == 0) then
            write(headerfid,"(I8)") tid
        end if  
        NumDumps = NumDumps + 1

        nullify(fieldsPhys)
    end subroutine 

    subroutine finalize_io
        if (nrank == 0) then
            write(headerfid,*) "--------------------------------------------------------------"
            write(headerfid,*) "------------------ END OF HEADER FILE ------------------------"
            close(headerfid)
        end if 
    end subroutine 
end module 
