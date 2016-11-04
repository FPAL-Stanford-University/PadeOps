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
        use IncompressibleGridWallM, only: igridWallM
        use mpi 
        
        class(igridWallM), target, intent(in) :: gp 
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
        !call dumpData4Matlab(gp)
        !call gp%dumpFullField(gp%u,'uVel')
        !call gp%dumpFullField(gp%v,'vVel')
        !call gp%dumpFullField(gp%wC,'wVel')
        !call output_tecplot(gp)
    end subroutine

    subroutine dumpData4Matlab(gp)
        use IncompressibleGridWallM, only: igridWallM
        use gridtools,          only: alloc_buffs
        use decomp_2d,        only: transpose_y_to_x
        
        class(igridWallM), target, intent(in) :: gp 
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

    !subroutine output_tecplot(gp)
    !    use IncompressibleGridWallM, only: igridWallM
    !    use gridtools,          only: alloc_buffs
    !    use decomp_2d,        only: transpose_y_to_x
    !    use turbineMod, only: turbineArray
 
    !    class(igridWallM), target, intent(in) :: gp 
    !    integer :: tid, runIDX
    !    character(len=clen) :: fname
    !    character(len=clen) :: tempname
    !    character(len=clen) :: OutputDir
    !    real(rkind), dimension(:,:,:,:), pointer :: fieldsPhys, xyz
    !    type(turbineArray), pointer :: turbarr
    !    integer :: fid = 1234, turbID, blID, ptID, i, j, k, T_indx

    !    OutputDir = gp%outputdir
    !    fieldsPhys => gp%PfieldsC
    !    xyz => gp%mesh
    !    runIDX = gp%runID
    !    tid = gp%step 

    !    ! output field variables
    !    if(gp%isStratified) then
    !      T_indx = 7
    !    else
    !      T_indx = 3
    !    endif
    !    write(tempname,"(A4,I2.2,A2,I4.4,A4)") "tec_", RunIDX, "_p",nrank,".dat"
    !    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
    !    if(tid==0) then
    !      open(fid,file=trim(fname),status='replace')
    !      write(fid,'(75a)') 'VARIABLES="x","y","z","u","v","wC","T"'
    !      write(fid,'(6(a,i7),a)') 'ZONE I=', gp%gpC%xsz(1), ' J=', gp%gpC%xsz(2), ' K=', gp%gpC%xsz(3), ' ZONETYPE=ORDERED'
    !      write(fid,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', gp%tsim
    !      do k = 1, gp%gpC%xsz(3)
    !       do j = 1, gp%gpC%xsz(2)
    !        do i = 1, gp%gpC%xsz(1)
    !            write(fid,'(7ES26.16)') xyz(i,j,k,1:3), fieldsPhys(i,j,k,1:3), fieldsPhys(i,j,k,T_indx)
    !        enddo
    !       enddo
    !      enddo
    !      close(fid)
    !    else
    !      open(fid,file=trim(fname),status='old',action='write',position='append')
    !      write(fid,'(6(a,i7),a)') 'ZONE I=', gp%gpC%xsz(1), ' J=', gp%gpC%xsz(2), ' K=', gp%gpC%xsz(3), ' ZONETYPE=ORDERED'
    !      write(fid,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', gp%tsim
    !      write(fid,'(a)') ' VARSHARELIST=([1, 2, 3]=1)'
    !      do k = 1, gp%gpC%xsz(3)
    !       do j = 1, gp%gpC%xsz(2)
    !        do i = 1, gp%gpC%xsz(1)
    !            write(fid,'(4ES26.16)') fieldsPhys(i,j,k,1:3), fieldsPhys(i,j,k,T_indx)
    !        enddo
    !       enddo
    !      enddo
    !      close(fid)
    !    endif

    !    ! output field variables
    !    if(gp%useWindTurbines) then
    !      turbarr => gp%WindTurbineArr
    !      do turbID = 1, turbarr%nTurbines
    !        if(turbarr%num_cells_cloud(turbID) > 0) then
    !          write(tempname,"(A4,I2.2,A2,I4.4,A3,I2.2,A4)") "tec_", RunIDX,"_p",nrank,"_wt",turbID,".dat"
    !          fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
    !          if(tid==0) then
    !            open(fid,file=trim(fname),status='replace')
    !            write(fid,'(110a)') 'VARIABLES="x","y","z","blade_forcex", "blade_forcey" "blade_forcez"'
    !          else
    !            open(fid,file=trim(fname),status='old',action='write',position='append')
    !          endif
    !          write(fid,'(6(a,i7),a)') 'ZONE I=', turbarr%num_blades(turbID)*turbarr%num_blade_points(turbID), ' J=', 1, ' K=', 1, ' ZONETYPE=ORDERED'
    !          write(fid,'(a,ES26.16)') 'DATAPACKING=POINT, SOLUTIONTIME=', gp%tsim
    !          do blID = 1, turbarr%num_blades(turbID)
    !            do ptID = 1, turbarr%num_blade_points(turbID)
    !               write(fid,'(7ES26.16)') turbarr%blade_points(:, ptID, blID, turbID), turbarr%blade_forces(:, ptID, blID, turbID)
    !            enddo
    !          enddo
    !          close(fid)
    !        endif
    !      enddo
    !      nullify(turbarr)
    !    endif

    !    nullify(fieldsPhys,xyz)

    !end subroutine

    subroutine finalize_io
        if (nrank == 0) then
            write(headerfid,*) "--------------------------------------------------------------"
            write(headerfid,*) "------------------ END OF HEADER FILE ------------------------"
            close(headerfid)
        end if 
    end subroutine 
end module 
