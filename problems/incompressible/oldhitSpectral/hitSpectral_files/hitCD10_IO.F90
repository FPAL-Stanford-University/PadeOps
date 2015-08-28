module hitCD_IO

    use kind_parameters, only: rkind
    use decomp_2d,       only: decomp_info,nrank,nproc
    use exits,           only: GracefulExit, message
    implicit none

    integer, dimension(:,:), allocatable        :: xst,xen,xsz

contains

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

end module 
