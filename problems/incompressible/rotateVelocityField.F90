!module RoutinesUpsampling
!    use kind_parameters, only: rkind
!    
!    implicit none
!    contains
!    subroutine upsampleX(arrIn, arrOut)
!        real(rkind), dimension(:,:,:), intent(in) :: arrIn
!        real(rkind), dimension(:,:,:), intent(out) :: arrOut
!        integer :: i, j, k, nx
!
!        nx = size(arrIn,1)
!        do k = 1,size(arrIn,3)
!            do j = 1,size(arrIn,2)
!                do i = 1,nx
!                    arrOut(2*i-1,j,k) = arrIn(i,j,k)
!                end do 
!            end do 
!        end do 
!
!        do k = 1,size(arrIn,3)
!            do j = 1,size(arrIn,2)
!                do i = 2,2*nx-2,2
!                    arrOut(i,j,k) = 0.5d0*(arrOut(i - 1,j,k) + arrOut(i + 1,j,k))
!                end do 
!                arrOut(2*nx,j,k) = 0.5d0*(arrOut(2*nx-1,j,k) + arrOut(1,j,k))
!            end do 
!        end do 
!
!    end subroutine
!
!    subroutine upsampleY(arrIn, arrOut)
!        real(rkind), dimension(:,:,:), intent(in) :: arrIn
!        real(rkind), dimension(:,:,:), intent(out) :: arrOut
!        integer :: j, k, ny
!
!        ny = size(arrIn,2)
!        do k = 1,size(arrIn,3)
!            do j = 1,ny
!                    arrOut(:,2*j-1,k) = arrIn(:,j,k)
!            end do 
!        end do 
!
!        do k = 1,size(arrIn,3)
!            do j = 2,2*ny-2,2
!                    arrOut(:,j,k) = 0.5d0*(arrOut(:,j - 1,k) + arrOut(:,j + 1,k))
!            end do 
!            arrOut(:,2*ny,k) = 0.5d0*(arrOut(:,2*ny-1,k) + arrOut(:,1,k))
!        end do 
!
!    end subroutine
!
!    subroutine upsampleZ_cells(arrIn,arrOut)
!        real(rkind), dimension(:,:,:), intent(in) :: arrIn
!        real(rkind), dimension(:,:,:), intent(out) :: arrOut
!        integer :: k, nz, idx
!        real(rkind), parameter :: ratEven = 0.25d0, ratOdd = 0.75d0
!
!        nz = size(arrIn,3)
!        do k = 1,nz-1
!            arrOut(:,:,2*k+1) = arrIn(:,:,k) + (arrIn(:,:,k+1) - arrIn(:,:,k))*ratOdd 
!        end do 
!        
!        idx = 1
!        do k = 2,2*nz-1,2
!            arrOut(:,:,k) = arrIn(:,:,idx) + (arrIn(:,:,idx + 1) - arrIn(:,:,idx))*ratEven
!            idx = idx + 1
!        end do
!   
!        arrOut(:,:,1) = 2.d0*arrOut(:,:,2) - arrOut(:,:,3)
!        arrOut(:,:,2*nz) = 2.d0*arrOut(:,:,2*nz - 1) - arrOut(:,:,2*nz - 2)
!
!    end subroutine
!
!
!    subroutine upsampleZ_edges(arrIn,arrOut)
!        real(rkind), dimension(:,:,:), intent(in) :: arrIn
!        real(rkind), dimension(:,:,:), intent(out) :: arrOut
!        integer :: k, nz
!
!        nz = size(arrIn,3) - 1
!
!        do k = 1,nz
!            arrOut(:,:,2*k-1) = arrIn(:,:,k)
!        end do 
!        arrOut(:,:,2*nz+1) = arrOut(:,:,nz+1)
!
!        do k = 1,nz
!            arrOut(:,:,2*k) = 0.5d0*(arrOut(:,:,2*k-1) + arrOut(:,:,2*k+1))
!        end do 
!
!        arrOut(:,:,2*nz+1) = 0.d0
!        arrOut(:,:,1) = 0.d0
!    end subroutine
!
!end module



program rotateVelocityField
    use kind_parameters, only: rkind, clen
    use gridtools, only: alloc_buffs, destroy_buffs
    use decomp_2d
    use exits, only: GracefulExit, message
    use mpi 
    use timer, only: tic, toc
    use decomp_2d_io
    use constants, only: one, half, pi, zero
    !use RoutinesUpsampling
   
    implicit none

    character(len=clen) :: inputfile 
    character(len=clen) :: outputdir, inputdir
    integer :: ioUnit, nx, ny, nz, ierr, statsFile_TID, inputFile_TID, inputFile_RID, outputFile_TID, outputFile_RID, i, j, k, kindex
    type(decomp_info) :: gpC
    real(rkind), dimension(:,:,:), allocatable :: u, v
    real(rkind), dimension(:),     allocatable :: ustats, vstats, zC
    character(len=clen) :: tempname, fname
    real(rkind) :: tsim, utmp, vtmp, rotmatrix(2,2), theta, uhub, vhub, HubHeight, dz, utmp1, utmp2, vtmp1, vtmp2
    namelist /INPUT/ nx, ny, nz, inputdir, outputdir, inputFile_TID, inputFile_RID, &
    outputFile_TID, outputFile_RID, statsFile_TID, HubHeight 

    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file

    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=INPUT)
    close(ioUnit)

    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(gpC)
    
    !!!!!! READ HEADER !!!!!!!!!!!
    write(tempname,"(A7,A4,I2.2,A6,I6.6)") "RESTART", "_Run",inputFile_RID, "_info.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    open(unit=10,file=fname,access='sequential',form='formatted')
    read (10, *)  tsim
    close(10)

    !!!!!!!!!!!!! CELL FIELDS !!!!!!!!!!!!!!!
    allocate(u(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(v(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
    allocate(zC(gpC%xsz(3)))
    allocate(ustats(gpC%xsz(3)))
    allocate(vstats(gpC%xsz(3)))
    
    ! Read u - field
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_u.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,u,fname, gpC)

    ! Read v - field
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",inputFile_RID, "_v.",inputFile_TID
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    call decomp_2d_read_one(1,v,fname, gpC)

    ! determine the angle by which to rotate
    dz = one/real(nz,rkind)
    do k = 1, nz
      zC(k) = (real(k,rkind)-half)*dz
    enddo
    kindex = -1
    do k = 1, nz-1
      if(HubHeight >= zC(k) .and. HubHeight < zC(k+1))  then
        kindex = k
        exit
      endif
    enddo
    if(HubHeight >= zero   .and. HubHeight < zC(1)) kindex = 1
    if(HubHeight >= zC(nz) .and. HubHeight <= one)  kindex = nz-1
    if(kindex==-1) then
      write(*,*) 'HubHeight = ', HubHeight, 'not bracketed in 0 and 1'
      stop
    endif
    
    ! Read u, v stats
    write(tempname,"(A3,I2.2,A2,I6.6,A4)") "Run",inputFile_RID, "_t",statsFile_TID,".stt"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
    open(unit=12,file=fname,status='old',action='read')
    do k = 1, nz
      read(12,*) ustats(k), vstats(k)
    enddo
    close(12)

    uhub = ustats(kindex) + (ustats(kindex+1)-ustats(kindex))/dz*(HubHeight-zC(kindex))
    vhub = vstats(kindex) + (vstats(kindex+1)-vstats(kindex))/dz*(HubHeight-zC(kindex))
    theta = atan2(vhub, uhub)
    write(*,*) 'Based on time and horizontal average'
    write(*,*) 'Velocity components at HubHeight = ', uhub, vhub
    write(*,*) 'Velocity field to be rotated by ', -theta*real(180,rkind)/pi

    utmp1 = sum(u(:,:,kindex))/real(nx*ny, rkind)
    utmp2 = sum(u(:,:,kindex+1))/real(nx*ny, rkind)
    vtmp1 = sum(v(:,:,kindex))/real(nx*ny, rkind)
    vtmp2 = sum(v(:,:,kindex+1))/real(nx*ny, rkind)
    uhub = utmp1 + (utmp2-utmp1)/dz*(HubHeight-zC(kindex))
    vhub = vtmp1 + (vtmp2-vtmp1)/dz*(HubHeight-zC(kindex))
    theta = atan2(vhub, uhub)
    write(*,*) 'Based on instantaneous horizontal average'
    write(*,*) 'Velocity components at HubHeight = ', uhub, vhub
    write(*,*) 'Velocity field to be rotated by ', -theta*real(180,rkind)/pi


    ! setup rotation matrix and rotate velocity vector
    rotmatrix(1,1) = cos(-theta); rotmatrix(1,2) = -sin(-theta)!; rotmatrix(1,3) = zero
    rotmatrix(2,1) = sin(-theta); rotmatrix(2,2) =  cos(-theta)!; rotmatrix(2,3) = zero
   !rotmatrix(3,1) = cos(-theta); rotmatrix(3,2) = -sin(-theta)!; rotmatrix(3,3) = zero

    do k = 1, nz
     do j = 1, ny
      do i = 1, nx
        utmp = rotmatrix(1,1)*u(i,j,k) + rotmatrix(1,2)*v(i,j,k)
        vtmp = rotmatrix(2,1)*u(i,j,k) + rotmatrix(2,2)*v(i,j,k)
        u(i,j,k) = utmp
        v(i,j,k) = vtmp
      enddo
     enddo
    enddo

    ! Write out new restart files
    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_u.",outputFile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
    call decomp_2d_write_one(1,u,fname, gpC)

    write(tempname,"(A7,A4,I2.2,A3,I6.6)") "RESTART", "_Run",outputFile_RID, "_v.",outputFile_TID
    fname = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
    call decomp_2d_write_one(1,v,fname, gpC)

    utmp1 = sum(u(:,:,kindex))/real(nx*ny, rkind)
    utmp2 = sum(u(:,:,kindex+1))/real(nx*ny, rkind)
    vtmp1 = sum(v(:,:,kindex))/real(nx*ny, rkind)
    vtmp2 = sum(v(:,:,kindex+1))/real(nx*ny, rkind)
    uhub = utmp1 + (utmp2-utmp1)/dz*(HubHeight-zC(kindex))
    vhub = vtmp1 + (vtmp2-vtmp1)/dz*(HubHeight-zC(kindex))
    theta = atan2(vhub, uhub)
    write(*,*) 'After rotation'
    write(*,*) 'Velocity components at HubHeight = ', uhub, vhub
    write(*,*) 'Velocity vector angle with x axis: ', -theta*real(180,rkind)/pi

    deallocate(u,v, zC, ustats, vstats)

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
