subroutine getNeighbor(ist,ien,jst,isPeriodic,whichNeighbor,neighbor)
  use findMod, only: find
  integer, dimension(:), intent(in) :: ist, ien, jst ! Start and end
                                                          ! indices for all processes
  integer, intent(out) :: neighbor
  integer, intent(in) :: whichNeighbor
  logical, intent(in) :: isPeriodic
  integer, dimension(:), allocatable :: idx1, idx2
  integer, dimension(:), allocatable :: PEpossible
  integer, dimension(nproc) :: allPE
  integer :: ierr, i, j, n
  
  ! Array containing ID of all MPI ranks
  do n = 1,nproc
    allPE(n) = n-1
  end do

  case1: select case (whichNeighbor)
  case (1) ! X-neighbor to right (i.e. "High" neighbor)
    !------------------------------------------------------------------------
    if (ien(nrank+1) == maxval(ien) .and. ist(nrank+1) == minval(ist)) then ! No neighbor
      neighbor = -1
      exit case1
    elseif (ien(nrank+1) == maxval(ien) .and. isPeriodic) then ! Periodic neighbor
      call find(ist,1,idx1)
    elseif (ien(nrank+1) == maxval(ien)) then              ! There is no neighbor
      neighbor = -1
      exit case1
    else                                                   ! Sequential neighbor
      call find(ist,ien(nrank+1)+1,idx1)
    end if

    allocate(PEpossible(size(idx1)))
    PEpossible = allPE(idx1)

    call find(jst,jst(nrank+1),idx2)
    do i = 1,size(idx2)
      do j = 1,size(PEpossible)
        if (PEpossible(j) == allPE(idx2(i))) then
          neighbor = PEpossible(j)
          exit
        end if
      end do
    end do

  case (-1) ! X-neighbor to left (i.e. "Low" neighbor)
    !----------------------------------------------------------------------------
    if (ien(nrank+1) == maxval(ien) .and. ist(nrank+1) == minval(ist)) then ! No neighbor
      neighbor = -1
      exit case1
    elseif (ist(nrank+1) == minval(ist) .and. isPeriodic) then ! Periodic neighbor
      call find(ien,maxval(ien),idx1)
    elseif (ist(nrank+1) == minval(ist)) then              ! There is no neighbor
      neighbor = -1
      exit case1
    else                                          ! Sequential neighbor
      call find(ien,ist(nrank+1)-1,idx1)
    end if

    allocate(PEpossible(size(idx1)))
    PEpossible = allPE(idx1)

    call find(jst,jst(nrank+1),idx2)
    do i = 1,size(idx2)
      do j = 1,size(PEpossible)
        if (PEpossible(j) == allPE(idx2(i))) then
          neighbor = PEpossible(j)
          exit case1
        end if
      end do
    end do

  case default
    !--------------------------------------------------------------------------
    call gracefulExit("Must specify upwind (1) or downwind (-1) neighbor --"//&
      " GaborModeRoutines.F90",ierr)
  end select case1
  
  if (allocated(idx1)) deallocate(idx1)
  if (allocated(idx2)) deallocate(idx2)
  if (allocated(PEpossible)) deallocate(PEpossible)
end subroutine

subroutine mpiIsendIrecv3D(sendBuff, recvBuff, sendNeigh, recvNeigh, sendReq,&
    recvReq, tag)
  use kind_parameters, only: mpirkind
  real(rkind), dimension(:,:,:), intent(in) :: sendBuff
  real(rkind), dimension(:,:,:), intent(out) :: recvBuff
  integer, intent(in) :: sendNeigh, recvNeigh
  integer, intent(out) :: sendReq, recvReq
  integer, intent(in), optional :: tag
  integer :: ierr, mytag

  mytag = 1
  if (present(tag)) mytag = tag
  call MPI_Irecv(recvBuff, size(recvBuff), mpirkind, recvNeigh, mytag, &
    MPI_COMM_WORLD, recvReq, ierr)
  call assert(ierr == MPI_SUCCESS)
  call MPI_Isend(sendBuff, size(sendBuff), mpirkind, sendNeigh, mytag, &
    MPI_COMM_WORLD, sendReq, ierr)
  call assert(ierr == MPI_SUCCESS)
end subroutine

subroutine mpiIsendIrecv2D(sendBuff, recvBuff, sendNeigh, recvNeigh, sendReq,&
    recvReq, tag)
  use kind_parameters, only: mpirkind
  real(rkind), dimension(:,:), intent(in) :: sendBuff
  real(rkind), dimension(:,:), intent(out) :: recvBuff
  integer, intent(in) :: sendNeigh, recvNeigh
  integer, intent(out) :: sendReq, recvReq
  integer, intent(in), optional :: tag
  integer :: ierr, mytag

  mytag = 1
  if (present(tag)) mytag = tag

  call MPI_Irecv(recvBuff, size(recvBuff), mpirkind, recvNeigh, mytag, &
    MPI_COMM_WORLD, recvReq, ierr)
  call assert(ierr == MPI_SUCCESS)
  call MPI_Isend(sendBuff, size(sendBuff), mpirkind, sendNeigh, mytag, &
    MPI_COMM_WORLD, sendReq, ierr)
  call assert(ierr == MPI_SUCCESS)
end subroutine

subroutine haloExchangeMPI(u)
  use domainSetup, only: decomp2Dpencil
  real(rkind), dimension(:,:,:), intent(inout) :: u
  integer :: ierr
  integer :: isz, jsz, ksz, neighHi, neighLo, tag

  isz = size(u,1)
  jsz = size(u,2)
  ksz = size(u,3)

  select case(decomp2Dpencil)
  case('x')

    ! Now do data exchange in x-direction
    if (periodic(1)) then
      u(1+nxsupp/2:nxsupp,:,:) = u(1+nxsupp/2:nxsupp,:,:) + &
        u(isz-nxsupp/2+1:isz,:,:)
      u(isz-nxsupp+1:isz-nxsupp/2,:,:) = u(isz-nxsupp+1:isz-nxsupp/2,:,:) + &
        u(1:nxsupp/2,:,:)
    end if
    
    ! Now do the same in y-direction
    call getNeighbor(jstAll,jenAll,kstAll,periodic(2),1,neighHi)   
    call getNeighbor(jstAll,jenAll,kstAll,periodic(2),-1,neighLo)
    
    if (neighHi == -1 .and. periodic(2)) then
      u(:,1+nysupp/2:nysupp,:) = u(:,1+nysupp/2:nysupp,:) + &
        u(:,jsz-nysupp/2+1:jsz,:)
      u(:,jsz-nysupp+1:jsz-nysupp/2,:) = u(:,jsz-nysupp+1:jsz-nysupp/2,:) + &
        u(:,1:nysupp/2,:)
    elseif (neighHi == -1 .and. neighLo == -1) then
      continue
    elseif (neighHi == -1) then ! Not periodic in y
      call assert(.false., "Hasn't been implemented yet! -- haloExchangeMPI()"//&
        "in GaborModeRoutines.F90")
    elseif (neighLo == -1) then ! Not periodic in y
      call assert(.false., "Hasn't been implemented yet! -- haloExchangeMPI()"//&
        "in GaborModeRoutines.F90")
    else
      ! Loop over each z-station so MPI send/recv contiguous chunks of memory
      do tag = 1,ksz
        call mpiIsendIrecv(u(:,jsz-nysupp/2+1:jsz,tag), buff3(:,:,tag), &
          neighHi, neighLo, sendReqHi(tag), recvReqLo(tag), tag)
        call mpiIsendIrecv(u(:,1:nysupp/2,tag), buff4(:,:,tag), &
         neighLo, neighHi, sendReqLo(tag), recvReqHi(tag), tag)
      end do

      ! Wait for data traffic and modify velocity
      call MPI_Waitall(2*ksz,(/sendReqHi, recvReqLo/),MPI_STATUSES_IGNORE,ierr)
      u(:,1+nysupp/2:nysupp,:) = u(:,1+nysupp/2:nysupp,:) + buff3
      
      call MPI_Waitall(2*ksz,(/sendReqLo, recvReqHi/),MPI_STATUSES_IGNORE,ierr)
      u(:,jsz-nysupp+1:jsz-nysupp/2,:) = u(:,jsz-nysupp+1:jsz-nysupp/2,:) + buff4 
    end if

    ! Exchange velocity with z-neighbors
    call getNeighbor(kstAll,kenAll,jstAll,periodic(3),1,neighHi)   
    call getNeighbor(kstAll,kenAll,jstAll,periodic(3),-1,neighLo)
    
    if (neighHi == -1 .and. periodic(3)) then ! Periodic halo exchange in local memory
      u(:,:,1+nzsupp/2:nzsupp) = u(:,:,1+nzsupp/2:nzsupp) + &
        u(:,:,ksz-nzsupp/2+1:ksz)
      u(:,:,ksz-nzsupp+1:ksz-nzsupp/2) = u(:,:,ksz-nzsupp+1:ksz-nzsupp/2) + &
        u(:,:,1:nzsupp/2)
    elseif (neighHi == -1 .and. neighLo == -1) then ! No halo exchange
      continue
    elseif (neighHi == -1) then ! Send and receive data to lower neighbor only
      call mpiIsendIrecv(u(:,:,1:nzsupp/2),buff2,neighLo,neighLo,sendReqLo(1), recvReqHi(1))
      call MPI_Waitall(2,(/sendReqLo(1), recvReqHi(1)/),MPI_STATUSES_IGNORE,ierr)
      u(:,:,1+nzsupp/2:nzsupp) = u(:,:,1+nzsupp/2:nzsupp) + buff2

    elseif (neighLo == -1) then ! Send and receive data to higher neighbor only
      call mpiIsendIrecv(u(:,:,ksz-nzsupp/2+1:ksz),buff1,neighHi,neighHi,sendReqHi(1), recvReqLo(1))
      call MPI_Waitall(2,(/sendReqHi(1), recvReqLo(1)/),MPI_STATUSES_IGNORE,ierr)
      u(:,:,ksz-nzsupp+1:ksz-nzsupp/2) = u(:,:,ksz-nzsupp+1:ksz-nzsupp/2) + buff1
      
    else ! Send data to neighboring processes
      ! Data movement from low to high: lo -> me -> hi
      call mpiIsendIrecv(u(:,:,ksz-nzsupp/2+1:ksz),buff1,neighHi,neighLo,sendReqHi(1), recvReqLo(1))
      
      ! Data movement from high to low: hi -> me -> lo
      call mpiIsendIrecv(u(:,:,1:nzsupp/2),buff2,neighLo,neighHi,sendReqLo(1), recvReqHi(1))
 
      ! Wait for data traffic and modify velocity
      call MPI_Waitall(2,(/sendReqLo(1), recvReqLo(1)/),MPI_STATUSES_IGNORE,ierr)
      u(:,:,1+nzsupp/2:nzsupp) = u(:,:,1+nzsupp/2:nzsupp) + buff1
      
      call MPI_Waitall(2,(/sendReqHi(1), recvReqHi(1)/),MPI_STATUSES_IGNORE,ierr)
      u(:,:,ksz-nzsupp+1:ksz-nzsupp/2) = u(:,:,ksz-nzsupp+1:ksz-nzsupp/2) + buff2
    end if
  case('y')
    call assert(.false.,"Haven't implemented y-pencil yet")
  case('z')
    call assert(.false.,"Haven't implemented z-pencil yet")
  case default
    call gracefulExit('decomp2Dpencil must be one of "x", "y", or "z" --"//&
      " GaborModeRoutines.F90',ierr)
  end select

end subroutine

