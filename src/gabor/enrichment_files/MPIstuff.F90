subroutine sendRecvHaloModes(this,recvBuff)
  class(enrichmentOperator), intent(inout) :: this
  real(rkind), dimension(:,:), allocatable, intent(inout) :: recvBuff
  integer :: Nhalo

  call this%howManyHaloModes(Nhalo)

  if (allocated(recvBuff)) then
    if (size(recvBuff,1) < Nhalo) then
      deallocate(recvBuff)
      allocate(recvBuff(Nhalo + 100,12))
    end if
  else
    allocate(recvBuff(Nhalo+100,12))
  end if



end subroutine

subroutine howManyHaloModes(this,Nhalo)
  use decomp_2d, only: DECOMP_2D_COMM_CART_X
  class(enrichmentOperator), intent(inout), target :: this
  integer, intent(out) :: Nhalo
  integer :: n
  real(rkind) :: yhaloMin, yhaloMax, zhaloMin, zhaloMax
  real(rkind), dimension(:,:), pointer :: ySendBuff1,  ySendBuff2, &
    zSendBuff1, zSendBuff2, yRecvBuff1, yRecvBuff2, zRecvBuff1, zRecvBuff2
  logical, dimension(3) :: periodic
  integer :: y1SendCount, y2SendCount, z1SendCount, z2SendCount
  integer :: y1RecvCount, y2RecvCount, z1RecvCount, z2RecvCount
  integer :: y1idx, y2idx, z1idx, z2idx
  integer :: st, en, ierr
  integer :: tagY1, tagY2, tagZ1, tagZ2
  integer, dimension(96) :: sendReq, recvReq
  real(rkind) :: Ly, Lz
  real(rkind) :: tol = 1.d-12
  real(rkind) :: PEymin, PEymax, PEzmin, PEzmax

  Ly = yDom(2) - yDom(1)
  Lz = zDom(2) - zDom(1)

  PEymin = this%QHgrid%yE(1)
  PEymax = this%QHgrid%yE(this%QHgrid%gpC%xsz(2)+1)
  PEzmin = this%QHgrid%zE(1)
  PEzmax = this%QHgrid%zE(this%QHgrid%gpC%xsz(3)+1)

  ! Define the cutoff location for halo modes
  yhaloMin = this%Qhgrid%yE(1) + &
    & real(this%nysupp/2,rkind)*this%smallScales%dy
  yhaloMax = this%Qhgrid%yE(this%QHgrid%gpC%xsz(2)+1) - &
    & real(this%nysupp/2,rkind)*this%smallScales%dy
  
  zhaloMin = this%Qhgrid%zE(1) + &
    & real(this%nzsupp/2,rkind)*this%smallScales%dz
  zhaloMax = this%Qhgrid%zE(this%QHgrid%gpC%xsz(3)+1) - &
    & real(this%nzsupp/2,rkind)*this%smallScales%dz

  ! Find all modes inside halo regions
  do n = 1,this%nmodes
    if (this%y(n) < yhaloMin) then
      y1SendCount = y1SendCount + 1
    else if (this%y(n) > yhaloMax) then
      y2SendCount = y2SendCount + 1
    end if

    if (this%z(n) < zhaloMin) then
      z1SendCount = z1SendCount + 1
    else if (this%z(n) > zhaloMax) then
      z2SendCount = z2SendCount + 1
    end if
  end do

  if (allocated(this%sendBuff)) then
    if (size(this%sendBuff,1) < y1SendCount + y2SendCount + z1SendCount + &
      z2SendCount) then
      deallocate(this%sendBuff)
      allocate(this%sendBuff(y1SendCount + y2SendCount + z1SendCount + &
        z2SendCount + 100,12))
    end if
  else
    allocate(this%sendBuff(y1SendCount + y2SendCount + z1SendCount + &
      z2SendCount + 100,12))
  end if

  st = 1
  en = y1SendCount
  ySendBuff1 => this%sendBuff(st:en,:)
  
  st = en + 1
  en = en + y2SendCount
  ySendBuff2 => this%sendBuff(st:en,:)
  
  st = en + 1
  en = en + z1SendCount
  zSendBuff1 => this%sendBuff(st:en,:)
  
  st = en + 1
  en = en + z2SendCount
  zSendBuff2 => this%sendBuff(st:en,:)

  ! Pack the send buffers
  ! We must check for two different cases:
  ! --> Cases 1: is the mode within the halo region?
  ! --> Cases 2: If so, does the halo-region correspond to a periodic
  !              boundary?
  ! If Case 2 is true then we must shift the mode location by one
  ! domain length for the velocity rendering routine to work
  y1idx = 1; y2idx = 1; z1idx = 1; z2idx = 1
  do n = 1,this%nmodes
    if (this%y(n) < yhaloMin .and. abs(PEymin - yDom(1)) < tol) then
      
      call this%copyMode(ySendBuff1(y1idx,:),n,[0.d0,Ly,0.d0])
      y1idx = y1idx + 1
    else if (this%y(n) < yhaloMin) then
      call this%copyMode(ySendBuff1(y1idx,:),n)
      y1idx = y1idx + 1
    else if (this%y(n) > yhaloMax .and. abs(PEymax - yDom(2)) < tol) then
      call this%copyMode(ySendBuff2(y2idx,:),n,[0.d0,-Ly,0.d0])
      y2idx = y2idx + 1
    else if (this%y(n) > yhaloMax) then
      call this%copyMode(ySendBuff2(y2idx,:),n)
      y2idx = y2idx + 1
    end if

    if (this%z(n) < zhaloMin .and. abs(PEzmin - zDom(1)) < tol) then
      call this%copyMode(zSendBuff1(z1idx,:),n,[0.d0,0.d0,Lz])
      z1idx = z1idx + 1
    else if (this%z(n) < zhaloMin) then
      call this%copyMode(zSendBuff1(z1idx,:),n)
      z1idx = z1idx + 1
    else if (this%z(n) > zhaloMax .and. abs(PEzmax - zDom(2)) < tol) then
      call this%copyMode(zSendBuff2(z2idx,:),n,[0.d0,0.d0,-Lz])
      z2idx = z2idx + 1
    else if (this%z(n) > zhaloMax) then
      call this%copyMode(zSendBuff2(z2idx,:),n)
      z2idx = z2idx + 1
    end if
  end do

  ! Let your neighbors know how much data you will be sending
  tagY1 = coords(1)
  if (coords(1) == dims(1)-1 .and. periodicBCs(2)) then
    tagY2 = 0
  else
    tagY2 = coords(1) + 1
  end if
  
  tagZ1 = coords(2)
  if (coords(2) == dims(2)-1 .and. periodicBCs(3)) then
    tagZ2 = 0
  else
    tagZ2 = coords(2) + 1
  end if
  
  ! Receive from lower y-rank
  call MPI_IRecv(y1RecvCount,1,MPI_INTEGER,neighbour(1),tagY1,&
    DECOMP_2D_COMM_CART_X,recvReq(1),ierr)
  ! Receive from higher y-rank
  call MPI_IRecv(y2RecvCount,1,MPI_INTEGER,neighbour(2),tagY2,&
    DECOMP_2D_COMM_CART_X,recvReq(2),ierr)
  ! Receive from lower z-rank
  call MPI_IRecv(z1RecvCount,1,MPI_INTEGER,neighbour(3),tagZ1,&
    DECOMP_2D_COMM_CART_X,recvReq(3),ierr)
  ! Receive from higher z-rank
  call MPI_IRecv(z2RecvCount,1,MPI_INTEGER,neighbour(4),tagZ2,&
    DECOMP_2D_COMM_CART_X,recvReq(4),ierr)
  
  ! Send to lower y-rank
  call MPI_ISend(y1SendCount,1,MPI_INTEGER,neighbour(1),tagY1,&
    DECOMP_2D_COMM_CART_X,sendReq(1),ierr)
  ! Send to higher y-rank
  call MPI_ISend(y2SendCount,1,MPI_INTEGER,neighbour(2),tagY2,&
    DECOMP_2D_COMM_CART_X,sendReq(2),ierr)
  ! Send to lower z-rank
  call MPI_ISend(z1SendCount,1,MPI_INTEGER,neighbour(3),tagZ1,&
    DECOMP_2D_COMM_CART_X,sendReq(3),ierr)
  ! Send to higher z-rank
  call MPI_ISend(z2SendCount,1,MPI_INTEGER,neighbour(4),tagZ2,&
    DECOMP_2D_COMM_CART_X,sendReq(4),ierr)

  call MPI_WaitAll(8,[recvReq(1:4), sendReq(1:4)],MPI_STATUSES_IGNORE,ierr)

  st = 1
  en = y1RecvCount
  yRecvBuff1 => this%haloBuff(st:en,:)

  st = en + 1
  en = en + y2RecvCount
  yRecvBuff2 => this%haloBuff(st:en,:)

  st = en + 1
  en = en + z1RecvCount
  zRecvBuff1 => this%haloBuff(st:en,:)

  st = en + 1
  en = en + z2RecvCount
  zRecvBuff2 => this%haloBuff(st:en,:)
  
  ! Send the mode info to neighbors
  do n = 1,12
    call MPI_IRecv(yRecvBuff1(:,n),y1RecvCount,mpirkind,neighbour(1),tagY1,&
      DECOMP_2D_COMM_CART_X,recvReq(n),ierr)
    call MPI_IRecv(yRecvBuff2(:,n),y2RecvCount,mpirkind,neighbour(2),tagY2,&
      DECOMP_2D_COMM_CART_X,recvReq(n+12*1),ierr)
    call MPI_IRecv(zRecvBuff1(:,n),z1RecvCount,mpirkind,neighbour(3),tagZ1,&
      DECOMP_2D_COMM_CART_X,recvReq(n+12*2),ierr)
    call MPI_IRecv(zRecvBuff2(:,n),z2RecvCount,mpirkind,neighbour(4),tagZ2,&
      DECOMP_2D_COMM_CART_X,recvReq(n+12*3),ierr)
  
    call MPI_ISend(ySendBuff1(:,n),y1SendCount,mpirkind,neighbour(1),tagY1,&
     DECOMP_2D_COMM_CART_X,sendReq(n),ierr) 
    call MPI_ISend(ySendBuff2(:,n),y2SendCount,mpirkind,neighbour(2),tagY2,&
      DECOMP_2D_COMM_CART_X,sendReq(n+12*1),ierr)
    call MPI_ISend(zSendBuff1(:,n),z1SendCount,mpirkind,neighbour(3),tagZ1,&
      DECOMP_2D_COMM_CART_X,sendReq(n+12*2),ierr)
    call MPI_ISend(zSendBuff2(:,n),z2SendCount,mpirkind,neighbour(4),tagZ2,&
      DECOMP_2D_COMM_CART_X,sendReq(n+12*3),ierr)
  end do
  call MPI_WaitAll(96,[recvReq, sendReq],MPI_STATUSES_IGNORE,ierr)
end subroutine

subroutine copyMode(this,cpyArr,modeID,periodicCorrection)
  class(enrichmentOperator), intent(inout) :: this
  real(rkind), dimension(12), intent(inout) :: cpyArr
  real(rkind), dimension(3), intent(in), optional :: periodicCorrection
  real(rkind), dimension(3) :: PC
  integer, intent(in) :: modeID

  PC = 0.d0
  if (present(periodicCorrection)) PC = periodicCorrection

  cpyArr(1)  = this%x(modeID) + PC(1)
  cpyArr(2)  = this%y(modeID) + PC(2)
  cpyArr(3)  = this%z(modeID) + PC(3)
  
  cpyArr(4)  = this%kx(modeID)
  cpyArr(5)  = this%ky(modeID)
  cpyArr(6)  = this%kz(modeID)
  
  cpyArr(7)  = this%uhatR(modeID)
  cpyArr(8)  = this%uhatI(modeID)
  cpyArr(9)  = this%vhatR(modeID)
  cpyArr(10) = this%vhatI(modeID)
  cpyArr(11) = this%whatR(modeID)
  cpyArr(12) = this%whatI(modeID)
end subroutine
  
subroutine getNeighbours(neighbour)
  integer, dimension(4), intent(out) :: neighbour
  integer :: ierr

  ! For X-pencil
  call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 0, 1, &
       neighbour(1), neighbour(2), ierr) ! north & south
  call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 1, 1, &
       neighbour(3), neighbour(4), ierr) ! top & bottom

end subroutine

subroutine getMPIcartCoords(coords)
  integer, dimension(2), intent(out) :: coords
  integer :: ierr

  call MPI_Cart_Coords(DECOMP_2D_COMM_CART_X, nrank, 2, coords, ierr)
end subroutine

subroutine testMPIsendRecv()
  integer :: tagY1, tagY2, tagZ1, tagZ2
  integer :: rankY1, rankY2, rankZ1, rankZ2
  integer, dimension(4) :: recvReq, sendReq
  integer :: n, ierr
  
  do n = 1,nproc
    if (nrank == n-1) then
      print*, "Rank", nrank, ". coords:", coords
    end if
    call MPI_Barrier(DECOMP_2D_COMM_CART_X,ierr)
  end do
  
  tagY1 = coords(1)
  if (coords(1) == dims(1)-1 .and. periodicBCs(2)) then
    tagY2 = 0
  else
    tagY2 = coords(1) + 1
  end if
  
  tagZ1 = coords(2)
  if (coords(2) == dims(2)-1 .and. periodicBCs(3)) then
    tagZ2 = 0
  else
    tagZ2 = coords(2) + 1
  end if
  
  ! Receive from lower y-rank
  call MPI_IRecv(rankY1,1,MPI_INTEGER,neighbour(1),tagY1,&
    DECOMP_2D_COMM_CART_X,recvReq(1),ierr)
  ! Receive from higher y-rank
  call MPI_IRecv(rankY2,1,MPI_INTEGER,neighbour(2),tagY2,&
    DECOMP_2D_COMM_CART_X,recvReq(2),ierr)
  ! Receive from lower z-rank
  call MPI_IRecv(rankZ1,1,MPI_INTEGER,neighbour(3),tagZ1,&
    DECOMP_2D_COMM_CART_X,recvReq(3),ierr)
  ! Receive from higher z-rank
  call MPI_IRecv(rankZ2,1,MPI_INTEGER,neighbour(4),tagZ2,&
    DECOMP_2D_COMM_CART_X,recvReq(4),ierr)
  
  ! Send to lower y-rank
  call MPI_ISend(nrank,1,MPI_INTEGER,neighbour(1),tagY1,&
    DECOMP_2D_COMM_CART_X,sendReq(1),ierr)
  ! Send to higher y-rank
  call MPI_ISend(nrank,1,MPI_INTEGER,neighbour(2),tagY2,&
    DECOMP_2D_COMM_CART_X,sendReq(2),ierr)
  ! Send to lower z-rank
  call MPI_ISend(nrank,1,MPI_INTEGER,neighbour(3),tagZ1,&
    DECOMP_2D_COMM_CART_X,sendReq(3),ierr)
  ! Send to higher z-rank
  call MPI_ISend(nrank,1,MPI_INTEGER,neighbour(4),tagZ2,&
    DECOMP_2D_COMM_CART_X,sendReq(4),ierr)

  call MPI_WaitAll(8,[recvReq, sendReq],MPI_STATUSES_IGNORE,ierr)

  do n = 1,nproc
    if (nrank == n-1) then
      print*, "Rank", nrank, ". Y neighbours:", rankY1, rankY2
      print*, "Rank", nrank, ". Z neighbours:", rankZ1, rankZ2
    end if
    call MPI_Barrier(DECOMP_2D_COMM_CART_X,ierr)
  end do
  call assert(rankY1 == neighbour(1),'rankY1 == neighbour(1)')
  call assert(rankY2 == neighbour(2),'rankY2 == neighbour(2)')
  call assert(rankZ1 == neighbour(3),'rankZ1 == neighbour(3)')
  call assert(rankZ2 == neighbour(4),'rankZ2 == neighbour(4)')

  call MPI_Barrier(DECOMP_2D_COMM_CART_X,ierr)
  call message('Test PASSED!')
end subroutine

