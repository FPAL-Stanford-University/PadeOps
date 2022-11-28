subroutine sendRecvHaloModes(this,modeData,coordinate,haloBuff,last,extraModeData)
  ! Send and recieve Gabor modes in the halo region of the MPI ranks
  ! Inputs:
  !     modeData --> all mode attributes, i.e. physical location, wave-vecotr
  !                  components, and vector component amplitudes
  !     coordinate --> "y" or "z" specifying which component of the halo is
  !                    being transferred
  !     extraModeData (optional) --> same info as modeData, but for addtional
  !                                  modes not native to the MPI rank
  ! In/Out:
  !     haloBuff --> buffer array to store the halo mode info
  !
  ! Outputs:
  !     last --> The last index of haloBuff with relevant mode data. This is 
  !              because the input haloBuff may be larger than required for the
  !              output data and so any data beyond index "last" is erroneous
  !              and should not be used when rendering the velocity field
  class(enrichmentOperator), intent(inout), target :: this
  real(rkind), dimension(:,:), intent(in), target :: modeData
  character(len=1), intent(in) :: coordinate
  real(rkind), dimension(:,:), allocatable, target, intent(inout) :: haloBuff
  integer, intent(out) :: last
  real(rkind), dimension(:,:), intent(in), target, optional :: extraModeData

  integer :: coordID
  real(rkind), dimension(:), pointer :: loc
  integer :: n
  real(rkind) :: haloMin, haloMax
  real(rkind), dimension(:,:), pointer :: sendBuff1,  sendBuff2, &
    recvBuff1, recvBuff2
  integer :: sendCount1, sendCount2, recvCount1, recvCount2
  integer :: idx1, idx2
  integer :: st, en, ierr
  integer :: tag1, tag2
  integer :: sendID, recvID
  integer, dimension(48) :: sendReq, recvReq
  real(rkind) :: Ly, Lz
  real(rkind) :: PEmin, PEmax
  real(rkind), dimension(3) :: periodicCorrection
  real(rkind), dimension(2) :: Dom

  ! The following are a subset of the MPI cartesion communicator attributes
  real(rkind) :: coord, dim
  logical :: periodic
  integer, dimension(2) :: neigh

  call assert(coordinate == 'y' .or. coordinate == 'z','Must select'//&
    ' coordinate to be "y" or "z"')

  sendBuff1 => null(); sendBuff2 => null()
  recvBuff1 => null(); recvBuff2 => null()

  Ly = yDom(2) - yDom(1)
  Lz = zDom(2) - zDom(1)
    
  ! Periodic correction if using periodic boundaries
  periodicCorrection = 0.d0

  select case (coordinate)
  case ('y')
    coordID = 2

    PEmin = this%QHgrid%yE(1)
    PEmax = this%QHgrid%yE(this%QHgrid%gpC%xsz(2)+1)
  
    ! Define the cutoff location for halo modes
    haloMin = this%Qhgrid%yE(1) + &
      & real(this%nysupp/2,rkind)*this%smallScales%dy
    haloMax = this%Qhgrid%yE(this%QHgrid%gpC%xsz(2)+1) - &
      & real(this%nysupp/2,rkind)*this%smallScales%dy
    
    if (periodicBCs(2)) periodicCorrection  = [0.d0, Ly, 0.d0]
    
    Dom      = yDom
    coord    = coords(1)
    dim      = dims(1)
    periodic = periodicBCs(2)
    neigh    = [neighbour(1), neighbour(2)]
    loc      => modeData(:,2)
  case ('z')
    coordID = 3

    PEmin = this%QHgrid%zE(1)
    PEmax = this%QHgrid%zE(this%QHgrid%gpC%xsz(3)+1)
  
    ! Define the cutoff location for halo modes
    haloMin = this%Qhgrid%zE(1) + &
      & real(this%nzsupp/2,rkind)*this%smallScales%dz
    haloMax = this%Qhgrid%zE(this%QHgrid%gpC%xsz(3)+1) - &
      & real(this%nzsupp/2,rkind)*this%smallScales%dz

    if (periodicBCs(3)) periodicCorrection  = [0.d0, 0.d0, Lz]
    
    Dom      = zDom
    coord    = coords(2)
    dim      = dims(2)
    periodic = periodicBCs(3)
    neigh    = [neighbour(3), neighbour(4)]
    loc      => modeData(:,3)
  end select

  ! Find all modes inside halo regions
  sendCount1 = 0
  sendCount2 = 0
  call getSendCount(loc,haloMin,haloMax,sendCount1,sendCount2)
  if (present(extraModeData)) then
    loc => extraModeData(:,coordID)
    call getSendCount(loc,haloMin,haloMax,sendCount1,sendCount2)
    loc => modeData(:,coordID)
  end if

  if (allocated(this%sendBuff)) then
    if (size(this%sendBuff,1) < sendCount1 + sendCount2) then
      deallocate(this%sendBuff)
      allocate(this%sendBuff(sendCount1 + sendCount2,12))
    end if
  else
    allocate(this%sendBuff(sendCount1 + sendCount2,12))
  end if

  en = 0
  if (sendCount1 > 0) then
    st = en + 1
    en = en + sendCount1
    sendBuff1 => this%sendBuff(st:en,:)
  end if
  
  if (sendCount2 > 0) then
    st = en + 1
    en = en + sendCount2
    sendBuff2 => this%sendBuff(st:en,:)
  end if

  ! Pack the send buffers
  ! We must check for two different cases:
  ! --> Case 1: is the mode within the halo region?
  ! --> Case 2: If so, does the halo-region correspond to a periodic
  !             boundary?
  ! If Case 2 is true then we must shift the mode location by one
  ! domain length for the velocity rendering routine to work
  idx1 = 1
  idx2 = 1
  call packSendBuffer(loc,haloMin,haloMax,PEmin,PEmax,Dom,modeData,&
    periodicCorrection,sendBuff1,sendBuff2,idx1,idx2)
  if (present(extraModeData)) then
    loc => extraModeData(:,coordID)
    call packSendBuffer(loc,haloMin,haloMax,PEmin,PEmax,Dom,extraModeData,&
      periodicCorrection,sendBuff1,sendBuff2,idx1,idx2)
    loc => modeData(:,coordID)
  end if
  
  ! Let your neighbors know how much data you will be sending
  tag1 = coord
  if (coord == dim - 1 .and. periodic) then
    tag2 = 0
  else
    tag2 = coord + 1
  end if
  
  ! Receive from lower rank
  call MPI_IRecv(recvCount1,1,MPI_INTEGER,neigh(1),tag1,&
    DECOMP_2D_COMM_CART_X,recvReq(1),ierr)
  ! Receive from higher rank
  call MPI_IRecv(recvCount2,1,MPI_INTEGER,neigh(2),tag2,&
    DECOMP_2D_COMM_CART_X,recvReq(2),ierr)
  
  ! Send to lower rank
  call MPI_ISend(sendCount1,1,MPI_INTEGER,neigh(1),tag1,&
    DECOMP_2D_COMM_CART_X,sendReq(1),ierr)
  ! Send to higher rank
  call MPI_ISend(sendCount2,1,MPI_INTEGER,neigh(2),tag2,&
    DECOMP_2D_COMM_CART_X,sendReq(2),ierr)

  call MPI_WaitAll(4,[recvReq(1:2), sendReq(1:2)],MPI_STATUSES_IGNORE,ierr)

  ! Allocate memory for the recieve buffers
  if (allocated(haloBuff)) then
    if (size(haloBuff,1) < recvCount1 + recvCount2) then
      deallocate(haloBuff)
    end if
  end if
  if (.not. allocated(haloBuff)) then
    if (recvCount1 + recvCount2 > 0) then
      allocate(haloBuff(recvCount1 + recvCount2,12))
    end if
  end if
  if (allocated(haloBuff)) halobuff = 0.d0

  en = 0
  if (recvCount1 > 0) then
    st = en + 1
    en = en + recvCount1
    recvBuff1 => haloBuff(st:en,:)
  end if

  if (recvCount2 > 0) then
    st = en + 1
    en = en + recvCount2
    recvBuff2 => haloBuff(st:en,:)
  end if
  last = en

  call message(3,'Sending/Receiving modes')
  sendID = 1
  recvID = 1
  do n = 1,12
    if (sendCount1 > 0) then
      call MPI_ISend(sendBuff1(:,n),sendCount1,mpirkind,neigh(1),tag1,&
        DECOMP_2D_COMM_CART_X,sendReq(sendID),ierr) 
      sendID = sendID + 1
    end if
    if (recvCount2 > 0) then
      call MPI_IRecv(recvBuff2(:,n),recvCount2,mpirkind,neigh(2),tag2,&
        DECOMP_2D_COMM_CART_X,recvReq(recvID),ierr)
      recvID = recvID + 1
    end if

    if (sendCount2 > 0) then
      call MPI_ISend(sendBuff2(:,n),sendCount2,mpirkind,neigh(2),tag2,&
        DECOMP_2D_COMM_CART_X,sendReq(sendID),ierr)
      sendID = sendID + 1
    end if
    if (recvCount1 > 0) then
      call MPI_IRecv(recvBuff1(:,n),recvCount1,mpirkind,neigh(1),tag1,&
        DECOMP_2D_COMM_CART_X,recvReq(recvID),ierr)
      recvID = recvID + 1
    end if
  end do
  call MPI_WaitAll(recvID + sendID - 2,[recvReq(1:recvID-1), sendReq(1:sendID-1)], &
    MPI_STATUSES_IGNORE,ierr)
  
  call message(3,'All modes sent/received!')
  
  nullify(sendBuff1, sendBuff2)
  nullify(recvBuff1, recvBuff2)
  nullify(loc)
end subroutine

subroutine copyMode(modeData, cpyArr, periodicCorrection)
  real(rkind), dimension(:), intent(in) :: modeData
  real(rkind), dimension(:), intent(inout) :: cpyArr
  real(rkind), dimension(3), intent(in), optional :: periodicCorrection
  real(rkind), dimension(3) :: PC

  PC = 0.d0
  if (present(periodicCorrection)) PC = periodicCorrection

  cpyArr(1)  = modeData(1)  + PC(1)
  cpyArr(2)  = modeData(2)  + PC(2)
  cpyArr(3)  = modeData(3)  + PC(3)
  
  cpyArr(4)  = modeData(4) 
  cpyArr(5)  = modeData(5) 
  cpyArr(6)  = modeData(6) 

  cpyArr(7)  = modeData(7) 
  cpyArr(8)  = modeData(8) 
  cpyArr(9)  = modeData(9) 
  cpyArr(10) = modeData(10) 
  cpyArr(11) = modeData(11) 
  cpyArr(12) = modeData(12)

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

  call MPI_WaitAll(8,[recvReq(1:4), sendReq(1:4)],MPI_STATUSES_IGNORE,ierr)

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

subroutine getSendCount(loc,haloMin,haloMax,sendCount1,sendCount2)
  real(rkind), dimension(:), intent(in) :: loc
  real(rkind), intent(in) :: haloMin, haloMax
  integer, intent(inout) :: sendCount1, sendCount2
  integer :: n

  do n = 1,size(loc)
    if (loc(n) < haloMin) then
      sendCount1 = sendCount1 + 1
    else if (loc(n) > haloMax) then
      sendCount2 = sendCount2 + 1
    end if
  end do
end subroutine

subroutine packSendBuffer(loc,haloMin,haloMax,PEmin,PEmax,Dom,modeData,&
    periodicCorrection,sendBuff1,sendBuff2,idx1,idx2)
  real(rkind), dimension(:), intent(in) :: loc
  real(rkind), intent(in) :: haloMin, haloMax, PEmin, PEmax
  real(rkind), dimension(2), intent(in) :: Dom
  real(rkind), dimension(:,:), intent(in) :: modeData
  real(rkind), dimension(3), intent(in) :: periodicCorrection
  real(rkind), dimension(:,:), intent(inout) :: sendBuff1, sendBuff2
  integer, intent(inout) :: idx1, idx2
  real(rkind) :: tol = 1.d-12
  integer :: n

  do n = 1,size(loc)
    if (loc(n) < haloMin .and. abs(PEmin - Dom(1)) < tol .and. &
      periodicBCs(2)) then
      call copyMode(modeData(n,:), sendBuff1(idx1,:), periodicCorrection)
      idx1 = idx1 + 1
    else if (loc(n) < haloMin .and. abs(PEmin - Dom(1)) >= tol) then
      call copyMode(modeData(n,:), sendBuff1(idx1,:))
      idx1 = idx1 + 1
    else if (loc(n) > haloMax .and. abs(PEmax - Dom(2)) < tol .and. &
      periodicBCs(2)) then
      call copyMode(modeData(n,:), sendBuff2(idx2,:), -periodicCorrection)
      idx2 = idx2 + 1
    else if (loc(n) > haloMax .and. abs(PEmax - Dom(2)) >= tol) then
      call copyMode(modeData(n,:), sendBuff2(idx2,:))
      idx2 = idx2 + 1
    end if
  end do
end subroutine
