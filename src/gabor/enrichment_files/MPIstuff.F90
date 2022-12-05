subroutine sendRecvHaloModes(this,modeData,coordinate,haloBuff)
  ! Send and recieve Gabor modes in the halo region of the MPI ranks
  ! Inputs:
  !     modeData --> all mode attributes, i.e. physical location, wave-vecotr
  !                  components, and vector component amplitudes
  !     coordinate --> "y" or "z" specifying which component of the halo is
  !                    being transferred
  !                                  modes not native to the MPI rank
  ! In/Out:
  !     haloBuff --> buffer array to store the halo mode info
  !
  class(enrichmentOperator), intent(inout) :: this
  real(rkind), dimension(:,:), intent(inout), allocatable, target :: modeData
  character(len=1), intent(in) :: coordinate
  real(rkind), dimension(:,:), allocatable, intent(inout), optional :: haloBuff

  integer :: coordID
  real(rkind), dimension(:), pointer :: loc
  real(rkind) :: haloMin, haloMax
  integer :: sendCount1, sendCount2, recvCount1, recvCount2
  integer :: ierr
  integer :: tag1, tag2
  integer, dimension(2) :: sendReq, recvReq
  real(rkind) :: Ly, Lz
  real(rkind) :: PEmin, PEmax
  real(rkind), dimension(3) :: periodicCorrection
  real(rkind), dimension(2) :: Dom
  real(rkind), dimension(:,:), allocatable :: sendBuff1, sendBuff2, recvBuff1, recvBuff2, tmpData

  ! The following are a subset of the MPI cartesion communicator attributes
  real(rkind) :: coord, dim
  logical :: periodic
  integer, dimension(2) :: neigh

  call assert(coordinate == 'y' .or. coordinate == 'z','Must select'//&
    ' coordinate to be "y" or "z"')

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
    neigh    = [neighbor(1), neighbor(2)]
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
    neigh    = [neighbor(3), neighbor(4)]
    loc      => modeData(:,3)
  end select

  ! Find all modes inside halo regions
  sendCount1 = 0
  sendCount2 = 0
  call getSendCount(loc,haloMin,haloMax,sendCount1,sendCount2)

  allocate(sendBuff1(sendCount1,12))
  allocate(sendBuff2(sendCount2,12))

  ! Pack the send buffers
  ! We must check for two different cases:
  ! --> Case 1: is the mode within the halo region?
  ! --> Case 2: If so, does the halo-region correspond to a periodic
  !             boundary?
  ! If Case 2 is true then we must shift the mode location by one
  ! domain length for the velocity rendering routine to work
  call packSendBuffer(loc,haloMin,haloMax,PEmin,PEmax,Dom,modeData,&
    periodicCorrection,sendBuff1,sendBuff2)
  
  ! Let your neighbors know how much data you will be sending
  tag1 = 0
  tag2 = 0
  
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

  call MPI_WaitAll(4,[recvReq, sendReq],MPI_STATUSES_IGNORE,ierr)

  ! Allocate memory for the recieve buffers
  allocate(recvBuff1(recvCount1,12))
  allocate(recvBuff2(recvCount2,12))
  recvBuff1 = 0.d0
  recvBuff2 = 0.d0

  call MPI_IRecv(recvBuff1,size(recvBuff1),mpirkind,neigh(1),tag1,&
    DECOMP_2D_COMM_CART_X,recvReq(1),ierr)
  call MPI_IRecv(recvBuff2,size(recvBuff2),mpirkind,neigh(2),tag2,&
    DECOMP_2D_COMM_CART_X,recvReq(2),ierr)
  call MPI_ISend(sendBuff1,size(sendBuff1),mpirkind,neigh(1),tag1,&
      DECOMP_2D_COMM_CART_X,sendReq(1),ierr) 
  call MPI_ISend(sendBuff2,size(sendBuff2),mpirkind,neigh(2),tag2,&
      DECOMP_2D_COMM_CART_X,sendReq(2),ierr) 

  call MPI_WaitAll(4,[recvReq, sendReq], MPI_STATUSES_IGNORE,ierr)

  if (present(haloBuff)) then 
    if (allocated(haloBuff)) deallocate(haloBuff) 
    call mergeModes(modeData, recvBuff1, recvBuff2, haloBuff) 
  else ! We will append modeData
    allocate(tmpData(size(modeData,1),12))
    tmpData = modeData
    deallocate(modeData)
    call mergeModes(tmpData, recvBuff1, recvBuff2, modeData)
  end if 

  deallocate(recvBuff1,recvBuff2,sendBuff1,sendBuff2)
  nullify(loc)
end subroutine

subroutine mergeModes(Modes1, Modes2, Modes3, MergedModes)
  real(rkind), dimension(:,:), intent(in) :: Modes1, Modes2, Modes3
  real(rkind), dimension(:,:), allocatable, intent(out) :: MergedModes
  integer :: n1, n2, n3 

  n1 = size(Modes1,1)
  n2 = size(Modes2,1)
  n3 = size(Modes3,1)

  allocate(MergedModes(n1 + n2 + n3,12))
  MergedModes(1:n1,            :) = Modes1
  MergedModes(n1+1:n1+n2,      :) = Modes2
  MergedModes(n1+n2+1:n1+n2+n3,:) = Modes3

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
  
subroutine getneighbors(neighbor)
  integer, dimension(4), intent(out) :: neighbor
  integer :: ierr

  ! For X-pencil
  call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 0, 1, &
       neighbor(1), neighbor(2), ierr) ! north & south
  call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 1, 1, &
       neighbor(3), neighbor(4), ierr) ! top & bottom

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
  call MPI_IRecv(rankY1,1,MPI_INTEGER,neighbor(1),tagY1,&
    DECOMP_2D_COMM_CART_X,recvReq(1),ierr)
  ! Receive from higher y-rank
  call MPI_IRecv(rankY2,1,MPI_INTEGER,neighbor(2),tagY2,&
    DECOMP_2D_COMM_CART_X,recvReq(2),ierr)
  ! Receive from lower z-rank
  call MPI_IRecv(rankZ1,1,MPI_INTEGER,neighbor(3),tagZ1,&
    DECOMP_2D_COMM_CART_X,recvReq(3),ierr)
  ! Receive from higher z-rank
  call MPI_IRecv(rankZ2,1,MPI_INTEGER,neighbor(4),tagZ2,&
    DECOMP_2D_COMM_CART_X,recvReq(4),ierr)
  
  ! Send to lower y-rank
  call MPI_ISend(nrank,1,MPI_INTEGER,neighbor(1),tagY1,&
    DECOMP_2D_COMM_CART_X,sendReq(1),ierr)
  ! Send to higher y-rank
  call MPI_ISend(nrank,1,MPI_INTEGER,neighbor(2),tagY2,&
    DECOMP_2D_COMM_CART_X,sendReq(2),ierr)
  ! Send to lower z-rank
  call MPI_ISend(nrank,1,MPI_INTEGER,neighbor(3),tagZ1,&
    DECOMP_2D_COMM_CART_X,sendReq(3),ierr)
  ! Send to higher z-rank
  call MPI_ISend(nrank,1,MPI_INTEGER,neighbor(4),tagZ2,&
    DECOMP_2D_COMM_CART_X,sendReq(4),ierr)

  call MPI_WaitAll(8,[recvReq(1:4), sendReq(1:4)],MPI_STATUSES_IGNORE,ierr)

  do n = 1,nproc
    if (nrank == n-1) then
      print*, "Rank", nrank, ". Y neighbors:", rankY1, rankY2
      print*, "Rank", nrank, ". Z neighbors:", rankZ1, rankZ2
    end if
    call MPI_Barrier(DECOMP_2D_COMM_CART_X,ierr)
  end do
  call assert(rankY1 == neighbor(1),'rankY1 == neighbor(1)')
  call assert(rankY2 == neighbor(2),'rankY2 == neighbor(2)')
  call assert(rankZ1 == neighbor(3),'rankZ1 == neighbor(3)')
  call assert(rankZ2 == neighbor(4),'rankZ2 == neighbor(4)')

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
    periodicCorrection,sendBuff1,sendBuff2)
  real(rkind), dimension(:), intent(in) :: loc
  real(rkind), intent(in) :: haloMin, haloMax, PEmin, PEmax
  real(rkind), dimension(2), intent(in) :: Dom
  real(rkind), dimension(:,:), intent(in) :: modeData
  real(rkind), dimension(3), intent(in) :: periodicCorrection
  real(rkind), dimension(:,:), intent(inout) :: sendBuff1, sendBuff2
  integer :: idx1, idx2
  real(rkind) :: tol = 1.d-12
  integer :: n

  idx1 = 1
  idx2 = 1
  do n = 1,size(loc)
    if (loc(n) < haloMin) then
      if (abs(PEmin - Dom(1)) < tol .and. periodicBCs(2)) then
        call copyMode(modeData(n,:), sendBuff1(idx1,:), periodicCorrection)
        idx1 = idx1 + 1
      else
        call copyMode(modeData(n,:), sendBuff1(idx1,:))
        idx1 = idx1 + 1
      end if
    else if (loc(n) > haloMax) then
      if (abs(PEmax - Dom(2)) < tol .and. periodicBCs(2)) then
        call copyMode(modeData(n,:), sendBuff2(idx2,:), -periodicCorrection)
        idx2 = idx2 + 1
      else
        call copyMode(modeData(n,:), sendBuff2(idx2,:))
        idx2 = idx2 + 1
      end if
    end if
  end do
end subroutine

