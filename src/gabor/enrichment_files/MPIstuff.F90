subroutine sendRecvHaloModes(this,modeData)
  ! Send and recieve Gabor modes in the halo region of the MPI ranks
  ! In/Out:
  !     modeData --> all mode attributes, i.e. physical location, wave-vector
  !                  components, and vector component amplitudes
  
  class(enrichmentOperator), intent(inout) :: this
  real(rkind), dimension(:,:), intent(inout), allocatable, target :: modeData

  real(rkind), dimension(:), pointer :: yloc, zloc
  real(rkind) :: haloMinY, haloMaxY, haloMinZ, haloMaxZ
  integer :: sendCount
  integer, dimension(nproc) :: recvCount, displs
  integer :: ierr
  real(rkind) :: Ly, Lz
  real(rkind) :: PEminY, PEmaxY, PEminZ, PEmaxZ
  real(rkind), dimension(3) :: periodicShiftY, periodicShiftZ
  real(rkind), dimension(:,:), allocatable :: sendBuff, recvBuff, newModes
  integer :: numNewModes, st, stg, en, eng, n

  !PEminY = this%sd%PEybound(1)
  !PEmaxY = this%sd%PEybound(2)
  !PEminZ = this%sd%PEzbound(1)
  !PEmaxZ = this%sd%PEzbound(2)
  PEminY = this%PEybound(1)
  PEmaxY = this%PEybound(2)
  PEminZ = this%PEzbound(1)
  PEmaxZ = this%PEzbound(2)
    
  ! Define the cutoff location for halo modes
  haloMinY = PEminY + real(this%nysupp/2 + this%haloPad,rkind)*this%smallScales%dy
  haloMaxY = PEmaxY - real(this%nysupp/2 + this%haloPad,rkind)*this%smallScales%dy
  haloMinZ = PEminZ + real(this%nzsupp/2 + this%haloPad,rkind)*this%smallScales%dz
  haloMaxZ = PEmaxZ - real(this%nzsupp/2 + this%haloPad,rkind)*this%smallScales%dz
  
  yloc => modeData(:,2)
  zloc => modeData(:,3)
    
  ! Periodic correction if using periodic boundaries
  periodicShiftY = 0.d0
  periodicShiftZ = 0.d0
  Ly = yDom(2) - yDom(1)
  Lz = zDom(2) - zDom(1)
  if (periodicBCs(2)) periodicShiftY  = [0.d0, Ly, 0.d0]
  if (periodicBCs(3)) periodicShiftZ  = [0.d0, 0.d0, Lz]

  ! Find all modes inside halo regions
  sendCount = 0
  call getSendCount(yloc,zloc,haloMinY,haloMaxY,haloMinZ,haloMaxZ,&
    PEminY, PEmaxY, PEminZ, PEmaxZ, sendCount)

  allocate(sendBuff(sendCount,this%nvars))
  sendBuff = 0.d0

  call packSendBuffer(yloc,zloc,haloMinY,haloMaxY,haloMinZ,haloMaxZ,&
    PEminY,PEmaxY,PEminZ,PEmaxZ,yDom,zDom,modeData,periodicShiftY,&
    periodicShiftZ,sendBuff)
 
  ! Get the number of modes to be sent by each rank 
  recvCount = 0
  call MPI_Allgather(sendCount,1,MPI_INTEGER,recvCount,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

  ! Allocate memory for the recieve buffer
  allocate(recvBuff(sum(recvCount),this%nvars))
  recvBuff = 0.d0

  ! Compute the displacements for the receive buffer
  displs = 0
  do n = 2,nproc
    displs(n) = displs(n-1) + recvCount(n-1)
  end do

  ! Get all the modes from all other ranks
  do n = 1,this%nvars
    call MPI_AllgatherV(sendBuff(:,n),size(sendBuff,1),mpirkind,recvBuff(:,n),recvCount,&
      displs,mpirkind,MPI_COMM_WORLD,ierr)
  end do

  numNewModes = sum(recvCount(1:nrank)) + sum(recvCount(nrank+2:nproc))
  allocate(newModes(numNewModes,this%nvars))
  st = 1
  stg = 1
  do n = 1,nproc
    eng = stg + recvCount(n) - 1
    if (nrank == n-1) then
      continue
    else
      en = st + recvCount(n) - 1
      newModes(st:en,:) = recvBuff(stg:eng,:)
      st = en + 1
    end if
    stg = eng + 1
  end do

  nullify(yloc, zloc)
  ! Append modeData
  call mergeModes(modeData, newModes, this%nvars)

  deallocate(recvBuff,sendBuff,newModes)
end subroutine

subroutine mergeModes(mergedModes, newModes, nvars)
  real(rkind), dimension(:,:), intent(in) :: newModes
  integer, intent(in) :: nvars
  real(rkind), dimension(:,:), allocatable, intent(inout) :: mergedModes
  integer :: n1, n2
  real(rkind), dimension(:,:), allocatable :: tmpData
  
  allocate(tmpData(size(mergedModes,1),nvars))
  tmpData = mergedModes
  deallocate(mergedModes)

  n1 = size(tmpData,1)
  n2 = size(newModes,1)

  allocate(MergedModes(n1 + n2,nvars))
  MergedModes(1:n1,      :) = tmpData
  MergedModes(n1+1:n1+n2,:) = newModes

  deallocate(tmpData)

end subroutine 
  
subroutine getneighbors(neighbor,PEybound,PEzbound,periodic)
  integer, dimension(4), intent(out) :: neighbor
  real(rkind), dimension(2), intent(in) :: PEybound, PEzbound
  logical, dimension(3), intent(in) :: periodic
  integer :: ierr, i

  ! For X-pencil
  call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 0, 1, &
       neighbor(1), neighbor(2), ierr) ! north & south
  call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 1, 1, &
       neighbor(3), neighbor(4), ierr) ! top & bottom

  do i = 1,2
    if ((abs(PEybound(i) - yDom(i)) < tol) .and.(.not. periodic(2)) ) then
      neighbor(i) = MPI_PROC_NULL
    end if
    if ((abs(PEzbound(i) - zDom(i)) < tol) .and.(.not. periodic(3)) ) then
      neighbor(i+2) = MPI_PROC_NULL
    end if
  end do

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
  integer, dimension(8) :: requests 
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

  requests(1:4) = recvReq
  requests(5:8) = sendReq
  call MPI_WaitAll(8,requests,MPI_STATUSES_IGNORE,ierr)

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

subroutine getSendCount(yloc,zloc,haloMinY,haloMaxY,haloMinZ,haloMaxZ,&
    PEminY, PEmaxY, PEminZ, PEmaxZ, sendCount)
  real(rkind), dimension(:), intent(in) :: yloc, zloc
  real(rkind), intent(in) :: haloMinY, haloMaxY, haloMinZ, haloMaxZ, &
    PEminY, PEmaxY, PEminZ, PEmaxZ
  integer, intent(inout) :: sendCount
  integer :: n, sendCountAdd = 0
  
  do n = 1,size(yloc)
    if (insideHaloRegion(yloc(n),zloc(n),haloMinY,haloMaxY,haloMinZ,haloMaxZ)) then 
    !if (yloc(n) < haloMinY .or. yloc(n) > haloMaxY .or. &
    !    zloc(n) < haloMinZ .or. zloc(n) > haloMaxZ) then
      sendCount = sendCount + 1
    end if
  end do

  ! Need to account for periodic shifts as well. For example, if the problem is 
  ! periodic in Y, each mode must be copied twice at y = y+Ly and y = y-Ly.
  ! If it is periodic in Z there is another factor of two, and if both are 
  ! periodic then another factor of 2 to account for (y,z) --> (y+/-Ly,z+/-Lz)
  ! This is inefficient in terms of communication and storage, but simple
  ! to implement. A future task can be to optimize this
  if (periodicBCs(2))                      sendCountAdd = sendCountAdd + 2*sendCount
  if (periodicBCs(3))                      sendCountAdd = sendCountAdd + 2*sendCount
  if (periodicBCs(2) .and. periodicBCs(3)) sendCountAdd = sendCountAdd + 4*sendCount
  sendCount = sendCount + sendCountAdd
end subroutine

subroutine packSendBuffer(yloc,zloc,haloMinY,haloMaxY,haloMinZ,haloMaxZ,&
    PEminY,PEmaxY,PEminZ,PEmaxZ,yDom,zDom,modeData,&
    periodicShiftY,periodicShiftZ,sendBuff)
  ! Pack the send buffers
  ! We must check for two different cases:
  ! --> Case 1: is the mode within the halo region?
  ! --> Case 2: If so, does the halo-region correspond to a periodic
  !             boundary?
  ! If Case 2 is true then we must shift the mode location by one
  ! domain length for the velocity rendering routine to work
  use GaborModeRoutines, only: copyMode
  real(rkind), dimension(:), intent(in) :: yloc, zloc
  real(rkind), intent(in) :: haloMinY, haloMaxY, PEminY, PEmaxY
  real(rkind), intent(in) :: haloMinZ, haloMaxZ, PEminZ, PEmaxZ
  real(rkind), dimension(2), intent(in) :: yDom, zDom
  real(rkind), dimension(:,:), intent(in) :: modeData
  real(rkind), dimension(3), intent(in) :: periodicShiftY, periodicShiftZ
  real(rkind), dimension(:,:), intent(inout) :: sendBuff
  integer :: idx
  integer :: n

  idx = 0
  do n = 1,size(yloc)
    if (insideHaloRegion(yloc(n),zloc(n),haloMinY,haloMaxY,haloMinZ,haloMaxZ)) then 
    !if (yloc(n) < haloMinY .or. yloc(n) > haloMaxY .or. &
    !    zloc(n) < haloMinZ .or. zloc(n) > haloMaxZ) then
      idx = idx + 1
      call copyMode(modeData(n,:), sendBuff(idx,:), [0.d0, 0.d0, 0.d0])
    end if
  end do
  if (periodicBCs(2)) then
    call copyMode(sendBuff(1:idx,:), sendBuff(  idx+1:2*idx,:),  periodicShiftY)
    call copyMode(sendbuff(1:idx,:), sendBuff(2*idx+1:3*idx,:), -periodicShiftY)
    idx = 3*idx
  end if
  if (periodicBCs(3)) then
    call copyMode(sendBuff(1:idx,:), sendBuff(  idx+1:2*idx,:),  periodicShiftZ)
    call copyMode(sendbuff(1:idx,:), sendBuff(2*idx+1:3*idx,:), -periodicShiftZ)
    idx = 3*idx
  end if
end subroutine

subroutine sortAndExchangeModes(this, coor, coorMin, coorMax, neighLo, neighHi)
  class(enrichmentOperator), intent(inout), target :: this
  integer :: n, inactiveIndex, iterSelf, iterLo, iterHi
  integer :: howmanyLo, howmanyHi
  real(rkind), dimension(:), intent(in) :: coor
  real(rkind), intent(in) :: coorMin, coorMax
  integer, intent(in) :: neighLo, neighHi
  real(rkind), dimension(:,:), allocatable :: sendBufferLo, sendBufferHi
  real(rkind), dimension(:,:), allocatable :: recvBufferLo, recvBufferHi
  real(rkind), dimension(:,:,:), allocatable :: tmp
  integer, dimension(4) :: sendReq, recvReq, requests
  integer :: tag, ierr, sz1

  inactiveIndex = mod(this%activeIndex+1,2)
  
  allocate(sendBufferLo(this%nmodes,this%nvars))
  allocate(sendBufferHi(this%nmodes,this%nvars))

  iterSelf = 0
  iterLo = 0
  iterHi = 0

  do n = 1,this%nModes
      
      if ((coor(n) <= coorMax) .and. (coor(n) > coorMin)) then 
          iterSelf = iterSelf + 1
          this%rawModedata(iterSelf,:,inactiveIndex) = this%rawModedata(n,:,this%activeIndex) 
      else
          if (coor(n) > coorMax) then 
              iterHi = iterHi + 1
              sendBufferHi(iterHi,:) = this%rawModedata(n,:,this%activeIndex) 

          else if (coor(n) <= coorMin) then 
              iterLo = iterLo + 1
              sendBufferLo(iterLo,:) = this%rawModedata(n,:,this%activeIndex) 
          
          else
              print*, coor(n), coorMax, coorMin
              print*, "This should not happen"
              stop 
          end if 
      end if 

  end do


  ! isend howmany to hi
  ! isend howmany to lo
  tag = 0 
  ! Send to lower rank
  call MPI_ISend(iterLo,1,MPI_INTEGER,neighLo,tag,&
    DECOMP_2D_COMM_CART_X,sendReq(1),ierr)
  ! Send to higher rank
  call MPI_ISend(iterHi,1,MPI_INTEGER,neighHi,tag,&
    DECOMP_2D_COMM_CART_X,sendReq(2),ierr)
  
  ! irecv howmany from lo
  ! irecv howmany from hi
  call MPI_IRecv(howManyLo,1,MPI_INTEGER,neighLo,tag,&
    DECOMP_2D_COMM_CART_X,recvReq(1),ierr)
  ! Receive from higher rank
  call MPI_IRecv(howManyHi,1,MPI_INTEGER,neighHi,tag,&
    DECOMP_2D_COMM_CART_X,recvReq(2),ierr)
  
  ! isend data to lo
  ! isend data to hi 
  call MPI_ISend(sendBufferLo(1:iterLo,:),this%nvars*iterLo,mpirkind,neighLo,tag,&
      DECOMP_2D_COMM_CART_X,sendReq(3),ierr) 
  call MPI_ISend(sendBufferHi(1:iterHi,:),this%nvars*iterHi,mpirkind,neighHi,tag,&
      DECOMP_2D_COMM_CART_X,sendReq(4),ierr) 
  
  requests(1:2) = recvReq(1:2)
  requests(3:4) = sendReq(1:2)
  call MPI_WaitAll(4,requests,MPI_STATUSES_IGNORE,ierr)

  allocate(recvBufferLo(howmanyLo,this%nvars))
  allocate(recvBufferHi(howmanyHi,this%nvars))

  ! irecv data from lo
  ! irecv data from hi 
  call MPI_IRecv(recvBufferLo,size(recvBufferLo),mpirkind,neighLo,tag,&
    DECOMP_2D_COMM_CART_X,recvReq(3),ierr)
  call MPI_IRecv(recvBufferHi,size(recvBufferHi),mpirkind,neighHi,tag,&
    DECOMP_2D_COMM_CART_X,recvReq(4),ierr)

  requests(1:2) = recvReq(3:4)
  requests(3:4) = sendReq(3:4)
  call MPI_WaitAll(4,requests, MPI_STATUSES_IGNORE,ierr)

  if (size(this%rawModedata,1) < iterSelf+howmanyLo+howmanyHi) then
      sz1 = size(this%rawModeData,1)
      
      ! Temporarily copy rawModeData to tmp array
      allocate(tmp(sz1,size(this%rawModeData,2),size(this%rawModeData,3)))
      tmp = this%rawModeData

      ! Deallocate and reallocate rawModeData
      deallocate(this%rawModeData)
      allocate(this%rawModeData(iterSelf+howManyLo+howManyHi,this%nvars,0:1))

      ! Copy tmp data back to rawModeData
      this%rawModeData(1:sz1,:,:) = tmp

      ! Clear memory
      deallocate(tmp)
  end if 

  this%rawModedata(iterSelf+1:iterSelf+howmanyLo,:,inactiveIndex) = recvBufferLo 
  this%rawModedata(iterSelf+howmanyLo+1:iterSelf+howmanyLo+howmanyHi,:,inactiveIndex) = recvBufferHi 

  this%nModes = iterSelf+howmanyLo+howmanyHi

  call this%togglePointer()

  deallocate(sendBufferLo, sendBufferHi, recvBufferLo, recvBufferHi)

end subroutine 
function insideHaloRegion(yloc,zloc,haloMinY,haloMaxY,haloMinZ,haloMaxZ) result (modeInside)
  real(rkind), intent(in) :: yloc, zloc, haloMinY, haloMaxY, haloMinZ, haloMaxZ
  logical :: modeInside
  
  modeInside = .false.  
  if (yloc < haloMinY .or. yloc > haloMaxY .or. &
      zloc < haloMinZ .or. zloc > haloMaxZ) modeInside = .true.

end function
