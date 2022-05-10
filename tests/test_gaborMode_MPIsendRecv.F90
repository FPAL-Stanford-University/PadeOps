program test_gaborMode_MPIsendRecv
  use mpi
  use hdf5
  use decomp_2D
  use kind_parameters, only: rkind
  use GaborModeRoutines, only: mpiIsendIrecv
  implicit none

  real(rkind), dimension(3,3,3) :: sendbuff, recvbuff
  integer :: sendReq, recvReq, tag
  integer :: sendNeigh, recvNeigh
  integer :: ierr

  call MPI_Init(ierr)
  call decomp_2d_init(100,100,100,0,0)
  sendbuff = real(nrank,rkind)
  if (nrank == 0) then
    sendNeigh = nrank + 1
    recvNeigh = nproc - 1
  elseif (nrank == nproc-1) then
    sendNeigh = 0
    recvNeigh = nrank - 1
  else
    sendNeigh = nrank + 1
    recvNeigh = nrank - 1
  end if 

  print*, nrank, sendNeigh, recvNeigh
  tag = 1
  call mpiIsendIrecv(sendBuff, recvBuff, sendNeigh, recvNeigh, sendReq,&
    recvReq, tag)
  call MPI_Waitall(2,[sendReq, recvReq],MPI_STATUSES_IGNORE,ierr)
  print*, nrank, maxval(recvBuff), minval(recvBuff) 
 
  call decomp_2d_finalize() 
  call MPI_Finalize(ierr)
end program
