program test_randomnumbers
   use kind_parameters, only: rkind, clen
   use decomp_2d 
   use mpi
   use random
   use exits, only: message
   
   implicit none 

   real(rkind), dimension(:,:,:), allocatable :: randArr 
   type(decomp_info) :: gpC
   integer :: nx=64, ny=64, nz=64
   integer :: seed = 134438, ierr 
   call MPI_Init(ierr)
   call decomp_2d_init(nx, ny, nz, 0, 0)
   call get_decomp_info(gpC)
      

   ! Same size array for all procs (testing)
   allocate(randArr(nx,ny,nz))

   ! Same random numers for each proc

   call message(0,"Gaussian random variables(same across procs)")
   call gaussian_random(randArr, 0._rkind, 1._rkind, seed)
   call sleep(nrank)
   print*, "Rank:", nrank
   print*, "First 6:", randArr(1:6,1,1)
   print*, "Mean:", sum(randArr)/real(nx*ny*nz,rkind)
   print*, "Variance:", sum((randArr - sum(randArr)/real(nx*ny*nz,rkind))**2)/real(nx*ny*nz,rkind)
   print*, "---------------"
  
   call mpi_barrier(mpi_comm_world, ierr)
   call message("=====================================")
   call message("=====================================")
   
   
   call message(0,"Gaussian random variables(varying across procs)")
   call gaussian_random(randArr, 0._rkind, 1._rkind, seed*nrank)
   call sleep(nrank)
   print*, "Rank:", nrank
   print*, "First 6:", randArr(1:6,1,1)
   print*, "Mean:", sum(randArr)/real(nx*ny*nz,rkind)
   print*, "Variance:", sum((randArr - sum(randArr)/real(nx*ny*nz,rkind))**2)/real(nx*ny*nz,rkind)
   print*, "---------------"


   call MPI_Finalize(ierr)

end program 
