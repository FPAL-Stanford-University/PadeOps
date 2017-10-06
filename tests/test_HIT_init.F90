program test_HIT_init
   use kind_parameters, only: rkind, clen
   use mpi
   use basic_io, only: read_2d_ascii 
   use decomp_2d
   use decomp_2d_io
   
   implicit none
   
   type(decomp_info) :: gp, gpE
   integer, parameter :: nx = 128, ny = 128, nz = 128
   real(rkind), dimension(:,:,:), allocatable :: u, v, wE
   real(rkind), dimension(:,:), allocatable :: tmp1, tmp2
   
   integer :: ierr 

   call mpi_init(ierr)
   call decomp_2d_init(nx,ny,nz,1,1)
   call get_decomp_info(gp)
   call decomp_info_init(nx,ny,nz+1, gpE)

   allocate(u(nx,ny,nz), v(nx,ny,nz), wE(nx,ny,nz+1))

   call read_2d_ascii(tmp1,"/home/aditya90/Codes/PadeOps/data/PadeOps_HIT_input_postProj.txt") 
   call read_2d_ascii(tmp2,"/home/aditya90/Codes/PadeOps/data/PadeOps_HIT_input_postProj_edge.txt") 

   u  = reshape(tmp1(:,1),[nx,ny,nz  ])
   v  = reshape(tmp1(:,2),[nx,ny,nz  ])
   wE = reshape(tmp2(:,1),[nx,ny,nz+1])


   deallocate(tmp1, tmp2)
   print*, u(43, 23, 32)
   print*, v(43, 23, 32)
   print*, wE(43, 23, 32)

   call decomp_2d_write_one(1,u ,'/home/aditya90/Codes/PadeOps/data/HIT_u_PostProj.dat',gp)
   call decomp_2d_write_one(1,v ,'/home/aditya90/Codes/PadeOps/data/HIT_v_PostProj.dat',gp)
   call decomp_2d_write_one(1,wE,'/home/aditya90/Codes/PadeOps/data/HIT_w_PostProj.dat',gpE)
  
   print*, "Done writing data"
   deallocate(u, v, wE)


   call MPI_Finalize(ierr)

end program 
