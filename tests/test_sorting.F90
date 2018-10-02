program test_sorting
   use kind_parameters, only: rkind
   use sorting_mod, only: sortgroup, Qsort    
   use mpi
   use decomp_2d
   use constants, only: zero, two, pi
   use gridtools, only: linspace
   use reductions, only: p_sum
   use decomp_2d_io
   use timer, only: tic, toc

   real(rkind), dimension(:,:,:), allocatable :: zCoords, f, fsort, buffz
   real(rkind), dimension(:), allocatable :: z, y, x
   real(rkind) :: dx, dy, dz, PE, PE_loc, PE_old = 1.E10

   integer, parameter :: nx = 256, ny = 256, nz = 256
   real(rkind) :: xFull(nx), yFull(ny), zFull(nz)
   integer :: ierr, np, i, j, k, idx
   type(decomp_info) :: gp
   type(sortgroup), dimension(:), allocatable :: sort_z, sort_y
   integer :: iiter, niter = 8

   call mpi_init(ierr)
   call MPI_Comm_size ( MPI_COMM_WORLD, np, ierr )

   call decomp_2d_init(nx,ny,nz,1,np)
   call get_decomp_info(gp)

   allocate(f(gp%ysz(1),gp%ysz(2),gp%ysz(3)))
   allocate(fsort(gp%ysz(1),gp%ysz(2),gp%ysz(3)))
   allocate(buffz(gp%zsz(1),gp%zsz(2),gp%zsz(3)))
   allocate(z(gp%ysz(3)))
   allocate(y(gp%ysz(2)))
   allocate(x(gp%ysz(1)))
   allocate(zCoords(gp%ysz(1),gp%ysz(2),gp%ysz(3)))
   

   dx = two*pi/real(nx,rkind); dy = two*pi/real(ny,rkind); dz = two*pi/real(nz,rkind)
   
   xFull = linspace(zero, two*pi-dx, nx)
   yFull = linspace(zero, two*pi-dy, ny)
   zFull = linspace(zero, two*pi-dz, nz)

   x = xFull(gp%yst(1):gp%yen(1))
   y = yFull(gp%yst(2):gp%yen(2))
   z = zFull(gp%yst(3):gp%yen(3))

   do k = 1,gp%ysz(3)
      do j = 1,gp%ysz(2)
         do i = 1,gp%ysz(1)
            f(i,j,k) = 2.d0*cos(x(i))*sin(y(j))*cos(4*z(k)) + 0.3d0*cos(5.d0*(z(k)+2.d0*y(j)+6.d0*x(i)))+0.2d0*sin(8.d0*(x(i)*y(j)*z(k))**2)+0.25d0*sin(8.d0*(x(i)**2)*((y(j)+z(k))**2))
            zCoords(i,j,k) = z(k)
         end do 
      end do
   end do 

   PE_loc = zero
   do k = 1,gp%ysz(3)
      PE_loc = PE_loc - sum(f(:,:,k)*z(k))
   end do 
   PE = p_sum(PE_loc)
   if (nrank == 0) print*, "PE:", PE

   allocate(sort_z(gp%zsz(1)*gp%zsz(2)*gp%zsz(3))) 
   allocate(sort_y(gp%ysz(1)*gp%ysz(2)*gp%ysz(3))) 

   fsort = f
   

   do iiter = 1,niter!while (PE < (PE_old + 1.d-10))
      PE_old = PE
      call tic()
      call transpose_y_to_z(fsort,buffz,gp)
      idx = 1
      do k = 1,gp%zsz(3)
         do j = 1,gp%zsz(2)
            do i = 1,gp%zsz(1)
               sort_z(idx)%value = buffz(i,j,k)
               idx = idx + 1
            end do 
         end do 
      end do   
      call transpose_y_to_z(zCoords,buffz,gp)
      idx = 1
      do k = 1,gp%zsz(3)
         do j = 1,gp%zsz(2)
            do i = 1,gp%zsz(1)
               sort_z(idx)%zpos = buffz(i,j,k)
               idx = idx + 1
            end do 
         end do 
      end do   
  
      call Qsort(sort_z, gp%zsz(1)*gp%zsz(2)*gp%zsz(3)) 

      idx = 1
      do k = 1,gp%zsz(3)
         do j = 1,gp%zsz(2)
            do i = 1,gp%zsz(1)
               buffz(i,j,k) = sort_z(idx)%value
               idx = idx + 1
            end do 
         end do 
      end do
      call transpose_z_to_y(buffz,fsort,gp)
      
      idx = 1
      do k = 1,gp%zsz(3)
         do j = 1,gp%zsz(2)
            do i = 1,gp%zsz(1)
               buffz(i,j,k) = sort_z(idx)%zpos
               idx = idx + 1
            end do 
         end do 
      end do
      call transpose_z_to_y(buffz,zCoords)
  
      idx = 1
      do k = 1,gp%ysz(3)
         do j = 1,gp%ysz(2)
            do i = 1,gp%ysz(1)
               sort_y(idx)%value = fsort(i,j,k)
               idx = idx + 1
            end do 
         end do 
      end do 
      
      idx = 1
      do k = 1,gp%ysz(3)
         do j = 1,gp%ysz(2)
            do i = 1,gp%ysz(1)
               sort_y(idx)%zpos = zCoords(i,j,k)
               idx = idx + 1
            end do 
         end do 
      end do 

      if (iiter .ne. niter) then
         call Qsort(sort_y, gp%ysz(1)*gp%ysz(2)*gp%ysz(3))
         idx = 1
         do k = 1,gp%ysz(3)
            do j = 1,gp%ysz(2)
               do i = 1,gp%ysz(1)
                  fsort(i,j,k) = sort_y(idx)%value
                  idx = idx + 1
               end do 
            end do 
         end do

         idx = 1
         do k = 1,gp%ysz(3)
            do j = 1,gp%ysz(2)
               do i = 1,gp%ysz(1)
                  zCoords(i,j,k) = sort_y(idx)%zpos
                  idx = idx + 1
               end do 
            end do 
         end do
      end if 
      
      call toc()
   
      PE_loc = zero
      do k = 1,gp%ysz(3)
         PE_loc = PE_loc - sum(fsort(:,:,k)*z(k))
      end do 
      PE = p_sum(PE_loc)
  
      if (nrank == 0) print*, "PE:", PE
   end do 

   call decomp_2d_write_one(2,fsort,"/scratch/globus/aghate/dump/fsort.dat", gp)
   call decomp_2d_write_one(2,zCoords,"/scratch/globus/aghate/dump/zcoords.dat", gp)
   deallocate(sort_y, sort_z) 
   deallocate(f, fsort, zCoords, buffz)
   deallocate(x, y, z)
   call MPI_Finalize(ierr)

end program 
