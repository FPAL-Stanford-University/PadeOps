program test_KS_inflow
   use kind_parameters, only: rkind
   use decomp_2d
   use mpi

   implicit none
   !real(rkind), parameter :: Lx = two*pi, Ly = two*pi, Lz = two*pi
   !integer, parameter :: nxQH = 1, nyQH = 4, nzQH = 4
   !integer, parameter :: xLESpQH = 64, yLESpQH = 64, zLESpQH = 64
   !integer, parameter :: nx = nxQH*xLESpQH, ny = nyQH*yLESpQH, nz = nzQH*zLESpQH
   !real(rkind), parameter :: dx = Lx/real(nx,rkind), dy = Ly/real(ny,rkind), dz = Lz/real(nz,rkind) 
  
   !type(decomp_info) :: gp
   !real(rkind), dimension(:,:,:), allocatable :: u, v, w
   !real(rkind), dimension(:), allocatable :: xG, yG, zG
   !integer :: ierr, pcol = 0, prow = 0
   !
   !integer :: i, j, k

   !type(rdt), dimension(:,:,:), allocatable :: QH
   !real(rkind), dimension(2) :: xLims, yLim, zLims
   !integer, parameter :: nk = 100, ntheta = 100
   !real(rkind) :: kmin = 5.d0, kmax = 10.d0
   !integer :: seed1 = 31345, seed2 = 43245, seed3 = 78647, s1, s2, s3
   !real(rkind), parameter :: eps_nufft = 1.d-8

   !call MPI_Init(ierr)
   !call decomp_2d_init(nx, ny, nz, prow, pcol)
   !call get_decomp_info(gp)

   !allocate(u(gp%xsz(1), gp%xsz(2), gp%xsz(3)))
   !allocate(v(gp%xsz(1), gp%xsz(2), gp%xsz(3)))
   !allocate(w(gp%xsz(1), gp%xsz(2), gp%xsz(3)))

   !allocate(QH(nxQH, nyQH, nzQH))

   !do k = 1,nzQH
   !   do j = 1,nyQH
   !      do i = 1,nxQH
   !         !s1 = seed1 + 12*i + 123*j + 1234*k; s2 = seed2 + 12*i + 123*j + 1234*k; s3 = seed3 + 12*i + 123*j + 1234*k
   !         !nxin = nx; nyin = ; nzin = 
   !         !call QH(i,j,k)%init(nk, ntheta, kmin, kmax, .false., s1, s2, s3, nxin, nyin, nzin, &
   !         !                     & xLims, yLims, zLims, eps_nufft, .false.)
   !      end do 
   !   end do 
   !end do 


   !deallocate(u, v, w)
   !call MPI_Finalize(ierr)
end program
