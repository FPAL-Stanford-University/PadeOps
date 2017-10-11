module ModifiedWaveNumbers
   use kind_parameters, only: rkind
   implicit none
contains

subroutine GetKmod_Fourier(kinout)
   real(rkind), dimension(:), intent(inout) :: kinout

   kinout = kinout
end subroutine 

subroutine GetKmod_CD06_stagg(kinout)
   use constants, only: one, two, three, five 
   real(rkind), dimension(:), intent(inout) :: kinout
   real(rkind), parameter :: alpha = 9._rkind/62._rkind
   real(rkind), parameter :: beta = 0._rkind
   real(rkind), parameter :: a = 63._rkind/62._rkind
   real(rkind), parameter :: b = 17._rkind/62._rkind
   real(rkind), parameter :: c = 0._rkind
    
    kinout = (two*a*sin(kinout/two) + (two/three)*b*sin(three*kinout/two) + &
         (two/five)*c*sin(five*kinout/two))/(one + two*alpha*cos(kinout) + &
          two*beta*cos(two*kinout))

end subroutine 


end module 


program test_PoissonPeriodic
   use kind_parameters, only: rkind
   use constants, only: pi, two 
   use PoissonPeriodicMod, only: PoissonPeriodic
   use ModifiedWaveNumbers, only: GetKmod_Fourier, GetKmod_CD06_stagg
   use decomp_2d
   use mpi
   use reductions, only: p_maxval, p_sum
   use exits, only: GracefulExit

   implicit none
   integer, parameter  :: dir_id = 3                               ! Data in decomposition, x: 1, y: 2, z: 3
   integer, parameter :: nx = 64, ny = 32, nz = 16                 ! number of points in x, y, z
   real(rkind), parameter :: Lx = two*pi, Ly = two*pi, Lz = two*pi ! Domain size 
   integer, parameter :: l = 6, m = 3, n = 1                       ! number of waves in x, y, z directions
   integer, parameter :: numscheme = 0                             ! 0: Fourier collocation, 1: 6th order stagg compact FD
   logical, parameter :: useExhaustiveFFT = .false.                ! search exhaustively for the best FFTW algorithm?  

   real(rkind), dimension(:,:,:), allocatable :: x, y, z, f, rhs, ftrue
   type(decomp_info) :: gp
   type(PoissonPeriodic) :: poiss
   integer :: i, j, k, ii, jj, kk, ierr 
   real(rkind) :: dx, dy, dz, maxerror
   integer :: st(3), en(3), sz(3)

   call MPI_Init(ierr)
   
   call decomp_2d_init(nx, ny, nz, 0, 0)
   call get_decomp_info(gp)
    

   dx = Lx/nx
   dy = Ly/ny
   dz = Lz/nz

   select case (numscheme)
   case (0)
      call poiss%init(dx, dy, dz, gp, dir_id, useExhaustiveFFT, Get_ModKx=GetKmod_Fourier, &
                     & Get_ModKy=GetKMod_Fourier, Get_ModKz=GetKMod_Fourier) 
   case (1)
      call poiss%init(dx, dy, dz, gp, dir_id, useExhaustiveFFT, Get_ModKx=GetKmod_CD06_stagg, &
                     & Get_ModKy=GetKMod_CD06_stagg, Get_ModKz=GetKMod_CD06_stagg) 
   case default
      call GracefulExit("Incorrect choice of numerical scheme",312)
   end select 


   select case (dir_id)
   case(1)
      st = gp%xst; en = gp%xen; sz = gp%xsz
   case(2)
      st = gp%yst; en = gp%yen; sz = gp%ysz
   case(3)
      st = gp%zst; en = gp%zen; sz = gp%zsz
   end select 

   allocate(x(sz(1),sz(2),sz(3)))
   allocate(y(sz(1),sz(2),sz(3)))
   allocate(z(sz(1),sz(2),sz(3)))
   allocate(f(sz(1),sz(2),sz(3)))
   allocate(ftrue(sz(1),sz(2),sz(3)))
   allocate(rhs(sz(1),sz(2),sz(3)))
   kk = 1
   do k = st(3),en(3) 
       jj = 1
       do j = st(2),en(2) 
           ii = 1
           do i = st(1),en(1) 
               x(ii,jj,kk) = (i-1)*dx
               y(ii,jj,kk) = (j-1)*dy
               z(ii,jj,kk) = (k-1)*dz
               ii = ii + 1
           end do 
           jj = jj + 1
       end do 
       kk = kk + 1
   end do 

   ftrue = sin(two*pi*l*x/Lx)*cos(two*pi*m*y/Ly)*sin(two*pi*n*z/Lz)
   rhs = -(4*pi**2*cos((2*pi*m*y)/Ly)*sin((2*pi*l*x)/Lx)*sin((2*pi*n*z)/Lz)*(Lx**2*Ly**2*n**2 &
                 &  + Lx**2*Lz**2*m**2 + Ly**2*Lz**2*l**2))/(Lx**2*Ly**2*Lz**2);
   

   call poiss%poisson_solve(rhs, f)

   maxerror = p_maxval(maxval(abs(f - ftrue)))

   if (nrank == 0) print*, "Max error:", maxerror

   deallocate(x, y, z, f, ftrue, rhs)
   call poiss%destroy()
   call MPI_Finalize(ierr)

end program
