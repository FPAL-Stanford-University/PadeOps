module PoissonPeriodicMod
   use kind_parameters, only: rkind
   use decomp_2d
   use constants, only: pi, zero, one, three, five, two, imi
   use fft_3d_stuff, only: fft_3d
   use exits, only: GracefulExit 
   implicit none
   private 
   public :: PoissonPeriodic  

   interface 
      subroutine GetModifiedWavenum(kinout)
         import :: rkind
         real(rkind), dimension(:), intent(inout) :: kinout 
      end subroutine 
   end interface 

   type :: PoissonPeriodic
      type(decomp_info), pointer :: gp, sp_gp
      integer :: dir_id, nx_g, ny_g, nz_g
      integer :: nx_in, ny_in, nz_in
      integer :: nx_hat, ny_hat, nz_hat
      real(rkind), dimension(:), allocatable :: kx, ky, kz
      type(fft_3d) :: FT
      logical :: doIhaveZeroWaveNum 

      contains
         procedure :: init
         procedure :: destroy
         procedure, private :: poisson_solve_inplace
         procedure, private :: poisson_solve_outofplace
         procedure, private :: poisson3D_multiply
         generic :: poisson_solve => poisson_solve_inplace, poisson_solve_outofplace 
   end type
contains 

   subroutine poisson_solve_inplace(this, rhs)
      class(PoissonPeriodic), intent(inout) :: this
      real(rkind), dimension(this%nx_in,this%ny_in,this%nz_in), intent(inout) :: rhs
      
      real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)) :: rbuffinX
      complex(rkind), dimension(this%nx_hat, this%ny_hat, this%nz_hat) :: RHS_hat

      select case(this%dir_id) 
      case(1)
         call this%FT%fft3_x2z(rhs,RHS_hat)
         call this%poisson3D_multiply(RHS_hat)
         call this%FT%ifft3_z2x(RHS_hat,rhs)
      case(2)
         call transpose_y_to_x(rhs,rbuffinX,this%gp)
         call this%FT%fft3_x2z(rbuffinX,RHS_hat)
         call this%poisson3D_multiply(RHS_hat)
         call this%FT%ifft3_z2x(RHS_hat,rbuffinX)
         call transpose_x_to_y(rbuffinX, rhs, this%gp)
      case(3)
         call this%FT%fft3_z2x(rhs,RHS_hat)
         call this%poisson3D_multiply(RHS_hat)
         call this%FT%ifft3_x2z(RHS_hat,rhs)
      end select 
   end subroutine

   subroutine poisson_solve_outofplace(this, rhs, f)
      class(PoissonPeriodic), intent(inout) :: this
      real(rkind), dimension(this%nx_in,this%ny_in,this%nz_in), intent(in)  :: rhs
      real(rkind), dimension(this%nx_in,this%ny_in,this%nz_in), intent(out) :: f
      
      real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)) :: rbuffinX
      complex(rkind), dimension(this%nx_hat, this%ny_hat, this%nz_hat) :: RHS_hat

      select case(this%dir_id) 
      case(1)
         call this%FT%fft3_x2z(rhs,RHS_hat)
         call this%poisson3D_multiply(RHS_hat)
         call this%FT%ifft3_z2x(RHS_hat,f)
      case(2)
         call transpose_y_to_x(rhs,rbuffinX,this%gp)
         call this%FT%fft3_x2z(rbuffinX,RHS_hat)
         call this%poisson3D_multiply(RHS_hat)
         call this%FT%ifft3_z2x(RHS_hat,rbuffinX)
         call transpose_x_to_y(rbuffinX, f, this%gp)
      case(3)
         call this%FT%fft3_z2x(rhs,RHS_hat)
         call this%poisson3D_multiply(RHS_hat)
         call this%FT%ifft3_x2z(RHS_hat,f)
      end select 

   end subroutine 

   subroutine poisson3D_multiply(this, RHS_hat)
      use constants, only: imi
      class(PoissonPeriodic), intent(in) :: this
      complex(rkind), dimension(this%nx_hat, this%ny_hat, this%nz_hat), intent(inout) :: RHS_hat
      integer :: i, j, k
      real(rkind) :: Minus_One_by_kabs_sq, ky_sq, kz_sq
      complex(rkind), parameter :: czero = 0._rkind + imi*0._rkind

      do k = 1,this%nz_hat
         kz_sq = this%kz(k)*this%kz(k)
         do j = 1,this%ny_hat
            ky_sq = this%ky(j)*this%ky(j)
            !$omp simd 
            do i = 1,this%nx_hat
               Minus_One_by_kabs_sq = -1._rkind/(this%kx(i)*this%kx(i) + ky_sq + kz_sq + 1.d-20)
               RHS_hat(i,j,k) = RHS_hat(i,j,k)*Minus_One_by_kabs_sq 
            end do 
         end do 
      end do 

      if (this%doIhaveZeroWavenum) RHS_hat(1,1,1) = czero

   end subroutine 


   subroutine init(this, dx, dy, dz, gp, dir_id, useExhaustiveFFT, Get_ModKx, Get_ModKy, Get_ModKz)
      class(PoissonPeriodic), intent(inout) :: this
      real(rkind), intent(in) :: dx, dy, dz
      integer, intent(in) :: dir_id
      logical, intent(in), optional :: useExhaustiveFFT
      type(decomp_info), intent(in), target :: gp
      logical :: useExhaustiveFFT_
      procedure(GetModifiedWavenum), optional :: Get_ModKx, Get_ModKy, Get_ModKz

      real(rkind), dimension(:), allocatable :: ktmp 
      integer :: nx, ny, nz, ierr, xst, xen, yst, yen, zst, zen

      nx = gp%xsz(1)
      ny = gp%ysz(2)
      nz = gp%zsz(3)

      this%gp => gp

      if (present(useExhaustiveFFT)) then
         useExhaustiveFFT_ = useExhaustiveFFT
      else
         useExhaustiveFFT_ = .true.
      end if 

      ierr = 0
      select case (dir_id)
      case (1)
         ierr = this%FT%init(nx,ny,nz, "x",dx,dy,dz,useExhaustiveFFT_, &
                  & fixOddball_=.false., allocK=.false., dodealiasing=.false.)   
         this%nx_in = gp%xsz(1); this%ny_in = gp%xsz(2); this%nz_in = gp%xsz(3)
      
      case(2)
         ! NOTE: the "x" used for input here is NOT a typo - you need to
         ! transpose inputs to x decomposition before calling the 3d fft
         ierr = this%FT%init(nx,ny,nz, "x",dx,dy,dz,useExhaustiveFFT_, &
                 & fixOddball_=.false., allocK=.false., dodealiasing=.false.)  
         this%nx_in = gp%ysz(1); this%ny_in = gp%ysz(2); this%nz_in = gp%ysz(3)
      
      case(3)
         ierr = this%FT%init(nx,ny,nz, "z",dx,dy,dz,useExhaustiveFFT_, &
                 & fixOddball_=.false., allocK=.false., dodealiasing=.false.)   
         this%nx_in = gp%zsz(1); this%ny_in = gp%zsz(2); this%nz_in = gp%zsz(3)
      
      case default
         call GracefulExit("Incorrect option for DIR_ID", 31243)
      end select 
      
      if (ierr .ne. 0) then
         call GracefulExit("Couldn't initialize 3d FFT inside SPECTRAL derived type",123)
      end if 

      this%dir_id = dir_id
      call this%FT%get_complex_output_size(this%nx_hat, this%ny_hat, this%nz_hat)
      call this%FT%get_complex_output_start_end_indices(xst,xen,yst,yen,zst,zen)
      call this%FT%link_spectral_gp(this%sp_gp)

      allocate(this%kx(this%nx_hat))
      allocate(this%ky(this%ny_hat))
      allocate(this%kz(this%nz_hat))
      
      allocate(ktmp(nx))
      ktmp = GetWaveNums(nx, dx)
      this%kx = ktmp(xst:xen)
      deallocate(ktmp)

      allocate(ktmp(ny))
      ktmp = GetWaveNums(ny, dy)
      this%ky = ktmp(yst:yen)
      deallocate(ktmp)
      
      allocate(ktmp(nz))
      ktmp = GetWaveNums(nz, dz)
      this%kz = ktmp(zst:zen)
      deallocate(ktmp)

      if (present(Get_ModKx)) then
         this%kx = this%kx*dx
         call Get_modKx(this%kx)
         this%kx = this%kx/dx
      end if
      
      if (present(Get_ModKy)) then
         this%ky = this%ky*dy
         call Get_modKy(this%ky)
         this%ky = this%ky/dy
      end if

      if (present(Get_ModKz)) then
         this%kz = this%kz*dz
         call Get_modKz(this%kz)
         this%kz = this%kz/dz
      end if

      if ((xst == 1) .and. (yst == 1) .and. (zst == 1)) then
         this%doIhaveZeroWaveNum = .true.
      else
         this%doIhaveZeroWaveNum = .false.
      end if 

   end subroutine 


   subroutine destroy(this)
      class(PoissonPeriodic), intent(out) :: this
      if (allocated(this%kx)) deallocate(this%kx)
      if (allocated(this%ky)) deallocate(this%ky)
      if (allocated(this%kz)) deallocate(this%kz)
      call this%FT%destroy()

   end subroutine 



   pure function GetWaveNums(nx,dx) result(k)
      use constants, only: pi, two
      integer, intent(in) :: nx
      real(rkind), intent(in) :: dx
      real(rkind), dimension(nx) :: k

      integer :: i,dummy

      dummy = nx - MOD(nx,2)

      do i = 1,nx
          k(i) = ( -pi + (i-1)*two*pi/real(dummy,rkind) ) / dx
      end do

      k = ifftshift(k)

   end function
   
   pure function ifftshift(k) result(kshift)
      
      real(rkind), dimension(:), intent(in) :: k
      real(rkind), dimension(SIZE(k)) :: kshift
      integer :: n

      n = SIZE(k)

      select case ( MOD(n,2) )
      case (0)
          kshift(1:n/2) = k(n/2+1:n)
          kshift(n/2+1:n) = k(1:n/2)
      case (1)
          kshift(1:(n+1)/2) = k((n+1)/2:n)
          kshift((n+1)/2+1:n) = k(1:(n-1)/2)
      end select

   end function
   

end module 

