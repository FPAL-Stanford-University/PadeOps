module igrid_Operators
   use kind_parameters, only: rkind, clen
   use spectralMod, only: spectral  
   use decomp_2d
   use cd06staggstuff, only: cd06stagg
   use reductions, only: p_sum
   use PadeDerOps, only: Pade6stagg 
   use PoissonPeriodicMod, only: PoissonPeriodic
   use exits, only: GracefulExit, message
   use gaussianstuff, only: gaussian  
   use PadePoissonMod, only: padepoisson
   implicit none
   
   private
   public :: igrid_ops

   type :: igrid_ops
      private
      complex(rkind), dimension(:,:,:), allocatable, public :: cbuffy1, cbuffy2, cbuffy3
      real(rkind),    dimension(:,:,:), allocatable :: rbuffx, rbuffy, rbuffz1, rbuffz2
      type(decomp_info), public :: gp, gpE
      type(spectral), public  :: spect, spectE
      type(Pade6stagg), public :: derZ
      type(cd06stagg) :: derZ1d
      real(rkind), dimension(:,:,:), allocatable :: zarr1d_1, zarr1d_2

      type(PoissonPeriodic) :: poiss_periodic
      real(rkind) :: dx, dy, dz, Lx, Ly, Lz
      character(len=clen) ::  inputdir, outputdir, RestartDir
      real(rkind) :: mfact_xy
    
      logical :: PeriodicInZ, PoissonSolverInitiatized = .false. 
      integer :: RunID, vfilt_times 
    
      real(rkind), dimension(:), allocatable :: gxfilt, gyfilt
      type(gaussian) :: gfilt 
      contains
         procedure :: init
         procedure :: destroy
         procedure :: ddx
         procedure :: ddy
         procedure :: ddz
         procedure :: ddz_1d
         procedure :: d2dz2_1d 
         procedure :: d2dz2
         procedure :: TakeMean_xy
         procedure :: TakeMean_y
         procedure :: getFluct_from_MeanZ
         procedure :: ReadField3D
         procedure :: WriteField3D
         procedure :: WriteSummingRestart 
         procedure :: ReadSummingRestart 
         procedure :: WriteSummingRestartInfo
         procedure :: ReadSummingRestartInfo
         procedure :: allocate3Dfield
         procedure :: getCenterlineQuantity
         procedure :: WriteASCII_2D
         procedure :: dump_plane
         procedure :: getSimTime
         procedure :: getVolumeIntegral
         procedure :: getGradient
         procedure :: getCurl
         procedure :: check_dump_existence
         procedure :: initPoissonSolver
         procedure :: PoissonSolvePeriodic_inplace
         procedure :: dealias
         procedure :: softdealias
         procedure :: alloc_cbuffz
         procedure :: Read_VizSummary
         procedure :: alloc_zvec
         procedure :: initFilter
         procedure :: FilterField
         procedure :: Project_DivergenceFree_BC
     end type 

contains

subroutine project_DivergenceFree_BC(this, u, v, w, poiss)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(inout) :: u, v, w 
   class(padepoisson), intent(inout) :: poiss
    
   complex(rkind), dimension(this%spect%spectdecomp%ysz(1),this%spect%spectdecomp%ysz(2),this%spect%spectdecomp%ysz(3)) :: uhat, vhat
   complex(rkind), dimension(this%spectE%spectdecomp%ysz(1),this%spectE%spectdecomp%ysz(2),this%spectE%spectdecomp%ysz(3)) :: what
   complex(rkind), dimension(this%spect%spectdecomp%zsz(1),this%spect%spectdecomp%zsz(2),this%spect%spectdecomp%zsz(3)) :: tmp1
   complex(rkind), dimension(this%spectE%spectdecomp%zsz(1),this%spectE%spectdecomp%zsz(2),this%spectE%spectdecomp%zsz(3)) :: tmp2

   call this%spect%fft(u, uhat)
   call this%spect%fft(w, vhat) !< this is intentional to save memory
   call transpose_y_to_z(vhat, tmp1, this%spect%spectdecomp)
   call this%derZ%interpz_C2E(tmp1,tmp2,0,0)
   call transpose_z_to_y(tmp2, what, this%spectE%spectdecomp)
   call this%spect%fft(v, vhat)
   call poiss%PressureProjection(uhat, vhat, what)
   call this%spect%ifft(uhat, u)
   call this%spect%ifft(vhat, v)
   call transpose_y_to_z(what, tmp2, this%spectE%spectdecomp)
   call this%derZ%interpz_E2C(tmp2,tmp1,-1,-1)
   call transpose_z_to_y(tmp1, vhat, this%spect%spectdecomp)
   call this%spect%ifft(vhat, w)


end subroutine 

subroutine FilterField(this, f, fout)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in) :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: fout
   integer :: fid, j, k

   call this%spect%fft(f,this%cbuffy1)
   do k = 1,size(this%cbuffy1,3)
       do j = 1,size(this%cbuffy1,2)
          this%cbuffy1(:,j,k) = this%gyfilt(j)*this%gxfilt*this%cbuffy1(:,j,k)
       end do 
   end do 
    
   call this%spect%ifft(this%cbuffy1,fout)

   call transpose_x_to_y(fout,this%rbuffy,this%gp)
   call transpose_y_to_z(this%rbuffy,this%rbuffz1,this%gp)
   do fid = 1,this%vfilt_times
        call this%gfilt%filter3(this%rbuffz1,this%rbuffz2,size(this%rbuffz1,1),size(this%rbuffz1,2))
        this%rbuffz1 = this%rbuffz2
   end do 
   call transpose_z_to_y(this%rbuffz2,this%rbuffy,this%gp)
   call transpose_y_to_x(this%rbuffy,fout,this%gp)

   

end subroutine 


subroutine initFilter(this, nx_filt, ny_filt, vfilt_times) 
    use constants, only: pi
    class(igrid_ops), intent(inout) :: this
    integer, intent(in) :: nx_filt, ny_filt, vfilt_times
    real(rkind) :: kx_co, ky_co, dxf, dyf
    integer :: i, j, ierr

    allocate(this%gxfilt(this%spect%spectdecomp%ysz(1)))
    allocate(this%gyfilt(this%spect%spectdecomp%ysz(2)))
    

    dxf = this%Lx/nx_filt
    dyf = this%Ly/ny_filt
    kx_co = ((2.d0/3.d0)*pi/dxf)  
    ky_co = ((2.d0/3.d0)*pi/dyf)  

    do i = 1,size(this%gxfilt)
        if (abs(this%spect%k1(i,1,1)) < kx_co) then
            this%gxfilt(i) = 1.d0 
        else
            this%gxfilt(i) = 0.d0 
        end if
    end do 
    
    do j = 1,size(this%gyfilt)
        if (abs(this%spect%k2(1,j,1)) < ky_co) then
            this%gyfilt(j) = 1.d0 
        else
            this%gyfilt(j) = 0.d0 
        end if
    end do 

    this%vfilt_times = vfilt_times
    
    ierr = this%gfilt%init(this%gp%zsz(3),.false.) 
    
end subroutine 

subroutine alloc_zvec(this, vec)
   class(igrid_ops), intent(in) :: this
   real(rkind), dimension(:), allocatable, intent(out) :: vec

   allocate(vec(this%gp%zsz(3)))

end subroutine 

subroutine Read_VizSummary(this, times, timesteps)
   class(igrid_ops), intent(inout) :: this
   integer, intent(out), dimension(:), allocatable :: timesteps
   real(rkind), intent(out), dimension(:), allocatable :: times
   character(len=clen) :: tempname, fname 
   logical :: exists
   integer :: nr, i, ierr

   write(tempname,"(A3,I2.2,A12,A4)") "Run",this%runID, "_vis_summary",".smm"
   fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
   inquire(file=fname, exist=exists)
   if (exists) then
       open(unit=10,file=fname,access='sequential',form='formatted')
       nr = 0
       do 
           read(10,*, iostat = ierr)
           if (ierr == -1) exit       
           nr = nr + 1
       end do 
       rewind(10)
       
       allocate(times(nr), timesteps(nr))
       do i=1, nr
           read (10, *) times(i), timesteps(i) 
       end do     
    else
       call GracefulExit("Summary file not found.", 34)
    end if 


end subroutine 

subroutine dealias(this, field)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(inout) :: field 

   call this%spect%fft(field,this%cbuffy1)
   call this%spect%dealias(this%cbuffy1)
   call this%spect%ifft(this%cbuffy1,field)

end subroutine 

subroutine softdealias(this, field)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(inout) :: field 

   call this%spect%fft(field,this%cbuffy1)
   call this%spect%softdealias(this%cbuffy1)
   call this%spect%ifft(this%cbuffy1,field)

end subroutine 

subroutine initPoissonSolver(this, dx, dy, dz)
   class(igrid_ops), intent(inout) :: this
   real(rkind), intent(in) :: dx, dy, dz

   if (this%PeriodicInZ) then
      call this%poiss_periodic%init(dx, dy, dz, this%gp, 1, useExhaustiveFFT=.true.)
      call message(0,"Poisson solver initialized.")  
      this%PoissonSolverInitiatized = .true.  
   else
      call GracefulExit("Poisson solver only supported for periodic problems",134)
   end if

end subroutine 

subroutine PoissonSolvePeriodic_inplace(this,rhs, dealiasRHS)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(inout) :: rhs
   logical, intent(in), optional :: dealiasRHS

   if (this%PoissonSolverInitiatized) then
      if (present(dealiasRHS)) then
          if (dealiasRHS) then
              call this%dealias(rhs)
          end if 
      end if 
      call this%poiss_periodic%poisson_solve(rhs)
   else
      call message(0,"WARNING: Poisson solver called without initializing.")
   end if 
end subroutine 

subroutine GetGradient(this, f, dfdx, dfdy, dfdz, botBC, topBC)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: dfdx, dfdy, dfdz
   integer, intent(in) :: botBC, topBC
   
   call this%spect%fft(f,this%cbuffy1)
   call this%spect%mtimes_ik1_oop(this%cbuffy1, this%cbuffy2)
   call this%spect%ifft(this%cbuffy2,dfdx,setOddball=.True.)
   call this%spect%mtimes_ik2_ip(this%cbuffy1)
   this%cbuffy1(:,size(this%cbuffy1,2)/2+1,:) = 0.d0 ! need to set the oddball to zero 
   call this%spect%ifft(this%cbuffy1,dfdy)
   call this%ddz(f, dfdz, botBC, topBC)

end subroutine 

function getVolumeIntegral(this, arr) result(val)
   class(igrid_ops), intent(in) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in) :: arr
   real(rkind) :: val 

   val = this%dx*this%dy*this%dz*p_sum(arr)

end function 

subroutine getCurl(this, u, v, w, omegax, omegay, omegaz, botBCu, topBCu, botBCv, topBCv)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: u, v, w
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: omegax, omegay, omegaz
   integer, intent(in) :: botBCu, topBCu, botBCv, topBCv

   call this%ddy(w,omegax)
   call this%ddz(u,omegay,botBCu,topBCu)
   call this%ddx(v,omegaz)

   call this%ddz(v,this%rbuffx,botBCv,topBCv) 
   omegax = omegax - this%rbuffx

   call this%ddx(w,this%rbuffx)
   omegay = omegay - this%rbuffx

   call this%ddy(u,this%rbuffx)
   omegaz = omegaz - this%rbuffx

end subroutine
   
function getCenterlineQuantity(this, vec) result(val)
   class(igrid_ops), intent(in) :: this
   real(rkind), dimension(this%gp%zsz(3)), intent(in) :: vec
   integer :: nz
   real(rkind) :: val

   nz = size(vec)

   val = 0.5d0*(vec(nz/2) + vec(nz/2+1))

end function

subroutine init(this, nx, ny, nz, dx, dy, dz, InputDir, OutputDir, RunID, isPeriodicinZ, NumericalSchemeZ, RestartDir)
   class(igrid_ops), intent(out), target :: this
   integer, intent(in) :: nx, ny, nz
   real(rkind), intent(in) :: dx, dy, dz
   character(len=clen), intent(in) ::  inputdir, outputdir
   character(len=clen), intent(in), optional :: RestartDir
   logical, intent(in) :: isPeriodicinZ
   integer, intent(in) :: RunID, NumericalSchemeZ
   logical, dimension(3) :: periodicbcs

   periodicbcs(1) = .true.; periodicbcs(2) = .true.; periodicbcs(3) = isPeriodicinZ
   call decomp_2d_init(nx, ny, nz, 0, 0, periodicbcs)
   call get_decomp_info(this%gp)
   
   call decomp_info_init(nx,ny,nz+1,this%gpE)

   call this%spect%init("x",nx,ny,nz,dx, dy, dz, "FOUR", "2/3rd", 2 , fixOddball=.false., &
                  exhaustiveFFT=.TRUE., init_periodicInZ=isPeriodicinZ, dealiasF=(2.d0/3.d0))
   call this%spectE%init("x",nx,ny,nz+1,dx, dy, dz, "FOUR", "2/3rd", 2 , fixOddball=.false., &
                  exhaustiveFFT=.TRUE., init_periodicInZ=.FALSE., dealiasF=(2.d0/3.d0))
  
   call this%derZ%init(this%gp, this%spect%spectdecomp, this%gpE, this%spectE%spectdecomp, dz, NumericalSchemeZ, isPeriodicinZ, this%spect)
   
   call this%derZ1d%init(nz, dz, isTopEven = .false., isBotEven = .false., &
                                   isTopsided = .true., isBotSided = .true.)
 
   call this%spect%alloc_r2c_out(this%cbuffy1) 
   call this%spect%alloc_r2c_out(this%cbuffy2) 

   this%Lx = nx*dx
   this%Ly = ny*dy
   this%Lz = nz*dz

   this%inputdir  = inputdir
   this%outputdir = outputdir
   this%RunID = RunID
   if (present(RestartDir)) this%RestartDir = RestartDir

   this%dx = dx
   this%dy = dy
   this%dz = dz
    
   this%PeriodicInZ = isPeriodicInZ

   this%mfact_xy = 1.d0/(real(nx,rkind)*real(ny,rkind))
   allocate(this%rbuffy (this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3)))
   allocate(this%rbuffz1(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
   allocate(this%rbuffz2(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
   allocate(this%rbuffx (this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)))

   allocate(this%zarr1d_1(1,1,nz))
   allocate(this%zarr1d_2(1,1,nz))
end subroutine 

subroutine destroy(this)
   class(igrid_ops), intent(inout), target :: this
  
   deallocate(this%rbuffx, this%rbuffy, this%rbuffz1, this%rbuffz2)
   deallocate(this%zarr1d_1, this%zarr1d_2)
end subroutine 

subroutine ddx(this,f, dfdx)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: dfdx
  
   call this%spect%fft(f,this%cbuffy1)
   call this%spect%mtimes_ik1_ip(this%cbuffy1)
   call this%spect%ifft(this%cbuffy1,dfdx,setOddball=.true.)
end subroutine 

subroutine ddy(this,f, dfdy)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: dfdy
  
   call this%spect%fft(f,this%cbuffy1)
   call this%spect%mtimes_ik2_ip(this%cbuffy1)
   this%cbuffy1(:,size(this%cbuffy1,2)/2+1,:) = 0.d0 ! need to set the oddball to zero 
   call this%spect%ifft(this%cbuffy1,dfdy)
end subroutine 

subroutine ddz(this, f, dfdz, botBC, topBC)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: dfdz
   integer, intent(in) :: botBC, topBC

   call transpose_x_to_y(f,this%rbuffy,this%gp)
   call transpose_y_to_z(this%rbuffy,this%rbuffz1,this%gp)
   call this%derZ%ddz_C2C(this%rbuffz1,this%rbuffz2, botBC, topBC)
   call transpose_z_to_y(this%rbuffz2,this%rbuffy,this%gp)
   call transpose_y_to_x(this%rbuffy,dfdz,this%gp)
end subroutine 

subroutine d2dz2(this, f, d2fdz2, botBC, topBC)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: d2fdz2
   integer, intent(in) :: botBC, topBC

   call transpose_x_to_y(f,this%rbuffy,this%gp)
   call transpose_y_to_z(this%rbuffy,this%rbuffz1,this%gp)
   call this%derZ%d2dz2_C2C(this%rbuffz1,this%rbuffz2, botBC, topBC)
   call transpose_z_to_y(this%rbuffz2,this%rbuffy,this%gp)
   call transpose_y_to_x(this%rbuffy,d2fdz2,this%gp)
end subroutine 

subroutine ddz_1d(this, f1d, dfdz1d, botBC, topBC)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%zsz(3)), intent(in)  :: f1d
   real(rkind), dimension(this%gp%zsz(3)), intent(out) :: dfdz1d
   integer, intent(in), optional :: botBC, topBC

   this%zarr1d_1(1,1,:) = f1d
   if (present(botBC) .and. present(topBC)) then
        call this%derZ%ddz_1d_C2C(this%zarr1d_1,this%zarr1d_2,botBC,topBC)
   else
        call this%derZ%ddz_1d_C2C(this%zarr1d_1,this%zarr1d_2,0,0)
   end if 
   dfdz1d = this%zarr1d_2(1,1,:)
end subroutine 


subroutine d2dz2_1d(this, f1d, d2fdz21d, botBC, topBC)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%zsz(3)), intent(in)  :: f1d
   real(rkind), dimension(this%gp%zsz(3)), intent(out) :: d2fdz21d
   integer, intent(in) :: botBC, topBC

   this%zarr1d_1(1,1,:) = f1d
   call this%derZ%ddz_1d_C2C(this%zarr1d_1,this%zarr1d_2,botBC,topBC)
   call this%derZ%ddz_1d_C2C(this%zarr1d_2,this%zarr1d_1,-botBC,-topBC)
   d2fdz21d = this%zarr1d_1(1,1,:)
end subroutine 


subroutine getFluct_from_MeanZ(this, f, ffluct)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: ffluct 

   real(rkind) :: mnZ
   integer :: k

   call transpose_x_to_y(f,this%rbuffy,this%gp)
   call transpose_y_to_z(this%rbuffy,this%rbuffz1,this%gp)
   do k = 1,this%gp%zsz(3)
      mnZ = p_sum(sum(this%rbuffz1(:,:,k)))*this%mfact_xy
      this%rbuffz2(:,:,k) = this%rbuffz1(:,:,k) - mnZ
   end do 
   call transpose_z_to_y(this%rbuffz2,this%rbuffy,this%gp)
   call transpose_y_to_x(this%rbuffy,ffluct,this%gp)

end subroutine 

subroutine TakeMean_xy(this, f, fmean)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%zsz(3)), intent(out) :: fmean 
   integer :: k

   call transpose_x_to_y(f,this%rbuffy,this%gp)
   call transpose_y_to_z(this%rbuffy,this%rbuffz1,this%gp)
   do k = 1,this%gp%zsz(3)
      fmean(k) = p_sum(sum(this%rbuffz1(:,:,k)))*this%mfact_xy
   end do 
end subroutine

subroutine TakeMean_y(this, f, fmean)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: fmean
   integer :: i,k 

   call transpose_x_to_y(f,this%rbuffy,this%gp)
   do k = 1,this%gp%ysz(3)
      do i = 1,this%gp%ysz(1)
        this%rbuffy(i,:,k) = sum(this%rbuffy(i,:,k))/this%gp%ysz(2) 
      end do
   end do
   call transpose_y_to_x(this%rbuffy,fmean,this%gp)

end subroutine

subroutine ReadField3D(this, field, label, tidx)
   use decomp_2d_io
   use exits, only: gracefulExit
   use mpi
   class(igrid_ops), intent(in) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out)  :: field
   character(len=clen) :: tempname, fname
   character(len=4), intent(in) :: label
   integer, intent(in) :: tidx
   integer :: ierr 

   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",tidx,".out"
   fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
   
   open(777,file=trim(fname),status='old',iostat=ierr)
   if(ierr .ne. 0) then
    print*, "Rank:", nrank, ". File:", fname, " not found"
    call mpi_barrier(mpi_comm_world, ierr)
    call gracefulExit("File I/O issue.",44)
   end if 
   close(777)
   call decomp_2d_read_one(1,field,fname,this%gp)
end subroutine  


function check_dump_existence(this, label, tidx) result(file_found)
   class(igrid_ops), intent(in) :: this
   character(len=clen) :: tempname, fname
   character(len=4), intent(in) :: label
   integer, intent(in) :: tidx
   integer :: ierr 
   logical :: file_found 

   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",tidx,".out"
   fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
   
   open(777,file=trim(fname),status='old',iostat=ierr)
   if (ierr == 0) then
       file_found = .true. 
   else
       file_found = .false. 
   end if 

end function


subroutine WriteField3D(this, field, label, tidx)
   use decomp_2d_io
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: field
   character(len=clen) :: tempname, fname
   character(len=4), intent(in) :: label
   integer, intent(in) :: tidx
         
   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",tidx,".out"
   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
   
   call decomp_2d_write_one(1,field,fname,this%gp)
end subroutine  


subroutine WriteSummingRestart(this, field, label, tidx)
   use decomp_2d_io
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: field
   character(len=clen) :: tempname, fname
   character(len=4), intent(in) :: label
   integer, intent(in) :: tidx
         
   write(tempname,"(A14,I2.2,A1,A4,A2,I6.6,A4)") "SummingRESTART",this%runID, "_",label,"_t",tidx,".rss"
   fname = this%RestartDir(:len_trim(this%RestartDir))//"/"//trim(tempname)
   call decomp_2d_write_one(1,field,fname,this%gp)
end subroutine  

subroutine ReadSummingRestart(this, field, label, tidx)
   use decomp_2d_io
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out)  :: field
   character(len=clen) :: tempname, fname
   character(len=4), intent(in) :: label
   integer, intent(in) :: tidx
   integer :: ierr
         
   write(tempname,"(A14,I2.2,A1,A4,A2,I6.6,A4)") "SummingRESTART",this%runID, "_",label,"_t",tidx,".rss"
   fname = this%RestartDir(:len_trim(this%RestartDir))//"/"//trim(tempname)
   
   open(777,file=trim(fname),status='old',iostat=ierr)
   if(ierr .ne. 0) then
    print*, "Rank:", nrank, ". File:", fname, " not found"
    call gracefulExit("File I/O issue.",44)
   end if 
   close(777)
   
   call decomp_2d_read_one(1,field,fname,this%gp)
end subroutine  

subroutine WriteSummingRestartInfo(this,tidx,nsum)
   class(igrid_ops), intent(inout) :: this
   integer, intent(in) :: tidx, nsum
   character(len=clen) :: tempname, fname
   if (nrank == 0) then
       write(tempname,"(A14,I2.2,A1,A4,A2,I6.6,A4)") "SummingRESTART",this%runID, "_","info","_t",tidx,".rss"
       fname = this%RestartDir(:len_trim(this%RestartDir))//"/"//trim(tempname)
       OPEN(UNIT=10, FILE=trim(fname))
       write(10,"(I7.7)") nsum
       close(10)
   end if 
end subroutine

subroutine ReadSummingRestartInfo(this,tidx,nsum)
   class(igrid_ops), intent(inout) :: this
   integer, intent(in) :: tidx
   integer, intent(out) :: nsum
   character(len=clen) :: tempname, fname
   integer :: ierr

   if (nrank == 0) then
       write(tempname,"(A14,I2.2,A1,A4,A2,I6.6,A4)") "SummingRESTART",this%runID, "_","info","_t",tidx,".rss"
       fname = this%RestartDir(:len_trim(this%RestartDir))//"/"//trim(tempname)
       OPEN(UNIT=10, FILE=trim(fname),status='old',iostat=ierr)
       if (ierr .ne. 0) then
           print*, "Rank:", nrank, ". File:", fname, " not found"
           call gracefulExit("File I/O issue.",44)
       else
           read(10,"(I7.7)") nsum
       end if 
       close(10)
   end if 
end subroutine

subroutine allocate3Dfield(this, field)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(:,:,:), allocatable, intent(out) :: field
   allocate(field(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)))
end subroutine 

subroutine alloc_cbuffz(this, field)
   class(igrid_ops), intent(inout) :: this
   complex(rkind), dimension(:,:,:), allocatable, intent(out) :: field
   allocate(field(this%spect%spectdecomp%zsz(1),this%spect%spectdecomp%zsz(2),this%spect%spectdecomp%zsz(3)))
end subroutine 


function getSimTime(this, tidx) result(time)
   class(igrid_ops), intent(in) :: this
   character(len=clen) :: tempname, fname
   integer, intent(in) :: tidx
   real(rkind) :: time

   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_","info","_t",tidx,".out"
   fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
   open(unit=10,file=fname,access='sequential',form='formatted')
   read(10,*) time
   close(10)

end function

subroutine WriteASCII_2D(this, field, flabel)
   use basic_io, only: write_2d_ascii 
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(:,:), intent(in) :: field
   character(len=4), intent(in) :: flabel
   character(len=clen) :: tempname, fname
   
   write(tempname,"(A3,I2.2,A1,A4,A4)") "Run",this%runID, "_",flabel,".stt"
   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
   
   call write_2d_ascii(field, fname)

end subroutine 


subroutine dump_plane(this, field, dir_id, plane_id, tid, label)
   use decomp_2d_io
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: field
   integer, intent(in) :: dir_id, plane_id, tid
   character(len=4), intent(in) :: label
   character(len=clen) :: fname
   character(len=clen) :: tempname
  
   select case(dir_id) 
   case(1)
       write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A1,A4,A4)") "Run", this%RunID,"_t",tid,"_x",plane_id,"_",label,".pln"
   case(2)
       write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A1,A4,A4)") "Run", this%RunID,"_t",tid,"_y",plane_id,"_",label,".pln"
   case(3)
       write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A1,A4,A4)") "Run", this%RunID,"_t",tid,"_z",plane_id,"_",label,".pln"
   end select 
   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
   call decomp_2d_write_plane(1,field,dir_id, plane_id, fname, this%gp)

end subroutine 


end module 
