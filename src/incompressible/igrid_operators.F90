module igrid_Operators
   use kind_parameters, only: rkind, clen
   use spectralMod, only: spectral  
   use decomp_2d
   use cd06staggstuff, only: cd06stagg
   use reductions, only: p_sum
   use PadeDerOps, only: Pade6stagg 

   implicit none
   
   private
   public :: igrid_ops

   type :: igrid_ops
      private
      complex(rkind), dimension(:,:,:), allocatable :: cbuffy1, cbuffy2
      real(rkind),    dimension(:,:,:), allocatable :: rbuffx, rbuffy, rbuffz1, rbuffz2
      type(decomp_info) :: gp, gpE
      type(spectral)  :: spect, spectE
      type(Pade6stagg) :: derZ
      type(cd06stagg) :: derZ1d
      real(rkind), dimension(:,:,:), allocatable :: zarr1d_1, zarr1d_2

      real(rkind) :: dx, dy, dz
      character(len=clen) ::  inputdir, outputdir
      real(rkind) :: mfact_xy

      integer :: RunID

      contains
         procedure :: init
         procedure :: destroy
         procedure :: ddx
         procedure :: ddy
         procedure :: ddz
         procedure :: ddz_1d
         procedure :: TakeMean_xy
         procedure :: getFluct_from_MeanZ
         procedure :: ReadField3D
         procedure :: WriteField3D
         procedure :: allocate3Dfield
         procedure :: getCenterlineQuantity
         procedure :: WriteASCII_2D
         procedure :: getSimTime
         procedure :: getVolumeIntegral
         procedure :: getGradient
         procedure :: getCurl
     end type 

contains

subroutine GetGradient(this, f, dfdx, dfdy, dfdz, botBC, topBC)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: dfdx, dfdy, dfdz
   integer, intent(in) :: botBC, topBC
   
   call this%spect%fft(f,this%cbuffy1)
   call this%spect%mtimes_ik1_oop(this%cbuffy1, this%cbuffy2)
   call this%spect%ifft(this%cbuffy2,dfdx)
   call this%spect%mtimes_ik2_ip(this%cbuffy1)
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

subroutine init(this, nx, ny, nz, dx, dy, dz, InputDir, OutputDir, RunID, isPeriodicinZ, NumericalSchemeZ)
   class(igrid_ops), intent(out), target :: this
   integer, intent(in) :: nx, ny, nz
   real(rkind), intent(in) :: dx, dy, dz
   character(len=clen), intent(in) ::  inputdir, outputdir
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
   
   this%inputdir  = inputdir
   this%outputdir = outputdir
   this%RunID = RunID

   this%dx = dx
   this%dy = dy
   this%dz = dz

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
   call this%spect%ifft(this%cbuffy1,dfdx)
end subroutine 

subroutine ddy(this,f, dfdy)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: dfdy
  
   call this%spect%fft(f,this%cbuffy1)
   call this%spect%mtimes_ik2_ip(this%cbuffy1)
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

subroutine ddz_1d(this, f1d, dfdz1d)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%zsz(3)), intent(in)  :: f1d
   real(rkind), dimension(this%gp%zsz(3)), intent(out) :: dfdz1d

   this%zarr1d_1(1,1,:) = f1d
   call this%derZ1d%ddz_C2C(this%zarr1d_1,this%zarr1d_2,1,1)
   dfdz1d = this%zarr1d_2(1,1,:)
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

subroutine ReadField3D(this, field, label, tidx)
   use decomp_2d_io
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out)  :: field
   character(len=clen) :: tempname, fname
   character(len=4), intent(in) :: label
   integer, intent(in) :: tidx
         
   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",tidx,".out"
   fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
   
   call decomp_2d_read_one(1,field,fname,this%gp)
end subroutine  

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

subroutine allocate3Dfield(this, field)
   class(igrid_ops), intent(inout) :: this
   real(rkind), dimension(:,:,:), allocatable, intent(out) :: field
   allocate(field(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)))
end subroutine 


function getSimTime(this, tidx) result(time)
   class(igrid_ops), intent(in) :: this
   character(len=clen) :: tempname, fname
   integer, intent(in) :: tidx
   real(rkind) :: time

   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_","info","_t",tidx,".out"
   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
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


end module 
