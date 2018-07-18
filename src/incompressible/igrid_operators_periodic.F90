module igrid_Operators_Periodic
   use kind_parameters, only: rkind, clen
   use spectralMod, only: spectral  
   use decomp_2d
   use cd06staggstuff, only: cd06stagg
   use reductions, only: p_sum
   use PoissonPeriodicMod, only: PoissonPeriodic 
   implicit none
   
   private
   public :: Ops_Periodic

   type :: Ops_Periodic
      private
      complex(rkind), dimension(:,:,:), allocatable :: cbuffy1, cbuffy2, cbuffz
      real(rkind),    dimension(:,:,:), allocatable :: rbuffy, rbuffz1, rbuffz2
      type(decomp_info), pointer :: gp
      type(spectral)  :: spect
      type(cd06stagg) :: derZ
      real(rkind), dimension(:,:,:), allocatable :: zarr1d_1, zarr1d_2

      character(len=clen) ::  inputdir, outputdir
      real(rkind) :: mfact_xy

      type(PoissonPeriodic) :: poiss

      contains
         procedure :: init
         procedure :: destroy
         procedure :: ddx
         procedure :: ddy
         procedure :: ddz
         procedure :: ddz_cmplx2cmplx
         procedure :: ReadField3D
         procedure :: WriteField3D
         procedure :: allocate3Dfield
         procedure :: SolvePoisson_oop
         procedure :: SolvePoisson_ip
         procedure :: dealiasField
         procedure :: link_spect
   end type 

contains



subroutine link_spect(this, spectout)
   class(Ops_Periodic), intent(in), target :: this
   type(spectral), pointer, intent(out) :: spectout

   spectout => this%spect

end subroutine


subroutine dealiasField(this, f)
   class(Ops_Periodic), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(inout)  :: f
   call this%spect%fft(f,this%cbuffy1)
   call this%spect%dealias(this%cbuffy1)
   call this%spect%ifft(this%cbuffy1,f)

end subroutine

subroutine GetKmod_Fourier(kinout)
   real(rkind), dimension(:), intent(inout) :: kinout

   kinout = kinout
end subroutine 

subroutine SolvePoisson_oop(this, rhs, p)
   class(Ops_Periodic), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: rhs
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: p

   call this%poiss%poisson_solve(rhs,p)
end subroutine 

subroutine SolvePoisson_ip(this, rhs)
   class(Ops_Periodic), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(inout)  :: rhs

   call this%poiss%poisson_solve(rhs)
end subroutine 

subroutine init(this, nx, ny, nz, dx, dy, dz, gp, InputDir, OutputDir)
   class(Ops_Periodic), intent(out), target :: this
   integer, intent(in) :: nx, ny, nz
   real(rkind), intent(in) :: dx, dy, dz
   class(decomp_info), intent(in), target :: gp
   character(len=clen), intent(in) ::  inputdir, outputdir

   this%gp => gp
   call this%spect%init("x",nx,ny,nz,dx, dy, dz, "four", "2/3rd", 2 , fixOddball=.false., &
                  exhaustiveFFT=.TRUE., init_periodicInZ=.TRUE., dealiasF=(2.d0/3.d0))
   
   call this%spect%alloc_r2c_out(this%cbuffy1) 
   allocate(this%cbuffz(this%spect%spectdecomp%zsz(1),this%spect%spectdecomp%zsz(2),this%spect%spectdecomp%zsz(3)))
   !call this%spect%alloc_r2c(this%cbuffy2) 
   
   this%inputdir  = inputdir
   this%outputdir = outputdir

   this%mfact_xy = 1.d0/(real(nx,rkind)*real(ny,rkind))
   allocate(this%rbuffy(gp%ysz(1),gp%ysz(2),gp%ysz(3)))
   allocate(this%rbuffz1(gp%zsz(1),gp%zsz(2),gp%zsz(3)))
   !allocate(this%rbuffz2(gp%zsz(1),gp%zsz(2),gp%zsz(3)))

   call this%poiss%init(dx, dy, dz, gp, 1, .true., Get_ModKx=GetKmod_Fourier, &
                     & Get_ModKy=GetKMod_Fourier, Get_ModKz=GetKMod_Fourier) 
end subroutine 

subroutine destroy(this)
   class(Ops_Periodic), intent(inout), target :: this
  
   deallocate(this%rbuffy, this%rbuffz1)
end subroutine 

subroutine ddx(this,f, dfdx)
   class(Ops_Periodic), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: dfdx
  
   call this%spect%fft(f,this%cbuffy1)
   call this%spect%mtimes_ik1_ip(this%cbuffy1)
   call this%spect%ifft(this%cbuffy1,dfdx)
end subroutine 

subroutine ddy(this,f, dfdy)
   class(Ops_Periodic), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: dfdy
  
   call this%spect%fft(f,this%cbuffy1)
   call this%spect%mtimes_ik2_ip(this%cbuffy1)
   call this%spect%ifft(this%cbuffy1,dfdy)
end subroutine 

subroutine ddz_cmplx2cmplx(this, fhat)
   class(Ops_Periodic), intent(inout) :: this
   complex(rkind), dimension(this%spect%spectdecomp%ysz(1),this%spect%spectdecomp%ysz(2),this%spect%spectdecomp%ysz(3)), intent(inout)  :: fhat

   call transpose_y_to_z(fhat, this%cbuffz,this%spect%spectdecomp)
   !call this%spect%ddz_C2C_complex_inplace(this%cbuffz)
   call this%spect%ddz_C2C_spect(this%cbuffz)
   call transpose_z_to_y(this%cbuffz, fhat, this%spect%spectdecomp)
end subroutine 



subroutine ddz(this, f, dfdz)
   class(Ops_Periodic), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: dfdz


   call transpose_x_to_y(f,this%rbuffy,this%gp)
   call transpose_y_to_z(this%rbuffy,this%rbuffz1,this%gp)
   !call this%spect%ddz_C2C_real_inplace(this%rbuffz1)
   call this%spect%ddz_C2C_spect(this%rbuffz1)
   call transpose_z_to_y(this%rbuffz1,this%rbuffy,this%gp)
   call transpose_y_to_x(this%rbuffy,dfdz,this%gp)
end subroutine 

subroutine ReadField3D(this, field, label, tidx, runID, newinputdir)
   use decomp_2d_io
   use exits, only: GracefulExit
   class(Ops_Periodic), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out)  :: field
   character(len=clen) :: tempname, fname
   character(len=4), intent(in) :: label
   integer, intent(in) :: tidx, runID
   character(len=clen), intent(in), optional :: newinputdir
   integer :: ierr

   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",runID, "_",label,"_t",tidx,".out"
   if (present(newinputdir)) then
      fname = newinputdir(:len_trim(newinputdir))//"/"//trim(tempname)
   else
      fname = this%InputDir(:len_trim(this%InputDir))//"/"//trim(tempname)
   end if
   open(777,file=fname,status='old',iostat=ierr)
   if (ierr .ne. 0) then
      if (nrank == 0) print*, "FILE NOT FOUND"
      if (nrank == 0) print*, fname
      call gracefulExit("File not found", 321)
   end if
   close(777)
   call decomp_2d_read_one(1,field,fname,this%gp)
end subroutine  

subroutine WriteField3D(this, field, label, tidx, runID, newOutputDir)
   use decomp_2d_io
   class(Ops_Periodic), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: field
   character(len=clen) :: tempname, fname
   character(len=4), intent(in) :: label
   integer, intent(in) :: tidx, runID
   character(len=clen), intent(in), optional :: newOutputdir
         
   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",runID, "_",label,"_t",tidx,".out"
   if (present(newOutputDir)) then
      fname = newOutputDir(:len_trim(newOutputDir))//"/"//trim(tempname)
   else
      fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
   end if
   call decomp_2d_write_one(1,field,fname,this%gp)
end subroutine  

subroutine allocate3Dfield(this, field)
   class(Ops_Periodic), intent(inout) :: this
   real(rkind), dimension(:,:,:), allocatable, intent(out) :: field
   allocate(field(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)))
end subroutine 


end module 
