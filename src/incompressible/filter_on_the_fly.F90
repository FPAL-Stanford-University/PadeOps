module fof_mod
   use kind_parameters, only: rkind, clen
   use spectralMod, only: spectral
   use decomp_2d
   use decomp_2d_io
   use exits, only: message, gracefulExit


   implicit none
   private
   public :: fof


   type fof
      private
      character(len=clen) :: outputdir
      logical :: PeriodicInZ
      type(decomp_info), pointer :: gpC, sp_gpC
      integer :: fof_id, runID
      type(spectral), pointer :: spectC 
      real(rkind), dimension(:), allocatable :: Gfilt_x, Gfilt_y, Gfilt_z
      complex(rkind), dimension(:,:,:), pointer :: fhaty, fhatz
     
      integer :: sz_x, sz_y, sz_z
      contains 
         procedure :: init
         procedure :: destroy
         procedure :: filter_Complex2Real
         procedure :: filter_Real2Real
         procedure :: dumpFullField
         procedure :: dumpYPlanes 
         !procedure :: filter_uvw
         !procedure :: filter_field
         !procedure :: dump_filter_fields
         !procedure :: dump_planes 
   end type 

contains

subroutine dumpFullField(this, f, label, step)
   class(fof), intent(in) :: this
   character(len=4), intent(in) :: label
   character(len=clen) :: tempname, fname
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in ) :: f
   integer, intent(in) :: step 

   write(tempname,"(A3,I2.2,A7,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID,"_FILTER",this%fof_id, "_",label,"_t",step,".out"
   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
   call decomp_2d_write_one(1,f,fname,this%gpC)

end subroutine 

subroutine dumpYPlanes(this, f, label, pid, step)
   class(fof), intent(in) :: this
   character(len=4), intent(in) :: label
   character(len=clen) :: tempname, fname
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in ) :: f
   integer, dimension(:), intent(in) :: pid
   integer, intent(in) :: step 
   integer :: j 

   do j = 1,size(pid)
      write(tempname,"(A3,I2.2,A7,I2.2,A2,I4.4,A1,A4,A2,I6.6,A4)") "Run",this%runID,"_FILTER",this%fof_id,"_y",pid(j),"_",label,"_t",step,".pln"
      fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
      call decomp_2d_write_plane(1,f,2, pid(j), fname, this%gpC)
   end do 
end subroutine


subroutine init(this, runID, inputdir, outputdir, fof_id, spectC, cbuffyC, cbuffzC, PeriodicInZ, gpC)
   use reductions, only: p_maxval
   class(fof), intent(inout) :: this
   character(len=clen), intent(in) :: outputdir, inputdir
   integer, intent(in) :: fof_id, runID
   class(spectral), intent(in), target :: spectC
   complex(rkind), dimension(:,:,:), intent(in), target :: cbuffyC, cbuffzC
   logical, intent(in) :: PeriodicInZ
   type(decomp_info), intent(in), target :: gpC

   real(rkind) :: kx_max, ky_max, kz_max
   real(rkind) :: filterfact_x, filterfact_y, filterfact_z
   character(len=clen) :: fname, tempname
   integer :: ioUnit

   namelist /NLFOF/ filterfact_x, filterfact_y, filterfact_z 

   write(tempname,"(A10,I3.3,A10)") "FilterInfo", fof_id, "_input.inp"
   fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)

   ioUnit = 11
   open(unit=ioUnit, file=trim(fname), form='FORMATTED')
   read(unit=ioUnit, NML=NLFOF)
   close(ioUnit)


   this%spectC => spectC
   this%sp_gpC => spectC%spectdecomp

   this%outputdir = outputdir

   this%PeriodicInZ = PeriodicInZ

   this%fhaty => cbuffyC(:,:,:)
   this%fhatz => cbuffzC(:,:,:)

   this%fof_id = fof_id
   this%runID = runID 

   this%gpC => gpC

   allocate(this%Gfilt_x(size(this%fhatz,1)))
   allocate(this%Gfilt_y(size(this%fhatz,2)))
   allocate(this%Gfilt_z(size(this%fhatz,3)))

   kx_max = p_maxval(maxval(abs(this%spectC%k1)))
   ky_max = p_maxval(maxval(abs(this%spectC%k2)))
   kz_max = p_maxval(maxval(abs(this%spectC%k3)))

   this%Gfilt_x = 1.d0
   where(abs(this%spectC%k1inZ) > filterfact_x*kx_max)
      this%Gfilt_x = 0.d0
   end where

   this%Gfilt_y = 1.d0
   where(abs(this%spectC%k2inZ) > filterfact_y*ky_max)
      this%Gfilt_y = 0.d0
   end where
   
   this%Gfilt_z = 1.d0
   where(abs(this%spectC%k3inZ) > filterfact_z*kz_max)
      this%Gfilt_z = 0.d0
   end where
 
   if (.not. this%PeriodicInZ) then
      call gracefulExit("The problem must be periodic in z, if you use fof derived type", 1324)
   end if 

   call message(1,"Initialized on-the-fly filter number:", fof_id)
   call message(2,"Filter factor in x:", filterfact_x)
   call message(2,"Filter factor in y:", filterfact_y)
   call message(2,"Filter factor in z:", filterfact_z)
   call message(2,"Max allowable wave in x:", filterfact_x*kx_max)
   call message(2,"Max allowable wave in y:", filterfact_x*ky_max)
   call message(2,"Max allowable wave in z:", filterfact_x*kz_max)

end subroutine 

subroutine destroy(this)
   class(fof), intent(inout) :: this

   nullify(this%fhaty, this%fhatz)
   deallocate(this%Gfilt_x, this%Gfilt_y, this%Gfilt_z)
end subroutine 

subroutine filter_Complex2Real(this, fhat, ffilt)
   class(fof), intent(inout) :: this
   complex(rkind), dimension(size(this%fhaty,1),size(this%fhaty,2),size(this%fhaty,3)), intent(in) :: fhat
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(out) :: ffilt 
   integer :: i, j, k

   call transpose_y_to_z(fhat, this%fhatz, this%sp_gpC)
   call this%spectC%take_fft1d_z2z_ip(this%fhatz)
   do k = 1,size(this%fhatz,3)
      do j = 1,size(this%fhatz,2)
         !$omp simd
         do i = 1,size(this%fhatz,1)
            this%fhatz(i,j,k) = this%fhatz(i,j,k)*this%Gfilt_x(i)
            this%fhatz(i,j,k) = this%fhatz(i,j,k)*this%Gfilt_y(j)
            this%fhatz(i,j,k) = this%fhatz(i,j,k)*this%Gfilt_z(k)
         end do 
      end do 
   end do

   call this%spectC%take_ifft1d_z2z_ip(this%fhatz)
   call transpose_z_to_y(this%fhatz,this%fhaty,this%sp_gpC)
   call this%spectC%ifft(this%fhaty, ffilt)
end subroutine 

subroutine filter_Real2Real(this, f, ffilt)
   class(fof), intent(inout) :: this
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in ) :: f
   real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(out) :: ffilt 
   integer :: i, j, k

   call this%spectC%fft(f,this%fhaty)
   call transpose_y_to_z(this%fhaty, this%fhatz, this%sp_gpC)
   call this%spectC%take_fft1d_z2z_ip(this%fhatz)
   do k = 1,size(this%fhatz,3)
      do j = 1,size(this%fhatz,2)
         !$omp simd
         do i = 1,size(this%fhatz,1)
            this%fhatz(i,j,k) = this%fhatz(i,j,k)*this%Gfilt_x(i)
            this%fhatz(i,j,k) = this%fhatz(i,j,k)*this%Gfilt_y(j)
            this%fhatz(i,j,k) = this%fhatz(i,j,k)*this%Gfilt_z(k)
         end do 
      end do 
   end do

   call this%spectC%take_ifft1d_z2z_ip(this%fhatz)
   call transpose_z_to_y(this%fhatz,this%fhaty,this%sp_gpC)
   call this%spectC%ifft(this%fhaty, ffilt)
end subroutine 

end module
