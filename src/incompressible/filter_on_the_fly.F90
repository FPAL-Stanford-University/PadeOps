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
      class(decomp_info), pointer :: gpC, sp_gpC
      integer :: fof_id, runID
      class(spectral), pointer :: spectC 
      real(rkind), dimension(:), allocatable :: Gfilt_x, Gfilt_y, Gfilt_z
      complex(rkind), dimension(:,:,:), pointer :: fhaty, fhatz
     
      integer :: sz_x, sz_y, sz_z, filterType
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
   use constants, only: pi
   use reductions, only: p_maxval
   class(fof), intent(inout), target :: this
   character(len=clen), intent(in) :: outputdir, inputdir
   integer, intent(in) :: fof_id, runID
   class(spectral), intent(in), target :: spectC
   complex(rkind), dimension(:,:,:), intent(in), target :: cbuffyC, cbuffzC
   logical, intent(in) :: PeriodicInZ
   type(decomp_info), intent(in), target :: gpC

   real(rkind) :: kx_max, ky_max, kz_max, dx, dy, dz, fact_x, fact_y, fact_z
   real(rkind) :: filterfact_x, filterfact_y, filterfact_z
   character(len=clen) :: fname, tempname
   integer :: filterType = 1
   integer :: ioUnit

   !real(rkind) :: a0 = 0.5d0, a1 = 0.8056734d0, a2 = 0.4684092d0, a3 = 0.1627358d0
   real(rkind) :: a0 = 0.5d0, a1 = 0.75d0, a2 = 0.3d0, a3 = 0.05d0
   !real(rkind) :: alpha1 = 0.d0, alpha2 = 0.9368185d0
   real(rkind) :: alpha1 = 0.d0, alpha2 = 0.6d0
   real(rkind), dimension(:), pointer :: k

   namelist /FOF/ filterType, filterfact_x, filterfact_y, filterfact_z 

   write(tempname,"(A10,I3.3,A10)") "FilterInfo", fof_id, "_input.inp"
   fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)

   ioUnit = 11
   open(unit=ioUnit, file=trim(fname), form='FORMATTED')
   read(unit=ioUnit, NML=FOF)
   close(ioUnit)


   this%spectC => spectC
   this%sp_gpC => spectC%spectdecomp

   this%outputdir = outputdir

   this%PeriodicInZ = PeriodicInZ

   this%filtertype = filtertype
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

   dx = this%spectC%dx
   dy = this%spectC%dy
   dz = this%spectC%dz
   
   select case(this%filtertype)
   case(1)
      this%Gfilt_x = 1.d0
      where(abs(this%spectC%k1inZ) > filterfact_x*kx_max)
         this%Gfilt_x = 0.d0
      end where
   case(2)
      this%Gfilt_x = exp(-((this%spectC%k1inZ*dx/filterfact_x)**2.d0)/4.d0)
   case(3)
      k => this%spectC%k1inZ
      fact_x = 1.d0/filterfact_x
      this%Gfilt_x = a0*cos(0.d0) + a1*cos(fact_x*k*dx) + a2*cos(2*fact_x*k*dx) + a3*cos(3*fact_x*k*dx)
      this%Gfilt_x = this%Gfilt_x/(1 + alpha1*cos(fact_x*k*dx) + alpha2*cos(2*fact_x*k*dx))
      where(abs(fact_x*k*dx)>pi) this%Gfilt_x = 0.d0
   case default 
      this%Gfilt_x = 1.d0
   end select

   select case(this%filtertype)
   case(1)
      this%Gfilt_y = 1.d0
      where(abs(this%spectC%k2inZ) > filterfact_y*ky_max)
         this%Gfilt_y = 0.d0
      end where
   case(2)
      this%Gfilt_y = exp(-((this%spectC%k2inZ*dy/filterfact_y)**2.d0)/4.d0)
   case(3)
      k => this%spectC%k2inZ
      fact_y = 1.d0/filterfact_y
      this%Gfilt_y = a0*cos(0.d0) + a1*cos(fact_y*k*dy) + a2*cos(2*fact_y*k*dy) + a3*cos(3*fact_y*k*dy)
      this%Gfilt_y = this%Gfilt_y/(1 + alpha1*cos(fact_y*k*dy) + alpha2*cos(2*fact_y*k*dy))
      where(abs(fact_y*k*dy)>pi) this%Gfilt_y = 0.d0
   case default 
      this%Gfilt_y = 1.d0
   end select
   
   select case(this%filtertype)
      case(1)
      this%Gfilt_z = 1.d0
      where(abs(this%spectC%k3inZ) > filterfact_z*kz_max)
         this%Gfilt_z = 0.d0
      end where
   case(2)
      this%Gfilt_z = exp(-((this%spectC%k3inZ*dz/filterfact_z)**2.d0)/4.d0)
   case(3)
      k => this%spectC%k3inZ
      fact_z = 1.d0/filterfact_z
      this%Gfilt_z = a0*cos(0.d0) + a1*cos(fact_z*k*dz) + a2*cos(2*fact_z*k*dz) + a3*cos(3*fact_z*k*dz)
      this%Gfilt_z = this%Gfilt_z/(1 + alpha1*cos(fact_z*k*dz) + alpha2*cos(2*fact_z*k*dz))
      where(abs(fact_z*k*dz)>pi) this%Gfilt_z = 0.d0
   case default 
      this%Gfilt_z = 1.d0
   end select
 
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
