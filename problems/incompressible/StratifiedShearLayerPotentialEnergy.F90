module igrid_Operators_NonPeriodicZ
   use kind_parameters, only: rkind, clen
   use spectralMod, only: spectral  
   use decomp_2d
   use cd06staggstuff, only: cd06stagg
   use reductions, only: p_sum
   use sorting_mod, only: sortgroup, Qsort    

   implicit none
   
   private
   public :: Ops_NP_Z

   integer :: niter = 4
   type :: Ops_NP_Z
      private
      complex(rkind), dimension(:,:,:), allocatable :: cbuffy1, cbuffy2
      real(rkind),    dimension(:,:,:), allocatable :: rbuffy, rbuffz1, rbuffz2
      type(decomp_info), pointer :: gp
      type(spectral)  :: spect
      type(cd06stagg) :: derZ
      real(rkind), dimension(:,:,:), allocatable :: zarr1d_1, zarr1d_2

      character(len=clen) ::  inputdir, outputdir
      real(rkind) :: mfact_xy, dx, dy, dz, mfact_xyz
      type(sortgroup), dimension(:), allocatable :: sort_z, sort_y
      
      real(rkind), dimension(:), allocatable :: zline
      integer :: RunID

      real(rkind), dimension(:), allocatable :: zold 
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
         procedure :: getPotentialEnergy
         procedure :: getAPE
         procedure :: sortPotT
   end type 

contains

subroutine sortPotT(this, T, Tsort) 
   class(Ops_NP_Z), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in) :: T
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: Tsort
   integer :: iiter, idx, i, j, k 

   call transpose_y_to_z(T, this%rbuffy, this%gp)
   do iiter = 1,niter
      call transpose_y_to_z(this%rbuffy,this%rbuffz1,this%gp)
      idx = 1
      do k = 1,this%gp%zsz(3)
         do j = 1,this%gp%zsz(2)
            do i = 1,this%gp%zsz(1)
               this%sort_z(idx)%value = this%rbuffz1(i,j,k)
               idx = idx + 1
            end do 
         end do 
      end do   
      
      call Qsort(this%sort_z, this%gp%zsz(1)*this%gp%zsz(2)*this%gp%zsz(3)) 
      
      idx = 1
      do k = 1,this%gp%zsz(3)
         do j = 1,this%gp%zsz(2)
            do i = 1,this%gp%zsz(1)
               this%rbuffz1(i,j,k) = this%sort_z(idx)%value
               idx = idx + 1
            end do 
         end do 
      end do
      
   
      call transpose_z_to_y(this%rbuffz1,this%rbuffy,this%gp)
      
      idx = 1
      do k = 1,this%gp%ysz(3)
         do j = 1,this%gp%ysz(2)
            do i = 1,this%gp%ysz(1)
               this%sort_y(idx)%value = this%rbuffy(i,j,k)
               idx = idx + 1
            end do 
         end do 
      end do 

      call Qsort(this%sort_y, this%gp%ysz(1)*this%gp%ysz(2)*this%gp%ysz(3)) 

      idx = 1
      do k = 1,this%gp%ysz(3)
         do j = 1,this%gp%ysz(2)
            do i = 1,this%gp%ysz(1)
               this%rbuffy(i,j,k) = this%sort_y(idx)%value
               idx = idx + 1
            end do 
         end do 
      end do 

   end do

   call transpose_y_to_x(this%rbuffy, Tsort, this%gp)
end subroutine 

function getCenterlineQuantity(this, vec) result(val)
   class(Ops_NP_Z), intent(in), target :: this
   real(rkind), dimension(this%gp%zsz(3)), intent(in) :: vec
   integer :: nz
   real(rkind) :: val

   nz = size(vec)

   val = 0.5d0*(vec(nz/2) + vec(nz/2+1))

end function

subroutine init(this, nx, ny, nz, dx, dy, dz, Lz, gp, InputDir, OutputDir, RunID)
   use gridtools, only: linspace
   class(Ops_NP_Z), intent(out), target :: this
   integer, intent(in) :: nx, ny, nz
   real(rkind), intent(in) :: dx, dy, dz, Lz
   class(decomp_info), intent(in), target :: gp
   character(len=clen), intent(in) ::  inputdir, outputdir
   integer :: RunID, idx, k

   this%gp => gp
   call this%spect%init("x",nx,ny,nz,dx, dy, dz, "four", "2/3rd", 2 , fixOddball=.false., &
                  exhaustiveFFT=.TRUE., init_periodicInZ=.FALSE., dealiasF=(2.d0/3.d0))
   
   call this%derZ%init(nz, dz, isTopEven = .false., isBotEven = .false., &
                                   isTopSided = .true., isBotSided = .true.)

   call this%spect%alloc_r2c_out(this%cbuffy1) 
   
   allocate(this%zline(nz))
   this%zline(1) = -Lz + dz/2.d0
   do idx = 2,nz
      this%zline(idx) = this%zline(idx-1) + dz
   end do 

   allocate(this%sort_z(gp%zsz(1)*gp%zsz(2)*gp%zsz(3))) 
   allocate(this%sort_y(gp%ysz(1)*gp%ysz(2)*gp%ysz(3))) 
   allocate(this%zold(gp%ysz(1)*gp%ysz(2)*gp%ysz(3))) 

   allocate(this%rbuffy(gp%ysz(1),gp%ysz(2),gp%ysz(3)))
   allocate(this%rbuffz1(gp%zsz(1),gp%zsz(2),gp%zsz(3)))
   allocate(this%rbuffz2(gp%zsz(1),gp%zsz(2),gp%zsz(3)))
  
   idx = 1
   do k = gp%yst(3),gp%yen(3)
      this%rbuffy(:,:,idx) = this%zline(k)
      idx = idx + 1
   end do 
   this%zold = reshape(this%rbuffy,[size(this%zold)])

   this%dx = dx
   this%dy = dy
   this%dz = dz
   this%inputdir  = inputdir
   this%outputdir = outputdir
   this%RunID = RunID

   this%mfact_xy = 1.d0/(real(nx,rkind)*real(ny,rkind))
   this%mfact_xyz = 1.d0/(real(nx,rkind)*real(ny,rkind)*real(nz,rkind))

   allocate(this%zarr1d_1(1,1,nz))
   allocate(this%zarr1d_2(1,1,nz))

end subroutine 

subroutine destroy(this)
   class(Ops_NP_Z), intent(inout), target :: this
  
   deallocate(this%rbuffy, this%rbuffz1, this%rbuffz2)
   deallocate(this%zarr1d_1, this%zarr1d_2, this%zold, this%sort_z, this%sort_y)
end subroutine 

subroutine ddx(this,f, dfdx)
   class(Ops_NP_Z), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: dfdx
  
   call this%spect%fft(f,this%cbuffy1)
   call this%spect%mtimes_ik1_ip(this%cbuffy1)
   call this%spect%ifft(this%cbuffy1,dfdx)
end subroutine 

subroutine ddy(this,f, dfdy)
   class(Ops_NP_Z), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: dfdy
  
   call this%spect%fft(f,this%cbuffy1)
   call this%spect%mtimes_ik2_ip(this%cbuffy1)
   call this%spect%ifft(this%cbuffy1,dfdy)
end subroutine 

subroutine ddz(this, f, dfdz)
   class(Ops_NP_Z), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: f
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: dfdz


   call transpose_x_to_y(f,this%rbuffy,this%gp)
   call transpose_y_to_z(this%rbuffy,this%rbuffz1,this%gp)
   call this%derZ%ddz_C2C(this%rbuffz1,this%rbuffz2,this%gp%zsz(1),this%gp%zsz(2))
   call transpose_z_to_y(this%rbuffz2,this%rbuffy,this%gp)
   call transpose_y_to_x(this%rbuffy,dfdz,this%gp)
end subroutine 

subroutine ddz_1d(this, f1d, dfdz1d)
   class(Ops_NP_Z), intent(inout) :: this
   real(rkind), dimension(this%gp%zsz(3)), intent(in)  :: f1d
   real(rkind), dimension(this%gp%zsz(3)), intent(out) :: dfdz1d

   this%zarr1d_1(1,1,:) = f1d
   call this%derZ%ddz_C2C(this%zarr1d_1,this%zarr1d_2,1,1)
   dfdz1d = this%zarr1d_2(1,1,:)
end subroutine 

subroutine getFluct_from_MeanZ(this, f, ffluct)
   class(Ops_NP_Z), intent(inout) :: this
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
   class(Ops_NP_Z), intent(inout) :: this
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
   class(Ops_NP_Z), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out)  :: field
   character(len=clen) :: tempname, fname
   character(len=4), intent(in) :: label
   integer, intent(in) :: tidx
         
   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",tidx,".out"
   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
   
   call decomp_2d_read_one(1,field,fname,this%gp)
end subroutine  

subroutine WriteField3D(this, field, label, tidx)
   use decomp_2d_io
   class(Ops_NP_Z), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: field
   character(len=clen) :: tempname, fname
   character(len=4), intent(in) :: label
   integer, intent(in) :: tidx
         
   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",tidx,".out"
   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
   
   call decomp_2d_write_one(1,field,fname,this%gp)
end subroutine  

subroutine allocate3Dfield(this, field)
   class(Ops_NP_Z), intent(inout) :: this
   real(rkind), dimension(:,:,:), allocatable, intent(out) :: field
   allocate(field(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)))
end subroutine 


function getSimTime(this, tidx) result(time)
   class(Ops_NP_Z), intent(in) :: this
   character(len=clen) :: tempname, fname
   integer, intent(in) :: tidx
   real(rkind) :: time

   write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_","info","_t",tidx,".out"
   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
   open(unit=10,file=fname,access='sequential',form='formatted')
   read(10,*) time
   close(10)

end function

subroutine getPotentialEnergy(this, T, TPE)
   class(Ops_NP_Z), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: T
   real(rkind), intent(out) :: TPE
   integer :: k 

   call transpose_x_to_y(T,this%rbuffy,this%gp)
   call transpose_y_to_z(this%rbuffy,this%rbuffz1,this%gp)
   do k = 1,this%gp%zsz(3)
      this%rbuffz1(:,:,k) = -this%zline(k)*this%rbuffz1(:,:,k)
   end do 

   TPE = p_sum(sum(this%rbuffz1))*this%dx*this%dy*this%dz

end subroutine 

subroutine getAPE(this, Tfluct, dTdz, APE)
   class(Ops_NP_Z), intent(inout) :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in)  :: Tfluct
   real(rkind), dimension(this%gp%zsz(3)), intent(in) :: dTdz
   real(rkind), intent(out) :: APE
   integer :: k 

   call transpose_x_to_y(Tfluct,this%rbuffy,this%gp)
   call transpose_y_to_z(this%rbuffy,this%rbuffz1,this%gp)

   this%rbuffz1 = this%rbuffz1*this%rbuffz1

   do k = 1,this%gp%zsz(3)
      this%rbuffz1(:,:,k) = dTdz(k)*this%rbuffz1(:,:,k)
   end do 

   APE = p_sum(sum(this%rbuffz1))*this%dx*this%dy*this%dz


end subroutine 

subroutine WriteASCII_2D(this, field, flabel)
   use basic_io, only: write_2d_ascii 
   class(Ops_NP_Z), intent(inout) :: this
   real(rkind), dimension(:,:), intent(in) :: field
   character(len=4), intent(in) :: flabel
   character(len=clen) :: tempname, fname
   
   write(tempname,"(A3,I2.2,A1,A4,A4)") "Run",this%runID, "_",flabel,".stt"
   fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
   
   call write_2d_ascii(field, fname)

end subroutine 


end module 

program StratifiedShearLayerPotentialEnergy
   use kind_parameters, only: rkind, clen
   use igrid_Operators_NonPeriodicZ, only: Ops_NP_Z
   use constants, only: pi, two
   use mpi
   use decomp_2d
   use timer, only: tic, toc
   use exits, only: message
   use decomp_2d_io
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: T, Tsort!, uFluct, vFluct, pFluct, Tfluct
   !real(rkind), dimension(:),     allocatable :: Prod, Transp_Conv, Transp_Press, Transp_Visc, Dissp, DisspSGS, DisspTheta, DisspThetaSGS, Buoy

   real(rkind), dimension(:), allocatable :: Tmn, dTdz
   real(rkind), parameter :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0
   real(rkind) :: dx, dy, dz, Re = 3000.d0, Rib = 0.2d0, Pr = 1.d0
   real(rkind) :: T0 = 100.d0
   integer :: nx, ny, nz, RunID, TIDX, np
   type(decomp_info) :: gp
   type(Ops_NP_Z) :: ops
   logical :: periodicbcs(3)
   integer :: ierr
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile

   real(rkind), dimension(10000) :: time, TPE, APE, BPE 

   real(rkind), dimension(:,:), allocatable :: data2write
   integer :: idx, tstart, tstop, tstep

   namelist /INPUT/ InputDir, OutputDir, RunID, tstart, tstop, T0, tstep, nx, ny, nz, Re, Rib, Pr
    
   call MPI_Init(ierr)               
   call MPI_Comm_size ( MPI_COMM_WORLD, np, ierr )
   call GETARG(1,inputfile)          
   open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=99, NML=INPUT)
   close(unit=99)
   periodicbcs(1) = .true.; periodicbcs(2) = .true.; periodicbcs(3) = .false.  
   
   call decomp_2d_init(nx, ny, nz, 1, np, periodicbcs)
   call get_decomp_info(gp)

   dx =     Lx/real(nx,rkind) 
   dy =     Ly/real(ny,rkind) 
   dz = two*Lz/real(nz,rkind)

   call ops%init(nx, ny, nz, dx, dy, dz, Lz, gp, InputDir, OutputDir, RunID)

   call ops%allocate3DField(T)
   call ops%allocate3DField(Tsort)


   allocate(Tmn(nz))
   allocate(dTdz(nz))


   tidx = tstart
   idx = 1
   do while(tidx <= tstop)
      call message(0, "Reading fields for tid:", TIDX)
      call tic()
      call ops%ReadField3D(T,"potT",TIDX)
      T = Rib*(T - T0)
      time(idx) = ops%getSimTime(tidx)
      call message(0, "Read simulation data at time:", time(idx))

      call ops%getPotentialEnergy(T, TPE(idx))
      call ops%sortPotT(T,Tsort)
      !print*, maxval(abs(Tsort - T))
      call ops%getPotentialEnergy(Tsort, BPE(idx))
      APE(idx) = TPE(idx) - BPE(idx)
  
      call message(1, "TPE:", TPE(idx))
      call message(1, "APE:", APE(idx))
      call message(1, "BPE:", BPE(idx))
      call toc()

      tidx = tidx + tstep
      idx = idx + 1
   end do 

   idx = idx - 1
   allocate(data2write(idx,4))
   data2write(:,1) = time(1:idx)
   data2write(:,2) = TPE(1:idx) 
   data2write(:,3) = APE(1:idx) 
   data2write(:,4) = BPE(1:idx)
   if (nrank == 0) then
      call ops%WriteASCII_2D(data2write, "TPEd")
   end if 

   call mpi_barrier(mpi_comm_world, ierr)
   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

