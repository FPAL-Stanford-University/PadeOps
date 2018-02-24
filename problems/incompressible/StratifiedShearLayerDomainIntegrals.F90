module igrid_Operators_NonPeriodicZ
   use kind_parameters, only: rkind, clen
   use spectralMod, only: spectral  
   use decomp_2d
   use cd06staggstuff, only: cd06stagg
   use reductions, only: p_sum

   implicit none
   
   private
   public :: Ops_NP_Z

   type :: Ops_NP_Z
      private
      complex(rkind), dimension(:,:,:), allocatable :: cbuffy1, cbuffy2
      real(rkind),    dimension(:,:,:), allocatable :: rbuffy, rbuffz1, rbuffz2
      type(decomp_info), pointer :: gp
      type(spectral)  :: spect
      type(cd06stagg) :: derZ
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
   end type 

contains

function getVolumeIntegral(this, arr) result(val)
   class(Ops_NP_Z), intent(in), target :: this
   real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in) :: arr
   real(rkind) :: val 

   val = this%dx*this%dy*this%dz*p_sum(arr)

end function 

function getCenterlineQuantity(this, vec) result(val)
   class(Ops_NP_Z), intent(in), target :: this
   real(rkind), dimension(this%gp%zsz(3)), intent(in) :: vec
   integer :: nz
   real(rkind) :: val

   nz = size(vec)

   val = 0.5d0*(vec(nz/2) + vec(nz/2+1))

end function

subroutine init(this, nx, ny, nz, dx, dy, dz, gp, InputDir, OutputDir, RunID)
   class(Ops_NP_Z), intent(out), target :: this
   integer, intent(in) :: nx, ny, nz
   real(rkind), intent(in) :: dx, dy, dz
   class(decomp_info), intent(in), target :: gp
   character(len=clen), intent(in) ::  inputdir, outputdir
   integer :: RunID

   this%gp => gp
   call this%spect%init("x",nx,ny,nz,dx, dy, dz, "FOUR", "2/3rd", 2 , fixOddball=.false., &
                  exhaustiveFFT=.TRUE., init_periodicInZ=.FALSE., dealiasF=(2.d0/3.d0))
  
   call this%derZ%init(nz, dz, isTopEven = .false., isBotEven = .false., &
                                   isTopSided = .true., isBotSided = .true.)

   call this%spect%alloc_r2c_out(this%cbuffy1) 
   !call this%spect%alloc_r2c(this%cbuffy2) 
   
   this%inputdir  = inputdir
   this%outputdir = outputdir
   this%RunID = RunID

   this%dx = dx
   this%dy = dy
   this%dz = dz

   this%mfact_xy = 1.d0/(real(nx,rkind)*real(ny,rkind))
   allocate(this%rbuffy(gp%ysz(1),gp%ysz(2),gp%ysz(3)))
   allocate(this%rbuffz1(gp%zsz(1),gp%zsz(2),gp%zsz(3)))
   allocate(this%rbuffz2(gp%zsz(1),gp%zsz(2),gp%zsz(3)))

   allocate(this%zarr1d_1(1,1,nz))
   allocate(this%zarr1d_2(1,1,nz))
end subroutine 

subroutine destroy(this)
   class(Ops_NP_Z), intent(inout), target :: this
  
   deallocate(this%rbuffy, this%rbuffz1, this%rbuffz2)
   deallocate(this%zarr1d_1, this%zarr1d_2)
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

program StratifiedShearLayerCenterline
   use kind_parameters, only: rkind, clen
   use igrid_Operators_NonPeriodicZ, only: Ops_NP_Z
   use constants, only: pi, two
   use mpi
   use decomp_2d
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: buff2, buff3, buff4
   real(rkind), dimension(:,:,:), allocatable :: u, v, w, ufluct, vfluct, T, nuSGS, Tfluct!, uFluct, vFluct, pFluct, Tfluct
   !real(rkind), dimension(:),     allocatable :: Prod, Transp_Conv, Transp_Press, Transp_Visc, Dissp, DisspSGS, DisspTheta, DisspThetaSGS, Buoy

   real(rkind) :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0
   real(rkind) :: dx, dy, dz, Re = 3000.d0, Rib = 0.05d0, Pr = 1.d0, Tref = 100.d0 
   integer :: nx, ny, nz, RunID, TIDX
   type(decomp_info) :: gp
   type(Ops_NP_Z) :: ops
   logical :: periodicbcs(3)
   integer :: ierr
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile

   real(rkind), dimension(10000) :: P, B, D, Dv, Dsgs, time, IEL, MKE, TKE 

   real(rkind), dimension(:,:), allocatable :: data2write
   integer :: idx, tstart, tstop, tstep

   namelist /INPUT/ Lx, Ly, Lz, InputDir, OutputDir, RunID, tstart, tstop, tstep, nx, ny, nz, Re, Rib, Pr, Tref 
    
   call MPI_Init(ierr)               
   call GETARG(1,inputfile)          
   open(unit=99, file=trim(inputfile), form='FORMATTED', iostat=ierr)
   read(unit=99, NML=INPUT)
   close(unit=99)
   periodicbcs(1) = .true.; periodicbcs(2) = .true.; periodicbcs(3) = .false.  
   call decomp_2d_init(nx, ny, nz, 0, 0, periodicbcs)
   call get_decomp_info(gp)

   dx =     Lx/real(nx,rkind) 
   dy =     Ly/real(ny,rkind) 
   dz = two*Lz/real(nz,rkind)

   call ops%init(nx, ny, nz, dx, dy, dz, gp, InputDir, OutputDir, RunID)

   call ops%allocate3DField(u)
   call ops%allocate3DField(ufluct)
   call ops%allocate3DField(v)
   call ops%allocate3DField(vfluct)
   call ops%allocate3DField(w)
   call ops%allocate3DField(T)
   call ops%allocate3DField(Tfluct)
   call ops%allocate3DField(nuSGS)
   call ops%allocate3DField(buff2)
   call ops%allocate3DField(buff3)
   call ops%allocate3DField(buff4)

   tidx = tstart
   idx = 1
   do while(tidx <= tstop)
      call message(0, "Reading fields for tid:", TIDX)
      call tic()
      call ops%ReadField3D(u,"uVel",TIDX)
      call ops%ReadField3D(v,"vVel",TIDX)
      call ops%ReadField3D(w,"wVel",TIDX)
      call ops%ReadField3D(T,"potT",TIDX)
      call ops%ReadField3D(nuSGS,"nSGS",TIDX)

      T = Rib*(T - Tref)

      time(idx) = ops%getSimTime(tidx)
      call message(0, "Read simulation data at time:", time(idx))

      ! STEP 1: Compute the gain in PE from IE
      call ops%ddz(T,buff2)
      buff2 = (1.d0/(Pr*Re))*buff2
      IEL(idx) = ops%getVolumeIntegral(buff2)
      
  
      ! STEP 4: Compute the production
      call ops%getFluct_from_MeanZ(u,ufluct)
      buff3 = u - ufluct ! mean
      call ops%ddz(buff3,buff2) ! dUdz
      buff3 = 0.5d0*buff3*buff3
      MKE(idx) = ops%getVolumeIntegral(buff3)
      buff3 = -ufluct*w      ! u'w' (since wmean = 0)
      buff2 = buff2*buff3 
      P(idx) = ops%getVolumeIntegral(buff2)
      

      call ops%getFluct_from_MeanZ(v,vfluct)
      buff3 = v - vfluct ! mean
      call ops%ddz(buff3,buff2) ! dVdz
      buff3 = 0.5d0*buff3*buff3
      MKE(idx) = MKE(idx) + ops%getVolumeIntegral(buff3)
      buff3 = -vfluct*w      ! u'w' (since wmean = 0)
      buff2 = buff2*buff3 
      P(idx) = P(idx) + ops%getVolumeIntegral(buff2)

      ! STEP 5a: Compute TKE
      buff2 = (ufluct*ufluct + vfluct*vfluct + w*w)
      TKE(idx) = 0.5d0*ops%getVolumeIntegral(buff2)

      ! STEP 5: Compute the Buoyancy term
      call ops%getFluct_from_MeanZ(T,buff2)
      buff2 = buff2*w
      B(idx) = ops%getVolumeIntegral(buff2)


      ! STEP 6: Compute dissipation rate
      call ops%ddx(ufluct,buff2)
      buff3 = buff2*buff2
      call ops%ddy(ufluct,buff2)
      buff3 = buff3 + buff2*buff2
      call ops%ddz(ufluct,buff2)
      buff3 = buff3 + buff2*buff2

      call ops%ddx(vfluct,buff2)
      buff3 = buff3 + buff2*buff2
      call ops%ddy(vfluct,buff2)
      buff3 = buff3 + buff2*buff2
      call ops%ddz(vfluct,buff2)
      buff3 = buff3 + buff2*buff2

      call ops%ddx(w,buff2)
      buff3 = buff3 + buff2*buff2
      call ops%ddy(w,buff2)
      buff3 = buff3 + buff2*buff2
      call ops%ddz(w,buff2)
      buff3 = buff3 + buff2*buff2
      
      buff3 = (1.d0/Re)*buff3
      Dv(idx) = ops%getVolumeIntegral(buff3)

      ! STEP 7: SGS sink term
      ! s11*s11
      call ops%ddx(ufluct,buff3)
      buff2 = buff3*buff3
      
      ! 2*s12*s12
      call ops%ddy(ufluct,buff3)
      call ops%ddx(vfluct,buff4)
      buff3 = 0.5d0*(buff3 + buff4)
      buff2 = buff2 + 2.d0*buff3*buff3 

      ! 2*s13*s13 
      call ops%ddz(ufluct,buff3)
      call ops%ddx(w     ,buff4)
      buff3 = 0.5d0*(buff3 + buff4)
      buff2 = buff2 + 2.d0*buff3*buff3 

      ! 2*s23*s23 
      call ops%ddz(vfluct,buff3)
      call ops%ddy(w     ,buff4)
      buff3 = 0.5d0*(buff3 + buff4)
      buff2 = buff2 + 2.d0*buff3*buff3 

      ! s22*s22
      call ops%ddy(vfluct,buff3)
      buff2 = buff2 + buff3*buff3
      
      ! s33*s33
      call ops%ddz(w,buff3)
      buff2 = buff2 + buff3*buff3
      
      buff2 = 2.d0*nuSGS*buff2
      Dsgs(idx) = ops%getVolumeIntegral(buff2)

      D(idx) = Dsgs(idx) + Dv(idx)
      call toc()

      if (nrank == 0) then
         print*, time(idx), P(idx), B(idx), D(idx), Dv(idx), Dsgs(idx), IEL(idx)
         print*, "ddt_TKE:", P(idx) + B(idx) - D(idx)
      end if 
     
      tidx = tidx + tstep
      idx = idx + 1
   end do 

   idx = idx - 1
   allocate(data2write(idx,9))
   data2write(:,1) = time(1:idx)
   data2write(:,2) = IEL(1:idx) 
   data2write(:,3) = P(1:idx) 
   data2write(:,4) = B(1:idx) 
   data2write(:,5) = D(1:idx) 
   data2write(:,6) = Dv(1:idx) 
   data2write(:,7) = Dsgs(1:idx)
   data2write(:,8) = MKE(1:idx)
   data2write(:,9) = TKE(1:idx)

   if (nrank == 0) then
      call ops%WriteASCII_2D(data2write, "avgV")
   end if 

   call mpi_barrier(mpi_comm_world, ierr)
   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

