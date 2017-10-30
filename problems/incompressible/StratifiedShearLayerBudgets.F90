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
   end type 

contains

subroutine init(this, nx, ny, nz, dx, dy, dz, gp, InputDir, OutputDir, RunID)
   class(Ops_NP_Z), intent(out), target :: this
   integer, intent(in) :: nx, ny, nz
   real(rkind), intent(in) :: dx, dy, dz
   class(decomp_info), intent(in), target :: gp
   character(len=clen), intent(in) ::  inputdir, outputdir
   integer :: RunID

   this%gp => gp
   call this%spect%init("x",nx,ny,nz,dx, dy, dz, "four", "2/3rd", 2 , fixOddball=.false., &
                  exhaustiveFFT=.TRUE., init_periodicInZ=.FALSE., dealiasF=(2.d0/3.d0))
   
   call this%derZ%init(nz, dz, isTopEven = .false., isBotEven = .false., &
                                   isTopSided = .true., isBotSided = .true.)

   call this%spect%alloc_r2c_out(this%cbuffy1) 
   !call this%spect%alloc_r2c(this%cbuffy2) 
   
   this%inputdir  = inputdir
   this%outputdir = outputdir
   this%RunID = RunID

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
   call this%spect%mtimes_ik1_ip(this%cbuffy1)
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


end module 


program StratifiedShearLayerBudgets
   use kind_parameters, only: rkind, clen
   use igrid_Operators_NonPeriodicZ, only: Ops_NP_Z
   use constants, only: pi, two
   use mpi
   use decomp_2d
   use timer, only: tic, toc
   use exits, only: message
   implicit none

   real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, buff3, buff4, buff5
   real(rkind), dimension(:,:,:), allocatable :: u, v, w, p, T, nuSGS, kappaSGS!, uFluct, vFluct, pFluct, Tfluct
   !real(rkind), dimension(:),     allocatable :: Prod, Transp_Conv, Transp_Press, Transp_Visc, Dissp, DisspSGS, DisspTheta, DisspThetaSGS, Buoy

   real(rkind), parameter :: Lx = 9.d0*pi, Ly = 9.d0*pi, Lz = 8.d0
   real(rkind) :: dx, dy, dz, Re = 3000.d0, Rib = 0.05d0, Pr = 1.d0
   integer :: nx, ny, nz, RunID, TIDX
   type(decomp_info) :: gp
   type(Ops_NP_Z) :: ops
   logical :: periodicbcs(3)
   integer :: ierr
   character(len=clen) ::  inputdir, outputdir
   character(len=clen) :: inputfile

   namelist /INPUT/ InputDir, OutputDir, RunID, TIDX, nx, ny, nz, Re, Rib, Pr
    
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
   call ops%allocate3DField(v)
   call ops%allocate3DField(w)
   call ops%allocate3DField(p)
   call ops%allocate3DField(T)
   call ops%allocate3DField(nuSGS)
   call ops%allocate3DField(kappaSGS)
   call ops%allocate3DField(buff1)
   call ops%allocate3DField(buff2)
   call ops%allocate3DField(buff3)
   call ops%allocate3DField(buff4)
   call ops%allocate3DField(buff5)



   call message(0, "Reading fields..")
   call tic()
   call ops%ReadField3D(u,"uVel",TIDX)
   call ops%ReadField3D(v,"vVel",TIDX)
   call ops%ReadField3D(w,"wVel",TIDX)
   call ops%ReadField3D(T,"potT",TIDX)
   call ops%ReadField3D(P,"prss",TIDX)
   call ops%ReadField3D(nuSGS,"nSGS",TIDX)
   call ops%ReadField3D(kappaSGS,"kSGS",TIDX)
   call toc()


   ! STEP 1: Compute TKE
   call ops%getFluct_from_MeanZ(u,buff1)
   buff2 = buff1*buff1
   call ops%getFluct_from_MeanZ(v,buff1)
   buff2 = buff2 + buff1*buff1
   buff2 = buff2 + w*w        ! w mean is exactly 0
   buff2 = buff2*0.5d0
   call ops%WriteField3d(buff2,"itke",tidx)

   ! STEP 2: Compute the convective flux
   buff1 = w*buff2
   call ops%ddz(buff1,buff3)
   call ops%WriteField3d(buff1,"cflx",tidx)
   call ops%WriteField3d(buff3,"ctrn",tidx)

   ! STEP 3: Compute the pressure flux
   call ops%getFluct_from_MeanZ(p,buff1)
   buff1 = buff1*buff2
   call ops%ddz(buff1,buff3)
   call ops%WriteField3d(buff1,"pflx",tidx)
   call ops%WriteField3d(buff3,"ptrn",tidx)

   ! STEP 4: Compute the viscous
   call ops%ddz(buff2,buff1)
   buff1 = (1.d0/Re)*buff1
   call ops%ddz(buff1,buff3)
   call ops%WriteField3d(buff1,"vflx",tidx)
   call ops%WriteField3d(buff3,"vtrn",tidx)
  
   ! STEP 5: Compute the production
   call ops%getFluct_from_MeanZ(u,buff1)
   buff3 = u - buff1 ! mean
   call ops%ddz(buff3,buff2) ! dUdz
   buff3 = -buff1*w      ! u'w' (since wmean = 0)
   buff2 = buff2*buff3 
   call ops%WriteField3d(buff2,"prod",tidx)
  
   ! STEP 6: Compute the Buoyancy term
   call ops%getFluct_from_MeanZ(T,buff1)
   buff2 = Rib*buff1*w
   call ops%WriteField3d(buff2,"buoy",tidx)


   ! STEP 7: Compute dissipation rate
   buff3 = 0.d0
   call ops%getFluct_from_MeanZ(u,buff1)
   call ops%ddx(buff1,buff2)
   buff3 = buff3 + buff2*buff2
   call ops%ddy(buff1,buff2)
   buff3 = buff3 + buff2*buff2
   call ops%ddz(buff1,buff2)
   buff3 = buff3 + buff2*buff2

   call ops%getFluct_from_MeanZ(v,buff1)
   call ops%ddx(buff1,buff2)
   buff3 = buff3 + buff2*buff2
   call ops%ddy(buff1,buff2)
   buff3 = buff3 + buff2*buff2
   call ops%ddz(buff1,buff2)
   buff3 = buff3 + buff2*buff2

   buff1 = w
   call ops%ddx(buff1,buff2)
   buff3 = buff3 + buff2*buff2
   call ops%ddy(buff1,buff2)
   buff3 = buff3 + buff2*buff2
   call ops%ddz(buff1,buff2)
   buff3 = buff3 + buff2*buff2
   
   buff3 = (1.d0/Re)*buff3
   call ops%WriteField3d(buff3,"diss",tidx)

   ! STEP 8: SGS sink term
   call ops%getFluct_from_MeanZ(u,buff1)
   call ops%getFluct_from_MeanZ(v,buff2)
   
   ! s11*s11
   call ops%ddx(buff1,buff3)
   buff5 = buff3*buff3
   
   ! 2*s12*s12
   call ops%ddy(buff1,buff3)
   call ops%ddx(buff2,buff4)
   buff3 = 0.5d0*(buff3 + buff4)
   buff5 = buff5 + 2.d0*buff3*buff3 

   ! 2*s13*s13 
   call ops%ddz(buff1,buff3)
   call ops%ddx(w    ,buff4)
   buff3 = 0.5d0*(buff3 + buff4)
   buff5 = buff5 + 2.d0*buff3*buff3 

   ! 2*s23*s23 
   call ops%ddz(buff2,buff3)
   call ops%ddy(w    ,buff4)
   buff3 = 0.5d0*(buff3 + buff4)
   buff5 = buff5 + 2.d0*buff3*buff3 

   ! s22*s22
   call ops%ddy(buff2,buff3)
   buff5 = buff5 + buff3*buff3
   
   ! s33*s33
   call ops%ddz(w,buff3)
   buff5 = buff5 + buff3*buff3
   buff5 = 2.d0*nuSGS*buff5
   call ops%WriteField3d(buff5,"sgsD",tidx)

   ! STEP 9: Mean SGS term
   call ops%getFluct_from_MeanZ(u,buff1)
   buff2 = u - buff1 ! Mean u
   call ops%ddz(buff2,buff3)  ! dUdz
   
   call ops%ddz(buff1,buff2)
   call ops%ddx(w    ,buff4)
   buff4 = 0.5d0*(buff2 + buff4)
   buff4 = 2.d0*nuSGS*buff4*buff3
   call ops%WriteField3d(buff4,"SGSd",tidx)

   ! STEP 10: SGS transport 
   call ops%getFluct_from_MeanZ(u,buff1)
   call ops%getFluct_from_MeanZ(v,buff2)
   
   call ops%ddz(u,buff3)
   call ops%ddx(w,buff4)
   buff3 = 0.5d0*(buff3 + buff4)
   buff5 = buff1*buff3

   call ops%ddz(v,buff3)
   call ops%ddy(w,buff4)
   buff3 = 0.5d0*(buff3 + buff4)
   buff5 = buff5 + buff2*buff3

   call ops%ddz(w,buff3)
   buff5 = buff5  + w*buff3
   
   buff5 = -2.d0*nuSGS
   call ops%ddz(buff5,buff3)
   call ops%WriteField3d(buff3,"Tsgs",tidx)


   deallocate(u, v, w, p, T, nuSGS, kappaSGS)
   call ops%destroy()
   call MPI_Finalize(ierr)           


end program 

