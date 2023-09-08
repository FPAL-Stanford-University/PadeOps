module scalar_igridMod
   use kind_parameters, only: rkind, clen
   use constants, only: imi, pi, zero,one,two,three,half, four,eight, nine, six, kappa 
   use decomp_2d
   use PadeDerOps, only: Pade6stagg
   use sgsmod_igrid, only: sgs_igrid
   use spectralMod, only: spectral
   use igrid_hooks!, only: setDirichletBC_Temp, set_Reference_Temperature, meshgen_WallM, initfields_wallM, set_planes_io, set_KS_planes_io 
   use fringeMethod, only: fringe 
   use io_hdf5_stuff, only: io_hdf5 
   implicit none
   
   private
   public :: scalar_igrid
   
   
   
   type :: scalar_igrid
      private 
      type(decomp_info), pointer :: gpC, gpE
      type(spectral), pointer :: spectC, spectE
      type(decomp_info), pointer :: sp_gpC, sp_gpE
      type(Pade6stagg), pointer  :: der
      
      real(rkind), dimension(:,:,:), allocatable, public :: F
      complex(rkind), dimension(:,:,:), pointer, public :: Fhat

      real(rkind), dimension(:,:,:), allocatable :: d2Fdz2, dFdxC, dFdyC, dFdzC
      real(rkind), dimension(:,:,:), allocatable :: dfdzE
      
      complex(rkind), dimension(:,:,:), allocatable, public :: source_hat

      real(rkind), dimension(:,:,:,:), pointer :: rbuffxC, rbuffyC, rbuffzC
      real(rkind), dimension(:,:,:,:), pointer :: rbuffxE, rbuffyE, rbuffzE

      complex(rkind), dimension(:,:,:,:), pointer :: cbuffyC, cbuffzC
      complex(rkind), dimension(:,:,:,:), pointer :: cbuffyE, cbuffzE
      complex(rkind), dimension(:,:,:), pointer, public :: rhs
      complex(rkind), dimension(:,:,:,:), allocatable, public :: rhs_storage
      
      type(fringe), pointer :: fringe_x, fringe_x1, fringe_x2

      complex(rkind), dimension(:,:,:), pointer, public :: Fhat1, Fhat2, Fhat3, Fhat4
      complex(rkind), dimension(:,:,:,:), allocatable, public   :: Sfields

      real(rkind) :: Re, PrandtlNum, TurbPrandtlNum, Cy 
      character(len=clen) :: inputDataDir, outputDataDir
      real(rkind), dimension(:,:,:), pointer :: u, v, w, wC
      real(rkind), dimension(:,:,:,:), pointer :: duidxj

      integer :: RunID, scalar_number
      integer, public :: bc_bottom, bc_top
      type(sgs_igrid), pointer :: sgsmodel
      logical :: useSource, isinviscid, useSGS, usefringe, usedoublefringe

      real(rkind) :: lowbound, highbound 
      contains
         procedure :: init
         procedure :: destroy
         procedure :: dealias
         procedure :: populateRHS
         procedure :: prep_scalar
         procedure :: reset_pointers
         procedure, private :: addAdvectionTerm
         procedure, private :: addViscousTerm
         procedure, private :: get_derivatives
         procedure :: dumpRestart
         procedure :: readRestart
         procedure :: dump_planes
         procedure :: dumpScalarField
   end type 
contains


subroutine dealias(this)
   class(scalar_igrid), intent(inout) :: this
   call this%spectC%dealias(this%Fhat)
end subroutine 

subroutine prep_scalar(this)
   class(scalar_igrid), intent(inout) :: this

   call this%spectC%ifft(this%Fhat,this%F)
   call this%get_derivatives()

end subroutine 

subroutine reset_pointers(this,reset_rhs)
   class(scalar_igrid), intent(inout), target :: this
   logical, intent(in), optional :: reset_rhs
   this%Fhat  => this%Sfields(:,:,:,1)
   if (present(reset_rhs)) then
      if (reset_rhs) then  
         this%rhs => this%rhs_storage(:,:,:,1) 
      end if
   end if

end subroutine 

subroutine addAdvectionTerm(this)
   class(scalar_igrid), intent(inout) :: this

   this%rbuffxC(:,:,:,1) = -this%dFdxC*this%u
   call this%spectC%fft(this%rbuffxC(:,:,:,1), this%rhs)

   this%rbuffxC(:,:,:,1) = -this%dFdyC*this%v
   call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))

   this%rbuffxC(:,:,:,1) = -this%dFdzC*this%wC
   call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,2))

   this%rhs = this%rhs + this%cbuffyC(:,:,:,1) + this%cbuffyC(:,:,:,2)


end subroutine 

subroutine addViscousTerm(this)
   class(scalar_igrid), intent(inout) :: this
   integer :: i, j, k
   complex(rkind) :: tmp
   real(rkind) :: molecularDiff

   call this%spectC%fft(this%d2Fdz2, this%cbuffyC(:,:,:,1))
  
   molecularDiff = one/(this%Re*this%PrandtlNum)
   do k = 1,size(this%rhs,3)
      do j = 1,size(this%rhs,2)
         !$omp simd
         do i = 1,size(this%rhs,1)
            tmp = -this%spectC%kabs_sq(i,j,k)*this%Fhat(i,j,k) + this%cbuffyC(i,j,k,1)
            this%rhs(i,j,k) = this%rhs(i,j,k) + molecularDiff*tmp
         end do 
      end do 
   end do 

end subroutine

subroutine populateRHS(this, dt)
   class(scalar_igrid), intent(inout) :: this
   real(rkind), intent(in) :: dt

   call this%addAdvectionTerm()

   if (this%useSource) this%rhs = this%rhs + this%source_hat

   if (this%useSGS) then
      call this%sgsmodel%getRHS_SGS_Scalar(this%rhs, this%dFdxC, this%dFdyC, this%dFdzC, this%dFdzE, &
         this%u, this%v, this%wC, this%F, this%Fhat, this%duidxj, this%TurbPrandtlNum, this%Cy, this%lowbound, this%highbound)
   end if

   if (.not. this%isinviscid) then
      call this%addViscousTerm()
   end if 
   
   ! Add the fringe contributions
   if (this%usedoublefringe) then
      call this%fringe_x1%addFringeRHS_scalar(dt, this%rhs, this%F)
      call this%fringe_x2%addFringeRHS_scalar(dt, this%rhs, this%F)
   else
      if (this%useFringe) then
         call this%fringe_x%addFringeRHS_scalar(dt, this%rhs, this%F)
      end if
   end if 

end subroutine 


subroutine destroy(this)
   class(scalar_igrid), intent(inout) :: this

   deallocate(this%F, this%Fhat, this%dFdxC, this%dFdyC, this%dFdzC, this%dFdzE, this%rhs)
   deallocate(this%fhat1, this%fhat2, this%fhat3)

   if (allocated(this%source_hat)) deallocate(this%source_hat)
   if (allocated(this%d2Fdz2)) deallocate(this%d2Fdz2)

end subroutine


subroutine init(this,gpC,gpE,spectC,spectE,sgsmodel,der,inputFile, inputDir,mesh,u,v,w,wC,duidxj,    &
                  & rbuffxC,rbuffyC,rbuffzC,rbuffxE,rbuffyE,rbuffzE, cbuffyC,&
                  & cbuffzC,cbuffyE,cbuffzE, Re, isinviscid, useSGS, scalar_number, &
                  & InputDataDir, OutputDataDir, RunID, restartSim, tid_restart, &
                  & usefringe, usedoublefringe, fringe_x, fringe_x1, fringe_x2)
         
   use kind_parameters, only: clen

   class(scalar_igrid), intent(inout), target :: this
   type(decomp_info), target, intent(in) :: gpC, gpE
   type(spectral), target, intent(in) :: spectC, spectE
   type(sgs_igrid), target, intent(in) :: sgsmodel
   type(Pade6stagg), target, intent(in) :: der
   real(rkind), dimension(:,:,:), target, intent(in) :: u, v, w, wC
   real(rkind), dimension(:,:,:,:), target, intent(in) :: rbuffxC, rbuffyC, rbuffzC
   real(rkind), dimension(:,:,:,:), target, intent(in) :: rbuffxE, rbuffyE, rbuffzE
   real(rkind), dimension(:,:,:,:), intent(in) :: mesh
   real(rkind), dimension(:,:,:,:), target, intent(in) :: duidxj 
   complex(rkind), dimension(:,:,:,:), target, intent(in) ::  cbuffyC, cbuffzC
   complex(rkind), dimension(:,:,:,:), target, intent(in) ::  cbuffyE, cbuffzE
   logical, intent(in) :: isinviscid, useSGS, restartSim, usefringe, usedoublefringe
   type(fringe), target, intent(in) :: fringe_x, fringe_x1, fringe_x2
   real(rkind), intent(in) :: Re
   character(len=*), intent(in) :: inputFile, inputDir, InputDataDir, OutputDataDir
   integer, intent(in) :: scalar_number, RunID, tid_restart 
   integer :: ierr
   logical :: useSource = .false., RejectScalarRestart = .false. 
   real(rkind) :: PrandtlNum = 1.d0, TurbPrandtlNum = 1.d0, Cy = 100.d0 
   integer ::  bc_bottom = 1, bc_top = 1 
   character(len=clen) :: tempname, fname
   real(rkind) :: lowbound = 0.d0, highbound = 1.d0


   namelist /SCALAR_INFO/ useSource, PrandtlNum, bc_bottom, bc_top,TurbPrandtlNum, Cy, RejectScalarRestart, lowbound, highbound

   
   this%InputDataDir = InputDataDir
   this%OutputDataDir = OutputDataDir
   
    write(tempname,"(A7,I2.2,A10)") "SCALAR_", scalar_number, "_input.inp"
    fname = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
   
   open(unit=123, file=trim(fname), form='FORMATTED', iostat=ierr)
   read(unit=123, NML=SCALAR_INFO)
   close(123)

   this%PrandtlNum = PrandtlNum
   this%TurbPrandtlNum = TurbPrandtlNum
   this%Re = Re

   ! Scalar bounding
   this%Cy = Cy 
   this%lowbound = lowbound
   this%highbound = highbound

   this%der => der

   this%isinviscid = isinviscid
   this%useSGS = useSGS

   this%usefringe = usefringe
   this%usedoublefringe = usedoublefringe


   this%fringe_x => fringe_x
   this%fringe_x1 => fringe_x1
   this%fringe_x2 => fringe_x2

   this%RunID = RunID
   this%spectC => spectC
   this%spectE => spectE

   this%bc_bottom = bc_bottom
   this%bc_top = bc_top

   this%gpC    => gpC
   this%gpE    => gpE
   this%sp_gpC => spectC%spectdecomp
   this%sp_gpE => spectE%spectdecomp

   this%useSource = useSource

   this%scalar_number = scalar_number

   this%rbuffxC => rbuffxC
   this%rbuffyC => rbuffyC
   this%rbuffzC => rbuffzC

   this%rbuffxE => rbuffxE
   this%rbuffyE => rbuffyE
   this%rbuffzE => rbuffzE

   this%cbuffyC => cbuffyC
   this%cbuffzC => cbuffzC

   this%cbuffyE => cbuffyE
   this%cbuffzE => cbuffzE

   this%sgsmodel => sgsmodel

   this%u => u
   this%v => v
   this%w => w
   this%wC => wC
   this%duidxj => duidxj


   allocate(this%F (gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))

   allocate(this%dFdxC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
   allocate(this%dFdyC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
   allocate(this%dFdzC(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
   allocate(this%dFdzE(gpE%xsz(1),gpE%xsz(2),gpE%xsz(3)))
  
   allocate(this%Sfields(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3),4))
   allocate(this%rhs_storage(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3),2))

   this%Fhat  => this%Sfields(:,:,:,1)
   this%Fhat1 => this%Sfields(:,:,:,2)
   this%Fhat2 => this%Sfields(:,:,:,3)
   this%Fhat3 => this%Sfields(:,:,:,4)
   this%rhs => this%rhs_storage(:,:,:,1)


   if (.not. this%isInviscid) then
      allocate(this%d2Fdz2(gpC%xsz(1),gpC%xsz(2),gpC%xsz(3)))
   end if

   ! initialize the scalar
   if (restartSim) then
       if (RejectScalarRestart) then
        call initScalar(this%gpC, inputFile, mesh, this%scalar_number, this%F) 
       else
        call this%readRestart(tid_restart)
       end if
   else
      call initScalar(this%gpC, inputFile, mesh, this%scalar_number, this%F) 
   end if 
   call this%spectC%fft(this%F, this%Fhat)
   
   ! initialize the scalar source term
   if (useSource) then
      allocate(this%source_hat(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)))
      call SetScalar_source(this%gpC, inputFile, mesh, this%scalar_number, this%rbuffxC(:,:,:,1)) 
      call this%spectC%fft(this%rbuffxC(:,:,:,1),this%source_hat)
   end if

   call this%dealias()
   call this%prep_scalar()
end subroutine

subroutine get_derivatives(this)
   class(scalar_igrid), intent(inout) :: this
   
   call transpose_x_to_y(this%F, this%rbuffyC(:,:,:,1), this%gpC)
   call transpose_y_to_z(this%rbuffyC(:,:,:,1),this%rbuffzC(:,:,:,1), this%gpC)
   call this%der%ddz_C2E(this%rbuffzC(:,:,:,1),this%rbuffzE(:,:,:,1), this%bc_bottom, this%bc_top)
   call this%der%interpz_C2E(this%rbuffzC(:,:,:,1),this%rbuffzE(:,:,:,2), this%bc_bottom, this%bc_top)
   call this%der%ddz_E2C(this%rbuffzE(:,:,:,2),this%rbuffzC(:,:,:,2), this%bc_bottom, this%bc_top)
   
   call transpose_z_to_y(this%rbuffzE(:,:,:,1), this%rbuffyE(:,:,:,1), this%gpE)
   call transpose_y_to_x(this%rbuffyE(:,:,:,1), this%dFdzE, this%gpE)

   call transpose_z_to_y(this%rbuffzC(:,:,:,2), this%rbuffyC(:,:,:,2), this%gpC)
   call transpose_y_to_x(this%rbuffyC(:,:,:,2), this%dFdzC, this%gpC)
   
   call this%spectC%mtimes_ik1_oop(this%Fhat, this%cbuffyC(:,:,:,1))
   call this%spectC%ifft(this%cbuffyC(:,:,:,1),this%dFdxC)

   call this%spectC%mtimes_ik2_oop(this%Fhat, this%cbuffyC(:,:,:,1))
   call this%spectC%ifft(this%cbuffyC(:,:,:,1),this%dFdyC)

   if (.not. this%isInviscid) then
      call this%der%d2dz2_C2C(this%rbuffzC(:,:,:,1), this%rbuffzC(:,:,:,2), this%bc_bottom, this%bc_top)
      call transpose_z_to_y(this%rbuffzC(:,:,:,2), this%rbuffyC(:,:,:,2), this%gpC)
      call transpose_y_to_x(this%rbuffyC(:,:,:,2), this%d2Fdz2, this%gpC)
   end if

end subroutine 

subroutine readRestart(this, tid)
   use decomp_2d_io
   use reductions, only: p_minval, p_maxval
   use exits, only: message, message_min_max
   class(scalar_igrid), intent(inout) :: this
   integer, intent(in) :: tid
   character(len=clen) :: tempname, fname


   write(tempname,"(A7,A4,I2.2,A7,I2.2,A1,I6.6)") "RESTART", "_Run", this%RunID, "_SCALAR",this%scalar_number,".",tid
   fname = this%InputDataDir(:len_trim(this%InputDataDir))//"/"//trim(tempname)
   call decomp_2d_read_one(1,this%F,fname, this%gpC)
   
   call message(1,"Restart data successfully read for scalar number", this%scalar_number)
   call message_min_max(2,"Bounds for F:", p_minval(minval(this%F)), p_maxval(maxval(this%F)))

end subroutine 

subroutine dumpRestart(this, tid)
   use decomp_2d_io
   use exits, only: message
   class(scalar_igrid), intent(in) :: this
   integer, intent(in) :: tid
   character(len=clen) :: tempname, fname

   write(tempname,"(A7,A4,I2.2,A7,I2.2,A1,I6.6)") "RESTART", "_Run", this%RunID, "_SCALAR",this%scalar_number,".",tid
   fname = this%OutputDataDir(:len_trim(this%OutputDataDir))//"/"//trim(tempname)

   call message(0,"Dumping restart for scalar number", this%scalar_number)
   call decomp_2d_write_one(1,this%F,fname, this%gpC)

end subroutine 

subroutine dump_planes(this, tid, pid, dirid, dirlabel)
   use decomp_2d_io
   class(scalar_igrid), intent(in) :: this
   integer, intent(in) :: pid, tid, dirid
   character(len=2), intent(in) :: dirlabel

   character(len=clen) :: tempname, fname

   write(tempname,"(A3,I2.2,A2,I6.6,A2,I5.5,A2,I2.2)") "Run", this%RunID,"_t",tid,dirlabel,pid,".p", this%scalar_number
   fname = this%OutputDataDir(:len_trim(this%OutputDataDir))//"/"//trim(tempname)
   call decomp_2d_write_plane(1,this%F,dirid, pid, fname, this%gpC)

end subroutine

subroutine dumpScalarField(this, tid, viz_hdf5 )
   use decomp_2d_io
   class(scalar_igrid), intent(in) :: this
   integer, intent(in) :: tid
   character(len=clen) :: tempname, fname
   type(io_hdf5), intent(inout), optional :: viz_hdf5
   character(len=4) :: scalar_label

   if (present(viz_hdf5)) then
      write(scalar_label,"(A2,I2.2)") "sc", this%scalar_number
      call viz_hdf5%write_variable(this%F, scalar_label)
   else
      write(tempname,"(A3,I2.2,A3,I2.2,A2,I6.6,A4)") "Run", this%RunID, "_sc",this%scalar_number,"_t",tid,".out" 
      fname = this%OutputDataDir(:len_trim(this%OutputDataDir))//"/"//trim(tempname)
      call decomp_2d_write_one(1,this%F,fname, this%gpC)
   end if 

end subroutine

end module 
