module fringeMethod
   use kind_parameters, only: rkind, clen
   use decomp_2d
   use spectralMod, only: spectral  
   use exits, only: message
   implicit none
   private
   public :: fringe

   type :: fringe 
      private
      logical                                       :: TargetsAssociated = .false. 
      real(rkind), dimension(:,:,:), pointer        :: u_target, v_target, w_target, T_target
      real(rkind), dimension(:,:,:), allocatable    :: Fringe_kernel_cells_x, Fringe_kernel_edges_x
      real(rkind)                                   :: Fringe_Lambda_x
      type(spectral),    pointer                    :: spectC, spectE
      type(decomp_info), pointer                    :: gpC, gpE, sp_gpC, sp_gpE
      real(rkind),    dimension(:,:,:,:), pointer   :: rbuffxC, rbuffxE
      complex(rkind), dimension(:,:,:,:), pointer   :: cbuffyC, cbuffyE
      real(rkind)                                   :: LambdaFact
      contains
         procedure :: init
         procedure :: destroy
         procedure :: addFringeRHS
         procedure :: associateFringeTargets
         procedure :: allocateTargetArray_Cells
         procedure :: allocateTargetArray_Edges
   end type
    
contains
   
   subroutine allocateTargetArray_Cells(this, array)
      class(fringe), intent(in)  :: this
      real(rkind), dimension(:,:,:), allocatable, intent(out) :: array

      allocate(array(this%gpC%xsz(1), this%gpC%xsz(2), this%gpC%xsz(3)))

   end subroutine
   
   subroutine allocateTargetArray_Edges(this, array)
      class(fringe), intent(in)  :: this
      real(rkind), dimension(:,:,:), allocatable, intent(out) :: array

      allocate(array(this%gpE%xsz(1), this%gpE%xsz(2), this%gpE%xsz(3)))

   end subroutine

   subroutine addFringeRHS(this, dt, urhs, vrhs, wrhs, uC, vC, wE)
      class(fringe),                                                                        intent(inout)  :: this
      real(rkind),                                                                         intent(in)     :: dt
      real(rkind),    dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)),          intent(in)     :: uC, vC 
      real(rkind),    dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)),          intent(in)     :: wE
      complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout)  :: urhs, vrhs
      complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout)  :: wrhs


      if (this%targetsAssociated) then
         ! u velocity source term 
         this%rbuffxC(:,:,:,1) = (this%Lambdafact/dt)*(this%Fringe_kernel_cells_x)*(this%u_target - uC)
         call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))      
         urhs = urhs + this%cbuffyC(:,:,:,1)

         ! v velocity source term 
         this%rbuffxC(:,:,:,1) = (this%Lambdafact/dt)*(this%Fringe_kernel_cells_x)*(this%v_target - vC)
         call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))      
         vrhs = vrhs + this%cbuffyC(:,:,:,1)
         
         ! w velocity source term 
         this%rbuffxE(:,:,:,1) = (this%Lambdafact/dt)*(this%Fringe_kernel_edges_x)*(this%w_target - wE)
         call this%spectE%fft(this%rbuffxE(:,:,:,1), this%cbuffyE(:,:,:,1))      
         wrhs = wrhs + this%cbuffyE(:,:,:,1)
      end if 

   end subroutine


   subroutine destroy(this)
      class(fringe), intent(inout) :: this

      deallocate(this%Fringe_kernel_cells_x)
      deallocate(this%Fringe_kernel_edges_x)
      this%TargetsAssociated = .false.
   end subroutine


   subroutine associateFringeTargets(this, utarget, vtarget, wtarget, Ttarget)
      class(fringe), intent(inout) :: this
      real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in), target           :: utarget, vtarget
      real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)), intent(in), target           :: wtarget 
      real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in), optional, target :: Ttarget 
  
      this%u_target => utarget
      this%v_target => vtarget
      this%w_target => wtarget

      if (present(Ttarget)) then
         this%T_target => Ttarget
      end if
      this%TargetsAssociated = .true.

      call message(0, "Fringe targets successfully associated.")
   end subroutine

   subroutine init(this, inputfile, dx, x, spectC, spectE, gpC, gpE, rbuffxC, rbuffxE, cbuffyC, cbuffyE)
      class(fringe), intent(inout) :: this
      character(len=clen), intent(in) :: inputfile 
      type(decomp_info), intent(in), target :: gpC, gpE
      real(rkind), dimension(gpC%xsz(1)), intent(in) :: x
      real(rkind), intent(in) :: dx
      type(spectral), intent(in), target :: spectC, spectE
      real(rkind),    dimension(:,:,:,:), target, intent(in) :: rbuffxC, rbuffxE
      complex(rkind), dimension(:,:,:,:), target, intent(in) :: cbuffyC, cbuffyE

      real(rkind) :: Lx, LambdaFact = 2.45d0
      real(rkind) :: Fringe_xst = 0.75d0, Fringe_xen = 1.d0, Fringe_delta_st = 0.2d0, Fringe_delta_en = 0.05d0
      integer :: ioUnit = 10, j, k, nx, ierr
      real(rkind), dimension(:), allocatable :: x1, x2, Fringe_func, S1, S2
      namelist /FRINGE/ Fringe_xst, Fringe_xen, Fringe_delta_st, Fringe_delta_en, LambdaFact 
     
      nx = gpC%xsz(1)
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=ioUnit, NML=FRINGE)
      close(ioUnit)

      Lx = maxval(x) + dx
      Fringe_xst      = Fringe_xst*Lx
      Fringe_xen      = Fringe_xen*Lx
      Fringe_delta_st = Fringe_delta_st*Lx
      Fringe_delta_en = Fringe_delta_en*Lx
      this%LambdaFact = LambdaFact
      
      allocate(this%Fringe_kernel_cells_x(nx, gpC%xsz(2), gpC%xsz(3)))
      allocate(this%Fringe_kernel_edges_x(nx, gpE%xsz(2), gpE%xsz(3)))
        
      allocate(x1         (nx))
      allocate(x2         (nx))
      allocate(S1         (nx))
      allocate(S2         (nx))
      allocate(Fringe_func(nx))
     

      this%gpC => gpC
      this%gpE => gpE
      this%spectC => spectC
      this%spectE => spectE
      this%sp_gpC => spectC%spectdecomp
      this%sp_gpE => spectE%spectdecomp
      this%rbuffxC => rbuffxC 
      this%rbuffxE => rbuffxE
      this%cbuffyC => cbuffyC 
      this%cbuffyE => cbuffyE
      
      x1 = ((x -  Fringe_xst)/Fringe_delta_st)
      x2 = ((x -  Fringe_xen)/Fringe_delta_en) + 1.d0
      call S_fringe(x1, S1)
      call S_fringe(x2, S2)
      Fringe_func = S1 - S2

      do k = 1,this%gpC%xsz(3)
         do j = 1,this%gpC%xsz(2)
             this%Fringe_kernel_cells_x(:,j,k) = Fringe_func    
         end do 
      end do

      do k = 1,this%gpE%xsz(3)
         do j = 1,this%gpE%xsz(2)
             this%Fringe_kernel_edges_x(:,j,k) = Fringe_func    
         end do 
      end do
      deallocate(x1, x2, S1, S2, Fringe_func)

      
      call message(0, "Fringe initialized successfully.")

   end subroutine


   pure subroutine S_fringe(x, output)
      real(rkind), dimension(:), intent(in)    :: x
      real(rkind), dimension(:), intent(out)   :: output
      integer :: i

      do i = 1,size(x)
        if (x(i) .le. 0.d0) then
           output(i) = 0.d0
        else if (x(i) .ge. 1.d0) then
           output(i) = 1.d0
        else
           output(i) = 1.d0/(1.d0 + exp((1.d0/(x(i) - 1.d0)) + (1.d0/(x(i)))))
        end if
      end do

   end subroutine


end module 
