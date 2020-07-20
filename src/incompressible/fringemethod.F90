module fringeMethod
   use kind_parameters, only: rkind, clen
   use decomp_2d
   use spectralMod, only: spectral  
   use exits, only: message, GracefulExit
   implicit none
   private
   public :: fringe

   type :: fringe 
      private
      logical                                       :: TargetsAssociated = .false. 
      real(rkind), dimension(:,:,:), pointer        :: u_target, v_target, w_target, T_target, F_target
      real(rkind), dimension(:,:,:), allocatable    :: Fringe_kernel_cells, Fringe_kernel_edges
      real(rkind)                                   :: Fringe_Lambda_x
      type(spectral),    pointer                    :: spectC, spectE
      type(decomp_info), pointer                    :: gpC, gpE, sp_gpC, sp_gpE
      real(rkind),    dimension(:,:,:,:), pointer   :: rbuffxC, rbuffxE
      complex(rkind), dimension(:,:,:,:), pointer   :: cbuffyC, cbuffyE
      real(rkind)                                   :: LambdaFact, LambdaFactPotTemp
      integer :: myFringeID = 1
      logical :: useTwoFringex = .false. 
      logical, public :: useFringeAsSponge_Scalar = .true. 
      logical :: firstCallComplete = .false.
      logical :: firstCallCompleteScalar = .false.
      logical :: useShiftedPeriodicBC = .false. 
      real(rkind) :: shift_distance, Lxeff
      integer :: ixst_1, ixen_1, ixst_2, ixen_2
      contains
         procedure :: init
         procedure :: destroy
         procedure :: addFringeRHS
         procedure :: associateFringeTargets
         procedure :: allocateTargetArray_Cells
         procedure :: allocateTargetArray_Edges
         procedure :: addFringeRHS_scalar
         procedure :: associateFringeTarget_scalar
         procedure :: getShiftedFields
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

   subroutine addFringeRHS(this, dt, urhs, vrhs, wrhs, uC, vC, wE, uhat, vhat, what, urhsF,vrhsF,wrhsF,addF)
      class(fringe),                                                                        intent(inout)  :: this
      real(rkind),                                                                         intent(in)     :: dt
      real(rkind),    dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)),          intent(in)     :: uC, vC 
      real(rkind),    dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)),          intent(in)     :: wE
      complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(in)     :: uhat, vhat
      complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(in)     :: what
      complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout)  :: urhs, vrhs
      complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout)  :: wrhs
      logical, intent(in), optional :: addF
      complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout), optional  :: urhsF, vrhsF
      complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout), optional  :: wrhsF
      
      logical :: AllOptionalsPresent

      if (present(urhsF) .and. present(vrhsF) .and. present(wrhsF) .and. present(addF)) then
         AllOptionalsPresent = .true. 
      else
         AllOptionalsPresent = .false. 
      end if

      if(this%useShiftedPeriodicBC) then
          call this%getShiftedFields(uhat, vhat, what) ! in uvw_target 
      endif

      if (this%targetsAssociated) then
         ! u velocity source term 
         this%rbuffxC(:,:,:,1) = (this%Lambdafact/dt)*(this%Fringe_kernel_cells)*(this%u_target - uC)
         call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))      
         urhs = urhs + this%cbuffyC(:,:,:,1)
         if (allOptionalsPresent) then
            if (addF) then
               urhsF = urhsF +  this%cbuffyC(:,:,:,1)
            else
               urhsF = this%cbuffyC(:,:,:,1)
            end if
         end if

         ! v velocity source term 
         this%rbuffxC(:,:,:,1) = (this%Lambdafact/dt)*(this%Fringe_kernel_cells)*(this%v_target - vC)
         call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))      
         vrhs = vrhs + this%cbuffyC(:,:,:,1)
         if (allOptionalsPresent) then
            if (addF) then
               vrhsF = vrhsF +  this%cbuffyC(:,:,:,1)
            else
               vrhsF = this%cbuffyC(:,:,:,1)
            end if
         end if
         
         ! w velocity source term 
         this%rbuffxE(:,:,:,1) = (this%Lambdafact/dt)*(this%Fringe_kernel_edges)*(this%w_target - wE)
         call this%spectE%fft(this%rbuffxE(:,:,:,1), this%cbuffyE(:,:,:,1))      
         wrhs = wrhs + this%cbuffyE(:,:,:,1)
         if (allOptionalsPresent) then
            if (addF) then
               wrhsF = wrhsF +  this%cbuffyE(:,:,:,1)
            else
               wrhsF = this%cbuffyE(:,:,:,1)
            end if
         end if
      end if 

   end subroutine

   subroutine addFringeRHS_scalar(this, dt, Frhs, F)
      class(fringe),                                                                      intent(inout)        :: this
      complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)),intent(inout)        :: Frhs  
      real(rkind),                                                                        intent(in)           :: dt
      real(rkind),    dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)),         intent(in)           :: F 

      if (this%firstCallCompleteScalar) then      
          if (this%useFringeAsSponge_Scalar) then
                this%rbuffxC(:,:,:,1) = -(this%LambdafactPotTemp/dt)*(this%Fringe_kernel_cells)*(F)
                call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))      
                Frhs = Frhs + this%cbuffyC(:,:,:,1)
          else
             if (associated(this%F_target)) then
                this%rbuffxC(:,:,:,1) = (this%LambdafactPotTemp/dt)*(this%Fringe_kernel_cells)*(this%F_target - F)
                call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))      
                Frhs = Frhs + this%cbuffyC(:,:,:,1)
             end if 
          end if 
      else
          this%firstCallCompleteScalar = .true.
      end if
   end subroutine

   subroutine destroy(this)
      class(fringe), intent(inout) :: this

      deallocate(this%Fringe_kernel_cells)
      deallocate(this%Fringe_kernel_edges)
      this%TargetsAssociated = .false.
   end subroutine

   subroutine associateFringeTarget_scalar(this, Ftarget)
      class(fringe), intent(inout) :: this
      real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in), target :: Ftarget

      this%F_target => Ftarget
      this%useFringeAsSponge_Scalar = .false. 

   end subroutine 

   subroutine associateFringeTargets(this, utarget, vtarget, wtarget, Ttarget)
      class(fringe), intent(inout) :: this
      real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in), target           :: utarget, vtarget
      real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)), intent(in), target           :: wtarget 
      real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)), intent(in), optional, target :: Ttarget 
 
      if(this%useShiftedPeriodicBC) then
          call GracefulExit("Fringe targets already associated at setup if &
            ShiftedPeriodicBC is true. Check input and problem files",999)
      endif
 
      this%u_target => utarget
      this%v_target => vtarget
      this%w_target => wtarget

      if (present(Ttarget)) then
         this%T_target => Ttarget
         this%useFringeAsSponge_Scalar = .false. 
      end if
      this%TargetsAssociated = .true.

      call message(0, "Fringe targets successfully associated.")
   end subroutine

   subroutine init(this, inputfile, dx, x, dy, y, spectC, spectE, gpC, gpE, rbuffxC, rbuffxE, cbuffyC, cbuffyE, fringeID)
      use reductions, only: p_minval, p_maxval
      use exits, only: message_min_max
       use decomp_2d_io
      use mpi
      class(fringe), intent(inout) :: this
      character(len=clen), intent(in) :: inputfile 
      type(decomp_info), intent(in), target :: gpC, gpE
      real(rkind), dimension(gpC%xsz(1)), intent(in) :: x
      real(rkind), dimension(gpC%xsz(2)), intent(in) :: y
      real(rkind), intent(in) :: dx, dy
      type(spectral), intent(in), target :: spectC, spectE
      real(rkind),    dimension(:,:,:,:), target, intent(in) :: rbuffxC, rbuffxE
      complex(rkind), dimension(:,:,:,:), target, intent(in) :: cbuffyC, cbuffyE
      integer, intent(in), optional :: fringeID

      real(rkind) :: Lx, Ly, LambdaFact = 2.45d0, LambdaFact2 = 2.45d0, LambdaFactPotTemp = 2.45d0
      real(rkind) :: Fringe_yst = 1.d0, Fringe_yen = 1.d0
      real(rkind) :: Fringe_xst = 0.75d0, Fringe_xen = 1.d0
      real(rkind) :: Fringe_delta_st_x = 1.d0, Fringe_delta_st_y = 1.d0, Fringe_delta_en_x = 1.d0, Fringe_delta_en_y = 1.d0
      
      real(rkind) :: Fringe1_delta_st_x = 1.d0, Fringe1_delta_en_x = 1.d0
      real(rkind) :: Fringe2_delta_st_x = 1.d0, Fringe2_delta_en_x = 1.d0
      real(rkind) :: Fringe1_xst = 0.75d0, Fringe1_xen = 1.d0
      real(rkind) :: Fringe2_xst = 0.75d0, Fringe2_xen = 1.d0
      
      integer :: ioUnit = 10, i, j, k, nx, ierr
      real(rkind), dimension(:), allocatable :: x1, x2, Fringe_func, S1, S2, y1, y2
      logical :: Apply_x_fringe = .true., Apply_y_fringe = .false.
      logical :: useShiftedPeriodicBC = .false.
      real(rkind) :: shift_distance = 0.5d0, Lxeff = 0.5d0

      namelist /FRINGE/ Apply_x_fringe, Apply_y_fringe, Fringe_xst, Fringe_xen, Fringe_delta_st_x, Fringe_delta_en_x, &
                        Fringe_delta_st_y, Fringe_delta_en_y, LambdaFact, LambdaFactPotTemp, LambdaFact2, Fringe_yen, Fringe_yst, Fringe1_delta_st_x, &
                        Fringe2_delta_st_x, Fringe1_delta_en_x, Fringe2_delta_en_x, Fringe1_xst, Fringe2_xst, Fringe1_xen, Fringe2_xen, &
                        useShiftedPeriodicBC, shift_distance, Lxeff
 
      if (present(fringeID)) then
         this%myFringeID = fringeID
         this%useTwoFringex = .true. 
      end if
      nx = gpC%xsz(1)
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=ioUnit, NML=FRINGE)
      close(ioUnit)

      if(useShiftedPeriodicBC .and. (Apply_y_fringe .or. this%useTwoFringex)) then
          call GracefulExit("Shifted Periodic BC supported only with one fringe in x direction",9999)
      endif

      if(useShiftedPeriodicBC .and. (.not. Apply_x_fringe)) then
          call GracefulExit("Shifted Periodic BC requires Apply_x_fringe to be true",9999)
      endif

      Lx = maxval(x) + dx
      Ly = p_maxval(maxval(y)) + dy
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

      this%useFringeAsSponge_Scalar = .true. 
      
      allocate(this%Fringe_kernel_cells(nx, gpC%xsz(2), gpC%xsz(3)))
      allocate(this%Fringe_kernel_edges(nx, gpE%xsz(2), gpE%xsz(3)))
      
      this%Fringe_kernel_cells = 0.d0
      this%Fringe_kernel_edges = 0.d0
      ! Maybe need a flag to ensure that it is define and the problem is
      ! stratified
      this%LambdaFactPotTemp = LambdaFactPotTemp
   
      if (this%usetwoFringex) then
         select case (this%myFringeID)
         case(1)
            this%LambdaFact   = LambdaFact
            Fringe_xst        = Fringe1_xst*Lx
            Fringe_xen        = Fringe1_xen*Lx
            Fringe_delta_st_x = Fringe1_delta_st_x*Lx
            Fringe_delta_en_x = Fringe1_delta_en_x*Lx
         case(2)
            this%LambdaFact   = LambdaFact2
            Fringe_xst        = Fringe2_xst*Lx
            Fringe_xen        = Fringe2_xen*Lx
            Fringe_delta_st_x = Fringe2_delta_st_x*Lx
            Fringe_delta_en_x = Fringe2_delta_en_x*Lx
         end select
      else
            this%LambdaFact    = LambdaFact
            Fringe_xst        = Fringe_xst*Lx
            Fringe_xen        = Fringe_xen*Lx
            Fringe_delta_st_x = Fringe_delta_st_x*Lx
            Fringe_delta_en_x = Fringe_delta_en_x*Lx
      end if 
      
      if (Apply_x_fringe) then
         ! x - direction fringe
         allocate(x1         (nx))
         allocate(x2         (nx))
         allocate(S1         (nx))
         allocate(S2         (nx))
         allocate(Fringe_func(nx))
     
         x1 = ((x -  Fringe_xst)/Fringe_delta_st_x)
         x2 = ((x -  Fringe_xen)/Fringe_delta_en_x) + 1.d0
         call S_fringe(x1, S1)
         call S_fringe(x2, S2)
         Fringe_func = S1 - S2

         do k = 1,this%gpC%xsz(3)
            do j = 1,this%gpC%xsz(2)
                this%Fringe_kernel_cells(:,j,k) = Fringe_func    
            end do 
         end do

         do k = 1,this%gpE%xsz(3)
            do j = 1,this%gpE%xsz(2)
                this%Fringe_kernel_edges(:,j,k) = Fringe_func    
            end do 
         end do
         deallocate(x1, x2, S1, S2, Fringe_func)
      end if 
       
      call message_min_max(1,"Bounds for Fringe_funcC:", p_minval(minval(this%Fringe_kernel_cells)), p_maxval(maxval(this%Fringe_kernel_cells)))
      call message_min_max(1,"Bounds for Fringe_funcE:", p_minval(minval(this%Fringe_kernel_edges)), p_maxval(maxval(this%Fringe_kernel_edges)))

      if (Apply_y_fringe) then
         Fringe_yst        = Fringe_yst*Ly
         Fringe_yen        = Fringe_yen*Ly
         Fringe_delta_st_y = Fringe_delta_st_y*Ly
         Fringe_delta_en_y = Fringe_delta_en_y*Ly
         ! y direction fringe 1
         allocate(y1         (this%gpC%xsz(2)))
         allocate(y2         (this%gpC%xsz(2)))
         allocate(S1         (this%gpC%xsz(2)))
         allocate(S2         (this%gpC%xsz(2)))
         allocate(Fringe_func(this%gpC%xsz(2)))
     
         y1 = ((y -  Fringe_yst)/Fringe_delta_st_y)
         y2 = ((y -  Fringe_yen)/Fringe_delta_en_y) + 1.d0
         call S_fringe(y1, S1)
         call S_fringe(y2, S2)
         Fringe_func = S1 - S2

         do k = 1,this%gpC%xsz(3)
            do j = 1,this%gpC%xsz(2)
               do i = 1,nx
                  this%Fringe_kernel_cells(i,j,k) = this%Fringe_kernel_cells(i,j,k) + Fringe_func(j)
               end do 
            end do 
         end do

         do k = 1,this%gpE%xsz(3)
            do j = 1,this%gpE%xsz(2)
               do i = 1,nx
                  this%Fringe_kernel_edges(i,j,k) = this%Fringe_kernel_edges(i,j,k) + Fringe_func(j)
               end do 
            end do 
         end do
         deallocate(y1, y2, S1, S2, Fringe_func)
      
         call message_min_max(1,"Bounds for Fringe_funcC:", p_minval(minval(this%Fringe_kernel_cells)), p_maxval(maxval(this%Fringe_kernel_cells)))
      
         call message_min_max(1,"Bounds for Fringe_funcE:", p_minval(minval(this%Fringe_kernel_edges)), p_maxval(maxval(this%Fringe_kernel_edges)))

      end if 
      
      where (this%Fringe_kernel_cells > 1) this%Fringe_kernel_cells = 1.d0 
      where (this%Fringe_kernel_edges > 1) this%Fringe_kernel_edges = 1.d0 

      this%firstCallComplete = .false.
      this%firstCallCompleteScalar = .false.

      if(useShiftedPeriodicBC) then
          this%useShiftedPeriodicBC = .true.
          this%shift_distance = shift_distance

          this%u_target => rbuffxC(:,:,:,2)
          this%v_target => rbuffxC(:,:,:,3)
          this%w_target => rbuffxE(:,:,:,2)
          !if (someLogical) then
          !   this%T_target => rbuffxC(:,:,:,4) !shifted_T
          !end if
          this%TargetsAssociated = .true.

          ! set indices so target can be in the interior in the x-direction
          ! (target ends at Lxeff)
          Lxeff = Lxeff*Lx
          this%ixen_1 = minloc(abs(x-Fringe_xen),1);   this%ixen_1 = min(this%ixen_1+1, nx)
          this%ixst_1 = minloc(abs(x-Fringe_xst),1);   this%ixst_1 = max(this%ixst_1-1, 1)

          this%ixen_2 = minloc(abs(x-Lxeff),1);   this%ixen_2 = min(this%ixen_2+1, nx)
          this%ixst_2 = this%ixen_2 - (this%ixen_1 - this%ixst_1); this%ixst_2 = max(this%ixst_2, 1)

          call message(1, "ixst_1 = ", this%ixst_1)
          call message(1, "ixen_1 = ", this%ixen_1)
          call message(1, "ixst_2 = ", this%ixst_2)
          call message(1, "ixen_2 = ", this%ixen_2)
          

          if( (this%ixen_2-this%ixst_2) .ne. (this%ixen_1-this%ixst_1)) then
              call GracefulExit("Lxeff, Fringe_xst and Fringe_xen incompatible; Check details for Shifted periodic BC to work.", 11)
          endif 

          call message(0, "Shifted periodic BC successfully set up in fringe init.")
      endif

      call message(0, "Fringe initialized successfully.")

   end subroutine

   subroutine getShiftedFields(this, uhat, vhat, what)
      !use reductions, only: p_minval, p_maxval
      !use exits, only: message_min_max
      !use decomp_2d_io
      class(fringe), intent(inout) :: this
      complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(in)  :: uhat, vhat
      complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(in)  :: what

      ! ----------Step 1 :: u_target----------------------
      call this%spectC%shift_in_y(uhat, this%cbuffyC(:,:,:,1), this%shift_distance)
      call this%spectC%ifft(this%cbuffyC(:,:,:,1), this%u_target)
      this%u_target(this%ixst_1:this%ixen_1,:,:) = this%u_target(this%ixst_2:this%ixen_2,:,:)
      this%u_target(1:this%ixst_1-1,:,:) = 0.0_rkind

      ! ----------Step 2 :: v_target----------------------
      call this%spectC%shift_in_y(vhat, this%cbuffyC(:,:,:,1), this%shift_distance)
      call this%spectC%ifft(this%cbuffyC(:,:,:,1), this%v_target)
      this%v_target(this%ixst_1:this%ixen_1,:,:) = this%v_target(this%ixst_2:this%ixen_2,:,:)
      this%v_target(1:this%ixst_1-1,:,:) = 0.0_rkind

      ! ----------Step 1 :: w_target----------------------
      call this%spectE%shift_in_y(what, this%cbuffyE(:,:,:,1), this%shift_distance)
      call this%spectE%ifft(this%cbuffyE(:,:,:,1), this%w_target)
      this%w_target(this%ixst_1:this%ixen_1,:,:) = this%w_target(this%ixst_2:this%ixen_2,:,:)
      this%w_target(1:this%ixst_1-1,:,:) = 0.0_rkind

   end subroutine


   pure subroutine S_fringe(x, output)
      real(rkind), dimension(:), intent(in)    :: x
      real(rkind), dimension(:), intent(out)   :: output
      integer :: i
      real(rkind) :: exparg

      do i = 1,size(x)
        if (x(i) .le. 0.d0) then
           output(i) = 0.d0
        else if (x(i) .ge. 1.d0) then
           output(i) = 1.d0
        else
           exparg = 1.d0/(x(i) - 1.d0 + 1.0D-32) + 1.d0/(x(i) + 1.0D-32)
           exparg = min(exparg,708.0d0) ! overflows if exparg > 709. need a better fix for this
           output(i) = 1.d0/(1.d0 + exp(exparg))
        end if
      end do

   end subroutine

end module 
