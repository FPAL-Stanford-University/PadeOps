module angleControl
   use kind_parameters, only: rkind, clen
   use decomp_2d
   use spectralMod, only: spectral  
   use exits, only: message
   use constants, only: pi
   implicit none


   type :: angCont
      private
      !logical                                       :: TargetsAssociated = .false. 
      !real(rkind), dimension(:,:,:), pointer        :: u_target, v_target, w_target, T_target
      !real(rkind), dimension(:,:,:), allocatable    :: Fringe_kernel_cells, Fringe_kernel_edges
      !real(rkind)                                   :: Fringe_Lambda_x
      type(spectral),    pointer                    :: spectC, spectE
      type(decomp_info), pointer                    :: gpC, gpE, sp_gpC, sp_gpE
      real(rkind),    dimension(:,:,:,:), pointer   :: rbuffxC, rbuffxE
      complex(rkind), dimension(:,:,:,:), pointer   :: cbuffyC, cbuffyE
      !real(rkind)                                   :: LambdaFact
      integer :: z_ref !myFringeID = 1
      !logical :: useTwoFringex = .false. 
      real(rkind)                          :: phi, beta, phi_ref, sigma
      contains
         procedure :: init
         procedure :: destroy
         procedure :: update_RHS_control
   end type
    
contains
   
    subroutine update_RHS_control(this, dt, urhs, vrhs, wrhs, uC, vC, newTimestep, phi_n)
      class(angCont),                                                                        intent(inout)  :: this
      real(rkind),                                                                         intent(in)     :: dt
      real(rkind),    dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)),          intent(in)     :: uC, vC 
      !real(rkind),    dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)),          intent(in)     :: wE
      complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout)  :: urhs, vrhs
      complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout)  :: wrhs
      logical, intent(in) :: newTimestep
      integer :: nx, ny, i, j
      ! PID tuning parameters
      real(rkind) :: wControl
      real(rkind), intent(out) :: phi_n 
      nx = this%gpC%xsz(1)
      ny = this%gpC%xsz(2)

      ! PID controller
      if (newTimestep) then
      this%rbuffxC(:,:,:,1) = atan((vC / uC)) !* 180.d0 / pi
      phi_n = 0.d0 
      do j = 1, nx
              do i = 1, ny
                      phi_n = phi_n + this%rbuffxC(j,i,this%z_ref,1)
              enddo
      enddo
      phi_n = phi_n / (float(nx) * float(ny)) 
      wControl = (phi_n - this%phi) / dt
      this%phi = phi_n
      wControl = wControl + this%beta * (phi_n - this%phi_ref)

      !if (this%targetsAssociated) then
         ! u velocity source term 
         this%rbuffxC(:,:,:,1) =  vC * wControl
         call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))      
         urhs = urhs + this%cbuffyC(:,:,:,1)

         ! v velocity source term 
         this%rbuffxC(:,:,:,1) = - uC * wControl
         call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))      
         vrhs = vrhs + this%cbuffyC(:,:,:,1)
         
         ! w velocity source term 
         ! w source term is zero
         !this%rbuffxE(:,:,:,1) = (this%Lambdafact/dt)*(this%Fringe_kernel_edges)*(this%w_target - wE)
         !call this%spectE%fft(this%rbuffxE(:,:,:,1), this%cbuffyE(:,:,:,1))      
         !wrhs = wrhs + this%cbuffyE(:,:,:,1)
      end if 

   end subroutine


   subroutine destroy(this)
      class(angCont), intent(inout) :: this

      !deallocate(this%Fringe_kernel_cells)
      !deallocate(this%Fringe_kernel_edges)
      !this%TargetsAssociated = .false.
   end subroutine

   subroutine init(this, inputfile, spectC, spectE, gpC, gpE, rbuffxC, rbuffxE, cbuffyC, cbuffyE)
      use reductions, only: p_maxval
      use mpi
      class(angCont), intent(inout) :: this
      character(len=clen), intent(in) :: inputfile 
      type(decomp_info), intent(in), target :: gpC, gpE
      !real(rkind), dimension(gpC%xsz(1)), intent(in) :: x
      !real(rkind), dimension(gpC%xsz(2)), intent(in) :: y
      !real(rkind), intent(in) :: dx, dy
      type(spectral), intent(in), target :: spectC, spectE
      real(rkind),    dimension(:,:,:,:), target, intent(in) :: rbuffxC, rbuffxE
      complex(rkind), dimension(:,:,:,:), target, intent(in) :: cbuffyC, cbuffyE
      !integer, intent(in), optional :: fringeID
      real(rkind) :: phi_ref, beta, sigma, phi 

      !real(rkind) :: Lx, Ly, LambdaFact = 2.45d0, LambdaFact2 = 2.45d0
      !real(rkind) :: Fringe_yst = 1.d0, Fringe_yen = 1.d0
      !real(rkind) :: Fringe_xst = 0.75d0, Fringe_xen = 1.d0
      !real(rkind) :: Fringe_delta_st_x = 1.d0, Fringe_delta_st_y = 1.d0, Fringe_delta_en_x = 1.d0, Fringe_delta_en_y = 1.d0
      
      !real(rkind) :: Fringe1_delta_st_x = 1.d0, Fringe1_delta_en_x = 1.d0
      !real(rkind) :: Fringe2_delta_st_x = 1.d0, Fringe2_delta_en_x = 1.d0
      !real(rkind) :: Fringe1_xst = 0.75d0, Fringe1_xen = 1.d0
      !real(rkind) :: Fringe2_xst = 0.75d0, Fringe2_xen = 1.d0
      
      integer :: ioUnit = 10, i, j, k, nx, ierr, z_ref
      !real(rkind), dimension(:), allocatable :: x1, x2, Fringe_func, S1, S2, y1, y2
      !logical :: Apply_x_fringe = .true., Apply_y_fringe = .false.
      !namelist /FRINGE/ Apply_x_fringe, Apply_y_fringe, Fringe_xst, Fringe_xen, Fringe_delta_st_x, Fringe_delta_en_x, &
      !                  Fringe_delta_st_y, Fringe_delta_en_y, LambdaFact, LambdaFact2, Fringe_yen, Fringe_yst, Fringe1_delta_st_x, &
      !                  Fringe2_delta_st_x, Fringe1_delta_en_x, Fringe2_delta_en_x, Fringe1_xst, Fringe2_xst, Fringe1_xen, Fringe2_xen
    
      !if (present(fringeID)) then
      !   this%myFringeID = fringeID
      !   this%useTwoFringex = .true. 
      !end if
      !nx = gpC%xsz(1)
      real(rkind)  :: Lx = 1.d0, Ly = 1.d0, Lz = 1.d0, Tref = 0.d0, Tsurf0 = 1.d0, dTsurf_dt = -0.05d0, z0init = 1.d-4, frameAngle = 0.d0
      namelist /PROBLEM_INPUT/ Lx, Ly, Lz, Tref, Tsurf0, dTsurf_dt, z0init, frameAngle, beta, sigma, phi_ref, z_ref
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=ioUnit, NML=PROBLEM_INPUT)
      close(ioUnit)

      !Lx = maxval(x) + dx
      !Ly = p_maxval(maxval(y)) + dy
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
      this%beta = beta
      this%sigma = sigma
      this%phi_ref = phi_ref
      this%z_ref = z_ref
      this%phi = phi_ref
      
      !allocate(this%Fringe_kernel_cells(nx, gpC%xsz(2), gpC%xsz(3)))
      !allocate(this%Fringe_kernel_edges(nx, gpE%xsz(2), gpE%xsz(3)))
      
      !this%Fringe_kernel_cells = 0.d0
      !this%Fringe_kernel_edges = 0.d0

!   
!      if (this%usetwoFringex) then
!         select case (this%myFringeID)
!         case(1)
!            this%LambdaFact   = LambdaFact
!            Fringe_xst        = Fringe1_xst*Lx
!            Fringe_xen        = Fringe1_xen*Lx
!            Fringe_delta_st_x = Fringe1_delta_st_x*Lx
!            Fringe_delta_en_x = Fringe1_delta_en_x*Lx
!         case(2)
!            this%LambdaFact   = LambdaFact2
!            Fringe_xst        = Fringe2_xst*Lx
!            Fringe_xen        = Fringe2_xen*Lx
!            Fringe_delta_st_x = Fringe2_delta_st_x*Lx
!            Fringe_delta_en_x = Fringe2_delta_en_x*Lx
!         end select
!      else
!            this%LambdaFact    = LambdaFact
!            Fringe_xst        = Fringe_xst*Lx
!            Fringe_xen        = Fringe_xen*Lx
!            Fringe_delta_st_x = Fringe_delta_st_x*Lx
!            Fringe_delta_en_x = Fringe_delta_en_x*Lx
!      end if 
!      
!      if (Apply_x_fringe) then
!         ! x - direction fringe
!         allocate(x1         (nx))
!         allocate(x2         (nx))
!         allocate(S1         (nx))
!         allocate(S2         (nx))
!         allocate(Fringe_func(nx))
!     
!         x1 = ((x -  Fringe_xst)/Fringe_delta_st_x)
!         x2 = ((x -  Fringe_xen)/Fringe_delta_en_x) + 1.d0
!         call S_fringe(x1, S1)
!         call S_fringe(x2, S2)
!         Fringe_func = S1 - S2
!
!         do k = 1,this%gpC%xsz(3)
!            do j = 1,this%gpC%xsz(2)
!                this%Fringe_kernel_cells(:,j,k) = Fringe_func    
!            end do 
!         end do
!
!         do k = 1,this%gpE%xsz(3)
!            do j = 1,this%gpE%xsz(2)
!                this%Fringe_kernel_edges(:,j,k) = Fringe_func    
!            end do 
!         end do
!         deallocate(x1, x2, S1, S2, Fringe_func)
!      end if 
!
!      if (Apply_y_fringe) then
!         Fringe_yst        = Fringe_yst*Ly
!         Fringe_yen        = Fringe_yen*Ly
!         Fringe_delta_st_y = Fringe_delta_st_y*Ly
!         Fringe_delta_en_y = Fringe_delta_en_y*Ly
!         ! y direction fringe 1
!         allocate(y1         (this%gpC%xsz(2)))
!         allocate(y2         (this%gpC%xsz(2)))
!         allocate(S1         (this%gpC%xsz(2)))
!         allocate(S2         (this%gpC%xsz(2)))
!         allocate(Fringe_func(this%gpC%xsz(2)))
!     
!         y1 = ((y -  Fringe_yst)/Fringe_delta_st_y)
!         y2 = ((y -  Fringe_yen)/Fringe_delta_en_y) + 1.d0
!         call S_fringe(y1, S1)
!         call S_fringe(y2, S2)
!         Fringe_func = S1 - S2
!
!         do k = 1,this%gpC%xsz(3)
!            do j = 1,this%gpC%xsz(2)
!               do i = 1,nx
!                  this%Fringe_kernel_cells(i,j,k) = this%Fringe_kernel_cells(i,j,k) + Fringe_func(j)
!               end do 
!            end do 
!         end do
!
!         do k = 1,this%gpE%xsz(3)
!            do j = 1,this%gpE%xsz(2)
!               do i = 1,nx
!                  this%Fringe_kernel_edges(i,j,k) = this%Fringe_kernel_edges(i,j,k) + Fringe_func(j)
!               end do 
!            end do 
!         end do
!         deallocate(y1, y2, S1, S2, Fringe_func)
!      end if 
!
      call message(0, "Control initialized successfully.")

   end subroutine

end module 
