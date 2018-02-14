module angleControl
   use kind_parameters, only: rkind, clen
   use decomp_2d
   use spectralMod, only: spectral  
   use exits, only: message
   use constants, only: pi
   use reductions, only: p_sum
   implicit none


   type :: angCont
      private
      !logical                                       :: TargetsAssociated = .false. 
      !real(rkind), dimension(:,:,:), pointer        :: u_target, v_target, w_target, T_target
      !real(rkind), dimension(:,:,:), allocatable    :: Fringe_kernel_cells, Fringe_kernel_edges
      !real(rkind)                                   :: Fringe_Lambda_x
      type(spectral),    pointer                    :: spectC, spectE
      type(decomp_info), pointer                    :: gpC, gpE, sp_gpC, sp_gpE
      real(rkind),    dimension(:,:,:,:), pointer   :: rbuffxC, rbuffxE, rbuffyC, rbuffzC
      complex(rkind), dimension(:,:,:,:), pointer   :: cbuffyC, cbuffyE
      !real(rkind)                                   :: LambdaFact
      integer :: z_ref, controlType !myFringeID = 1
      !logical :: useTwoFringex = .false. 
      real(rkind)                          :: phi, beta, phi_ref, sigma, wFilt, alpha, wFilt_n, angleTrigger
      contains
         procedure :: init
         procedure :: destroy
         procedure :: update_RHS_control
         procedure :: getPhi
    end type
    
contains
   
    pure function getPhi(this) result (val) 
      class(angCont), intent(in) :: this
      real(rkind) :: val
      val = this%phi
    end function  

    subroutine update_RHS_control(this, dt, urhs, vrhs, wrhs, uC, vC, newTimestep, phi_n, wFilt_n, deltaGalpha, z_hub, trigger)
      class(angCont),                                                                        intent(inout)  :: this
      real(rkind),                                                                         intent(in)     :: dt
      real(rkind),    dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)),          intent(in)     :: uC, vC 
      !real(rkind),    dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)),          intent(in)     :: wE
      complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout)  :: urhs, vrhs
      complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout)  :: wrhs
      logical, intent(in) :: newTimestep
      integer :: nx, ny, i, j
      ! PID tuning parameters
      real(rkind) :: wControl_n, vM, uM
      real(rkind), intent(out) :: wFilt_n
      real(rkind), intent(out) :: phi_n, deltaGalpha, trigger
      integer, intent(out) :: z_hub
      !real(rkind), intent(out) :: totalAngle
      nx = this%gpC%xsz(1)
      ny = this%gpC%ysz(2)

      ! PID controller
         !this%rbuffxC(:,:,:,1) = atan2(vC, uC) !* 180.d0 / pi
         !call transpose_x_to_y(this%rbuffxC(:,:,:,1),this%rbuffyC(:,:,:,1),this%gpC)
         !call transpose_y_to_z(this%rbuffyC(:,:,:,1),this%rbuffzC(:,:,:,1),this%gpC)
         !phi_n = p_sum(sum(this%rbuffzC(:,:,this%z_ref,1))) / (real(nx,rkind) * real(ny,rkind))      
         !this%rbuffxC(:,:,:,1) = atan2(vC, uC) !* 180.d0 / pi
         ! Compute the angle at hub height
         call transpose_x_to_y(uC,this%rbuffyC(:,:,:,1),this%gpC)
         call transpose_y_to_z(this%rbuffyC(:,:,:,1),this%rbuffzC(:,:,:,1),this%gpC)
         uM = p_sum(sum(this%rbuffzC(:,:,this%z_ref,1))) / (real(nx,rkind) * real(ny,rkind))
         call transpose_x_to_y(vC,this%rbuffyC(:,:,:,1),this%gpC)
         call transpose_y_to_z(this%rbuffyC(:,:,:,1),this%rbuffzC(:,:,:,2),this%gpC)
         !phi_n = p_sum(sum(this%rbuffzC(:,:,this%z_ref,1))) / (real(nx,rkind) * real(ny,rkind))      
         vM = p_sum(sum(this%rbuffzC(:,:,this%z_ref,2))) / (real(nx,rkind) * real(ny,rkind))        
         phi_n = atan2(vM, uM)
         z_hub = this%z_ref
         trigger = this%angleTrigger 

      if (newTimestep .AND. abs(phi_n * 180.d0 / pi) > this%angleTrigger) then
         
         if (this%controlType == 1) then
            ! Meneveau 2014 psuedo force paper
            ! Rotation rate
            wControl_n = (phi_n - this%phi) / dt 
            ! Update the total angle and stored angles
            !totalAngle = totalAngle + (phi_n - this%phi)
            this%phi = phi_n
            ! First order time filter         
            !wFilt_n = (dt*wControl_n + dt*this%wControl + this%wFilt*(2.d0*this%sigma - dt)) / (this%sigma*2.d0 + dt) 
            wFilt_n = (1.d0 - (dt/(this%sigma+dt))) * this%wFilt + dt/(this%sigma + dt) * wControl_n
            this%wFilt = wFilt_n
            this%wFilt_n = this%alpha * wFilt_n + this%beta * (phi_n - this%phi_ref)
            wFilt_n = this%wFilt_n 
         elseif (this%controlType == 2) then
            ! Control geostrophic velocity directly
            wControl_n = (phi_n - this%phi)
            this%phi = phi_n
            deltaGalpha = (1.d0 - (dt/(this%sigma+dt))) * this%wFilt + dt/(this%sigma + dt) * wControl_n            
            this%wFilt = deltaGalpha 
            deltaGalpha = this%alpha * deltaGalpha + this%beta * (phi_n - this%phi_ref)
            deltaGalpha = deltaGalpha * pi / 180.d0
            wFilt_n = 0.d0    

         endif        
     end if
         ! Update the RHS 
         this%rbuffxC(:,:,:,1) =  2.d0 * vC * this%wFilt_n 
         call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))      
         urhs = urhs + this%cbuffyC(:,:,:,1)    
         this%rbuffxC(:,:,:,1) = -2.d0 * uC * this%wFilt_n   
         call this%spectC%fft(this%rbuffxC(:,:,:,1), this%cbuffyC(:,:,:,1))      
         vrhs = vrhs + this%cbuffyC(:,:,:,1)
         !!!!!!!!!!!!!!!!!!!!!!!
         ! Here I added the factor of 2 to deltaGalpha
         !!!!!!!!!!!!!!!!!!!!!!!
         deltaGalpha = 2.d0 * this%wFilt_n * dt * 180.d0 / pi
          

   end subroutine


   subroutine destroy(this)
      class(angCont), intent(inout) :: this

      !deallocate(this%Fringe_kernel_cells)
      !deallocate(this%Fringe_kernel_edges)
      !this%TargetsAssociated = .false.
   end subroutine

   subroutine init(this, inputfile, spectC, spectE, gpC, gpE, rbuffxC, rbuffxE, cbuffyC, cbuffyE, rbuffyC, rbuffzC, phiRestart)
      use reductions, only: p_maxval
      use mpi
      class(angCont), intent(inout) :: this
      character(len=clen), intent(in) :: inputfile 
      type(decomp_info), intent(in), target :: gpC, gpE
      !real(rkind), dimension(gpC%xsz(1)), intent(in) :: x
      !real(rkind), dimension(gpC%xsz(2)), intent(in) :: y
      !real(rkind), intent(in) :: dx, dy
      type(spectral), intent(in), target :: spectC, spectE
      real(rkind),    dimension(:,:,:,:), target, intent(in) :: rbuffxC, rbuffxE, rbuffyC, rbuffzC
      complex(rkind), dimension(:,:,:,:), target, intent(in) :: cbuffyC, cbuffyE
      !integer, intent(in), optional :: fringeID
      real(rkind) :: phi_ref, beta, sigma, phi, alpha , angleTrigger
      integer :: controlType
      real(rkind), intent(in) :: phiRestart
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
      !real(rkind)  :: Lx = 1.d0, Ly = 1.d0, Lz = 1.d0, Tref = 0.d0, Tsurf0 = 1.d0, dTsurf_dt = -0.05d0, z0init = 1.d-4, frameAngle = 0.d0
      !namelist /PROBLEM_INPUT/ Lx, Ly, Lz, Tref, Tsurf0, dTsurf_dt, z0init, frameAngle, beta, sigma, phi_ref, z_ref
      namelist /CONTROL/ beta, sigma, phi_ref, z_ref, alpha, controlType, angleTrigger
      !open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      !read(unit=ioUnit, NML=CONTROL)
      !close(ioUnit)
      ioUnit = 11
      open(unit=ioUnit, file=trim(inputfile), form='FORMATTED', iostat=ierr)
      read(unit=ioUnit, NML=CONTROL)
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
      this%rbuffyC => rbuffyC
      this%rbuffzC => rbuffzC
      this%beta = beta
      this%sigma = sigma
      this%phi_ref = phi_ref
      this%z_ref = z_ref
      this%wFilt = 0.d0     
      this%alpha = alpha
      this%controlType = controlType
      this%phi = phiRestart
      this%wFilt_n = 0.d0
      this%angleTrigger = angleTrigger
      call message(0, "Control initialized successfully.")

   end subroutine

end module 
