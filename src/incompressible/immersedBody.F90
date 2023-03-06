module immersedBodyMod
    use kind_parameters, only: rkind
    use decomp_2d
    use spectralMod, only: spectral  
    use constants, only: zero, one 
    implicit none

    type :: immersedBody 
        private 
        real(rkind), dimension(:,:,:), allocatable, public :: RbodyC, RbodyE
        real(rkind), dimension(:,:,:), allocatable, public :: utarget, vtarget, wtarget     
        type(decomp_info), pointer :: gpC, gpE, sp_gpC, sp_gpE 
        type(spectral), pointer :: spectC, spectE  
        real(rkind) :: taufact
        real(rkind), dimension(:,:,:), pointer :: rbuffxC, rbuffxE 
        complex(rkind), dimension(:,:,:), pointer :: cbuffyE, cbuffyC 

        contains 
            procedure :: init 
            procedure :: destroy
            procedure :: updateRHS 

    end type 

contains 

    subroutine init(this, spectC, spectE, gpC, gpE, sp_gpC, sp_gpE, taufact, rbuffxC, rbuffxE, cbuffyC, cbuffyE)
        class(immersedBody), intent(inout) :: this 
        class(spectral), target, intent(inout) :: spectC, spectE
        class(decomp_info), target, intent(inout) :: gpC, gpE, sp_gpC, sp_gpE 
        real(rkind), intent(in) :: taufact 
        complex(rkind), intent(inout), dimension(:,:,:), target :: cbuffyC, cbuffyE 
        real(rkind), intent(inout), dimension(:,:,:), target :: rbuffxC, rbuffxE 

        this%spectC => spectC 
        this%spectE => spectE

        this%gpC => gpC 
        this%gpE => gpE
        this%sp_gpC => sp_gpC 
        this%sp_gpE => sp_gpE

        this%taufact = taufact
        allocate(this%RbodyC(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
        allocate(this%RbodyE(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))

        this%RbodyC = zero 
        this%RbodyE = zero

        allocate(this%utarget(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
        allocate(this%vtarget(gpC%xsz(1), gpC%xsz(2), gpC%xsz(3)))
        allocate(this%wtarget(gpE%xsz(1), gpE%xsz(2), gpE%xsz(3)))

        this%rbuffxC => rbuffxC 
        this%rbuffxE => rbuffxE
        this%cbuffyC => cbuffyC 
        this%cbuffyE => cbuffyE

    end subroutine 


    subroutine destroy(this)
        class(immersedBody), intent(inout) :: this 
        
<<<<<<< HEAD
        if (allocated(this%RbodyC)) deallocate(this%RbodyC)
        ! TO DO (Ryan) : finish deallocating everything 
=======
        if (allocated(this%RbodyC))   deallocate(this%RbodyC)
        if (allocated(this%RbodyE))   deallocate(this%RbodyE)
        if (allocated(this%utarget))  deallocate(this%utarget)
        if (allocated(this%vtarget))  deallocate(this%vtarget)
        if (allocated(this%wtarget))  deallocate(this%wtarget)
        if (associated(this%gpC))     nullify(this%gpC)
        if (associated(this%gpE))     nullify(this%gpE)
        if (associated(this%sp_gpC))  nullify(this%sp_gpC)
        if (associated(this%sp_gpE))  nullify(this%sp_gpE)
        if (associated(this%spectC))  nullify(this%spectC)
        if (associated(this%spectE))  nullify(this%spectE)
        if (associated(this%rbuffxC)) nullify(this%rbuffxC)
        if (associated(this%rbuffxE)) nullify(this%rbuffxE)
        if (associated(this%cbuffyE)) nullify(this%cbuffyE)
        if (associated(this%cbuffyC)) nullify(this%cbuffyC)
>>>>>>> bca0e56daec689d11e59fe4531fb699163924c45
    end subroutine

    subroutine updateRHS(this, uC, vC, wE, urhs, vrhs, wrhs, dt) 
        class(immersedBody), intent(inout) :: this
        real(rkind), dimension(this%gpC%xsz(1),this%gpC%xsz(2),this%gpC%xsz(3)),intent(in) :: uC, vC
        real(rkind), dimension(this%gpE%xsz(1),this%gpE%xsz(2),this%gpE%xsz(3)),intent(in) :: wE

        complex(rkind), dimension(this%sp_gpC%ysz(1),this%sp_gpC%ysz(2),this%sp_gpC%ysz(3)), intent(inout) :: urhs, vrhs    
        complex(rkind), dimension(this%sp_gpE%ysz(1),this%sp_gpE%ysz(2),this%sp_gpE%ysz(3)), intent(inout) :: wrhs  
        real(rkind), intent(in) :: dt 
        real(rkind) :: tau 

        tau = one/(this%taufact*dt)

        this%rbuffxC =  tau*this%rBodyC*(this%uTarget - uC)
        call this%spectC%fft(this%rbuffxC, this%cbuffyC) 
        urhs = urhs + this%cbuffyC 

        this%rbuffxC =  tau*this%rBodyC*(this%vTarget - vC)
        call this%spectC%fft(this%rbuffxC, this%cbuffyC) 
        vrhs = vrhs + this%cbuffyC 

        this%rbuffxE =  tau*this%rBodyE*(this%wTarget - wE)
        call this%spectE%fft(this%rbuffxE, this%cbuffyE) 
        wrhs = wrhs + this%cbuffyE
    end subroutine 

end module 
