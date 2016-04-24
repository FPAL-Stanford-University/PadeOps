module wallmodelMod
    use kind_parameters, only: rkind, clen
    use constants, only: imi, zero,one,two,three,half,fourth 
    use decomp_2d

    implicit none
    private
    public :: wallmodel

    real(rkind), parameter :: kappa = 0.4_rkind

    type wallmodel
        integer :: modelID 
        real(rkind) :: mfactor
        real(rkind) :: meanSpeed 
        logical :: isInitialized = .false. 
        real(rkind), dimension(:,:,:), pointer :: xbuff, ybuff, zbuff1, zbuff2
        type(decomp_info), pointer :: decomp
        integer :: nxX, nyX, nzX
    contains
        procedure :: init
        procedure :: destroy
        procedure :: updateWallStress       
    end type

contains

    subroutine init(this,dz,inputfile, decompC, xbuff, ybuff, zbuff)
        class(wallmodel), intent(inout) :: this
        real(rkind), intent(in) :: dz
        character(len=*), intent(in) :: inputfile
        real(rkind), dimension(:,:,:,:), intent(in), target :: xbuff, ybuff, zbuff
        type(decomp_info), target, intent(in) :: decompC
        real(rkind) :: z0 = 0.1_rkind, H = 1000._rkind, G = 15._rkind, f = 1.454d-4, dpdx
        integer :: iounit

        namelist /PBLINPUT/ H, z0, G, f, dpdx

        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=PBLINPUT)
        close(ioUnit)    

        ! Non - dimensionalize z0
        z0 = z0/H
        this%mfactor = (kappa/log((dz/two)/z0))**2
    
        this%decomp => decompC   
        this%xbuff => xbuff(:,:,:,1) 
        this%ybuff => ybuff(:,:,:,1) 
        this%zbuff1 => zbuff(:,:,:,1) 
        this%zbuff2 => zbuff(:,:,:,2) 
        
        this%nxX = decompC%xsz(1)
        this%nyX = decompC%xsz(2)
        this%nzX = decompC%xsz(3)

        this%IsInitialized = .true. 
    end subroutine
   
    subroutine destroy(this)
        class(wallmodel), intent(inout) :: this

        nullify(this%xbuff, this%ybuff, this%zbuff1, this%zbuff2)
        this%IsInitialized = .false. 
    end subroutine 

    subroutine updateWallStress(this, tau13, tau23, u, v, umn)
        class(wallmodel), intent(inout) :: this
        real(rkind), dimension(this%nxX,this%nyX,this%nzX), intent(inout) :: tau13, tau23
        real(rkind), dimension(this%nxX,this%nyX,this%nzX), intent(in) :: u, v
        real(rkind), intent(in) :: umn

        !print*, "Tau13 before:"
        !print*, tau13(3,2,1:3)
        call transpose_x_to_y(tau13,this%ybuff,this%decomp)
        call transpose_y_to_z(this%ybuff,this%zbuff1,this%decomp)
        call transpose_x_to_y(u,this%ybuff,this%decomp)
        call transpose_y_to_z(this%ybuff,this%zbuff2,this%decomp)
        this%zbuff1(:,:,1) =  this%zbuff2(:,:,1)*this%mfactor*umn
        !if (nrank == 0) then
        !        print*, this%zbuff1(1,1,1:3)
        !end if 
        call transpose_z_to_y(this%zbuff1,this%ybuff,this%decomp)
        call transpose_y_to_x(this%ybuff,tau13,this%decomp)
        !print*, "Tau13 after:"
        !print*, tau13(3,2,1:3)
        !if(nrank==0) write(*,'(a,3(1x,e19.12))') "===", umn, sqrt(this%mfactor)*umn


        !print*, "Tau23 before:"
        !print*, tau23(3,2,1:3)
        call transpose_x_to_y(tau23,this%ybuff,this%decomp)
        call transpose_y_to_z(this%ybuff,this%zbuff1,this%decomp)
        call transpose_x_to_y(v,this%ybuff,this%decomp)
        call transpose_y_to_z(this%ybuff,this%zbuff2,this%decomp)
        this%zbuff1(:,:,1) = this%zbuff2(:,:,1)*this%mfactor*umn
        !if (nrank == 0) then
        !        print*, this%zbuff1(1,1,1:3)
        !end if 
        call transpose_z_to_y(this%zbuff1,this%ybuff,this%decomp)
        call transpose_y_to_x(this%ybuff,tau23,this%decomp)
        !print*, "Tau23 after:"
        !print*, tau23(3,2,1:3)
        !print*, "============================"

    end subroutine
     

end module 
