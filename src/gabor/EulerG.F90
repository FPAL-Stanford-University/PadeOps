module EulerG_mod
    use kind_parameters, only: rkind, clen  
    use decomp_2d  
    use Efield_mod, only: EField 
    use exits, only: GracefulExit 
    implicit none 

    type :: EulerG
        private
        integer :: nlevels 
        type(Efield), dimension(:), allocatable :: fields
        type(decomp_info), pointer :: gp 
        
        contains 
            procedure :: init 
            procedure :: agglomerate
            procedure :: writeFields 
            procedure :: destroy
            procedure :: setField 
            procedure :: getPointersToFields 
            procedure :: allocateFieldForLevel 
            procedure :: getDeltaForLevel 
            procedure :: getDomainRangeForLevel
            procedure :: ResetVelocityToZero 
    end type 

contains 

subroutine init(this, gpLES, xLim, yLim, zLim, nlevels)
    class(EulerG), intent(inout), target :: this
    type(decomp_info), intent(in), target :: gpLES 
    real(rkind), dimension(2) :: xLim, yLim, zLim 
    integer, intent(in) :: nLevels 
    integer :: iLevel
    type(decomp_info), pointer :: myGP

    if (nlevels < 1) call GracefulExit("nlevels in EulerG needs to be at least 1", 124)
    this%nlevels = nlevels    
    allocate(this%fields(nlevels))
    this%gp => gpLES 

    call this%fields(1)%init(this%gp, xLim, yLim, zLim)

    do iLevel = 2,nLevels
        myGP => this%fields(ilevel-1)%gpFinerLevel 
        call this%fields(ilevel)%init(myGP, xLim, yLim, zLim)
    end do 

    call this%ResetVelocityToZero()

end subroutine

subroutine setField(this, u, v, w, level)
    class(EulerG), intent(inout) :: this 
    real(rkind), dimension(:,:,:), intent(in) :: u, v, w
    integer, intent(in) :: level 

    this%fields(level)%u = u 
    this%fields(level)%v = v
    this%fields(level)%w = w 
end subroutine 

subroutine agglomerate(this)
    class(EulerG), intent(inout) :: this 
    integer :: iLevel 

    do iLevel = 1,this%nLevels-1
        call this%fields(ilevel)%aggDown(this%fields(ilevel+1)%u, this%fields(ilevel+1)%v, this%fields(ilevel+1)%w, this%fields(ilevel+1)%buffer)
    end do 
    
end subroutine 

subroutine getDeltaForLevel(this, delta, level)
    class(EulerG), intent(in) :: this 
    integer, intent(in) :: level 
    real(rkind), dimension(3), intent(out) :: delta 

    delta = this%fields(level)%delta 

end subroutine 

subroutine getDomainRangeForLevel(this, xRange, yRange, zRange, level)
    class(EulerG), intent(in) :: this 
    integer, intent(in) :: level 
    real(rkind), dimension(2), intent(out) :: xRange, yRange, zRange 

    xRange = this%fields(level)%xRange 
    yRange = this%fields(level)%yRange 
    zRange = this%fields(level)%zRange 
end subroutine 

subroutine destroy(this)
    class(EulerG), intent(inout) :: this 
    integer :: iLevel
    if (allocated(this%fields)) then 
        do iLevel = 1,this%nLevels
            call this%fields(iLevel)%destroy()
        end do
        deallocate(this%fields) 
    end if 

end subroutine 

subroutine allocateFieldForLevel(this, f, level)
    class(EulerG), intent(in) :: this 
    integer, intent(in) :: level 
    real(rkind), dimension(:,:,:), allocatable, intent(out) :: f

    if (allocated(f)) deallocate(f)
    allocate(f(this%fields(level)%gp%xsz(1), this%fields(level)%gp%xsz(2), this%fields(level)%gp%xsz(3)))
end subroutine 

subroutine getPointersToFields(this, u, v, w, level)
    class(EulerG), intent(inout), target :: this 
    real(rkind), dimension(:,:,:), pointer, intent(out) :: u, v, w 
    integer, intent(in) :: level 

    u => this%fields(level)%u
    v => this%fields(level)%v
    w => this%fields(level)%w

end subroutine 

subroutine writeFields(this, fname, level)
    use decomp_2d_io
    class(EulerG), intent(inout) :: this 
    integer :: level 
    character(len=*), intent(in) :: fname 
    character(len=clen) :: writefname
    character(len=17) :: tmp 

    write(tmp,"(A6,I2.2,A5,A4)") "_Level",level, "_uVEL", ".out"
    writefname = fname // tmp 
    call decomp_2d_write_one(1, this%fields(level)%u, trim(writefname), this%fields(level)%gp)

    write(tmp,"(A6,I2.2,A5,A4)") "_Level",level, "_vVEL", ".out"
    writefname = fname // tmp 
    call decomp_2d_write_one(1, this%fields(level)%v, trim(writefname), this%fields(level)%gp)

    write(tmp,"(A6,I2.2,A5,A4)") "_Level",level, "_wVEL", ".out"
    writefname = fname // tmp 
    call decomp_2d_write_one(1, this%fields(level)%w, trim(writefname), this%fields(level)%gp)

end subroutine 

subroutine ResetVelocityToZero(this)
    class(EulerG), intent(inout) :: this 
    integer :: ilevel 

    do ilevel = 1,this%nlevels
        this%fields(ilevel)%u = 0.d0 
        this%fields(ilevel)%v = 0.d0 
        this%fields(ilevel)%w = 0.d0 
    end do 

end subroutine 

end module 