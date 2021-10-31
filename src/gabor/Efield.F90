module Efield_mod
    use kind_parameters, only: rkind 
    use decomp_2d  
    use interpolatorMod, only: interpolator  
    use gridTools, only: linspace 

    implicit none 

    type :: Efield 
        real(rkind), dimension(:,:,:), allocatable :: u, v, w, buffer
        type(interpolator), allocatable  :: interpToFinerLevel 
        type(decomp_info), pointer :: gp
        type(decomp_info), allocatable  :: gpFinerLevel  

        contains 
        procedure :: init
        procedure :: AggDown    
        procedure :: destroy 
    end type 

contains 

subroutine init(this, gp, xDom, yDom, zDom)
    class(Efield), intent(inout) :: this 
    type(decomp_info), intent(in), target :: gp 
    real(rkind), dimension(2), intent(in) :: xDom, yDom, zDom 
    real(rkind), dimension(:), allocatable :: x, y, z, xFine, yFine, zFine
    integer :: nx, ny, nz
    real(rkind), parameter :: small = 1.d-13 
    
    this%gp => gp 
    nx = this%gp%xsz(1)-1 
    ny = this%gp%ysz(2)-1 
    nz = this%gp%zsz(3)-1 

    allocate(x(nx+1), y(ny+1), z(nz+1))
    allocate(xFine(2*nx+1), yFine(2*ny+1), zFine(2*nz+1))

    x      = linspace(xDom(1), xDom(2),   nx+1)
    xFine  = linspace(xDom(1)+small, xDom(2)-small, 2*nx+1)
    
    y      = linspace(yDom(1), yDom(2),   ny+1)
    yFine  = linspace(yDom(1)+small, yDom(2)-small, 2*ny+1)
    
    z      = linspace(zDom(1), zDom(2),   nz+1)
    zFine  = linspace(zDom(1)+small, zDom(2)-small, 2*nz+1)

    allocate(this%gpFinerLevel)
    call decomp_info_init(2*nx+1,2*ny+1,2*nz+1,this%gpFinerLevel)

    allocate(this%interpToFinerLevel)
    call this%interpToFinerLevel%init(this%gp, this%gpFinerLevel, x, y, z, xFine, yFine, zFine) 
   
    allocate(this%u(this%gp%xsz(1), this%gp%xsz(2), this%gp%xsz(3)))
    allocate(this%v(this%gp%xsz(1), this%gp%xsz(2), this%gp%xsz(3)))
    allocate(this%w(this%gp%xsz(1), this%gp%xsz(2), this%gp%xsz(3)))
    allocate(this%buffer(this%gp%xsz(1), this%gp%xsz(2), this%gp%xsz(3)))
    this%u = 0.d0 
    this%v = 0.d0 
    this%w = 0.d0 
end subroutine 

subroutine AggDown(this, uFiner, vFiner, wFiner, buffer)
    class(Efield), intent(inout) :: this 
    real(rkind), dimension(:,:,:), intent(inout) :: uFiner, vFiner, wFiner, buffer

    call this%interpToFinerLevel%LinInterp3D(this%u, uFiner, buffer)
    call this%interpToFinerLevel%LinInterp3D(this%v, vFiner, buffer)
    call this%interpToFinerLevel%LinInterp3D(this%w, wFiner, buffer)
end subroutine 

subroutine destroy(this)
    class(Efield), intent(inout) :: this 
    if (allocated(this%u)) deallocate(this%u)
    if (allocated(this%v)) deallocate(this%v)
    if (allocated(this%w)) deallocate(this%w)
    if (allocated(this%gpFinerLevel)) deallocate(this%gpFinerLevel)
    if (associated(this%gp)) nullify(this%gp)
    if (allocated(this%interpToFinerLevel)) then 
        call this%interpToFinerLevel%destroy()
        deallocate(this%interpToFinerLevel)
    end if 
end subroutine 

end module 