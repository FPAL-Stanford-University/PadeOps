subroutine initfields(decomp,dx,dy,dz,inpDirectory,mesh,fields)
    use kind_parameters,    only: rkind
    use constants,          only: two,pi
    use GridMod,            only: alloc_buffs
    use IncompressibleGrid, only: u_index,v_index,w_index
    use hitCD_IO 
    use decomp_2d,          only: decomp_info, transpose_x_to_y
    implicit none
    type(decomp_info),               intent(in)    :: decomp
    character(len=*),                intent(in)    :: inpDirectory
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    real(rkind), dimension(:,:,:,:), allocatable   :: xfields
    integer :: nx, ny, nz

    ! Get global size
    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! Allocate buffer in x decomp to read in data
    call alloc_buffs(xfields,3,'x',decomp)
    call getHit3d_uvw(nx,ny,nz,xfields,decomp,inpDirectory)

    fields(:,:,:,1:3) = xfields
    deallocate(xfields)

end subroutine

