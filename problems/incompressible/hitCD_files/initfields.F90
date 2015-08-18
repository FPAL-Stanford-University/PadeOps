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
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn

    ! Get global size
    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition in in Y
    ix1 = decomp%yst(1); iy1 = decomp%yst(2); iz1 = decomp%yst(3)
    ixn = decomp%yen(1); iyn = decomp%yen(2); izn = decomp%yen(3)

    ! Allocate buffer in x decomp to read in data
    call alloc_buffs(xfields,3,'x',decomp)
    call getHit3d_uvw(nx,ny,nz,xfields,decomp,inpDirectory)

    ! Transpose the data to Y decomp
    call transpose_x_to_y(xfields(:,:,:,1),fields(:,:,:,1),decomp)
    call transpose_x_to_y(xfields(:,:,:,2),fields(:,:,:,2),decomp)
    call transpose_x_to_y(xfields(:,:,:,3),fields(:,:,:,3),decomp)

    deallocate(xfields)

end subroutine

