subroutine initfields(decomp,dx,dy,dz,inpDirectory,mesh,fields)
    use kind_parameters,    only: rkind, ckind
    use constants,          only: two,pi
    use IncompressibleGrid, only: u_index,v_index,w_index
    use hitCD_IO 
    use decomp_2d,          only: decomp_info
    implicit none
    type(decomp_info),               intent(in)    :: decomp
    character(len=*),                intent(in)    :: inpDirectory
    real(rkind),                     intent(in)    :: dx,dy,dz
    real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
    real(rkind), dimension(:,:,:,:), intent(inout) :: fields
    integer :: nx, ny, nz, ix1, ixn, iy1, iyn, iz1, izn

    ! Get global size
    nx = decomp%xsz(1); ny = decomp%ysz(2); nz = decomp%zsz(3)

    ! If base decomposition in in X
    ix1 = decomp%xst(1); iy1 = decomp%xst(2); iz1 = decomp%xst(3)
    ixn = decomp%xen(1); iyn = decomp%xen(2); izn = decomp%xen(3)


    call getHit3d_uvw(nx,ny,nz,fields,decomp,inpDirectory)


end subroutine

