subroutine initfields(nx,ny,nz,proc_st,proc_en,proc_sz,dx,dy,dz,nvars,mesh,fields)
    use kind_parameters,  only: rkind
    use constants,        only: two,pi
    use IncompressibleGrid, only: u_index,v_index,w_index
    implicit none

    integer,                                                        intent(in)    :: nx,ny,nz,nvars
    integer, dimension(3),                                          intent(in)    :: proc_st,proc_en,proc_sz
    real(rkind),                                                    intent(in)    :: dx,dy,dz
    real(rkind), dimension(proc_sz(1),proc_sz(2),proc_sz(3),3),     intent(in)    :: mesh
    real(rkind), dimension(proc_sz(1),proc_sz(2),proc_sz(3),nvars), intent(inout) :: fields

    associate( u => fields(:,:,:,u_index), v => fields(:,:,:,v_index), w => fields(:,:,:,w_index))

    end associate

end subroutine

