subroutine meshgen(nx,ny,nz,proc_st,proc_en,proc_sz,dx,dy,dz,mesh)
    use kind_parameters,  only: rkind
    use constants,        only: two,pi
    implicit none

    integer,                                                    intent(in)    :: nx,ny,nz
    integer, dimension(3),                                      intent(in)    :: proc_st,proc_en,proc_sz
    real(rkind),                                                intent(inout) :: dx,dy,dz
    real(rkind), dimension(proc_sz(1),proc_sz(2),proc_sz(3),3), intent(inout) :: mesh

    integer :: i,j,k

    ! Create mesh from [0,2*pi)x[0,2*pi)x[0,2*pi) using nx, ny, nz points in x, y and z respectively
    ! Need to set x, y and z as well as  dx, dy and dz

    associate( x => mesh(:,:,:,1), y => mesh(:,:,:,2), z => mesh(:,:,:,3) )

        dx = two*pi/real(nx,rkind)
        dy = two*pi/real(ny,rkind)
        dz = two*pi/real(nz,rkind)

        do k=1,proc_sz(3)
            do j=1,proc_sz(2)
                do i=1,proc_sz(1)
                    x(i,j,k) = real( proc_st(1) + i - 1, rkind ) * dx
                    y(i,j,k) = real( proc_st(2) + j - 1, rkind ) * dy
                    z(i,j,k) = real( proc_st(3) + k - 1, rkind ) * dz
                end do
            end do
        end do

    end associate

end subroutine

