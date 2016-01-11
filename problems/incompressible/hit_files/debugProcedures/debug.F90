module debugMod
    use IncompressibleGrid, only: igrid
    use kind_parameters, only: rkind
    use decomp_2d, only: transpose_y_to_x
    use gridtools,            only: alloc_buffs

    implicit none 
contains
    subroutine debug(igp)
        class(igrid), intent(in) :: igp
        real(rkind), dimension(:,:,:,:), allocatable :: fieldinX

        call alloc_buffs(fieldinX,1,'x',igp%decomp)
        call transpose_y_to_x(igp%v,fieldinX(:,:,:,1),igp%decomp)
        print*, fieldinX(1:5,1,1,1)
        deallocate(fieldinX)
    end subroutine 
end module
