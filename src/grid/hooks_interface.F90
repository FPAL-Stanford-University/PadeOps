module hooks
    use kind_parameters, only: rkind
    use decomp_2d, only: decomp_info 
    implicit none

    interface meshgen
        subroutine meshgen(decomp, dx, dy, dz, mesh)
            import :: rkind
            type(decomp_info), intent(in) :: decomp
            real(rkind), intent(inout) :: dx, dy, dz
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh

        end subroutine 
    end interface

    interface initfields
        subroutine initfields(decomp, dx, dy, dz, inpDirectory, mesh, fields)
            import :: rkind
            import :: decomp_info
            type(decomp_info), intent(in) :: decomp
            real(rkind), intent(inout) :: dx, dy, dz
            character(len=*), intent(in) :: inpDirectory
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(inout) :: fields

        end subroutine 
    end interface


end module 
