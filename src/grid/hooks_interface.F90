module hooks
    use kind_parameters, only: rkind
    use decomp_2d, only: decomp_info 
    implicit none

    interface meshgen
        subroutine meshgen(nx     , ny     , nz     , &
                           proc_st, proc_en, proc_sz, &
                           dx     , dy     , dz     , &
                           mesh)
            import :: rkind
            integer, intent(in) :: nx, ny, nz
            real(rkind), intent(inout) :: dx, dy, dz
            integer, dimension(3), intent(in) :: proc_st, proc_en, proc_sz
            real(rkind), dimension(proc_sz(1),proc_sz(2),proc_sz(3),3), intent(in) :: mesh

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
