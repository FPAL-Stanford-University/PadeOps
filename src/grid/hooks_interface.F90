module hooks
    use kind_parameters, only: rkind 
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
        subroutine initfields(nx     , ny     , nz     , &
                              proc_st, proc_en, proc_sz, &
                              dx     , dy     , dz     , &
                              nvars, mesh, fields)
            import :: rkind
            integer, intent(in) :: nx, ny, nz
            real(rkind), intent(inout) :: dx, dy, dz
            integer, dimension(3), intent(in) :: proc_st, proc_en, proc_sz
            integer, intent(in) :: nvars 
            real(rkind), dimension(proc_sz(1),proc_sz(2),proc_sz(3),3), intent(in) :: mesh
            real(rkind), dimension(proc_sz(1),proc_sz(2),proc_sz(3),nvars), intent(inout) :: fields

        end subroutine 
    end interface


end module 
