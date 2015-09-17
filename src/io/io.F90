module io_stuff

    use kind_parameters, only: rkind, clen
    use decomp_2d,       only: decomp_info
    implicit none

    type, abstract :: io

        integer :: nprimary
        integer :: vizcount
        character(len=clen) :: file_prefix
        character(len=clen) :: vizdir
        character(len=clen), dimension(:), allocatable :: primary_names

    contains
    
        procedure(init_interface),              deferred :: init
        procedure(destroy_interface),           deferred :: destroy
        procedure(WriteViz_interface),          deferred :: WriteViz

    end type

    abstract interface

        subroutine init_interface(this, vizdir_, file_prefix_, nprimary_, primary_names_)
            import :: io
            import :: clen
            class(io), intent(inout) :: this
            character(len=*),    intent(in) :: vizdir_
            character(len=*),    intent(in) :: file_prefix_
            integer,             intent(in) :: nprimary_
            character(len=*), dimension(nprimary_), intent(in) :: primary_names_
        end subroutine

        subroutine destroy_interface(this)
            import :: io
            class(io), intent(inout) :: this
        end subroutine

        subroutine WriteViz_interface(this, gp, mesh, primary, secondary, secondary_names)
            import :: io
            import :: decomp_info
            import :: rkind
            import :: clen
            class(io), intent(inout) :: this
            class(decomp_info), intent(in) :: gp
            real(rkind), dimension(gp%ysz(1),gp%ysz(2),gp%ysz(3),3), intent(in) :: mesh
            real(rkind), dimension(gp%ysz(1),gp%ysz(2),gp%ysz(3),this%nprimary), intent(in) :: primary
            real(rkind), dimension(:,:,:,:), intent(in), optional :: secondary
            character(len=*), dimension(:), intent(in), optional :: secondary_names
        end subroutine

    end interface

end module
