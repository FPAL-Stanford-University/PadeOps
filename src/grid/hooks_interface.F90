module hooks
    use kind_parameters, only: rkind
    use decomp_2d, only: decomp_info
    implicit none

    interface meshgen
        subroutine meshgen(decomp, dx, dy, dz, mesh)
            import :: rkind
            import :: decomp_info
            type(decomp_info), intent(in) :: decomp
            real(rkind), intent(inout) :: dx, dy, dz
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh

        end subroutine 
    end interface

    interface initfields
        subroutine initfields(decomp, dx, dy, dz, inpDirectory, mesh, fields, rho0, mu, gam, PInf)
            import :: rkind
            import :: decomp_info
            type(decomp_info), intent(in) :: decomp
            real(rkind), intent(inout) :: dx, dy, dz
            character(len=*), intent(in) :: inpDirectory
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(inout) :: fields
            real(rkind),           optional, intent(inout) :: rho0, mu, gam, PInf

        end subroutine 
    end interface


    interface initfields_stagg
        subroutine initfields_stagg(decompC, decompE, dx, dy, dz, inpDirectory, mesh, fieldsC, fieldsE, u_g)
            import :: rkind
            import :: decomp_info
            type(decomp_info), intent(in) :: decompC
            type(decomp_info), intent(in) :: decompE
            real(rkind), intent(inout) :: dx, dy, dz
            character(len=*), intent(in) :: inpDirectory
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(inout) :: fieldsC
            real(rkind), dimension(:,:,:,:), intent(inout) :: fieldsE
            real(rkind), intent(out), optional :: u_g

        end subroutine 
    end interface


    interface hook_output
        subroutine hook_output(decomp,dx,dy,dz,outputdir,mesh,fields,tsim,vizcount)
            import :: rkind
            import :: decomp_info
            character(len=*),                intent(in) :: outputdir
            type(decomp_info),               intent(in) :: decomp
            real(rkind),                     intent(in) :: dx,dy,dz,tsim
            integer,                         intent(in) :: vizcount
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(in) :: fields

        end subroutine
    end interface

    interface hook_bc
        subroutine hook_bc(decomp,mesh,fields,tsim)
            import :: rkind
            import :: decomp_info
            type(decomp_info),               intent(in)    :: decomp
            real(rkind),                     intent(in)    :: tsim
            real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
            real(rkind), dimension(:,:,:,:), intent(inout) :: fields

        end subroutine
    end interface

    interface hook_timestep
        subroutine hook_timestep(decomp,mesh,fields,tsim)
            import :: rkind
            import :: decomp_info
            type(decomp_info),               intent(in) :: decomp
            real(rkind),                     intent(in) :: tsim
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(in) :: fields

        end subroutine
    end interface

    interface hook_source
        subroutine hook_source(decomp,mesh,fields,tsim,rhs,rhsg)
            import :: rkind
            import :: decomp_info
            type(decomp_info),               intent(in)    :: decomp
            real(rkind),                     intent(in)    :: tsim
            real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
            real(rkind), dimension(:,:,:,:), intent(in)    :: fields
            real(rkind), dimension(:,:,:,:), intent(inout) :: rhs, rhsg

        end subroutine
    end interface

end module 
