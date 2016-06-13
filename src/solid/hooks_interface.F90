module sgrid_hooks
    use kind_parameters, only: rkind
    use decomp_2d, only: decomp_info
    use DerivativesMod, only: derivatives
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
        subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,rho0,mu,yield,gam,PInf,tau0,tstop,dt,tviz)
            import :: rkind
            import :: decomp_info
            type(decomp_info), intent(in) :: decomp
            real(rkind), intent(inout) :: dx, dy, dz
            character(len=*), intent(in) :: inputfile
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(inout) :: fields
            real(rkind),           optional, intent(inout) :: rho0, mu, gam, PInf, tstop, dt, tviz
            real(rkind),           optional, intent(inout) :: yield, tau0

        end subroutine 
    end interface


    interface initfields_stagg
        subroutine initfields_stagg(decompC, decompE, dx, dy, dz, inpDirectory, mesh, fieldsC, fieldsE, u_g, fcorr)
            import :: rkind
            import :: decomp_info
            type(decomp_info), intent(in) :: decompC
            type(decomp_info), intent(in) :: decompE
            real(rkind), intent(inout) :: dx, dy, dz
            character(len=*), intent(in) :: inpDirectory
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(inout) :: fieldsC
            real(rkind), dimension(:,:,:,:), intent(inout) :: fieldsE
            real(rkind), intent(out) :: u_g, fcorr

        end subroutine 
    end interface

    interface hook_output
        subroutine hook_output(decomp,dx,dy,dz,outputdir,mesh,fields,tsim,vizcount,der)
            import :: rkind
            import :: decomp_info
            import :: derivatives
            character(len=*),                intent(in) :: outputdir
            type(decomp_info),               intent(in) :: decomp
            real(rkind),                     intent(in) :: dx,dy,dz,tsim
            integer,                         intent(in) :: vizcount
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(in) :: fields
            type(derivatives),               intent(in) :: der

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
        subroutine hook_timestep(decomp,mesh,fields,step,tsim,hookcond)
            import :: rkind
            import :: decomp_info
            type(decomp_info),               intent(in)    :: decomp
            integer,                         intent(in)    :: step
            real(rkind),                     intent(in)    :: tsim
            logical,               optional, intent(inout) :: hookcond
            real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
            real(rkind), dimension(:,:,:,:), intent(in)    :: fields

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
            real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
            real(rkind), dimension(:,:,:,:), optional, intent(inout) ::rhsg

        end subroutine
    end interface

end module 
