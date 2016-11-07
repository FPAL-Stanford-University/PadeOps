module sgrid_hooks
    use kind_parameters, only: rkind
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use SolidMixtureMod, only: solid_mixture
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
        subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tstop,dt,tviz)
            import :: rkind
            import :: decomp_info
            import :: solid_mixture
            type(decomp_info),               intent(in)    :: decomp
            real(rkind),                     intent(in)    :: dx, dy, dz
            character(len=*),                intent(in)    :: inputfile
            real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
            real(rkind), dimension(:,:,:,:), intent(inout) :: fields
            type(solid_mixture),             intent(inout) :: mix
            real(rkind),                     intent(inout) :: tstop, dt, tviz
        end subroutine 
    end interface


    interface hook_output
        ! subroutine hook_output(decomp,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
        subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount,x_bc,y_bc,z_bc)
            import :: rkind
            import :: decomp_info
            import :: derivatives
            import :: solid_mixture
            character(len=*),                intent(in) :: outputdir
            type(decomp_info),               intent(in) :: decomp
            type(derivatives),               intent(in) :: der
            type(solid_mixture),             intent(in) :: mix
            real(rkind),                     intent(in) :: dx,dy,dz,tsim
            integer,                         intent(in) :: vizcount
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(in) :: fields
            integer, dimension(2),           intent(in) :: x_bc,y_bc,z_bc
        end subroutine
    end interface

    interface hook_postproc
        subroutine hook_postproc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
            import :: rkind
            import :: decomp_info
            import :: solid_mixture
            type(decomp_info),               intent(in)    :: decomp
            real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
            real(rkind), dimension(:,:,:,:), intent(inout) :: fields
            type(solid_mixture),             intent(inout) :: mix
            real(rkind),                     intent(in)    :: tsim
            integer, dimension(2),           intent(in)    :: x_bc,y_bc,z_bc
        end subroutine
    end interface

    interface hook_bc
        subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc)
            import :: rkind
            import :: decomp_info
            import :: solid_mixture
            type(decomp_info),               intent(in)    :: decomp
            real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
            real(rkind), dimension(:,:,:,:), intent(inout) :: fields
            type(solid_mixture),             intent(inout) :: mix
            real(rkind),                     intent(in)    :: tsim
            integer, dimension(2),           intent(in)    :: x_bc,y_bc,z_bc
        end subroutine
    end interface

    interface hook_timestep
        subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim)
            import :: rkind
            import :: decomp_info
            import :: solid_mixture
            type(decomp_info),               intent(in) :: decomp
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(in) :: fields
            type(solid_mixture),             intent(in) :: mix
            integer,                         intent(in) :: step
            real(rkind),                     intent(in) :: tsim
        end subroutine
    end interface

    interface hook_mixture_source
        subroutine hook_mixture_source(decomp,mesh,fields,mix,tsim,rhs)
            import :: rkind
            import :: decomp_info
            import :: solid_mixture
            type(decomp_info),               intent(in)    :: decomp
            type(solid_mixture),             intent(in)    :: mix
            real(rkind),                     intent(in)    :: tsim
            real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
            real(rkind), dimension(:,:,:,:), intent(in)    :: fields
            real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
        end subroutine
    end interface

end module 
