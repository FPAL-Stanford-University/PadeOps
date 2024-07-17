module cgrid_hooks
    use kind_parameters, only: rkind
    use decomp_2d,       only: decomp_info
    use DerivativesMod,  only: derivatives
    use MixtureEOSMod,   only: mixture
    use sgsmod_cgrid,    only: sgs_cgrid
    implicit none

    interface meshgen
        subroutine meshgen(decomp, dx, dy, dz, mesh, inputfile, xmetric, ymetric, zmetric, xi, eta, zeta, dxs, dys, dzs, xbuf, zbuf)
            import :: rkind
            import :: decomp_info
            type(decomp_info),               intent(in)    :: decomp
            real(rkind),                     intent(inout) :: dx, dy, dz
            real(rkind), dimension(:,:,:,:), intent(inout) :: mesh
            character(len=*),                intent(in)    :: inputfile
            logical    ,                     intent(in)    :: xmetric, ymetric, zmetric
            real(rkind), dimension(:,:,:),   intent(inout) :: xi, eta, zeta
            real(rkind), dimension(:,:,:),   intent(inout) :: dxs, dys, dzs
            real(rkind), dimension(:,:,:,:), intent(in)    :: xbuf, zbuf

        end subroutine 
    end interface

    interface initfields
        subroutine initfields(decomp,dx,dy,dz,inputfile,mesh,fields,mix,tsim,tstop,dt,tviz)
            import :: rkind
            import :: decomp_info
            import :: mixture
            type(decomp_info),               intent(in)    :: decomp
            type(mixture),                   intent(inout) :: mix
            real(rkind),                     intent(in)    :: dx, dy, dz
            character(len=*),                intent(in)    :: inputfile
            real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
            real(rkind), dimension(:,:,:,:), intent(inout) :: fields
            real(rkind),                     intent(inout) :: tsim, tstop, dt, tviz

        end subroutine 
    end interface


    interface hook_output
        subroutine hook_output(decomp,der,dx,dy,dz,outputdir,mesh,fields,mix,tsim,vizcount)
            import :: rkind
            import :: decomp_info
            import :: derivatives
            import :: mixture
            character(len=*),                intent(in) :: outputdir
            type(decomp_info),               intent(in) :: decomp
            type(mixture),                   intent(in) :: mix
            type(derivatives),               intent(in) :: der
            real(rkind),                     intent(in) :: dx,dy,dz,tsim
            integer,                         intent(in) :: vizcount
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(in) :: fields

        end subroutine
    end interface

    interface hook_bc
        subroutine hook_bc(decomp,mesh,fields,mix,tsim,x_bc,y_bc,z_bc,newTimeStep, time_step)
            import :: rkind
            import :: decomp_info
            import :: mixture
            type(decomp_info),               intent(in)    :: decomp
            type(mixture),                   intent(in)    :: mix
            real(rkind),                     intent(in)    :: tsim
            real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
            real(rkind), dimension(:,:,:,:), intent(inout) :: fields
            integer, dimension(2),           intent(in)    :: x_bc, y_bc, z_bc
            logical,                         intent(in)    :: newTimeStep
            integer,                         intent(in)    :: time_step

        end subroutine
    end interface

    interface hook_timestep
        subroutine hook_timestep(decomp,mesh,fields,mix,step,tsim,sgsmodel)
            import :: rkind
            import :: decomp_info
            import :: mixture
            import :: sgs_cgrid
            type(decomp_info),               intent(in) :: decomp
            type(mixture),                   intent(in) :: mix
            integer,                         intent(in) :: step
            real(rkind),                     intent(in) :: tsim
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(in) :: fields
            type(sgs_cgrid), optional,       intent(in) :: sgsmodel

        end subroutine
    end interface

    interface hook_source
        !subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs,Wcnsrv,tkeb,tsim_0,dtheta_0,der,dt)
        subroutine hook_source(decomp,mesh,fields,mix,tsim,rhs,der,dt,step,dys)
            use kind_parameters, only: rkind
            use decomp_2d,       only: decomp_info
            use MixtureEOSMod,   only: mixture
            use TKEBudgetMod,    only: tkeBudget
            use DerivativesMod,  only: derivatives

            implicit none
            type(decomp_info),               intent(in)    :: decomp
            type(mixture),                   intent(in)    :: mix
            real(rkind),                     intent(in)    :: tsim
            real(rkind), dimension(:,:,:,:), intent(in)    :: mesh
            real(rkind), dimension(:,:,:,:), intent(in)    :: fields
            real(rkind), dimension(:,:,:,:), intent(inout) :: rhs
            ! KVM 2021
            !real(rkind), dimension(:,:,:,:), optional,intent(in)    :: Wcnsrv
            !type(tkeBudget), optional,       intent(inout) :: tkeb
            !real(rkind), optional,           intent(in)    :: tsim_0 ! the previous time 
            !real(rkind), optional,           intent(in)    :: dtheta_0 ! the previous L99
            ! channel flow, vishwaja 2024 
            type(derivatives), optional,     intent(in)    :: der
            real(rkind), optional,           intent(in)    :: dt
            integer, optional,               intent(in)    :: step
            real(rkind),dimension(:,:,:),optional,  intent(in) :: dys
        end subroutine
    end interface

end module 
