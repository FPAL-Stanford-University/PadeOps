module igrid_hooks
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
    
    interface meshgen_WallM
        subroutine meshgen_WallM(decomp, dx, dy, dz, mesh, inputfile)
            import :: rkind
            import :: decomp_info
            type(decomp_info), intent(in) :: decomp
            real(rkind), intent(inout) :: dx, dy, dz
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            character(len=*), intent(in) :: inputfile

        end subroutine 
    end interface 

    interface initfields
        subroutine initfields(decomp, dx, dy, dz, inpDirectory, mesh, fields, rho0, mu, gam, PInf, tstop, dt, tviz)
            import :: rkind
            import :: decomp_info
            type(decomp_info), intent(in) :: decomp
            real(rkind), intent(inout) :: dx, dy, dz
            character(len=*), intent(in) :: inpDirectory
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(inout) :: fields
            real(rkind),           optional, intent(inout) :: rho0, mu, gam, PInf, tstop, dt, tviz

        end subroutine 
    end interface

    interface getForcing
        subroutine getForcing(inputfile, dpdx)
            import :: rkind
            real(rkind), intent(out) :: dpdx 
            character(len=*), intent(in) :: inputfile
        end subroutine
    end interface

    interface initfields_stagg
        subroutine initfields_stagg(decompC, decompE, inpDirectory, mesh, fieldsC, fieldsE, u_g, Ro)
            import :: rkind
            import :: decomp_info
            type(decomp_info), intent(in) :: decompC
            type(decomp_info), intent(in) :: decompE
            character(len=*), intent(in) :: inpDirectory
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(inout) :: fieldsC
            real(rkind), dimension(:,:,:,:), intent(inout) :: fieldsE
            real(rkind), intent(out) :: u_g, Ro

        end subroutine 
    end interface


   interface initScalar
      subroutine initScalar(decompC, inpDirectory, mesh, scalar_id, scalarField)
            import :: rkind
            import :: decomp_info
            type(decomp_info), intent(in) :: decompC
            character(len=*), intent(in) :: inpDirectory
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:), intent(out) :: scalarField
            integer, intent(in) :: scalar_id
      end subroutine 
   end interface


   interface SetScalar_source
      subroutine SetScalar_source(decompC, inpDirectory, mesh, scalar_id, scalarSource)
            import :: rkind
            import :: decomp_info
            type(decomp_info), intent(in) :: decompC
            character(len=*), intent(in) :: inpDirectory
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:), intent(out) :: scalarSource
            integer, intent(in) :: scalar_id
      end subroutine 
   end interface


    interface setInhomogeneousNeumannBC_Temp
        subroutine setInhomogeneousNeumannBC_Temp(inpDirectory, wTh_surf)
            import :: rkind
            character(len=*), intent(in) :: inpDirectory
            real(rkind), intent(out) :: wTh_surf

        end subroutine 
    end interface

    interface setDirichletBC_Temp
        subroutine setDirichletBC_Temp(inpDirectory, Tsurf, dTsurfdt)
            import :: rkind
            character(len=*), intent(in) :: inpDirectory
            real(rkind), intent(out) :: Tsurf, dTsurfdt

        end subroutine 
    end interface

    interface set_Reference_Temperature
        subroutine set_Reference_Temperature(inputfile, Tref)
            import :: rkind
            character(len=*), intent(in) :: inputfile
            real(rkind), intent(out) :: Tref

        end subroutine 
    end interface

    interface initfields_wallM
        subroutine initfields_wallM(decompC, decompE, inpDirectory, mesh, fieldsC, fieldsE)
            import :: rkind
            import :: decomp_info
            type(decomp_info), intent(in) :: decompC
            type(decomp_info), intent(in) :: decompE
            character(len=*), intent(in) :: inpDirectory
            real(rkind), dimension(:,:,:,:), intent(in) :: mesh
            real(rkind), dimension(:,:,:,:), intent(inout) :: fieldsC
            real(rkind), dimension(:,:,:,:), intent(inout) :: fieldsE

        end subroutine 
    end interface
    
    interface set_planes_io
        subroutine set_planes_io(xplanes, yplanes, zplanes)
            integer, dimension(:), allocatable,  intent(inout) :: xplanes
            integer, dimension(:), allocatable,  intent(inout) :: yplanes
            integer, dimension(:), allocatable,  intent(inout) :: zplanes
        end subroutine
    end interface

    interface set_KS_planes_io
        subroutine set_KS_planes_io(planesCourseGrid, planesFineGrid)
            integer, dimension(:), allocatable,  intent(inout) :: planesFineGrid
            integer, dimension(:), allocatable,  intent(inout) :: planesCourseGrid
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

    interface hook_probes
        subroutine hook_probes(inputfile, probe_locs)
            import :: rkind
            character(len=*),                intent(in)    :: inputfile
            real(rkind), dimension(:,:), allocatable, intent(inout) :: probe_locs

        end subroutine
    end interface
end module 
