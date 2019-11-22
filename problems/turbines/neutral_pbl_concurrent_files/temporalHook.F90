module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi

    implicit none 

    integer :: nt_print2screen = 1
    real(rkind) :: maxDiv, DomMaxDiv, time0
    integer :: ierr, stepper=0 
    
contains

    subroutine doTemporalStuff(gp, simid)
        class(igrid), intent(inout) :: gp 
        integer, intent(in) :: simid

        if (mod(gp%step,nt_print2screen) == 0) then
            maxDiv = maxval(gp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            select case (simid)
            case (1)
               call message(0,"Primary Simulation Info:")
            case (2)
               call message(0,"Concurrent Simulation Info:")
            end select 
            call message(0,"Time",gp%tsim)
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message(1,"u_star:",gp%sgsmodel%get_ustar())
            call message(1,"Inv. Ob. Len:",gp%sgsmodel%get_InvObLength())
            call message(1,"Surface Flux (K*nd velocity):",gp%wTh_surf)
            call message_min_max(1,"Bounds for u:", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))
            if (simid == 1) then
                call message(0,"Wind direction hub height", gp%WindTurbineArr%windAngle)
            end if 
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            call message(0,"------------------------------------------")
            if (simid == 1) then
               if (allocated(gp%scalars)) then
                  call message_min_max(1,"Bounds for SCALAR 1:", p_minval(minval(gp%scalars(1)%F)), p_maxval(maxval(gp%scalars(1)%F)))
                  call message_min_max(1,"Bounds for SCALAR 2:", p_minval(minval(gp%scalars(2)%F)), p_maxval(maxval(gp%scalars(2)%F)))
                  call message_min_max(1,"Bounds for SCALAR 3:", p_minval(minval(gp%scalars(3)%F)), p_maxval(maxval(gp%scalars(3)%F)))
               end if
            elseif (simid == 2) then
               call toc()
               call tic()
            end if 
        end if 
        

        if (stepper==0) then
            time0 = gp%Tsim
            stepper=1
        end if
        
        call update_surface_flux(gp%Tsim-time0, gp%wTh_surf)

    end subroutine

    subroutine update_surface_flux(time, surfaceFlux)
        real(rkind), intent(in) :: time
        real(rkind), intent(inout) :: surfaceFlux

        ! ADITYA TO MIKE: Set your function here
        ! Nocturnal surface flux slope Kumar et al. 2006, -0.003 (K*m/s) / hour
        ! -0.003 (K*m/s)/hour * (1 nd speed / G m/s) * (1hour/3600s) * (80s/T
        ! nd) = -1.3333e-05

        surfaceFlux = -1.3333D-5 * time

    end subroutine

end module 
