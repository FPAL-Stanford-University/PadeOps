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
    integer :: tid_statsDump = 5000
    integer :: tid_compStats = 100
    real(rkind) :: maxDiv, DomMaxDiv
    
    integer :: tid_start_planes = 1
    integer :: tid_stop_planes = 100000
    integer :: tid_dump_plane_every = 10000

contains

    subroutine doTemporalStuff(gp, simid)
        class(igrid), intent(inout) :: gp 
        integer, intent(in) :: simid

        if (mod(gp%step,nt_print2screen) == 0) then
            maxDiv = maxval(gp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            select case (simid)
            case (1)
               call message(0,"Actuator Disk Simulation Info:")
            case (2)
               call message(0,"HIT Simulation Info:")
            end select 
            call message(0,"Time",gp%tsim)
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message_min_max(1,"Bounds for u:", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))
            call message_min_max(1,"Bounds for T:", p_minval(minval(gp%T)), p_maxval(maxval(gp%T)))
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            call message(0,"------------------------------------------")
            if (simid == 1) then
               if (allocated(gp%scalars)) then
                  call message_min_max(1,"Bounds for SCALAR 1:", p_minval(minval(gp%scalars(1)%F)), p_maxval(maxval(gp%scalars(1)%F)))
                  call message_min_max(1,"Bounds for SCALAR 2:", p_minval(minval(gp%scalars(2)%F)), p_maxval(maxval(gp%scalars(2)%F)))
               end if
            elseif (simid == 2) then
               call message(1,"Mean TKE for HIT:", gp%getMeanKE())
               call message(1,"Mean  uu for HIT:", gp%getMeanuu())
               call message(1,"Mean  uv for HIT:", gp%getMeanuv())
               call message(1,"Mean  uw for HIT:", gp%getMeanuw())
               call message(1,"Mean  vv for HIT:", gp%getMeanvv())
               call message(1,"Mean  vw for HIT:", gp%getMeanvw())
               call message(1,"Mean  ww for HIT:", gp%getMeanww())
               call toc()
               call tic()
            end if 
        end if 

    end subroutine


end module 
