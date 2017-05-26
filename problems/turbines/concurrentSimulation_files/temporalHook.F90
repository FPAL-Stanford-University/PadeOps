module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi

    implicit none 

    integer :: nt_print2screen = 20
    integer :: tid_statsDump = 5000
    integer :: tid_compStats = 100
    real(rkind) :: maxDiv, DomMaxDiv
    
    integer :: tid_start_planes = 1
    integer :: tid_stop_planes = 100000
    integer :: tid_dump_plane_every = 10000

contains

    subroutine doTemporalStuff(prec, gp)
        class(igrid), intent(inout) :: prec, gp 

        if (mod(gp%step,nt_print2screen) == 0) then
            maxDiv = maxval(gp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            call message(0,"Time",gp%tsim)
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message(1,"Prec u_star:",prec%sgsmodel%get_ustar())
            call message_min_max(1,"Bounds for u:", p_minval(minval(prec%u)), p_maxval(maxval(prec%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(prec%v)), p_maxval(maxval(prec%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(prec%w)), p_maxval(maxval(prec%w)))
            call message(1,"Main u_star:",gp%sgsmodel%get_ustar())
            call message_min_max(1,"Bounds for u:", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            call toc()
            call tic()
        end if 

    end subroutine


end module 
