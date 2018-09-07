module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi
    use decomp_2d

    implicit none 

    integer :: nt_print2screen = 1
    integer :: tid_statsDump = 5000
    integer :: tid_compStats = 100
    real(rkind) :: time_startDumping = 10.0_rkind, maxDiv, DomMaxDiv
    integer :: ierr 
    
    integer :: tid_start_planes = 1
    integer :: tid_stop_planes = 100000
    integer :: tid_dump_plane_every = 10000

contains

    subroutine doTemporalStuff(gp)
        class(igrid), intent(inout) :: gp 
      
        if (mod(gp%step,nt_print2screen) == 0) then
            maxDiv = maxval(gp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            call message(0,"Time",gp%tsim)
            if (gp%useSGS) call message(1,"u_star:",gp%sgsmodel%get_ustar())
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message_min_max(1,"Bounds for u:", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))
            if (gp%initspinup) then
               call message_min_max(1,"Bounds for T:", p_minval(minval(gp%T)), p_maxval(maxval(gp%T)))
            end if
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if
            call message("==========================================================")
            call toc()
            call tic()
        end if 

    end subroutine


end module 
