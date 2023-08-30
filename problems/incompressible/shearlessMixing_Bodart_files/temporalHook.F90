module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval, p_sum
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

    subroutine doTemporalStuff(gp)
        class(igrid), intent(inout) :: gp

        if (mod(gp%step,nt_print2screen) == 0) then
            maxDiv = maxval(gp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            call message(0,"Time",gp%tsim)
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message_min_max(1,"Bounds for u:" , p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v:" , p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w:" , p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))
            call message_min_max(1,"Bounds for T:" , p_minval(minval(gp%T)), p_maxval(maxval(gp%T)))
            if (gp%localizedForceLayer == 1) then
                call message(1,"Forcing ampFact:",gp%forcelayer%ampFact)
                call message_min_max(1,"Bounds for fx:", &
                  gp%forcelayer%ampFact*p_minval(minval(gp%forcelayer%fx)), &
                  gp%forcelayer%ampFact*p_maxval(maxval(gp%forcelayer%fx)))
                call message_min_max(1,"Bounds for fy:", &
                  gp%forcelayer%ampFact*p_minval(minval(gp%forcelayer%fy)), &
                  gp%forcelayer%ampFact*p_maxval(maxval(gp%forcelayer%fy)))
                call message_min_max(1,"Bounds for fz:", &
                  gp%forcelayer%ampFact*p_minval(minval(gp%forcelayer%fz)), &
                  gp%forcelayer%ampFact*p_maxval(maxval(gp%forcelayer%fz)))
            elseif (gp%localizedForceLayer == 2) then
                call message(1,"Forcing layer divergence",gp%spectForceLayer%maxDiv)
                call message(1,"Forcing layer max all-time divergence",gp%spectForceLayer%maxDivAllTime)
                call message_min_max(1,"Bounds for fx:", &
                  gp%spectForceLayer%ampFact*p_minval(minval(gp%spectForceLayer%fx)), &
                  gp%spectForceLayer%ampFact*p_maxval(maxval(gp%spectForceLayer%fx)))
                call message_min_max(1,"Bounds for fy:", &
                  gp%spectForceLayer%ampFact*p_minval(minval(gp%spectForceLayer%fy)), &
                  gp%spectForceLayer%ampFact*p_maxval(maxval(gp%spectForceLayer%fy)))
                call message_min_max(1,"Bounds for fz:", &
                  gp%spectForceLayer%ampFact*p_minval(minval(gp%spectForceLayer%fz)), &
                  gp%spectForceLayer%ampFact*p_maxval(maxval(gp%spectForceLayer%fz)))

            end if
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            call message(0,"------------------------------------------")
            if (allocated(gp%scalars)) then
               call message_min_max(1,"Bounds for SCALAR 1:", p_minval(minval(gp%scalars(1)%F)), p_maxval(maxval(gp%scalars(1)%F)))
               call message_min_max(1,"Bounds for SCALAR 2:", p_minval(minval(gp%scalars(2)%F)), p_maxval(maxval(gp%scalars(2)%F)))
            end if
        end if 

    end subroutine


end module 
