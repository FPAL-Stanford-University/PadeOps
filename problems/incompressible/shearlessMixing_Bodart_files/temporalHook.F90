module temporalHook
    use kind_parameters,    only: rkind, clen
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval, p_sum
    use exits,              only: message, message_min_max
    use constants,          only: half
    use mpi
    use fortran_assert,     only: assert

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
        integer :: n
        character(len=clen) :: mssg
        real(rkind) :: umax, umin

        if (mod(gp%step,nt_print2screen) == 0) then
            maxDiv = maxval(gp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            call message(0,"Time",gp%tsim)
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            umax = p_maxval(maxval(gp%u))
            umin = p_minval(minval(gp%u))
            if (umax > 100.d0 .or. umin < -100.d0) then
                call assert(.false.,'Simulation is blowing up')
            end if
            call message_min_max(1,"Bounds for u:" , umin, umax)
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
                  gp%spectForceLayer%ampFact_x*p_minval(minval(gp%spectForceLayer%fx)), &
                  gp%spectForceLayer%ampFact_x*p_maxval(maxval(gp%spectForceLayer%fx)))
                call message_min_max(1,"Bounds for fy:", &
                  gp%spectForceLayer%ampFact_y*p_minval(minval(gp%spectForceLayer%fy)), &
                  gp%spectForceLayer%ampFact_y*p_maxval(maxval(gp%spectForceLayer%fy)))
                call message_min_max(1,"Bounds for fz:", &
                  gp%spectForceLayer%ampFact_z*p_minval(minval(gp%spectForceLayer%fz)), &
                  gp%spectForceLayer%ampFact_z*p_maxval(maxval(gp%spectForceLayer%fz)))
                call message(1,"Forcing throttling factor",gp%spectForceLayer%ampFact_time)

            end if
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            call message(0,"------------------------------------------")
            if (allocated(gp%scalars)) then
               do n = 1,gp%n_scalars
                 write(mssg,'(A18,I2.2,A1)')"Bounds for SCALAR ",n,":"
                 call message_min_max(1,trim(mssg), p_minval(minval(gp%scalars(n)%F)), p_maxval(maxval(gp%scalars(n)%F)))
               end do
            end if
        end if 

    end subroutine


end module 
