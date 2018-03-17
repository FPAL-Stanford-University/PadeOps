module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max, GracefulExit
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi

    implicit none 

    integer :: nt_print2screen = 1
    real(rkind) ::  maxDiv, DomMaxDiv
    integer :: ierr 

contains

    subroutine doTemporalStuff(igp)
        class(igrid), intent(inout) :: igp 
      
        if (mod(igp%step,nt_print2screen) == 0) then
            maxDiv = maxval(igp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            call message(0,"Time",igp%tsim)
            call message(1,"u_star:",igp%sgsmodel%get_ustar())
            call message(1,"TIDX:",igp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message(1,"Inv. Ob. Len:",igp%sgsmodel%get_InvObLength())
            call message_min_max(1,"Bounds for u:", p_minval(minval(igp%u)), p_maxval(maxval(igp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(igp%v)), p_maxval(maxval(igp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(igp%w)), p_maxval(maxval(igp%w)))
            call message_min_max(1,"Bounds for T:", p_minval(minval(igp%T)), p_maxval(maxval(igp%T)))
            if (igp%useCFL) then
                call message(1,"Current dt:",igp%dt)
            end if
            call message("==========================================================")
            call toc()
            call tic()
            ! Added by MH, debugging blowup in neutral BL simulation
            !if (p_maxval(maxval(igp%u))>2.) then
            !    call message(1, "this step has blown up", igp%tsim)
            !    call igp%dumpFullField(igp%u,"uVel")
            !    call igp%dumpFullField(igp%v,"vVel")
            !    call igp%dumpFullField(igp%wC,"wVel")
            !    call igp%dumpFullField(igp%T, "potT")
            !    call GracefulExit("u-velocity has blown up",1)
            !end if
        end if 

    end subroutine


end module 
