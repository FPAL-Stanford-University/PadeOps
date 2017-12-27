module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_maxval, p_minval
    use exits,              only: message, message_min_max, GracefulExit
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi

    implicit none 

    integer :: nt_print2screen = 1
    real(rkind) ::  maxDiv, DomMaxDiv, maxnusgs, maxkappasgs
    integer :: ierr 

contains

    subroutine doTemporalStuff(igp)
        class(igrid), intent(inout) :: igp 
      
        if (mod(igp%step,nt_print2screen) == 0) then
            maxDiv = maxval(igp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            call message(0,"Time",igp%tsim)
            call message(1,"TIDX:",igp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message_min_max(1,"Bounds for u:", p_minval(minval(igp%u)), p_maxval(maxval(igp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(igp%v)), p_maxval(maxval(igp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(igp%w)), p_maxval(maxval(igp%w)))
            call message_min_max(1,"Bounds for T:", p_minval(minval(igp%T)), p_maxval(maxval(igp%T)))
            if (igp%useSGS) then
                maxnusgs = p_maxval(igp%nu_SGS)
                maxkappasgs = p_maxval(igp%kappaSGS)
                call message(1,"Maximum SGS viscosity:", maxnusgs)
                call message(1,"Maximum SGS scalar kappa:", maxkappasgs)
               
                if (igp%sgsModel%usingDynProc()) then
                  call message(1,"Maximum lambda_dynamic:", igp%sgsModel%getMax_DynSmagConst())
                  call message(1,"Maximum beta_dynamic:", igp%sgsModel%getMax_DynPrandtl())
                end if 

            end if 
            if (igp%useCFL) then
                call message(1,"Current dt:",igp%dt)
            end if
            call message("==========================================================")
            call toc()
            call tic()
        end if 

    end subroutine


end module 
