module temporalHook
    use kind_parameters,    only: rkind
    !use IncompressibleGridWallM, only: igridWallM
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max
    !use pblwt_IO,           only: output_tecplot!dumpData4Matlab 
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi

    implicit none 

    integer :: nt_print2screen = 20
    integer :: tid_statsDump = 5000
    integer :: tid_compStats = 100
    real(rkind) :: time_startDumping = 10.0_rkind, maxDiv, DomMaxDiv
    integer :: ierr 
    
    integer :: tid_start_planes = 1
    integer :: tid_stop_planes = 100000
    integer :: tid_dump_plane_every = 10000

contains

    subroutine doTemporalStuff(gp)
        !class(igridWallM), intent(inout) :: gp 
        class(igrid), intent(inout) :: gp 
        real(rkind) :: ubulk_f(3), ubulk_s(3)

        if (mod(gp%step,nt_print2screen) == 0) then
            maxDiv = maxval(gp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            call message(0,"Time",gp%tsim)
            call message(1,"u_star:",gp%sgsmodel%get_ustar())
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message(1,"ubulk:",gp%ubulk)
            call message(1,"dpFdx:",gp%dpFdx)
            call message_min_max(1,"Bounds for u:", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            !call message(1,gp%dtlimit)
            if(gp%useIBM) then
                ! from ghost points
                call message(1,         "Average of utauIB-Gpts:", gp%ibm%get_ibwm_utau_avg())
                call message_min_max(1, "Bounds for utauIB-Gpts:", gp%ibm%get_ibwm_utau_min(), gp%ibm%get_ibwm_utau_max())

                ! from stressC points
                call message(1,         "Average of utauIB StrC:", gp%ibm%get_ibwm_utauStrsC_avg())
                call message_min_max(1, "Bounds for utauIB StrC:", gp%ibm%get_ibwm_utauStrsC_min(), gp%ibm%get_ibwm_utauStrsC_max())

                ! from stressE points
                call message(1,         "Average of utauIB StrE:", gp%ibm%get_ibwm_utauStrsE_avg())
                call message_min_max(1, "Bounds for utauIB StrE:", gp%ibm%get_ibwm_utauStrsE_min(), gp%ibm%get_ibwm_utauStrsE_max())

                ! from stressE points
                call gp%ibm%get_ubulk_f(ubulk_f);          call gp%ibm%get_ubulk_s(ubulk_s);
                call message_min_max(1, "ubulk fluid-solid     :", ubulk_f(1), ubulk_s(1))
                call message_min_max(1, "vbulk fluid-solid     :", ubulk_f(2), ubulk_s(2))
                call message_min_max(1, "wbulk fluid-solid     :", ubulk_f(3), ubulk_s(3))
                call message(1,         "dpFdxs                :", gp%dpFdxs)
            endif
            call toc()
            call tic()
        end if 

    end subroutine

end module 
