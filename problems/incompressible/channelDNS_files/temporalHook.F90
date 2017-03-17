module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max
    !use channelDNS_IO,           only: output_tecplot!dumpData4Matlab 
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
        class(igrid), intent(inout) :: gp 
        real(rkind) :: utau, gdcoeff

        if (mod(gp%step,nt_print2screen) == 0) then
            maxDiv = maxval(gp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            call message(0,"Time",gp%tsim)
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message_min_max(1,"Bounds for u:", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            ! compute utau
            utau = sqrt(2.0d0*sqrt(real(gp%uhat(1,1,1),rkind)**2 + real(gp%vhat(1,1,1),rkind)**2)*gp%meanFact/(gp%dz*gp%Re))
            call message(1,"utau_lo:", utau)
            if (gp%useSGS) then
               if ((gp%sgsModel%DynamicProcedureType==2)) then
                  gdcoeff = gp%sgsModel%cmodel_global/(1.5d0*gp%dx*1.5d0*gp%dy*gp%dz)**(1.d0/3.d0)
                  call message(1,"Glob Dyn Coeff:",gdcoeff)
               end if
            end if 
            call toc()
            call tic()
        end if 

    end subroutine


end module 
