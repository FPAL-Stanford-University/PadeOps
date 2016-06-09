module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGridWallM, only: igridWallM
    use reductions,         only: P_MAXVAL
    use exits,              only: message
    use gabls_IO,           only: dumpData4Matlab 
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi

    implicit none 

    integer :: nt_print2screen = 20
    integer :: tid_statsDump = 5000
    integer :: tid_compStats = 100
    real(rkind) :: time_startDumping = 10.0_rkind
    integer :: ierr 
    
    integer :: tid_start_planes = 1
    integer :: tid_stop_planes = 100000
    integer :: tid_dump_plane_every = 100

contains

    subroutine doTemporalStuff(gp)
        class(igridWallM), intent(inout) :: gp 
      
        if (mod(gp%step,nt_print2screen) == 0) then
            call message(0,"Time",gp%tsim)
            call message(1,"Max KE:",gp%getMaxKE())
            call message(1,"Max nuSGS:",gp%max_nuSGS)
            call message(1,"u_star:",gp%ustar)
            call message(1,"Inv. Ob. Length:",gp%InvObLength)
            call message(1,"Max u:", p_maxval(maxval(gp%u)))
            call message(1,"Max v:", p_maxval(maxval(gp%v)))
            call message(1,"Max w:", p_maxval(maxval(gp%w)))
            call message(1,"Max T:", p_maxval(maxval(gp%T)))
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            if ((gp%useDynamicProcedure) .and. (gp%useSGS)) then
                call message(1,"Max cSGS:",p_maxval(maxval(gp%c_SGS(1,1,:))))
            end if 
            call toc()
            call tic()
        end if 

        if (mod(gp%step,gp%t_dataDump)==0) then
           call message(0,"Data dump!")
           call gp%dumpFullField(gp%u,'uVel')
           call gp%dumpFullField(gp%v,'vVel')
           call gp%dumpFullField(gp%wC,'wVel')
           call gp%dumpFullField(gp%T,'potT')
        end if 

    end subroutine


end module 
