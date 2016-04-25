module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGridWallM, only: igridWallM
    use reductions,         only: P_MAXVAL
    use exits,              only: message,GracefulExit
    use pbl_IO,           only: dumpData4Matlab 
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi

    implicit none 

    integer :: nt_print2screen = 20
    integer :: nt_getMaxKE = 20
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
        end if

        if (mod(gp%step,nt_getMaxKE) == 0) then
            call message(1,"Max KE:",gp%getMaxKE())
            call message(1,"Max nuSGS:",gp%max_nuSGS)
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            if ((gp%useDynamicProcedure) .and. (gp%useSGS)) then
                call message(1,"Max cSGS:",p_maxval(maxval(gp%c_SGS(1,1,:))))
            end if 
            !call message(1,"Mean ustar:",gp%mean_ustar)
            !call message(1,"Max ustar:",gp%max_ustar)
            !call message(1,"Min ustar:",gp%min_ustar)
            call toc()
            call tic()
        end if 

        if (mod(gp%step,100*nt_print2screen)==0) then
           call message("Dumped Planes")
           call gp%dump_planes()
        end if 

        if (mod(gp%step,gp%t_dataDump)==0) then
           call message("Data dump!")
           !call dumpData4Matlab(gp) 
           call gp%dumpFullField(gp%u,'uVel')
           call gp%dumpFullField(gp%v,'vVel')
           call gp%dumpFullField(gp%wC,'wVel')
           call gp%dumpFullField(gp%nu_SGS,'nuSG')
        end if 

        if ((mod(gp%step,tid_compStats)==0) .and. (gp%tsim > time_startDumping)) then
            call gp%compute_stats()
        end if 

        if ((mod(gp%step,tid_statsDump) == 0) .and. (gp%tsim > time_startDumping)) then
            call gp%compute_stats()
            call gp%dump_stats()
        end if 
        
        if (mod(gp%step,gp%t_restartDump) == 0) then
            call gp%dumpRestartfile()
        end if
        
        if ((gp%dumpPlanes) .and. (mod(gp%step,tid_dump_plane_every) == 0) .and. &
                 (gp%step .ge. tid_start_planes) .and. (gp%step .le. tid_stop_planes)) then
            call gp%dump_planes()
        end if 

        open(unit=20,file='exit',status='old',iostat=ierr)
        if(ierr==0) call GracefulExit("exit file found in working directory. Stopping run.",123)

    end subroutine


end module 
