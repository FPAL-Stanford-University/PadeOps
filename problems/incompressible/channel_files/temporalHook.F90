module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGridNP, only: igrid
    use reductions,         only: P_MAXVAL
    use exits,              only: GracefulExit, message
    use channel_IO,           only: dumpData4Matlab 
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi
    integer :: nt_print2screen = 10
    integer :: nt_getMaxKE = 10
    integer :: tid_statsDump = 40000
    integer :: tid_compStats = 10
    real(rkind) :: time_startDumping = 2000.0_rkind
    integer :: ierr 
contains

    subroutine doTemporalStuff(gp)
        class(igrid), intent(inout) :: gp 
        integer :: ierr
      
        if (mod(gp%step,nt_print2screen) == 0) then
            call message(0,"Time",gp%tsim)
        end if

        if (mod(gp%step,nt_getMaxKE) == 0) then
            call message(1,"Max KE:",gp%getMaxKE())
            call message(1,"Max nuSGS:",gp%max_nuSGS)
            call toc()
            call tic()
        end if 

        if (mod(gp%step,gp%t_dataDump)==0) then
           call message("Data dump!")
           call dumpData4Matlab(gp) 
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

        open(unit=20,file='exit',status='old',iostat=ierr)
        if(ierr==0) call GracefulExit("exit file found in working directory. Stopping run.",123)

    end subroutine


end module 
