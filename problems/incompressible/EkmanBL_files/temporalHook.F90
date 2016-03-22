module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGridNP, only: igrid
    use reductions,         only: P_MAXVAL
    use exits,              only: message
    use EkmanBL_IO,           only: dumpData4Matlab 
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi
    use AllStatistics 
    integer :: nt_print2screen = 10
    integer :: nt_getMaxKE = 10
    integer :: tid_statsDump = 100
    integer :: tid_compStats = 20
    real(rkind) :: time_startDumping = 2._rkind
    integer :: ierr 
contains

    subroutine doTemporalStuff(gp)
        class(igrid), intent(inout) :: gp 
      
        if (mod(gp%step,nt_print2screen) == 0) then
            call message(0,"Time",gp%tsim)
        end if

        if (mod(gp%step,nt_getMaxKE) == 0) then
            call message(1,"Max KE:",P_MAXVAL(gp%getMaxKE()))
            call toc()
            call tic()
        end if 

        if (mod(gp%step,gp%t_dataDump)==0) then
           call message("Data dump!")
           call dumpData4Matlab(gp) 
        end if 

        if ((mod(gp%step,tid_compStats)==0) .and. (gp%tsim > time_startDumping)) then
            call compute_stats(gp)
        end if 

        if ((mod(gp%step,tid_statsDump) == 0) .and. (gp%tsim > time_startDumping)) then
            call compute_stats(gp)
            call dump_stats(gp)
        end if 
        
        
        if (mod(gp%step,gp%t_restartDump) == 0) then
            call gp%dumpRestartfile()
        end if

        call mpi_barrier(mpi_comm_world,ierr)
        !if (mod(gp%step,gp%t_restartDump) == 0) then
        !    ! Incomplete 
        !end if 

    end subroutine


end module 
