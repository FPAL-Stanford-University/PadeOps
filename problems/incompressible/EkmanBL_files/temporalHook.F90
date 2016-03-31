module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGridNP, only: igrid
    use reductions,         only: P_MAXVAL
    use exits,              only: message
    use EkmanBL_IO,           only: dumpData4Matlab 
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi
    integer :: nt_print2screen = 20
    integer :: nt_getMaxKE = 20
    integer :: tid_statsDump = 2000
    integer :: tid_compStats = 20
    real(rkind) :: time_startDumping = 75._rkind
    integer :: ierr 
contains

    subroutine doTemporalStuff(gp)
        class(igrid), intent(inout) :: gp 
      
        if (mod(gp%step,nt_print2screen) == 0) then
            call message(0,"Time",gp%tsim)
        end if

        if (mod(gp%step,nt_getMaxKE) == 0) then
            call message(1,"Max KE:",gp%getMaxKE())
            call message(1,"Max nuSGS:",gp%max_nuSGS)
            call message(1,"Max cSGS:",p_maxval(maxval(gp%c_SGS(1,1,:))))
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

    end subroutine


end module 
