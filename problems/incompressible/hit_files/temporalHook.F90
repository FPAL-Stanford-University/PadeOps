module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL
    use exits,              only: message
    use hitCD_IO,           only: dumpData4Matlab 
    use constants,          only: half 
    integer :: nt_print2screen = 20
    integer :: nt_getMaxKE = 20
contains

    subroutine doTemporalStuff(gp)
        class(igrid), intent(in) :: gp 
      
        if (mod(gp%step,nt_print2screen) == 0) then
            call message(0,"Time",gp%tsim)
        end if

        if (mod(gp%step,nt_getMaxKE) == 0) then
            call message(1,"Max KE:",P_MAXVAL(half*(gp%u**2 + gp%v**2 + gp%w**2)))
        end if 

        if (mod(gp%step,gp%t_dataDump)==0) then
            call message("Data dump!")
            call dumpData4Matlab(gp) 
        end if 

        if (mod(gp%step,gp%t_restartDump) == 0) then
            ! Incomplete 
        end if 

    end subroutine


end module 
