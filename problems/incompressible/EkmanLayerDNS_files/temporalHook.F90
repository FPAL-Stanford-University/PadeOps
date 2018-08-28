module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max
    !use EkmanDNS_IO,           only: output_tecplot!dumpData4Matlab 
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi

    implicit none 

    integer :: nt_print2screen = 1
    real(rkind) :: maxDiv, DomMaxDiv, maxnusgs
    integer :: ierr 
   
    real(rkind), dimension(:), allocatable :: tmp
    real(rkind), dimension(:,:), allocatable :: ubudget, vbudget, AllBudgets


    integer :: compute_frequency = 9999999, write_frequency = 999999999
    logical :: dobudgetcalcs = .false.
    integer :: NumEnsembles

contains
    

    subroutine doTemporalStuff(gp)
        class(igrid), intent(inout) :: gp 
      
        if (mod(gp%step,nt_print2screen) == 0) then
            maxDiv = p_maxval(maxval(gp%divergence))
            call message(0,"Time",gp%tsim)
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",maxDiv)
            call message_min_max(1,"Bounds for u:", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))
            if (gp%initspinup .or. gp%isStratified) then
                call message_min_max(1,"Bounds for T:", p_minval(minval(gp%T)), p_maxval(maxval(gp%T)))
            end if 
            if (gp%useSGS) then
                maxnusgs = p_maxval(gp%nu_SGS)
                call message(1,"Maximum SGS viscosity:", maxnusgs)
                if (gp%sgsModel%usingDynProc()) then
                  call message(1,"Maximum lambda_dynamic:", gp%sgsModel%getMax_DynSmagConst())
                end if 
            end if 
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            call toc()
            call tic()
        end if 

        if (mod(gp%step,compute_frequency) == 0) then
            call process_stats(gp)
        end if 
        
        if (mod(gp%step,write_frequency) == 0) then
            call dump_problem_stats(gp)
        end if 

    end subroutine

    subroutine initialize_processing(nz,inputfile)
        integer, intent(in) :: nz
        character(len=*), intent(in) :: inputfile
        integer :: ioUnit
        namelist /ProcessingStats/doBudgetcalcs, write_frequency, compute_frequency

        ioUnit = 11
        open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
        read(unit=ioUnit, NML=ProcessingStats)
        close(ioUnit)    

        if (dobudgetcalcs) then
            allocate(tmp(nz))
            allocate(ubudget(nz,4))
            allocate(vbudget(nz,4))
            
            allocate(AllBudgets(nz,8))

            ubudget = 0.d0 
            vbudget = 0.d0 
            numEnsembles = 0
        end if 

    end subroutine 

    subroutine dump_problem_stats(gp)
        use kind_parameters, only: clen 
        use basic_io, only: write_2d_ascii
        class(igrid), intent(in) :: gp 
        character(len=clen) :: tempname, fname

        if (dobudgetcalcs) then
           AllBudgets(:,1:4) = ubudget/NumEnsembles 
           AllBudgets(:,5:8) = vbudget/NumEnsembles 

         write(tempname,"(A3,I2.2,A15,A2,I6.6,A2,I6.6,A4)") "Run",gp%runID,"_MomentumBudget","_t",gp%step,"_n",NumEnsembles,".stt"
         fname = gp%OutputDir(:len_trim(gp%OutputDir))//"/"//trim(tempname)
         call write_2d_ascii(AllBudgets, fname)
        end if 

    end subroutine 

    subroutine process_stats(gp)
        class(igrid), intent(inout) :: gp 
        
        if (dobudgetcalcs) then

        end if 

    end subroutine 

    subroutine finalize_processing()
        if (dobudgetcalcs) deallocate(tmp, ubudget, vbudget)
    end subroutine 
end module 
