module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi
    use HIT_Periodic_parameters, only: uTarget, vTarget, wTarget, x_shift, uadvect, useBandpassFilter

    implicit none 

    integer :: nt_print2screen = 1
    real(rkind) :: DomMaxDiv
    integer :: ierr 

contains

    subroutine doTemporalStuff(gp)
        class(igrid), intent(inout), target :: gp 
        real(rkind), dimension(:,:,:), pointer ::x, y, z
         
         z => gp%mesh(:,:,:,3)
         y => gp%mesh(:,:,:,2)
         x => gp%mesh(:,:,:,1)
   
         x_shift = x_shift + uadvect*gp%get_dt()

         if (mod(gp%step,gp%t_dataDump) == 0) then
            if (useBandpassFilter) then
               call message(0,"Shifting solution by", x_shift)
               call gp%spectC%bandpassFilter_and_PhaseShift(gp%uhat , uTarget, x_shift)
               call gp%spectC%bandpassFilter_and_PhaseShift(gp%vhat , vTarget, x_shift)
               call gp%spectC%bandpassFilter_and_PhaseShift(gp%whatC, wTarget, x_shift)
               call toc()
               call gp%dumpFullField(uTarget,'uBPF')
               call gp%dumpFullField(vTarget,'vBPF')
               call gp%dumpFullField(wTarget,'wBPF')
               call message(0, "Just dumped bandpass filtered fields")
            end if 
         end if 
        
         if (mod(gp%step,nt_print2screen) == 0) then
            DomMaxDiv = p_maxval(gp%divergence)
            call message(0,"Time",gp%tsim)
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message_min_max(1,"Bounds for u:", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))
            if (allocated(gp%pressure)) then
               call message_min_max(1,"Bounds for P:", p_minval(minval(gp%pressure)), p_maxval(maxval(gp%pressure)))
            end if
            if (gp%useSGS) then
                call message_min_max(1,"Bounds for nuSGS:", p_minval(minval(gp%nu_sgs)), p_maxval(maxval(gp%nu_sgs)))
                if (gp%sgsModel%usingDynProc()) then
                    call message(2,"Model multiplier, (c*delta)^2:",gp%sgsModel%getMax_DynSmagConst())
                end if 
            end if

            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
           call message(1,"Mean TKE for HIT:",gp%getMeanKE())
            call toc()
            call tic()
        end if 

    end subroutine

end module 
