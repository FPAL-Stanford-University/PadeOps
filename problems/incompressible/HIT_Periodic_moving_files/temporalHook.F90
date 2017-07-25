module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi
    use HIT_Periodic_parameters

    implicit none 

    integer :: nt_print2screen = 1
    integer :: tid_statsDump = 5000
    integer :: tid_compStats = 100
    real(rkind) :: time_startDumping = 10.0_rkind, maxDiv, DomMaxDiv
    integer :: ierr 
    
    integer :: tid_start_planes = 1
    integer :: tid_stop_planes = 100000
    integer :: tid_dump_plane_every = 10000

contains

    subroutine doTemporalStuff(gp)
        class(igrid), intent(inout), target :: gp 
        real(rkind), dimension(:,:,:), pointer ::x, y, z
         
         z => gp%mesh(:,:,:,3)
         y => gp%mesh(:,:,:,2)
         x => gp%mesh(:,:,:,1)
   
        ! case (1)
        !    uexact  =  sin(x)*cos(y)*exp(-(2.d0/gp%Re)*gp%tsim)
        !    vexact  = -cos(x)*sin(y)*exp(-(2.d0/gp%Re)*gp%tsim)
        !    wexact  = 0.d0
        !    pexact  = 0.25d0*(cos(2.d0*x) + cos(2.d0*y))*((exp(-(2.d0/gp%Re)*gp%tsim))**2)
        ! case (2)
        !    uexact  =  sin(x)*cos(z)*exp(-(2.d0/gp%Re)*gp%tsim)
        !    vexact  = 0.d0
        !    wexact  = -cos(x)*sin(z)*exp(-(2.d0/gp%Re)*gp%tsim)
        !    pexact  = 0.25d0*(cos(2.d0*x) + cos(2.d0*z))*((exp(-(2.d0/gp%Re)*gp%tsim))**2)
        ! case (3)
        !    uexact  = 0.d0
        !    vexact  =  sin(y)*cos(z)*exp(-(2.d0/gp%Re)*gp%tsim)
        !    wexact  = -cos(y)*sin(z)*exp(-(2.d0/gp%Re)*gp%tsim)
        !    pexact  = 0.25d0*(cos(2.d0*y) + cos(2.d0*z))*((exp(-(2.d0/gp%Re)*gp%tsim))**2)
        ! end select
        
         if (mod(gp%step,nt_print2screen) == 0) then
            maxDiv = maxval(gp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            call message(0,"Time",gp%tsim)
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message(1,"Mean TKE:",gp%getMeanKE())
            call message_min_max(1,"Bounds for u:", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))
            if (allocated(gp%pressure)) then
               call message_min_max(1,"Bounds for P:", p_minval(minval(gp%pressure)), p_maxval(maxval(gp%pressure)))
            end if 
            !call message(1,"Max Error u:",p_maxval(maxval(abs(gp%u  - uexact))))
            !call message(1,"Max Error v:",p_maxval(maxval(abs(gp%v  - vexact))))
            !call message(1,"Max Error w:",p_maxval(maxval(abs(gp%wC - wexact))))
            !call message(1,"Max Error P:",p_maxval(maxval(abs(gp%pressure - pexact))))
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            call toc()
            call tic()
        end if 

    end subroutine


end module 
