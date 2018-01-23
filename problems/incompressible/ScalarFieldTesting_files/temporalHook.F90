module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max
    use constants,          only: half, one
    use timer,              only: tic, toc 
    use mpi
    use ScalarFieldTesting_parameters

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
   
         select case(direction) 
         case (1)
            scalarExact = cos(4.d0*(x - one*gp%tsim))*sin(2.d0*y)*cos(3.d0*z)  
         case (2)
            scalarExact = cos(4.d0*x)*sin(2.d0*(y - one*gp%tsim))*cos(3.d0*z)  
         case (3)
            scalarExact = cos(4.d0*x)*sin(2.d0*y)*cos(3.d0*(z - one*gp%tsim))  
         end select
        
         if (mod(gp%step,nt_print2screen) == 0) then
            maxDiv = maxval(gp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            call message(0,"Time",gp%tsim)
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message_min_max(1,"Bounds for u:", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))
            call message_min_max(1,"Bounds for S:", p_minval(minval(gp%scalars(1)%F)), p_maxval(maxval(gp%scalars(1)%F)))
            call message(1,"Max Error S:", p_maxval(maxval(abs(gp%scalars(1)%F - scalarExact))))
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            call toc()
            call tic()
        end if 

    end subroutine


end module 
