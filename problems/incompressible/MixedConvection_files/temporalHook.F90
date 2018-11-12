module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval, p_sum
    use exits,              only: message, message_min_max
    !use EkmanDNS_IO,           only: output_tecplot!dumpData4Matlab 
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi
    use decomp_2d
    use reductions, only: p_sum

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

    subroutine doTemporalStuff(igp)
        class(igrid), intent(inout), target :: igp 
        real(rkind) :: wvar, uvar, kappa, E, tke
        real(rkind), dimension(:,:,:), pointer :: z

        if (mod(igp%step,nt_print2screen) == 0) then
            maxDiv = maxval(igp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            igp%rbuffxC(:,:,:,1) = igp%wC*igp%wC
            wvar = p_sum(sum(igp%rbuffxC(:,:,:,1)))/(igp%nx*igp%ny*igp%nz)
            
            ! get base state
            kappa = (igp%Ra**(1.d0/4.d0))/sqrt(2.d0)
            E = -2.d0*kappa*sqrt(igp%Ra)*(cosh(2*kappa)  + cos(2*kappa))/(sinh(2*kappa) - sin(2*kappa))
            z => igp%mesh(:,:,:,3)
            igp%rbuffxC(:,:,:,1) = -(E/sqrt(igp%Ra))*(sinh(kappa*(1+z))*sin(kappa*(1-z)) &
                                 & + sin(kappa*(1+z))*sinh(kappa*(1-z)))/(cosh(2*kappa) +  cos(2*kappa))
            igp%rbuffxC(:,:,:,1) = igp%u - igp%rbuffxC(:,:,:,1)
            igp%rbuffxC(:,:,:,1) = igp%rbuffxC(:,:,:,1)*igp%rbuffxC(:,:,:,1)
            uvar = p_sum(sum(igp%rbuffxC(:,:,:,1)))/(igp%nx*igp%ny*igp%nz)
            
            tke = 0.5d0*(uvar + wvar)

            call message(0,"Time",igp%tsim)
            call message(1,"TIDX:",igp%step)
            call message(1,"MaxDivergence:",DomMaxDiv)
            call message(1,"TKE:", tke)
            call message_min_max(1,"Bounds for u:", p_minval(minval(igp%u)),p_maxval(maxval(igp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(igp%v)),p_maxval(maxval(igp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(igp%w)),p_maxval(maxval(igp%w)))
            call message_min_max(1,"Bounds for T:", p_minval(minval(igp%T)),p_maxval(maxval(igp%T)))
            if (igp%useCFL) then
                call message(1,"Current dt:",igp%dt)
            end if
            call message("==========================================================")
            call toc()
            call tic()
        end if      
    end subroutine


end module 
