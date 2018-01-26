module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
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
        class(igrid), intent(inout) :: igp 
        real(rkind) :: speedTop, um, vm, speedHub
        igp%rbuffxC(:,:,:,1) = atan2(igp%v, igp%u) !* 180.d0 / 3.14d0
        call transpose_x_to_y(igp%rbuffxC(:,:,:,1),igp%rbuffyC(:,:,:,1),igp%gpC)
        call transpose_y_to_z(igp%rbuffyC(:,:,:,1),igp%rbuffzC(:,:,:,1),igp%gpC)
        igp%angleHubHeight = p_sum(sum(igp%rbuffzC(:,:,igp%zHubIndex,1))) / &
                        (real(igp%gpC%xsz(1),rkind) * real(igp%gpC%ysz(2),rkind))
        !igp%rbuffxC(:,:,:,1) = (igp%u*igp%u + igp%v * igp%v)**0.5
        !call transpose_x_to_y(igp%rbuffxC(:,:,:,1),igp%rbuffyC(:,:,:,1),igp%gpC)
        !call transpose_y_to_z(igp%rbuffyC(:,:,:,1),igp%rbuffzC(:,:,:,1),igp%gpC)
        !speed = p_sum(sum(igp%rbuffzC(:,:,igp%gpC%zsz(3),1))) / &
        !        (real(igp%gpC%xsz(1),rkind) * real(igp%gpC%ysz(2),rkind))

        call transpose_x_to_y(igp%u,igp%rbuffyC(:,:,:,1),igp%gpC)
        call transpose_y_to_z(igp%rbuffyC(:,:,:,1),igp%rbuffzC(:,:,:,1),igp%gpC) 
        call transpose_x_to_y(igp%v,igp%rbuffyC(:,:,:,1),igp%gpC)
        call transpose_y_to_z(igp%rbuffyC(:,:,:,1),igp%rbuffzC(:,:,:,2),igp%gpC)
        um = p_sum(sum(igp%rbuffzC(:,:,igp%gpC%zsz(3),1))) / &
                (real(igp%gpC%xsz(1),rkind) * real(igp%gpC%ysz(2),rkind))
        vm = p_sum(sum(igp%rbuffzC(:,:,igp%gpC%zsz(3),2))) / &
                (real(igp%gpC%xsz(1),rkind) * real(igp%gpC%ysz(2),rkind))
        speedTop = (um*um + vm*vm)**0.5
        um = p_sum(sum(igp%rbuffzC(:,:,igp%zHubIndex,1))) / &
                (real(igp%gpC%xsz(1),rkind) * real(igp%gpC%ysz(2),rkind))
        vm = p_sum(sum(igp%rbuffzC(:,:,igp%zHubIndex,2))) / &
                (real(igp%gpC%xsz(1),rkind) * real(igp%gpC%ysz(2),rkind))
        speedHub = (um*um + vm*vm)**0.5
        igp%angleHubHeight = atan2(vm,um)

        if (mod(igp%step,nt_print2screen) == 0) then
            maxDiv = maxval(igp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            call message(0,"Time",igp%tsim)
            call message(1,"TIDX:",igp%step)
            !call message(1,"MaxDiv:",DomMaxDiv)
            call message_min_max(1,"Bounds for u:", p_minval(minval(igp%u)),p_maxval(maxval(igp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(igp%v)),p_maxval(maxval(igp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(igp%w)),p_maxval(maxval(igp%w)))
            !call message_min_max(1,"Bounds for T:", p_minval(minval(igp%T)),p_maxval(maxval(igp%T)))
                       call message(1,"hub angle, degrees:",igp%angleHubHeight *180.d0/3.14d0)
            call message(1,"frameAngle:",igp%frameAngle)
            call message(1,"Control w, rad/time:",igp%wFilt)
            call message(1,"Control Galpha:", igp%G_alpha)
            call message(1,"speed at the top:", speedTop)
            call message(1,"u at the top:", um)
            call message(1,"v at the top:", vm)
            call message(1,"speed at the hub:",speedHub)
            call message(1,"Hub",igp%zHubIndex)
            call message(1,"Top",igp%gpC%zsz(3))
            if (igp%useCFL) then
                call message(1,"Current dt:",igp%dt)
            end if
            call message("==========================================================")
            call toc()
            call tic()
        end if      
    end subroutine


end module 
