module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max, GracefulExit
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi
    use decomp_2d
    use reductions, only: p_sum

    implicit none 

    integer :: i, j, nt_print2screen = 1
    real(rkind) ::  maxDiv, DomMaxDiv, angle
    integer :: ierr 

contains
    subroutine initialize_controller_location(igp, filename)
        character(len=*), intent(in) :: filename
        class(igrid), intent(inout) :: igp
        real(rkind) :: beta, sigma, phi_ref, alpha, angleTrigger
        integer :: ioUnit,  controlType, z_ref = 16
        namelist /CONTROL/ beta, sigma, phi_ref, z_ref, alpha, controlType, angleTrigger

        ioUnit = 11
        open(unit=ioUnit, file=trim(filename), form='FORMATTED', iostat=ierr)
        read(unit=ioUnit, NML=CONTROL)
        close(ioUnit)
        
        igp%zHubindex = z_ref


    end subroutine 

    subroutine doTemporalStuff(igp)
        class(igrid), intent(inout) :: igp  
        real(rkind) :: speedTop, um, vm, speedHub, utop, vtop
        igp%rbuffxC(:,:,:,1) = atan2(igp%v, igp%u) !* 180.d0 / 3.14d0
        call transpose_x_to_y(igp%rbuffxC(:,:,:,1),igp%rbuffyC(:,:,:,1),igp%gpC)
        call transpose_y_to_z(igp%rbuffyC(:,:,:,1),igp%rbuffzC(:,:,:,1),igp%gpC)
        igp%angleHubHeight = p_sum(sum(igp%rbuffzC(:,:,igp%zHubIndex,1))) / &
                        (real(igp%gpC%xsz(1),rkind) * real(igp%gpC%ysz(2),rkind))
        call transpose_x_to_y(igp%u,igp%rbuffyC(:,:,:,1),igp%gpC)
        call transpose_y_to_z(igp%rbuffyC(:,:,:,1),igp%rbuffzC(:,:,:,1),igp%gpC)
        call transpose_x_to_y(igp%v,igp%rbuffyC(:,:,:,1),igp%gpC)
        call transpose_y_to_z(igp%rbuffyC(:,:,:,1),igp%rbuffzC(:,:,:,2),igp%gpC)
        utop = p_sum(sum(igp%rbuffzC(:,:,igp%gpC%zsz(3),1))) / &
                (real(igp%gpC%xsz(1),rkind) * real(igp%gpC%ysz(2),rkind))
        vtop = p_sum(sum(igp%rbuffzC(:,:,igp%gpC%zsz(3),2))) / &
                (real(igp%gpC%xsz(1),rkind) * real(igp%gpC%ysz(2),rkind))
        speedTop = (utop*utop + vtop*vtop)**0.5
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
            call message_min_max(1,"Bounds for u:", p_minval(minval(igp%u)), p_maxval(maxval(igp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(igp%v)), p_maxval(maxval(igp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(igp%w)), p_maxval(maxval(igp%w)))
            call message_min_max(1,"Bounds for T:", p_minval(minval(igp%T)), p_maxval(maxval(igp%T)))
            call message(1,"u_star:",igp%sgsmodel%get_ustar())
            call message(1,"Inv. Ob. Length:",igp%sgsmodel%get_InvObLength())
            call message(1,"wTh_surf:",igp%sgsmodel%get_wTh_surf())
            call message(1,"hub angle, degrees:",igp%angleHubHeight * 180.d0/3.14d0)
            call message(1,"frameAngle:",igp%frameAngle)
            call message(1,"Control w, rad/time:",igp%wFilt)
            call message(1,"Control Galpha:", igp%G_alpha)
            call message(1,"speed at the top:", speedTop)
            call message(1,"u at the top:", utop)
            call message(1,"v at the top:", vtop)
            call message(1,"speed at the hub:",speedHub)
            call message(1,"u at the hub:", um)
            call message(1,"v at the hub:", vm)
            call message(1,"Hub",igp%zHubIndex)
            if (igp%useCFL) then
                call message(1,"Current dt:",igp%dt)
            end if
            call message("==========================================================")
            call toc()
            call tic()
        end if 

    end subroutine


end module 
