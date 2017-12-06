module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGrid, only: igrid
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max
    !use pbl_IO,           only: output_tecplot!dumpData4Matlab 
    use constants,          only: half
    use timer,              only: tic, toc 
    use mpi
    use decomp_2d

    implicit none 

    integer :: nt_print2screen = 20
    integer :: tid_statsDump = 5000
    integer :: tid_compStats = 100
    real(rkind) :: time_startDumping = 10.0_rkind, maxDiv, DomMaxDiv
    integer :: ierr 
    
    integer :: tid_start_planes = 1
    integer :: tid_stop_planes = 100000
    integer :: tid_dump_plane_every = 10000

contains

    subroutine doTemporalStuff(gp)
        class(igrid), intent(inout) :: gp 
      
        if (mod(gp%step,nt_print2screen) == 0) then
            maxDiv = maxval(gp%divergence)
            DomMaxDiv = p_maxval(maxDiv)
            call message(0,"Time",gp%tsim)
            if (gp%useSGS) call message(1,"u_star:",gp%sgsmodel%get_ustar())
            call message(1,"TIDX:",gp%step)
            call message(1,"MaxDiv:",DomMaxDiv)
            call message_min_max(1,"Bounds for u:", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message_min_max(1,"Bounds for v:", p_minval(minval(gp%v)), p_maxval(maxval(gp%v)))
            call message_min_max(1,"Bounds for w:", p_minval(minval(gp%w)), p_maxval(maxval(gp%w)))
            if (gp%initspinup) then
               call message_min_max(1,"Bounds for T:", p_minval(minval(gp%T)), p_maxval(maxval(gp%T)))
            end if
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if
            call message("==========================================================")
            call toc()
            call tic()
        end if 

        if (mod(gp%step,gp%t_dataDump) == 0) then
            call interpTau(gp)
        endif

    end subroutine


    !!! Added by Mike to interpolate tau13 and multiply and interpolate u'w'
    subroutine interpTau(this)
        class(igrid), intent(inout), target :: this
        complex(rkind), dimension(:,:,:), pointer :: ybuffC, zbuffC, zbuffE, ybuffE
        real(rkind), dimension(:,:,:), pointer :: rbuff, rbuffE

        zbuffE => this%cbuffzE(:,:,:,1)
        zbuffC => this%cbuffzC(:,:,:,1)
        ybuffC => this%cbuffyC(:,:,:,1)
        ybuffE => this%cbuffyE(:,:,:,1)
        rbuff => this%rbuffxC(:,:,:,1)
        rbuffE => this%rbuffxE(:,:,:,1)
        !rbuffv => this%rbuffxC(1,1,:,1)
        !rbuffEv => this%rbuffxE(1,1,:,1)

        call this%spectE%fft(this%tau13,ybuffE)
        call transpose_y_to_z(ybuffE,zbuffE,this%sp_gpE)
        call this%Pade6opZ%interpz_E2C(zbuffE,zbuffC, 0, 0)
        call transpose_z_to_y(zbuffC,ybuffC,this%sp_gpC)
        call this%spectC%ifft(ybuffC,rbuff)
        ! Dump the interpolated and non-interpolated tau13

        call this%dumpFullField(rbuff,'tauC')
        call this%dumpFullField(this%tau13,'tauE',this%gpE)

        ! Create u'w'
        call this%spectC%fft(this%u,ybuffC)
        call transpose_y_to_z(ybuffC,zbuffC,this%sp_gpC)
        call this%Pade6opZ%interpz_C2E(zbuffC,zbuffE,0,1)
        call transpose_z_to_y(zbuffE,ybuffE,this%sp_gpE)
        call this%spectE%ifft(ybuffE,rbuffE)
        rbuffE = rbuffE * this%w
        call this%spectE%fft(rbuffE,ybuffE)
        call transpose_y_to_z(ybuffE,zbuffE,this%sp_gpE)
        call this%Pade6opZ%interpz_E2C(zbuffE,zbuffC,0,0)
        call transpose_z_to_y(zbuffC,ybuffC,this%sp_gpC)
        call this%spectC%ifft(ybuffC,rbuff)
        ! Get spatial mean
        !rbuffv = zero
        !rbuffEv = zero
        !do i = 1, this%nx
        !   do j = 1, this%ny
        !       rbuffv = rbuffv + this%u(i,j,:)
        !       rbuffEv = rbuffEv + this%w(i,j,:)
        !   end do
        !end do 
        !call this%spect 
        ! Print u'w'
        call this%dumpFullField(rbuff,'upwp')




    end subroutine


end module 
