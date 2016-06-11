module temporalHook
    use kind_parameters,    only: rkind
    use IncompressibleGridWallM, only: igridWallM
    use reductions,         only: P_MAXVAL, p_minval
    use exits,              only: message, message_min_max
    use gabls_IO,           only: dumpData4Matlab 
    use constants,          only: half
    use timer,              only: tic, toc 
    use decomp_2d,          only: nrank
    use mpi

    implicit none 

    integer :: nt_print2screen = 20
    integer :: tid_statsDump = 5000
    integer :: tid_compStats = 100
    real(rkind) :: time_startDumping = 10.0_rkind
    integer :: ierr 
    
    integer :: tid_start_planes = 1
    integer :: tid_stop_planes = 100000
    integer :: tid_dump_plane_every = 100

contains

    subroutine doTemporalStuff(gp)
        class(igridWallM), intent(inout) :: gp 
      
        if (mod(gp%step,nt_print2screen) == 0) then
            call message(0,"Time Step",gp%step)
            call message(0,"Time",gp%tsim)
            call message_min_max(1,"Bounds for u:", p_minval(minval(gp%u)), p_maxval(maxval(gp%u)))
            call message(1,"T_surf:", gp%Tsurf)
            call message(1,"u_star:",gp%ustar)
            call message(1,"Inv. Ob. Length:",gp%InvObLength)
            if (gp%useCFL) then
                call message(1,"Current dt:",gp%dt)
            end if 
            if (nrank == 0) then
                print*, "==================================="
            end if 
            call toc()
            call tic()
        end if 

        if (mod(gp%step,gp%t_dataDump)==0) then
           call message(0,"NOW PERFORMING DATA DUMP")
           call gp%dumpFullField(gp%u,'uVel')
           call gp%dumpFullField(gp%v,'vVel')
           call gp%dumpFullField(gp%wC,'wVel')
           call gp%dumpFullField(gp%T,'potT')
        end if 

    end subroutine


end module 
